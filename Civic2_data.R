#### LOAD PACKAGES and SET DIRECTORY 
library(corrplot)
library(viridis)
library(tidyverse)
library(tidycensus) 
library(sf)
library(terra)
library(stars)
library(raster) #make sure ncdf4 package is installed 
library(lubridate)

setwd("~/Documents/01_VECTOR.nosync/Sandia")

##################################################################################################
#### CREATE CENSUS MAP ###########################################################################
##################################################################################################
mycrs = 5070 #chose projected coordinate system: EPSG 5070 NAD83 Conus Albers
year=2019 # year for county boundaries 

# Get map
options(tigris_use_cache = TRUE) #cache shapefiles for future sessions
state = "Texas"
county = c("Harris County")
census_map = get_acs(geography = "tract", state = state, county = county,
                     variables=c("B01003_001"), year = year, geometry = TRUE, 
                     cache_table = TRUE)
census_map = census_map %>%
  mutate(POPULATION = estimate) %>%
  dplyr::select(GEOID, NAME, POPULATION) %>% 
  st_transform(mycrs) # project to Conic Equal Area Albers, EPSG:5070 

# Calculate area and population density of each county 
census_map_area = census_map %>%  
  mutate(AREA = as.vector(st_area(census_map))) %>% #sq-meters; as.vector removes units suffix 
  mutate(DENSITY = POPULATION / AREA * 1000^2) #population per sq-km 

################################################################################################
#### STATIC ENVIRONMENTAL FACTORS ##############################################################
################################################################################################
## NLCD
nlcd_path = "./Data/nlcd_2019_land_cover_l48_20210604/nlcd_2019_land_cover_l48_20210604.img"
nlcd_ras = terra::rast(nlcd_path) #read in with terra package
nlcd_ras_100 = terra::aggregate(nlcd_ras, fact = 100/30, fun = "modal") #downsample to 100m resolution
nlcd_levels_100 = levels(nlcd_ras_100) #get levels
nlcd_proj_100 = terra::project(nlcd_ras_100, paste("EPSG:", mycrs)) #project raster to desired crs
nlcd_proj_100_crop = crop(nlcd_proj_100, vect(census_map_area)) #crop
nlcd_proj_100_mask = mask(nlcd_proj_100_crop, vect(census_map_area)) #mask
nlcd_proj_100_extract = terra::extract(x = nlcd_proj_100_mask, y = vect(census_map_area)) #extract raster values
nlcd_stack  = nlcd_proj_100_extract %>% # Assign land class based on NLCD land type then group by ID and NLCD, county number 
  mutate(LANDCLASS = case_when( # https://www.mrlc.gov/data/legends/national-land-cover-database-2011-nlcd2011-legend
    `NLCD Land Cover Class` %in% c("Open Water", "Perennial Snow/Ice") ~ "Water",
    `NLCD Land Cover Class` %in% c("Developed, Open Space", "Developed, Low Intensity", "Developed, Medium Intensity", "Developed, High Intensity") ~ "Developed",
    `NLCD Land Cover Class` %in% c("Barren Land") ~ "Barren",
    `NLCD Land Cover Class` %in% c("Deciduous Forest", "Evergreen Forest", "Mixed Forest") ~ "Forest",
    `NLCD Land Cover Class` %in% c("Shrub/Scrub") ~ "Shrub",
    `NLCD Land Cover Class` %in% c("Herbaceous")~ "Herbaceous",
    `NLCD Land Cover Class` %in% c("Hay/Pasture", "Cultivated Crops") ~ "Cultivated",
    `NLCD Land Cover Class` %in% c("Woody Wetlands", "Emergent Herbaceous Wetlands") ~ "Wetlands",
    TRUE ~ "Other"
  )) %>%
  group_by(ID, LANDCLASS) %>%
  count()
nlcd_wide = nlcd_stack %>% # Expand LANDCLASS into wide format 
  pivot_wider(names_from = LANDCLASS, values_from = n)
land_names = colnames(nlcd_wide[2:ncol(nlcd_wide)]) #get column names
nlcd_wide2 = nlcd_wide %>% 
  replace(is.na(.), 0) %>%
  mutate(Total = rowSums(across(all_of(land_names)))) # Column totals of land type 
nlcd_perc = nlcd_wide2 %>%
  mutate(across(all_of(land_names), ~ .x / Total)) # Percentages of land type
census_map_nlcd = census_map_area %>% # join to census shapefile
  bind_cols(nlcd_perc)
#save(census_map_nlcd, file = "./Data/census_map_nlcd.Rda")
#load(file = "./Data/census_map_nlcd.Rda")

## DEM
dem_se_path = "./Data/gt30w100n40_dem/gt30w100n40.dem" 
dem_ne_path = "./Data/gt30w100n90_dem/gt30w100n90.dem"
dem_sw_path = "./Data/gt30w140n40_dem/gt30w140n40.dem"
dem_nw_path = "./Data/gt30w140n90_dem/gt30w140n90.dem"
dem_se = rast(dem_se_path) #read in raster
dem_ne = rast(dem_ne_path)
dem_sw = rast(dem_sw_path)
dem_nw = rast(dem_nw_path)
dem_merge1 = terra::merge(dem_se, dem_ne) #merge rasters
dem_merge2 = terra::merge(dem_merge1, dem_sw)
dem_merge3 = terra::merge(dem_merge2, dem_nw)
dem_proj = terra::project(dem_merge3, paste("EPSG:", mycrs)) # re-project to crs
dem_crop = crop(dem_proj, vect(census_map_area)) #crop
dem_final = mask(dem_crop, vect(census_map_area)) #mask
dem_extract = terra::extract(x = dem_final, y = vect(census_map_area)) #extract
dem_extract_group = dem_extract %>% #group by ID
  group_by(ID) %>%
  summarize(DEM_mean = mean(gt30w100n40, na.rm = TRUE),
            DEM_med = median(gt30w100n40, na.rm = TRUE),
            DEM_sd = sd(gt30w100n40, na.rm = TRUE),
            DEM_min = min(gt30w100n40, na.rm = TRUE),
            DEM_max = max(gt30w100n40, na.rm = TRUE)
  )
census_map_dem = census_map_area %>% #join DEM to counties
  bind_cols(dem_extract_group)
#save(census_map_dem, file = "./Data/census_map_dem.Rda")
#load(file = "./Data/census_map_dem.Rda")

## ROOT ZONE
tx_gdb = "./Data/July2020_gSSURGO_by_State/gSSURGO_TX/gSSURGO_TX.gdb"
tx_sf = sf::st_read(dsn = tx_gdb, layer = "MUPOLYGON") # read in MUPOLYGON layer 
tx_group = tx_sf %>% group_by(MUKEY) %>% summarise(n = n()) #group by MUKEY 
tx_val = sf::st_read(dsn = tx_gdb, layer = "Valu1") # read in Value1 table
tx_join = tx_group %>% left_join(tx_val, by = c("MUKEY" = "mukey")) # join values
tx_ras = st_rasterize(tx_join["rootznemc"], dx = 100, dy = 100) #100m x 100km resolution raster
tx_write = paste("TX", "_rootznemc100.tif", sep = "") #name of raster
write_stars(tx_ras, dsn = tx_write) # write raster in working directory 
tx_terra = rast(tx_ras)
tx_proj = terra::project(tx_terra, paste("EPSG:", mycrs)) #project to crs
tx_crop = crop(tx_proj, vect(census_map_area)) #crop
tx_mask = mask(tx_crop, vect(census_map_area)) #mask
tx_extract = terra::extract(x = tx_mask, y = vect(census_map_area)) # extract values
mode = function(x) { # function for mode
  ux = na.omit(unique(x) )
  tab = tabulate(match(x, ux)); ux[tab == max(tab)]
}
tx_extract_group = tx_extract %>% # group extracted data 
  group_by(ID) %>%
  summarize(RZ_mean = mean(filed5f1314bf64, na.rm = TRUE),
            RZ_med = median(filed5f1314bf64, na.rm = TRUE),
            RZ_mode = mode(filed5f1314bf64)
  )
census_map_rz = census_map_area %>% #join rz values
  bind_cols(tx_extract_group)
#save(census_map_rz, file = "./Data/census_map_rz.Rda")
#load(file = "./Data/census_map_rz.Rda")

## SOCIO-ECONOMIC (from feature selection) 
tx_vars_use = c("B01003_001", #total population
                "B01001_026", #female population
                "B01001C_001", #native american population
                "B25024_010", "B25001_001", #mobile homes and total housing units 
                "B25077_001", #median housing value (Factor 1)
                "B01001D_001", #asian population (part of Factor 1)
                "B17001_002", #poverty population (Factor 2) 
                "B01001B_001", #black population (part of Factor 2)
                "B19055_002", "B19055_001", #households receiving social security and total households (Factor 4) 
                "B01001I_001", #hispanic/latino population (Factor 5)
                "B06009_004", "B06009_005", "B06009_006", "B06009_002", #education equity: -abs((B06009_004 + B06009_005 + B06009_006) - B06009_002)
                "B25034_007", "B25034_008", "B25034_009", "B25034_001", #housing stock construction quality: (B25034_007 + B25034_008 + B25034_009) / B25034_001
                "B05007_002", #%population not foreign-born persons who came to US within previous five years: (B01003_001 - B05007_002) / B01003_001
                "B06001_013", #%population born in state of residence B06001_013 / B01003_001
                "B06009_002", #%population less than 12-th grade degree B06009_002 / B01003_001
                "B23025_004" #%population in labor force: B23025_004 / B01003_001 
                )
tx_acs_data = get_acs(geography = "tract", state = state, county = county, variables=tx_vars_use, year = year, geometry = FALSE)
tx_acs_data_w = tx_acs_data %>%
  dplyr::select(-moe, -NAME) %>% #remove errors of estimates because we're spreading by variable and don't want duplicates, and don't need name b/c joining by GEOID
  spread(key=variable, value = estimate)
tx_sec_dat = tx_acs_data_w %>%
  mutate(QFEMALE = B01001_026 / B01003_001) %>% # %female
  mutate(QNATIVE = B01001C_001 / B01003_001) %>% # %native american 
  mutate(QMOHO = B25024_010 / B25001_001) %>% # %mobile homes
  mutate(MDHVAL = B25077_001) %>% # median housing value 
  mutate(QASIAN = B01001D_001 / B01003_001) %>% # %asian 
  mutate(QPOVERTY = B17001_002 / B01003_001) %>% # %below poverty line
  mutate(QBLACK = B01001B_001 / B01003_001) %>% # %black 
  mutate(QSSBEN =  B19055_002 / B19055_001) %>% #%households receiving social security 
  mutate(QSPANISH = B01001I_001 / B01003_001) %>% # %hispanic or latino 
  mutate(EDEQUITY = -abs((B06009_004 + B06009_005 + B06009_006) - B06009_002)) %>% #education equity 
  mutate(HOUSEQUAL = (B25034_007 + B25034_008 + B25034_009) / B25034_001) %>% #housing quality 
  mutate(PATTACHIM = (B01003_001 - B05007_002) / B01003_001) %>% #long-time immigrants
  mutate(PATTACHRES = B06001_013 / B01003_001) %>% #native state residents 
  mutate(QED12 =  B06009_002 / B01003_001) %>% #percent less than 12th grade education  
  mutate(QEMPL = B23025_004 / B01003_001) %>% #percent employed
  dplyr::select(c(GEOID, QFEMALE:QEMPL))
census_map_social = census_map_area %>%
  inner_join(tx_sec_dat, by = "GEOID")
#save(census_map_social, file = "./Data/census_map_social.Rda")
#load(file = "./Data/census_map_social.Rda")

#########################################################################################################
#### DYNAMIC ENVIRONMENTAL FACTORS ######################################################################
#########################################################################################################
## SPI - monthly lagged variable 
spi03_name = "./Data/SPI_ncdf/nclimgrid-spi-pearson-03.nc"
spi03_brick = raster::brick(spi03_name) #read in ncdf4 file as raster brick
spi03_subset = raster::subset(spi03_brick, dim(spi03_brick)[3]) #subset most recent layer 
spi03_proj = raster::projectRaster(spi03_subset, crs = crs(census_map_area)) #project raster into my crs
spi03_spat = terra::rast(spi03_proj) #convert to SpatRaster (makes for easier projection and extraction) 
spi03_crop = terra::crop(spi03_spat, vect(census_map_area)) #crop to map
spi03_mask = terra::mask(spi03_crop, vect(census_map_area)) #mask to map (make anything outside NA)
spi03_ext = terra::extract(spi03_mask, y = vect(census_map_area)) #extract SPI layer for each tract 
spi03_ext_g = spi03_ext %>% #group by ID 
  group_by(ID) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
colnames(spi03_ext_g) = c("ID","spi03_mean")
spi12_name = "./Data/SPI_ncdf/nclimgrid-spi-pearson-12.nc"
spi12_brick = raster::brick(spi12_name) #read in ncdf4 file as raster brick
spi12_subset = raster::subset(spi12_brick, dim(spi12_brick)[3]) #subset most recent layer 
spi12_proj = raster::projectRaster(spi12_subset, crs = crs(census_map_area)) #project raster into my crs
spi12_spat = terra::rast(spi12_proj) #convert to SpatRaster (makes for easier projection and extraction) 
spi12_crop = terra::crop(spi12_spat, vect(census_map_area)) #crop to map
spi12_mask = terra::mask(spi12_crop, vect(census_map_area)) #mask to map (make anything outside NA)
spi12_ext = terra::extract(spi12_mask, y = vect(census_map_area)) #extract SPI layer for each tract 
spi12_ext_g = spi12_ext %>% #group by ID 
  group_by(ID) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
colnames(spi12_ext_g) = c("ID","spi12_mean")
spi24_name = "./Data/SPI_ncdf/nclimgrid-spi-pearson-24.nc"
spi24_brick = raster::brick(spi24_name) #read in ncdf4 file as raster brick
spi24_subset = raster::subset(spi24_brick, dim(spi24_brick)[3]) #subset most recent layer 
spi24_proj = raster::projectRaster(spi24_subset, crs = crs(census_map_area)) #project raster into my crs
spi24_spat = terra::rast(spi24_proj) #convert to SpatRaster (makes for easier projection and extraction) 
spi24_crop = terra::crop(spi24_spat, vect(census_map_area)) #crop to map
spi24_mask = terra::mask(spi24_crop, vect(census_map_area)) #mask to map (make anything outside NA)
spi24_ext = terra::extract(spi24_mask, y = vect(census_map_area)) #extract SPI layer for each tract 
spi24_ext_g = spi24_ext %>% #group by ID 
  group_by(ID) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
colnames(spi24_ext_g) = c("ID","spi24_mean")
spi_stack = spi03_ext_g %>%
  inner_join(spi12_ext_g, by = "ID") %>%
  inner_join(spi24_ext_g, by = "ID")
census_map_spi = census_map_area %>%
  bind_cols(spi_stack)
#save(census_map_spi, file = "./Data/census_map_spi.Rda")
#load(file = "./Data/census_map_spi.Rda")

# ## SOIL MOISTURE - 3 day lag 
# soil.names = list.files(path = "./Data/Soil_current", full.names = T, pattern = "\\.grb$") # get most recent day of hourly forecasts (3 day lag)
# for (i in 1:length(soil.names)){
#   temp_i = soil.names[i]
#   print(temp_i)
#   temp_rast = raster::brick(temp_i) #read in via raster brick
#   temp_proj = raster::projectRaster(temp_rast, crs = crs(census_map_area)) #project into crs
#   temp_terra = terra::rast(temp_proj) #convert to Spat Raster for easier masking and extraction
#   temp_crop = terra::crop(temp_terra, vect(census_map_area)) #crop to map
#   temp_mask = terra::mask(temp_crop, vect(census_map_area)) #mask to map (make anything outside NA)
#   temp_s26 = temp_mask[[26]] # soil moisture 0-10cm (kg/m^2)
#   temp_s27 = temp_mask[[27]] # soil moisture 10-40cm (kg/m^2)
#   temp_s28 = temp_mask[[28]] # soil moisture 40-100cm (kg/m^2)
#   # Extract soil layer values by tract 
#   temp_extract26 = terra::extract(x = temp_s26, y = vect(census_map_area)) %>% 
#     group_by(ID) %>%
#     summarise(across(everything(), mean, na.rm = TRUE))
#   temp_extract27 = terra::extract(x = temp_s27, y = vect(census_map_area)) %>%
#     group_by(ID) %>%
#     summarise(across(everything(), mean, na.rm = TRUE))
#   temp_extract28 = terra::extract(x = temp_s28, y = vect(census_map_area)) %>%
#     group_by(ID) %>%
#     summarise(across(everything(), mean, na.rm = TRUE))
#   # Join temp extractions by ID
#   temp_extract = temp_extract26 %>%
#     inner_join(temp_extract27, by = c("ID")) %>%
#     inner_join(temp_extract28, by = c("ID"))
#   colnames(temp_extract) = c("ID", "soil0_10", "soil10_40", "soil40_100")
#   # Add date_hour (character) variable whose value is the file name, formatted to match outage data ("2018-01-01T00:00:00Z")
#   temp_name1 = sub(".002.grb", "", sub(".*NLDAS_NOAH0125_H.A", "", temp_i)) #"YYYYMMDD.HHMM"
#   temp_name2 = paste(substr(temp_name1, 1, 4), "-", substr(temp_name1, 5, 6), "-", substr(temp_name1, 7, 8), "T", substr(temp_name1, 10, 11), ":", substr(temp_name1, 12, 13), ":00Z", sep = "")
#   temp_day = substr(temp_name2, 1,10)
#   temp_out = temp_extract
#   temp_out$date_hour = temp_name2
#   temp_out$date_day = temp_day
#   # Data-frame of GEOID, date-time (hourly and daily), and 3 soil moisture layer variables 
#   temp_map_soil = census_map_area %>%
#     bind_cols(temp_out) %>%
#     dplyr::select(-ID, -POPULATION) %>%
#     st_set_geometry(NULL)
#   if (i == 1) {
#     census_map_soil = temp_map_soil
#   }
#   else {
#     census_map_soil = census_map_soil %>%
#       bind_rows(temp_map_soil)
#   }
# }
# census_map_soil_day = census_map_soil %>%  
#   group_by(GEOID, date_day) %>%
#   summarise(soil10_3dLAG = mean(soil0_10, na.rm = T), 
#             soil40_3dLAG = mean(soil10_40, na.rm = T),
#             soil100_3dLAG = mean(soil40_100, na.rm = T)
#             )
# #save(census_map_soil_day, file = "./Data/census_map_soil_day.Rda")
# #load(file = "./Data/census_map_soil_day.Rda")

## WIND SPEED - forecasts from StormGeo
hurricane.frcst.files = list.files(path = "./Data/StormGEO/HarveyNC", full.names = T, pattern = "\\.nc$") #get forecast (NCDF4 format) 
hurricane.frcst.raster = lapply(hurricane.frcst.files, function(i){raster::raster(i, varname = "wspd")}) #read list of files as rasters
hurricane.frcst.proj = lapply(hurricane.frcst.raster, function(i){raster::projectRaster(i, crs = crs(census_map_area))}) #project into crs  
hurricane.frcst.terra = lapply(hurricane.frcst.proj, function(i){terra::rast(i)}) #convert to Spat Raster for easier masking and extraction
hurricane.frcst.crop = lapply(hurricane.frcst.terra, 
                              function(i){
                                if("try-error" %in% class(try(terra::crop(i, vect(census_map_area)), silent = T))){NULL}
                                else {terra::crop(i, vect(census_map_area))}
                                }
                              ) #crop to map and assign NULL for when the extents do not overlap 
hurricane.frcst.mask = lapply(hurricane.frcst.crop, 
                              function(i){
                                if("try-error" %in% class(try(terra::mask(i, vect(census_map_area)), silent = T))){NULL}
                                else {terra::mask(i, vect(census_map_area))}
                                }
                              ) #mask map and assign NULL for when the extents do not overlap 
hurricane.frcst.extract = lapply(hurricane.frcst.mask, 
                              function(i){
                                if("try-error" %in% class(try(terra::extract(i, vect(census_map_area)), silent = T))){NULL}
                                else {terra::extract(i, vect(census_map_area))}
                                }
                              ) #extract values to census tracts
hurricane.frcst.group = lapply(hurricane.frcst.extract, 
                                function(i){
                                  if(is.null(i)){rep(NA, nrow(census_map))} #want same rows as number of tracts even if null to make unlisting easier
                                  else {i %>% group_by(ID) %>% summarize(WIND_mean = mean(X10.meter.Windspeed, na.rm = TRUE))}
                                }
                              ) #group extracted values by census tract
hurricane.frcst.df = as.data.frame(hurricane.frcst.group) %>% dplyr::select(ID, contains("WIND_mean")) #put list into data.frame
hurricane.frcst.maps = census_map_area %>% bind_cols(hurricane.frcst.df) #map of forceasted wind speeds
max_wind_names = hurricane.frcst.df %>% dplyr::select(contains("WIND_mean")) %>% names()
census_map_WINDmax = hurricane.frcst.maps %>%
  mutate(WIND_max = pmax(!!!rlang::syms(max_wind_names), na.rm = T)) %>% # https://stackoverflow.com/questions/32978458/dplyr-mutate-rowwise-max-of-range-of-columns
  dplyr::select(GEOID, WIND_max)
#save(census_map_WINDmax, file = "./Data/census_map_WINDmax.Rda")
#load(file = "./Data/census_map_WINDmax.Rda")  


gg2 = ggplot(hurricane.frcst.max)+
  geom_sf(aes(fill = WIND_max), color = NA) + 
  #scale_fill_viridis_c(option="plasma", na.value = "grey50") +
  scale_fill_viridis_c(option="plasma", na.value = "skyblue") +
  #geom_sf(fill = NA, show.legend = F, color = "black", lwd = 0.005)+
  #coord_sf(datum = NA) + #removes gridlines 
  #guides(fill = "none") + #removes legend
  #theme_minimal()  #removes background
  theme_dark() +
  #labs(title = "NLCD Land Use", fill = "Developed\n(prop.)") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
  )
gg2


# ############ COUNTY HURICANE 
# county_map = get_acs(geography = "county", state = state,
#                      variables=c("B01003_001"), year = year, geometry = TRUE, 
#                      cache_table = TRUE)
# county_map = county_map %>%
#   mutate(POPULATION = estimate) %>%
#   dplyr::select(GEOID, NAME, POPULATION) %>% 
#   st_transform(mycrs) # project to Conic Equal Area Albers, EPSG:5070 
# 
# county_map_area = county_map %>%  
#   mutate(AREA = as.vector(st_area(county_map))) %>% #sq-meters; as.vector removes units suffix 
#   mutate(DENSITY = POPULATION / AREA * 1000^2) #population per sq-km 
# 
# i = 60
# csd = raster::raster(hurricane.frcst.files[i], varname = "wspd")
# csd2 = raster::projectRaster(csd, crs = crs(county_map_area))
# csd3 = terra::rast(csd2) #convert to Spat Raster for easier masking and extraction
# csd4 = terra::crop(csd3, vect(county_map_area)) #crop to map
# csd5 = terra::mask(csd4, vect(county_map_area)) #mask to map (make anything outside NA)
# csd6 = terra::extract(x = csd5, y = vect(county_map_area)) # extract values
# csd_group = csd6 %>% # group extracted data 
#   group_by(ID) %>%
#   summarize(WIND_mean = mean(X10.meter.Windspeed, na.rm = TRUE))
# county_map_wind = county_map_area %>% #join rz values
#   bind_cols(csd_group)
# 
# county_map_Harris = county_map_area %>%
#   dplyr::filter(GEOID == "48201")
# rast_df = as.data.frame(csd3, xy= TRUE)
# title = substr(hurricane.frcst.files[i], 26, 100)
# 
# gg3 = ggplot() +
#   geom_tile(data = rast_df, aes(x = x, y = y, fill = X10.meter.Windspeed), alpha = 0.9) +
#   geom_sf(data = county_map_area, color = "black", fill = "NA", lwd = 0.4) +
#   geom_sf(data = county_map_Harris, color = "black", fill = "Skyblue", lwd = 0.4) +
#   scale_fill_viridis_c(option="plasma", na.value = "grey10") +
#   theme_dark() +
#   ggtitle(paste(title)) + 
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank()
#   )
# gg3

