#### LOAD PACKAGES and SET DIRECTORY 
options(java.parameters = "-Xmx7g")
library(bartMachine)
set_bart_machine_num_cores(4)
library(lme4)
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

##############################################################################################
## WHICH VARIABLES TO USE? ###################################################################
##############################################################################################
#### LOAD COUNTY VARIABLES
load("./Data/SandiaVariables.Rda") #from Sandia2_vars.R

#### POWER DATA
#outages_csv = read.csv("./Data/SE_states_outage_merra_2018.csv", header = T)
#save(outages_csv, file = "./Data/outages.Rda")
load(file = "Data/outages.Rda")

# Format FIPS to character and include starting zero
outages_csv$fips_code = as.character(outages_csv$fips_code)
fips = outages_csv$fips_code
fips0 = str_pad(fips, 5, pad = "0")
outages_csv$fips_code = fips0

# Add day and month to data
county_map_outages = outages_csv %>%
  mutate(rowID = row_number()) %>% # add row number to index events 
  mutate(date_day = str_extract(date_hour, "^.{10}")) %>% #returns first 10 characters of date_time (yyyy-mm-dd)
  mutate(date_month = str_extract(date_hour, "^.{7}")) # returns yyyy-mm

#### PROCESS POWER DATA
# Filter dates that are flagged as outages 
county_map_outages_FILTER = county_map_outages %>%
  filter(outage_status %in% c("pre", "start", "during", "end")) #filter events flagged as outages

# # Join power outages by GEOID and time-stamp 
# county_map_outages_JOIN = county_map_area %>%
#   inner_join(county_map_outages_FILTER, by = c("GEOID" = "fips_code")) %>% # join outage/hourly weather data 
#   left_join(county_map_static_CLEAN, by = c("GEOID")) %>% #join static environmental variables
#   left_join(county_scores_CLEAN, by = c("GEOID")) %>% #join socio-economic variables  
#   left_join(county_map_soil_CLEAN, by = c("GEOID", "date_hour")) %>% #join soil moisture by county and hour time-stamp 
#   left_join(county_map_spi_CLEAN, by = c("GEOID", "date_month")) #join SPI by GEOID and month 

# Daily soil - 3 day lag for forecasting 
county_map_soil_CLEAN$date_day = substr(county_map_soil_CLEAN$date_hour, 1,10)
daily_soil = county_map_soil_CLEAN %>%
  group_by(GEOID, date_day) %>%
  summarise(soil10_DAY = mean(soil0_10),
            soil100_DAY = mean(soil40_100))
daily_soil$ts_day = lubridate::ymd(daily_soil$date_day)
daily_soil$ts_3daylag = daily_soil$ts_day + 3 #take 3-day lag (what's available for forecasting)
daily_soil_CLEAN = daily_soil %>% 
  mutate(date_3daylag = as.character(ts_3daylag)) %>% # back to character format to join later on
  dplyr::select(-date_day, -ts_day, -ts_3daylag) %>%
  rename(soil10_3dLAG = soil10_DAY, soil100_3dLAG = soil100_DAY) 

# Monthly SPI - 1 month lag for forecasting 
monthly_spi = county_map_spi_CLEAN
monthly_spi$ts_month = lubridate::ym(monthly_spi$date_month)
monthly_spi$ts_1mthlag = monthly_spi$ts_month %m+% months(1) 
monthly_spi_CLEAN = monthly_spi %>%
  mutate(date_1mthlag = substr(as.character(ts_1mthlag), 1, 7)) %>%
  dplyr::select(-date_month, -ts_month, -ts_1mthlag) %>%
  rename(spi03_1mLAG = spi_03, spi12_1mLAG = spi_12, spi24_1mLAG = spi_24) 

county_map_outages_CIVIC = county_map_area %>%
  inner_join(county_map_outages_FILTER, by = c("GEOID" = "fips_code")) %>% # join outage/hourly weather data
  left_join(county_map_static_CLEAN, by = c("GEOID")) %>% #join static environmental variables
  left_join(county_scores_CLEAN, by = c("GEOID")) %>% #join socio-economic variables
  inner_join(daily_soil_CLEAN, by = c("GEOID", "date_day" = "date_3daylag")) %>% #join soil moisture by county and 3day lag
  inner_join(monthly_spi_CLEAN, by = c("GEOID", "date_month" = "date_1mthlag")) #join SPI by GEOID and month

# Group by event and aggregate dynamic variables 
outages_group = county_map_outages_CIVIC %>%
  st_set_geometry(NULL) %>%
  group_by(outage_number, GEOID) %>%
  summarise(out_hrs= n(), out_maxcust = max(hr_mean_customers_out), out_percust = sum(hr_mean_customers_out) / sum(POPULATION),
            Density = mean(DENSITY), 
            PS_mean = mean(PS), PS_sd = sd(PS), #mean and sd of surface pressure
            SLP_mean = mean(SLP), PS_sd = sd(SLP), # same for sea level pressure
            QV_max = max(QV10M), # max specific humidity
            U_max = max(U10M),  
            V_max = max(V10M),
            WIND_max = max(WIND10M), 
            T_mean = mean(T10M), T_sd = sd(T10M), 
            TQI_mean = mean(TQI), TQL_mean = mean(TQL), TQV_mean = mean(TQV),
            soil10_lag = first(soil10_3dLAG), soil100_lag = first(soil100_3dLAG),
            spi03_lag = first(spi03_1mLAG), spi12_lag = first(spi12_1mLAG), spi24_lag = first(spi24_1mLAG)
  )

# Join back in the static environmental and socio-economic variables
county_outages_GROUP = outages_group %>%
  left_join(county_map_static_CLEAN, by = c("GEOID")) %>% #join static environmental variables
  left_join(county_scores_CLEAN, by = c("GEOID")) #join socio-economic variables 

# Select only data we'll have available for prediction at the census tract level  
county_outages_GROUP_CLEAN = county_outages_GROUP %>%
  dplyr::select(-c(PS_mean:TQV_mean, AMBULANCES:DEATHS, EMPBLDG:HOSPBEDS, INTERNET:NURSHOMES, PHYSICIANS, PROFORGS:PROXCAP, RADIO:WATEFF, QEXTRCT, QNRRES)) %>%
  dplyr::select(-c(RZ_mode, QGROUPHSE, QSERVIND)) #take out variables where effect has high uncertaintly, likely from lingering multicollinearity

##############################################################################################################
#### MODELING ################################################################################################
##############################################################################################################
# Filter to large events (90th percentile in duration, > 12 hrs)
df_bart = data.frame(county_outages_GROUP_CLEAN) %>%
  drop_na %>% 
  dplyr::filter(out_hrs > 12) #filter to big events 
# y = log(df_bart$out_hrs) #take log to help deal with extreme values 
y = log(df_bart$out_maxcust)
X = df_bart %>%
  dplyr::select(-outage_number, -GEOID, -out_hrs, -out_maxcust, -out_percust)
model_name = paste("Predicting Max Outages - No Weather Data")

bart = bartMachine(X, y) 
predictions = predict(bart, X)
CI = round(calc_credible_intervals(bart, X, ci_conf = 0.95), 2)
gg = dplyr::tibble(predictions = predictions,
                   lower = CI[,1],
                   upper = CI[,2],
                   actual = y
)
gg = arrange(gg, actual)
gg$index = seq.int(nrow(gg))
gg$Ymean = mean(gg$actual, na.rm = T)

rmse = sqrt(mean((gg$predictions - gg$actual)^2)) 
rsq = 1 - sum((gg$actual - gg$predictions)^2) / sum((gg$actual - mean(gg$actual))^2) 
lb1 = paste("R^2 == ", round(rsq, 3))


plot_filtering_estimates2 <- function(df) {
  p <- ggplot(data = gg, aes(x = index)) +
    theme_classic() +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = "indianred"), alpha = 0.5) + #credible intervals 
    geom_line(aes(y = predictions, colour = "red"), size = 0.25, alpha = .9) + #prediction point estimate
    geom_point(aes(y = actual, colour = "black"), size = 0.55, shape = 16, alpha = 0.9) + #actual observation points
    geom_line(aes(y = Ymean, colour = "blue"), size = 0.55, lty = "solid", alpha = 0.9) + #null model (mean only) 
    #ylab("Outage Duration (log hours)") + 
    ylab("Max Cust. Outages (log)") + 
    scale_y_continuous(labels = function(x) paste0(x)) +
    xlab("Outage Index (event x county)") +
    ggtitle(model_name) +
    scale_fill_identity(name = "fill", guide = 'legend', labels = c('95% CI')) + 
    scale_colour_manual(name = 'colour',
                        values = c('black',"blue", "red"),
                        labels = c('Actual', 'Mean', 'BART'),
                        guide = guide_legend(
                          reverse = T,
                          override.aes = list(
                            linetype = c("solid", "solid","blank"),
                            shape = c(NA, NA, 16))
                        )) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(),
          legend.direction = "vertical",
          legend.box = "horizontal",
          legend.position = c(.25, .75) #x and y percentage
    ) +
    annotate("text", x = quantile(gg$index, 0.8), y = quantile(gg$actual, .05), label=lb1, parse=T, color="red", size = 3) 
  print(p)
}
plot_filtering_estimates2(gg)

# Feature selection (Kapelner and Bleich, 2016)
vs = var_selection_by_permute(bart, bottom_margin = 10, num_reps_for_avg = 20, num_permute_samples = 20, plot = T)
vs$important_vars_local_names
vs$important_vars_global_se_names
# > vs$important_vars_local_names
# [1] "spi24_lag" "Density"   "Developed" "QMOHO"     "QFEMALE"   "QNATIVE"  
# [7] "RZ_mean"  
# > vs$important_vars_global_se_names
# [1] "spi24_lag" "Density"   "Developed"

#feature selection
X = df_bart %>%
  dplyr::select(-outage_number, -GEOID, -out_hrs, -out_maxcust, -out_percust) %>%
  dplyr::select(c("spi24_lag","Density", "Developed", "QMOHO", "QFEMALE", "QNATIVE", "RZ_mean"))
model_name = paste("Predicting Max Outages - No Weather Data")
###############################################################################################
###############################################################################################




#### CREATE CENSUS MAP
mycrs = 5070 #chose projected coordinate system: EPSG 5070 NAD83 Conus Albers
year=2019 # year for county boundaries 

# Get map
options(tigris_use_cache = TRUE) #cache shapefiles for future sessions
soptions(tigris_use_cache = TRUE) #to cache shapefiles for future sessions
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

#### STATIC ENVIRONMENTAL FACTORS
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

#### DYNAMIC ENVIRONMENTAL FACTORS
## SPI 
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





gg2 = ggplot(census_map_spi)+
  geom_sf(aes(fill = spi24_mean), color = NA) + 
  scale_fill_viridis_c(option="plasma", na.value = "grey50") +
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

