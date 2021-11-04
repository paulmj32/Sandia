## NLCD land use (30m resolution) 
# https://www.mrlc.gov/data/nlcd-2019-land-cover-conus

nlcd_path = "./Data/nlcd_2019_land_cover_l48_20210604/nlcd_2019_land_cover_l48_20210604.img"

## Read in raster with Terra package 
nlcd_ras = terra::rast(nlcd_path)
# methods(class = class(nlcd_ras))

# nrow(nlcd_ras) #104,424
# ncol(nlcd_ras) #161,190 
# nrow(nlcd_ras) * ncol(nlcd_ras) # 16.8 billion 
# # N=45000; memorytestmatrix <- matrix(nrow=N, ncol=N) # need to downsample raster

# Downsample to 100m resolution
# https://www.patrickbaylis.com/blog/2021-03-13-downsampling-magic/
nlcd_ras_100 = terra::aggregate(nlcd_ras, fact = 100/30, fun = "modal") 
nlcd_levels_100 = levels(nlcd_ras_100)

# project raster to desired crs
nlcd_proj_100 = terra::project(nlcd_ras_100, paste("EPSG:", mycrs))

# crop and mask raster 
nlcd_proj_100_crop = crop(nlcd_proj_100, vect(county_map_proj))
nlcd_proj_100_mask = mask(nlcd_proj_100_crop, vect(county_map_proj))

# extract values to each county 
nlcd_proj_100_extract = terra::extract(x = nlcd_proj_100_mask, y = vect(county_map_proj))
#save(nlcd_proj_100_extract, file = "./Data/nlcd_proj_100_extract.Rda")
#load(file = "./Data/nlcd_proj_100_extract.Rda")

try1 = nlcd_proj_100_extract %>% filter(ID < 301)

# Assign land classes based on NLCD legend 
# https://www.mrlc.gov/data/legends/national-land-cover-database-2011-nlcd2011-legend
nlcd_extract_join = try %>%
  mutate(LANDCLASS = case_when(
    `NLCD Land Cover Class` %in% c("Open Water", "Perennial Snow/Ice") ~ "Water",
    `NLCD Land Cover Class` %in% c("Developed, Open Space", "Developed, Low Intensity", "Developed, Medium Intensity", "Developed, High Intensity") ~ "Developed",
    `NLCD Land Cover Class` %in% c("Barren Land") ~ "Barren",
    `NLCD Land Cover Class` %in% c("Deciduous Forest", "Evergreen Forest", "Mixed Forest") ~ "Forest",
    `NLCD Land Cover Class` %in% c("Shrub/Scrub") ~ "Shrub",
    `NLCD Land Cover Class` %in% c("Herbaceous")~ "Herbaceous",
    `NLCD Land Cover Class` %in% c("Hay/Pasture", "Cultivated Crops") ~ "Cultivated",
    `NLCD Land Cover Class` %in% c("Woody Wetlands", "Emergent Herbaceous Wetlands") ~ "Wetlands",
    TRUE ~ "Other"
  ))



#################################################################################################

asd = crop(nlcd_proj_100, vect(census_map_proj))
asdf = mask(asd, vect(census_map_proj))
asdf_extract = terra::extract(x = asdf, y = vect(census_map_proj))


# https://www.mrlc.gov/data/legends/national-land-cover-database-2011-nlcd2011-legend
water = which(nlcd_levels_100[[1]] %in% c("Open Water", "Perennial Snow/Ice"))
developed = which(nlcd_levels_100[[1]] %in% c("Developed, Open Space", "Developed, Low Intensity", "Developed, Medium Intensity", "Developed, High Intensity"))
barren = which(nlcd_levels_100[[1]] %in% c("Barren Land"))
forest = which(nlcd_levels_100[[1]] %in% c("Deciduous Forest", "Evergreen Forest", "Mixed Forest"))
shrub = which(nlcd_levels_100[[1]] %in% c("Shrub/Scrub"))
herbaceous = which(nlcd_levels_100[[1]] %in% c("Herbaceous"))
cultivated = which(nlcd_levels_100[[1]] %in% c("Hay/Pasture", "Cultivated Crops"))
wetland = which(nlcd_levels_100[[1]] %in% c("Woody Wetlands", "Emergent Herbaceous Wetlands"))






state = "Texas"
county = c("Harris County")
census_map = get_acs(geography = "tract", state = state, county = county,
                     variables=c("B01003_001"), year = year, geometry = TRUE, 
                     cache_table = TRUE)
census_map = census_map %>% select(GEOID, NAME) 
census_map_proj = census_map %>% 
  st_transform(mycrs) # project to Conic Equal Area Albers, EPSG:5070 

asd = crop(nlcd_proj_100, vect(census_map_proj))
asdf = mask(asd, vect(census_map_proj))
asdf_extract = terra::extract(x = asdf, y = vect(census_map_proj))

# Assign land classes based on NLCD legend 
# https://www.mrlc.gov/data/legends/national-land-cover-database-2011-nlcd2011-legend
asdf_extract_join = asdf_extract %>%
  mutate(LANDCLASS = case_when(
    `NLCD Land Cover Class` %in% c("Open Water", "Perennial Snow/Ice") ~ "Water",
    `NLCD Land Cover Class` %in% c("Developed, Open Space", "Developed, Low Intensity", "Developed, Medium Intensity", "Developed, High Intensity") ~ "Developed",
    `NLCD Land Cover Class` %in% c("Barren Land") ~ "Barren",
    `NLCD Land Cover Class` %in% c("Deciduous Forest", "Evergreen Forest", "Mixed Forest") ~ "Forest",
    `NLCD Land Cover Class` %in% c("Shrub/Scrub") ~ "Shrub",
    `NLCD Land Cover Class` %in% c("Herbaceous")~ "Herbaceous",
    `NLCD Land Cover Class` %in% c("Hay/Pasture", "Cultivated Crops") ~ "Cultivated",
    `NLCD Land Cover Class` %in% c("Woody Wetlands", "Emergent Herbaceous Wetlands") ~ "Wetlands",
    TRUE ~ "Other"
    ))

# Group based on shapefile polygon (e.g., county) and LANDCLASS
asdf_summ = asdf_extract_join %>%
  group_by(ID, LANDCLASS) %>%
  count()

# expand LANDCLASS into wide format 
asdf_wide = asdf_summ %>%
  spread(key = LANDCLASS, value = n)

# column totals of land type 
asdf_wide2 = asdf_wide %>%
  replace(is.na(.), 0) %>%
  mutate(Total = rowSums(across(Barren:Wetlands)))

# percentages of land type
asdf_perc = asdf_wide2 %>%
  mutate(across(Barren:Total, ~ .x / Total))

# join to census shapefile
census_map_static = census_map_proj %>%
  bind_cols(asdf_perc)

gg2 = ggplot(census_map_static)+
  geom_sf(aes(fill = Developed), color = NA) + 
  scale_fill_viridis_c(option="plasma", na.value = "grey50") +
  #geom_sf(fill = NA, show.legend = F, color = "black", lwd = 0.005)+
  #coord_sf(datum = NA) + #removes gridlines 
  #guides(fill = "none") + #removes legend
  #theme_minimal()  #removes background
  theme_dark() +
  labs(title = "Harris County: Land Use", fill = "Developed\n(prop.)") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
  )
gg2

# mask raster based on land use type 
water_mask = developed_mask = barren_mask = forest_mask = shrub_mask = herbaceous_mask = cultivated_mask = wetland_mask = nlcd_proj_100

water_mask[!(water_mask %in% water)] = NA #set non-water pixels to NA 
developed_mask[!(developed_mask %in% developed)] = NA #set non-water pixels to NA 

# extract values for each county
nlcd_extract = terra::extract(x = nlcd_proj_100, y = vect(county_map_area))



# # https://cyberhelp.sesync.org/geospatial-packages-in-R-lesson/2017/07/18/
# library(raster)
# nlcd_ras = raster(nlcd_path)
# 
# df = nlcd_ras@data@attributes[[1]]
# land_class = df$NLCD.Land.Cover.Class
# levels(land_class)
# crs(nlcd_ras)
# 
# 
# pasture = mask(nlcd_ras, nlcd_ras == 81, maskvalue = FALSE)

# Star Proxy
# https://cyberhelp.sesync.org/geospatial-packages-in-R-lesson/index.html
nlcd_stars = read_stars(nlcd_path)
class(nlcd_stars)
methods(class = "stars_proxy")
nlcd_stars = droplevels(nlcd_stars)


TX = get_acs(geography = "county", state = "TX", county = "Harris County",
             variables=c("B01003_001"), year = year, geometry = TRUE, 
             cache_table = TRUE)
TX_proj = st_transform(TX, crs = st_crs(nlcd_stars))

census_map = get_acs(geography = "tract", state = "TX", county = "Harris County",
                     variables=c("B01003_001"), year = year, geometry = TRUE, 
                     cache_table = TRUE)
census_proj = st_transform(TX, crs = st_crs(nlcd_stars)) 


nlcd_stars_crop = st_crop(nlcd_stars, TX_proj) #crop 
nlcd_stars_mask = st_crop(nlcd_stars_crop, TX_proj, crop = FALSE) #mask 

asdf = st_as_stars(nlcd_stars_mask)
levels(asdf$nlcd_2019_land_cover_l48_20210604.img)
asdf = droplevels(asdf)

forest_types <- c('Evergreen Forest', 'Deciduous Forest', 'Mixed Forest')
forest_mask <- asdf
forest_mask[!(forest_mask %in% forest_types)] <- NA

#just sum entire table 
try = forest_mask %>% pull %>% table
ncells = dim(forest_mask)[1] * dim(forest_mask)[2]

mymode <- function(x) names(which.max(table(x)))
modal_lc <- aggregate(forest_mask, TX_proj, FUN = mymode)

try2 = aggregate(forest_mask, TX_proj, FUN = count)