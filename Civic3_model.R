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
### LOAD VARIABLES ###########################################################################
##############################################################################################
# Get census tract map 
mycrs = 5070 #chose projected coordinate system: EPSG 5070 NAD83 Conus Albers
year=2019 # year for county boundaries 
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
census_map_area = census_map %>%  
  mutate(AREA = as.vector(st_area(census_map))) %>% #sq-meters; as.vector removes units suffix 
  mutate(DENSITY = POPULATION / AREA * 1000^2) #population per sq-km 

load(file = "./Data/census_map_nlcd.Rda")  #NLCD
census_map_nlcd = census_map_nlcd %>% st_set_geometry(NULL)
load(file = "./Data/census_map_dem.Rda") #DEM
census_map_dem = census_map_dem %>% st_set_geometry(NULL)
load(file = "./Data/census_map_rz.Rda") #Root Zone
census_map_rz = census_map_rz %>% st_set_geometry(NULL)
load(file = "./Data/census_map_social.Rda") #Socio-economic
census_map_social = census_map_social %>% st_set_geometry(NULL)
load(file = "./Data/census_map_spi.Rda") #SPI 
census_map_spi = census_map_spi %>% st_set_geometry(NULL)
load(file = "./Data/census_map_soil_day.Rda") #Soil Moisture

census_map_CLEAN = census_map_area %>%
  dplyr::select(GEOID, DENSITY) %>%
  inner_join(dplyr::select(census_map_nlcd, c(GEOID, Developed:Wetlands)), by = "GEOID") %>%
  inner_join(dplyr::select(census_map_dem, c(GEOID, DEM_mean:DEM_max)), by = "GEOID") %>%
  inner_join(dplyr::select(census_map_rz, c(GEOID, RZ_mean:RZ_mode)), by = "GEOID") %>%
  inner_join(dplyr::select(census_map_social, c(GEOID, 1, 6:ncol(census_map_social))), by = "GEOID") %>%
  inner_join(dplyr::select(census_map_spi, c(GEOID, spi03_mean:spi24_mean)), by = "GEOID") %>%
  inner_join(dplyr::select(census_map_soil_day, c(GEOID, soil10_3dLAG:soil100_3dLAG)), by = "GEOID") %>%
  rename(spi03_lag = spi03_mean, spi12_lag = spi12_mean, spi24_lag = spi24_mean, 
         soil10_lag = soil10_3dLAG, soil40_lag = soil40_3dLAG, soil100_lag = soil100_3dLAG,
         Density = DENSITY)

#impute missing values based on neighboring polygons
which(is.na(census_map_CLEAN), arr.ind = T)
index = st_touches(census_map_CLEAN, census_map_CLEAN) %>%
impute = function(x)

##############################################################################################
### BART MODEL ###############################################################################
##############################################################################################
load(file = "bart_civic.Rda")  #load trained model 
X_2 = census_map_CLEAN %>%
  #st_set_geometry(NULL) %>%
  dplyr::select(c(GEOID, "spi24_lag","Density", "Developed", "QMOHO", "soil100_lag" , "QFEMALE", "QNATIVE", "spi12_lag", "Wetlands", "RZ_mean")) %>%
  drop_na()
X = X_2 %>%
  dplyr::select(-GEOID) %>%
  st_set_geometry(NULL)
model_name = paste("Predicting Outages - No Weather Data - Harris County, TX")

predictions = predict(bart, X)

X_2$predictions = predictions


gg2 = ggplot(X_2)+
  geom_sf(aes(fill = predictions), color = NA) + 
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


