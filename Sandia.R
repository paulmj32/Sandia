library(tidyverse)
library(tidycensus) 
library(sf)
library(terra)
library(stars)

# library(lme4)
# options(java.parameters = "-Xmx5g")
# library(bartMachine)
# set_bart_machine_num_cores(4)

setwd("~/Documents/01_VECTOR.nosync/Sandia")
mycrs = 5070

##### CONTIGUOUS US COUNTY MAP ############################################
year=2019
options(tigris_use_cache = TRUE) #to cache shapefiles for future sessions
state_list = c("AL", "AR", "AZ", "CA", "CO", "CT", "DE", "FL", "GA", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")
county_map = get_acs(geography = "county", state = state_list,
                     variables=c("B01003_001"), year = year, geometry = TRUE, 
                     cache_table = TRUE)
county_map = county_map %>%
  mutate(POPULATION = estimate) %>%
  dplyr::select(GEOID, NAME, POPULATION) 

# # county map info 
# st_crs(county_map)
# st_crs(county_map)$IsGeographic  
# st_crs(county_map)$units_gdal
# st_crs(county_map)$proj4string

## Project county map and calculate area and population density of each county 
county_map_proj = county_map %>% 
  st_transform(mycrs) # project to Conic Equal Area Albers, EPSG:5070 

county_map_area = county_map_proj %>%  
  mutate(AREA = as.vector(st_area(county_map_proj))) %>% ## calculate area of each county (sq-meters); as.vector removes units suffix 
  mutate(DENSITY = POPULATION / AREA * 1000^2) #population per sq-km 

gg = ggplot(county_map_area)+
  geom_sf(aes(fill = DENSITY), color = NA) + 
  scale_fill_viridis_c(option="plasma", na.value = "grey50") +
  #geom_sf(fill = NA, show.legend = F, color = "black", lwd = 0.005)+
  #coord_sf(datum = NA) + #removes gridlines 
  #guides(fill = "none") + #removes legend
  theme_minimal()  #removes background

pdf("figure_county.pdf", width = 7.48, height = 4.5)
gg
dev.off()

##### STATIC VARIABLES ##############################################################
## NLCD

## DEM

## Root Zone 



##### OUTAGE DATA #############################################
# Read in data
outages_csv = read.csv("./Data/SE_states_outage_merra_2018.csv", header = T)
save(outages_csv, file = "./Data/outages.Rda")
load(file = "Data/outages.Rda")

# Format FIPS to character and include starting zero
outages_csv$fips_code = as.character(outages_csv$fips_code)
fips = outages_csv$fips_code
fips0 = str_pad(fips, 5, pad = "0")
outages_csv$fips_code = fips0

view(head(outages_csv, n = 100))



## Exploratory analysis

#autocorrelation of one county through time (make sure to find example where outagedata) 





#boxplot of hr_mean_frac_out to see where outliers may be
summary(outages_csv$hr_mean_frac_out) #highly skewed right 
boxplot(outages_csv$hr_mean_frac_out)

#filter only outage events
outages_only = outages_csv %>%
  filter(hr_mean_frac_out >= 0.01) #could potentially use status instead? 
summary(outages_only$hr_mean_frac_out)
boxplot(outages_only$hr_mean_frac_out)

wind_10m = outages_csv %>%
  filter(WIND10M >= 4.69) #90th percentile

#order data with respect to mean frac_out
outages_only1 = outages_only %>% arrange(hr_mean_frac_out)
plot(outages_only1$hr_mean_frac_out)

outages_bad = outages_csv %>%
  filter(hr_mean_frac_out > .8)


## Try some statistical models
# BART
df_bart = outages_csv %>% 
  filter(hr_mean_frac_out >= 0.3) %>% #filter out non-outages
  dplyr::select(-c(state_fips, state, county, fips_code, year, date_hour, county_population, hr_mean_customers_out, outage_status, outage_number))
y = df_bart$hr_mean_frac_out
X = df_bart; X$hr_mean_frac_out = NULL

set.seed(32) #for some reason BART doesn't like set.seed() / the command doesn't afffect stochastic elements (just repeat)
bart_cv = bartMachineCV(X, y) #bartMachine CV win: k: 2 nu, q: 5, 0.99 m: 200 
print(bart_cv)

#predictions
predictions = predict(bart_cv, X)
rmse = sqrt(mean((predictions - df_bart$hr_mean_frac_out)^2)) 
rsq = 1 - sum((df_bart$hr_mean_frac_out - predictions)^2) / sum((df_bart$hr_mean_frac_out - mean(df_bart$hr_mean_frac_out))^2)

#check BART assumptions
check_bart_error_assumptions(bart_cv)
plot_convergence_diagnostics(bart_cv)

CI = round(calc_credible_intervals(bart_cv, X, ci_conf = 0.95), 2)
gg = dplyr::tibble(x_mean = predictions,
                   lower = CI[,1],
                   upper = CI[,2],
                   actual = df_bart$hr_mean_frac_out
)
gg = arrange(gg, actual)
gg$index = seq.int(nrow(gg))

plot_filtering_estimates <- function(df) {
  p <- ggplot(data = gg, aes(x = index)) +
    theme_classic() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4, fill = "red") + #quantiles 
    geom_point(aes(y = actual), colour = "gray32", #actual observation points
               size = 0.9, shape = 16, alpha = 0.9) +
    geom_line(aes(y = x_mean), colour = "red", size = 0.4) + #mean estimate
    ylab("hr_mean_frac_out") + 
    scale_y_continuous(labels = function(x) paste0(x)) +
    xlab("Index") +
    ggtitle("Bayesian Additive Regression Tree") +
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
}
plot_filtering_estimates(gg)

summary = investigate_var_importance(bart_cv, num_replicates_for_avg = 5)
pd_plot(bart_cv, j = "PS")
pd_plot(bart_cv, j = "WIND2M")
pd_plot(bart_cv, j = "SLP")
pd_plot(bart_cv, j = "V2M")
pd_plot(bart_cv, j = "TQV")

library(corrplot)
library(viridis)
mycor = cor(X)
col = viridis(100, direction = -1, option = "C")
corrplot(mycor, method = "circle", tl.col="black", tl.srt=45, tl.cex = 0.7, col = col, cl.cex = 0.7,
         order = "hclust", type = "upper", diag = T, mar = c(1, 1, 1, 1))



## Get US county sf file
year=2018
options(tigris_use_cache = TRUE) #to cache shapefiles for future sessions
county_map = get_acs(geography = "county", 
                     variables=c("B01003_001"), year = year, geometry = TRUE, 
                     cache_table = TRUE)

# #join data with sf file
# county_outages = county_map %>%
#   left_join(outages_csv, by = c("GEOID" = "fips_code"))





#### DATA ###################################################################################


## Root zone data
TX_gdb =  "/Users/paulmj/Downloads/gSSURGO_TX/gSSURGO_TX.gdb"
# TX = sf::st_read(dsn = TX_gbd_ext, layer = "MUPOLYGON")
# TX_group = TX %>%
#   group_by(MUKEY) %>%
#   summarise(n = n())
# save(TX_group, file = "TX_group.Rda")
load("Data/TX_group.Rda")
TX_Valu1 = sf::st_read(dsn = TX_gdb, layer = "Valu1")
TX_group_val1 = TX_group %>% left_join(TX_Valu1, by = c("MUKEY" = "mukey"))

TX_rast_1k = st_rasterize(TX_group_val1["rootznemc"], dx = 1000, dy = 1000) #1km x 1km resolution 
TX_rast_100m = st_rasterize(TX_group_val1["rootznemc"], dx = 100, dy = 100) #100m x 100km resolution 

rr = ggplot() + 
  geom_stars(data = TX_rast_1k, aes(x = x, y = y, fill = rootznemc)) + 
  scale_fill_viridis_c(direction = -1, na.value = "gray") +
  theme_minimal()
rr

st_crs(TX_group_val1)
TX_map = st_transform(TX_group_val1, crs_map)

### Geometry base - US Counties 
year=2019
options(tigris_use_cache = TRUE) #to cache shapefiles for future sessions
state = "Texas"
county = c("Harris County")

county_map = get_acs(geography = "county", state = state,
                     variables=c("B01003_001"), year = year, geometry = TRUE, 
                     cache_table = TRUE)
county_map = county_map %>% select(GEOID, NAME) 

census_map = get_acs(geography = "tract", state = state, county = county,
                     variables=c("B01003_001"), year = year, geometry = TRUE, 
                     cache_table = TRUE)
census_map = census_map %>% select(GEOID, NAME) 


## PROJECTION (USA Equal Area Conic: EPSG 5070) 
Houston_census = st_transform(census_map, crs = st_crs(TX_rast_100m))
TX_rast_100m_crop = st_crop(TX_rast_100m, Houston_census) #crop to Harris county 
TX_rast_100m_crop_mask = st_crop(TX_rast_100m_crop, Houston_census, crop = FALSE) #mask to Harris county

hh = ggplot() + 
  geom_stars(data = TX_rast_100m_crop_mask, aes(x = x, y = y, fill = rootznemc)) + 
  scale_fill_viridis_c(direction = -1, na.value = "transparent") +
  geom_sf(data = Houston_census, colour = alpha("white", 0.5), fill = NA, size = 0.45) +
  theme_dark() +
  labs(title = "Harris County: Root Zone Depth", fill = "Root Zone\n(cm)") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
        )
hh

## TREE DATA
tree_db = "/Users/paulmj/Downloads/L48_Totals/L48_Totals.gdb"
st_layers(tree_db)


crs_map = st_crs(census_map)
st_crs(census_map)$IsGeographic  
st_crs(census_map)$units_gdal
st_crs(census_map)$proj4string

## SOIL MOISTURE
soil_db = "/Users/paulmj/Downloads/NDLAS_NOAH_2018_monthly/NLDAS_NOAH0125_M.A201801.002.grb"
asd = read_stars(soil_db, along = "band")
st_dimensions(asd)

library(rNOMADS)
GribInfo(soil_db, file.type = "grib1")
soil0_10 = asd[ , , , 26]
soil10_40 = asd[ , , , 27]
soil40_100 = asd[ , , , 28]


