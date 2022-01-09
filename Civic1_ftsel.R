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
vs = var_selection_by_permute(bart, bottom_margin = 10, num_reps_for_avg = 20, num_permute_samples = 20, plot = F)
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

#feature selection2 - before we get drop off in importance 
X = df_bart %>%
  dplyr::select(-outage_number, -GEOID, -out_hrs, -out_maxcust, -out_percust) %>%
  dplyr::select(c("spi24_lag","Density", "Developed", "QMOHO", "soil100_lag" , "QFEMALE", "QNATIVE", "spi12_lag", "Wetlands", "RZ_mean")) 
model_name = paste("Predicting Max Outages - No Weather Data")
bart = bartMachine(X, y, serialize = T) 
save(bart, file = "bart_civic.Rda")

