#### LOAD PACKAGES and SET DIRECTORY 
options(java.parameters = "-Xmx7g")
library(bartMachine)
set_bart_machine_num_cores(4)
library(lme4)
library(corrplot)
library(viridis)
library(tidyverse)
library(tidymodels)
library(tidycensus) 
library(sf)
library(terra)
library(stars)
library(raster) #make sure ncdf4 package is installed 
library(lubridate)
library(doParallel)
library(xgboost)
library(vip)
library(pdp)

setwd("~/Documents/01_VECTOR.nosync/Sandia")
tidymodels_prefer()

##############################################################################################
## LOAD  VARIABLES TO USE ####################################################################
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

# Join in socio-economic variables that were included in the Five Factors 
factor_variables = c("B01003_001", #total population
                     "B25077_001", #median housing value (Factor 1)
                     "B01001D_001", #asian population (part of Factor 1)
                     "B17001_002", #poverty population (Factor 2) 
                     "B01001B_001", #black population (part of Factor 2)
                     "B19055_002", "B19055_001", #households receiving social security and total households (Factor 4) 
                     "B01001I_001" #hispanic/latino population (Factor 5)
                     )
mycrs = 5070 #chose projected coordinate system: EPSG 5070 NAD83 Conus Albers
year=2019 # year for county boundaries 
options(tigris_use_cache = TRUE) #cache shapefiles for future sessions
state_list = c("AL", "AR", "AZ", "CA", "CO", "CT", "DE", "FL", "GA", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")
factor_acs_data = get_acs(geography = "county",state = state_list, variables=factor_variables, year = year, geometry = FALSE)
factor_acs_data_w = factor_acs_data %>%
  dplyr::select(-moe, -NAME) %>% #remove errors of estimates because we're spreading by variable and don't want duplicates, and don't need name b/c joining by GEOID
  spread(key=variable, value = estimate)
factor_dat = factor_acs_data_w %>%
  mutate(MDHVAL = B25077_001) %>% # median housing value 
  mutate(QASIAN = B01001D_001 / B01003_001) %>% # %asian 
  mutate(QPOVERTY = B17001_002 / B01003_001) %>% # %below poverty line
  mutate(QBLACK = B01001B_001 / B01003_001) %>% # %black 
  mutate(QSSBEN =  B19055_002 / B19055_001) %>% #%households receiving social security 
  mutate(QSPANISH = B01001I_001 / B01003_001) %>% # %hispanic or latino 
  dplyr::select(c(GEOID, MDHVAL:QSPANISH))

# Select only data we'll have available for prediction at the census tract level  
county_outages_GROUP_CLEAN = county_outages_GROUP %>%
  dplyr::left_join(factor_dat, by = c("GEOID")) %>% #join in factor variables 
  dplyr::select(-c(PS_mean:V_max, T_mean:TQV_mean, AMBULANCES:DEATHS, EMPBLDG:HOSPBEDS, INTERNET:NURSHOMES, PHYSICIANS, PROFORGS:PROXCAP, RADIO:WATEFF, QEXTRCT, QNRRES, POPSTAB)) %>%
  dplyr::select(-c(Factor1, Factor2, Factor3, Factor4, Factor5)) %>% #instead use variables within factors
  dplyr::select(-c(soil10_lag, soil100_lag)) %>% #don't contribute enough to warrant the hassle of getting them
  dplyr::select(-c(RZ_mode, QGROUPHSE, QSERVIND, QAGEWORK, QFEMLBR)) #take out variables where effect has high uncertainty, likely from lingering multicollinearity

# Final data frame 
df_data = data.frame(county_outages_GROUP_CLEAN) %>%
  dplyr::filter(out_hrs > quantile(county_outages_GROUP$out_hrs, .9)) %>% #filter to big events 
  dplyr::select(-outage_number, -GEOID, -out_percust) %>% #get rid of variables not used 
  mutate(out_hrs = log(out_hrs), out_maxcust = log(out_maxcust)) %>% #take log of DVs
  rename(ln_hrs = out_hrs, ln_cust = out_maxcust) #rename to emphasize ln 
rm(list=setdiff(ls(), "df_data")) 
gc() 

##############################################################################################################
#### MACHINE LEARNING ########################################################################################
##############################################################################################################
num_cores = detectCores() - 1
unregister_dopar = function() {env <- foreach:::.foreachGlobals; rm(list=ls(name=env), pos=env)}

## Split into training vs testing
set.seed(23)
df_split = initial_split(df_data, prop = 0.80, strata = "ln_hrs")
df_train = training(df_split)
df_test = testing(df_split)  
df_cv = vfold_cv(df_train, v = 10, repeats = 1)

## Pre-processing (recipe)
cust_recipe = recipe(ln_cust ~ . , data = df_data) %>%
  step_rm(ln_hrs) %>% 
  step_impute_knn(all_predictors()) %>% #knn impute missing predictors (if any)
  # step_normalize(all_predictors()) %>% #z-score standardize all predictors (important for PLS or NN)
  step_zv(all_predictors()) %>% #removes predictors of single value 
  step_corr(all_predictors())  #removes highly correlated 
prep(cust_recipe) #shows changes 
cust_juice = prep(cust_recipe) %>% juice() #view prepared dataset 

## BART 
X = cust_juice %>% slice(df_split$in_id) %>% dplyr::select(-ln_cust) %>% as.data.frame()
y = cust_juice %>% dplyr::select(ln_cust) %>% slice(df_split$in_id) %>% pull()
bart_train = bartMachineCV(X, y, k_folds = 10)  #bartMachine CV win: k: 2 nu, q: 3, 0.9 m: 200 

X_test = cust_juice %>% slice(-df_split$in_id) %>% dplyr::select(-ln_cust) %>% as.data.frame()
y_test = cust_juice %>% slice(-df_split$in_id) %>% dplyr::select(ln_cust) %>% pull()
bart_test = predict(bart_train, X_test)
bart_CI = round(calc_credible_intervals(bart_train, X_test, ci_conf = 0.95), 2)
rsq_bart = 1 - sum((y_test - bart_test)^2) / sum((y_test - mean(y_test))^2) 
cverror_bart = paste(data.frame(bart_train$cv_stats) %>% dplyr::slice(1) %>% pull(oos_error) %>% round(3) %>% format(nsmall = 3))

vi = investigate_var_importance(bart_train, num_replicates_for_avg = 20)
vs = var_selection_by_permute(bart_train, bottom_margin = 10, num_reps_for_avg = 20, num_permute_samples = 20, plot = F)
vs$important_vars_local_names
vs$important_vars_global_se_names

## XGB 
show_model_info("boost_tree")
xgb_model = boost_tree(mode = "regression",
                       trees = tune(), min_n = tune(), tree_depth = tune(), learn_rate = tune(), loss_reduction = tune(), mtry = tune()) %>%
  set_engine(engine = "xgboost") %>%
  translate()
xgb_work = workflow() %>%
  add_recipe(cust_recipe) %>%
  add_model(xgb_model)
xgb_grid = dials::grid_max_entropy(parameters(trees(), min_n(), tree_depth(), learn_rate(), loss_reduction(), finalize(mtry(), cust_juice)), size = 100)
cl = makeCluster(num_cores, type = "FORK")
registerDoParallel(cl, cores = num_cores)
xgb_tune = xgb_work %>%
  tune_grid(resamples = df_cv, 
            grid = xgb_grid,
            metrics = metric_set(yardstick::rmse, yardstick::rsq),
            control = tune::control_grid(verbose = T, allow_par = T, parallel_over = "resamples")) #parallel processing turns off verbose
stopCluster(cl) 
unregister_dopar()
show_best(xgb_tune, metric = "rmse")
xgb_tune_results = xgb_tune %>% collect_metrics()
xgb_best = xgb_tune %>% select_best(metric = "rmse")
xgb_fit = xgb_work %>%
  finalize_workflow(xgb_best) %>%
  last_fit(df_split)
xgb_test = xgb_fit %>% collect_metrics() #metrics evaluated on test sample (b/c last_fit() function) 
xgb_predictions = xgb_fit %>% collect_predictions() #predictions for test sample (b/c last_fit() function)
rsq_xgb = 1 - sum((y_test - xgb_predictions$.pred)^2) / sum((y_test - mean(y_test))^2) 
cverror_xgb = paste(show_best(xgb_tune, metric = "rmse") %>% dplyr::slice(1) %>% pull(mean) %>% round(3) %>% format(nsmall = 3))

xgb_train = xgb_work %>%
  finalize_workflow(xgb_best) %>%
  fit(df_train) %>%
  extract_fit_parsnip()
importance = xgb.importance(model = xgb_train$fit)
head(importance, n = 15)
vip(xgb_train$fit, n = 20)


## Lasso/Ridge/ElasticNet 
show_model_info("linear_reg")
lre_model = linear_reg(penalty = tune(), mixture = tune()) %>% #lambda (penalty) and alpha/mixture (1 lasso, 0 ridge)
  set_engine("glmnet") %>%
  translate()
lre_work = workflow() %>% 
  add_recipe(cust_recipe) %>%
  add_model(lre_model)
set.seed(32); lre_grid = dials::grid_max_entropy(parameters(penalty(), mixture()), size = 40)
cl = makeCluster(num_cores, type = "FORK")
registerDoParallel(cl, cores = num_cores)
lre_tune = lre_work %>%
  tune_grid(resamples = df_cv,
            grid = lre_grid,
            metrics = metric_set(yardstick::rmse, yardstick::rsq),
            control = tune::control_grid(verbose = T, allow_par = T, parallel_over = "resamples")
  ) #parallel processing turns off verbose
stopCluster(cl) 
unregister_dopar()
show_best(lre_tune, metric = "rmse")
lre_tune_results = lre_tune %>% collect_metrics()
lre_best = lre_tune %>% select_best(metric = "rmse")
lre_fit = lre_work %>%
  finalize_workflow(lre_best) %>%
  last_fit(df_split)
lre_test = lre_fit %>% collect_metrics() #metrics evaluated on test sample (b/c last_fit() function) 
lre_predictions = lre_fit %>% collect_predictions() #predictions for test sample (b/c last_fit() function)



##############################################################################################################
#### PLOTTING ################################################################################################
##############################################################################################################
gg_test = dplyr::tibble(actual = y_test,
                        #bart = bart_test,
                        eNet = as.vector(lre_predictions$.pred), 
                        xgb = as.vector(xgb_predictions$.pred)
                        )
gg_test = arrange(gg_test, actual)
gg_test$index = seq.int(nrow(gg_test))
gg_test_l = gg_test %>% pivot_longer(!index, names_to = "Model", values_to = "ypred")
plot_filtering_estimates2 <- function() {
  p = ggplot() + 
    theme_classic() + 
    geom_hline(yintercept = mean(gg_test$actual, na.rm = T), linetype="dashed", color = "gray50") +
    geom_line(data = gg_test_l, aes(x = index, y = ypred, color = Model), alpha = 0.8) +
    ylab("Max Cust. Outages (log)") + 
    #ylab("Max Cust. Outages (log)") + 
    scale_y_continuous(labels = function(x) paste0(x)) +
    xlab("Outage Index (event x county)") +
    ggtitle("County-level Test") + 
    theme(legend.spacing.y = unit(-0.25, "cm"),
          legend.direction = "vertical",
          legend.box = "vertical",
          legend.position = c(.215, .8),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)
    ) 
  print(p)
}
plot_filtering_estimates2()

##### EXPORT MODEL
# BART does a tiny bit better and WIND_max is very intuitive/we could just impute missing entries as sub-cyclone winds b/c they're effect is level
bart_final = bartMachine(as.data.frame(cust_juice %>% dplyr::select(-ln_cust)), as.vector(cust_juice %>% pull(ln_cust)), k = 2, nu = 3, q = 0.99, num_trees = 200, serialize = T)
pdp_wind = pd_plot(bart_final, j = "WIND_max")
save(bart_final, file = "bart_final.Rda")

# But XGB handles missing data much better and pretty much same fit
xgb_final = xgb_work %>%
  finalize_workflow(xgb_best) %>%
  fit(df_data) 
xgb_final_model = extract_fit_parsnip(xgb_final)$fit
save(xgb_final_model, file = "xgb_final_model.Rda")
xgb.plot.tree(model = xgb_final_model, trees = 1:5) 
p_wind =  pdp::partial(xgb_final_model, pred.var = "WIND_max", ice = T, center = F,
                       plot = T, rug= T, alpha = 0.1, plot.engine = "ggplot2",
                       train = as.data.frame(cust_juice %>% dplyr::select(-ln_cust)), type = "regression") +
  geom_vline(xintercept = 8.9, color = "blue") + 
  ggtitle("PDP and ICE - Max Wind Speed (m/s)") + 
  annotate("text", x =12.5, y = 6, label="Tropical Storm Wind Forecast", color="blue", size = 3) 

p_oth =  pdp::partial(xgb_final_model, pred.var = "spi03_lag", ice = T, center = F,
                       plot = T, rug= T, alpha = 0.1, plot.engine = "ggplot2",
                       train = as.data.frame(cust_juice %>% dplyr::select(-ln_cust)), type = "regression") +
  ggtitle("PDP and ICE Plot") 

##### PLOT FINAL MODEL (ON TRAINING SET) WITH AND WITHOUT DYNAMIC FORECAST
xgb_final_train = xgb_work %>%
  finalize_workflow(xgb_best) %>%
  fit(df_train) 
xgb_final_train_model = extract_fit_parsnip(xgb_final_train)$fit
X_static = cust_juice %>% 
  slice(-df_split$in_id) %>% 
  dplyr::select(-ln_cust) %>% 
  mutate(WIND_max = NA) %>% 
  as.matrix()
X_dynamic = cust_juice %>% 
  slice(-df_split$in_id) %>% 
  dplyr::select(-ln_cust) %>% 
  #mutate(WIND_max = replace(WIND_max, which(WIND_max < 10), NA)) %>% #strip out non-hurricane winds 
  as.matrix()
xgb_pred_static = predict(xgb_final_train_model, X_static)
xgb_pred_dynamic = predict(xgb_final_train_model, X_dynamic)
rsq_xgb_static = round(1 - sum((y_test - xgb_pred_static)^2) / sum((y_test - mean(y_test))^2), 3)
rsq_xgb_dynamic = round(1 - sum((y_test - xgb_pred_dynamic)^2) / sum((y_test - mean(y_test))^2), 3)

color_vec = c("black", "#eb8055ff", "#253582ff") 
gg_civic = dplyr::tibble(actual = y_test,
                        xgb_dynamic = as.vector(xgb_pred_dynamic),
                        xgb_static = as.vector(xgb_pred_static)
)
gg_civic = arrange(gg_civic, actual)
gg_civic$index = seq.int(nrow(gg_civic))
gg_civic_l = gg_civic %>% pivot_longer(!index, names_to = "Model", values_to = "ypred")
plot_filtering_estimates2 <- function() {
  p = ggplot() + 
    theme_classic() + 
    geom_hline(yintercept = mean(gg_civic$actual, na.rm = T), linetype="dashed", color = "gray50") +
    geom_line(data = gg_civic_l, aes(x = index, y = ypred, color = Model), alpha = 0.8) +
    scale_color_manual(
      values = color_vec, 
      labels = c("Actual",
                 bquote("XGB Dynamic (" * R^2 ~ "=" ~ .(rsq_xgb_dynamic) * ")"), 
                 bquote("XGB Baseline (" * R^2 ~ "=" ~ .(rsq_xgb_static) * ")")
      ),
      name = element_blank()) +
    ylab("Max Cust. Outages (log)") + 
    #ylab("Max Cust. Outages (log)") + 
    scale_y_continuous(labels = function(x) paste0(x)) +
    xlab("Outage Index (event x county)") +
    ggtitle("County-level Test Sample") + 
    theme(legend.spacing.y = unit(-0.25, "cm"),
          legend.direction = "vertical",
          legend.box = "vertical",
          legend.position = c(.215, .8),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)
    ) 
  print(p)
}
plot_filtering_estimates2()

