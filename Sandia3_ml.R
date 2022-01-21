#### LOAD PACKAGES and SET DIRECTORY 
options(java.parameters = "-Xmx7g")
library(bartMachine)
set_bart_machine_num_cores(4)
library(tidyverse)
library(tidymodels)
library(tidycensus) 
library(sf)
library(lme4)
library(corrplot)
library(viridis)
#library(sqldf)

tidymodels_prefer()
setwd("~/Documents/01_VECTOR.nosync/Sandia")

#### LOAD VARIABLES ##############################################################################
load("./Data/SandiaVariables.Rda") #from Sandia1_vars.R
load(file = "Data/outages.Rda")
outages_csv$fips_code = as.character(outages_csv$fips_code)
fips = outages_csv$fips_code
fips0 = str_pad(fips, 5, pad = "0")
outages_csv$fips_code = fips0
county_map_outages = outages_csv %>%
  mutate(rowID = row_number()) %>% # add row number to index events 
  mutate(date_day = str_extract(date_hour, "^.{10}")) %>% #returns first 10 characters of date_time (yyyy-mm-dd)
  mutate(date_month = str_extract(date_hour, "^.{7}")) # returns yyyy-mm
county_map_outages_FILTER = county_map_outages %>%
  filter(outage_status %in% c("pre", "start", "during", "end")) #filter events flagged as outages
county_map_outages_JOIN = county_map_area %>%
  inner_join(county_map_outages_FILTER, by = c("GEOID" = "fips_code")) %>% # join outage/hourly weather data 
  left_join(county_map_static_CLEAN, by = c("GEOID")) %>% #join static environmental variables
  left_join(county_scores_CLEAN, by = c("GEOID")) %>% #join socio-economic variables  
  left_join(county_map_soil_CLEAN, by = c("GEOID", "date_hour")) %>% #join soil moisture by county and hour time-stamp 
  left_join(county_map_spi_CLEAN, by = c("GEOID", "date_month")) #join SPI by GEOID and month 
outages_group = county_map_outages_JOIN %>%
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
            soil10_mean = mean(soil0_10), soil100_mean = mean(soil40_100),
            spi03_mean = mean(spi_03), spi12_mean = mean(spi_12),spi24_mean = mean(spi_24)
  )
county_outages_GROUP = outages_group %>%
  left_join(county_map_static_CLEAN, by = c("GEOID")) %>% #join static environmental variables
  left_join(county_scores_CLEAN, by = c("GEOID")) #join socio-economic variables 

#### MACHINE LEARNING ##########################################################################
df_data = data.frame(county_outages_GROUP) %>%
  dplyr::filter(out_hrs > quantile(county_outages_GROUP$out_hrs, .9)) %>% #filter to big events 
  dplyr::select(-outage_number, -GEOID, -out_percust) %>% #get rid of variables not used 
  mutate(out_hrs = log(out_hrs), out_maxcust = log(out_maxcust)) %>% #take log of DVs
  rename(ln_hrs = out_hrs, ln_cust = out_maxcust) #rename to emphasize ln 
rm(list=setdiff(ls(), "df_data")) 
gc() 

## Split into training vs testing
set.seed(32)
df_split = initial_split(df_data, prop = 0.80, strata = "ln_hrs")
df_train = training(df_split)
df_test = testing(df_split)  
df_cv = vfold_cv(df_train, v = 10, repeats = 1)

## Recipe for models
hours_recipe = recipe(ln_hrs ~ . , data = df_data) %>%
  step_rm(ln_cust) %>% 
  step_impute_knn(all_predictors()) %>% #knn impute missing predictors (if any)
  # step_normalize(all_predictors()) %>% #z-score standardize all predictors (important for PLS or NN)
  step_zv(all_predictors()) %>% #removes predictors of single value 
  step_corr(all_predictors())  #removes highly correlated 
prep(hours_recipe) #shows changes 
#hours_juice = prep(hours_recipe) %>% juice() #view prepared dataset 

## Specify models
model_name = paste("Large Events - Static, Socio-economic, and Dynamic Variables")

# BART (not part of tidymodels yet) 
df_bart = prep(hours_recipe) %>% juice()
df_bart_train = df_bart %>% slice(df_split$in_id) %>% dplyr::select(-ln_hrs)
X = data.frame(df_bart_train)
y = df_bart %>% dplyr::select(ln_hrs) %>% slice(df_split$in_id) %>% pull()
bart_fit = bartMachineCV(X, y, k_folds = 10) #bartMachine CV win: k: 2, nu: 3, q: 0.99, num_trees: 50 

df_bart_test = df_bart %>% slice(-df_split$in_id) %>% dplyr::select(-ln_hrs)
X_test = data.frame(df_bart_test)
y_test = df_bart %>% slice(-df_split$in_id) %>% dplyr::select(ln_hrs) %>% pull()
bart_pred = predict(bart_fit, X_test)
bart_rmse = sqrt(mean((bart_pred - y_test)^2)) 
bart_rsq = 1 - sum((y_test - bart_pred)^2) / sum((y_test - mean(y_test))^2) 
bart_CI = round(calc_credible_intervals(bart_fit, X_test, ci_conf = 0.95), 2)

# Random Forest
#https://www.rebeccabarter.com/blog/2020-03-25_machine_learning/
show_model_info("rand_forest")
rf_model = rand_forest(mtry = tune(), trees = tune(), min_n = tune()) %>%
  set_engine("ranger", importance = "permutation") %>%
  set_mode("regression") 
rf_work = workflow() %>%
  add_recipe(hours_recipe) %>%
  add_model(rf_model)
rf_grid = expand.grid(mtry = c(1, 3, 5), trees = c(100, 500, 1000), min_n = c(3, 5, 10))
rf_tune = rf_work %>%
  tune_grid(resamples = df_cv,
            grid = rf_grid,
            metrics = metric_set(yardstick::rmse, yardstick::rsq))
show_best(rf_tune, metric = "rmse")
rf_tune_results = rf_tune %>% collect_metrics()
rf_best = rf_tune %>% select_best(metric = "rmse")
rf_fit = rf_work %>%
  finalize_workflow(rf_best) %>%
  last_fit(df_split)
rf_test = rf_fit %>% collect_metrics() #metrics evaluated on test sample (b/c last_fit() function) 
rf_predictions = rf_fit %>% collect_predictions() #predictions for test sample (b/c last_fit() function)

# Lasso/Ridge/ElasticNet 
#https://dnield.com/posts/tidymodels-intro/ 
show_model_info("linear_reg")
lre_model = linear_reg(penalty = tune(), mixture = tune()) %>%
  set_engine("glmnet") 
lre_work = workflow() %>% 
  add_recipe(hours_recipe) %>%
  add_model(lre_model)
lre_grid = grid_regular(parameters(penalty(), mixture()), levels = c(5, 5))
lre_tune = lre_work %>%
  tune_grid(resamples = df_cv,
            grid = lre_grid,
            metrics = metric_set(yardstick::rmse, yardstick::rsq))
show_best(lre_tune, metric = "rmse")
lre_tune_results = lre_tune %>% collect_metrics()
lre_best = lre_tune %>% select_best(metric = "rmse")
lre_fit = lre_work %>%
  finalize_workflow(lre_best) %>%
  last_fit(df_split)
lre_test = lre_fit %>% collect_metrics() #metrics evaluated on test sample (b/c last_fit() function) 
lre_predictions = lre_fit %>% collect_predictions() #predictions for test sample (b/c last_fit() function)

# GBM 
#https://www.r-bloggers.com/2020/05/using-xgboost-with-tidymodels/
show_model_info("boost_tree")
gb_model = boost_tree(mode = "regression", trees = 1000, 
                      min_n = tune(), tree_depth = tune(), learn_rate = tune(), loss_reduction = tune()) %>%
  set_engine(engine = "xgboost")
gb_work = workflow() %>%
  add_recipe(hours_recipe) %>%
  add_model(gb_model)
gb_grid = dials::grid_max_entropy(parameters(min_n(), tree_depth(), learn_rate(), loss_reduction()), size = 100)
gb_tune = gb_work %>%
  tune_grid(resamples = df_cv,
            grid = gb_grid,
            metrics = metric_set(yardstick::rmse, yardstick::rsq),
            control = tune::control_grid(verbose = T)) 
show_best(gb_tune, metric = "rmse")
gb_tune_results = gb_tune %>% collect_metrics()
gb_best = gb_tune %>% select_best(metric = "rmse")
gb_fit = gb_work %>%
  finalize_workflow(gb_best) %>%
  last_fit(df_split)
gb_test = gb_fit %>% collect_metrics() #metrics evaluated on test sample (b/c last_fit() function) 
gb_predictions = gb_fit %>% collect_predictions() #predictions for test sample (b/c last_fit() function)








## Final model fit on all data 
rf_final = rf_work %>%
  finalize_workflow(rf_best) %>%
  fit(df_data)
rf_final_predictions = rf_final %>% predict(df_data)





gg = dplyr::tibble(predictions = bart_pred,
                   lower = bart_CI[,1],
                   upper = bart_CI[,2],
                   actual = y_test,
                   rf = as.vector(rf_predictions$.pred),
                   lre = as.vector(lre_predictions$.pred)
                   
)
gg = arrange(gg, actual)
gg$index = seq.int(nrow(gg))
gg$Ymean = mean(gg$actual, na.rm = T)


lb1 = paste("R^2 == ", round(rsq, 3))


plot_filtering_estimates2 <- function(df) {
  p <- ggplot(data = gg, aes(x = index)) +
    theme_classic() +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = "indianred"), alpha = 0.5) + #credible intervals 
    geom_line(aes(y = predictions, colour = "red"), size = 0.25, alpha = .9) + #prediction point estimate
    geom_point(aes(y = actual, colour = "black"), size = 0.55, shape = 16, alpha = 0.9) + #actual observation points
    geom_line(aes(y = Ymean, colour = "blue"), size = 0.55, lty = "solid", alpha = 0.9) + #null model (mean only) 
    geom_line(aes(y = rf, colour = "orange"), size = 0.55, lty = "solid", alpha = 0.9) + #null model (mean only) 
    geom_line(aes(y = lre, colour = "green"), size = 0.55, lty = "solid", alpha = 0.9) + #null model (mean only) 
    ylab("Outage Duration (log hours)") + 
    #ylab("Max Cust. Outages (log)") + 
    scale_y_continuous(labels = function(x) paste0(x)) +
    xlab("Outage Index (event x county)") +
    ggtitle(model_name) +
    scale_fill_identity(name = "fill", guide = 'legend', labels = c('95% CI')) + 
    scale_colour_manual(name = 'colour',
                        values = c("green", "orange", 'black',"blue", "red"),
                        #labels = c("RF",'Actual', 'Mean', 'BART'),
                        guide = guide_legend(
                          reverse = T,
                          override.aes = list(
                            linetype = c("solid", "solid", "solid", "solid","solid"),
                            shape = c(NA, NA, NA, NA, 16))
                        )) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(),
          legend.direction = "vertical",
          legend.box = "horizontal",
          legend.position = c(.25, .75) #x and y percentage
    ) 
    #annotate("text", x = quantile(gg$index, 0.8), y = quantile(gg$actual, .05), label=lb1, parse=T, color="red", size = 3) 
  print(p)
}
plot_filtering_estimates2(gg)


