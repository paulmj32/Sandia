library(tidyverse)
library(tidycensus) 
library(sf)
library(terra)

library(lme4)
options(java.parameters = "-Xmx5g")
library(bartMachine)
set_bart_machine_num_cores(4)

setwd("~/Documents/01_VECTOR.nosync/Sandia")

##### OUTAGE DATA #############################################
# Read in data
#outages_csv = read.csv("SE_states_outage_merra_2018.csv", header = T)
#save(outages_csv, file = "outages.Rda")
load(file = "outages.Rda")

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
## DEM Data
dem_path = "./Data/gt30w100n40_dem/gt30w100n40.dem"
dem = rast(dem_path)
plot(dem)

## NLCD land use
nlcd_path = "./Data/nlcd_2019_land_cover_l48_20210604/nlcd_2019_land_cover_l48_20210604.img"
nlcd = rast(nlcd_path)
plot(nlcd)
