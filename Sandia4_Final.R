#### LOAD PACKAGES and SET DIRECTORY 
options(java.parameters = "-Xmx7g")
library(bartMachine)
set_bart_machine_num_cores(4)
library(tidyverse)
library(tidymodels)
library(tidycensus) 
library(sf)
library(corrplot)
library(viridis)
library(doParallel)
library(terra)
library(stars)
library(raster) #make sure ncdf4 package is installed 
library(lubridate)
library(spdep)
library(xgboost)
#library(sqldf)
library(DBI)
tidymodels_prefer()

setwd("~/Documents/01_VECTOR.nosync/Sandia")

num_cores = detectCores() - 1
#function to un-register parallel processing in doParallel
unregister_dopar = function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

## Import SQL data 
# https://stackoverflow.com/questions/9802680/importing-files-with-extension-sqlite-into-r
# https://cran.r-project.org/web/packages/RSQLite/vignettes/RSQLite.html
mydb = RSQLite::dbConnect(drv = RSQLite::SQLite(), 
                          dbname = "./Data/Outages sqlite/AGM_full_CONUS_outages_merra_noaa_summary_Proc31Jan2022.sqlite")
tables = dbListTables(mydb) #see tables
dbGetQuery(mydb, 'SELECT * FROM outages LIMIT 5')
dbDisconnect(con)
