## Root zone data
# https://www.nrcs.usda.gov/wps/portal/nrcs/detail/soils/survey/geo/?cid=nrcs142p2_053628
# https://gdg.sc.egov.usda.gov/GDGHome_DirectDownLoad.aspx --> soil geographic database 

# NOTE: Too large to read in CONUS dataset (8GB Mac) ... need to run for each state 
# RZ_conus_gdb =  "/Users/paulmj/Downloads/gSSURGO_CONUS/gSSURGO_CONUS.gdb"
# RZ_conus = sf::st_read(dsn = RZ_conus_gdb, layer = "MUPOLYGON")
# 
# RZ_conus_sql = sf::st_read(dsn = RZ_conus_gdb,
#                            query = 'SELECT * FROM "MUPOLYGON" WHERE FID = 1',
#                            layer = "MUPOLYGON")


try_list = state_list[3:32]
# Iterate: read in state shapefile from state_list, group by key, join value table, rasterize, save raster 
for (i in try_list){
  print(i)
  temp_wd = paste("/Users/paulmj/Downloads/July2020_gSSURGO_by_State/gSSURGO_", i, sep = "")   # change working directory to file 
  setwd(temp_wd) #set working directory 
  temp_gdb = paste("gSSURGO_", i, ".gdb", sep = "") # define geodatabase 
  temp_sf = sf::st_read(dsn = temp_gdb, layer = "MUPOLYGON") # read in MUPOLYGON layer 
  temp_group = temp_sf %>% group_by(MUKEY) %>% summarise(n = n()) #group by MUKEY 
  temp_val = sf::st_read(dsn = temp_gdb, layer = "Valu1") # read in Value1 table
  temp_join = temp_group %>% left_join(temp_val, by = c("MUKEY" = "mukey")) # join values
  temp_ras = st_rasterize(temp_join["rootznemc"], dx = 100, dy = 100) #100m x 100km resolution raster
  temp_write = paste(i, "_rootznemc100.tif", sep = "") #name of raster
  #temp_ras= st_rasterize(temp_join["rootznemc"], dx = 1000, dy = 1000) #1000m x 1000km resolution raster
  #temp_write = paste(i, "_rootznemc1000.tif", sep = "") #name of raster
  write_stars(temp_ras, dsn = temp_write) # write raster in working directory 
}


