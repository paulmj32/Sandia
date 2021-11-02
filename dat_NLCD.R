# 
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