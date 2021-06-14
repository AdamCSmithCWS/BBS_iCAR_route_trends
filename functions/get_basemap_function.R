get_basemap <- function(strata_type = NULL,
                        transform_laea = TRUE,
                        append_area_weights = FALSE){
  
  ## returns an sf (simple features) geometry object
  ## of one of hte bbsBayes strata options
  ## optionally reprojected into a Lambert Equal Area coordinate reference system
  
  if (isFALSE(is.element(strata_type, c("state", "bcr", 
                                        "latlong", "bbs_cws", "bbs_usgs")))) {
    stop("Invalid stratification specified, choose one of state, bcr, latlong, bbs_cws, or bbs_usgs")
    return(NULL)
  }
  require(bbsBayes)
  require(sf)
  require(dplyr)
  
  strat_opts <- c("BBS_ProvState_strata",
                  "BBS_BCR_strata",
                  "BBS_LatLong_strata",
                  "BBS_CWS_strata",
                  "BBS_USGS_strata")
  names(strat_opts) <- c("state", 
                         "bcr", 
                         "latlong", 
                         "bbs_cws", 
                         "bbs_usgs")
  
  
  laea = st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area coord reference system
  
  
  locat = system.file("maps",
                      package = "bbsBayes")
  
  map.file = strat_opts[strata_type] # the most standard BBS stratification, used here just to provide some spatial limits
  
  # reading in the strata spatial file (GIS shapefile)
  strata_map = read_sf(dsn = locat,
                       layer = map.file)
  
  if(transform_laea){
    laea = st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area coord reference system
    
    strata_map = st_transform(strata_map,crs = laea) #reprojecting the geographic coordinate file to an equal area projection
  }
  
  if(append_area_weights){
    strat_df <- get_composite_regions(strata_type)
    strata_map <- left_join(strata_map,strat_df,by = c("ST_12" = "region"))
  }
  
  
  return(strata_map)
}

