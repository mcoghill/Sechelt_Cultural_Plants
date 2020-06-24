#' Predict Landscape
#'
#' Takes a caret model object and associated covariate rasters to generate a thematic map.
#' In order to process the landscape level prediction the input co-variated are tiled
#' and then mosaic'd together at the end.
#'
#' @param model Location of a `mlr` model object (i.e. path to it.)
#' @param cov   A list of raster files.  These will be loaded as a `stars` object
#' @param tilesize Specify the number of pixels in the `x` and `y` directions for the tiles to be generated.  If your computer is mememory limited use a smaller tile (e.g. 500).
#' @param outDir directory for the output file.
#' @keywords predict, landscape
#' @export
#' @examples
#' 

predict_landscape <- function(model, covariates, tilesize = 500,
                              outDir = "./predicted", type = "raw") {
  
  cov <- sapply(covariates@layers, function(x) x@file@name)
  
  ## create output dir -----------
  dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
  
  ## alternate -- outside of pemgeneratr
  source('./_functions/tile_index.R')
  tiles <- tile_index(cov[1], tilesize)
  
  ## begin loop through tiles -----
  
  ## set up progress messaging
  a <- 0 ## running total of area complete
  ta <- sum(as.numeric(sf::st_area(tiles)))
  wkey <- 0
  
  for (i in 1:nrow(tiles)) {   ## testing first 2 tiles       ##nrow(tiles)) {e
    
    t <- tiles[i, ]  ## get tile
    cat(paste("\nWorking on tile", i, "of", nrow(tiles)))
    
    ## * load tile area---------
    cat("\n...loading new data (from rasters)...")
    r <- stars::read_stars(cov,
                           RasterIO = list(nXOff  = t$offset.x[1] + 1, ## hack -- stars and tile_maker issue??
                                           nYOff  = t$offset.y[1] + 1,
                                           nXSize = t$region.dim.x[1],
                                           nYSize = t$region.dim.y[1]))
    
    ## * update names ---------
    names(r) <- tools::file_path_sans_ext(basename(cov))
    
    ## * convert tile to dataframe ---------
    rsf <- sf::st_as_sf(r, as_points = TRUE)
    
    ## * Test if tile is empty -------------
    na_table <- as.data.frame(sapply(rsf, function(x) all(is.na(x))))
    ## * determine na counts -- this will be used to restore NA values if tile is run.
    na_table$count <- as.data.frame(sapply(rsf, function(x) sum(is.na(x))))[, 1]
    
    ## If any attribute is empty skip the tile
    if(any(na_table[, 1] == TRUE)) {
      
      cat("\nSome variables with all NA values, skipping tile...")
      
    } else {
      
      ## * predict ---------
      
      ## * * Managing NA values ----------------
      ## When some of the values are NA change them to zero
      ## * Identify the attribute with the highest number of NA values.
      na_max <- na_table[na_table$count ==  max(na_table$count), ]
      na_max <- row.names(na_max[1, ]) ## if multiple attributes -- take the first
      ## make a copy of the attribute with the highest number of na values
      rsf_bk <- rsf[, na_max]  ## -- this will be used to restore NA values
      
      ## When some of the values are NA change them to zero -- facilitates pred()
      rsf[is.na(rsf)] <- 0 ## convert NA to zero as the predict function cannot handle NA
      rsf.df <- sf::st_drop_geometry(rsf)
      
      # Default prediction types for most modelling packages use the response/
      # classification type as the default model type. Difficult to hard code all
      # of the options, but by removing the variable it should work out properly
      cat("\n...modelling outcomes (predicting)...")
      if(type != "prob" || class(model) == "ranger") {
        pred <- predict(model, rsf.df)
      } else {
        pred <- predict(model, rsf.df, type = type)
      }
      
      ## Restore NA values
      if(class(model) == "ranger") {
        pred_dat <- pred$predictions
        pred_dat[is.na(st_drop_geometry(rsf_bk)[, 1])] <- NA
        pred <- pred_dat
        
      } else {
        pred_dat <- pred
        pred_dat[is.na(st_drop_geometry(rsf_bk)[, 1])] <- NA
        pred <- pred_dat
      }
      
      ## * geo-link predicted values ---------
      r_out <- cbind(rsf, pred)
      
      ## layers to keep (i.e. newly predicted layers)
      keep <- setdiff(names(r_out), names(r))
      keep <- keep[-length(keep)] ## drops the last entry (geometry field, not a name)
      
      r_out <- r_out %>% dplyr::select(all_of(keep))
      
      ## Save the names of the model response -----
      ## The levels are in the multiclass 'response'
      if(wkey == 0) {
        if(type != "prob") {
          respNames <- levels(r_out$pred) ## this becomes the dictionary to describe the raster values
        } else {
          respNames <- keep ## this becomes the dictionary to describe the raster values
        }
        wkey <- 1
        write.csv(respNames, file.path(outDir, "response_names.csv"), row.names = TRUE)
      }
      
      ## change the text values to numeric values.
      if(type != "prob") {
        r_out$pred <- as.numeric(r_out$pred)
      }
      
      ## Set up subdirectories for rastertile outputs
      cat("\n...exporting raster tiles...")
      
      ## * save tile (each pred item saved) ---------
      for (j in 1:length(keep)) {
        dir.create(file.path(outDir, keep[j]), showWarnings = FALSE)
        out <- stars::st_rasterize(r_out[j],
                                   template = r[1])
        update <- ifelse(file.exists(paste0(outDir,"/",
                                            keep[j], "/",             #sub-directoy
                                            keep[j], "_", i, ".tif")), TRUE, FALSE)
        stars::write_stars(out,
                           paste0(outDir,"/",
                                  keep[j], "/",             #sub-directoy
                                  keep[j], "_", i, ".tif"), 
                           update = update) #tile name
      }
      
      ## * report progress -----
      a <- a + as.numeric(sf::st_area(t))
      cat(paste0("\n", round(a / ta * 100, 1), "% completed at ", format(Sys.time(), "%X %b %d %Y")))
      
    } ## end if statement -- for when tile is empty
  } ## END LOOP -------------
  
  cat("\nAll predicted tiles generated")
  
  ## Mosaic Tiles ---------------
  
  cat("\nGenerating raster mosaics")
  for (k in keep) {
    # get list of tiles
    #k = "response" # testing
    r_tiles <- list.files(paste(outDir, k, sep = "/"),
                          pattern = ".tif$",
                          full.names = TRUE)
    
    ## mosaic
    mos <- gdalUtils::mosaic_rasters(
      gdalfile = r_tiles, ## list of rasters to mosaic
      dst_dataset = paste0(outDir, "/", k, ".tif"),  #output: dir and filename
      output_Raster = TRUE) %>% ## saves the raster (not just a virtual raster)
      rast() %>% 
      mask(rast(subset(covariates, 1)))
  }
  
  if(length(keep) == 1) {
    return(mos)
  } else {
    return(list.files(file.path(outDir), pattern = paste0(keep, ".tif", collapse = "|"), full.names = TRUE))
  }
  
} ### end function
