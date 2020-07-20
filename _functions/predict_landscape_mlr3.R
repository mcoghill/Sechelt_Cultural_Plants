#' Predict Landscape
#'
#' Takes a mlr3 model object and associated covariate rasters to generate a 
#' thematic map.In order to process the landscape level prediction the input 
#' covariates are tiled and then mosaicked together at the end.
#'
#' @param model Location of a `mlr` model object (i.e. path to it.)
#' @param cov   A list of raster files.  These will be loaded as a `stars` object
#' @param tilesize Specify the number of pixels in the `x` and `y` directions for the tiles to be generated.  If your computer is mememory limited use a smaller tile (e.g. 500).
#' @param outDir directory for the output file.
#' @keywords predict, landscape, mlr3, tile
#' @export
#' @examples
#' 

predict_landscape_mlr3 <- function(
  learner, 
  task,
  covariates, 
  tilesize = 500,
  outDir = "./predicted", 
  type = "raw", 
  mask_layer) {
  
  source('./_functions/tile_index.R')
  
  if(!is.character(mask_layer) || missing(mask_layer)) 
    stop("The mask_layer variable needs to be a character vector of length 1 
    \ridentifying the layer to subset from the 'covariates' object")
  
  # Get directories for each raster layer
  cov <- covariates@ptr$filenames
  
  ## create output dir -----------
  dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
  
  ## Use the tile_index function to generate all tiles
  tiles <- tile_index(cov[1], tilesize)
  
  ## set up progress messaging
  a <- 0 ## running total of area complete
  ta <- sum(as.numeric(sf::st_area(tiles)))
  tiles_keep <- NULL
  
  ## begin loop through tiles -----
  tile_files <- foreach(i = 1:nrow(tiles), .combine = c) %do% {
    
    t <- tiles[i, ]  ## get tile
    a <- a + as.numeric(sf::st_area(t))
    cat(paste("\nWorking on tile", i, "of", nrow(tiles)))
    
    ## Do a test run (work in progress, may not speed things up)
    r <- stars::read_stars(cov[1],
                           RasterIO = list(nXOff  = t$offset.x[1] + 1, 
                                           nYOff  = t$offset.y[1] + 1,
                                           nXSize = t$region.dim.x[1],
                                           nYSize = t$region.dim.y[1]))
    
    if(!any(sapply(r, function(x) all(is.na(x))))) {
      
      ## load tile raster data ---------
      cat("\n...Loading new data (from rasters)...")
      r <- stars::read_stars(cov,
                             RasterIO = list(nXOff  = t$offset.x[1] + 1, 
                                             nYOff  = t$offset.y[1] + 1,
                                             nXSize = t$region.dim.x[1],
                                             nYSize = t$region.dim.y[1])) %>% 
        magrittr::set_names(tools::file_path_sans_ext(names(.)))
      cat("done!")
      
    }
    
    ## If any attribute is empty skip the tile, else continue (backwards here)
    if(!any(sapply(r, function(x) all(is.na(x))))) {
      
      ## * update names ---------
      tiles_keep <- c(tiles_keep, i)
      
      ## * convert tile to dataframe ---------
      ## Note: can also simply do sf::st_as_sf() call, but was running into memory 
      ## issues and it was also slower 
      rsf <- as.data.frame(r) %>%
        sf::st_as_sf(coords = c("x", "y")) %>%
        filter_at(vars(names(r)), any_vars(!is.na(.))) %>%
        replace(is.na(.), 0)
      
      ## * predict ---------
      cat("\n...Predicting outcomes...")
      pred <- learner$predict_newdata(newdata = sf::st_drop_geometry(rsf), task = task)
      
      if(type == "prob") {
        pred_dat <- as.data.frame(pred$prob)
      } else {
        pred_dat <- data.frame(pred = pred$response)
      }
      
      ## * geo-link predicted values ---------
      r_out <- sf::st_sf(pred_dat, sf::st_geometry(rsf))
      
      ## layers to keep (i.e. newly predicted layers)
      keep <- names(pred_dat)
      
      ## Save the names of the model response -----
      ## Only needs to be completed for the first iteration
      ## The levels are in the multiclass 'response'
      if(length(tiles_keep) == 1) {
        if(type != "prob") {
          if(is.numeric(r_out$pred)) {
            respNames <- "pred"
          } else respNames <- levels(r_out$pred)
        } else {
          respNames <- keep
        }
        ## this becomes the dictionary to describe the raster values
        write.csv(respNames, file.path(outDir, "response_names.csv"), row.names = TRUE)
      }
      
      ## change the text values to numeric values.
      if(type != "prob") {
        r_out$pred <- as.numeric(r_out$pred)
      }
      
      ## Set up subdirectories for rastertile outputs
      cat("\n...Exporting raster tiles...")
      
      ## * save tile (each pred item saved) ---------
      out_files <- foreach(j = 1:length(keep), .combine = c) %do% {
        dir.create(file.path(outDir, keep[j]), showWarnings = FALSE)
        write_path <- file.path(outDir, keep[j], paste0(keep[j], "_", i, ".tif"))
        if(file.exists(write_path)) unlink(write_path)
        out <- stars::st_rasterize(r_out[j], template = r[1], file = write_path)
        return(write_path)
      }
      
      ## end if statement -- for when tile is not empty
      
    } else {
      
      cat("\nSome variables with all NA values, skipping tile...")
      out_files <- NULL
      
    }
    
    ## * report progress -----
    cat(paste0("\n", round(a / ta * 100, 1), "% completed at ", format(Sys.time(), "%X %b %d %Y"), "\n"))
    return(out_files)
    
  } ## END LOOP -------------
  
  cat("\nAll predicted tiles generated")
  
  ## Mosaic Tiles ---------------
  
  cat("\nGenerating raster mosaics\n")
  
  def_ops <- capture.output(terraOptions())
  terraOptions(progress = 0)
  
  resamp_method <- if(learner$predict_type == "prob" || learner$task_type == "regr") {
    "bilinear"
  } else if(learner$predict_type != "prob" || learner$task_type != "regr") {
    "near"
  } else "bilinear"
  
  for(k in unique(dirname(tile_files))) {
    
    cat(paste("\nMosaicking", basename(k), "tiles"))
    
    r_tiles <- list.files(
      k, pattern = paste0("^", basename(k), "_", tiles_keep, ".tif$", collapse = "|"),
      full.names = TRUE)
    
    # This will overwrite temp files, saving storage space on a PC
    temp <- if(file.exists(grep(".tif$", tmpFiles(), value = TRUE)[1])) {
      grep(".tif$", tmpFiles(), value = TRUE)[1]
    } else ""
    
    ## mosaic
    mos <- foreach(i = r_tiles, .final = function(x) 
      do.call(merge, c(x, filename = paste0(k, ".tif"), overwrite = TRUE))) %do% {
        rast(i)
      } %>% 
      terra::resample(subset(covariates, mask_layer$layer), method = resamp_method, filename = temp, overwrite = TRUE) %>% 
      magrittr::set_names(basename(k)) %>% 
      mask(subset(covariates, mask_layer$layer), filename = paste0(k, ".tif"), overwrite = TRUE)
  }
  
  terraOptions(progress = readr::parse_number(grep("progress", def_ops, value = TRUE)))
  
  if(length(keep) == 1) {
    return(mos)
  } else {
    return(c(file.path(outDir, "response_names.csv"), list.files(
      file.path(outDir), pattern = paste0(keep, ".tif", collapse = "|"), 
      full.names = TRUE)))
  }
} ### end function
