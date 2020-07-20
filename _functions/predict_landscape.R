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

predict_landscape <- function(
  model, 
  covariates, 
  tilesize = 500,
  outDir = "./predicted", 
  type = "raw") {
  
  source('./_functions/tile_index.R')
  
  # Count NA values in the covariate data to determine best layer to use for masking later on
  cov <- covariates@ptr$filenames
  
  mask_layer <- foreach(i = 1:nlyr(covariates), .combine = rbind) %do% {
    cat(paste0("Counting NA values in ", names(covariates[[i]]), 
               " [", i, " of ", nlyr(covariates), "]\n"))
    new <- subset(covariates, i) * 0
    data.frame(layer = names(new), 
               data_cells = data.frame(freq(new))$count)
  } %>% mutate(data_cells = ifelse(data_cells == ncell(new), 0, data_cells)) %>% 
    dplyr::slice(which.max(data_cells))
  
  ## create output dir -----------
  dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
  
  tiles <- tile_index(cov[1], tilesize)
  
  ## begin loop through tiles -----
  
  ## set up progress messaging
  a <- 0 ## running total of area complete
  ta <- sum(as.numeric(sf::st_area(tiles)))
  tiles_keep <- NULL
  
  tile_files <- foreach(i = 74:nrow(tiles), .combine = c) %do% {
    
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
      
      ## * load tile area---------
      cat("\n...loading new data (from rasters)...")
      r <- stars::read_stars(cov,
                             RasterIO = list(nXOff  = t$offset.x[1] + 1, 
                                             nYOff  = t$offset.y[1] + 1,
                                             nXSize = t$region.dim.x[1],
                                             nYSize = t$region.dim.y[1])) %>% 
        magrittr::set_names(tools::file_path_sans_ext(names(.)))
      cat("done!")
      
    }
    
    if(any(sapply(r, function(x) all(is.na(x))))) {
      
      cat("\nSome variables with all NA values, skipping tile...")
      out_files <- NULL
      
    } else {
      ## * update names ---------
      tiles_keep <- c(tiles_keep, i)
      
      ## * convert tile to dataframe ---------
      rsf <- as.data.frame(r) %>% 
        sf::st_as_sf(coords = c("x", "y")) %>%
        filter_at(vars(names(r)), any_vars(!is.na(.))) %>%
        replace(is.na(.), 0)

      # Default prediction types for most modelling packages use the response/
      # classification type as the default model type. Difficult to hard code all
      # of the options, but by removing the variable it should work out properly
      cat("\n...modelling outcomes (predicting)...")
      if(type != "prob") {
        pred <- predict(model, st_drop_geometry(rsf))
      } else {
        pred <- cbind(predict(model, st_drop_geometry(rsf), type = type), 
                      pred = predict(model, st_drop_geometry(rsf)))
      }
      
      ## * geo-link predicted values ---------
      r_out <- sf::st_sf(pred, sf::st_geometry(rsf))
      
      ## layers to keep (i.e. newly predicted layers)
      keep <- names(pred)
      
      r_out <- r_out %>% dplyr::select(all_of(keep))
      
      ## Save the names of the model response -----
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
      if("pred" %in% names(r_out)) {
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
    } ## end if statement -- for when tile is empty
    
    ## * report progress -----
    cat(paste0("\n", round(a / ta * 100, 1), "% completed at ", format(Sys.time(), "%X %b %d %Y"), "\n"))
    return(out_files)
      
    } ## END LOOP -------------
  
  cat("\nAll predicted tiles generated")
  
  ## Mosaic Tiles ---------------
  
  cat("\nGenerating raster mosaics")
  
  def_ops <- capture.output(terraOptions())
  terraOptions(progress = 0)
  
  for(k in unique(dirname(tile_files))) {
    
    cat(paste("\nMosaicking", basename(k), "tiles"))
    
    resamp_method <- if(basename(k) == "pred") {
      "near"
    } else "bilinear"
    
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
    return(list.files(file.path(outDir), pattern = paste0(keep, ".tif", collapse = "|"), full.names = TRUE))
  }
  
} ### end function
