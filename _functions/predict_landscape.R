#' Predict Landscape
#'
#' Takes a caret model object and associated covariate rasters to generate a 
#' thematic map. In order to process the landscape level prediction the input 
#' covariates are tiled and then mosaicked together at the end.
#'
#' @param model A `caret` model object
#' @param covariates A \code{SpatRaster} object, containing the same variables as
#' the \code{model}
#' @param tilesize Numerical vector of length 1 to specify the number of pixels 
#' in the `x` and `y` directions for the tiles to be generated.  If your computer 
#' is mememory limited use a smaller tile (e.g. 500).
#' @param outDir Directory for the output file(s).
#' @param type The type of predictions to generate (either "prob" for probabilty
#' estimation, or "raw" for raw predictions).
#' @keywords predict, landscape
#' @export
#' @examples
#' 

predict_landscape <- function(
  model, 
  covariates, 
  tilesize = 500,
  outDir = file.path(getwd(), "predicted_tiles"),
  type = "raw") {
  
  source('./_functions/tile_index.R')
  
  if(missing(model) || class(model) != "train") 
    stop("A valid 'model' object from the caret package is required to proceed.")
  
  if(missing(model) || (!class(covariates) %in% "SpatRaster"))
    stop("A SpatRaster object of the raster covariates is required to run predictions.
       Additionally, the layers should be the same as the model variables.")
  
  # Get the raster file paths for the covariates
  cov <- sources(covariates)[, "source"]
  
  # Count NA values in rasters data to determine best layer to use for masking
  mask_layer <- foreach(i = 1:nlyr(covariates), .combine = rbind) %do% {
    cat(paste0("Counting NA values in ", names(covariates[[i]]), 
               " [", i, " of ", nlyr(covariates), "]\n"))
    new <- subset(covariates, i) * 0
    data.frame(layer = names(new), 
               data_cells = data.frame(freq(new))$count)
  } %>% mutate(data_cells = ifelse(data_cells == ncell(new), 0, data_cells)) %>% 
    dplyr::slice(which.max(data_cells))
  
  # Create directories to store results in
  dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
  
  # Generate tile index
  tiles <- tile_index(cov[1], tilesize)
  
  # Set up progress messaging
  # Get total area of all tiles
  a <- 0
  ta <- sum(as.numeric(sf::st_area(tiles)))
  
  # Create index of which tiles have data in them
  tiles_keep <- NULL
  
  # Begin loop through tiles
  tile_files <- foreach(i = 1:nrow(tiles), .combine = c) %do% {
    
    # Subset a tile
    t <- tiles[i, ]
    
    # Get total area of processed area + new tile area
    a <- a + as.numeric(sf::st_area(t))
    cat(paste("\nWorking on tile", i, "of", nrow(tiles)))
    
    # Do a test run on a single layer, if any variable is all NA then return to
    # top of loop
    r <- stars::read_stars(cov[1],
                           RasterIO = list(nXOff  = t$offset.x[1] + 1, 
                                           nYOff  = t$offset.y[1] + 1,
                                           nXSize = t$region.dim.x[1],
                                           nYSize = t$region.dim.y[1]))
    
    if(!any(sapply(r, function(x) all(is.na(x))))) {
      # Load all tile data from each raster, if any variable is all NA then
      # return to top of loop
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
      # If tile has data in all variables, then add this to the "keep" index
      tiles_keep <- c(tiles_keep, i)
      
      # Convert tile to sf dataframe, only keep rows that are not all NA, 
      # and replace any NA values with 0
      rowAny <- function(x) rowSums(!is.na(x)) > 0
      rsf <- as.data.frame(r) %>% 
        sf::st_as_sf(coords = c("x", "y")) %>%
        dplyr::filter(rowAny(across(names(r), ~ .x > 0))) %>%
        replace(is.na(.), 0)
      
      # Old dplyr::filter_at method below
      # rsf <- as.data.frame(r) %>% 
      #   sf::st_as_sf(coords = c("x", "y")) %>%
      #   dplyr::filter_at(vars(names(r)), any_vars(!is.na(.))) %>%
      #   replace(is.na(.), 0)

      # Carry out model prediction and format depending on predict type
      cat("\n...modelling outcomes (predicting)...")
      if(type != "prob") {
        pred <- predict(model, st_drop_geometry(rsf))
      } else {
        pred <- cbind(predict(model, st_drop_geometry(rsf), type = type), 
                      pred = predict(model, st_drop_geometry(rsf)))
      }
      
      # Geo-link predicted values
      r_out <- sf::st_sf(pred, sf::st_geometry(rsf))
      
      # Layers to keep (i.e. newly predicted layers)
      keep <- names(r_out)[-ncol(r_out)]
      r_out <- r_out %>% dplyr::select(all_of(keep))
      
      # Save the names of the model response
      # The levels are in the multiclass 'response'
      if(type != "prob") {
        names(r_out)[1] <- "pred"
        keep <- names(r_out)[-ncol(r_out)]
      }
      
      # This becomes the dictionary to describe the raster values
      if(length(tiles_keep) == 1) {
        if(type != "prob") {
          if(is.numeric(r_out$pred)) {
            respNames <- "pred"
          } else respNames <- levels(r_out$pred)
        } else respNames <- keep
        write.csv(respNames, file.path(outDir, "response_names.csv"), row.names = TRUE)
      }
      
      # Change the text values to numeric values.
      if("pred" %in% names(r_out)) {
        r_out$pred <- as.numeric(r_out$pred)
      }
      
      # Set up subdirectories for rastertile outputs
      cat("\n...Exporting raster tiles...")
      
      # Save tile (each pred item saved)
      out_files <- foreach(j = 1:length(keep), .combine = c) %do% {
        dir.create(file.path(outDir, keep[j]), showWarnings = FALSE)
        write_dir <- file.path(outDir, keep[j], paste0(keep[j], "_", i, ".tif"))
        if(file.exists(write_dir)) unlink(write_dir)
        out <- stars::st_rasterize(r_out[j], template = r[1], file = write_dir)
        return(write_dir)
      }
    }
    
    # Report progress
    cat(paste0("\n", round(a / ta * 100, 1), "% completed at ", 
               format(Sys.time(), "%X %b %d %Y"), "\n"))
    return(out_files)
      
    }
  
  cat("\nAll predicted tiles generated")
  
  # Mosaic Tiles
  cat("\nGenerating raster mosaics")
  
  # Don't want to display a bunch of progress bars here, so turn that off for
  # the time being
  def_ops <- capture.output(terraOptions())
  terraOptions(progress = 0)
  
  pred_out <- foreach(k = unique(dirname(tile_files)), .combine = c) %do% {
    
    cat(paste("\nMosaicking", basename(k), "tiles"))
    
    # In order to properly mask the layer, the CRS and extents need to match 
    # perfectly, hence the resampling step
    resamp_method <- if(basename(k) == "pred") {
      "near"
    } else "bilinear"
    
    r_tiles <- list.files(
      k, pattern = paste0("^", basename(k), "_", tiles_keep, ".tif$", collapse = "|"),
      full.names = TRUE)
    
    # Mosaic, resample, mask and save as temp file
    mos <- foreach(i = r_tiles, .final = function(x) do.call(terra::merge, x)) %do% {
        rast(i)
      } %>% 
      terra::resample(subset(covariates, mask_layer$layer), method = resamp_method) %>% 
      magrittr::set_names(basename(k)) %>% 
      terra::mask(subset(covariates, mask_layer$layer), 
                  filename = tempfile(pattern = basename(k), fileext = ".tif"), 
                  overwrite = TRUE)
    
    # Depending on prediction type, output either a single SpatRaster layer or
    # a list of files
    if(length(keep) == 1) {
      return(mos)
    } else {
      return(terra::sources(mos)[, "source"])
    }
  }
  
  # Reapply default progress bar options
  terraOptions(progress = readr::parse_number(grep("progress", def_ops, value = TRUE)))
  
  return(pred_out)
  
}

#END