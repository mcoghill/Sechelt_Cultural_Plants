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
  type = "raw", mask = NULL) {
  
  source('./_functions/tile_index.R')
  
  lapply(c("tidyverse", "stars", "terra"), 
         require, character.only = TRUE)
  
  if(missing(model) || !class(model) %in% c("train", "C5.0")) 
    stop("A valid 'model' object from the caret package is required to proceed.")
  
  if(missing(model) || (!class(covariates) %in% "SpatRaster"))
    stop("A SpatRaster object of the raster covariates is required to run predictions.
       Additionally, the layers should be the same as the model variables.")
  
  if(is.null(mask)) mask <- 1L
  
  # Create a polygon of the masking layer
  mask_poly <- terra::as.polygons(subset(covariates, mask) * 0) %>% 
    as("Spatial") %>% 
    sf::st_as_sfc() %>% 
    sf::st_make_valid()
  
  # Create directories to store results in
  dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
  
  # Get the raster file paths for the covariates
  cov <- sources(covariates)[, "source"]
  
  # If rasters are loaded in memory, save them as .tif files elsewhere
  if(any(cov == "")) {
    cov_names <- names(covariates)
    covariates <- terra::writeRaster(
      covariates, filename = tempfile(pattern = names(covariates), fileext = ".tif")) %>% 
      stats::setNames(cov_names)
    cov <- sources(covariates)[, "source"]
  }
  
  # Generate tile index and intersect with polygon to generate area for sampling
  tiles <- tile_index(covariates, tilesize) %>% 
    sf::st_intersection(mask_poly)
  
  # Set up progress messaging
  # Get total area of all tiles
  ta <- sum(as.numeric(sf::st_area(tiles)))
  
  # Create index of which tiles have data in them
  tiles_keep <- NULL
  
  # Begin loop through tiles
  tile_files <- sapply(1:nrow(tiles), function(i) {
    
    # Subset a tile
    t <- tiles[i, ]
    
    # Get total area of processed area + new tile area
    a <- sum(as.numeric(sf::st_area(tiles[i:1, ])))
    cat(paste("\nWorking on tile", i, "of", nrow(tiles), "\n"))
    
    # Do a test run on a single layer, if any variable is all NA then return to
    # top of loop
    r <- stars::read_stars(
      cov[1], RasterIO = list(nXOff  = t$offset.x[1] + 1, 
                              nYOff  = t$offset.y[1] + 1,
                              nXSize = t$region.dim.x[1],
                              nYSize = t$region.dim.y[1]),
      proxy = FALSE)
    
    if(!any(sapply(r, function(x) all(is.na(x))))) {
      # Load all tile data from each raster, if any variable is all NA then
      # return to top of loop. Handle errors where necessary
      cat("\r...Loading raster data...")
      r <- tryCatch({
        stars::read_stars(
          cov, RasterIO = list(nXOff  = t$offset.x[1] + 1,
                               nYOff  = t$offset.y[1] + 1,
                               nXSize = t$region.dim.x[1],
                               nYSize = t$region.dim.y[1]))
      }, error = function(e) {
        lapply(cov, stars::read_stars, RasterIO = list(
          nXOff  = t$offset.x[1] + 1, 
          nYOff  = t$offset.y[1] + 1,
          nXSize = t$region.dim.x[1],
          nYSize = t$region.dim.y[1])) %>% 
          lapply("[[", 1) %>%
          stars::st_as_stars(dimensions = stars::st_dimensions(r), 
                             coordinates = st_coordinates(r))}) %>% 
        stats::setNames(names(covariates))
      
      cat("done!\n")
      
    }
    
    if(any(sapply(r, function(x) all(is.na(x))))) {
      
      cat("\rSome variables with all NA values, skipping tile...\n")
      out_files <- NULL
      
    } else {
      # If tile has data in all variables, then add this to the "keep" index
      tiles_keep <<- c(tiles_keep, i)
      
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
      cat("\r...Predicting outcomes...")
      if(type != "prob") {
        pred <- predict(model, sf::st_drop_geometry(rsf))
      } else {
        pred <- cbind(predict(model, sf::st_drop_geometry(rsf), type = type), 
                      pred = predict(model, sf::st_drop_geometry(rsf)))
      }
      cat("done!\n")
      
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
      cat("\r...Exporting raster tiles...\n")
      
      # Save tile (each pred item saved)
      out_files <- sapply(1:length(keep), function(j) {
        dir.create(file.path(outDir, keep[j]), showWarnings = FALSE)
        write_dir <- file.path(outDir, keep[j], paste0(keep[j], "_", i, ".tif"))
        if(file.exists(write_dir)) unlink(write_dir)
        out <- stars::st_rasterize(r_out[j], template = r[1], file = write_dir)
        return(write_dir)
      })
    }
    
    # Report progress
    cat(paste0("\r", round(a / ta * 100, 1), "% completed at ", 
               format(Sys.time(), "%X %b %d %Y"), "\n"))
    return(out_files)
    
  }) %>% as.character()
  
  cat("\nAll predicted tiles generated\n")
  
  # Mosaic Tiles
  cat("\rGenerating raster mosaics\n\n")
  
  # Don't want to display a bunch of progress bars here, so turn that off for
  # the time being
  def_ops <- terra:::spatOptions()$progress
  terraOptions(progress = 0)
  
  pred_out <- do.call(c, lapply(unique(dirname(tile_files)), function(k) {
    
    cat(paste("\rMosaicking", basename(k), "tiles", "\n"))
    
    # In order to properly mask the layer, the CRS and extents need to match 
    # perfectly, hence the resampling step
    if(basename(k) == "pred") {
      resamp_method <- "near"
      wopt <- list(datatype = "INT2S")
    } else {
      resamp_method <- "bilinear"
      wopt <- list(datatype = "FLT4S")
    }
    
    r_tiles <- list.files(
      k, pattern = paste0("^", basename(k), "_", tiles_keep, ".tif$", collapse = "|"),
      full.names = TRUE)
    
    # Mosaic, resample, mask and save as temp file
    mos <- {if(length(r_tiles) > 1) {
      do.call(mosaic, lapply(r_tiles, terra::rast))
    } else terra::rast(r_tiles)} %>% 
      # terra::resample(subset(covariates, mask), method = resamp_method) %>% 
      stats::setNames(basename(k)) %>% 
      # terra::mask(subset(covariates, mask), 
      #             filename = tempfile(pattern = basename(k), fileext = ".tif"), 
      #             overwrite = TRUE, wopt = wopt)
      terra::writeRaster(tempfile(pattern = basename(k), fileext = ".tif"),
                         overwrite = TRUE, wopt = wopt)
  }))
  
  # Reapply default progress bar options
  terraOptions(progress = def_ops)
  
  return(pred_out)
  
}

#END