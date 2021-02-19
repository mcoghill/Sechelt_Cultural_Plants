#' Summarise the results of a DSMART spatial disaggregation.
#' 
#' \code{summarise} summarises the results of the spatial disaggregation of a 
#' polygon soil map in several ways. First, it computes the probabilities of 
#' occurrence of the soil classes that occur across all the map's map units. 
#' Second, it ranks the soil class predictions according to their probabilities 
#' of occurrence and maps the \emph{n}-most-probable soil classes at each grid 
#' cell, and their probabilities. Finally, it computes Shannon's entropy on the
#' class probabilities, and the degree of confusion between the most probable 
#' and second-most-probable soil classes.
#' 
#' @param realisations A \code{SpatRaster} where each layer contains one 
#'   realisation of the soil class distribution across the soil map area, as 
#'   produced by \code{\link{disaggregate}}. If probabilistic predictions are
#'   used (\code{type = "prob"}), a list of SpatRaster objects with predicted
#'   class probabilities must be passed.
#' @param lookup A two-column \code{data.frame} containing a mapping between the
#'   integer soil class codes in the layers of \code{realisations}, and the soil
#'   class codes defined by the map unit composition \code{data.frame} used as
#'   an argument to \code{disaggregate} and \code{dsmart}. \code{lookup} is the
#'   same lookup table that is produced by \code{dsmart}. First column is the
#'   soil class code of the map unit composition; second column is the integer
#'   soil class code.
#' @param n.realisations An integer that identifies the number of realisations
#'   of the soil class distribution that were computed by \code{disaggregate}.
#'   Default value is \code{terra::nlyr(realisations)}.
#' @param nprob At any location, disaggregated soil class predictions can be 
#'   ranked according to their probabilities of occurence. \code{rdsmart} can 
#'   map the class predictions, and their probabilities, at any rank. 
#'   \code{nprob} is an integer that identifies the number of probability ranks 
#'   to map. For example, if \code{n = 3}, DSMART will map the first-, second- 
#'   and third-most-probable soil classes and their probabilities of occurrence.
#' @param outputdir A character string that identifies the location of the main 
#'   output directory. The folder \code{output} and its subfolders will be 
#'   placed here. Default is the current working directory, \code{getwd()}.
#' @param stub \emph{optional} A character string that identifies a short name
#'   that will be prepended to all output.
#' @param prob A character vector to specify the type of the predictions to be 
#'   summarised. By default, raw class predictions are used. If set to "prob",
#'   probabilistic predictions are used.
#'
#' @return A list that contains metadata about the current run of
#'   \code{summarise}.
#'   
#' @references McBratney, A.B., Mendonca Santos, M. de L., Minasny, B., 2003. On
#'   digital soil mapping. Geoderma 117, 3--52. doi: 
#'   \href{http://dx.doi.org/10.1016/S0016-7061(03)00223-4}{10.1016/S0016-7061(03)00223-4}
#'   
#'   Odgers, N.P., McBratney, A.B., Minasny, B., Sun, W., Clifford, D., 2014. 
#'   DSMART: An algorithm to spatially disaggregate soil map units, \emph{in:} 
#'   Arrouays, D., McKenzie, N.J., Hempel, J.W., Richer de Forges, A., 
#'   McBratney, A.B. (Eds.), GlobalSoilMap: Basis of the Global Spatial Soil 
#'   Information System. Taylor & Francis, London, pp. 261--266.
#'   
#'   Odgers, N.P., Sun, W., McBratney, A.B., Minasny, B., Clifford, D., 2014. 
#'   Disaggregating and harmonising soil map units through resampled 
#'   classification trees. Geoderma 214, 91--100. doi: 
#'   \href{http://dx.doi.org/10.1016/j.geoderma.2013.09.024}{10.1016/j.geoderma.2013.09.024}
#'   
#' @examples
#' # Load datasets
#' data(dalrymple_lookup)
#' data(dalrymple_realisations)
#' 
#' # Summarise
#' summarise(dalrymple_realisations, dalrymple_lookup, nprob = 5)
#' 
#' @export

summarise <- function(
  realisations, lookup, 
  n.realisations = ifelse(is.list(realisations), length(realisations), terra::nlyr(realisations)), 
  nprob = 3, outputdir = getwd(), stub = NULL, type = "raw") {
  
  # Requires the following packages:
  sapply(c("tidyverse", "terra"), 
         require, character.only = TRUE)[0]
  
  # Create list to store output
  output <- base::list()
  
  # Save start time
  output$timing <- base::list(start = Sys.time())
  
  # Check arguments before proceeding
  messages <- c("Attention is required with the following arguments:\n")
  
  if(type != "prob") {
    if(!(class(realisations) == "SpatRaster")) {
      messages <- append(messages, "'realisations': Not a valid SpatRaster\n")
    } 
  } else {
    if(is.list(realisations) == FALSE) {
      messages <- append(messages, "'realisations' must be a list of SpatRaster objects when probabilistic predictions are used.'.\n")
    } else if(sum(unlist(lapply(realisations, function(x) class(x) != "SpatRaster"))) > 0) {
      messages <- append(messages, "'realisations' must be a list of SpatRaster objects when probabilistic predictions are used.'.\n")
    }
  }
  if(!(class(lookup) == "data.frame")) {
    messages <- append(messages, "'lookup': Not a valid data.frame.\n")
  }
  if(n.realisations <= 0) {
    messages <- append(messages, "'n.realisations': Value must be greater than 0.\n")
  }
  if(nprob <= 0) {
    messages <- append(messages, "'nprob': Value must be greater than 0.\n")
  }
  if(!(file.exists(outputdir))) {
    messages <- append(messages, "'outputdir': Output directory does not exist.")
  }
  if(length(messages) > 1) {
    stop(messages)
  }
  if(is.null(stub)) {
    stub <- ""
  } else if(stub == "") {
    stub <- ""
  } else if(!(substr(stub, nchar(stub), nchar(stub)) == "_")) {
    stub <- paste0(stub, "_")
  }
  
  # Save function call
  output$call <- base::match.call()
  
  # Save parameters
  output$parameters <- base::list(
    n.realisations = n.realisations, 
    nprob = nprob, stub = stub, type = type)
  
  # Set up output directories
  outputdir <- file.path(outputdir)
  ndirs <- length(dir(outputdir, pattern = "output"))
  if(ndirs > 1) {
    last <- max(as.numeric(
      gsub("\\D", "", dir(outputdir, pattern = "output"))), na.rm = TRUE)
    subdir <- paste0("output_", formatC(last, width = 3, format = "d", flag = 0))
  } else {
    subdir <- "output"
  }
  dir.create(file.path(outputdir, subdir), showWarnings = FALSE)
  dir.create(file.path(outputdir, subdir, "probabilities"), showWarnings = FALSE)
  dir.create(file.path(outputdir, subdir, "mostprobable"), showWarnings = FALSE)
  
  # Save output locations
  output$locations <- base::list(
    root = file.path(outputdir, subdir), 
    probabilities = file.path(outputdir, subdir, "probabilities"), 
    mostprobable = file.path(outputdir, subdir, "mostprobable"))
  
  # Don't want to display a bunch of progress bars here, so turn that off for
  # the time being
  def_ops <- terra:::spatOptions()$progress
  terraOptions(progress = 0)
  
  # Make sure lookup table column names are correct
  names(lookup) <- c("name", "code")
  
  # If raw predictions are used, calculate class probabilities by counting.
  cat("\nComputing class probabilities from realisations")
  if(type != "prob") {
    
    # Compute counts
    counts <- terra::app(realisations, function(x) {
      if(is.na(sum(x))) {
        rep(NA, nrow(lookup))
      } else {
        tabulate(x, nbins = nrow(lookup))
      }
    }, wopt = list(names = lookup$name))
    
    # Compute probabilities and write probabilities to raster files
    probs <- (counts / n.realisations) %>% 
      terra::writeRaster(filename = file.path(
        outputdir, subdir, "probabilities", 
        paste0(stub, "prob_", names(.), ".tif")), 
        overwrite = TRUE)
  } else {
    # If probabilistic predictions are used, calculate class probabilities by 
    # averaging the predicted probabilities across the realisations.
    if(length(realisations) == 1 | n.realisations == 1) {
      # If only one realisation is used, no averaging is needed.
      probs <- terra::writeRaster(realisations[[1]], filename = file.path(
        outputdir, subdir, "probabilities", 
        paste0(stub, "prob_", names(realisations[[1]]), ".tif")), 
        overwrite = TRUE)
    } else {
      probs <- rast(lapply(1:nrow(lookup), function(i) {
        rlist <- mean(rast(lapply(realisations, "[[", i)), na.rm = TRUE) %>% 
          terra::writeRaster(file.path(
            outputdir, subdir, "probabilities", 
            paste0(stub, "prob_", lookup$name[which(lookup$code == i)], ".tif")),
            overwrite = TRUE, wopt = list(names = lookup$name[which(lookup$code == i)]))
        # return(terra::app(rlist, fun = "mean", na.rm = TRUE, filename = file.path(
        #   outputdir, subdir, "probabilities", 
        #   paste0(stub, "prob_", lookup$name[which(lookup$code == i)], ".tif")), 
        #   overwrite = TRUE))
      })) %>% stats::setNames(lookup$name)
    }
  }
  
  # Compute the class indices of the n-most-probable soil classes
  cat("\nGenerating most probable classification rasters")
  ind.names <- paste0(stub, "mostprob_", 
                      formatC(1:nlyr(probs), width = nchar(nrow(lookup)), 
                              format = "d", flag = "0"), "_class")
  ordered.indices <- if(type != "prob") {
    # If raw class predictions are used, use "counts" for indicing.
    terra::app(counts, function(x) {
      if(is.na(sum(x))) {
        rep(NA, nlyr(probs))
      } else {
        order(x, decreasing = TRUE, na.last = TRUE)
      }
    }, wopt = list(names = ind.names, datatype = "INT2S"))[[1:nprob]]
  } else {
    # If probabilistic predictions are used, use "probs" for indicing.
    terra::app(probs, function(x) {
      if(is.na(sum(x))) {
        rep(NA, nlyr(probs))
      } else {
        order(x, decreasing = TRUE, na.last = TRUE)
      }
    }, wopt = list(names = ind.names, datatype = "INT2S"))[[1:nprob]]
  } 
  
  # Write ith-most-probable soil class raster to file. Ordering will fill areas
  # that are supposed to be NA, so mask NA areas with the mask layer
  ordered.indices <- terra::writeRaster(
    ordered.indices, filename = file.path(
      outputdir, subdir, "mostprobable", 
      paste0(ind.names[1:nprob], ".tif")), 
    overwrite = TRUE, wopt = list(datatype = "INT2S"))
  
  # Compute the class probabilities of the n-most-probable soil classes
  cat("\nCalculating probabilities of classification rasters being accurate")
  prob.names <- paste0(stub, "mostprob_", formatC(
    1:nlyr(probs), width = nchar(nrow(lookup)), format = "d", flag = "0"), 
    "_probs")
  ordered.probs <- terra::app(probs, function(x) {
    if(is.na(sum(x))) {
      rep(NA, nlyr(probs))
    } else {
      sort(x, decreasing = TRUE, na.last = TRUE)
    }
  }, wopt = list(names = prob.names))[[1:max(2, nprob)]]
  
  # Compute the confusion index on the class probabilities
  cat("\nGenerating confusion index raster from the class probabilities")
  confusion <- (1 - (ordered.probs[[1]] - ordered.probs[[2]])) %>% 
    terra::writeRaster(
      filename = file.path(outputdir, subdir, "mostprobable", 
                           paste0(stub, "confusion.tif")), overwrite = TRUE,
      wopt = list(names = "confusion"))
  
  # Compute Shannon's entropy on the class probabilities
  # See rationale in section 2.5 of Kempen et al. (2009) "Updating the 1:50,000
  # Dutch soil map"
  cat("\nGenerating Shannon diversity raster (entropy) from the class probabilities")
  shannon <- terra::app(ordered.probs, function(x) {
    -sum(x * (log(x, base = length(x))), na.rm = TRUE)}) %>% 
    terra::mask(confusion, filename = file.path(
      outputdir, subdir, "mostprobable", 
      paste0(stub, "shannon.tif")), overwrite = TRUE, 
    wopt = list(names = "shannon"))
  
  # This will hopefully work in a newer "terra" version, will be faster:
  # shannon <- -sum(
  #   ordered.probs * (log(ordered.probs, base = nlyr(ordered.probs))), 
  #   na.rm = TRUE) %>% 
  #   terra::writeRaster(file.path(
  #     outputdir, subdir, "mostprobable", 
  #     paste0(stub, "shannon.tif")), overwrite = TRUE, 
  #     wopt = list(names = "shannon"))
  
  # Write ith-most-probable soil class probability raster to file
  cat("\nWriting probability raster(s)")
  ordered.probs <- terra::writeRaster(
    ordered.probs[[1:nprob]], filename = file.path(
      outputdir, subdir, "mostprobable", 
      paste0(prob.names[1:nprob], ".tif")), 
    overwrite = TRUE)
  
  # Remove temporary files from system
  suppressWarnings(terra::tmpFiles(old = TRUE, remove = TRUE))
  
  # Save finish time
  output$timing$finish <- Sys.time()
  
  # Reapply default progress bar options
  terraOptions(progress = def_ops)
  
  # Return output
  return(output)
}

#END