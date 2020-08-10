summarise <- function(
  realisations, lookup, 
  n.realisations = ifelse(is.list(realisations), length(realisations), terra::nlyr(realisations)), 
  nprob = 3, outputdir = getwd(), stub = NULL, type = "raw") {
  
  output <- base::list()
  output$timing <- base::list(start = base::date())
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
  
  output$call <- base::match.call()
  output$parameters <- base::list(
    n.realisations = n.realisations, 
    nprob = nprob, stub = stub, type = type)
  outputdir <- file.path(outputdir)
  dir.create(file.path(outputdir, "output"), showWarnings = FALSE)
  dir.create(file.path(outputdir, "output", "probabilities"), showWarnings = FALSE)
  dir.create(file.path(outputdir, "output", "mostprobable"), showWarnings = FALSE)
  output$locations <- base::list(
    root = file.path(outputdir, "output"), 
    probabilities = file.path(outputdir, "output", "probabilities"), 
    mostprobable = file.path(outputdir, "output", "mostprobable"))
  names(lookup) <- c("name", "code")
  param <- nrow(lookup)
  if(type != "prob") {
    mask_layer <- realisations[[1]]
    counts <- terra::app(realisations, tabulate, nbins = param) %>% 
      terra::mask(mask_layer) %>% 
      magrittr::set_names(lookup$name)
    probs <- (counts / n.realisations) %>% 
      terra::writeRaster(filename = file.path(
        outputdir, "output", "probabilities", 
        paste0(stub, "prob_", names(.), ".tif")), 
        overwrite = TRUE)
  } else {
    mask_layer <- realisations[[1]][[1]]
    if(length(realisations) == 1 | n.realisations == 1) {
      probs <- terra::writeRaster(realisations[[1]], filename = file.path(
        outputdir, "output", "probabilities", 
        paste0(stub, "prob_", names(realisations[[1]]), ".tif")), 
        overwrite = TRUE)
    } else {
      probs <- foreach(i = 1:param, .final = function(x) do.call(c, x)) %do% {
        rlist <- foreach(j = 1:n.realisations, .final = function(y) do.call(c, y)) %do% {
          return(realisations[[j]][[i]])
        }
        return(terra::app(rlist, fun = "mean", na.rm = TRUE, filename = file.path(
          outputdir, "output", "probabilities", 
          paste0(stub, "prob_", lookup$name[which(lookup$code == i)], ".tif")), 
          overwrite = TRUE))
      }
    }
  }
  
  ordered.indices <- if(type != "prob") {
    if(nprob > 1) {
      terra::app(counts, order, decreasing = TRUE, na.last = TRUE)[[1:nprob]]
    } else {
      terra::modal(realisations, ties = "first", na.rm = TRUE)
    }
  } else {
    terra::app(probs, order, decreasing = TRUE, na.last = TRUE)[[1:nprob]]
  } 
  
  ordered.indices <- ordered.indices %>% 
    terra::mask(mask_layer) %>% 
    magrittr::set_names(paste0(stub, "mostprob_", 
                               formatC(1:nprob, width = nchar(nrow(lookup)), 
                                       format = "d", flag = "0"), "_class")) %>% 
    terra::writeRaster(filename = file.path(
      outputdir, "output", "mostprobable", 
      paste0(names(.), ".tif")), 
      overwrite = TRUE, wopt = list(datatype = "INT2S"))
  
  ordered.probs <- terra::app(probs, sort, decreasing = TRUE, 
                              na.last = TRUE)[[1:max(2, nprob)]] %>% 
    terra::mask(mask_layer) %>% 
    magrittr::set_names(
      paste0(stub, "mostprob_", formatC(1:max(2, nprob), width = nchar(nrow(lookup)), 
                                        format = "d", flag = "0"), "_probs"))
  
  confusion <- terra::app(ordered.probs, function(x) {
    (1 - (x[[1]] - x[[2]]))
  }, filename = file.path(outputdir, "output", "mostprobable", 
                          paste0(stub, "confusion.tif")), overwrite = TRUE)
  shannon <- terra::app(ordered.probs, function(x) {
    -sum(x * (log(x, base = length(x))), na.rm = TRUE)}) %>% 
    terra::mask(mask_layer, filename = file.path(
      outputdir, "output", "mostprobable", 
      paste0(stub, "shannon.tif")), overwrite = TRUE)
  
  ordered.probs <- terra::subset(ordered.probs, 1:nprob) %>% 
    terra::writeRaster(filename = file.path(
      outputdir, "output", "mostprobable", 
      paste0(names(.), ".tif")), 
      overwrite = TRUE)
  
  suppressWarnings(terra::tmpFiles(old = TRUE, remove = TRUE))
  output$timing$finish <- base::date()
  return(output)
}
