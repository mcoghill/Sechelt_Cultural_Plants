summarise <-
function (realisations, lookup, n.realisations = raster::nlayers(realisations), 
    nprob = 3, cpus = 1, outputdir = getwd(), stub = NULL, type = "raw") 
{
    source("./_functions/dsmart_custom/order_stack_values.R")
    
    output <- base::list()
    output$timing <- base::list(start = base::date())
    messages <- c("Attention is required with the following arguments:\n")
    if (type != "prob") {
        if (!(class(realisations) == "RasterStack")) {
            messages <- append(messages, "'realisations': Not a valid RasterStack.\n")
        }
    } else {
        if (is.list(realisations) == FALSE) {
            messages <- append(messages, "'realisations' must be a list of RasterBrick objects when probabilistic predictions are used.'.\n")
        } else {
            if (sum(unlist(lapply(realisations, function(x) class(x) != 
                "RasterBrick"))) > 0) {
                messages <- append(messages, "'realisations' must be a list of RasterBrick objects when probabilistic predictions are used.'.\n")
            }
        }
    }
    if (!(class(lookup) == "data.frame")) {
        messages <- append(messages, "'lookup': Not a valid data.frame.\n")
    }
    if (n.realisations <= 0) {
        messages <- append(messages, "'n.realisations': Value must be greater than 0.\n")
    }
    if (nprob <= 0) {
        messages <- append(messages, "'nprob': Value must be greater than 0.\n")
    }
    if (cpus <= 0) {
        messages <- append(messages, "'cpus': Value must be greater than 0.\n")
    }
    if (!(file.exists(outputdir))) {
        messages <- append(messages, "'outputdir': Output directory does not exist.")
    }
    if (length(messages) > 1) {
        stop(messages)
    }
    if (is.null(stub)) {
        stub <- ""
    } else if (stub == "") {
        stub <- ""
    } else if (!(substr(stub, nchar(stub), nchar(stub)) == "_")) {
        stub <- paste0(stub, "_")
    }
    
    output$call <- base::match.call()
    output$parameters <- base::list(n.realisations = n.realisations, 
        nprob = nprob, cpus = cpus, stub = stub, type = type)
    outputdir <- file.path(outputdir)
    dir.create(file.path(outputdir, "output"), showWarnings = FALSE)
    dir.create(file.path(outputdir, "output", "probabilities"), 
        showWarnings = FALSE)
    dir.create(file.path(outputdir, "output", "mostprobable"), 
        showWarnings = FALSE)
    output$locations <- base::list(root = file.path(outputdir, 
        "output"), probabilities = file.path(outputdir, "output", 
        "probabilities"), mostprobable = file.path(outputdir, 
        "output", "mostprobable"))
    names(lookup) <- c("name", "code")
    param <- nrow(lookup)
    assign("param", param, envir = .GlobalEnv)
    if (type != "prob") {
        
        counts <- app(rast(realisations), tabulate, nbins = param)
        names(counts) <- lookup$name
        probs <- counts / n.realisations
        names(probs) <- lookup$name
        
        # terra::modal(realisations, ties = "first", na.rm = TRUE)
        
        # raster::beginCluster(cpus)
        # counts <- raster::clusterR(realisations, calc, args = list(fun = function(x) {
        #     if (is.na(sum(x))) {
        #         rep(NA, param)
        #     } else {
        #         tabulate(x, nbins = param)
        #     }
        # }), export = "param")
        # raster::endCluster()
        # assign("n.realisations", n.realisations, envir = .GlobalEnv)
        # raster::beginCluster(cpus)
        # probs <- raster::clusterR(counts, calc, args = list(fun = function(x) {
        #     x/n.realisations
        # }), export = "n.realisations")
        # raster::endCluster()
    } else {
        if (length(realisations) == 1 | n.realisations == 1) {
            probs <- realisations[[1]]
        } else {
            # Code for probability estimation not edited in this loop
            raster::beginCluster(cpus)
            probs <- list()
            for (i in 1:param) {
                rlist <- list()
                for (j in 1:n.realisations) {
                  rlist[[j]] <- realisations[[j]][[i]]
                }
                rlist <- stack(rlist)
                probs[[i]] <- raster::clusterR(rlist, calc, args = list(fun = mean))
            }
            raster::endCluster()
            probs <- stack(probs)
        }
    }
    for (i in 1:terra::nlyr(probs)) {
        terra::writeRaster((probs[[i]]), filename = file.path(outputdir, 
            "output", "probabilities", paste0(stub, "prob_", 
                lookup$name[which(lookup$code == i)], ".tif")), 
            overwrite = TRUE)
    }
    if(nprob > 1) {
        if (type != "prob") {
            ordered.indices <- app(counts, order, decreasing = TRUE, na.last = TRUE)[[1:nlyr(counts)]]
        } else {
            ordered.indices <- app(probs, order, decreasing = TRUE, na.last = TRUE)[[1:nlyr(counts)]]
        }
    } else {
        ordered.indices <- terra::modal(realisations, ties = "first", na.rm = TRUE)
    }
    ordered.probs <- app(probs, sort, decreasing = TRUE, na.last = TRUE)[[1:max(2, nprob)]]
    
    # raster::beginCluster(cpus)
    # ordered.probs = raster::clusterR(probs, calc, args = list(fun = function(x) {
    #     if (is.na(sum(x))) {
    #         rep(NA, max(2, nprob))
    #     } else {
    #         sort(x, decreasing = TRUE, na.last = TRUE)[1:max(2, 
    #             nprob)]
    #     }
    # }))
    # raster::endCluster()
    for (i in 1:nprob) {
        terra::writeRaster(ordered.indices[[i]], filename = file.path(outputdir, 
            "output", "mostprobable", paste0(stub, "mostprob_", 
                formatC(i, width = nchar(nrow(lookup)), format = "d", 
                  flag = "0"), "_class.tif")),
            overwrite = TRUE)
        terra::writeRaster(ordered.probs[[i]], filename = file.path(outputdir, 
            "output", "mostprobable", paste0(stub, "mostprob_", 
                formatC(i, width = nchar(nrow(lookup)), format = "d", 
                  flag = "0"), "_probs.tif")),
            overwrite = TRUE)
    }
    confusion <- app(ordered.probs, function(x) {
        (1 - (x[[1]] - x[[2]]))
    }, filename = file.path(outputdir, "output", "mostprobable", 
       paste0(stub, "confusion.tif")), overwrite = TRUE)
    shannon <- app(ordered.probs, function(x) {
        x %>% magrittr::multiply_by(log(x, base = length(x))) %>% 
            sum(na.rm = TRUE) %>% magrittr::multiply_by(-1)
    }, filename = file.path(outputdir, "output", "mostprobable", 
       paste0(stub, "shannon.tif")), overwrite = TRUE)
    
    # raster::beginCluster(cpus)
    # confusion <- raster::clusterR(ordered.probs, fun = function(x) {
    #     (1 - (x[[1]] - x[[2]]))
    # }, filename = file.path(outputdir, "output", "mostprobable", 
    #     paste0(stub, "confusion.tif")), format = "GTiff", overwrite = TRUE, 
    #     NAflag = -9999)
    # shannon <- raster::clusterR(ordered.probs, fun = function(x) {
    #     x %>% magrittr::multiply_by(log(x, base = length(x))) %>% 
    #         sum(na.rm = TRUE) %>% magrittr::multiply_by(-1)
    # }, filename = file.path(outputdir, "output", "mostprobable", 
    #     paste0(stub, "shannon.tif")), format = "GTiff", overwrite = TRUE, 
    #     NAflag = -9999)
    # raster::endCluster()
    output$timing$finish <- base::date()
    return(output)
}
