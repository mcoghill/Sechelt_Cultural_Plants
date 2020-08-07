disaggregate <-
function (covariates, polygons, composition, rate = 15, reals = 100, 
    observations = NULL, method.sample = "by_polygon", method.allocate = "weighted", 
    method.model = NULL, args.model = NULL, strata = NULL, outputdir = getwd(), 
    stub = NULL, factors = NULL, type = "raw", predict = TRUE) 
{
    source("./_functions/dsmart_custom/.getVirtualSamples.R")
    source("./_functions/dsmart_custom/.getStratifiedVirtualSamples.R")
    source("./_functions/dsmart_custom/.observations.R")
    source("./_functions/predict_landscape.R")
    output <- base::list()
    output$timing <- base::list(start = base::date())
    messages <- c("Attention is required with the following arguments:\n")
    if (!(class(covariates) == "SpatRaster")) {
        messages <- append(messages, "'covariates': Not a valid SpatRaster.\n")
    }
    if (!(class(polygons) == "SpatVector")) {
        messages <- append(messages, "'polygons': Not a valid SpatVector\n")
    }
    if (!(class(composition) == "data.frame")) {
        messages <- append(messages, "'composition': Not a valid data.frame.\n")
    }
    if (rate <= 0) {
        messages <- append(messages, "'n': Value must be greater than 0.\n")
    }
    if (reals <= 0) {
        messages <- append(messages, "'reals': Value be greater than 0.\n")
    }
    if (!(is.null(observations))) {
        if (!(class(observations) == "data.frame")) {
            messages <- append(messages, "'observations': Not a valid data.frame.\n")
        }
    }
    if (!(file.exists(outputdir))) {
        messages <- append(messages, "'outputdir': Output directory does not exist.\n")
    }
    if (!(is.null(strata))) {
        if (!(class(strata) == "SpatRaster")) {
            messages <- append(messages, "'strata': Not a valid SpatRaster\n")
        }
    }
    if (!(is.null(method.model))) {
        if (!(is.character(method.model) & length(method.model) ==
            1)) {
            messages <- append(messages, "'method.model' must be NULL or a single character value.\n")
        }
    }
    if (length(messages) > 1) {
        stop(messages)
    }
    if (is.character(method.model) & length(method.model) ==
        1) {
        require(caret)
    }
    if (is.null(stub)) {
        stub <- ""
    }  else if (stub == "") {
        stub <- ""
    }  else if (!(substr(stub, nchar(stub), nchar(stub)) == "_")) {
        stub <- paste0(stub, "_")
    }
    # Repeated cross validation generates the best model anyway, no need to
    # repeat the model for multiple realisations if this is the case
    if(!is.null(method.model) && 
       args.model$trControl$method == "repeatedcv" &&
       is.numeric(args.model$trControl$repeats)) {
        rate <- rate * reals
        reals <- 1
    }
    # Override probability method in the trainControl function
    if(!is.null(method.model) && type == "prob") {
        args.model$trControl$classProbs <- TRUE
    } else {
        args.model$trControl$classProbs <- FALSE
    }
    
    output$call <- base::match.call()
    output$parameters <- base::list(rate = rate, reals = reals, 
        method.sample = method.sample, method.allocate = method.allocate, 
        method.model = method.model, args.model = args.model, 
        stub = stub, factors = factors, type = type)
    outputdir <- file.path(outputdir)
    dir.create(file.path(outputdir, "output"), showWarnings = FALSE)
    dir.create(file.path(outputdir, "output", "realisations"), 
        showWarnings = FALSE)
    dir.create(file.path(outputdir, "output", "models"), showWarnings = FALSE)
    dir.create(file.path(outputdir, "output", "probabilities"), 
               showWarnings = FALSE)
    output$locations <- list(root = file.path(outputdir, "output"), 
        realisations = file.path(outputdir, "output", "realisations"), 
        models = file.path(outputdir, "output", "models"))
    if (!(is.null(strata))) {
        names(composition) <- c("poly", "mapunit", "stratum", 
            "soil_class", "proportion")
    } else {
        names(composition) <- c("poly", "mapunit", "soil_class", 
            "proportion")
    }
    
    # Deal with polygons
    polygons <- polygons[, 1] %>% magrittr::set_names("poly")
    polys_to_remove <- composition %>% stats::complete.cases() %>% 
        magrittr::equals(FALSE) %>% base::which() %>% composition$poly[.] %>% 
        base::unique()
    if (length(polys_to_remove) > 0) {
        polygons <- polygons[!polygons$poly %in% polys_to_remove]
        composition <- composition %>% dplyr::filter(!(poly %in% 
            polys_to_remove))
        warning(base::paste0("The following polygons were removed from further analysis because they have incomplete or undefined map unit compositions: ", 
            base::paste(base::as.character(polys_to_remove), 
                collapse = ", ")))
    }
    
    write.table(names(covariates), file.path(outputdir, "output", 
        paste0(stub, "covariate_names.txt")), quote = FALSE, 
        sep = ",", row.names = FALSE, col.names = FALSE)
    message(Sys.time(), " Generating samples for ", reals, " realisations")
    
    if (is.null(strata)) {
        samples <- .getVirtualSamples(covariates, polygons, composition, 
            n.realisations = reals, rate = rate, method.sample = method.sample, 
            method.allocate = method.allocate)
    } else {
        samples <- .getStratifiedVirtualSamples(covariates, polygons, 
            composition, strata = rast(strata), n.realisations = reals, rate = rate, 
            method.sample = method.sample, method.allocate = method.allocate)
    }
    samples <- samples %>% complete.cases() %>% dplyr::filter(samples, .)
    if (!(is.null(observations))) {
        names(observations) <- c("x", "y", "class")
        observations <- .observations(observations, covariates)
        write.table(observations, file.path(outputdir, "output", 
            paste0(stub, "observations_with_covariates.txt")), 
            sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)
    }
    write.table(samples, file.path(outputdir, "output", paste0(stub, 
        "virtual_samples.txt")), sep = ",", quote = FALSE, col.names = TRUE, 
        row.names = FALSE)
    levs <- as.character(unique(composition$soil_class))
    if (!(is.null(observations))) {
        levs <- base::union(levs, as.character(unique(observations$soil_class)))
    }
    levs <- sort(levs)
    lookup <- data.frame(name = levs, code = 1:length(levs), stringsAsFactors = FALSE)
    write.table(lookup, file.path(outputdir, "output", paste0(stub, 
        "lookup.txt")), sep = ",", quote = FALSE, col.names = TRUE, 
        row.names = FALSE)
    output$locations$lookup <- file.path(outputdir, "output", 
        paste0(stub, "lookup.txt"))
    
    # Modelling for each realisation
    for (j in 1:reals) {
        message(Sys.time(), " Realisation ", j)
        startcol <- numeric(0)
        if (is.null(strata)) {
            startcol = 8
        } else {
            startcol = 9
        }
        s <- samples[which(samples$realisation == j), startcol:ncol(samples)]
        soil_class <- as.character(samples$soil_class[which(samples$realisation == 
            j)])
        if (!(is.null(observations))) {
            s <- rbind(s, observations[, 8:ncol(observations)])
            soil_class <- append(soil_class, as.character(observations$soil_class))
        }
        soil_class <- factor(soil_class, levels = levs)
        if (!is.null(factors)) {
            s <- s %>% dplyr::mutate(across(factors, as.factor))
        }
        zeroes <- sum(table(soil_class) == 0) > 0
        rclt <- cbind(c(1:sum(table(soil_class) != 0)), c(1:length(table(soil_class)))[table(soil_class) != 
            0])
        if (is.null(method.model)) {
            model <- C50::C5.0(s, y = soil_class)
            out <- utils::capture.output(summary(model))
        } else {
            soil_class <- base::droplevels(soil_class)
            
            # Note: As of July 2020, mlr/mlr3 methods do NOT increase the speed here
            model <- base::do.call(caret::train, c(list(x = s, 
                y = soil_class, method = method.model), args.model))
            out <- utils::capture.output(model$finalModel)
        }
        
        cat(out, file = paste0(outputdir, "/output/models/", 
            stub, "model_", formatC(j, width = nchar(reals), 
                format = "d", flag = "0"), ".txt"), sep = "\n", 
            append = TRUE)
        save(model, file = file.path(outputdir, "output", "models", 
            paste0(stub, "model_", formatC(j, width = nchar(reals), 
                format = "d", flag = "0"), ".RData")))
        
        # Model prediction using tiling method
        if (predict) {
            if (type != "prob") {
                
                r1 <- predict_landscape(model, covariates, tilesize = 500,
                                        outDir = file.path(outputdir, "tiles"), 
                                        type = type)
                
                r1 <- writeRaster(r1, file.path(outputdir, "output",
                            "realisations", paste0(stub, "realisation_",
                            formatC(j, width = nchar(reals), format = "d",
                            flag = "0"), ".tif")), 
                            overwrite = TRUE, wopt = list(datatype = "INT2S"))
                
                if (zeroes == TRUE & is.null(method.model) == FALSE) {
                    r1 <- classify(r1, rcl = rclt, filename = file.path(outputdir, 
                                   "output", "realisations", paste0(stub, 
                                   "realisation_", formatC(j, width = nchar(reals), 
                                    format = "d", flag = "0"), ".tif")), 
                                   overwrite = TRUE, wopt = list(datatype = "INT2S"))
                }
            } else { # else if type == prob
                if (zeroes == FALSE | is.null(method.model)) {
                    r1 <- predict_landscape(model, covariates, tilesize = 500,
                                            outDir = file.path(outputdir, "tiles"), 
                                            type = "prob")
                    
                    # output is file list for each ss call
                    r1 <- rast(grep("pred.tif", r1, value = TRUE, invert = TRUE)) %>% 
                        writeRaster(file.path(outputdir, "output", "realisations",
                                              paste0(stub, "realisation_",
                                                     formatC(j, width = nchar(reals), format = "d", flag = "0"), ".tif")),
                                    overwrite = TRUE)
                  
                } else {
                    tmp1 <- predict_landscape(model, covariates, tilesize = 500,
                                      outDir = file.path(outputdir, "tiles"), 
                                      type = "prob")
                    tmp2 <- (rast(tmp1[1]) * 0) %>% 
                        writeRaster(tempfile(pattern = "spat_", fileext = ".tif"))
                    
                    r1 <- foreach(i = 1:nrow(lookup), .combine = c) %do% {
                        if (i %in% rclt[, 2]) {
                            rast(tmp1[which(rclt[, 2] == i)])
                        } else {
                            tmp2
                        }
                    } %>% magrittr::set_names(lookup$name) %>% 
                        writeRaster(file.path(outputdir, "output", "realisations",
                                              paste0(stub, "realisation_", formatC(j, width = nchar(reals), format = "d", flag = "0"),
                                                     ".tif")), 
                                    overwrite = TRUE)
                  rm(tmp1, tmp2)
                }
            }
        }
    }
    output$timing$finish <- base::date()
    return(output)
}
