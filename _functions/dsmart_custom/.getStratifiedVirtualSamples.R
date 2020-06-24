.getStratifiedVirtualSamples <-
function (covariates, polygons, composition, strata, n.realisations = 100, 
    rate = 15, method.sample = "by_polygon", method.allocate = "weighted") 
{
    source("./_functions/dsmart_custom/.allocate.R")
    if (ncol(composition) == 4) {
        names(composition) <- c("poly", "mapunit", "soil_class", 
            "proportion")
    } else if (ncol(composition) == 5) {
        names(composition) <- c("poly", "mapunit", "stratum", 
            "soil_class", "proportion")
    } else stop("Map unit composition in unknown format.")
    
    cell <- covariates[[1]] %>% magrittr::set_names("cell")
    values(cell) <- 1:ncell(cell)
    covariates <- c(cell, covariates)
    
    samples <- foreach::foreach(x = 1:length(polygons)) %do% {
        poly.id = as.data.frame(polygons)[, 1][x]
        cat(paste0(
            "\nGenerating stratified samples for polygon ", 
            poly.id, " [", x, "/", length(polygons), "]"))
        poly <- subset(polygons, as.data.frame(polygons)[, 1] == poly.id)
        n.samples <- 0
        if (method.sample == "by_area") {
            area <- terra::area(poly)/1000000
            if (area < 1) {
                area <- 1
            }
            n.samples <- rate * base::trunc(area)
        } else if (method.sample == "by_polygon") {
            n.samples <- rate
        }  else stop("Sampling method \"", method.sample, "\" is unknown")
        
        poly.samples <- as.data.frame(terra::extract(covariates, poly)) %>% 
            dplyr::rename_at(vars(names(.)), ~names(covariates)) %>% 
            dplyr::filter(complete.cases(.)) %>%
            {if(nrow(.) > 0) {
                dplyr::sample_n(., size = (n.samples * n.realisations), replace = TRUE)
            } else . }
        xy <- as.data.frame(terra::xyFromCell(covariates, poly.samples$cell))
        poly.samples.strata <- terra::extract(strata, xy)
        poly.samples <- cbind(poly.samples.strata, poly.samples[, 
            2:ncol(poly.samples)])
        names(poly.samples)[names(poly.samples) == colnames(poly.samples.strata)] <- "stratum"
        poly.samples <- poly.samples[order(poly.samples$stratum), ]
        soil_class <- character()
        for (stratum in unique(poly.samples[, 1])) {
            stratum.n <- length(which(poly.samples[, 1] == stratum))
            stratum.classes <- as.character(composition[which((composition$poly == 
                poly.id) & (composition$stratum == stratum)), 
                4])
            if (length(stratum.classes) == 0) {
                stop(paste0("No classes available in stratum ", 
                  stratum, " of polygon ", poly.id))
            } else if (method.allocate == "weighted") {
                if (!(ncol(composition) == 5)) {
                  stop("Weighted-random allocation specified but no stratum weights available.")
                } else {
                  stratum.weights <- composition[which((composition$poly == 
                    poly.id) & (composition$stratum == stratum)), 
                    5]
                  alloc <- .allocate(stratum.classes, n = stratum.n, 
                    method = "weighted", weights = stratum.weights)
                  soil_class <- append(soil_class, alloc)
                }
            }
            else if (method.allocate == "random") {
                alloc <- .allocate(stratum.classes, n = stratum.n, 
                  method = "random")
                soil_class <- append(soil_class, alloc)
            }
        }
        poly.samples <- cbind(xy, soil_class, poly.samples[, 
            1:ncol(poly.samples)])
        poly.samples <- poly.samples[base::sample(nrow(poly.samples)), ]
        meta <- list(realisation = rep(1:n.realisations, times = n.samples)[0:nrow(poly.samples)], 
            type = base::rep("virtual", nrow(xy)), sampling = base::rep(method.sample, 
                nrow(xy)), allocation = base::rep(method.allocate, 
                nrow(xy)))
        poly.samples <- cbind(as.data.frame(meta), poly.samples)
        return(poly.samples)
    }
    samples <- data.table::rbindlist(samples)
    return(samples)
}
