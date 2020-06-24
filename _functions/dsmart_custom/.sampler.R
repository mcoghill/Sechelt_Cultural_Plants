## NOTE: THIS FUNCTION IS NO LONGER NEEDED

# .sampler <-
# function (covariates, polygons, composition, poly.id, n.realisations = 100, 
#     rate = 15, method.sample = "by_polygon", method.allocate = "weighted") 
# {
#     source("./_functions/.allocate.R")
#     
#     poly <- base::subset(polygons, as.data.frame(polygons)[, 1] == poly.id)
#     n.samples <- 0
#     if (method.sample == "by_area") {
#         area <- terra::area(poly)/1000000
#         if (area < 1) {
#             area <- 1
#         }
#         n.samples <- rate * base::trunc(area)
#     } else if (method.sample == "by_polygon") {
#         n.samples <- rate
#     } else stop("Sampling method \"", method.sample, "\" is unknown")
#     
#     # Need to add cell number to ouput. Easiest way is to create a raster of
#     # cell values
#     cells <- covariates[[1]] %>% magrittr::set_names("cells")
#     values(cells) <- 1:ncell(cells)
#     covariates <- c(covariates, cells)
#     
#     poly.samples <- as.data.frame(terra::extract(covariates, poly)) %>% 
#         dplyr::rename_at(vars(names(.)), ~names(covariates)) %>% 
#         dplyr::filter(complete.cases(.)) %>% 
#         dplyr::sample_n(size = (n.samples * n.realisations), replace = TRUE)
#     
#     soil_class <- character()
#     if (method.allocate == "weighted") {
#         poly.classes <- base::as.character(composition[base::which(composition[, 
#             1] == poly.id), 3])
#         poly.weights <- composition[base::which(composition[, 
#             1] == poly.id), 4]
#         if (length(poly.classes) == 0) {
#             stop(paste0("No map unit composition for polygon ", 
#                 poly.id))
#         } else {
#             soil_class <- .allocate(poly.classes, n = n.samples * 
#                 n.realisations, method = "weighted", weights = poly.weights)
#         }
#     } else if (method.allocate == "random-mapunit") {
#         poly.classes <- as.character(composition[which(composition[, 
#             1] == poly.id), 3])
#         if (length(poly.classes) == 0) {
#             stop(paste0("No map unit composition for polygon ", 
#                 poly.id))
#         } else {
#             soil_class <- .allocate(poly.classes, n = n.samples * 
#                 n.realisations, method = "random", weights = NULL)
#         }
#     } else if (method.allocate == "random-all") {
#         poly.classes <- as.character(unique(composition[, 3]))
#         if (length(poly.classes) == 0) {
#             stop("No soil classes available to allocate to")
#         } else {
#             soil_class <- .allocate(poly.classes, n = n.samples * 
#                 n.realisations, method = "random", weights = NULL)
#         }
#     } else stop("Allocation method is unknown")
#     
#     xy <- as.data.frame(terra::xyFromCell(covariates, poly.samples$cell))
#     meta <- list(realisation = base::rep(1:n.realisations, times = n.samples), 
#         type = base::rep("virtual", nrow(xy)), sampling = base::rep(method.sample, 
#             base::nrow(xy)), allocation = base::rep(method.allocate, 
#             base::nrow(xy)))
#     poly.samples <- cbind(as.data.frame(meta), xy, soil_class, 
#         poly.samples[, 2:base::ncol(poly.samples)])
#     return(poly.samples)
# }
