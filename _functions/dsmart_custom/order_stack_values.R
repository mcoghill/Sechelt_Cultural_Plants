# order_stack_values <-
# function (r, cpus, n = nlyr(r)) 
# {
#     # tuning <- .blocks_per_node(terra::nrow(r), terra::ncol(r), 
#     #     cpus = cpus)
#     output <- order(r, decreasing = TRUE, na.last = TRUE)[1:n]
#     
#     # raster::beginCluster(cpus)
#     # output = raster::clusterR(r, calc, args = list(fun = function(x) {
#     #     if (is.na(sum(x))) {
#     #         rep(NA, n)
#     #     } else {
#     #         order(x, decreasing = TRUE, na.last = TRUE)[1:n]
#     #     }
#     # }), m = tuning)
#     # raster::endCluster()
#     return(output)
# }
