# sort_stack_values <-
# function (r, cpus, n = nlayers(r), decreasing = TRUE) 
# {
#     tuning <- .blocks_per_node(raster::nrow(r), raster::ncol(r), 
#         cpus = cpus)
#     raster::beginCluster(cpus)
#     output = raster::clusterR(r, calc, args = list(fun = function(x) {
#         if (is.na(sum(x))) {
#             rep(NA, n)
#         } else {
#             sort(x, decreasing = decreasing, na.last = TRUE)[1:n]
#         }
#     }), m = tuning)
#     raster::endCluster()
#     return(output)
# }
