order_stack_values <-
function (r, cpus, n = nlayers(r)) 
{
    tuning <- .blocks_per_node(raster::nrow(r), raster::ncol(r), 
        cpus = cpus)
    raster::beginCluster(cpus)
    output = raster::clusterR(r, calc, args = list(fun = function(x) {
        if (is.na(sum(x))) {
            rep(NA, n)
        } else {
            order(x, decreasing = TRUE, na.last = TRUE)[1:n]
        }
    }), m = tuning)
    raster::endCluster()
    return(output)
}
