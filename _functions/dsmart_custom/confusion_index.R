confusion_index <-
function (r, cpus) 
{
    tuning <- .blocks_per_node(raster::nrow(r), raster::ncol(r), 
        cpus = cpus)
    sorted <- sort_stack_values(r, cpus)
    raster::beginCluster(cpus)
    output <- raster::clusterR(sorted, fun = function(x) {
        (1 - (x[[1]] - x[[2]]))
    }, filename = tempfile(fileext = ".tif"), format = "GTiff", 
        overwrite = TRUE, NAflag = -9999, datatype = "FLT4S", 
        m = tuning)
    raster::endCluster()
    return(output)
}
