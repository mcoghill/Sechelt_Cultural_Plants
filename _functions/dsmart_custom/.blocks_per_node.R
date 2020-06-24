.blocks_per_node <-
function (rows, cols, cpus, max_cells = 1000000) 
{
    m <- NA
    if (max_cells < (rows * cols)) {
        block_rows <- max_cells/cols
        total_blocks <- rows/block_rows
        m <- ceiling(total_blocks/cpus)
    }
    else {
        m <- 2
    }
    return(m)
}
