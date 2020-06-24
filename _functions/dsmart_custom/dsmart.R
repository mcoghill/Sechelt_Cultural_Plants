dsmart <-
function (covariates, polygons, composition, rate = 15, reals = 100, 
    observations = NULL, method.sample = "by_polygon", method.allocate = "weighted", 
    method.model = NULL, args.model = NULL, strata = NULL, nprob = 3, 
    outputdir = getwd(), stub = NULL, cpus = 1, factors = NULL, 
    type = "raw") 
{
    output <- base::list()
    output$timing <- base::list(start = base::date())
    if (is.null(stub)) {
        stub <- ""
    }
    else if (stub == "") {
        stub <- ""
    }
    else if (!(substr(stub, nchar(stub), nchar(stub)) == "_")) {
        stub <- paste0(stub, "_")
    }
    outputdir <- file.path(outputdir)
    dir.create(file.path(outputdir, "output"), showWarnings = FALSE)
    base::match.call() %>% base::deparse() %>% base::write(file = file.path(outputdir, 
        "output", "dsmart_function_call.txt"))
    output$disaggregate <- disaggregate(covariates, polygons, 
        composition, rate = rate, reals = reals, cpus = cpus, 
        observations = observations, method.sample = method.sample, 
        method.allocate = method.allocate, method.model = method.model, 
        args.model = args.model, strata = strata, outputdir = outputdir, 
        stub = stub, factors = factors, type = type)
    if (type != "prob") {
        realisations <- raster::stack()
        for (filename in base::list.files(path = file.path(outputdir, 
            "output", "realisations"), pattern = ".tif$", full.names = TRUE)) {
            r <- raster::raster(filename)
            realisations <- raster::stack(realisations, r)
        }
    }
    else {
        realisations <- list()
        i <- 1
        for (filename in base::list.files(path = file.path(outputdir, 
            "output", "realisations"), pattern = ".tif$", full.names = TRUE)) {
            realisations[[i]] <- raster::brick(filename)
            i <- i + 1
        }
    }
    lookup <- read.table(file.path(outputdir, "output", paste0(stub, 
        "lookup.txt")), header = TRUE, sep = ",")
    output$summarise <- summarise(realisations, lookup, n.realisations = reals, 
        nprob = nprob, cpus = cpus, outputdir = outputdir, stub = stub, 
        type = type)
    output$timing$finish <- base::date()
    return(output)
}
