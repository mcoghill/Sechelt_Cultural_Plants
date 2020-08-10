dsmart <- function(
  covariates, polygons, composition, rate = 15, reals = 100, 
  observations = NULL, method.sample = "by_polygon", method.allocate = "weighted", 
  method.model = NULL, args.model = NULL, strata = NULL, nprob = 3, 
  outputdir = getwd(), stub = NULL, factors = NULL, type = "raw") {
  
  source("./_functions/dsmart_custom/disaggregate.R")
  source("./_functions/dsmart_custom/summarise.R")
  
  output <- base::list()
  output$timing <- base::list(start = base::date())
  if(is.null(stub)) {
    stub <- ""
  } else if(stub == "") {
    stub <- ""
  } else if(!(substr(stub, nchar(stub), nchar(stub)) == "_")) {
    stub <- paste0(stub, "_")
  }
  outputdir <- file.path(outputdir)
  dir.create(file.path(outputdir, "output"), showWarnings = FALSE)
  base::match.call() %>% base::deparse() %>% 
    base::write(file = file.path(outputdir, "output", "dsmart_function_call.txt"))
  output$disaggregate <- disaggregate(
    covariates, polygons, 
    composition, rate = rate, reals = reals, 
    observations = observations, method.sample = method.sample, 
    method.allocate = method.allocate, method.model = method.model, 
    args.model = args.model, strata = strata, outputdir = outputdir, 
    stub = stub, factors = factors, type = type)
  
  if(type != "prob") {
    realisations <- terra::rast(output$disaggregate$locations$realisations)
  } else {
    realisations <- lapply(output$disaggregate$locations$realisations, terra::rast)
  }
  lookup <- read.table(output$disaggregate$locations$lookup, header = TRUE, sep = ",")
  output$summarise <- summarise(
    realisations, lookup, n.realisations = reals, 
    nprob = nprob, outputdir = outputdir, stub = stub, 
    type = type)
  output$timing$finish <- base::date()
  return(output)
}
