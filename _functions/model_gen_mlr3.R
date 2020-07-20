# Please note that I have tried many different things to try to get probability estimation to work
# using the ranger package in mlr and mlr3 however it never seems to work when there are
# few data points for a class, hence why predict_type is response. Other learners have no
# issue using the prob predict_type

# Essentially the way I have this working is as follows: First, create the pipe
# operators which allows for subsetting certain percentages of the top variables
# (ex: a subset of the top 10%, 20%, etc. of variables) and then creates single
# models from the subsetted features. The model performance is evaluated for each
# and then the best model is fitted. This is the "tuning" step and it is performed
# on each of the input "traindat" dataframes. The other option for
# tuning is recursive feature elimination, however this would take much longer, 
# especially if there are a lot fo variables in a dataset. 
# After the best feature subset is found, 10 fold cross validation repeated 5 times 
# is performed in order to generate the model metrics

model_gen_mlr3 <- function(
  traindat, 
  target, 
  feature_selection = "filter", 
  type = "response", 
  folds = 10, repeats = 5) {
  
  library(mlr3verse)
  #traindat <- pres_abs
  #target <- "Pres"
  #feature_selection = "filter"
  
  # perform input checks
  if(!is.list(traindat)) 
    stop("Please provide a named list of dataframes for the 'traindat' argument.")
  
  if(!all(sapply(traindat, function(x) is.data.frame(x))))
    stop("List argument provided, but not all list objects are dataframes.")
  
  if(!is.character(target))
    stop("'target' argument is not a character vector. Additionally, 'target' should
         be of length 1 and should be the same for each listed dataframe.")
  
  if(feature_selection == "rfe") {
    if(!("mlr3fselect" %in% installed.packages()[, "Package"]))
      devtools::install_github("mlr-org/mlr3fselect", force = TRUE)
    
    library(mlr3fselect)
    
  } else if(feature_selection != "filter") {
    stop("Feature selection variable is not properly set. Choose one of 'filter' or 'rfe'")
  }
  
  # Check target data column for determining regression vs classification
  if(all(sapply(traindat, function(x) 
    class(x[, target])) == "factor")) {
    mod_type <- "classif.ranger"
  } else {
    mod_type <- "regr.ranger"
  }
  
  if(all(sapply(traindat, function(x) 
    nlevels(x[, target])) <= 2)) {
    if(any(c(TRUE, FALSE) %in% levels(traindat[[1]][, target]))) {
      positive <- "TRUE"
    } else {
      positive <- NULL
    }
  }
  
  measures <- if(mod_type == "regr.ranger") {
    "regr.mse"
  } else {
    "classif.ce"
  }
  
  # Check for X and Y columns. If present, a spatial cv will be performed
  if(all(sapply(traindat, function(x) 
    length(which(names(x) %in% c("X", "x"))) == 1 && length(which(names(x) %in% c("Y", "y"))) == 1))) {
    spatial <- TRUE
    coordinate_cols <- lapply(traindat, function(x)
      names(x)[names(x) %in% c("X", "x", "Y", "y")])
    
    if(!("mlr3spatiotempcv" %in% installed.packages()[, "Package"]))
      devtools::install_github("mlr-org/mlr3spatiotempcv", force = TRUE)
    
    library(mlr3spatiotempcv)
    message("Performing spatial sampling")
    
    # Fixes issue descrived above
    mlr_reflections$task_types <- mlr_reflections$task_types[package == "mlr3spatiotempcv", ]
    
    resampling_outer <- rsmp("repeated-spcv-coords", folds = folds, repeats = repeats)
    
    if(mod_type == "regr.ranger") {
      tasks <- lapply(seq_along(traindat), function(x) {
        TaskRegrST$new(id = gsub(".gpkg$", "", basename(names(traindat[x]))), 
                       backend = traindat[[x]], 
                       target = target, 
                       coords_as_features = FALSE, 
                       crs = st_crs(3005)$proj4string, 
                       coordinate_names = coordinate_cols[[x]])
      })
    } else {
      tasks <- lapply(seq_along(traindat), function(x) {
        TaskClassifST$new(id = gsub(".gpkg$", "", basename(names(traindat[x]))), 
                          backend = traindat[[x]], 
                          target = target, 
                          coords_as_features = FALSE, 
                          crs = st_crs(3005)$proj4string, 
                          positive = positive,
                          coordinate_names = coordinate_cols[[x]])
      })
    }
    
  } else {
    # remove any x and y columns if it isn't spatial
    spatial <- FALSE
    
    message("Performing non-spatial sampling")
    
    sp_index <- which(sapply(traindat, function(x) 
      length(which(names(x) %in% c("X", "Y"))) == 2 || length(which(names(x) %in% c("x", "y")))))
    
    traindat <- lapply(traindat, function(x) x[, names(x)[!names(x) %in% c("X", "Y", "x", "y")]])
    
    resampling_outer <- rsmp("repeated_cv", folds = folds, repeats = repeats)
    
    if(mod_type == "regr.ranger") {
      tasks <- lapply(seq_along(traindat), function(x) {
        TaskRegr$new(id = gsub(".gpkg$", "", basename(names(traindat[x]))), 
                     backend = traindat[[x]], 
                     target = target)
      })
    } else {
      tasks <- lapply(seq_along(traindat), function(x) {
        TaskClassif$new(id = gsub(".gpkg$", "", basename(names(traindat[x]))), 
                        backend = traindat[[x]], 
                        positive = positive,
                        target = target)
      })
    }
  }
  names(tasks) <- names(traindat)
  
  # Create pipeop learner, note I had some issues with this, error included duplicated value
  # in "name" column, but if the mlr3spatiaotempcv package is not loaded here, this 
  # works fine. I submitted this as an issue to the mlr3pipelines page: 
  # https://github.com/mlr-org/mlr3pipelines/issues/368
  # The issue is still ongoing as of June 26...I appear to have opened a can of
  # worms with this!
  
  # Note: don't use the mtry argument here, or problems will arrive later on
  lrn <- lrn(
    mod_type, 
    num.threads = parallel::detectCores(),
    importance = "impurity",
    predict_type = type
  )
  
  if(feature_selection == "filter") {
    po_lrn <- po("learner", lrn)
    
    # Create feature filter based on variable importance
    filter <- po("filter", filter = mlr3filters::flt("importance", learner = lrn))
    
    # Create process (new learner) for filtering the task
    glrn <- GraphLearner$new(filter %>>% po_lrn)
    glrn$predict_type <- type
    
    # Create filter parameters (i.e.: filter the top features of importance using
    # a given percentage of total features)
    param_set <- ParamSet$new(
      params = list(ParamDbl$new("importance.filter.frac", 
                                 lower = min(
                                   10 / length(tasks[[1]]$feature_names), 
                                   0.1), 
                                 upper = 1))
    )
    
    # Create the autotuner which will generate models using the top features of 
    # importance, using 2 features and numbers between 2 and the total number of 
    # features (ex: if tuner has resolution 10, there will be 9 other values of 
    # features used in model creation. The best model will be used to create the 
    # resampling task.)
    resolution <- round(max(sapply(tasks, function(x) length(x$feature_names))) / 10)
    feature_sel <- AutoTuner$new(
      learner = glrn, 
      resampling = resampling_outer, 
      measures = msr(measures),
      tune_ps = param_set, 
      terminator = term("none"), 
      tuner = tnr("grid_search", resolution = resolution)
    )
    feature_sel$store_tuning_instance <- TRUE
    
    design <- benchmark_grid(tasks, feature_sel, resamplings = rsmp("holdout", ratio = 1))
    bmr <- benchmark(design, store_models = TRUE)
    all_results <- lapply(1:bmr$n_resample_results, function(x)
      bmr$resample_result(x)$learners[[1]])
    names(all_results) <- names(traindat)
    
  } else if(feature_selection == "rfe") {
    f_select <- fs("rfe", min_features = 2, recursive = TRUE) # Recalculates importance each iteration
    resampling_outer <- rsmp("repeated_cv", folds = folds, repeats = repeats) # Fixes problem for resampling
    
    # Create the feature selection learner. Models will be compared at each iteration
    # of a feature being removed. If the models don't improve after 5 iterations, then
    # feature selection stops, using the best evaluated model
    
    # PROBLEM WITH INTEGRATING THIS INTO A BENCHMARK: THERE IS NO PROPER WAY
    # TO STORE THE MODELS, SO I HAVE TO GET CREATIVE AND DO IT MYSELF
    # feature_sel <- AutoFSelect$new(
    #   learner = lrn,
    #   resampling = resampling_outer, 
    #   measure = msr(measures), 
    #   terminator = term("stagnation", iters = 5, threshold = 1e-5), 
    #   fselect = f_select
    # )
    # feature_sel$store_fselect_instance <- TRUE
    
    all_results <- lapply(tasks, function(x) {
      feature_sel <- FSelectInstance$new(
        task = x,
        learner = lrn,
        resampling = resampling_outer,
        measure = msr(measures),
        terminator = term("stagnation", iters = round(0.1 * x$ncol), threshold = 1e-5),
        store_models = TRUE
      )
      f_select$optimize(feature_sel)
      
      find_feats <- lapply(feature_sel$archive$data()$x_domain, unlist)
      feature_sel_id <- which(unlist(lapply(find_feats, function(y)
        all(y == unlist(feature_sel$result$x_domain[[1]])))))
      
      model_set <- feature_sel$archive$data()[feature_sel_id, ]$resample_result[[1]]
      best_learner <- model_set$learners[[which.min(lapply(model_set$learners, function(y) y$oob_error()))]]
      return(best_learner)
    })
    names(all_results) <- names(tasks)
  } else if (feature_selection == "none") {
    
    rr <- resample(tasks[[1]], lrn, resampling_outer, store_models = TRUE)
    bmr <- as_benchmark_result(rr)
    
  }
  return(list(learners = all_results, tasks = tasks))
}
