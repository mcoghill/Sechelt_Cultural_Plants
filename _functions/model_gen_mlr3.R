# Please note that I have tried many different things to try to get probability 
# estimation to work using the ranger package in mlr and mlr3 however it never 
# seems to work when there are few data points for a class, hence why predict_type 
# defaults to "response". Other learners have no issue using the prob predict_type

# Essentially the way I have this working is as follows: First, create the pipe
# operators which allows for subsetting certain percentages of the top variables
# (ex: a subset of the top 10%, 20%, etc. of variables) and then creates single
# models from the subsetted features. The model performance is evaluated for each
# and then the best model is fitted. This is the "tuning" step and it is performed
# on each of the input "traindat" dataframes. The other option for
# tuning is recursive feature elimination, however this would take much longer, 
# especially if there are a lot of variables in a dataset. 
# 10 fold cross validation repeated 5 times is performed on each of the data 
# subsets in order to generate the final models and metrics

# This function requires a list of sf dataframes to perform spatial analyses
# (note: will need to put stops in place that require that as an input)

model_gen_mlr3 <- function(
  traindat, 
  target, 
  filter_correlation = FALSE,
  corr_cutoff = 0.9,
  feature_selection = "filter", 
  type = "response", 
  folds = 10, repeats = 5) {
  
  require(tidyverse)
  require(mlr3verse)
  require(mlr3spatiotempcv)
  require(sf)
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
  
  if(feature_selection %in% c("rfe", "none")) {
    if(!("mlr3fselect" %in% installed.packages()[, "Package"]))
      install.packages("mlr3fselect")
    
    require(mlr3fselect)
    
  } else if(feature_selection != "filter") {
    stop("Feature selection variable is not properly set. Choose one of 'filter', 'rfe', or 'none'")
  }
  
  # Check target data column for determining regression vs classification
  if(all(sapply(traindat, function(x) 
    class(x[, target][[1]])) == "factor")) {
    mod_type <- "classif.ranger"
  } else {
    mod_type <- "regr.ranger"
  }
  
  # If there are only 2 levels that are either true or false, create a positive true result
  if(all(sapply(traindat, function(x) 
    nlevels(x[, target][[1]])) <= 2)) {
    if(any(c(TRUE, FALSE) %in% levels(traindat[[1]][, target][[1]]))) {
      positive <- "TRUE"
    } else {
      positive <- NULL
    }
  }
  
  measures <- ifelse(mod_type == "regr.ranger", "regr.mse", "classif.ce")
  
  # Check for sf inheritance. If true, a spatial cv will be performed
  if(all(sapply(traindat, inherits, "sf"))) {
    coordinate_cols <- lapply(traindat, function(x) colnames(st_coordinates(x)))
    
    message("Performing spatial sampling")
    
    resampling_outer <- rsmp("repeated_spcv_coords", folds = folds, repeats = repeats)
    
    traindat <- lapply(traindat, function(x) {
      fct <- names(Filter(is.factor, x))
      fct <- fct[!fct %in% c(target, attr(x, "sf_column"))]
      if(length(fct)) {
        message("Converting factor variables to their respective integer values")
        x <- dplyr::mutate(x, across(all_of(fct), ~as.numeric(levels(.x))[.x]))
      }
      return(x)
    })
    
    if(mod_type == "regr.ranger") {
      tasks <- sapply(names(traindat), function(x) {
        TaskRegrST$new(
          id = gsub(".gpkg$", "", basename(x)), 
          backend = cbind(sf::st_drop_geometry(traindat[[x]]), 
                          sf::st_coordinates(traindat[[x]])), 
          target = target, 
          extra_args = list(
            coords_as_features = FALSE, coordinate_names = coordinate_cols[[x]],
            crs = as.character(sf::st_crs(traindat[[x]])$epsg)
          ))
      }, simplify = FALSE, USE.NAMES = TRUE)
      
    } else {
      tasks <- sapply(names(traindat), function(x) {
        TaskClassifST$new(
          id = gsub(".gpkg$", "", basename(x)), 
          backend = cbind(sf::st_drop_geometry(traindat[[x]]), 
                          sf::st_coordinates(traindat[[x]])), 
          target = target, positive = positive,
          extra_args = list(
            coords_as_features = FALSE, coordinate_names = coordinate_cols[[x]],
            crs = as.character(sf::st_crs(traindat[[x]])$epsg)
          ))
      }, simplify = FALSE, USE.NAMES = TRUE)
    }
  } else {
    # Throw error if both X and Y columns are present
    stop("Spatial input is required.")
  }
  
  # Create correlation matrix (mandatory output)
  cor <- lapply(tasks, function(x) {
    cor(as.data.frame(as.data.table(x)) %>% 
          dplyr::select(-all_of(target)), method = "pearson")
  })
  
  # Apply a correlation filter based on input data, must remove variables 
  # with no variation first or calculations fail
  # Reference for cutoff value being inverse:
  # https://github.com/mlr-org/mlr3filters/blob/master/R/FilterFindCorrelation.R
  # Associated pdf:
  # https://cran.r-project.org/web/packages/mlr3filters/mlr3filters.pdf
  # Cutoff is 0.1 (opposite of caret::findCorrelation, where cutoff is 0.9)
  # due to the nature of the mlr3 algorithms:
  # https://github.com/topepo/caret/blob/master/pkg/caret/R/findCorrelation.R
  # Chose 0.1 because it is the default value - can adjust it if need be
  if(filter_correlation) {
    tasks <- lapply(tasks, function(x) {
      x <- x$select(
        names(which(apply(x$data() %>% dplyr::select(-Pres), 2, function(y) 
          var(y, na.rm = TRUE) != 0))))
      cor <- flt("find_correlation", method = "pearson")
      cor$calculate(x)
      x <- x$select(
        as.data.table(cor) %>% 
          dplyr::filter(score >= (1 - corr_cutoff)) %>% 
          dplyr::pull(feature))
      return(x)
    })
  }
  
  # Create pipeop learner, note I had some issues with this, error included 
  # duplicated value in "name" column, but if the mlr3spatiaotempcv package is 
  # not loaded here, this works fine. I submitted this as an issue to the 
  # mlr3pipelines page: 
  # https://github.com/mlr-org/mlr3pipelines/issues/368, which resulted in a PR
  # here: https://github.com/mlr-org/mlr3pipelines/pull/371, which was eventually
  # moved to here: https://github.com/mlr-org/mlr3/issues/470, and another PR 
  # here: https://github.com/mlr-org/mlr3/pull/529
  # The issue is still ongoing as of June 26...I appear to have opened a can of
  # worms with this!
  # Edit Sep. 23, 2020: The issue is now resolved and no more workarounds should
  # be necessary provided that a users packages are up to date.
  
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
    glrn <- GraphLearner$new(graph = filter %>>% po_lrn, 
                             predict_type = type)
    
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
      measure = msr(measures),
      search_space = param_set,
      terminator = trm("none"),
      tuner = tnr("grid_search", resolution = resolution)
    )
    
    # Create a benchmark grid which describes the models to be run
    design <- benchmark_grid(tasks = tasks, learners = feature_sel, 
                             resamplings = rsmp("holdout", ratio = 1))
    
  } else if(feature_selection == "rfe") {
    # Do recursive feature selection, recursive = TRUE will recalculate importance
    # at each iteration
    if(length(tasks[[1]]$feature_names) == 1) {
      f_select <- fs("sequential", max_features = 1) 
      trm <- trm("evals", n_evals = 1)
    } else {
      f_select <- fs("rfe", min_features = 1, feature_number = 1, recursive = TRUE) 
      trm <- trm("stagnation", iters = 5, threshold = 1e-5)
    }
    
    # Create the feature selection learner. Models will be compared at each iteration
    # of a feature being removed. If the models don't improve after 5 iterations, then
    # feature selection stops, using the best evaluated model
    feature_sel <- AutoFSelector$new(
      learner = lrn,
      resampling = resampling_outer,
      measure = msr(measures),
      terminator = trm,
      fselector = f_select,
      store_models = TRUE
    )
    
    design <- benchmark_grid(tasks = tasks, learners = feature_sel, 
                             resamplings = rsmp("holdout", ratio = 1))
    
  } else if (feature_selection == "none") {
    
    f_select <- fs("rfe", min_features = length(tasks[[1]]$feature_names)) 
    trm <- trm("evals", n_evals = 1)
    feature_sel <- AutoFSelector$new(
      learner = lrn,
      resampling = resampling_outer,
      measure = msr(measures),
      terminator = trm,
      fselector = f_select,
      store_models = TRUE
    )
    
    design <- benchmark_grid(tasks = tasks, learners = feature_sel, 
                             resamplings = rsmp("holdout", ratio = 1))
    
  }
  
  bmr <- benchmark(design, store_models = TRUE)
  all_results <- lapply(1:bmr$n_resample_results, function(x) {
      bmr$resample_result(x)$learners[[1]]
  })
  names(all_results) <- names(traindat)
  
  return(list(learners = all_results, tasks = tasks, cor_matrix = cor))
}
