# Rewrite of model generation for 08 script

model_gen_mlr3_binary_simple <- function(
  traindat, 
  target, 
  positive = NULL, # The ID of the positive class
  filter_correlation = FALSE,
  corr_cutoff = 0.9,
  folds = 10, repeats = 5, 
  covariates,
  outDir = "./predicted",
  id, rfe = FALSE) {
  
  # Load required packages
  invisible(
    sapply(c("tidyverse", "lgr", "mlr3verse", "mlr3fselect", "ggpubr",
             "mlr3spatiotempcv", "mlr3tuning"), require, character.only = TRUE)[0])
  
  
  # Use external parallelism with future package and set amount of logging to print
  future::plan("multisession")
  def_mlr3_lgr <- get_logger("mlr3")$threshold
  def_bbotk_lgr <- get_logger("bbotk")$threshold
  get_logger("mlr3")$set_threshold("warn")
  get_logger("bbotk")$set_threshold("warn")
  
  # Get the names of the coordinate columns
  coordinate_cols <- colnames(st_coordinates(traindat))
  
  # Create the mlr3 "task" from the sf dataframe (note: process may change in 
  # later versions of mlr3spatiotempcv package)
  task <- TaskClassifST$new(
    id = "cross_validation", 
    backend = cbind(sf::st_drop_geometry(traindat), 
                    sf::st_coordinates(traindat)) %>% 
      dplyr::rename(x = coordinate_cols[1], y = coordinate_cols[2]), 
    target = target, positive = positive,
    extra_args = list(
      coords_as_features = FALSE, coordinate_names = c("x", "y"),
      crs = as.character(sf::st_crs(traindat)$proj4string)))
  
  # Write input data to directory
  write.csv(cbind(st_drop_geometry(traindat), st_coordinates(traindat)), 
            file.path(outDir, "input_data.csv"), row.names = FALSE)
  
  # Create correlation matrix (mandatory output)
  cor_dat <- as.data.table(task) %>% 
    dplyr::select(-any_of(target)) %>% 
    dplyr::select(where(is.numeric))
  
  # Define the variable names from original dataset
  vars <- names(st_drop_geometry(traindat))[names(st_drop_geometry(traindat)) != target]
  
  # If there is data in the correlation matrix, produce correlation outputs
  if(ncol(cor_dat) > 0) {
    
    # Create correlation matrix
    cor_mat <- cor(cor_dat, method = "pearson")
    
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
    
    # Create similar task to the original but remove the x and y variables to 
    # prevent them from being used in the correlation calculations
    task_cor <- task$clone(deep = FALSE)
    task_cor <- task_cor$select(
      names(which(apply(
        task_cor$data() %>% dplyr::select(-any_of(target)), 
        2, function(y) var(y, na.rm = TRUE) != 0))))
    cor <- flt("find_correlation", method = "pearson")
    cor$calculate(task_cor)
    cor_mat_out <- as.data.frame(cor_mat) %>% 
      rownames_to_column("feature") %>% 
      merge(as.data.table(cor)) %>% 
      dplyr::mutate(score = 1 - score) %>% 
      dplyr::rename(pearson_score = score)
    
    # Remove correlated variables if necessary
    if(filter_correlation) {
      task <- task$select(
        as.data.table(cor) %>% 
            dplyr::filter(score >= (1 - corr_cutoff)) %>% 
            dplyr::pull(feature))
    }
    
    # Save correlation image and matrix
    cor_out <- function(cor_matrix) {
      pal <- colorRampPalette(c("green", "white", "red")) (20)
      png(file.path(outDir, "correlation_matrix.jpg"),
          height = 960, width = 960)
      corrplot(cor_matrix)
      dev.off()
      return(invisible())
    }
    cor_out(cor_mat)
    write.csv(cor_mat_out, file.path(outDir, "correlation_matrix.csv"), row.names = FALSE)
  }
  
  # Resampling for raw (cross validation) and class balancing models
  resamp_raw <- rsmp("repeated_cv", folds = folds, repeats = repeats)
  resamp_raw$id <- paste0(folds, "_fold_repeated_", repeats, "_times")
  
  # Resampling for spatial models
  resamp_sp <- rsmp("repeated_spcv_coords", repeats = repeats, folds = folds)
  resamp_sp$id <- paste0("spatial_", folds, "_fold_repeated_", repeats, "_times")
  
  # Follow tuneRF code from randomForest package
  mtryStart <- floor(sqrt(length(task$feature_names)))
  stepFactor <- 2
  
  # Create vector of mtrys to tune model with (from tuneRF in randomForest)
  mtrys <- sort(unique(c(mtryStart, unlist(lapply(c("left", "right"), function(x) {
    mtryCur <- mtryStart
    mtryOld <- Inf
    out <- c()
    while(mtryCur != mtryOld) {
      mtryOld <- mtryCur
      mtryCur <- if(x == "left") {
        max(1, ceiling(mtryCur / stepFactor))
      } else {
        min(length(task$feature_names), floor(mtryCur * stepFactor))
      }
      out <- c(out, mtryCur)
      if(mtryCur == mtryOld) break
    }
    return(unique(out))
  })))))
  
  # Create batch of learners from mtrys
  lrns <- lapply(mtrys, function(x) {
    lrn("classif.ranger", num.threads = 1L,
        importance = "impurity", predict_type = "prob", mtry = x,
        id = paste0("mtry_", x))
  })
  
  # Create non-spatial class balanced models
  design_raw <- benchmark_grid(task, lrns, resamp_raw)
  
  # Create spatial model
  task_sp <- task$clone(deep = FALSE)
  task_sp$id <- "spatial_cross_validation"
  design_sp <- benchmark_grid(task_sp, lrns, resamp_sp)
  
  # Combine spatial and non-spatial models into single object
  design <- rbind(design_raw, design_sp)
  tasks <- list(task, task_sp)
  
  # Generate appropriate performance measures - get all twoclass probability measures
  all_msrs <- suppressWarnings(as.data.table(mlr_measures)) %>% 
    dplyr::filter(task_type == "classif", predict_type == "prob",
                  task_properties == "twoclass")
  
  # Also generate other measurements - note: cannot use OOB becauase models
  # need to be saved in order to do that. Will calculate OOB in saved models later
  prob_measures <- all_msrs$key
  extra_measures <- c("classif.acc", "classif.mcc", 
                      "classif.sensitivity", "classif.specificity",
                      "classif.fbeta")
  
  measures <- lapply(c(prob_measures, extra_measures), function(x) 
    msr(x, id = sub('.*\\.', '', x)))
  
  design_tsks <- sapply(design$task, function(y) y$id)
  
  ### ITERATE THROUGH EACH TASK - MAKES MORE SENSE AND CAN SAVE MEMORY
  run <- lapply(tasks, function(x) {
    cat("Generating", x$id, "models\n")
    out_subdir <- file.path(outDir, x$id)
    dir.create(out_subdir, showWarnings = FALSE)
    
    # Filter design to keep just tasks in iteration, and run benchmark of each learner
    design_flt <- design[which(design_tsks == x$id), ]
    inst <- benchmark(design_flt, store_models = FALSE, store_backends = FALSE)
    
    # Get measures for each model run
    all_results <- inst$score(measures) %>% 
      dplyr::select(nr, sapply(measures, function(x) x$id))
    
    # Create index that defines the iteration and mtry
    mtry_index <- dplyr::select(inst$print(), nr, learner_id) %>%
      group_by(nr) %>%
      slice_head(n = 1) %>%
      ungroup()
    
    # Create a plot used for output later
    mtry_plot <- autoplot(inst, measure = msr("classif.fbeta"))
    
    # Aggregate with mean and variance scores. Write out this table to the subdir
    mtry_results <- group_by(all_results, nr) %>% 
      dplyr::summarise(across(sapply(measures, function(x) x$id), 
                              ~mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
                       across(sapply(measures, function(x) x$id),
                              ~var(.x, na.rm = TRUE), .names = "var_{.col}"), 
                       count_fbeta = sum(!is.na(fbeta)), .groups = "drop") %>% 
      merge(mtry_index) %>% 
      dplyr::relocate(learner_id)
    write.csv(mtry_results, file.path(out_subdir, "tune_mtry_results.csv"), row.names = FALSE)
    
    # Remove mtry runs that had any aggregated NA scores
    if(nrow(mtry_results %>% drop_na()) > 0) 
      mtry_results <- mtry_results %>% drop_na()
    
    # Select the best mtry run from what is left. Use counts of highest number
    # of reported fbeta scores, then select the ~2 mtry runs with the lowest 
    # summed variance 
    best_mtry_var <- mtry_results %>% 
      dplyr::mutate(var_sum = rowSums(.[, paste0("var_", sapply(measures, function(x) x$id))])) %>% 
      dplyr::slice_max(count_fbeta, prop = ifelse(nrow(.) > 2, 0.75, 1)) %>% 
      dplyr::slice_min(var_sum, n = max(2, ceiling(0.25 * nrow(.)))) %>% 
      dplyr::pull(learner_id)
    
    # From the above mtry values, select the one with the highest fbeta and
    # lowest bbrier scores. Outputs single mtry value
    best_mtry <- dplyr::filter(mtry_results, learner_id %in% best_mtry_var) %>% 
      dplyr::slice_max(mean_fbeta) %>% 
      dplyr::slice_min(mean_bbrier, with_ties = FALSE) %>% 
      dplyr::pull(learner_id)
    
    # Get the corresponding model run iteration from the mtry value
    best_inst <- mtry_index %>% 
      dplyr::filter(learner_id == best_mtry) %>% 
      dplyr::pull(nr)
    
    # Use the best mtry value to select the corresponding best learner
    best_lrn <- inst$learners[learner_id == best_mtry]$learner[[1]]
    
    # With the best tune defined, perform resampling and select best model
    # using similar parameters
    rr <- resample(x, best_lrn, design_flt$resampling[[1]], store_models = TRUE, store_backends = FALSE)
    all_rr_results <- rr$score(c(measures, msr("oob_error")))
    if(nrow(all_rr_results %>% drop_na()) > 0)
      all_rr_results <- all_rr_results %>% drop_na()
    
    best_rr_result <- all_rr_results %>% 
      dplyr::slice_max(fbeta, n = ceiling(nrow(.) / 10)) %>% 
      dplyr::slice_min(bbrier, with_ties = FALSE)
    
    mean_rr_results <- all_rr_results %>% 
      dplyr::select(which(apply(., 2, function(x) !inherits(x[[1]], "R6"))),
                    -any_of("iteration")) %>% 
      group_by(task_id, learner_id, resampling_id) %>% 
      dplyr::summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE), 
                              .names = "mean_{.col}"), .groups = "drop")
    
    # Extract model and iteration ID
    model <- best_rr_result$learner[[1]]$model
    iter <- best_rr_result$iteration
    
    # Create table of tuning results for the best model that can be written
    tune_results_model <- best_rr_result %>%
      dplyr::select(which(apply(., 2, function(x) !inherits(x[[1]], "R6"))))
    
    ### IF RFE NEEDS TO BE DONE, DO IT HERE ###
    
    # Some model outputs: Confusion matrices and importance values
    confusion_best <- rr$predictions()[[iter]]$confusion
    confusion_whole <- rr$prediction()$confusion
    importance <- data.frame(gini = best_rr_result$learner[[1]]$importance()) %>% 
      rownames_to_column("layer")
    
    conf_best_out <- format(ftable(confusion_best), quote = FALSE, method = "non.compact")
    conf_whole_out <- format(ftable(confusion_whole), quote = FALSE, method = "non.compact")
    
    # Best testing and training datasets
    rr_reps <- unique(rr$resampling$instance$rep)
    rr_folds <- unique(rr$resampling$instance$fold)
    rr_resampling <- data.frame(iters = 1:rr$resampling$iters) %>% 
      dplyr::mutate(reps = rr$resampling$repeats(iters)) %>% 
      group_by(reps) %>% 
      dplyr::mutate(folds = 1:n()) %>% 
      ungroup()
    
    train_whole <- bind_rows(lapply(rr_reps, function(y) {
      bind_rows(lapply(rr_folds, function(z) {
        train_set_ids <- rr$resampling$instance[rep == y & fold != z]
        train_set <- x$data(train_set_ids$row_id) %>% 
          dplyr::mutate(rep = y, fold = z, row_id = train_set_ids$row_id)
      }))
    }))
    
    test_whole <- bind_rows(lapply(rr_reps, function(y) {
      bind_rows(lapply(rr_folds, function(z) {
        test_set_ids <- rr$resampling$instance[rep == y & fold == z]
        test_set <- x$data(test_set_ids$row_id) %>% 
          dplyr::mutate(rep = y, fold = z, row_id = test_set_ids$row_id)
      }))
    }))
    
    train_best <- cbind(model_set = "training", 
                        x$data(rr$resampling$train_set(iter)),
                        row_id = rr$resampling$train_set(iter))
    test_best <- cbind(model_set = "testing", 
                       x$data(rr$resampling$test_set(iter)),
                       row_id = rr$resampling$test_set(iter))
    
    best_input <- bind_rows(train_best, test_best) %>% 
      dplyr::mutate(rep = rr_resampling[iter, ]$reps,
                    fold = rr_resampling[iter, ]$folds) %>% 
      dplyr::relocate(row_id, .after = last_col()) %>% 
      arrange(row_id)
    
    # Capture text of ranger model
    model_txt <- utils::capture.output(model)
    
    # Write remaining outputs to subdir
    write.table(conf_best_out, file.path(out_subdir, "confusion_best_fold.csv"), 
                row.names = FALSE, col.names = FALSE, sep = ",")
    write.table(conf_whole_out, file.path(out_subdir, "confusion_all.csv"), 
                row.names = FALSE, col.names = FALSE, sep = ",")
    cat(model_txt, file = file.path(out_subdir, "best_model.txt"), sep = "\n")
    write.csv(tune_results_model, file.path(out_subdir, "model_metrics_best.csv"), row.names = FALSE)
    write.csv(importance, file.path(out_subdir, "variable_importance.csv"), row.names = FALSE)
    write.csv(train_whole, file.path(out_subdir, "resampling_train_data.csv"), row.names = FALSE)
    write.csv(test_whole, file.path(out_subdir, "resampling_test_data.csv"), row.names = FALSE)
    write.csv(best_input, file.path(out_subdir, "best_input_data.csv"), row.names = FALSE)
    saveRDS(model, file.path(out_subdir, "best_model.rds"))
    
    return(list(mtry_plot = mtry_plot, model = model, metrics = mean_rr_results))
  })
  
  # Create image for visual comparison
  bmr_compare <- ggarrange(plotlist = lapply(run, function(x) {
    p <- x$mtry_plot
    p$coordinates$limits$y <- c(0, 1)
    return(p)})) +
    ggsave(file.path(outDir, "run_compare.jpg"), width = 10, height = 5, dpi = 300)
  
  # Extract models for output
  all_best <- bind_rows(lapply(run, "[[", "metrics"))
  m <- lapply(run, "[[", "model") %>% setNames(all_best$task_id)
  best <- all_best %>% 
    dplyr::slice_max(mean_fbeta, n = max(2, nrow(.) * (1/3))) %>% 
    dplyr::slice_min(mean_bbrier, with_ties = FALSE)
  
  write.csv(all_best, file.path(outDir, "run_compare.csv"), row.names = FALSE)
  
  best_m <- m[[best$task_id]]
  
  # Return original logger levels
  get_logger("mlr3")$set_threshold(def_mlr3_lgr)
  get_logger("bbotk")$set_threshold(def_bbotk_lgr)
  
  future::plan("default")
  return(best_m)
}
