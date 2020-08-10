.allocate <- function (classes, n = 15, method = "random", weights = NULL) {
  
  if((base::length(classes) == 0) | (base::is.null(classes))) {
    stop("Classes are not specified.")
  }
  if(n < 1) {
    stop("n must be greater than 0.")
  }
  if(method == "weighted") {
    if(base::is.null(weights)) {
      stop("Weighted random allocation specified but no weights supplied.")
    }
    else if(!(base::length(classes) == base::length(weights))) {
      stop("Number of classes is not the same as number of weights.")
    }
  }
  if(base::all(base::is.na(classes)) == TRUE) {
    return(base::rep(NA, times = n))
  }
  allocation <- character()
  if(method == "weighted") {
    s <- gtools::rdirichlet(1, weights)
    allocation = base::sample(classes, size = n, replace = TRUE, 
                              prob = s[1, ])
  }
  else if(method == "random") {
    allocation = base::sample(classes, size = n, replace = TRUE)
  }
  else stop("Allocation method \"", method, "\" is unknown")
  return(allocation)
}
