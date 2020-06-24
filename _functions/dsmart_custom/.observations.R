.observations <-
function (observations, covariates) 
{
    base::colnames(observations) <- c("x", "y", "soil_class")
    o.covariates <- observations %>% 
        cbind(terra::extract(covariates, observations[, c("x", "y")]))
    meta <- list(realisation = numeric(length = nrow(observations)), 
        type = base::rep("actual", nrow(observations)), sampling = base::rep("observed", 
            nrow(observations)), allocation = base::rep("observed", 
            nrow(observations)))
    obs <- cbind(as.data.frame(meta), o.covariates) %>% 
        dplyr::filter(complete.cases(.))
    return(obs)
}
