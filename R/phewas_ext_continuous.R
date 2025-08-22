#' PheWAS Extension for Continuous Predictors
#'
#' @export

phewas_ext_continuous <-
function(data, phenotypes, estimated_repeat_length, covariates=NULL,
         cores=1, min.records=20, return.models=FALSE, confint.level=NA,
         factor.contrasts=NA, ...) {
  
  library(parallel)
  
  # build phenotype-predictor-covariate combinations
  combos <- expand.grid(phenotypes, estimated_repeat_length, stringsAsFactors=FALSE)
  combos <- split(combos, seq(nrow(combos)))
  combos <- lapply(combos, function(x) list(x[,1], x[,2], covariates))
  
  # Run in parallel
  results <- mclapply(combos, phe_as_ext_continuous,
                      min.records=min.records,
                      return.models=return.models,
                      confint.level=confint.level,
                      factor.contrasts=factor.contrasts,
                      my.data=data,
                      mc.cores=cores)
  
  results <- do.call(rbind, results)
  results
}
