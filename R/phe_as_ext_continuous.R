phe_as_ext_continuous <-
function(phe.gen, min.records=20, return.models=FALSE, confint.level=NA,
         factor.contrasts=NA, my.data, ...) {
  if(!missing(my.data)) data=my.data
  # Retrieve the targets for this loop
  phe=phe.gen[[1]]
  estimated_repeat_length=phe.gen[[2]]
  cov=phe.gen[[3]]
  
  # Subset the data
  d=data %>% select(one_of(na.omit(unlist(c(phe,estimated_repeat_length,cov)))))
  # Turn covariates into a string, if not NA
  if(!is.na(cov[1])) {covariates=paste(cov,collapse=",")}
  else {covariates=NA_character_}

  # Exclude NA phenotype
  d=d[!is.na(d[[phe]]),]
  n_total=nrow(d)
  
  n_cases=NA_integer_
  n_controls=NA_integer_
  beta=NA_real_
  se=NA_real_
  p=NA_real_
  type=NA_character_
  note=""
  model=NA
  formula.string=NA_character_
  expanded_formula=NA_character_
  estimated_repeat_length_expansion=1:length(estimated_repeat_length)
  
  # Drop columns with no variability
  drop.cols = names(d)[sapply(d, function(col) length(unique(col)))<=1]
  if(length(drop.cols)>0) {
    note=paste(note,"[Note: Dropped due to no variability: ",
               paste0(drop.cols,collapse=", "),"]")
    d=select(d, -one_of(drop.cols))
    cov=setdiff(cov,drop.cols)
  }
  
  if(n_total<min.records) {
    note=paste(note,"[Error: <",min.records," complete records]")
  } else if(sum(c(phe,estimated_repeat_length) %in% names(d))!=length(c(phe,estimated_repeat_length))) {
    note=paste(note,"[Error: missing phenotype or estimated_repeat_length]")
  } else {
    # Alter factors to use special contrasts
    if(suppressWarnings(!is.na(factor.contrasts))) {
      d=data.frame(lapply(d,function(x){
        if("factor" %in% class(x)){
          x=droplevels(x)
          contrasts(x)=factor.contrasts(x)
        }
        x}),check.names=F)
    }
    # Create formula
    formula.string=paste0("`",phe,"` ~ `",
                          paste(na.omit(c(estimated_repeat_length,cov)),collapse = "` + `"),'`')
    my.formula = as.formula(formula.string)
    
    # Logistic vs linear
    if(class(d[[phe]]) %in% c("logical","factor") & length(unique(d[[phe]]))==2) {
      type = "logistic"
      # Create logistic model
      n_cases=sum(d[[phe]]==1, na.rm=TRUE)
      n_controls=n_total-n_cases
      if(n_cases<min.records|n_controls<min.records) {
        note=paste(note,"[Error: <",min.records," cases or controls]")
      } else {
        model = glm(my.formula, data=d, family=binomial)
        modsum= summary(model)
        if(model$converged) {
          estimated_repeat_length_expansion=attr(model.matrix(my.formula, data=d),"assign")
          pred_list=which(estimated_repeat_length_expansion %in% 1:length(estimated_repeat_length))
          estimated_repeat_length_expansion=estimated_repeat_length_expansion[pred_list]
          pred_list=grep(estimated_repeat_length,row.names(modsum$coef))
          est_terms=row.names(modsum$coef)[pred_list]
          beta=modsum$coef[pred_list,1]
          se=modsum$coef[pred_list,2]
          p=modsum$coef[pred_list,4]
          expanded_formula=paste0(names(model$coefficients),collapse=" + ")
        } else {
          note=paste(note,"[Error: The model did not converge]")
        }
      }
    } else {
      type = "linear"
      model = glm(my.formula, data=d)
      modsum= summary(model)
      if(model$converged) {
        estimated_repeat_length_expansion=attr(model.matrix(my.formula, data=d),"assign")
        pred_list=which(estimated_repeat_length_expansion %in% 1:length(estimated_repeat_length))
        estimated_repeat_length_expansion=estimated_repeat_length_expansion[pred_list]
        pred_list=grep(estimated_repeat_length,row.names(modsum$coef))
        est_terms=row.names(modsum$coef)[pred_list]
        beta=modsum$coef[pred_list,1]
        se=modsum$coef[pred_list,2]
        p=modsum$coef[pred_list,4]
        expanded_formula=paste0(names(model$coefficients),collapse=" + ")
      } else {
        note=paste(note,"[Error: The model did not converge]")
      }
    }
  }
  
  output=data.frame(phenotype=phe,
                    estimated_repeat_length=est_terms,
                    covariates=covariates,
                    beta=beta, SE=se,
                    p=p, type=type,
                    n_total=n_total, n_cases=n_cases, n_controls=n_controls,
                    formula=formula.string,
                    expanded_formula=expanded_formula,
                    note=note, stringsAsFactors=F)

  # Add confidence intervals if requested
  if(!is.na(confint.level)) {
    if(!is.na(model)[1]){
      suppressMessages(conf<-confint(model,c(1,pred_list),level=confint.level))
      lower=conf[-1,1]
      upper=conf[-1,2]
      if(type=="logistic") {
        lower=exp(lower)
        upper=exp(upper)
      }
    } else {
      lower=NA_real_
      upper=NA_real_
    }
    output$lower=lower
    output$upper=upper
    
    output=output[,c("phenotype","estimated_repeat_length","beta","SE",
                     "lower","upper","p","type",
                     "n_total","n_cases","n_controls",
                     "formula","expanded_formula","note")]
  }
  
  if(return.models) {attributes(output)$model=model}
  attributes(output)$successful.phenotype=ifelse(is.na(p),NA,phe)
  attributes(output)$successful.estimated_repeat_length=ifelse(is.na(p),NA,estimated_repeat_length)
  output
}
