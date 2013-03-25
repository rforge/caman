#setup data structure

setClass("CAMAN.object", representation(dat="matrix", family="character", 
                                        LL="numeric", num.k="numeric", p="numeric", t="numeric",  component.var = "numeric", 
                                        prob="matrix", classification = "numeric", num.obs="numeric", steps = "numeric", 
                                        otherParams = "numeric", BIC = "numeric", VEM_result = "matrix", finalacc = "numeric", cl =  "call", is_metaAnalysis = "numeric"))

setClass("CAMAN.VEM.object", representation(dat="matrix", family="character", 
                                            LL="numeric", num.k="numeric", startk="numeric", p="numeric", t="numeric",  
                                            num.obs="numeric", steps = "numeric", BIC = "numeric", finalacc = "numeric",
                                            otherParams = "numeric", cl =  "call", is_metaAnalysis = "numeric", grid="data.frame", totalgrid="data.frame"))


setClass("CAMAN.glm.object", representation(dat="data.frame", family="character", 
                                            LL="numeric", num.k="numeric", p="numeric", t="numeric",  hetvar = "numeric", prob="matrix", 
                                            classification = "numeric", num.obs="numeric", steps = "numeric", otherParams = "numeric", 
                                            BIC = "numeric", coefMatrix = "data.frame", commonEffect = "numeric", cl =  "call", fittedObs="numeric",
                                            numPara = "numeric", depVar = "character", fixedVar = "character", random = "character",
                                            form="formula", glmModel = "glm", mode="character", residVar= "numeric", idxControl="list", inputData = "data.frame"))

setMethod("show", "CAMAN.object", function(object){
  cat("Computer Assisted Mixture Analysis: \n \n")
  cat("Data consists of", object@num.obs, "observations (rows). \n")
  cat("The Mixture Analysis identified", object@num.k, "component")
  if (length(object@num.k) >0) cat("s")
  cat(" of a", object@family, "distribution: \n \n")
  
  n <- object@num.obs    
  details <- matrix(0, nrow=object@num.k, ncol=2)
  descr_var =""
  if (object@family == "gaussian") {
    tmp <- "mean"
    if (object@is_metaAnalysis == 0) descr_var <- paste("component variance:", object@component.var, "\n")
  }
  else if (object@family == "poisson") tmp <- "lambda"
  else if (object@family == "binomial") tmp <- "prob"
  colnames(details) = c("p", tmp)
  rownames(details) = 1:object@num.k
  details[,1] <- object@p
  details[,2] <- object@t
  cat("DETAILS:\n")
  print(details)
  cat(descr_var,"\n")
  cat("Log-Likelihood:",object@LL,"    ")
  cat("BIC:",object@BIC,"\n")      
})

setMethod("show", "CAMAN.VEM.object", function(object){
  cat("Computer Assisted Mixture Analysis (VEM): \n \n")
  cat("Data consists of", object@num.obs, "observations (rows). \n")
  cat("The VEM-algorithm identified", nrow(object@grid), "grid point")
  if (length(object@num.k) >0) cat("s")
  cat(" with positive support \n \n")
  
  n <- object@num.obs    
  print(object@grid)
  
  cat("Log-Likelihood:",object@LL,"    ")
  cat("BIC:",object@BIC,"\n")      
})


setMethod("show", "CAMAN.glm.object", function(object){
  cat("Computer Assisted Mixture Analysis with covariates: \n \n")
  cat("Data consists of", object@num.obs, "observations (rows). \n")
  cat("The Mixture Analysis identified", object@num.k, "component")
  if (length(object@num.k) >0) cat("s")
  cat(" of a", object@family, "distribution: \n \n")
  
  n <- object@num.obs    
  cat("mixing weights:\n")
  p_tmp <- object@p
  names(p_tmp) = paste("comp.", 1:object@num.k)
  print(p_tmp)
  
  cat("\n Coefficients :\n")
  coefPrint <- round(object@coefMatrix[,1],3)
  names(coefPrint) <- rownames(object@coefMatrix)
  print(coefPrint)
  if (object@family=="gaussian") cat("residual variance:", object@residVar)
  cat("\n")
  
  cat("Log-Likelihood:",object@LL,"    ")
  cat("BIC:",object@BIC,"\n")      
})

