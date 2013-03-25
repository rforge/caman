# TODO
# 
# Author: Johannes
###############################################################################

# TODO aspirin_metaanalyse
# TODO aspirin_metaanalyse




# C:\charite\codeForSpringer\actual caman files\CAMAN\R\CAMAN.r
mixalg = function(obs, weights=NULL, family="gaussian", data=NULL, 
                  pop.at.risk=NULL, var.lnOR=NULL, limit=0.01, 
                  acc=10^(-7), numiter=5000, startk=50){        
    # Performs CAMAN (computer-assisted analysis of mixtures)
    if (family == "gaussian") dens_i = 0
    else if (family == "poisson") dens_i = 1
    else if (family == "binomial") dens_i = 2
    else return("Please enter a valid family of distribution (normal, poisson, binomial)")    
    
    #check data
    if (is.null(data)) data <- data.frame() # no data was given 
    if (!((obs %in% colnames(data))||is.numeric(obs))) stop("obs must be a colname of 'data' or a numeric vector")
    if (!((weights %in% colnames(data))||is.numeric(weights) || is.null(weights))) stop("weights must be a colname of 'data', a numeric vector or 'NULL'")
    if (!((var.lnOR %in% colnames(data))||is.numeric(var.lnOR) || is.null(var.lnOR))) stop("variances must be a colname of 'data', a numeric vector or 'NULL'")

	if (is.null(var.lnOR) ) is_metaAnalysis <- 0 #if variances are not specified, we apply don't perform a meta analysis and estimate the variances
	else is_metaAnalysis <- 1 #variances are given --> meta analysis --> no estimation

    
    n= max(nrow(data), length(obs) )
    datmat <- matrix(1,ncol=4, nrow=n)
    
    if (is.numeric(obs)&&length(obs>1)) datmat[,1] <- obs
    else datmat[,1] <- data[,obs]
    
    #build matrix 'datmat' by reading out the command
    tmpdat <- list(obs, weights, pop.at.risk, var.lnOR)
    for (i in 1:4){
        if (is.character(tmpdat[[i]]) ) datmat[,i] <- data[,tmpdat[[i]]] #colname was given
        else if (is.null(tmpdat[[i]]) ) datmat[,i] <- rep(1,n)   #NULL was given
        else if (is.numeric(tmpdat[[i]]) ) datmat[,i] <- tmpdat[[i]] #a numeric vector was given
        else stop("Data initialization failed...")
    }
    rm(tmpdat)
  #estimate variances  
	if ((sum(rep(1,n) == datmat[,4])== n) && family=="gaussian" && is.null(var.lnOR)) datmat[,4] <- rep(var(datmat[,1]),n)
    
	
    res1 <- .C("caman_C", as.vector(as.double(datmat[,1])), as.vector(as.double(datmat[,2])), as.vector(as.double(datmat[,3])), 
        as.vector(as.double(datmat[,4])), as.integer(n), as.integer(startk), as.integer(dens_i), 
        as.integer(999), as.double(999), rep(as.double(999), 150), rep(as.double(999), 150), 
        as.double(limit), as.double(acc), as.integer(numiter), as.double(c(-999)), 
               rep(as.double(-999), (2*startk + 2)),rep(as.double(-999), 2) ,
               as.integer(is_metaAnalysis)  ,PACKAGE = "CAMAN")
    
    numObs <- sum(datmat[,2])
    k <- res1[[8]]
    p=res1[[10]][1:k]    
    bic <- -2 * res1[[9]] + (2*k - 1) * log(numObs)
    VEM_tmp <- res1[[16]]
    EM_tmp <- res1[[17]]
    finalacc <- c(VEM_tmp[2], EM_tmp[1]) # VEM, EM 
    vem_res <- matrix(VEM_tmp[3: (2*VEM_tmp[1] + 2)], ncol=2)
    totalsteps <- c(res1[[14]][1], EM_tmp[2]) #VEM, EM
    res <- new("CAMAN.object",dat=datmat, family=family, LL=res1[[9]], 
               num.k=k, p=p, t=res1[[11]][1:k], num.obs = numObs, steps=totalsteps, 
               otherParams = c(limit, numiter, acc, startk), BIC = bic, VEM_result = 
                 vem_res, finalacc= finalacc, is_metaAnalysis=is_metaAnalysis)
    if (dens_i == 0) res@component.var=res1[[15]]
    
	#compute posterior prob
	probs <- mix.densDistr(res)
    res@prob <- probs
    res@classification <- apply(probs, 1, which.max)
    
    
    if (res@steps[1] >= res@otherParams[2]) {
      warning("Warning: Solution does not satisfy convergence criterion:\n The last VEM-iteration had a accurency of ",
              res@finalacc[1],". You asked for (acc=)", acc,
              "\n Please increase numiter or decrease acc", sep="")}
if (res@steps[2] >= res@otherParams[2]) {
  warning("Warning: Solution does not satisfy convergence criterion:\n The last EM-iteration had a accurency of ",
          res@finalacc[2],". You asked for (acc=)", acc,
          "\n Please increase numiter or decrease acc", sep="")}
    return(res)
}


anova.CAMAN.object <- function(mix0, mix1, nboot=2500, limit=0.01, acc=10^(-7), numiter=5000, giveBootstrapData=FALSE, giveLikelihood=FALSE){ 
#compute LL-ratio:
# simulate data from mix0 and compare LL of mix0 and mix1 based on this data   
	cl <- match.call()	
    if (nboot<1) return("Please enter a valid number for nboot")
    
    if (mix0@family == "gaussian") dens_i0 = 0
        else if (mix0@family == "poisson") dens_i0 = 1
        else if (mix0@family == "binomial") dens_i0 = 2

    if (mix1@family == "gaussian") dens_i1 = 0
        else if (mix1@family == "poisson") dens_i1 = 1
        else if (mix1@family == "binomial") dens_i1 = 2
        
    datUnpack <- rep(mix0@dat[,1], mix0@dat[,2])
    len <- length(datUnpack)
    obs <- sapply(1:nboot, function(x){simDat(mix0)}) #returns a matrix with bootstrapped observations of m1
	
    #do we need to unpack the data...
    if (sum(mix0@dat[,2]) == nrow(mix0@dat) ){ #no weights --> datUnpacked == mix@dat[,1]
		dataUnpacked=FALSE
    }
    else {  #else: datUnpacked != mix@dat[,1] --> there were weights
		dataUnpacked=TRUE
#if the got unpacked before, we now pack them together again
		tmp_obs <- apply(obs,2, function(x){as.numeric(names(table(x)))})
		tmp_weights <- apply(obs,2, function(x){as.numeric(table(x))})
		SimulatedObs <- tmp_obs
		SimulatedWeights <- tmp_weights
	}	
	
    LL0 <- NULL
    LL1 <- NULL
    for (i in 1:nboot){    
    #perform EM for each bootstrap sample
    #perform EM for mix0
		if (dataUnpacked){ #if the data was packed, we 
#need to handle it in another way because there might be less rows in our datamatrix... 
			obs_i <- SimulatedObs[[i]]
			obs_weights <- SimulatedWeights [[i]]
			tmplen <- length(obs_weights)
			col3 <- rep(1, tmplen) #packing--> use ones as parameter
			col4 <- rep(1, tmplen) #packing--> use ones as parameter
		}
		else { #no packing--> use original parameters for var.lnOR and pop.at.risk
			obs_i <- obs[,i]
			obs_weights <- mix0@dat[,2] #=rep(1, len)
			col3 <- mix0@dat[,3]
			col4 <- mix0@dat[,4]
			tmplen <- len			
		}
        res0 <- .C("mixalg_sub", as.double(obs_i), as.double(obs_weights), as.double(col3), 
            as.double(col4), as.integer(tmplen), as.integer(mix0@num.k), as.integer(dens_i0), 
            as.integer(mix0@num.k), as.double(999), as.double(mix0@p), as.double(mix0@t), as.double(limit), as.double(acc), 
            as.integer(numiter), as.double(c(-999)), as.integer(1),  as.integer(mix0@is_metaAnalysis) , PACKAGE = "CAMAN")   
        res1 <- .C("mixalg_sub", as.double(obs_i), as.double(obs_weights), as.double(col3), 
            as.double(col4), as.integer(tmplen), as.integer(mix1@num.k), as.integer(dens_i1), 
            as.integer(mix1@num.k), as.double(999), as.double(mix1@p), as.double(mix1@t), as.double(limit), as.double(acc), 
            as.integer(numiter), as.double(c(-999)), as.integer(1) , as.integer(mix1@is_metaAnalysis), PACKAGE = "CAMAN")          

		LL0[i] <- res0[[9]]
        LL1[i] <- res1[[9]]
    }
res <- list() 
LL_ratios <- sort(-2*(LL0 - LL1))  #90, 95, 97.5, 99 - quartils
LL_ratio_quartils <- LL_ratios[floor(c(.9, .95, .975, .99)*nboot)]
names(LL_ratio_quartils) <- c(.9, .95, .975, .99)
res$overview <- data.frame(c(as.character(cl$mix0), as.character(cl$mix1)), c(mix0@num.k, mix1@num.k), c(mix0@BIC, mix1@BIC), c(mix0@LL, mix1@LL), c(NA, -2*(mix0@LL - mix1@LL)))
names(res$overview) = c("mixture model","k","BIC","LL", "LL-ratio")
res$"LL ratios in bootstrap-data" <- LL_ratio_quartils
res$"simulated p-value" <- sum(LL_ratios > (-2*(mix0@LL - mix1@LL)) )/nboot

if (giveBootstrapData) {
	if (dataUnpacked) res$BootStrapData <- SimulatedObs
	else res$BootStrapData <- obs
}
if (giveLikelihood) res$LL <- rbind(LL0,LL1) 

return(res)
}

mixalg.paraBoot <- function(mix.estim, nboot=50, limit=0.01, acc=10^(-7), numiter=5000, startk=50, giveBootstrapData=FALSE){
    # Performs a parametric bootstrap on data.
    # Returns the standard deviation of patameters t and p of the given mixture model! 
    if (nboot<1) return("Please enter a valid number for nboot")
    
    if (mix.estim@family == "gaussian") dens_i = 0
        else if (mix.estim@family == "poisson") dens_i = 1
        else if (mix.estim@family == "binomial") dens_i = 2

    datUnpack <- rep(mix.estim@dat[,1], mix.estim@dat[,2])
    len <- length(datUnpack)
    obs <- sapply(1:nboot, function(x){simDat(mix.estim)}) #returns a matrix with bootstrapped observations
    
    #other parameters for the EM-algorithm
    if (sum(mix.estim@dat[,2]) == nrow(mix.estim@dat) ){ #no weights --> datUnpacked == mix@dat[,1] ==> use original parameters for var.lnOR and pop.at.risk
        c2 <- mix.estim@dat[,2]
        c3 <- mix.estim@dat[,3]
        c4 <- mix.estim@dat[,4] 
    }
    else {  #else: there were weights --> datUnpacked != mix@dat[,1] ==> use ones as add. parameters!
        c2 <- rep(1, len)
        c3 <- c2
        c4 <- c2           
    }

    p_mat <- matrix(0, ncol=mix.estim@num.k, nrow=nboot)
    t_mat <- matrix(0, ncol=mix.estim@num.k, nrow=nboot)

cat("\nProgress:\n0.........50.......100% of Bootstraps done.\n")
    for (i in 1:nboot){    
    #perform EM for each bootstrap sample
        res1 <- .C("mixalg_sub", as.double(obs[,i]), as.double(c2), as.double(c3), 
        as.double(c4), as.integer(len), as.integer(mix.estim@num.k), as.integer(dens_i), 
        as.integer(mix.estim@num.k), as.double(999), as.double(mix.estim@p), as.double(mix.estim@t), as.double(limit), as.double(acc), 
        as.integer(numiter), as.double(c(-999)), as.integer(1) , as.integer(mix.estim@is_metaAnalysis), PACKAGE = "CAMAN")   
    if (i %in% seq(0,nboot,max(floor(nboot/20),1))  ) cat ("|")

    p_mat[i, ] <- res1[[10]][1:mix.estim@num.k]
    t_mat[i, ] <- res1[[11]][1:mix.estim@num.k]
    }
cat ("\n")
return(list(p_mat, t_mat, obs=obs))
    sd.p <- apply(p_mat, 2, sd)
    sd.t <- apply(t_mat, 2, sd)
    if (giveBootstrapData) return(list(sd.p = sd.p, sd.t = sd.t, giveBootstrapData = obs))
    return(list(sd.p = sd.p, sd.t = sd.t))
}


mixalg.boot <- function(mix, nboot=500, limit=0.01, acc=10^(-5), numiter=5000, startk=50, returnBootstrapRep= FALSE){
    # Performs a nonparametric bootstrap on data.
    # Returns the computed optimal number of component for each bootstrap replication 
    if (nboot<1) return("Please enter a valid number for nboot")
    
    if (mix@family == "gaussian") dens_i = 0
    else if (mix@family == "poisson") dens_i = 1
    else if (mix@family == "binomial") dens_i = 2

    datUnpack <- rep(mix@dat[,1], mix@dat[,2]) # If there are weights --> unpack data and do not use any weights/ -- freq = 1
    #generate bootstrap samples:

    x.sample <- matrix( rep(datUnpack,nboot), ncol=nboot, byrow=FALSE) #initialization
    xlen <- nrow(x.sample)
    x.sample <- apply(x.sample, 2, function(yy){sample(yy, xlen, replace=TRUE)}) #sample (bootstrap)
	
	#pack data
    tmp_obs <- apply(x.sample,2, function(x){as.numeric(names(table(x)))})
    tmp_weights <- apply(x.sample,2, function(x){as.numeric(table(x))})
    vec_n <- sapply(tmp_obs, length)
    bootSamples <- unlist(tmp_obs)
    bootWeights <- unlist(tmp_weights)
    
    tmpn <- nrow(mix@dat)
    bootVar <- rep(1, length(bootSamples))
    if (sum(mix@dat[,4] == rep(1,tmpn)) != tmpn)
    {  #variances were given --> extract them
        tmpvar <- mix@dat[,4]
        names(tmpvar) <- mix@dat[,1] 
        tmpBootVar <- sapply(tmp_obs, function(x){tmpvar[as.character(x)]})
        bootVar <-  as.numeric(unlist(tmpBootVar))
    }
    bootPopAtRisk <- rep(1, length(bootSamples))
    if (sum(mix@dat[,3] == rep(1,tmpn)) != tmpn)
    {  #popAtRisk were given --> extract them
        tmppop <- mix@dat[,3]
        names(tmppop) <- mix@dat[,1] 
        tmpBootPop <- sapply(tmp_obs, function(x){tmppop[as.character(x)]})
        bootPopAtRisk <- as.numeric(unlist(tmpBootPop))
    }
    
    #permutation step; write observations in a vector (columnwise)    
    res1 <- .C("caman_boot", as.double(bootSamples), as.double(bootWeights), as.double(bootPopAtRisk), 
        as.double(bootVar), as.vector(as.integer(vec_n)), as.integer(startk), as.integer(dens_i), 
        as.integer(999), rep(as.double(999.99), nboot), rep(as.double(999), 150), rep(as.double(999), 150), 
        as.double(limit), as.double(acc), as.integer(numiter), as.double(c(-999)), as.integer(nboot), rep(as.integer(999),nboot), rep(as.double(999.99), nboot), as.integer(mix@is_metaAnalysis),PACKAGE = "CAMAN")
    
    if (returnBootstrapRep) res <- list(dat.bootstrap=x.sample, LL=res1[[9]], numk.boot=res1[[17]], LL_k1 = res1[[18]]) #unpacked bootstrap data is returned! 
    else res <- list(LL=res1[[9]], numk.boot=res1[[17]], LL_k1 = res1[[18]])
    return(res)
}




mixalg.EM <- function(mix = NULL, p, t, obs=NULL, weights=NULL, family="gaussian", 
                      data=NULL, pop.at.risk=NULL, var.lnOR=NULL,  limit=0.01, 
                      acc=10^(-7), numiter=5000, giveProb=TRUE){
    # computes the seconde (EM-) part of the CAMAN algorithm: 
    # use manualy defined values for p, t &  k and 
    #   --> refiend soultion with EM algorithm 
    # return updated estimates for p, t, acc & number of iterations (numiter)
    if (length(p) == length(t) ) num.k = length(t)
    else stop("Please enter valid data for p and t")
    if (!is.null(mix)){
        family = mix@family
        datmat = mix@dat   
		is_metaAnalysis = mix@is_metaAnalysis
    }
    else{     
	    #check data
	    if (is.null(data)) data <- data.frame() # no data was given 
	    if (!((obs %in% colnames(data))||is.numeric(obs))) stop("obs must be a colname of 'data' or a numeric vector")
	    if (!((weights %in% colnames(data))||is.numeric(weights) || is.null(weights))) stop("weights must be a colname of 'data', a numeric vector or 'NULL'")
	    if (!((var.lnOR %in% colnames(data))||is.numeric(var.lnOR) || is.null(var.lnOR))) stop("weights must be a colname of 'data', a numeric vector or 'NULL'")
 
 	    if (is.null(var.lnOR) ) is_metaAnalysis <- 0 #if variances are not specified, we apply don't perform a meta analysis and estimate the variances
    	else is_metaAnalysis <- 1 #variances are given --> meta analysis --> no estimation

	    
	    n= max(nrow(data), length(obs) )
	    datmat <- matrix(1,ncol=4, nrow=n)
	    
	    if (is.numeric(obs)&&length(obs>1)) datmat[,1] <- obs
	    else datmat[,1] <- data[,obs]
		
		
	    
	    #build matrix 'datmat' by reading out the command
	    tmpdat <- list(obs, weights, pop.at.risk, var.lnOR)
	    for (i in 1:4){
	        if (is.character(tmpdat[[i]]) ) datmat[,i] <- data[,tmpdat[[i]]] #colname was given
	        else if (is.null(tmpdat[[i]]) ) datmat[,i] <- rep(1,n)   #NULL was given
	        else if (is.numeric(tmpdat[[i]]) ) datmat[,i] <- tmpdat[[i]] #a numeric vector was given
	        else stop("Data initialization failed...")
	    }
		#estimate variances  
		if ((sum(rep(1,n) == datmat[,4])== n) && family=="gaussian" && is.null(var.lnOR)) datmat[,4] <- rep(var(datmat[,1]),n)
		
		rm(tmpdat)
    }
	if (family == "gaussian") dens_i = 0
	else if (family == "poisson") dens_i = 1
	else if (family == "binomial") dens_i = 2    
	else stop("Please enter a valid density distribution (gaussian, poisson, binomial)")
	
	
    res1 <- .C("mixalg_sub", as.double(datmat[,1]), as.double(datmat[,2]), as.double(datmat[,3]), 
    as.double(datmat[,4]), as.integer(nrow(datmat)), as.integer(num.k), as.integer(dens_i), 
    as.integer(num.k), as.double(999), as.double(p), as.double(t), as.double(limit), as.double(acc), 
    as.integer(numiter), as.double(c(-999)), as.integer(1), as.integer(is_metaAnalysis) ,PACKAGE = "CAMAN")
    
    numObs <- sum(datmat[,2])
    bic <- -2 * res1[[9]] + (2*num.k - 1) * log(numObs)
    totalsteps <- c(NA, res1[[14]]) #VEM, EM
    finalacc <- c(NA, res1[[13]])
    
    res <- new("CAMAN.object",dat=datmat, family=family, LL=res1[[9]], num.k=num.k, p=res1[[10]], t=res1[[11]], 
			num.obs = numObs, steps=totalsteps, otherParams = c(limit, numiter, acc, startk=num.k), BIC = bic, 
			VEM_result = matrix(), finalacc= finalacc, is_metaAnalysis = is_metaAnalysis)
    if (dens_i == 0) {
		if  (num.k>1) res@component.var=res1[[ 15 ]]
		else res@component.var = var(res@dat[,1])
	}
	
    #compute posterior probabilities
    probs <- mix.densDistr(res)
    res@prob <- probs
    res@classification <- apply(probs, 1, which.max)
    
return(res)
}


mixalg.VEM <- function(mix = NULL, obs=NULL, weights=NULL, data=NULL, pop.at.risk=NULL, var.lnOR=NULL, family="gaussian", limit=0.01, acc=10^(-7), numiter=5000, startk=50){
    # computes the first part of the CAMAN algorithm:
    #    1. construct grid of potential subpopulation means 
    #    2. calculation of mxing kernel density
    #    3. VEM algorithm 
    # returns estimates for parameters t & weights p
	
    if (!is.null(mix)){
		family = mix@family
		datmat = mix@dat   
		is_metaAnalysis = mix@is_metaAnalysis
	}
	else{     
		#check data
		if (is.null(data)) data <- data.frame() # no data was given 
		if (!((obs %in% colnames(data))||is.numeric(obs))) stop("obs must be a colname of 'data' or a numeric vector")
		if (!((weights %in% colnames(data))||is.numeric(weights) || is.null(weights))) stop("weights must be a colname of 'data', a numeric vector or 'NULL'")
		if (!((var.lnOR %in% colnames(data))||is.numeric(var.lnOR) || is.null(var.lnOR))) stop("weights must be a colname of 'data', a numeric vector or 'NULL'")
		
		if (is.null(var.lnOR) ) {
			is_metaAnalysis <- 0 #if variances are not specified, we apply don't perform a meta analysis and estimate the variances
		}
		else {
			is_metaAnalysis <- 1 #variances are given --> meta analysis --> no estimation
		}
		
		n = max(nrow(data), length(obs) )
		datmat <- matrix(1,ncol=4, nrow=n)
			
		if (is.numeric(obs)&&length(obs>1)) datmat[,1] <- obs
		else datmat[,1] <- data[,obs]
		
		#build matrix 'datmat' by reading out the command
		tmpdat <- list(obs, weights, pop.at.risk, var.lnOR)
		for (i in 1:4){
			if (is.character(tmpdat[[i]]) ) datmat[,i] <- data[,tmpdat[[i]]] #colname was given
			else if (is.null(tmpdat[[i]]) ) datmat[,i] <- rep(1,n)   #NULL was given
			else if (is.numeric(tmpdat[[i]]) ) datmat[,i] <- tmpdat[[i]] #a numeric vector was given
			else stop("Data initialization failed...")
		}
		#estimate variances  
		if ((sum(rep(1,n) == datmat[,4])== n) && family=="gaussian" && is.null(var.lnOR)) datmat[,4] <- rep(var(datmat[,1]),n)
		
		rm(tmpdat)
	}
    if (family == "gaussian") dens_i = 0
    else if (family == "poisson") dens_i = 1
    else if (family == "binomial") dens_i = 2
    else stop("Please enter a valid density distribution (normal, poisson, binomial)")

    num.k <- startk #min(sum(datmat[,2]),startk)
    p <- rep(999, num.k)
    t <- rep(999, num.k)
	#perform the EM algorithm for the given (simulated) data
    res1 <- .C("mixalg_sub", as.double(datmat[,1]), as.double(datmat[,2]), as.double(datmat[,3]), 
    as.double(datmat[,4]), as.integer(nrow(datmat)), as.integer(num.k), as.integer(dens_i), 
    as.integer(num.k), as.double(999), as.double(p), as.double(t), as.double(limit), as.double(acc), 
    as.integer(numiter), as.double(c(-999)), as.integer(0) , as.integer(is_metaAnalysis), 
	as.double(rep(-999.9,num.k)),PACKAGE = "CAMAN")

    LL <- res1[[9]]
	numObs <- sum(datmat[,2])
    bic <- -2 * res1[[9]] + (2*num.k - 1) * log(numObs)
    grid <- data.frame(p=res1[[10]], t=res1[[11]])
	grid <- grid[grid$p>0,]
	rownames(grid) <- as.character(1:nrow(grid))
	finalacc <- res1[[13]][1]
	totalsteps <- res1[[14]][1]
    
	totalgrid <- data.frame(p=res1[[10]], t=res1[[11]], gradient=res1[[18]])

	res <- new("CAMAN.VEM.object",dat=datmat, family=family, LL=LL, 
			num.k=num.k, num.obs = numObs, steps=totalsteps, otherParams = c(limit, numiter, acc, startk), 
			BIC = bic, finalacc= finalacc, startk=startk, grid=grid, totalgrid=totalgrid)
	
	
	
	return(res)    
}


#
#mix.densDistr <- function(mix){
#    # computes the probability for each observation (1..n - row of mix@dat) belonging to each component (1..k)
#    # returns a matrix of dimension n x k
#    dat <- mix@dat[,1]
#    res <- matrix(ncol=mix@num.k, nrow=length(dat))
#    p <- mix@p
#   if (mix@family == "gaussian") {
#        mu <- mix@t
#        mix.sd <- sqrt(mix@component.var)
#        for (i in 1:mix@num.k) res[,i] <- sapply(dat,
#			function(x){p[i]*dnorm(x, mu[i], mix.sd ) / sum(p*dnorm(x, mu,
#			mix.sd ))})
#        }
#   if (mix@family == "binomial") {
#        prob <- mix@t
#        popAtRisk <- mix@dat[,3]
#        for (i in 1:mix@num.k) res[,i] <- apply(cbind(dat, popAtRisk), 1,
#			function(x){p[i]*dbinom(x[1], x[2], prob[i]) / sum(p*dbinom(x[1],
#			x[2], prob))})
#        }
#   if (mix@family == "poisson") {
#        lambda <-  mix@t
#        popAtRisk <- mix@dat[,3]
#        for (i in 1:mix@num.k) res[,i] <- apply(cbind(dat, popAtRisk), 1, 
#					function(x){p[i]*dpois(x[1], x[2]* lambda[i]) / 
#								sum(p*dpois(x[1], x[2] * lambda))})
#        }
#    return(res)
#}


getFDR <- function(dat, threshold=.7, idxNotDiff=1 ){
   #computes False Discovery Rate, etc. 
   tau0 <- dat@prob[,idxNotDiff] #p(not differentail genes)
   tau1 <- 1 - dat@prob[,idxNotDiff] #p(differentail genes)
   n <- nrow(dat@prob)
   fdr.hat <- sum(tau0 * (tau0 <= threshold)) / sum((tau0 <= threshold) ) 
   fndr.hat <- sum( (tau1) * (tau0 >= threshold))  / (n - sum(tau0 <= threshold) )
   fpr.hat <- sum(tau0 * (tau0 <= threshold)) / sum(tau0)
   fnr.hat <- sum(tau1 * (tau0> threshold) ) / sum(tau1)
return(list(FDR = fdr.hat, FNDR=fndr.hat, FPR = fpr.hat, FNR=fnr.hat) )
}


simDat <- function(mix){
    #simulate data for parametric bootstrap & compareMixModels
    #if weights are != 1, other parameters needs to be == 1!!! (--> otherwise, data cannot be unpacked reasonably!) 
     
    n= mix@num.obs #number of observations 

    k= mix@num.k
    myunif <- runif(n)
    z <- matrix(0, ncol=k, nrow=n)
    for (i in 1:k)
		if (i==1) z[,i] <- as.integer(sum(mix@p[1]) > myunif )
	    else z[,i] <- as.integer(sum(mix@p[1:(i-1)]) < myunif & (sum(mix@p[1:i]) >= myunif) )
    if (k==1) z <-matrix(1, ncol=k, nrow=n)
	#print (z)
    # z is a matrix that randomly assigns each observation (row) to a component (column) 
    # --> each row consists of 1x1 and (k-1)x0
    x = matrix(0,ncol=mix@num.k, nrow=n)

    if (mix@family == "gaussian"){
        if (sum(mix@dat[,4] == rep(1, n) ) == n) #dat[,4] only consists of ones
        tmpsd <-  sqrt(mix@component.var)
        else tmpsd <- mix@dat[,4]
        
        for (i in 1:k)
           x[,i] <- z[,i] * rnorm(n, mean=mix@t[i], sd= sqrt(tmpsd))
        res = apply(x,1,sum)  #a row consists of one number != 0 and (k-1) numbers == 0
    }
    else if (mix@family == "poisson"){
        Ei <- mix@dat[,3]
        for (i in 1:k)
           x[,i] <- z[,i] * rpois(n, mix@t[i]*Ei)
        res = apply(x,1,sum)  #a row consists of one number != 0 and (k-1) numbers == 0 
    }
    if (mix@family == "binomial"){
        tmpsize <- mix@dat[,3]
        for (i in 1:k)
           x[,i] <- z[,i] * rbinom(n, tmpsize, mix@t[i]) 
        res = apply(x,1,sum)  #a row consists of one number != 0 and (k-1) numbers == 0                 
    }
return(res)
}





summary.CAMAN.object <- function(object){
	cl <- match.call()
	cat("Summary of a Computer Assisted Mixture Analysis: \n \n")
	cat("Data consists of", object@num.obs, "observations (rows). \n")
	cat("The Mixture Analysis identified", object@num.k, "component")
	if (length(object@num.k) >0) cat("s")
	cat(" of a", object@family, "distribution: \n \n")
	
	details <- matrix(0, nrow=object@num.k, ncol=2)
	descr_var <- ""
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
	cat(descr_var)
	cat("\n \n")
	
	cat("Classification:\n")
    if (object@num.obs > 20 || object@num.k>8) cat("The classification matrix is too big to visualize here. \n type ",as.character(cl[2]),"@prob to watch the probability matrix \n or type ",as.character(cl[2]),"@classification to watch the \n class labeling (each row was assigned to its most likely component!)\n", sep="")
	else {
		cat("Classification matrix:\n")
		print(object@prob)
		cat("Class labeling:\n")
		print(object@classification)
	}
	cat("\n \nnumber of VEM-iterations done:",object@steps[1],"\n")        
	cat("alteration within the final VEM-iteration step:",object@finalacc[1],"\n")
	cat("number of EM-iterations done:",object@steps[2],"\n")        
	cat("alteration within the final EM-iteration step:",object@finalacc[2],"\n \n")
	cat("Log-Likelihood:",object@LL,"    ")
	cat("BIC:",object@BIC,"\n \n")
	
	# @otherParams = c(limit, numiter, acc, startk)
	
	cat("User-defined parameters:\n")
	cat("   max number of iterations:",object@otherParams[2],"\n")
	cat("   limit for combining components:",object@otherParams[1],"\n")
	cat("   threshold for converging:",object@otherParams[3],"\n")
	cat("   number of grid points (startk):",object@otherParams[4],"\n")
}


#some abbrevated commands
mixboot <- mixalg.boot
mixalg.Boot <- mixalg.boot
mixpboot <- mixalg.paraBoot
mix.anova <- anova.CAMAN.object
