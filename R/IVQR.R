ivqr <- function(formula, taus=0.5, data, grid, gridMethod="Default", ivqrMethod="iqr",
				qrMethod="br"){
	
	formula <- Formula(formula) 
	#!
	# warning("cannot exclude intercept")
	
	# --- 
	# Check if taus has appropriate range
	eps <- .Machine$double.eps ^ (2/3)
    if (length(taus) > 1){
		if (any(taus < 0) || any(taus > 1)) 
		stop("invalid taus:  taus should be >= 0 and <= 1")
	}
	if (any(taus == 0)) taus[taus == 0] <- eps
	if (any(taus == 1)) taus[taus == 1] <- 1 - eps

	# By default, we use D projection on (X,Z) as instrument
	# Add the instruments to the data and update the formula respectively
	XZ <- formula(formula, lhs = 1, rhs = c(2,3), collapse = TRUE)
	XZ <- model.matrix(XZ,data)
	D  <- update( formula(formula,lhs = 1, rhs = 1), . ~ . - 1)
	D  <- model.matrix(D,data)

	# If my code doesn't run for multiple endg var, uncomment below
	#! ---
	# if (dim(D)[2]==1){ # Make sure my code works at least for dim(D) == 1
	# 	D_hat <- XZ %*% solve(t(XZ) %*% XZ) %*% t(XZ) %*% D
	# }else{
	# 	D_hat <- D
	# 	for (i in 1:dim(D)[2]){
	# 		D_hat[,i] <- XZ %*% solve(t(XZ) %*% XZ) %*% t(XZ) %*% D[,i]
	# 	}
	# }	
	D_hat <- D
	for (i in 1:dim(D)[2]){
		D_hat[,i] <- XZ %*% solve(t(XZ) %*% XZ) %*% t(XZ) %*% D[,i]
	}
	#! ---
	if (any(grepl(".ivqr_dhat",colnames(data)))){
		stop("No names of variables in the data set should include .ivqr_dhat")
	}
	dhat_formula <- c("~")
	for (i in 1:dim(D)[2]){
		data[,paste(".ivqr_dhat",i,sep="")] <- D_hat[,i]
		dhat_formula <- cbind(dhat_formula,paste(".ivqr_dhat",i,sep=""))
		if (i < dim(D)[2]) dhat_formula <- cbind(dhat_formula,"+")
	}

	dhat_formula <- formula(paste(dhat_formula,collapse=""))
	
	iqr_formula <- as.Formula(formula(formula,lhs=1,rhs=1),
							  dhat_formula,formula(formula,lhs=0,rhs=3))

	# ---
	# Define variables to store coef estimates and names
	coef <- list()

	endg_varnames <- formula(iqr_formula, lhs = 1, rhs = 1)
	endg_varnames <- update(endg_varnames,.~.-1)
	inst_varnames <- formula(iqr_formula, lhs = 1, rhs = 2)
	inst_varnames <- update(inst_varnames,.~.-1)
	exog_varnames <- formula(iqr_formula, lhs = 1, rhs = 3)	

	D <- model.matrix(endg_varnames, data)
	PHI <- model.matrix(inst_varnames, data)
	X <- model.matrix(exog_varnames, data)

	coef$endg_var <- matrix(NA, ncol(D), length(taus))
	coef$inst_var <- matrix(NA, ncol(PHI), length(taus))
	coef$exog_var <- matrix(NA, ncol(X), length(taus))

	if (ncol(D) != ncol(PHI)){
		warning("Shouldn't ncol(D) always equal ncol(PHI)?")
	}
	# ---
	# Check the provided grid has appropriate dim
	if (ncol(D) == 1) {
		if (!is.vector(grid)) {
			stop("Dimension of the grid does not match the numbers of endogenous 
				variables. Grid should be a vector for dim(D) = 1")
	} else if (ncol(D) == 2) { # Grid search on a square when dim(D) == 2
		if (!is.matrix(grid) )
			stop("Dimension of the grid does not match the number of endogenous 
				variables. Grid should be a matrix with 2 rows for dim(D) == 2")
		} else if (ncol(D) >= 2) {
			stop("Complexity grows exponentially in the number of endogenous 
				variables. This version only deals with cases when dim(D) <= 2")
		}
	}
	# ---

	residuals <- matrix(0,nrow(X),length(taus))

	for (i in 1:length(taus)) {
		iqr_est <- ivqr.fit(iqr_formula, tau = taus[i], data, grid, gridMethod, 
			ivqrMethod, qrMethod)

		coef$endg_var[,i] <- iqr_est$coef_endg_var
		coef$inst_var[,i] <- iqr_est$coef_inst_var
		coef$exog_var[,i] <- iqr_est$coef_exog_var
		residuals[,i] <- iqr_est$residuals
	}
	
	# Preparing for standard erros
	fit <- list()
	class(fit) <- "ivqr"	
	fit$coef <- coef
	fit$residuals <- residuals
	fit$formula <- formula
	fit$taus <- taus
	fit$data <- data
	fit$dim_d_d_k <- c(ncol(D),ncol(D),ncol(X))
	PSI <- cbind(PHI, X)
	DX <- cbind(D,X)
	
	fit$DX <- DX
	fit$PSI <- PSI

	vc <- ivqr.vc(fit,covariance,"Silver")
	
	plot.ivqr(fit)
	fit$se <- vc$se
	return(fit)


	# return(list(coef=coef,vc=vc$covThetaMat,h=vc$bandwith_h,residual=residual,formula=formula(formula),taus=taus,
	#			data=data,ks=ks))

}

ivqr.fit <- function(iqr_formula, tau, data, grid, gridMethod, ivqrMethod, qrMethod){
	if (length(tau) > 1) {
		tau <- tau[1]
		warning("Multiple taus not allowed in rq.fit: solution restricted to first element")
	}
    fit <- switch(ivqrMethod,
		iqr = ivqr.fit.iqr(iqr_formula, tau, data, grid, gridMethod, qrMethod))
    #! should implement warining message about unimplemented method
	return(fit)
}

ivqr.fit.iqr <- function(iqr_formula, tau, data, grid, gridMethod, qrMethod){
	# Grid Search
	objFcn <- GetObjFcn(iqr_formula,tau,data,qrMethod)
	coef_endg_var <- GridSearch(objFcn,tau,grid)
	
	# reg again to get theta = (beta,gamma) and residuals
	inv_formula <- formula(iqr_formula,lhs=1,rhs=c(2,3),collapse=TRUE)
	response_varname <- as.character(inv_formula[[2]])

	Y <- c(data[[response_varname]])
	
	D <- formula(iqr_formula, lhs = 1, rhs = 1)
	D <- update(D,.~.-1)
	D <- model.matrix(D,data)
	
	X <- formula(iqr_formula, lhs = 1, rhs = 3)	
	X <- model.matrix(X, data)

	dim_inst_var <- formula(iqr_formula, lhs = 1, rhs = 2)
	dim_inst_var <- update(dim_inst_var,.~.-1)
	dim_inst_var <- ncol(model.matrix(dim_inst_var,data))

	data[[response_varname]] <- Y - D %*% coef_endg_var
	rq_fit <- rq(inv_formula,tau,data)
	coef_exog_var <- c(rq_fit$coef[1],rq_fit$coef[(2 + dim_inst_var) : length(rq_fit$coef)])
	coef_inst_var <- rq_fit$coef[2 : (2 + dim_inst_var - 1)] 
	residuals <- Y - D %*% coef_endg_var - X %*% coef_exog_var
	return(list(coef_endg_var = coef_endg_var, coef_exog_var = coef_exog_var,
		coef_inst_var = coef_inst_var, residuals = residuals))
}

GetObjFcn <- function(iqr_formula, taus, data, qrMethod) {
	
	inv_formula <- formula(iqr_formula,lhs=1,rhs=c(2,3),collapse=TRUE)
	response_varname <- as.character(inv_formula[[2]])

	Y <- c(data[[response_varname]])
	
	D <- formula(iqr_formula, lhs = 1, rhs = 1)
	D <- update(D,.~.-1)
	D <- model.matrix(D,data)

	dim_inst_var <- formula(iqr_formula, lhs = 1, rhs = 2)
	dim_inst_var <- update(dim_inst_var,.~.-1)
	dim_inst_var <- ncol(model.matrix(dim_inst_var,data))
	
	# Is not passing data, Y, D, response_varname Ok?
	objFcn <- function(tau, alpha){
		data[[response_varname]] <- Y - D %*% alpha
		rq_fit <- rq(inv_formula,tau,data,method=qrMethod)
		gamma <- rq_fit$coef[2 : (2 + dim_inst_var - 1)] 
		return(gamma)
	}
	return(objFcn)
}

GridSearch <- function(objFcn,tau,grid){ 
	if (is.vector(grid)) {
		gridValue <- rep(NA,length(grid))		
		for (i in 1:length(grid)) { 
			gridValue[i] <- objFcn(tau,grid[i])	
		}
		endg_var <- grid[which.min(abs(gridValue))]
		return(endg_var)	
	} else if (is.matrix(grid)) {
		gridValue <- marix(NA,length(grid[1,]),length(grid[2,]))				
		for (i in 1:length(grid[1])) {
			for (j in 1:length(grid[2])) {
				alpha <- c(grid[1,i],grid[2,i])
				gridValue[i] <- objFcn(tau,alpha)	
			}
		} 
		endg_var <- grid[which.min(abs(gridValue))]
		return(endg_var)
	}
}

ivqr.vc <- function(object,covariance,bd_rule="Silver") {
	taus <- object$taus
	DX <- object$DX
	PSI <- object$PSI
	residuals <- object$residuals
	tPSI_PSI <- t(PSI) %*% PSI # tPSI_PSI = sum_over_i(psi_i %*% t(psi_i))
	n <- nrow(PSI)
	kd <- ncol(PSI)

	se <- matrix(NA,kd,length(taus))
	cov_mats <- array(NA,dim = c(kd,kd,length(taus)))
	for( tau_index in 1:length(taus)){
		e <- residuals[,tau_index]

		# Silverman's rule of thumb
		if ( bd_rule == "Silver" ) { 
			h <- 1.364 * ( (2*sqrt(pi)) ^ (-1/5) ) * std(e) * ( n ^ (-1/5) )
		}

		S <- (taus[tau_index] - taus[tau_index] ^ 2) * (1 / n) * tPSI_PSI
		kernel_mat <- diag(as.numeric( abs(e) < h ))

		J <- (1 / (2 * n * h)) * t(kernel_mat %*% PSI) %*% DX
		invJ <- solve(J)

		cov_mats[,,tau_index] <- (1/n) * invJ %*% S %*% invJ
		se[,tau_index] <- diag(cov_mats[,,tau_index]) ^ (1/2)
	}
	
	vc <- list()
	class(vc) <- "ivqr.vc"
	vc$se <- se
	vc$cov_mats <- cov_mats
	return(vc)
}

print.ivqr <- function(x, ...) {
	d <- x$dim_d_d_k[1]
	k <- x$dim_d_d_k[3]
	taus <- x$taus
	coef <- x$coef
	se <- x$se

	endg_var <- cbind(coef$endg_var, se[1:d])
	exog_var <- cbind(coef$endg_var, se[(d + 1):(d + k)])

	# if ( d + k != length(x$se[,1])) {
	# 	warning("Dimension of se does not match that of ceof estimates")
	# }

	for (tau_index in 1:length(taus)){

		cat("\ntau:")
		print(taus[tau_index])
		cat("\nCoefficients of endogenous variables\n")
		table <- cbind(t(coef$endg_var)[tau_index,],se[1:d,tau_index])
		colnames(table) <- c("coef","se")
		print(table)
		cat("\nCoefficients of control variables\n")
		table <- cbind(t(coef$exog_var)[tau_index,],se[(d + 1):(d + k),tau_index])
		colnames(table) <- c("coef","se")
		print(table)			

	}

}

plot.ivqr(object){
	warning("plot.ivqr only implements for dim(D) == 1")
	up_bdd <- object$coef$endg_var + 1.96 * object$se[1,]
	lw_bdd <- object$coef$endg_var - 1.96 * object$se[1,]
	plot(taus,object$coef$endg_var,ylim=c(0,30000),type='n')
	polygon(c(taus,rev(taus)),c(up_bdd,rev(lw_bdd)),col="grey")
	lines(taus,object$coef$endg_var,ylim=c(0,30000))
}
