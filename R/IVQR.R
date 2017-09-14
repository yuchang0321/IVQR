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
		beta <- solve(crossprod(XZ), crossprod(XZ,D[,i]))
		D_hat[,i] <- XZ %*% beta
	}

	#! ---
	if (any(grepl(".ivqr_dhat",colnames(data)))){
		stop("No names of variables in the data set should include .ivqr_dhat")
	}
	dhat_formula <- c("~")

	copy_data <- data
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

	fitted <- matrix(NA,nrow(X),length(taus))
	residuals <- matrix(NA,nrow(X),length(taus))
	grid_value <- matrix(NA,length(grid),length(taus))
	for (i in 1:length(taus)) {
		print(paste("Now at tau=",taus[i]))
		ivqr_est <- ivqr.fit(iqr_formula, tau = taus[i], data, grid, gridMethod,
			ivqrMethod, qrMethod)

		if (is.na(ivqr_est)){

			if (i == 1 & length(taus) > 1){
				new_taus <- taus[2:length(taus)]
			} else if ( i < length(taus)){
				new_taus <- c(taus[1:(i-1)],taus[(i+1):length(taus)])
			} else if (i == length(taus)){
				new_taus <- taus[1:(length(taus) - 1)]
			}
			return(ivqr(formula = formula(formula), taus= new_taus,
				data = copy_data, grid = grid, gridMethod = gridMethod,
				ivqrMethod = ivqrMethod, qrMethod = qrMethod))
		}

		coef$endg_var[,i] <- ivqr_est$coef_endg_var
		coef$inst_var[,i] <- ivqr_est$coef_inst_var
		coef$exog_var[,i] <- ivqr_est$coef_exog_var
		residuals[,i] <- ivqr_est$residuals
		fitted[,i] <- ivqr_est$fitted
		grid_value[,i] <- ivqr_est$grid_value

	}
	# Preparing for standard erros
	fit <- list()
	class(fit) <- "ivqr"
	fit$coef <- coef
	fit$fitted <- fitted
	fit$residuals <- residuals
	fit$formula <- formula
	fit$taus <- taus
	fit$data <- data
	fit$dim_d_d_k <- c(ncol(D),ncol(D),ncol(X))
	fit$n <- nrow(X)
	fit$obj_fcn <- grid_value
	fit$grid <- grid
	fit$copy_data <- copy_data
	fit$gridMethod <- gridMethod
	fit$ivqrMethod <- ivqrMethod
	fit$qrMethod <- qrMethod
	PSI <- cbind(PHI, X)
	DX <- cbind(D,X)

	fit$DX <- DX
	fit$PSI <- PSI

	print("Estimating se")
	vc <- ivqr.vc(fit,covariance,"Silver")

	#plot.ivqr(fit)
	fit$se <- vc$se
	fit$vc <- vc

	return(fit)

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
	grid_search <- GridSearch(objFcn,tau,grid)
	if (is.na(grid_search)) return(NA)
	coef_endg_var <- grid_search$coef_endg_var
	grid_value <- grid_search$grid_value

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
	fitted <- D %*% coef_endg_var + X %*% coef_exog_var
	residuals <- Y - fitted

	return(list(coef_endg_var = coef_endg_var, coef_exog_var = coef_exog_var,
		coef_inst_var = coef_inst_var, fitted = fitted, residuals = residuals,
		grid_value = grid_value))
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
		  rq_fit <- tryCatch({
			rq(inv_formula,tau,data,method=qrMethod)
		  },
		  error=function(e){
		  	cat("ERROR :",conditionMessage(e), "\n")
		  	return(NA)})
		if (is.na(rq_fit)) return(NA)
		gamma <- rq_fit$coef[2 : (2 + dim_inst_var - 1)]
		cov_mat <- summary.rq(rq_fit, se = "ker", covariance = TRUE)$cov
		cov_mat <- cov_mat[(2 : (2 + dim_inst_var - 1)),(2 : (2 + dim_inst_var - 1))]
		wald_stat <- t(gamma) %*% solve(cov_mat) %*% gamma
		return(wald_stat)
		#return(gamma)
	}
	return(objFcn)
}

GridSearch <- function(objFcn,tau,grid){
	if (is.vector(grid)) {
		grid_value <- rep(NA,length(grid))
		for (i in 1:length(grid)) {
			grid_value[i] <- objFcn(tau,grid[i])
			if (is.na(grid_value[i])) {
				print(paste("Error happened at tau=", tau))
				return(NA)
			}

		}
		output <- list()
		output$coef_endg_var <- grid[which.min(abs(grid_value))]
		output$grid_value <- grid_value
		return(output)
	} else if (is.matrix(grid)) {
		grid_value <- marix(NA,length(grid[1,]),length(grid[2,]))
		for (i in 1:length(grid[1])) {
			for (j in 1:length(grid[2])) {
				alpha <- c(grid[1,i],grid[2,i])
				grid_value[i] <- objFcn(tau,alpha)
  			}
			endg_var <- grid[which.min(abs(grid_value))]
			return(endg_var)
		}
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
	J_array <- array(NA,dim = c(kd,kd,length(taus)))
	print("Done with Preparing")
	for(tau_index in 1:length(taus)){
		e <- residuals[,tau_index]
		# Silverman's rule of thumb

		if ( bd_rule == "Silver" ) {
			h <- 1.364 * ( (2*sqrt(pi)) ^ (-1/5) ) * std(e) * ( n ^ (-1/5) )
		}

		S <- (taus[tau_index] - taus[tau_index] ^ 2) * (1 / n) * tPSI_PSI
		print("Done with S")
		kernel <- c(as.numeric( abs(e) < h ))
		print("Done with kernel")

		J <- (1 / (2 * n * h)) * t(kernel * PSI) %*% DX
		J_array[,,tau_index] <- J
		print("Done with J")
		#invJ <- chol2inv(chol(J))
		invJ <- solve(J)
		print("Done with invJ")


		cov_mats[,,tau_index] <- (1/n) * invJ %*% S %*% invJ
		print("Done with multiplication")
		se[,tau_index] <- diag(cov_mats[,,tau_index]) ^ (1/2)
	}

	vc <- list()
	class(vc) <- "ivqr.vc"
	vc$se <- se
	vc$cov_mats <- cov_mats
	vc$J <- J_array
	return(vc)
}

ivqr.ks <- function(object, variable = NULL, trim = c(0.05,0.95), B = 2000,  b_scale = 1,
	nullH="No_Effect"){

	dim_d <- object$dim_d_d_k[1]

	if (dim_d > 1 & is.null(variable)) {
		warning("Please specify which (endogenous) variable to test. Default
			is the first endogenous variable.")
	}

	if (dim_d == 1 & is.null(variable)) variable = 1

	ks <- switch(nullH,
		No_Effect = ivqr.ks.no(object, trim, B, variable, b_scale),
		Dominance = ivqr.ks.dom(object, trim, B, variable, b_scale),
		Location_Shift = ivqr.ks.const(object, trim, B, variable, b_scale),
		"This test is not implemented")

	return(ks)
}

ivqr.ks.no <- function(object, trim, B, variable, b_scale) {
	taus <- object$taus
	data <- object$data
	coef <- object$coef
	fitted <- object$fitted
	residuals <- object$residuals
	endg_var_se <- object$se[variable]
	dim_dk <- object$dim_d_d_k[1] + object$dim_d_d_k[3]
	n <- object$n
	PSI <- object$PSI
	J <- object$vc$J

	# nrow = n, ncol = length(taus)
	L <- matrix(rep(taus,n), n, length(taus), byrow = TRUE) - as.numeric(residuals < 0)

	Z <- matrix(NA,n,length(taus))

	for (tau_index in 1:length(taus)) {
		invJ <- solve(J[,,tau_index])
		R_invJ <- invJ[variable,] # = R %*% invJ, R = c(0,0,,,0,1,0,..0)
		Z[,tau_index] <- as.vector(L[,tau_index]) * (PSI %*% R_invJ)
	}

	sd_z <- colMeans(Z ^ 2) ^ (1 / 2)

	block_size <- as.integer(5 * n ^ (2/5)) * b_scale
	V <- matrix(NA, B, length(taus))
	for (b in 1:B) {
		resample_indexes <- sample(n, block_size, replace=FALSE)
		V[b,] <- sum(Z[resample_indexes,]) / block_size
	}

	# Trim the extreme tails
	tl <- which(taus >= trim[1])[1]
	th <- which(taus <= trim[2])[length(which(taus <= trim[2]))]
	s <- block_size ^ (1 / 2) * apply(t(t(abs(V[,(tl:th)])) / sd_z[tl:th]), 1, max)
	s <- as.vector(s)
	critical_value <- c(quantile(s,0.90),quantile(s,0.95),quantile(s,0.99))

	process <- rbind(coef$endg_var,coef$exog_var)[variable,(tl:th)]
	process <- process / sd_z
	ks_stat <- max(abs(process)) * n ^ (1/2)

	ks <- list()
	class(ks) <- "ivqr_ks"
	ks$ks_stat <- ks_stat
	ks$critical_value <- critical_value
	ks$block_size <- block_size
	ks$B <- B
	ks$s <- s

	return(ks)
}

ivqr.ks.dom <- function(object, trim, B, variable, b_scale) {
	taus <- object$taus
	data <- object$data
	coef <- object$coef
	fitted <- object$fitted
	residuals <- object$residuals
	endg_var_se <- object$se[variable]
	dim_dk <- object$dim_d_d_k[1] + object$dim_d_d_k[3]
	n <- object$n
	PSI <- object$PSI
	J <- object$vc$J

	# nrow = n, ncol = length(taus)
	L <- matrix(rep(taus,n), n, length(taus), byrow=TRUE) - as.numeric(residuals < 0)

	Z <- matrix(NA,n,length(taus))

	for (tau_index in 1:length(taus)) {
		invJ <- solve(J[,,tau_index])
		R_invJ <- invJ[variable,] # = R %*% invJ, R = c(0,0,,,0,1,0,..0)
		Z[,tau_index] <- as.vector(L[,tau_index]) * (PSI %*% R_invJ)
	}

	sd_z <- colMeans(Z ^ 2) ^ (1 / 2)

	block_size <- as.integer(5 * n ^ (2/5)) * b_scale
	V <- matrix(NA, B, length(taus))
	for (b in 1:B) {
		resample_indexes <- sample(n, block_size, replace=FALSE)
		V[b,] <- sum(Z[resample_indexes,]) / block_size
	}

	# Trim the extreme tails
	tl <- which(taus >= trim[1])[1]
	th <- which(taus <= trim[2])[length(which(taus < trim[2]))]

	s <- block_size ^ (1 / 2) * apply(t(t(abs(V[,(tl:th)])) / sd_z[tl:th]), 1, max)
	s <- as.vector(s)
	critical_value <- c(quantile(s,0.90),quantile(s,0.95),quantile(s,0.99))

	process <- rbind(coef$endg_var,coef$exog_var)[variable,(tl:th)]
	ks_stat <- max(max(-1 * process / sd_z[tl:th]), 0) * n ^ (1/2)

	ks <- list()
	class(ks) <- "ivqr_ks"
	ks$ks_stat <- ks_stat
	ks$critical_value <- critical_value
	ks$block_size <- block_size
	ks$B <- B
	ks$s <- s

	return(ks)
}


ivqr.ks.const <- function(object, trim, B, variable, b_scale) {
	taus <- object$taus
	data <- object$copy_data
	coef <- object$coef
	grid <- object$grid
	fitted <- object$fitted
	gridMethod <- object$gridMethod
	ivqrMethod <- object$ivqrMethod
	qrMethod <- object$qrMethod
	residuals <- object$residuals
	formula <- formula(object$formula)
	# endg_var_se <- object$se[variable]
	dim_dk <- object$dim_d_d_k[1] + object$dim_d_d_k[3]
	n <- object$n
	PSI <- object$PSI
	J <- object$vc$J

	ivqr_median_fit <- ivqr(formula, taus = 0.5, data = data , grid = grid,
		gridMethod = gridMethod, ivqrMethod = ivqrMethod, qrMethod = qrMethod)

	coef_m <- ivqr_median_fit$coef
	coef_m <- rbind(coef_m$endg_var, coef_m$exog_var)[variable,1]
	residuals_m <- ivqr_median_fit$residuals[,1]
	J_m <- ivqr_median_fit$vc$J[,,1]
	L_m <- matrix(rep(0.5,n), n, length(0.5), byrow=TRUE) - as.numeric(residuals_m < 0)

	# nrow = n, ncol = length(taus)
	L <- matrix(rep(taus,n), n, length(taus), byrow=TRUE) - as.numeric(residuals < 0)

	Z <- matrix(NA,n,length(taus))

	for (tau_index in 1:length(taus)) {
		invJ <- solve(J[,,tau_index])
		R_invJ <- invJ[variable,] # = R %*% invJ, R = c(0,0,,,0,1,0,..0)
		Z[,tau_index] <- (as.vector(L[,tau_index]) - as.vector(L_m[,1])) %*% (PSI  %*% R_invJ)
	}

	sd_z <- colMeans(Z ^ 2) ^ (1 / 2)
	print(sd_z)
	stop()
	block_size <- as.integer(5 * n ^ (2/5)) * b_scale
	V <- matrix(NA, B, length(taus))
	for (b in 1:B) {
		resample_indexes <- sample(n, block_size, replace=FALSE)
		V[b,] <- sum(Z[resample_indexes,]) / block_size
	}

	# Trim the extreme tails cut out [1/2 - eps, 1/2 + eps]
	tl <- which(taus >= trim[1])[1]
	th <- which(taus <= trim[2])[length(which(taus < trim[2]))]
	left_to_median <- which(taus < 0.5)[length(which(taus < 0.5 ))]
	right_to_median <- which(taus > 0.5)[1]

	s <- (block_size) ^ (1 / 2) * apply(abs(
		V[,c(tl:left_to_median,right_to_median:th)]),1,max)
	s <- as.vector(s)
	critical_value <- c(quantile(s,0.90),quantile(s,0.95),quantile(s,0.99))

	process <- rbind(coef$endg_var,coef$exog_var)[variable,c(tl:left_to_median,
		right_to_median:th)]
	ks_stat <- max(abs(process - coef_m )) * n ^ (1/2)

	ks <- list()
	class(ks) <- "ivqr_ks"
	ks$ks_stat <- ks_stat
	ks$critical_value <- critical_value
	ks$block_size <- block_size
	ks$B <- B
	ks$s <- s

	return(ks)
}


print.ivqr <- function(x, ...) {
	d <- x$dim_d_d_k[1]
	k <- x$dim_d_d_k[3]
	taus <- x$taus
	coef <- x$coef
	se <- x$se

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

weakIVtest <- function(object, size = 0.05){

	grid <- object$grid
	dim_d <- object$dim_d_d_k[1]
	obj_fcn <- object$obj_fcn
	taus <- object$taus

	#!
	if (dim_d > 1) stop("weakIVtest() is only implented for single endogenous variable")
	critical_value <- qchisq((1 - size), dim_d)
	grid_mat <- replicate(length(taus),grid)
	grid_mat[obj_fcn > critical_value] <- NA
	result <- list()
	result$CI <- grid_mat
	result$taus <- taus
	class(result) <- "ivqr_weakIV"
	return(result)
}

print.ivqr_ks <- function(x, ...) {
	print(x$ks_stat)
	print(x$critical_value)
	print(x$block_size)
}

plot.ivqr <- function(object, trim = c(0.05,0.95), variable = 1){
	warning("plot.ivqr is only implemented for dim(D) == 1")
	taus <- object$taus
	tl <- which(taus >= trim[1])[1]
	th <- which(taus <= trim[2])[length(which(taus < trim[2]))]
	taus <- taus[tl:th]

	coef <- object$coef$endg_var[1,tl:th]
	se <- object$se[1,tl:th]
	up_bdd <- coef + 1.96 * se
	lw_bdd <- coef - 1.96 * se

	plot(taus,coef, ylim = c(min(lw_bdd) - max(se),
		max(up_bdd) + max(se)), type='n')
	polygon(c(taus, rev(taus)), c(up_bdd, rev(lw_bdd)), col = 'grey')
	lines(taus,coef, ylim = c(min(lw_bdd) - max(se),
		max(up_bdd) + max(se)))
}


plot.ivqr_weakIV <- function(object){
	warning("weakIVtest: CI can be not-convex")
	CI <- object$CI
	taus <- object$taus
	plot(rep(taus,nrow(CI)),c(t(CI)))
}
