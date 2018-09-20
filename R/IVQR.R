#' Instrumental Variable Qauntile Regression
#' @param fomula A Formula (not the built-in formula) object. A correct model 
#' specificaion should look like y ~ d | z | x, where y is the outcome variable,
#' d is the endogenous variable, z is the instrument, and x is the control 
#' variable'. If no control is needed, write y ~ d | z | 1. Constant has to be
#' included, and number of endogenous variable should be no more than two.
#' @param taus The quantile(s) to be estimated, which should be a vector that 
#' contains number strictly between 0 and 1.
#' @param data A data.frame in which to interpret the variables named in the 
#' formula.
#' @param grid A vector (resp. list) when number of endogenous variable is one 
#' (resp. two). 
#' The grid should be reasonably large 
#' and fine so that the objective function is properly minimized. 
#' When there are two endogenous variables, the list should contain two vectors,
#' each specifies a grid for one of the variables.
#' Users are encouraged to try different grids and use the function 
#' diagnose() to visually inspect the objective functions.
#' @param gridMethod The type of grid search. In the current version,
#' only the method "Default" is available, which is simply an exhaustive 
#' searching.
#' @param ivqrMethod The algorithm to compute the fit. In the current version,
#' only the "iqr" (instrumental quantile regression) algorithm is available.
#' @param qrMethod The algorithmic method used in the (conventional) quantile
#' regression, which is part of the instrumental quantile regression algorithm.
#' For more information, see the help file of function \code{rq()} from the 
#' package quantreg.
#' @return An ivqr object. Use functions such as summary() or plot() to see
#' the coefficient estimates and standard errors. Weak-IV robust inference and 
#' general inferences are implemented respectively by the fucntions 
#' \code{weakIVtest()} and \code{ivqr.ks()}.
#' @examples
#' data(ivqr_eg)
#' fit <- ivqr(y ~ d | z | x, 0.5, grid = seq(-2,2,0.2), data = ivqr_eg) # median
#' fit <- ivqr(y ~ d | z | x, 0.25, grid = seq(-2,2,0.2), data = ivqr_eg) # the first quartile
#' fit <- ivqr(y ~ d | z | x, seq(0.1,0.9,0.1), grid = seq(-2,2,0.2), data = ivqr_eg) # a process
#' plot(fit) # plot the coefficients of the endogenous variable along with 95% CI
#' 
#' data(ivqr_eg2) # An example with two endogenous variables
#' grid2 <- list(seq(-1,1,0.05),seq(1.5,2.5,0.05)) # Two grids
#' fit <- ivqr(y ~ d1 + d2 | z1 + z2 | x, seq(0.1,0.9,0.1), grid = grid2, data = ivqr_eg2) # median
#' plot(fit,1) # plot the coefficients of the first endogenous variable along with 95% CI
#' plot(fit,2) # plot the coefficients of the second endogenous variable along with 95% CI
#' summary(fit) # print the results
#' @import Formula
#' @export
ivqr <- function(formula, taus=0.5, data, grid, gridMethod="Default", ivqrMethod="iqr",
				qrMethod="br"){

	formula <- Formula(formula)
	if(length(formula)[2] < 3){
		stop("If there's no control variable, specify the model like y ~ d | z | 1")
	}
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

	if (!any(grepl(qrMethod,c("br","fn","fnb")))){
		stop("Please specify one of the quantreg solution method: br, fn, or fnb")
	}
	# By default, we use D projection on (X,Z) as instrument
	# Add the instruments to the data and update the formula respectively
	XZ <- formula(formula, lhs = 1, rhs = c(2,3), collapse = TRUE)
	XZ <- model.matrix(XZ,data)
	D  <- update(formula(formula,lhs = 1, rhs = 1), . ~ . - 1)
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
			stop("Dimension of the grid does not match the numbers of endogenous variables. Grid should be a vector for dim(D) = 1")
		}
	} else if (ncol(D) == 2) { # Grid search on a square when dim(D) == 2
		if (!is.list(grid) )
			stop("Dimension of the grid does not match the number of endogenous variables. Grid should be a matrix with 2 rows for dim(D) == 2")
	} else if (ncol(D) >= 2) {
			stop("Complexity grows exponentially in the number of endogenous variables. This version only deals with cases when dim(D) <= 2")
	}
	# ---

	fitted <- matrix(NA,nrow(X),length(taus))
	residuals <- matrix(NA,nrow(X),length(taus))
	if (!is.list(grid)){
		grid_value <- matrix(NA,length(grid),length(taus))
	} else {
		grid_value <- array(NA, dim = c(length(grid[[1]]),length(grid[[2]]),length(taus)))
	}
	error_tau_flag <- !logical(length(taus))
	for (i in 1:length(taus)) {
		# print(paste("Now at tau=",taus[i]))
		ivqr_est <- ivqr.fit(iqr_formula, tau = taus[i], data, grid, gridMethod,
			ivqrMethod, qrMethod)
		if (is.list(ivqr_est)){
			coef$endg_var[,i] <- ivqr_est$coef_endg_var
			coef$inst_var[,i] <- ivqr_est$coef_inst_var
			coef$exog_var[,i] <- ivqr_est$coef_exog_var

			# Name labeling
			rownames(coef$endg_var) <- names(ivqr_est$coef_endg_var)
			rownames(coef$inst_var) <- names(ivqr_est$coef_inst_var)
			rownames(coef$exog_var) <- names(ivqr_est$coef_exog_var)

			residuals[,i] <- ivqr_est$residuals
			fitted[,i] <- ivqr_est$fitted
			
			if (!is.list(grid)){
				grid_value[,i] <- ivqr_est$grid_value
			} else {
				grid_value[,,i] <- ivqr_est$grid_value
			}
			error_tau_flag[i] <- FALSE
		}
	}

	# Column names
	taulabs <- paste("tau=",format(round(taus,3)))
	colnames(coef$endg_var) <- colnames(coef$inst_var) <- colnames(coef$exog_var) <- taulabs

	# Preparing for standard erros
	fit <- list()
	class(fit) <- "ivqr"
	fit$coef <- coef
	fit$fitted <- fitted
	fit$residuals <- residuals
	fit$formula <- formula
	fit$taus <- taus
	fit$error_tau_flag <- error_tau_flag
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

	# print("Estimating se")
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
	coef_endg_var <- grid_search$coef_endg_var
	grid_value <- grid_search$grid_value

	# reg again to get theta = (beta,gamma) and residuals
	inv_formula <- formula(iqr_formula,lhs=1,rhs=c(2,3),collapse=TRUE)
	response_varname <- as.character(inv_formula[[2]])

	Y <- c(data[[response_varname]])

	D <- formula(iqr_formula, lhs = 1, rhs = 1)
	D <- update(D,.~.-1)
	D <- model.matrix(D,data)
	names(coef_endg_var) <- colnames(D)

	X <- formula(iqr_formula, lhs = 1, rhs = 3)
	X <- model.matrix(X, data)

	dim_inst_var <- formula(iqr_formula, lhs = 1, rhs = 2)
	dim_inst_var <- update(dim_inst_var,.~.-1)
	dim_inst_var <- ncol(model.matrix(dim_inst_var,data))

	data[[response_varname]] <- Y - D %*% coef_endg_var

	fit_rq <-
	withCallingHandlers(
	  	tryCatch(quantreg::rq(inv_formula,tau,data,method=qrMethod),
		  	error = function(e){
		  		cat("ERROR :",conditionMessage(e), "\n")
		  		#! Why X singular or not can depend on y & tau?
		  		cat(paste("Singular Design Matrix for whole grid", tau))
		  		cat(paste("This occurs at tau=", tau))
		  		return(Inf)
		  		}
		  	),
	  	warning = Suppress_Sol_not_Unique
  	)

	if (!grepl("rq",class(fit_rq))) return(Inf)

	#coef_exog_var <- c(fit_rq$coef[1])
	if(length(fit_rq$coef) >= (2 + dim_inst_var)){
		coef_exog_var <- c(fit_rq$coef[1],fit_rq$coef[(2 + dim_inst_var) : length(fit_rq$coef)])	
	} else {
		coef_exog_var <- c(fit_rq$coef[1])
	}
	
	coef_inst_var <- fit_rq$coef[2 : (2 + dim_inst_var - 1)]
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
		  rq_fit <- withCallingHandlers(
			  	tryCatch(quantreg::rq(inv_formula,tau,data,method=qrMethod),
			  	error = function(e){
			  		cat("ERROR :",conditionMessage(e), "\n")
			  		cat(paste("This occurs at tau=", tau,
			  			"when evaluating grid value alpha=",alpha),"\n")
			  		return(Inf)}),
			  	warning = Suppress_Sol_not_Unique
		  )

		if (!grepl("rq", class(rq_fit))) return(Inf)
		gamma <- rq_fit$coef[2 : (2 + dim_inst_var - 1)]
		cov_mat <- quantreg::summary.rq(rq_fit, se = "ker", covariance = TRUE)$cov
		cov_mat <- cov_mat[(2 : (2 + dim_inst_var - 1)),(2 : (2 + dim_inst_var - 1))]
		wald_stat <- t(gamma) %*% solve(cov_mat) %*% gamma
		return(wald_stat)
		#return(gamma)
	}
	return(objFcn)
}

GridSearch <- function(objFcn,tau,grid){
	if (!is.list(grid)) {
		grid_value <- rep(NA,length(grid))
		for (i in 1:length(grid)) {
			grid_value[i] <- objFcn(tau,grid[i])
		}
		output <- list()
		output$coef_endg_var <- grid[which.min(abs(grid_value))]
		output$grid_value <- grid_value
		return(output)
	} else if (is.list(grid)) {
		grid_value <- matrix(NA,length(grid[[1]]),length(grid[[2]]))
		for (i in 1:length(grid[[1]])) {
			for (j in 1:length(grid[[2]])) {
				alpha <- c(grid[[1]][i],grid[[2]][j])
				grid_value[i,j] <- objFcn(tau,alpha)
  			}
  		}
		idx <- which.min(abs(grid_value))
		row_idx <- idx %% length(grid[[1]])
		if (row_idx == 0) row_idx <- length(grid[[1]])	
		col_idx <- (idx - row_idx) / length(grid[[1]]) + 1

		output <- list()
		output$coef_endg_var <- c(grid[[1]][row_idx],grid[[2]][col_idx])
		output$grid_value <- grid_value
		return(output)
		
	}
}

ivqr.vc <- function(object, covariance, bd_rule="Silver", h_multi = 1) {
	taus <- object$taus
	error_tau_flag <- object$error_tau_flag
	DX <- object$DX
	PSI <- object$PSI
	residuals <- object$residuals
	tPSI_PSI <- t(PSI) %*% PSI # tPSI_PSI = sum_over_i(psi_i %*% t(psi_i))
	n <- nrow(PSI)
	kd <- ncol(PSI)

	se <- matrix(NA,kd,length(taus))
	cov_mats <- array(NA,dim = c(kd,kd,length(taus)))
	J_array <- array(NA,dim = c(kd,kd,length(taus)))

	Check_Invertible <- function(m) class(try(solve(m),silent=T))=="matrix"

	for (tau_index in 1:length(taus)){
		while (error_tau_flag[tau_index] & tau_index < length(taus)){
			tau_index <- tau_index + 1
		}

		if (error_tau_flag[tau_index] & tau_index == length(tau_index)){
			break
		}

		e <- residuals[,tau_index]
		# Silverman's rule of thumb

		if (bd_rule == "Silver") {
			h <- h_multi * 1.364 * ( (2*sqrt(pi)) ^ (-1/5) ) * sd(e) * ( n ^ (-1/5) )
		}

		S <- (taus[tau_index] - taus[tau_index] ^ 2) * (1 / n) * tPSI_PSI


		kernel <- c(as.numeric( abs(e) < h ))
		J <- (1 / (2 * n * h)) * t(kernel * PSI) %*% DX

		if (!Check_Invertible(J)) warning(paste("At tau=",taus[tau_index],":",
			"Matrix J in Powell kernel estimator is not invertible with default bandwith"))

		while (!Check_Invertible(J)){
			h <- h * 1.1
			kernel <- c(as.numeric( abs(e) < h ))
			J <- (1 / (2 * n * h)) * t(kernel * PSI) %*% DX
		}

		J_array[,,tau_index] <- J
		invJ <- solve(J)

		cov_mats[,,tau_index] <- (1/n) * invJ %*% S %*% t(invJ)
		se[,tau_index] <- diag(cov_mats[,,tau_index]) ^ (1/2)
	}

	vc <- list()
	class(vc) <- "ivqr.vc"
	vc$se <- se
	vc$cov_mats <- cov_mats
	vc$J <- J_array
	return(vc)
}

#' @export
print.ivqr <- function(x, ...){
	cat("\nCoefficients of endogenous variables:\n\n")
	print(x$coef$endg_var)
	cat("\nCoefficients of exogenous variables:\n\n")
	print(x$coef$exog_var)
}

#' @export
plot.ivqr <- function(object, variable = 1, trim = c(0.05,0.95), ...){
	taus <- object$taus
	tl <- which(taus >= trim[1])[1]
	th <- which(taus <= trim[2])[length(which(taus < trim[2]))]
	taus <- taus[tl:th]

	coef <- object$coef$endg_var[variable,tl:th]
	yname <- rownames(object$coef$endg_var)[variable]

	se <- object$se[1,tl:th]
	up_bdd <- coef + 1.96 * se
	lw_bdd <- coef - 1.96 * se

	plot(taus,coef, ylim = c(min(lw_bdd) - max(se),
		max(up_bdd) + max(se)), type='n', xlab = "tau", ylab = yname)
	polygon(c(taus, rev(taus)), c(up_bdd, rev(lw_bdd)), col = 'grey')
	lines(taus,coef, ylim = c(min(lw_bdd) - max(se),
		max(up_bdd) + max(se)))
}

#' @export
summary.ivqr <- function(x, i = NULL, ...) {
	d <- x$dim_d_d_k[1]
	k <- x$dim_d_d_k[3]
	taus <- x$taus
	coef <- x$coef
	se <- x$se

	# if ( d + k != length(x$se[,1])) {
	# 	warning("Dimension of se does not match that of ceof estimates")
	# }

	if(!is.null(i)){
		start <- end <- i 
	} else{
		start <- 1
		end <- length(taus)
	}

	for (tau_index in start:end){
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

#' Weak-IV Robust Test
#' @param object An ivqr object returned from the function \code{ivqr()}
#' @param size the size of the test. Default is 0.05.
#' package quantreg.
#' @return An ivqr_weakIV object, which will be directly inputted into the method
#' \code{plot.ivqr_weakIV()} to visualize the confidence 
#' region. It is also possible to print the confidence region by applying the method 
#' \code{print()} to the ivqr_weakIV object. 
#' Note that the confidence region is not guaranteed to be convex.
#' @examples 
#' data(ivqr_eg)
#' fit <- ivqr(y ~ d | z | x, seq(0.1,0.9,0.1), grid = seq(-2,2,0.2), data = ivqr_eg) # a process
#' weakIVtest(fit)
#' print(fit) # print the confidence region.
#' @export
weakIVtest <- function(object, size = 0.05){
	grid <- object$grid
	dim_d <- object$dim_d_d_k[1]
	obj_fcn <- object$obj_fcn
	taus <- object$taus

	if (dim_d > 1) stop("weakIVtest() is only implented for single endogenous variable")
	critical_value <- qchisq((1 - size), dim_d)
	grid_mat <- replicate(length(taus),grid)
	grid_mat[obj_fcn > critical_value] <- NA
	result <- list()
	result$CI <- grid_mat
	result$taus <- taus
	result$yname <- rownames(object$coef$endg_var)[1]
	class(result) <- "ivqr_weakIV"
	warning("weakIVtest: CI can be not-convex. Also, it the dots hits the boundary, the grid used in ivqr() should be widened")

	plot.ivqr_weakIV(result)
	return(result)
}



#' Visualizing the Weak-IV Robust Confidence Region
#' @param object An ivqr_weakIV object returned from the function \code{weakIVtest()}
#' @details This function will be called when evaluating \code{weakIVtest()}. 
#' The case of two endogenous varaibles is not supported.
#' @examples 
#' data(ivqr_eg)
#' fit <- ivqr(y ~ d | z | x, seq(0.1,0.9,0.1), grid = seq(-2,2,0.2), data = ivqr_eg) # a process
#' plot(weakIVtest(fit))
#' @export
plot.ivqr_weakIV <- function(object, ...){	
	CI <- object$CI
	taus <- object$taus
	yname <- object$yname
	plot(rep(taus,nrow(CI)), c(t(CI)), xlab = "tau", ylab = yname)
}
#' Printing the Weak-IV Robust Confidence Region
#' @param x An ivqr_weakIV object returned from the function \code{weakIVtest()}
#' @details The case of two endogenous varaibles is not supported.
#' @export 
print.ivqr_weakIV <- function(x,...){
	cat("\nBelow prints the confidence region of the endogenous variable. NA means the corresponding value (defined in the user-specified grid when calling ivqr()) is not in the confidence region. Note that CI can be not-convex.\n")
	print(x$CI)
}


Suppress_Sol_not_Unique <-function(w) {
	if( any( grepl( "Solution may be nonunique", w) ) ) invokeRestart( "muffleWarning" )
}

#' Plotting the Objective Function
#' @param object An ivqr object returned from the function \code{ivqr()}
#' @param i index of specific quantile to inspect. For example, if
#' argument \code{taus = c(0.25,0.5,0.75)} was used as an input of \code{ivqr()}
#' , specify \code{i=2} for the median.
#' @param size the size of the test. Default is 0.05.
#' @param trim a vector of two number indicating the lower and upper bounds 
#' of the quantiles to inspect.
#' @return A plot of the objective function as well as the (asymptotical) first
#' order conditions evaluated at the coefficient estimates.
#' @examples 
#' data(ivqr_eg)
#' fit <- ivqr(y ~ d | z | x, c(0.25,0.5,0.75), grid = seq(-2,2,0.2), data = ivqr_eg)
#' diagnose(fit,2) # Plot the objective function at the median.
#' @export
diagnose <- function(object, i = 1, size = 0.05, trim = NULL, GMM_only = 0){	
	dim_d <- object$dim_d_d_k[1]
	if (dim_d > 1) stop("diagnose() is only implented for single endogenous variable")
	grid <- object$grid
	dim_d <- object$dim_d_d_k[1]
	obj_fcn <- object$obj_fcn[,i]
	taus <- object$taus
	PSI <- object$PSI
	n <- object$n
	residuals <- object$residuals
	if (!is.null(trim)){
		if (min(trim) > max(grid) | max(trim) < min(grid)){
			stop("Parameter trim results in nothing to plot. Either
				min(trim) > max(grid) or max(trim) < min(grid)")
		}
		gl <- which(grid >= trim[1])[1]
		gh <- which(grid <= trim[2])[length(which(grid <= trim[2]))]
	} else {
		gl <- 1
		gh <- length(grid)
	}
	critical_value <- qchisq((1 - size), dim_d)

	if(length(i) > 1) {
		warning("Multiple taus not allowed in diagnose: plot restricted to first element")
	}


	# Plot the result from grid search
	if (!GMM_only){
		plot(grid[gl:gh], obj_fcn[gl:gh], type = 'l', col = "blue",
			ylab = "Objective Function", xlab = "Grid", main = paste0("Objective function at the ", taus[i],"-th quantile"))
		abline(h = critical_value, col = "green")
		legend("topright", legend = "Weak-IV Critical Value", col = "green", lty=1:2, cex=0.8)	
	}
	
	# Calculate the GMM FOC
	L <- rep(taus[i],n) - as.numeric(residuals[,i] < 0)
	GMM_criterion <- rowSums(L * t(PSI))
	print(paste("At tau=", taus[i], "GMM FOC has values:"))
	print("(These values should be closed to zero asymptotically)")
	print(cbind(GMM_criterion/(sqrt(n))))
}	
