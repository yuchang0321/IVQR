#' General Inference
#' Testing general hypothesis regarding the qunatile process of the endogenous
#' variable.
#' @param object An ivqr object returned from the function \code{ivqr()}
#' @param variable A number indicates which endogenous variable to test. Since 
#' at most two endongenous variables can be included in the function \code{ivqr},
#' this argument should be either 1 or 2.
#' @param trim A vector of two numbers indicating the lower and upper bounds 
#' of the quantiles to consider.
#' @param B Number of sub-sampling in the bootstrap. Default is 2000.
#' @param nullH The hypothesis to be tested. The four options are: No_Effect, 
#' Dominance, Location_Shift, and Exogeneity.
#' @return An ivqr_ks object which contains information regarding test statistics,
#' critical value, sub-sampling block size, ...etc.
#' @examples 
#' data(ivqr_eg)
#' fit <- ivqr(y ~ d | z | x, seq(0.15,0.85,0.02), grid = seq(-2,2,0.2), data = ivqr_eg) # taus should be a fine grid 
#' ivqr(fit,nullH=No_Effect) # Test of no effect.
#' ivqr(fit,nullH=Dominance) # Test of dominance.
#' ivqr(fit,nullH=Location_Shift) # Test of location shift.
#' ivqr(fit,nullH=Exogeneity) # Test of exogeneity.
#' @export
ivqr.ks <- function(object, variable = NULL, trim = c(0.05,0.95), B = 2000,  b_scale = 1,
	nullH="No_Effect", ...){
	if (any(object$error_tau_flag)) Stop("Error occurred for some tau. Re-specify taus
		and run ivqr() again")
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
		Exogeneity = ivqr.ks.exog(object, trim, B, variable, b_scale),
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
	process <- process / sd_z[tl:th]
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
		invJ_m <- solve(J_m)
		R_invJ_m <- invJ_m[variable,]
		R_invJ <- invJ[variable,] # = R %*% invJ, R = c(0,0,,,0,1,0,..0)
		Z[,tau_index] <- (as.vector(L[,tau_index]) * (PSI  %*% R_invJ)
			- as.vector(L_m[,1]) * (PSI  %*% R_invJ_m))
	}

	sd_z <- colMeans(Z ^ 2) ^ (1 / 2)
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
	indexes <- c(tl:left_to_median,right_to_median:th)

	s <- (block_size) ^ (1 / 2) * apply(t(t(abs(V[,indexes])) / sd_z[indexes]),1,max)
	s <- as.vector(s)

	critical_value <- c(quantile(s,0.90),quantile(s,0.95),quantile(s,0.99))

	process <- rbind(coef$endg_var,coef$exog_var)[variable,indexes]

	ks_stat <- n ^ (1/2) * max(abs(process - coef_m) / sd_z[indexes])
	# print(process)
	# print(coef_m)
	# print(sd_z[indexes])
	# print(ks_stat)
	# stop()
	ks <- list()
	class(ks) <- "ivqr_ks"
	ks$ks_stat <- ks_stat
	ks$critical_value <- critical_value
	ks$block_size <- block_size
	ks$B <- B
	ks$s <- s

	return(ks)
}

ivqr.ks.exog <- function(object, trim, B, variable, b_scale, bd_rule = "Silver") {
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
	dim_dk <- object$dim_d_d_k[1] + object$dim_d_d_k[3]
	dim_d <- object$dim_d_d_k[1]
	n <- object$n
	PSI <- object$PSI
	J <- object$vc$J


	formula_rq <- formula(Formula(formula), lhs = 1, rhs = c(1,3), collapse = TRUE)
	fit_rq <- withCallingHandlers(
			  	rq(formula_rq, tau = taus, data = data, method = qrMethod),
			  	warning = Suppress_Sol_not_Unique
		  )

	coef_rq <- fit_rq$coef
	rq_residuals <- fit_rq$residuals
	model_matrix <- model.matrix(formula_rq,data)
	X <- cbind(model_matrix[,1],model_matrix[,(2 + dim_d):(dim(model_matrix)[2])])
	D <- model_matrix[,2:(1 + dim_d)]
	X_tilde <- cbind(D,X)

	# nrow = n, ncol = length(taus)
	A <- matrix(rep(taus,n), n, length(taus), byrow=TRUE) - as.numeric(rq_residuals < 0)
	L <- matrix(rep(taus,n), n, length(taus), byrow=TRUE) - as.numeric(residuals < 0)
	Z <- matrix(NA,n,length(taus))
	for (tau_index in 1:length(taus)) {

		e <- rq_residuals[,tau_index]
		if ( bd_rule == "Silver" ) {
			h <- 1.364 * ( (2*sqrt(pi)) ^ (-1/5) ) * std(e) * ( n ^ (-1/5) )
		}
		kernel <- c(as.numeric( abs(e) < h ))

		H <- (1 / (2 * n * h)) * t(kernel * X_tilde) %*% X_tilde
		invH <- solve(H)
		invJ <- solve(J[,,tau_index])
		R_invH <- invH[variable,]
		R_invJ <- invJ[variable,] # = R %*% invJ, R = c(0,0,,,0,1,0,..0)
		Z[,tau_index] <- (as.vector(L[,tau_index]) * (PSI  %*% R_invJ)
			- as.vector(A[,tau_index]) * (X_tilde  %*% R_invH))
	}

	sd_z <- colMeans(Z ^ 2) ^ (1 / 2)
	block_size <- as.integer(5 * n ^ (2/5)) * b_scale

	V <- matrix(NA, B, length(taus))
	for (b in 1:B) {
		resample_indexes <- sample(n, block_size, replace=FALSE)
		V[b,] <- sum(Z[resample_indexes,]) / block_size
	}

	# Trim the extreme tails cut out [1/2 - eps, 1/2 + eps]
	tl <- which(taus >= trim[1])[1]
	th <- which(taus <= trim[2])[length(which(taus < trim[2]))]

	s <- block_size ^ (1 / 2) * apply(t(t(abs(V[,(tl:th)])) / sd_z[tl:th]), 1, max)
	s <- as.vector(s)

	critical_value <- c(quantile(s,0.90),quantile(s,0.95),quantile(s,0.99))
	process <- rbind(
		coef$endg_var - coef_rq[2:(1 + dim_d),],
		coef$exog_var[1,] - coef_rq[1,],
		coef$exog_var[2:dim(coef$exog_var)[1],] - coef_rq[(2 + dim_d):(dim(coef_rq))[1],]
		)[variable,(tl:th)]

	ks_stat <- n ^ (1/2) * max(abs(process) / sd_z[(tl:th)])

	ks <- list()
	class(ks) <- "ivqr_ks"
	ks$ks_stat <- ks_stat
	ks$critical_value <- critical_value
	ks$block_size <- block_size
	ks$B <- B
	ks$s <- s

	return(ks)
}
#' @export
print.ivqr_ks <- function(x, ...) {
	print(x$ks_stat)
	print(x$critical_value)
	print(paste("Block size:",x$block_size))
}