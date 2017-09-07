ivqr <- function(formula, taus=0.5, data, grid="Default"){
	
	formula <- Formula(formula) 

	# --- 
	# Check if taus has appropriate range
	eps <- .Machine$double.eps ^ (2/3)
    if(length(taus) > 1){
		if(any(taus < 0) || any(taus > 1)) 
		stop("invalid taus:  taus should be >= 0 and <= 1")
	}
	if(any(taus == 0)) taus[taus == 0] <- eps
	if(any(taus == 1)) taus[taus == 1] <- 1 - eps

	# ---
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
	if (any(grepl("ivqr_dhat",colnames(data)))){
		stop("No names of variables in the data set should include ivqr_dhat")
	}
	dhat_formula <- c("~")
	for (i in 1:dim(D)[2]){
		data[,paste("ivqr_dhat",i,sep="")] <- D_hat[,i]
		dhat_formula <- cbind(dhat_formula,paste("ivqr_dhat",i,sep=""))
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


	coef$endg_var <-  
	coef$exog_var <- 
	coef$inst_var <-
	iqr.fit(formula, taus, data, grid, gridMethod, maxIter)

	iqr.fit(formula, taus, data, grid, gridMethod, maxIter)
	
	residual <- GetResidual(formula, taus, data, coef)
	vc <- iqr.vc(formula, taus, data, coef, residual)
	ks <- iqr.ks(formula, taus, data, coef, residual)
	
	return(list(coef=coef,vc=vc$covThetaMat,h=vc$bandwith_h,residual=residual,formula=formula(formula),taus=taus,
				data=data))#,ks=ks))

}

ivqr.fit <- function(formula, tau, data, grid, gridMethod, maxIter){

	if(length(tau) > 1) {
		    tau <- tau[1]
		    warning("Multiple taus not allowed in rq.fit: solution restricted to first element")
		}

	dimEndgsVar <- length(attr(terms(formula(formula, lhs=0, rhs=c(1))),"term.labels"))
	dimCtrlVar <- 1 + length(attr(terms(formula(formula, lhs=0, rhs=c(3))),"term.labels"))
	
	estAlpha <- array(NA,dim=c(length(taus),dimEndgsVar))
	estGamma <- array(NA,dim=c(length(taus),dimEndgsVar)) 
	# We project treatment variable on instrument var so they have the same dimension

	estBeta <- array(NA,dim=c(length(taus),dimCtrlVar))

	# Project D on Z to get instruments
	# replace Z with the projected value of D
	data <- GetInstrument(formula, data)

	for (tau_index in 1:length(taus)) {
		print(paste("tau index is",tau_index))
		tau <- taus[tau_index]
		if(gridMethod == "Default"){
			initialGrid <- GetInitialGrid(formula, tau, data)
			value <- InvStep(formula,tau,data,initialGrid,gridMethod)
		} 

		if(gridMethod == "2SQR"){
			initialGrid <- GetInitialGrid2SQR(formula, tau, data)
			value <- InvStep(formula,tau,data,initialGrid,gridMethod)
		} 

		if(gridMethod == "Assigned"){
			value <- InvStep(formula,tau,data,grid,gridMethod)
		} 
		estAlpha[tau_index,] <- value$alpha
		estGamma[tau_index,] <- value$coef[2:(2+dimEndgsVar-1)]
		estBeta[tau_index,] <- c(value$coef[1],value$coef[(2+dimEndgsVar):(2+dimEndgsVar+dimCtrlVar-1-1)])		
	}
	

	# qr for the theta_hat
	

	
	# GetAsySE(formula, taus=0.5, data, residual)

	# format
	# S3 object

	return(list(alpha=estAlpha,beta=estBeta,gamma=estGamma))
}