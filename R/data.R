#' An example data created by simulation.
#'
#' A dataset illustrating the usage of the IVQR pacakge.
#'
#' @format A data frame with 10000 rows and 4 variables:
#' \describe{
#'   \item{y}{outcome variable}
#'   \item{d}{binary endogenous variable}
#' 	 \item{z}{binary instrumental variable}
#' 	 \item{x}{control variable}
#'   ...
#' }
#' @examples 
#' # The data generation process
#' sim_ivqr_eg <- function(n = 10 ^ 4){
#' 	u <- runif(n)
#' 	x <- rbinom(n, 1,0.2)
#' 	z <- rbinom(n, 1, 0.37)

#' 	v <- rnorm(n)
#' 	d <- z * (u > 0.5 * v)
#' 	y0 <- 0 + x * 2 + qnorm(u,0,1)
#' 	y1 <- (u - 0.5) + x * 2 + qnorm(u,0,1)
#' 	y <- d * y1 + (1 - d) * y0
	
#' 	value <- list()
#' 	value$y <- y
#' 	value$d <- d
#' 	value$z <- z
#' 	value$x <- x
#' 	value <- data.frame(value)
#' 	return(value)
#' }
"ivqr_eg"

sim_ivqr_eg <- function(n = 10 ^ 4){
	u <- runif(n)
	x <- rbinom(n, 1,0.2)
	z <- rbinom(n, 1, 0.37)

	v <- rnorm(n)
	d <- z * (u > 0.5 * v)
	y0 <- 0 + x * 2 + qnorm(u,0,1)
	y1 <- (u - 0.5) + x * 2 + qnorm(u,0,1)
	y <- d * y1 + (1 - d) * y0
	
	value <- list()
	value$y <- y
	value$d <- d
	value$z <- z
	value$x <- x
	value <- data.frame(value)
	return(value)
}

ivqr_eg <- sim_ivqr_eg(10 ^ 4)
devtools::use_data(ivqr_eg,overwrite=TRUE)
head(ivqr_eg)