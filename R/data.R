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

#' An (second) example data created by simulation.
#'
#' A dataset with 2 endogenous variables to illustrate the usage of the IVQR pacakge.
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
#' sim_ivqr_eg2 <- function(n = 10 ^ 4){
#'   u <- runif(n)
#'   x <- rbinom(n, 1,0.2)

#'   z1 <- rbinom(n, 1, 0.37)
#'   v1 <- rnorm(n)
#'   d1 <- z1 * (u > 0.5 * v1)

#'   z2 <- rbinom(n, 1, 0.37)
#'   v2 <- rnorm(n)
#'   d2 <- z2 * (u > 0.5 * v1)

#'   y00 <- 0 + x * 2 + qnorm(u,0,1)
#'   y10 <- u  + x * 2 + qnorm(u,0,1)
#'   y01 <- 2 + x * 2 + qnorm(u,0,1)
#'   y11 <- 2 + u + x * 2 + qnorm(u,0,1)
#'   y <- d1 * d2 * y11  + (1 - d1) * d2 * y01 + d1 * (1 - d2) * y10 + (1 - d1) * (1 - d2) * y00

#'   value <- list()
#'   value$y <- y
#'   value$d1 <- d1
#'   value$d2 <- d2
#'   value$z1 <- z1
#'   value$z2 <- z2
#'   value$x <- x
#'   value <- data.frame(value)
#'   return(value)
#' }
"ivqr_eg2"

# For my own reference
# ivqr_eg <- gen_ivqr_eg(10 ^ 4)
# devtools::use_data(ivqr_eg,overwrite=TRUE)

gen_ivqr_eg <- function(n = 10 ^ 4){
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

gen_ivqr_eg2 <- function(n = 10 ^ 4){
  u <- runif(n)
  x <- rbinom(n, 1,0.2)

  z1 <- rbinom(n, 1, 0.37)
  v1 <- rnorm(n)
  d1 <- z1 * (u > 0.5 * v1)

  z2 <- rbinom(n, 1, 0.37)
  v2 <- rnorm(n)
  d2 <- z2 * (u > 0.5 * v1)

  y00 <- 0 + x * 2 + qnorm(u,0,1)
  y10 <- u  + x * 2 + qnorm(u,0,1)
  y01 <- 2 + x * 2 + qnorm(u,0,1)
  y11 <- 2 + u + x * 2 + qnorm(u,0,1)
  y <- d1 * d2 * y11  + (1 - d1) * d2 * y01 + d1 * (1 - d2) * y10 + (1 - d1) * (1 - d2) * y00

  value <- list()
  value$y <- y
  value$d1 <- d1
  value$d2 <- d2
  value$z1 <- z1
  value$z2 <- z2
  value$x <- x
  value <- data.frame(value)
  return(value)
}