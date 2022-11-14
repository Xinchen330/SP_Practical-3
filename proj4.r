# Group 12 Practical 4
## Group Members: Xin Chen (s2340094)
##                Yicong Sun (s2445309)
##                Yihong Zhao (s2331659)

## Address of the github repo: 

## Contributions:

# Overview:

# Function fd to approximate the Hessian matrix by finite differencing the
# gradient vector
# Inputs:
## theta -- a vector of parameter values
## grad -- the gradient function
## eps -- the finite difference interval, with default value 1e-6
# The function fd returns the approximated Hessian matrix hfd
fd <- function(theta,grad,eps=1e-6) {
  ## Initialise a matrix for storing derivatives
  u_hfd <- matrix(0,length(theta),length(theta))
  ## Loop through elements to evaluate the upper triangular part of the 
  ## matrix
  for (j in 1:ncol(u_hfd)) {
    for (i in 1:j) {
      theta1 <- theta ## Make a copy of theta
      theta1[i] <- theta1[i]+eps ## Increase theta[i] by eps
      ## Approximate derivatives
      u_hfd[i,j] <- (grad(theta1)-grad(theta))[j]/eps
    }
  }
  ## Some pre-processing before copying the upper triangle part to lower
  ## triangle
  pre_hfd <- 2*u_hfd - diag(diag(u_hfd))
  ## Use the fact that (t(A)+A)/2 is symmetric for all square matrix A
  hfd <- (t(pre_hfd) + pre_hfd)/2
  return(hfd)
}

# Function hess_inv to find the inverse of a given Hessian matrix
# Inputs:
## h -- the Hessian matrix to be inverted
# The function hess_inv returns the inverse of the given Hessian matrix h,
# which is guaranteed to be positive definite
hess_inv <- function(h) {
  ## Test if the Hessian is positive definite with Cholesky decomposition
  hinv <- try(chol2inv(chol(h)),silent = TRUE)
  while (inherits(hinv,"try-error")) {
    ## Keep adding identity matrix to force the Hessian to be positive definite
    h <- h + diag(ncol(h))
    hinv <- try(chol2inv(chol(h)),silent=TRUE)
  }
  return(hinv)
}

newt(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,
     eps=1e-6) {
  ## Issue errors when the objective or derivatives are not finite at the
  ## initial theta
  if (func(theta)==Inf | any(grad(theta)==Inf)) {
    warning("The objective or derivatives are not finite at the initial 
            theta!")
  }
  iter <- 0 ## Initialise a counter for number of iterations
  ## The criteria for judging whether the gradient vector is zero
  threshold <- tol*(abs(func(theta)) + fscale)
  ## Case I: the Hessian matrix is provided
  if (!is.null(hess)) {
    while (any(abs(grad(theta)) >= threshold)) {
      h <- hess(theta) ## Hessian matrix
      ## The inverse of the Hessian matrix, guaranteed to be positive definite
      hinv <- hess_inv(theta,h)
      ## Compute the descent direction
      delta <- -hinv %*% grad(theta)
      ## Initialise a counter for the number of times the step was halved
      i_half <- 0
      ## Repeatedly halve step sizes until the objective decreases
      while (func(theta+delta) >= func(theta)) {
        delta <- delta/2
        i_half <- i_half+1 # Update counter
        ## Stop and print error message if the max.half is reached without 
        ## reducing the objective
        if (i_half > max.half) {
          stop("Failed to reduce the objective despite trying max.half step 
               halvings")
        }
      }
      theta <- theta + delta ## Update parameter values
    }
  }
  ## Case II: the Hessian matrix is not provided, approximate Hessian matrices 
  ## using finite difference
  else {
    hess <- fd(theta,grad,eps)
  }
}
