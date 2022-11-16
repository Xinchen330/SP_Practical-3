# Group 12 Practical 4
## Group Members: Xin Chen (s2340094)
##                Yicong Sun (s2445309)
##                Yihong Zhao (s2331659)

## Address of the github repo: https://github.com/Xinchen330/SP_Practical-4.git

## Contributions:

# Overview:
## This R script aims to give an independent implementation of Newton's method 
## for optimization. The minimum is judged by comparing the absolute value of 
## each element in the gradient vector with the threshold given by the user and 
## making sure that the Hessian matrix is positive definite at the minimum. The 
## step size to minimize the objective function is calculated by negating the 
## product of the inverse Hessian matrix and the gradient vector at a given set 
## of parameter values. If the user does not provide the Hessian matrix, it
## will be approximated by finite differencing the gradient vector. At each 
## step, two modifications are made to guarantee convergence. First, the 
## Hessian matrix might be perturbed by adding multiples of the identity 
## matrix until it is positive definite. Second, the step size might be 
## repeatedly halved if it fails to reduce the objective function. These two 
## modifications will ensure that the iteration converges to a minimum.

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

# Function newt to minimize the objective function func using Newton's 
# optimization method. The convergence is judged by comparing the absolute
# value of each element in the gradient vector with the threshold calculated by 
# multiplying the tolerance by the sum of the absolute value of the objective
# value and fscale
# Inputs:
## theta - the initial guess of parameter values
## func -- the objective function to minimize
## grad -- the gradient function
## hess -- the Hessian matrix, with default value NULL (not provided by the 
## user). If the Hessian matrix is not provided, the Hessian matrix will be 
## approximated by finite differencing the gradient vector
## ... -- pass extra arguments required by func, grad and hess
## tol -- convergence tolerence, with default value 1e-8
## fscale -- an estimate of the magnitude of the objective value near the 
## minimum, with default value 1
## maxit -- the maximum number of iterations to try before concluding that the
## Newton's method fails to converge, with default value 100
## max.half -- the maximum number of step halvings to try before concluding
## that the step fails to decrease the objective value, with default value 20
## eps -- the finite difference interval
# The function returns a list containing 5 elements:
## f -- the objective value at minimum
## theta -- the parameter values that minimizes the objective function
## iter -- the number of iterations to reach the minimum
## g -- the gradient vector at the minimum
## Hi -- the inverse Hessian matrix at the minimum
newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,
                 max.half=20,eps=1e-6) {
  ## Issue errors when the objective or derivatives are not finite at the
  ## initial theta
  if (abs(func(theta))==Inf | any(abs(grad(theta))==Inf) 
      | is.na(func(theta)) | any(is.na(grad(theta)))) {
    warning("The objective or derivatives are not finite at the initial 
            theta!")
  }
  iter <- 0 ## Initialise a counter for number of iterations
  ## Initialise a counter for the number of times the step was halved
  i_half <- 0
  ## The criteria for judging whether the gradient vector is zero
  threshold <- tol*(abs(func(theta)) + fscale)
  ## Case I: the Hessian matrix is provided
  if (!is.null(hess)) {
    while (any(abs(grad(theta)) >= threshold)) {
      h <- hess(theta) ## Hessian matrix
      ## The inverse of the Hessian matrix, guaranteed to be positive definite
      hinv <- hess_inv(h)
      ## Compute the descent direction
      delta <- -hinv %*% grad(theta)
      ## Repeatedly halve step sizes until the objective decreases
      while (func(theta+delta) >= func(theta) | is.na(func(theta+delta))) {
        delta <- delta/2
        i_half <- i_half+1 # Update counter
        ## Stop and print error message if the objective function is non-finite
        ## at max.half
        if ((is.na(func(theta+delta)) | abs(func(theta+delta))==Inf) & 
            i_half >= max.half) {
          stop("The objective function is non finite at max.half!")
        }
        ## Stop and print error message if the max.half is reached without 
        ## reducing the objective
        else if (func(theta+delta) >= func(theta) & i_half >= max.half) {
          stop("Failed to reduce the objective despite trying max.half step 
               halvings!")
        }
      }
      theta <- theta + delta ## Update parameter values
      iter <- iter + 1 ## Update counter
      ## Stop and print error message if the maxit is reached without 
      ## convergence
      if (any(abs(grad(theta)) >= threshold) & iter >= maxit) {
        stop("Maxit is reached without convergence!")
      }
    }
    ## Check if the Hessian is positive definite at convergence
    h <- hess(theta) ## Hessian matrix at convergence
    hinv <- try(chol2inv(chol(h)),silent=TRUE)
    ## Issue errors if the Hessian is not positive definite at convergence
    if (inherits(hinv,"try-error")) {
      stop("The Hessian is not positive definite at convergence")
    }
  }
  ## Case II: the Hessian matrix is not provided, approximate Hessian matrices 
  ## using finite difference
  else {
    while (any(abs(grad(theta)) >= threshold)) {
      ## Hesian matrix obtained by using finite difference
      hfd <- fd(theta,grad,eps)
      ## The inverse of the approximated Hessian matrix
      hinv <- hess_inv(hfd)
      ## Compute the descent direction
      delta <- -hinv %*% grad(theta)
      ## Repeatedly halve step sizes until the objective decreases
      while (func(theta+delta) >= func(theta) | is.na(func(theta+delta))) {
        delta <- delta/2
        i_half <- i_half+1 # Update counter
        ## Stop and print error message if the objective function is non-finite
        ## at max.half
        if ((is.na(func(theta+delta)) | abs(func(theta+delta))==Inf) & 
            i_half >= max.half) {
          stop("The objective function is non finite at max.half!")
        }
        ## Stop and print error message if the max.half is reached without 
        ## reducing the objective
        else if (func(theta+delta) >= func(theta) & i_half >= max.half) {
          stop("Failed to reduce the objective despite trying max.half step 
               halvings!")
        }
      }
      theta <- theta + delta ## Update parameter values
      iter <- iter + 1 ## Update counter
      ## Stop and print error message if the maxit is reached without 
      ## convergence
      if (any(abs(grad(theta)) >= threshold) & iter >= maxit) {
        stop("Maxit is reached without convergence!")
      }
    }
    ## Check if the Hessian is positive definite at convergence
    h <- fd(theta,grad,eps) ## Hessian matrix at convergence
    hinv <- try(chol2inv(chol(h)),silent=TRUE)
    ## Issue errors if the Hessian is not positive definite at convergence
    if (inherits(hinv,"try-error")) {
      stop("The Hessian is not positive definite at convergence")
    }
  }
  f <- func(theta) ## The value of the objective at the minimum
  g <- grad(theta) ## The gradient vector at the minimum
  Hi <- hinv ## Inverse Hessian matrix at the minimum
  ## A list containing the minimum objective function value f, value of 
  ## parameters theta, number of iterations iter, gradient vector at the 
  ## minimum g, the inverse Hessian matrix at the minimum Hi
  li_newt <- list(f=f,theta=theta,iter=iter,g=g,Hi=Hi)
}
