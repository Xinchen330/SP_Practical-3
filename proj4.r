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

newt(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,
     eps=1e-6) {
  ## Issue errors when the objective or derivatives are not finite at the
  ## initial theta
  if (func(theta)==Inf | any(grad(theta)==Inf)) {
    warning("The objective or derivatives are not finite at the initial 
            theta!")
  }
  ## If Hessian matrix is not provided, approximate the Hessian matrix by
  ## finite differencing the gradient vector
  if (hess=NULL) {
    hess <- fd(theta,grad,eps)
  }
}
