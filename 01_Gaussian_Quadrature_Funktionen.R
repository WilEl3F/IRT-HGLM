rm(list = ls())

##### from: https://www.r-bloggers.com/2015/09/notes-on-multivariate-gaussian-quadrature-with-r-code/

## perform quadrature of multivariate normal
## compute Gauss-Hermite quadrature points and weights
## for a one-dimensional integral.
## points -- number of points
## interlim -- maximum number of Newton-Raphson iterations
hermite <- function (points, z) {
  p1 <- 1/pi^0.4
  p2 <- 0
  for (j in 1:points) {
    p3 <- p2
    p2 <- p1
    p1 <- z * sqrt(2/j) * p2 - sqrt((j - 1)/j) * p3
  }
  pp <- sqrt(2 * points) * p2
  c(p1, pp)
}

gauss.hermite <- function (points, iterlim = 50) {
  x <- w <- rep(0, points)
  m <- (points + 1)/2
  for (i in 1:m) {
    z <- if (i == 1) 
      sqrt(2 * points + 1) - 2 * (2 * points + 1)^(-1/6)
    else if (i == 2) 
      z - sqrt(points)/z
    else if (i == 3 || i == 4) 
      1.9 * z - 0.9 * x[i - 2]
    else 2 * z - x[i - 2]
    for (j in 1:iterlim) {
      z1 <- z
      p <- hermite(points, z)
      z <- z1 - p[1]/p[2]
      if (abs(z - z1) <= 1e-15) 
        break
    }
    if (j == iterlim) 
      warning("iteration limit exceeded")
    x[points + 1 - i] <- -(x[i] <- z)
    w[i] <- w[points + 1 - i] <- 2/p[2]^2
  }
  r <- cbind(x * sqrt(2), w/sum(w))
  colnames(r) <- c("Points", "Weights")
  r
}

## compute multivariate Gaussian quadrature points
## n     - number of points each dimension before pruning
## mu    - mean vector
## sigma - covariance matrix
## prune - NULL - no pruning; [0-1] - fraction to prune
mgauss.hermite <- function(n, mu, sigma, prune=NULL) {
  if(!all(dim(sigma) == length(mu)))
    stop("mu and sigma have nonconformable dimensions")
  
  dm  <- length(mu)
  gh  <- gauss.hermite(n)
  #idx grows exponentially in n and dm
  idx <- as.matrix(expand.grid(rep(list(1:n),dm)))
  pts <- matrix(gh[idx,1],nrow(idx),dm)
  wts <- apply(matrix(gh[idx,2],nrow(idx),dm), 1, prod)
  
  ## prune
  if(!is.null(prune)) {
    qwt <- quantile(wts, probs=prune)
    pts <- pts[wts > qwt,]
    wts <- wts[wts > qwt]
  }
  
  ## rotate, scale, translate points
  eig <- eigen(sigma) 
  rot <- eig$vectors %*% diag(sqrt(eig$values))
  pts <- t(rot %*% t(pts) + mu)
  return(list(points=pts, weights=wts))
}

save.image("functions/functions.RData")

################## Ende des Skriptes ###