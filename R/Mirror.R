#' @import knitr
#' @import MASS
#' @import tree
#' @import boot
#' @import bootstrap
#' @import DAAG
#' @import parallel
#' @import microbenchmark
#' @import GM
#' @import DSfdr
#' @import glmnet
#' @importFrom Rcpp evalCpp
#' @importFrom stats lm rnorm rbeta rbinom
#' @useDynLib SA24204176
NULL



#' @title A dataset used for Homework 0.
#' @name hde
#' @description This is a part of the South African heart disease data, 
#'    we use classification tree method to provide a model for 
#'    determining whether a patient has heart disease (CHD).
NULL

Gaussian_Mirror <- function(y, X){
  n <- length(y); p <- ncol(X) 
  stopifnot(length(y) == nrow(X), n > p + 1)
  X <- scale(X)
  
  GM_j <- function(j){
    zj <- rnorm(n)
    X_j <- X[, -j, drop = FALSE] 
    Nj <- diag(1,n) - X_j %*% solve(crossprod(X_j)) %*% t(X_j) 
    cj <- as.vector(sqrt(crossprod(Nj %*% X[ ,j]) / crossprod(Nj %*% zj)))
    xjp <- X[ ,j] + cj * zj; xjn <- X[ ,j] - cj * zj
    bj <- lm(y ~ xjp + xjn + X_j)$coefficients[2:3]
    Mj <- abs(bj[1] + bj[2]) - abs(bj[1] - bj[2])
    return(Mj)
  }
  
  M <- sapply(1:p, GM_j)
  names(M) <- colnames(X)
  return(M)
}

Data_split <- function(y, X){
  n <- length(y); p <- ncol(X) 
  stopifnot(length(y) == nrow(X), round(n/2) > p)
  X <- scale(X)
  
  group1 <- sample(n, round(n/2)) 
  X1 <- X[ group1, ]; y1 <- y[ group1, ]
  X2 <- X[-group1, ]; y2 <- y[-group1, ]
  beta1 <- lm(y1 ~ X1)$coefficient[-1] 
  beta2 <- lm(y2 ~ X2)$coefficient[-1] 
  M <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
  names(M) <- colnames(X)
  return(M)
}


Mirror <- function(method){
  function(y, X, q){
    Mirror_stat <- method(y, X)
    feature <- Selection(Mirror_stat, q)
    names(feature) <- colnames(X)[feature]
    return(list(feature=feature, mirror_statistic = Mirror_stat))
  }
}


#' @title Gauss Mirrors
#' @description Use Gaussian Mirrors to select features while controlling FDR for linear regression model.
#' @param y vector of the response variable
#' @param X the design matrix
#' @param q the designated FDR level to be controlled at
#' @return a list with the selected features and the mirror statistic of all features. 
#' @examples
#' \dontrun{
#' n <- 100; p <- 5
#' beta <- c(0,1,1,2,0)
#' X <- matrix(rnorm(n*p),n,p)
#' y <- X %*% beta + rnorm(n)
#' GM(y,X,0.1)$feature
#' }
#' @export
GM <- Mirror(Gaussian_Mirror)

#' @title Data Splitting
#' @description Use Data Splitting to select features while controlling FDR for linear regression model.
#' @param y vector of the response variable
#' @param X the design matrix
#' @param q the designated FDR level to be controlled at
#' @return a list with the selected features and the mirror statistic of all features. 
#' @examples
#' \dontrun{
#' n <- 100; p <- 5
#' beta <- c(0,1,1,2,0)
#' X <- matrix(rnorm(n*p),n,p)
#' y <- X %*% beta + rnorm(n)
#' DS(y,X,0.1)$feature
#' }
#' @export
DS <- Mirror(Data_split)

#' @title Multiple Data Splitting
#' @description Use Multiple Data Splitting to select features while controlling FDR for linear regression model.
#' @param y vector of the response variable
#' @param X the design matrix
#' @param q the designated FDR level to be controlled at
#' @param m the number of independent random sample splits.
#' @return a list with the selected features and the inclusion rate of all features. 
#' @examples
#' \dontrun{
#' n <- 100; p <- 5
#' beta <- c(0,1,1,2,0)
#' X <- matrix(rnorm(n*p),n,p)
#' y <- X %*% beta + rnorm(n)
#' MDS(y,X,0.1)$feature
#' }
#' @export
MDS <- function(y, X, q, m=50){
  p <- ncol(X)
  inclusion_mat <- matrix(0, m, p)
  colnames(inclusion_mat) <- colnames(X)
  for(k in 1:m){
    inclusion_mat[k, DS(y, X, q)$feature] <- 1
  }
  IR <- colMeans(inclusion_mat / pmax(rowSums(inclusion_mat), 1))
  return(list(feature = which(IR >= sort(IR)[cumsum(sort(IR)) > q][1]),
              inclusion_rate = IR))
}
