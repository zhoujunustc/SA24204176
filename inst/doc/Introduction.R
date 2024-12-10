## ----eval=FALSE---------------------------------------------------------------
# Gaussian_Mirror <- function(y, X){
#   n <- length(y); p <- ncol(X)
#   stopifnot(length(y) == nrow(X), n > p + 1)
#   X <- scale(X)
# 
#   GM_j <- function(j){
#     zj <- rnorm(n)
#     X_j <- X[, -j, drop = FALSE]
#     Nj <- diag(1,n) - X_j %*% solve(crossprod(X_j)) %*% t(X_j)
#     cj <- as.vector(sqrt(crossprod(Nj %*% X[ ,j]) / crossprod(Nj %*% zj)))
#     xjp <- X[ ,j] + cj * zj; xjn <- X[ ,j] - cj * zj
#     bj <- lm(y ~ xjp + xjn + X_j)$coefficients[2:3]
#     Mj <- abs(bj[1] + bj[2]) - abs(bj[1] - bj[2])
#     return(Mj)
#   }
# 
#   M <- sapply(1:p, GM_j)
#   names(M) <- colnames(X)
#   return(M)
# }
# 
# GM <- function(y, X, q){
#   Mirror_stat <- Gaussian_Mirror(y, X)
#   feature <- Selection(Mirror_stat, q)
#   names(feature) <- colnames(X)[feature]
#   return(list(feature=feature, mirror_statistic = Mirror_stat))
# }

## ----eval=FALSE---------------------------------------------------------------
# Data_split <- function(y, X){
#   n <- length(y); p <- ncol(X)
#   stopifnot(length(y) == nrow(X), round(n/2) > p)
#   X <- scale(X)
# 
#   group1 <- sample(n, round(n/2))
#   X1 <- X[ group1, ]; y1 <- y[ group1, ]
#   X2 <- X[-group1, ]; y2 <- y[-group1, ]
#   beta1 <- lm(y1 ~ X1)$coefficient[-1]
#   beta2 <- lm(y2 ~ X2)$coefficient[-1]
#   M <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
#   names(M) <- colnames(X)
#   return(M)
# }
# 
# DS <- function(y, X, q){
#   Mirror_stat <- Data_split(y, X)
#   feature <- Selection(Mirror_stat, q)
#   names(feature) <- colnames(X)[feature]
#   return(list(feature=feature, mirror_statistic = Mirror_stat))
# }

## ----eval=FALSE---------------------------------------------------------------
# MDS <- function(y, X, q, m=50){
#   p <- ncol(X)
#   inclusion_mat <- matrix(0, m, p)
#   colnames(inclusion_mat) <- colnames(X)
#   for(k in 1:m){
#     inclusion_mat[k, DS(y, X, q)$feature] <- 1
#   }
#   IR <- colMeans(inclusion_mat / pmax(rowSums(inclusion_mat), 1))
#   return(list(feature = which(IR >= sort(IR)[cumsum(sort(IR)) > q][1]),
#               inclusion_rate = IR))
# }

## -----------------------------------------------------------------------------
n <- 1000; p <- 30
set.seed(0)
beta <- sample(-2:2, p, replace = TRUE, prob = c(.1, .1, .5, .2, .1))
X <- matrix(rnorm(n * p), n, p)
y <- X %*% beta + rnorm(n)

## ----include=FALSE------------------------------------------------------------
library(SA24204176)
library(parallel)
library(glmnet)

## ----message=FALSE------------------------------------------------------------
q <- 0.1
set.seed(0)
GM.fit <- GM(y, X, q)
DS.fit <- DS(y, X, q)
MDS.fit <- MDS(y, X, q, m = 50)

res <- matrix(0, 7, p)
colnames(res) <- paste0("x", 1:p)
rownames(res) <- c("true", "myGM", "theGM", "myDS", "theGS", "myMDS", "theMDS")

res[1, beta != 0] <- 1
res[2, GM.fit$feature] <- 1
res[3, GM::gm(y, X, q)$gm_selected] <- 1
res[4, DS.fit$feature] <- 1
res[5, DSfdr::DS(X, y, 1, q)$DS_feature] <- 1
res[6, MDS.fit$feature] <- 1
res[7, DSfdr::DS(X, y, 50, q)$MDS_feature] <- 1

knitr::kable(res)

## ----message=FALSE------------------------------------------------------------
library(microbenchmark)
ts <- microbenchmark(
  myGM   = GM(y, X, q),
  theGM  = GM::gm(y, X, q),
  myDS   = DS(y, X, q),
  theDS  = DSfdr::DS(X, y, 1, q),
  myMDS  = MDS(y, X, q),
  theMDS = DSfdr::DS(X, y, 50, q),
  times  = 20
)
summary(ts)[ , c(1, 3, 5, 6)]

