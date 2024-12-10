## ----example1, echo=FALSE-----------------------------------------------------
library(knitr)
data <- data.frame(
  Name = c("Alice", "Bob", "Charlie", "David", "Eve", "Frank", "Grace", "Henry", "Ivy", "Jack"),
  Age = c(25, 30, 28, 32, 27, 35, 33, 29, 31, 36),
  Grade = c("A", "B", "A-", "B+", "A", "C", "A", "B-", "A", "A+"),
  Department = c("Marketing", "Finance", "HR", "Operations", "IT", "Sales", "Finance", "HR", "Marketing", "Operations"),
  Experience = c(3, 5, 4, 7, 2, 8, 6, 4, 5, 9)
)
kable(data, caption = "Table of Employees")

## ----example2, echo=FALSE-----------------------------------------------------
theta <- seq(0, 2*pi, length.out = 100)
r <- 1 - sin(theta)
x <- r * cos(theta)
y <- r * sin(theta)
plot(x, y, type = "l", col = "red", xlab = "x", ylab = "y", main = "Cardioid Curve")

## ----example3a,echo=FALSE,warning=FALSE---------------------------------------
library(SA24204176)
library(MASS)
data("hde")
kable(head(hde), caption = "South African heart disease")
hde$famhist <- as.factor(hde$famhist)
hde$chd <- as.factor(hde$chd)

## ----example3b,echo=FALSE, out.width = '100%',warning=FALSE-------------------
library(tree)
attach(hde)
tree.hde=tree(chd~. , hde)
chd.pred=predict(tree.hde, hde, type="class")
plot(tree.hde)
text(tree.hde,pretty=0,cex=0.7)

## ----ex3.4, echo=FALSE--------------------------------------------------------
n <- 10000
sigma <- c(0.1,0.2,0.5,1)
par(mfrow=c(2,2), mar=c(3,3,1,1))
set.seed(233)
for(s in sigma){
  U <- runif(n)
  X <- s*sqrt(-2*log(U)) 
  hist(X,main=bquote(sigma~"="~.(s)))
  abline(v=s,col = "blue", lty = 2, lwd=3)
}

## ----ex3.11, echo=FALSE, fig.width=10-----------------------------------------
n <- 1000
mixnorm <- function(p){
  k <- sample(c(0,3),n,replace=TRUE,prob=c(p,1-p))
  X <- rnorm(n,k,1)
}
set.seed(233)
par(mfrow=c(2,5))
for(p in seq(0.05,0.95,0.1)){
  hist(mixnorm(p),prob=TRUE,xlab="x",ylim=c(0,0.4),main=bquote(p[1]~'='~.(p)))
  curve(p*exp(-x^2/2)/sqrt(2*pi)+(1-p)*exp(-(x-3)^2/2)/sqrt(2*pi),add=TRUE)
  abline(v=0,col="red",lwd=2,lty=2)
  abline(v=3,col="red",lwd=2,lty=2)
}

## ----ex3.20, echo=FALSE-------------------------------------------------------
lambda <- c(1,2,3)
shape <- c(1,2,3)
rate <- c(1,2,3)

result <- data.frame(matrix(ncol=7,nrow=0))
n <- 10000
t0 <- 10
set.seed(233)
for(l in lambda){
  for(s in shape){
    for(r in rate){
      Xt <- replicate(n,expr={
        Tn <- rexp(200,l)
        Sn <- cumsum(Tn)
        Nt <- min(which(Sn>t0))-1
        sum(rgamma(Nt,s,r))
      })
      result <- rbind(result,list(l,s,r,l*t0*s/r,mean(Xt),
                                  l*t0*s*(s+1)/r^2,var(Xt)))
    }
  }
}
colnames(result) <- c("lambda","shape","rate",
                      "Theoretical mean","Sample mean",
                      "Theoretical variance","Sample variance")
knitr::kable(result)

## ----ex5.4, echo=FALSE--------------------------------------------------------
pbeta.mc <- function(x,a=3,b=3,M=10000){
  g <- function(y) y^(a-1)*(1-y)^(b-1)/beta(a,b)
  cdf <- numeric(length(x))
  for(i in 1:length(x)){
    Y <- runif(M)
    cdf[i] <- mean(x[i]*g(x[i]*Y))
  }
  return(cdf)
}

cdf <- matrix(0,3,9)
rownames(cdf) <- c("x","pbeta","estimation")
cdf[1,] <- (1:9)/10
cdf[2,] <- round(pbeta((1:9)/10,3,3),5)
set.seed(233)
cdf[3,] <- round(pbeta.mc((1:9)/10,3,3),5)
knitr::kable(cdf)

## ----ex5.9, echo=FALSE--------------------------------------------------------
M <- 20000
sigma <-1
set.seed(233)
u <- runif(M/2)
v1 <- runif(M/2)
v2 <- 1-u
g <- function(u) sigma*sqrt(-2*log(u))
mean1 <- (g(u)+g(v1))/2
mean2 <- (g(u)+g(v2))/2

## ----ex5.13, echo=FALSE-------------------------------------------------------
set.seed(233)
M <- 10000
U <- runif(M)
### f1
X1 <-  qnorm(1-pnorm(-1)*U)
Y1 <- pnorm(-1)*X1^2
theta1 <- mean(Y1)
var1 <- var(Y1)/M

### f2
X2 <- sqrt(1-2*log(U))
Y2 <- X2*exp(-1/2)/sqrt(2*pi)
theta2 <- mean(Y2)
var2 <- var(Y2)/M

## ----FS, include=FALSE--------------------------------------------------------
Fast_sort <- function(array){
  d <- length(array)
  index <- sample(1:d,1)
  value <- array[index]
  smaller <- which(array < value)
  arrayl <- array[smaller]
  arrayr <- array[-union(smaller,index)]
  
  if(length(arrayl)!=0){
    arrayl <- Fast_sort(arrayl)    
  }
  if(length(arrayr)!=0){
    arrayr <- Fast_sort(arrayr)  
  } 
  array <- c(arrayl,value,arrayr)
  return(array)
}

n <- c(1,2,4,6,8)*1e4
an <- numeric(5)
M <- 100

## ----MC, eval=FALSE, include=FALSE--------------------------------------------
# set.seed(233)
# for(i in 1:5){
#   an[i] <- mean(replicate(M,{
#     t0 <- Sys.time()
#     Fast_sort(sample(1:n[i]))
#     difftime(Sys.time(),t0)
#   }))
# }

## ----echo=FALSE---------------------------------------------------------------
bn <- c(0.1887577,0.3752437,0.7572155,1.1584966,1.6326623)
tn <- n*log(n)

## ----plot, echo=FALSE---------------------------------------------------------
plot(bn~tn,xlab=expression(t[n]==n*log(n)),ylab=expression(a[n]))
abline(lm(bn~tn))

## ----ex6.6, echo=FALSE--------------------------------------------------------
generation_6.6 <- function(n=100,M=10000){
  data <- matrix(rnorm(M*n),c(M,n))
  return(data)
}

inference_6.6 <- function(data){
  skew <- function(x){      ## sample skewness
    xbar <- mean(x)
    m3 <- mean((x-xbar)^3)
    m2 <- mean((x-xbar)^2)
    return(m3/m2^1.5)
  }
  alpha <- c(.025,.05,.95,.975)
  n <- ncol(data)
  M <- nrow(data)
  stat <- apply(data,1,skew)
  result <- matrix(NA,3,4)
  rownames(result) <- c("Estimated","Normal","SE")
  colnames(result) <- alpha
  result[1,] <- quantile(stat,alpha) ## 估计分位数
  result[2,] <- qnorm(alpha,0,sqrt(6/n)) ## 正态分位数
  var.exact <- 6*(n-2)/((n+1)*(n+3)) ## 精确方差
  result[3,] <- sqrt(alpha*(1-alpha)/(M*(dnorm(result[1,],0,sqrt(var.exact))^2)))
  return(result)
}

report_6.6 <- function(result){
  knitr::kable(round(result,4))
} 

set.seed(233)
report_6.6(inference_6.6(generation_6.6()))

## ----ex6.B, echo=FALSE--------------------------------------------------------
generation_6.B <- function(mu1,mu2,sigma1,sigma2,r,n=100,M=10000){
  data <- array(NA,dim=c(n,2,M))
  for(i in 1:M){
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    x3 <- x1*r + x2*sqrt(1-r^2)
    data[,1,i] <- mu1 + sigma1*x1
    data[,2,i] <- mu2 + sigma2*x3
  }
  return(data)
}

inference_6.B <- function(data,alpha=0.05){
  M <- dim(data)[3]
  P <- S <- K <- numeric(M)
  for(i in 1:M){
    P[i] <- cor.test(data[,1,i],data[,2,i],method="pearson")$p.value < alpha
    S[i] <- cor.test(data[,1,i],data[,2,i],method="spearman")$p.value < alpha
    K[i] <- cor.test(data[,1,i],data[,2,i],method="kendall")$p.value < alpha
  }
  return(list(Pearson=mean(P),Spearman=mean(S),Kendall=mean(K)))
}

report_6.B <- function(power){
  report <- data.frame(method=c("Pearson","Spearman","Kendall"),
                       power=c(power$Pearson,power$Spearman,power$Kendall))
  knitr::kable(report)
} 

set.seed(0)
report_6.B(inference_6.B(generation_6.B(0,0,1,1,0.0001)))

## ----echo=FALSE---------------------------------------------------------------
    A <- matrix(c('V (FP)','U (TN)','m0=950','S (TP)','T (FN)','m1=50','R','m-R','m=1000'),3)
    colnames(A) <- c('H0 is true','Ha is true','Total')
    rownames(A) <- c('Positive(reject H0)','Negative (accept H0)','Total')
    knitr::kable(A)

## ----echo=FALSE---------------------------------------------------------------
N <- 1000
N0 <- 950; N1 <- 50
H0 <- c(rep(T,N0),rep(F,N1))

m <- 10000
result <- array(NA,dim=c(3,2,m))
dimnames(result) <- list(c("FWER","FDR","TPR"),
                         c("Bonferroni correction","B-H correction"),
                         NULL)

multitest <- function(H0,pval,alpha){
  m <- length(H0)
  m0 <- sum(H0==T); m1 <- m-m0
  V <- sum((pval<alpha)[H0==T])
  S <- sum((pval<alpha)[H0==F])
  R <- V+S
  return(c(FWER=(V!=0), FDR=V/R, TPR=S/m1))
}

set.seed(233)
for(i in 1:m){
  p <- c(runif(N0),rbeta(N1,0.1,1))
  p.bon <- p.adjust(p,method='bonferroni')
  p.bh <- p.adjust(p,method='fdr')
  
  result[,1,i] <- multitest(H0,pval=p.bon,alpha=0.1)
  result[,2,i] <- multitest(H0,pval=p.bh,alpha=0.1)
}

knitr::kable(round(apply(result,c(1,2),mean),4))

## ----echo=FALSE---------------------------------------------------------------
library(boot)
data <- aircondit$hours

B <- 1e4; set.seed(233); lambdastar <- numeric(B)
lambda <- 1/mean(data)
for(b in 1:B){
  datastar <- sample(data,replace=TRUE)
  lambdastar[b] <- 1/mean(datastar)
}

## ----echo=FALSE---------------------------------------------------------------
lambda.boot <- function(dat,ind){
  #function to compute the MLE of lambda
  return(mean(dat[ind]))
}
boot.obj <- boot(data,statistic=lambda.boot,R=1e4)
ci <- boot.ci(boot.obj,type=c("norm","basic","perc","bca"))
ci.norm<-ci$norm[2:3];ci.basic<-ci$basic[4:5]
ci.perc<-ci$percent[4:5];ci.bca<-ci$bca[4:5]
ci.mat <- rbind(ci.norm,ci.basic,ci.perc,ci.bca)
rownames(ci.mat) <- c("Standard Normal","Basic","Percentile","BCa")
colnames(ci.mat) <- c("Lower bound","Upper bound")
knitr::kable(round(ci.mat,1))

## ----ex7.8--------------------------------------------------------------------
library(boot)
library(bootstrap)
data <- scor
n <- nrow(data)

proportion1 <- function(data,index){
  ## calculate the first principal component of the given data
  S <- cov(data[index,])
  lambda <- eigen(S)$values
  return(lambda[1]/sum(lambda))
} 

theta.hat <- proportion1(data,1:n)
theta.jack <- numeric(n)
for(i in 1:n){
  theta.jack[i] <- proportion1(data,(1:n)[-i])
}
bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)
se.jack <- sqrt((n-1)*mean((theta.jack-theta.hat)^2))
round(c(original=theta.hat,bias.jack=bias.jack,se.jack=se.jack),4)

## ----ex7.10-------------------------------------------------------------------
library(DAAG); attach(ironslag)
n <- length(magnetic)  
e1 <- e2 <- e3 <- e4 <- numeric(n)

for (k in 1:n) {
  y <- magnetic[-k]
  x <- chemical[-k]
  
  J1 <- lm(y ~ x)
  yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
  e1[k] <- magnetic[k] - yhat1
  
  J2 <- lm(y ~ x + I(x^2))
  yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +
    J2$coef[3] * chemical[k]^2
  e2[k] <- magnetic[k] - yhat2
  
  J3 <- lm(log(y) ~ x)
  logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
  yhat3 <- exp(logyhat3)
  e3[k] <- magnetic[k] - yhat3
  
  J4 <- lm(y ~ x + I(x^2) + I(x^3))
  yhat4 <-  J4$coef[1] + J4$coef[2] * chemical[k] +
    J4$coef[3] * chemical[k]^2 + J4$coef[4] * chemical[k]^3
  e4[k] <- magnetic[k] - yhat4
}


## ----e2-----------------------------------------------------------------------
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))

## ----R2-----------------------------------------------------------------------
y <- magnetic
x <- chemical

M1 <- lm(y ~ x)
R21 <- summary(M1)$adj.r.squared
M2 <- lm(y ~ x + I(x^2))
R22 <- summary(M2)$adj.r.squared
M3 <- lm(log(y) ~ x)
R23 <- summary(M3)$adj.r.squared
M4 <- lm(y ~ x + I(x^2) + I(x^3))
R24 <- summary(M4)$adj.r.squared
c(max(R21),max(R22),max(R23),max(R24))

## ----ex8.1a-------------------------------------------------------------------
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))
z <- c(x, y);K <- 1:26
detach(chickwts)
R <- 999
D <- numeric(R);

CM.stat <- function(x,y){
  ## function to calculate the Cramer-von Mises statistic 
  n <- length(x); m <- length(y)
  Fn <- ecdf(x); Gm <- ecdf(y)
  W2 <- m*n/(m+n)^2*(sum(Fn(x)-Gm(x))^2+sum(Fn(y)-Gm(y))^2)
  return(W2)
} 

## implement the permutation test
set.seed(233)
D0 <- CM.stat(x,y)
for (i in 1:R) {
  #generate indices k for the first sample
  k <- sample(K, size = 14, replace = FALSE)
  x1 <- z[k]; y1 <- z[-k]
  D[i] <- CM.stat(x1,y1)
} 
p <- mean(c(D0, D) >= D0)

## ----ex8.1b-------------------------------------------------------------------
hist(D, main = "", freq = FALSE, xlab = "D",breaks = "scott")
text(6,0.5, paste0("p.perm = ", p))
abline(v=D0,col='red',lwd=2)  # observed D

## -----------------------------------------------------------------------------
library(boot)
set.seed(233)
# generate the sample
n <- 100
x <- rnorm(n)
y <- x^2 + rnorm(n)
z <- cbind(x,y)

rho.boot <- function(data,ind){
  # return the Spearman rank correlation of the data
  cor(data[,1],data[ind,2],method="spearman")
}

# permutation test
boot.obj <- boot(data=z,statistic=rho.boot,sim="permutation",R=999)
tb <- c(boot.obj$t0,boot.obj$t) 

## -----------------------------------------------------------------------------
mean(abs(tb) >= abs(boot.obj$t0))

## -----------------------------------------------------------------------------
cor.test(x,y,method="spearman")$p.value

## ----chain9.3-----------------------------------------------------------------
chain_9.3 <- function(scale=1,location=0, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], 1)
    if (u[i] <= dcauchy(y,location,scale) / dcauchy(x[i-1],location,scale))
      x[i] <- y  
    else {
      x[i] <- x[i-1]
      k <- k + 1
    }
  }
  return(list(x=x, k=k))
}

## ----Rhat, echo=FALSE---------------------------------------------------------
Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  
  psi.means <- rowMeans(psi)     #row means
  B <- n * var(psi.means)        #between variance est.
  psi.w <- apply(psi, 1, "var")  #within variances
  W <- mean(psi.w)               #within est.
  v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
  r.hat <- v.hat / W             #G-R statistic
  return(r.hat)
}

## ----GR9.3, echo=FALSE--------------------------------------------------------
k <- 5          #number of chains to generate
n <- 15000      #length of chains
b <- 1000       #burn-in length

#initial values
x0 <- c(-10, -5, 0, 5, 10)

#generate the chains
set.seed(233)
X <- matrix(0, nrow=k, ncol=n)
for(i in 1:k){
  X[i, ] <- chain_9.3(x0=x0[i],N=n)$x
}
  
#compute diagnostic statistics(the median)
psi <- matrix(0, nrow=k, ncol=n)
for(i in 1:k){
  for(j in 1:n){
    psi[i,j] <- median(X[i,1:j])
  }
}

#plot the sequence of R-hat statistics
rhat <- rep(0, n)
for (j in (b+1):n)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

## ----run9.3, echo=FALSE-------------------------------------------------------
N <- 7000
set.seed(233)
x0 <- 0
rw <- chain_9.3(location=0,scale=1,x0,N)
plot(rw$x,type="l",ylab='X')

## ----decile9.3, echo=FALSE----------------------------------------------------
a <- seq(0.1,0.9,0.1)
Q <- qcauchy(a)
mc <- rw$x[(b+1):N]
Qrw <- quantile(mc,a)
knitr::kable(round(cbind(Q,Qrw),3))

## ----chain9.8-----------------------------------------------------------------
chain_9.8 <- function(a,b,n,x0,y0,N){
  #initialize constants and parameters
  chain <- matrix(0, N, 2) #the chain, a bivariate sample
  chain[1,] <- c(x0,y0) #initialize
  for (i in 2:N) {
    chain[i,1] <- rbinom(1,n,chain[i-1,2])
    chain[i,2] <- rbeta(1,chain[i,1]+a, n-chain[i,1]+b)
  }
  return(chain)
}

## ----GR9.8, echo=FALSE--------------------------------------------------------
a<-2;b<-3;n<-10
k <- 9          #number of chains to generate
N <- 15000      #length of chains
burn <- 1000       #burn-in length


#initial values
x0 <- rep(round(c(n/3,n/2,2*n/3)),3)
y0 <- rep(c(1/3,1/2,2/3),each=3)

#generate the chains
set.seed(233)
X <- array(0, dim=c(N,2,k))
for(i in 1:k){
  X[,,i] <- chain_9.8(a,b,n,x0=x0[i],y0=y0[i],N=N)
}

#compute diagnostic statistics(the mean of X*Y)
Y <- t(X[,1,]*X[,2,])
psi <- t(apply(Y,1,cumsum))
for (i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi))

#plot the sequence of R-hat statistics
rhat <- rep(0, N)
for (j in (burn+1):N)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(burn+1):N], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

## ----ex11.3(a)----------------------------------------------------------------
xk <- function(k,d,a){
  (-1)^k/(2*k+1)/(2*k+2)*
    exp((k+1)*log(sum(a^2))+lgamma((d+1)/2)+lgamma(k+3/2)-lgamma(k+1)-k*log(2)-lgamma(k+d/2+1))
}

## ----ex11.3(b)----------------------------------------------------------------
S <- function(d,a){
  Sn <- 0; x <- 1
  k <- 0
  while(abs(x)>1e-8 & k<1000){
    x <- xk(k,d,a)
    Sn <- Sn + x
    k <- k+1
  }
  stopifnot(k<1000)
  return(Sn)
}

## ----ex11.3(c)----------------------------------------------------------------
S(d=2,a=c(1,2))

## ----ex11.4-------------------------------------------------------------------
ck <- function(k,a){
  sqrt(a^2*k/(k+1-a^2))
}

A <- function(k){
  Sk <- function(k,a){
    pt(ck(k,a), df=k, lower.tail=FALSE)
  }
  res <- uniroot(function(a){Sk(k-1,a)-Sk(k,a)}, c(0.1,0.8*sqrt(k)))  
  return(res$root)
}

## ----ex11.5-------------------------------------------------------------------
A2 <- function(k){
  Sk <- function(k,a){
    2*exp(lgamma((k+1)/2)-lgamma(k/2))/sqrt(pi*k)*
      integrate(function(u){(1+u^2/k)^(-(k+1)/2)},lower=0,upper=ck(k,a))$value
  }
  logit <- function(u){log(u/(1-u))}
  res <- uniroot(function(a){logit(Sk(k-1,a))-logit(Sk(k,a))}, c(0.1,(1-log10(k)/3.5)*sqrt(k)))  
  return(res$root)
}

## ----run11.5, echo=FALSE------------------------------------------------------
k <- c(4:25,100,500,1000)
result <- matrix(0,25,3)
colnames(result) <- c("k","A(k) in Ex 11.4","A(k) in Ex 11.5") 
for(i in 1:25){
  result[i,1] <- k[i]
  result[i,2] <- A(k[i])
  result[i,3] <- A2(k[i])
}

knitr::kable(result)

## ----ex11.7-------------------------------------------------------------------
library(boot)
A1 <- rbind(c(2, 1, 1), c(1, -1, 3))
b1 <- c(2, 3)
a <- c(4, 2, 9)
simplex(a = a, A1 = A1, b1 = b1)

## ----ex11.1.3, eval=FALSE-----------------------------------------------------
# formulas <- list(
#   mpg ~ disp,
#   mpg ~ I(1 / disp),
#   mpg ~ disp + wt,
#   mpg ~ I(1 / disp) + wt
# )
# 
# ## for loops
# models <- vector("list", length(formulas))
# for(i in seq_along(formulas)){
#   models[[i]] <- lm(formulas[[i]], data = mtcars)
# }
# 
# ## lapply()
# models <- lapply(formulas, lm, data = mtcars)

## ----ex11.1.4, eval=FALSE-----------------------------------------------------
# bootstraps <- lapply(1:10, function(i) {
#   rows <- sample(1:nrow(mtcars), rep = TRUE)
#   mtcars[rows, ]
# })
# 
# ## for loop
# fits <- vector("list", length(bootstraps))
# for(i in seq_along(bootstraps)){
#   fits[[i]] <- lm(mpg ~ disp, data = bootstraps[[i]])
# }
# 
# ## lapply()
# fits <- lapply(bootstraps, lm, formula = mpg ~ disp) # no anonymous function

## ----ex11.1.5, eval=FALSE-----------------------------------------------------
# rsq <- function(mod) summary(mod)$r.squared
# 
# rsq_3 <- lapply(models, rsq)
# rsq_4 <- lapply(fits, rsq)

## ----ex11.2.3, eval=FALSE-----------------------------------------------------
# trials <- replicate(
#   100,
#   t.test(rpois(10, 10), rpois(7, 10)),
#   simplify = FALSE
# )
# 
# ## use an anonymous function
# pvalues <- sapply(trials, function(res) res$p.value)
# 
# ## get rid of the anonymous function
# pvalues <- sapply(trials, `[[`, "p.value")

## ----ex11.2.6, eval=FALSE-----------------------------------------------------
# library(parallel)
# mc.Map.vapply <- function(FUN, ..., FUN.VALUE, mc.cores = 1){
#   FUN <- match.fun(FUN)
#   res <- mcMap(FUN, ..., mc.cores = mc.cores)
#   vapply(res, function(x) x, FUN.VALUE = FUN.VALUE)
# }

## ----ex17.5.4, warning=FALSE--------------------------------------------------
chisq.test2 <- function(x, y){
  stopifnot(length(x) == length(y))
  mat.o <- table(x, y)
  df <- prod(dim(mat.o) - 1)
  
  xsum <- rowSums(mat.o)
  ysum <- colSums(mat.o)
  n <- sum(xsum)
  mat.e <- outer(xsum, ysum, "*")/n
  X2 <- sum((mat.o-mat.e)^2/mat.e)
  
  return(list(test.statistic = X2, df = df, 
              p_value = 1 - pchisq(X2, df)))
}

set.seed(233)
x <- rpois(1000,3)
y <- rpois(1000,5)
chisq.test(x, y)
chisq.test2(x, y)

library(microbenchmark)
microbenchmark(chisq.test(x, y), chisq.test2(x, y))

## ----ex17.5.5 table2----------------------------------------------------------
table2 <- function(x, y){
  stopifnot(length(x) == length(y))
  fx <- as.factor(x); fy <- as.factor(y)
  lx <- levels(fx); ly <- levels(fy)
  bin <- (as.numeric(fx) - 1L) + length(lx) * (as.numeric(fy) - 1L) + 1L
  out <- matrix(tabulate(bin, length(lx) * length(ly)), length(lx), length(ly))
  dimnames(out) <- list(lx, ly)
  class(out) <- "table"
  out
}
table(x,y)
table2(x,y)
microbenchmark(table(x,y), table2(x,y))

## ----ex17.5.5 chisq.test3, warning=FALSE--------------------------------------
chisq.test3 <- function(x, y){
  stopifnot(length(x) == length(y))
  mat.o <- table2(x, y)
  df <- prod(dim(mat.o) - 1)
  
  xsum <- rowSums(mat.o)
  ysum <- colSums(mat.o)
  n <- sum(xsum)
  mat.e <- outer(xsum, ysum, "*")/n
  X2 <- sum((mat.o-mat.e)^2/mat.e)
  
  return(list(test.statistic = X2, df = df, 
              p_value = 1 - pchisq(X2, df)))
}
chisq.test3(x, y)
microbenchmark(chisq.test(x, y), chisq.test2(x, y), chisq.test3(x, y))

## -----------------------------------------------------------------------------
library(Rcpp)
cppFunction("
NumericMatrix gibbsC(int a, int b, int n, int N) {
  NumericMatrix mat(N, 2);
  double x = n / 2, y = 0.5;
  mat(0, 0) = x; mat(0, 1) = y;
  for(int i=1; i<N; i++){
 	  x = rbinom(1, n, y)[0];
  	y = rbeta(1, x + a, n - x + b)[0];
    mat(i, 0) = x; mat(i, 1) = y;
  }
  return(mat);
}
")

## -----------------------------------------------------------------------------
gibbsR <- function(a, b, n, N){
  chain <- matrix(0, N, 2)
  x <- round(n/2); y <- 1/2;
  chain[1, ] <- c(x, y) 
  for (i in 2:N) {
    x <- rbinom(1, n, y)
    y <- rbeta(1, x + a, n - x + b) 
    chain[i, ] <- c(x, y)
  }
  return(chain)
}

## -----------------------------------------------------------------------------
library(Rcpp)

burn <- 1000
N <- 3000
set.seed(233)
chainR <- gibbsR(a = 2, b = 3, n = 10, N = N)[(burn + 1) : N, ]
chainC <- gibbsC(a = 2, b = 3, n = 10, N = N)[(burn + 1) : N, ]

## -----------------------------------------------------------------------------
qqplot(chainR[ ,1], chainC[ ,1], xlab = "x(R)", ylab = "x(C)")

## -----------------------------------------------------------------------------
qqplot(chainR[ ,2], chainC[ ,2], xlab = "y(R)", ylab = "y(C)")

## -----------------------------------------------------------------------------
library(microbenchmark)
ts <- microbenchmark(gibbsR = gibbsR(a = 2, b = 3, n = 10, N = N),
                     gibbsC = gibbsC(a = 2, b = 3, n = 10, N = N))
summary(ts)[ , c(1, 3, 5, 6)]

