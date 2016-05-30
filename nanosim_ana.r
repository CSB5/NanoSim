##### Mix model functions
# Poisson-Geometric distribution
rpoisgeommix <- function(n, l, p, w) {
  ifelse(runif(n) < w,
         rpois(n, l),
         rgeom(n, p))
}

dpoisgeommix <- function(x, l, p, w,log=FALSE) {
  r <- w * dpois(x, l) + (1 - w) * dgeom(x, p)
  if (log) log(r) else r
}

ppoisgeommix <- function(x, l, p, w, log=FALSE){
  r <- w * ppois(x, l) + (1 - w) * pgeom(x, p)
  if (log) log(r) else r
}


# Wei-Geom distribution
rweigeommix <- function(n, shape, p, w, scale = 1) {
  ifelse(runif(n) < w,
         rweibull(n, shape, scale),
         rgeom(n, p))
}

dweigeommix <- function(x, shape, p, w, scale = 1, log=FALSE) {
  r <- w * dweibull(x, shape, scale) + (1 - w) * dgeom(x, p)
  if (log) log(r) else r
}

pweigeommix <- function(x, shape,p, w, scale = 1, log=FALSE) {
  r <- w * pweibull(x, shape, scale) + (1 - w) * pgeom(x, p)
  if (log) log(r) else r
}

weigeom.fit <- function(p, shape1, shape2, scale1, scale2, p1, p2, prob1, prob2, step){
  d <- 1
  for (i in seq(shape1, shape2, step)){
    for (j in seq(scale1, scale2, step)){
      for (n in seq(p1, p2, step)){
        for (m in seq(prob1, prob2, step)){
          e <- pweigeommix(1:length(p), i, n, m, j)
          d1 <- max(abs(e - p))
          if (d1 < d){
            d = d1
            I <- i
            J <- j
            N <- n
            P <- m
          }
        }
      }
    }
  }
  estimate <- c(I, J, N, P, d)
  estimate
}

library(fitdistrplus)

####################### mismatch ########################
lambdas <- c(0.5339, 0.4673, 0.3971, 0.4345)
ps <- c(0.7192, 0.7193, 0.7211, 0.6973)
alphas <- c(0.2325,0.2930,0.3705,0.2715)

for(i in 1:4){
    ## data generated with poisson geometric mixture
    dat <- (rpoisgeommix(10000,lambdas[i],ps[i],alphas[i]))
    ## data generated with simple geometric distribution
    dat.sim <- rgeom(10000, ps[i])
    print(ks.test(dat, dat.sim))
}

####################### insertion ########################
lambdas <- c(0.9571, 1.381, 1.1810, 1.180)
ks <- c(0.9797, 1.2183, 1.3406, 1.3021)
ps <- c(0.3955, 0.4704, 0.5267, 0.4816)
alphas <- c(0.9031,0.6023,0.6023,0.4880)

for(i in 1:4){
    ## data generated with poisson geometric mixture
    dat <- (rweigeommix(1000,ks[i], ps[i],alphas[i], lambdas[i]))
    ## data generated with simple geometric distribution
    dat.sim <- rexp(1000, 1/lambdas[i])
    print(ks.test(dat, dat.sim))
}


ks.test(dat, dat.g)

plot(fitdist(dat, "geom"))

