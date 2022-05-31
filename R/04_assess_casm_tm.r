####### CENTRAL AMERICA ASSESSMENT #######


############################# NON-SPATIAL ##############################

# Bayesian data-augmentation model in NIMBLE
library(coda)
source("Nclosed.pk.nimble.r")
### CAmSM
M <- 2000
dl <- list(Y = array2binom(apply(y3d.camsm.1921,MARGIN=c(1,3),sum)), K=3)
out.nim.CR0.camsm <- Nclosed.pk.fitfun(data = dl, M = M, ni = 22500, nb = 2500) 
gelman.diag(out.nim.CR0.camsm$fit1, multivariate = F)
effectiveSize(out.nim.CR0.camsm$fit1)
plot(out.nim.CR0.camsm$fit1)
mcmcplots::rmeanplot(out.nim.CR0.camsm$fit1)
summary(out.nim.CR0.camsm$fit1)
#### GOF
sumyi <- c(rowSums(dl$Y), rep(0, M-dim(dl$Y)[1]))
Ti <- var(sumyi)/mean(sumyi)
fit1 <- as.matrix(out.nim.CR0.camsm$fit1)
sum(fit1[,"Tk"] > fit1[,"Tknew"])/dim(fit1)[1]
sum(Ti > fit1[,"Tinew"])/dim(fit1)[1]
### save
save(out.nim.CR0.camsm, Ti, file="nim.CR0.camsm.1921.QAB.rdata")


################################################ SPATIAL ################################################

source("SCR0pjk.1D.nimble.r")
library(coda)
## spatial
### 2019-2021
#### model boundary at 19.2 N
dl <- list(Y = y3d.camsm.1921, traplocs = traplocs.1921, ylim = c(7.25, 19.2), K=3)
out.nim.1921.QAB.19 <- SCR0pjk.1D.fitfun(data = dl, M = 2000, ni = 57000, nb = 2000)
gelman.diag(out.nim.1921.QAB.19$fit1$samples, multivariate = F)
effectiveSize(out.nim.1921.QAB.19$fit1$samples)
plot(out.nim.1921.QAB.19$fit1$samples)
mcmcplots::rmeanplot(out.nim.1921.QAB.19$fit1$samples)
summary(out.nim.1921.QAB.19$fit1$samples)
#hist(as.matrix(out.nim.1921.QAB.19$fit1$samples2)[,23:2022][as.logical(as.matrix(out.nim.1921.QAB.19$fit1$samples2)[,2024:4023])])
save(out.nim.1921.QAB.19, file="nim.SCR0pjk1D.camsm.1921.QAB.19.2.rdata")
#### model boundary at 18.6 N
dl <- list(Y = y3d.camsm.1921, traplocs = traplocs.1921, ylim = c(7.25, 18.6), K=3)
out.nim.1921.QAB.18 <- SCR0pjk.1D.fitfun(data = dl, M = 2000, ni = 57000, nb = 2000)
gelman.diag(out.nim.1921.QAB.18$fit1$samples, multivariate = F)
effectiveSize(out.nim.1921.QAB.18$fit1$samples)
plot(out.nim.1921.QAB.18$fit1$samples)
mcmcplots::rmeanplot(out.nim.1921.QAB.18$fit1$samples)
summary(out.nim.1921.QAB.18$fit1$samples)
save(out.nim.1921.QAB.18, file="nim.SCR0pjk1D.camsm.1921.QAB.18.6.rdata")
#### model boundary at 20.4 N
dl <- list(Y = y3d.camsm.1921, traplocs = traplocs.1921, ylim = c(7.25, 20.4), K=3)
out.nim.1921.QAB.20 <- SCR0pjk.1D.fitfun(data = dl, M = 2000, ni = 57000, nb = 2000)
gelman.diag(out.nim.1921.QAB.20$fit1$samples, multivariate = F)
effectiveSize(out.nim.1921.QAB.20$fit1$samples)
plot(out.nim.1921.QAB.20$fit1$samples)
mcmcplots::rmeanplot(out.nim.1921.QAB.20$fit1$samples)
summary(out.nim.1921.QAB.20$fit1$samples)
save(out.nim.1921.QAB.20, file="nim.SCR0pjk1D.camsm.1921.QAB.20.4.rdata")
#### post-process abundance with boundary uncertainty
##### range of northern boundaries for 19.2 state-space
attach("nim.SCR0pjk1D.camsm.1921.QAB.19.2.rdata")
sigma <- as.matrix(out.nim.1921.QAB.19$fit1$samples)[,"sigma"]
N.1921.19.2 <- as.matrix(out.nim.1921.QAB.19$fit1$samples)[,"N"] # original
z <- as.matrix(out.nim.1921.QAB.19$fit1$samples2)[,2024:4023]
s <- as.matrix(out.nim.1921.QAB.19$fit1$samples2)[,23:2022]
nlive.casm.19 <- matrix(NA,dim(z)[1],20)
nbounds <- seq(18, 19.2, length.out=20)
for (i in 1:dim(z)[1]) {
  slive <- s[i,][as.logical(z[i,])]
  for (j in 1:length(nbounds)) {
    nlive.casm.19[i,j] <- sum(slive <= nbounds[j])
  }
}
detach(2)
##### lower and upper envelope simulations  
###### lower
attach("nim.SCR0pjk1D.camsm.1921.QAB.20.4.rdata")    
z <- as.matrix(out.nim.1921.QAB.20$fit1$samples2)[,2024:4023]  
s <- as.matrix(out.nim.1921.QAB.20$fit1$samples2)[,23:2022]
nlive.casm.20 <- rep(NA,dim(z)[1]) 
nbound <- 18.6
for (i in 1:dim(z)[1]) {
  slive <- s[i,][as.logical(z[i,])] 
  nlive.casm.20[i] <- sum(slive <= nbound)
}
detach(2)
###### upper
attach("nim.SCR0pjk1D.camsm.1921.QAB.18.6.rdata")
nlive.casm.18 <- as.matrix(out.nim.1921.QAB.18$fit1$samples)[,"N"]
detach(2)
##### correction factor
# CF = 1/(1+ %bias), with birthdeaths+nocalfs rnorm(-5.3,1.8), and sexhet rnorm(-19.0, 9.9)
cf <- 1/(1-rnorm(6600,.053,.018)-rnorm(6600,.190,.099))
##### save
save(nlive.casm.18, nlive.casm.19, nlive.casm.20, cf, file="Nprocessed.SCR0pjk.1921.rdata")


#### Change between 2004-6 and 2019-21? ####
##### current estimate
attach("Nprocessed.SCR0pjk.1921.rdata")
##### Wade estimate
n.w <- 755
sd.w <- 0.242*755
lx <- log((n.w^2)/sqrt(n.w^2+sd.w^2))
ls <- sqrt(log(1+sd.w^2/n.w^2))
n.0406 <- rlnorm(length(nlive.casm.19), meanlog=lx, sdlog=ls)
rm(lx, ls, n.w, sd.w)
##### direct comparison of estimates
n.1921.cf <- as.vector(cf * nlive.casm.19)
lambda.dir <- (n.1921.cf/n.0406)^(1/15)
##### add individuals up to 14.5
attach("nim.SCR0pjk1D.camsm.1921.QAB.19.2.rdata")
z <- as.matrix(out.nim.1921.QAB.19$fit1$samples2)[,2024:4023]
s <- as.matrix(out.nim.1921.QAB.19$fit1$samples2)[,23:2022]
nlive.casm.14.5 <-rep(NA,dim(z)[1])
nbound <- 14.5
for (i in 1:dim(z)[1]) {
  slive <- s[i,][as.logical(z[i,])]
  nlive.casm.14.5[i] <- sum(slive <= nbound)
}
detach(2)
n.1921.145.cf <- as.vector(cf * nlive.casm.14.5)
lambda <- (n.1921.145.cf/n.0406)^(1/15)
save(n.1921.145.cf, lambda, lambda.dir, file="lambda.0406Wade.1921SCR0pjk.rdata")


#### Goodness of Fit
#attach("nim.SCR0pjk1D.camsm.1921.QAB.19.2.rdata")
out.nim <- out.nim.1921.QAB.19
samples <- as.matrix(out.nim$fit1$samples2)
niter <- dim(samples)[1]
##### Ensure we have the nodes needed to simulate new datasets
M <- 2000
dl <- list(Y = y3d.camsm.1921, traplocs = traplocs.1921, ylim = c(7.25, 19.2), K=3)
X <- dl$traplocs; J <- length(X)
nind <- dim(dl$Y)[1]
y <- array(0, c(M, J, dl$K))
y[1:nind,,] <- dl$Y
constants <- list(X = X, K = dl$K, M = M, J = J, ylim = dl$ylim)
m.scr01d <- nimbleModel(SCR0pjk.1D.code, constants = constants, data = list(y = y))
cm.scr01d <- compileNimble(m.scr01d)
dataNodes <- m.scr01d$getNodeNames(dataOnly = TRUE, returnScalarComponents = TRUE)
parentNodes <- m.scr01d$getParents(dataNodes, stochOnly = TRUE)
simNodes <- m.scr01d$getDependencies(parentNodes, self = FALSE)   # includes p, N, y
# var name order from mcmc in out.nim: p0 (18), psi, s, sigma, z
vars <- dimnames(samples)[[2]]   # use this instead of parentNodes since also includes psi
nsbin <- 10
binlims <- dl$ylim[1] + (0:nsbin)*((dl$ylim[2]-dl$ylim[1])/nsbin)
ppc <- data.frame(Ts = rep(NA, niter), Tsnew = rep(NA, niter), 
                  Tij = rep(NA, niter), Tijnew = rep(NA, niter),
                  Ti = rep(NA, niter), Tinew = rep(NA, niter),
                  Tj = rep(NA, niter), Tjnew = rep(NA, niter))
set.seed(1)
system.time({
  for(ni in 1:niter) {
    values(cm.scr01d, vars) <- samples[ni, ]
    cm.scr01d$simulate(simNodes, includeData = TRUE)
    ysimflat <- values(cm.scr01d, dataNodes)   # Note data order different for ysimflat <- values(cm.scr01d, dataNodes) vs y2 <- values(cm.scr01d, "y"); see order of dataNodes
    ysim <- aperm(array(ysimflat, dim=c(dl$K, J, M)), perm=c(3,2,1))
    pflat <- values(cm.scr01d, "p")   # calculated p from scratch for two nodes to double check alignment
    p <- array(pflat, dim=c(M,J,dl$K))
    s <- values(cm.scr01d, "s")
    z <- values(cm.scr01d, "z")
    sran <- runif(M, dl$ylim[1], dl$ylim[2])
    # Goodness of fit components
    ## spatial randomness
    sbin <- rep(NA, nsbin)
    snewbin <- rep(NA, nsbin)
    for (b in 1:nsbin) {
     sbin[b] <- sum((z > 0) & (s > binlims[b]) & (s <= binlims[b+1]))
     snewbin[b] <- sum((z > 0) & (sran > binlims[b]) & (sran <= binlims[b+1]))
    }
    esbin <- sum(z)/nsbin
    ppc$Ts[ni] <- sum(pow(sqrt(sbin[1:nsbin]) - sqrt(esbin), 2))
    ppc$Tsnew[ni] <- sum(pow(sqrt(snewbin[1:nsbin]) - sqrt(esbin), 2))
    ## observations
    esumyij <- matrix(NA, nrow=M, ncol=J)
    sumyij <- matrix(NA, nrow=M, ncol=J)
    sumynewij <- matrix(NA, nrow=M, ncol=J)
    for (j in 1:J) {
      for (i in 1:M) {
        esumyij[i,j] <- sum(z[i]*p[i,j,])
        sumyij[i,j] <- sum(y[i,j,])
        sumynewij[i,j] <- sum(ysim[i,j,])
      }
    }
    esumyj <- colSums(esumyij)
    sumyj <- colSums(sumyij)
    sumynewj <- colSums(sumynewij)
    esumyi <- rowSums(esumyij)
    sumyi <- rowSums(sumyij)
    sumynewi <- rowSums(sumynewij)
    ppc$Tij[ni] <- sum(pow(sqrt(sumyij) - sqrt(esumyij), 2))
    ppc$Tijnew[ni] <- sum(pow(sqrt(sumynewij) - sqrt(esumyij), 2))
    ppc$Tj[ni] <- sum(pow(sqrt(sumyj) - sqrt(esumyj), 2))
    ppc$Tjnew[ni] <- sum(pow(sqrt(sumynewj) - sqrt(esumyj), 2))
    ppc$Ti[ni] <- sum(pow(sqrt(sumyi) - sqrt(esumyi), 2))
    ppc$Tinew[ni] <- sum(pow(sqrt(sumynewi) - sqrt(esumyi), 2))
  }
})
sum(ppc$Ts < ppc$Tsnew)/niter # 0.43
sum(ppc$Tj < ppc$Tjnew)/niter # 0.38
sum(ppc$Tij < ppc$Tijnew)/niter # 0.84
sum(ppc$Ti < ppc$Tinew)/niter  # 0.77
save(ppc, out.nim.1921.QAB.19, file = "nim.SCR0pjk1D.camsm.1921.QAB.19.2.rdata")
# note that need separate objects for separate chains if want to be able continue chains
# incorporating asymm in model leads to shift of model weight to less spatial "north" traps term, high N estim still 

# Mainland Mexico by subtraction
## estimate abundance via Chao Mth as in Calambokidis and Barlow 2020
## previous estimate from C & B 2020: x=4973, se=239
attach("Nprocessed.SCR0pjk.1921.rdata")
lx <- log((4973^2)/sqrt(4973^2+239^2))
ls <- sqrt(log(1+239^2/4973^2))
n.wc <- rlnorm(length(nlive.casm.19), meanlog=lx, sdlog=ls)
n.mmex <- n.wc - as.vector(cf * nlive.casm.19)
# x = 3478, cv = 0.100
# same as analytical answer, so use analytical:
4973-mean(nlive.casm.19*cf)   # 3477
sqrt(239^2 + sd(nlive.casm.19*cf)^2)   # 350; cv = 0.101
lx <- log((3477^2)/sqrt(3477^2+350^2))
ls <- sqrt(log(1+350^2/3477^2))
qlnorm(0.2, lx, ls)
detach(2)
