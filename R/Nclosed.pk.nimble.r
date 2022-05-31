# nimble data augmentation model for closed population estimate

library(nimble)
library(nimbleEcology)


# non-spatial closed model with pt
Nclosed.pk.code <- nimbleCode({
  psi ~ dunif(0, 1)     # DA parameter: E(N) = M*psi
  
  for (k in 1:K) {
    p[k] ~ dunif(0, 1)      # baseline encounter probability
  }

  for(i in 1:M) {
    #z[i] ~ dbern(psi)  # Is individual real?
    y[i,1:K] ~ dOcc_v(psi, p[1:K], len = K)
    ynew[i,1:K] ~ dOcc_v(psi, p[1:K], len = K)   # posterior predictive data
    sumyi[i] <- sum(y[i,1:K])   # calculate once in wrapper code
    sumynewi[i] <- sum(ynew[i,1:K])
  }
  
  # GOF metrics
  ## individual heterogeneity based on dispersion of individual captures (adapted from Royle et al. SCR book, used for spatial bins)
  Ti <- var(sumyi[1:M])/mean(sumyi[1:M])   # calculate once in wrapper code, because no "expected" component
  Tinew <- var(sumynewi[1:M])/mean(sumynewi[1:M])
  ## temporal heterogeneity
  for (k in 1:K) { 
    sumyk[k] <- sum(y[1:M,k])
    sumynewk[k] <- sum(ynew[1:M,k])
    esumyk[k] <- p[k] * psi * M
  }
  Tk <- sum(pow(sqrt(sumyk[1:K]) - sqrt(esumyk[1:K]), 2))
  Tknew <- sum(pow(sqrt(sumynewk[1:K]) - sqrt(esumyk[1:K]), 2))
  
  #N <- sum(z[1:M])                         # realized abundance (not available because z not directly estimated)
  EN <- psi*M
})



Nclosed.pk.fitfun <- function(data, M, rand = NULL, ni = 2000, nb = 1000, cm = NULL, cmcmc = NULL, 
                              model="Nclosed.pk.code") {
  
  if (!is.null(rand)) set.seed(rand)
  
  K <- data$K
  yobs <- data$Y
  nind <- dim(yobs)[1]
  y <- array(0, c(M,K))
  y[1:nind,] <- yobs
  # z <- c(rep(1, nind), rep(0, M - nind))
  
  constants <- list(K = K, M = M)
  
  inits <- function() {
    list(psi = runif(1, 0.5, 1),
         p = runif(K, 0, 1))
  }
  
  if (is.null(cm)) {
    m <- nimbleModel(get(model), constants = constants, data = list(y = y))
    cm <- compileNimble(m)
    conf <- configureMCMC(m, monitors = c("p", "EN","Ti","Tinew","Tk","Tknew"), thin=1)
    conf$addSampler(target = c("p","psi"), type = 'AF_slice')
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project=m)
  } else {
    cm$y <- y
  }
  
  fit1 = runMCMC(cmcmc, inits = inits, niter = ni, nburnin = nb, nchains = 3, samplesAsCodaMCMC = TRUE)
  
  return(list(fit1 = fit1, cm = cm, cmcmc = cmcmc))
}

