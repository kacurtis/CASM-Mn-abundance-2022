# nimble 1D SCR for closed populations

library(nimble)
library(nimbleEcology)

SCR0pjk.1D.code <- nimbleCode({
  for (j in 1:J) {
    for (k in 1:K) {
      p0[j,k] ~ dunif(0, 1)      # baseline encounter probability
    }
  }
  sigma ~ dunif(0,50)  # scale parameter of encounter function
  psi ~ dunif(0, 1)     # DA parameter: E(N) = M*psi
  
  for(i in 1:M) {
    z[i] ~ dbern(psi)  # Is individual real?
    s[i] ~ dunif(ylim[1], ylim[2]) # 1D spatial coordinate
    # # dist between activity center and trap
    # d[i,1:J] <- abs(s[i] - X[1:J])
    for(j in 1:J) {
      p[i,j,1:K] <- p0[j,1:K] * exp(-(s[i] - X[j])^2/(2*sigma^2))  # capture prob at trap j
      y[i,j,1:K] ~ dOcc_v(z[i], p[i,j,1:K], len = K)   # recode so vector is 1:M looped over J and K?
    }
  }
  
  N <- sum(z[1:M])   # realized abundance
  EN <- psi*M
  #D <- N/area   # realized density
})


SCR0pjk.1D.fitfun <- function(data, M, rand = NULL, ni = 2000, nb = 1000, cm.scr01d = NULL, cmcmc.scr01d = NULL, 
                              model="SCR0pjk.1D.code") {
  
  if (!is.null(rand)) set.seed(rand)
  
  X <- data$traplocs
  K <- data$K
  J <- length(X)
  ylim <- data$ylim
  # area <- ylim[2] - ylim[1]
  
  yobs <- data$Y
  nind <- dim(yobs)[1]
  y <- array(0, c(M,J,K))
  y[1:nind,,] <- yobs
  z <- c(rep(1, nind), rep(0, M - nind))
  sst <- runif(M, ylim[1], ylim[2])
  X2 <- array(X, dim = c(J, K))
  for (i in 1:nind) {
    if (sum(y[i, , ]) == 0) 
      next
    sst[i] <- mean(X2[y[i, , ] > 0])
  }

  constants <- list(X = X, K = K, M = M, J = J, ylim = ylim)   # area = area
  
  inits.scr01d <- function() {
    list(psi = runif(1, 0, 1),   # this seems to be ignored - maybe determined by z init? 
         p0 = array(runif(J*K, 0, 0.3), dim=c(J,K)), 
         sigma = runif(1, 0, 50), 
         s = sst,
         z=z)
  }
  
  if (is.null(cm.scr01d)) {
    m.scr01d <- nimbleModel(get(model), constants = constants, data = list(y = y))   #, calculate = FALSE)
    cm.scr01d <- compileNimble(m.scr01d)
    if (model=="SCR0pjk.1Dasymm.code") {
      conf.scr01d <- configureMCMC(m.scr01d, monitors = c("p0", "beta1", "sigma", "N"), thin=10)
      conf.scr01d$addSampler(target = c("beta1","sigma","psi"), type = 'AF_slice')
    } else {
      conf.scr01d <- configureMCMC(m.scr01d, monitors = c("p0", "sigma", "N"), thin=1,
                                   monitors2=c("s","z","p0","sigma","psi"), thin2=25)
      conf.scr01d$removeSamplers(c("psi","sigma")) 
      conf.scr01d$addSampler(target = c("psi","sigma"), type = 'AF_slice')  
    }
    mcmc.scr01d <- buildMCMC(conf.scr01d)
    cmcmc.scr01d <- compileNimble(mcmc.scr01d, project=m.scr01d)
  } else {
    cm.scr01d$y <- y
  }

  fit1 = runMCMC(cmcmc.scr01d, inits = inits.scr01d, niter = ni, nburnin = nb, nchains = 3, samplesAsCodaMCMC = TRUE)
  
  return(list(fit1 = fit1, cm.scr01d = cm.scr01d, cmcmc.scr01d = cmcmc.scr01d))
}
