#!/usr/bin/env Rscript

# New simulation study that explores the new 'informed' model that we've
# been discussing the past two weeks.  Massimo has created the code to
# produce the density of the pc prior on alpha1 | alpha2 and alpha2 | alpha1.
#
# In addition to the pc prior, we employ independent gamma priors on
# alpha1 and alpha2.
#
# We will compare methods based on the sparse FMM of silvia et al and also
# the Mixture of FMM of Miller et al.
#
# Our comparison will be carried out using five metrics.  They are as follows
#
# 1. true K+ - mode of posterior of K+ (this measures bias)
# 2. sum_{k-1}^K {(U - k)^2*pr(K+ = k | data)} (this measures )
# 3. cross-sectional MSE for density estimate
# 4. true K+ - K+ estimated from salso partition estimate
# 5. sum_i sum_{l<i} ((I[i~l] - Pr(zi=zl|y))^2
#
# We will generate data in two ways.
#
# 1. using the model based on fixed U and random alpha1=U and alpha2=1e-5?
# 2. Using approach we used before namely using
#    mu_k = c(0, 3, 6, 9)
#    sig2_k = 0.75*c(1, 1, 1, 1)



rm(list=ls())

#install.packages("TruncatedDistributions", repos="http://R-Forge.R-project.org")
library(TruncatedDistributions) # needed for the truncated exponential
library(MCMCpack) # needed for the rdirichlet function
library(miscPack) # needed to fit our models version
library(mixAK) # needed for NMixMCMC function (RJMCMC finite mixture)
library(BNPmix) # needed for PYdensity function
library(mclust) # needed for Mclust and densityMclust functions
library(AntMAN) # needed for the AM_mcmc_fit function
library(salso)

# install.packages("MCMCpack","mixAK","BNPmix","mclust","AntMAN","salso")

sessionInfo()

# set nobs in {100,  1000} - to begin
# set K in {25} - to begin
# set U in {1, 2, 5, 10}
# set datatype in {1, 2}
#       1 - from model - *just use this one to begin*
#       2 - as before
#
# There are 2 x 1 x 4 x 2  scenarios

# nohup ./SimulationStudy_asym.R 25 100 25 1 1 > out1.txt &

args <- commandArgs(TRUE)
# args <- c(100, 1000,  25, 2, 2)

ndata <- as.numeric(args[1])
nobs <- as.numeric(args[2])
K <- as.numeric(args[3])
U <- as.numeric(args[4])
datatype <- as.numeric(args[5])


cat("ndata = ", ndata, "\n")
cat("nobs = ", nobs, "\n")
cat("K = ", K, "\n")
cat("U = ", U, "\n")
cat("datatype = ", datatype, "\n")


nsamples <- 1e5 # monte carlo samples to get lambda and alpha density
agrid_length <- 1e5 # denseness of values of alpha

agrid <- exp(seq(log(1e-5), log(1), length.out = agrid_length))


Uprior <- c(1,2,5,10)

Kplus_var <- matrix(NA, nrow=ndata, ncol=7*length(Uprior)+5)
Kplus_bias <- matrix(NA, nrow=ndata, ncol=7*length(Uprior)+5)
density_mse <- matrix(NA, nrow=ndata, ncol=7*length(Uprior)+5)
ari <- matrix(NA, nrow=ndata, ncol=7*length(Uprior)+5)
Kplus_salso <- matrix(NA, nrow=ndata, ncol=7*length(Uprior)+5)
ccp_diff <- matrix(NA, nrow=ndata, ncol=7*length(Uprior)+5)
extra_w <- matrix(NA, nrow=ndata, ncol=7*length(Uprior)+5)

colnames(Kplus_var) <- colnames(Kplus_bias) <- colnames(density_mse) <-
colnames(ari) <- colnames(Kplus_salso) <- colnames(ccp_diff) <-
colnames(extra_w) <-
  c(apply(expand.grid(Uprior,c("Gam_alpha2F","Gam_alpha1F","Gam")),1,function(x) paste0("U_",as.numeric(x[1]),"_",as.character(x[2]))),
          "Sylvia_Sparse","FM_RevJump", "MBC", "DPM", "FM_Raffaele",
    apply(expand.grid(Uprior,c("PC_alpha2F_01")),1,function(x) paste0("U_",as.numeric(x[1]),"_",as.character(x[2]))),
    apply(expand.grid(Uprior,c("Gam_ProbKp_01")),1,function(x) paste0("U_",as.numeric(x[1]),"_",as.character(x[2]))),
    apply(expand.grid(Uprior,c("PC_alpha2F_09")),1,function(x) paste0("U_",as.numeric(x[1]),"_",as.character(x[2]))),
    apply(expand.grid(Uprior,c("Gam_ProbKp_09")),1,function(x) paste0("U_",as.numeric(x[1]),"_",as.character(x[2]))))

ndens_y <- 1000

niter <- 50000
nburn <- 40000
nthin <- 10
nout <- (niter - nburn)/nthin

# Function that relabels component indicators with using
# the Blocked GIBBS sampler of weighted DDP and DP models
relabel = function(z){
  as.numeric(factor(z, levels=unique(z), labels=1:length(unique(z))))
}

# These values came from Massimo, except the first one which I just inserted
if(nobs==100){
  shape <- c(0.1, 0.794, 1.804, 22.303)
  rate <- c(0.05, 0.265, 0.601, 7.434)
}
if(nobs==1000){
  shape <- c(0.1, 0.659, 1.041, 1.499)
  rate <- c(0.05, 0.220, 0.347, 0.500)
}
# These values came from Massimo, except the first one which I just inserted
if(nobs==100){
  shape2 <- c(0.009, 0.015, 0.021, 0.027)
  rate2 <- c(0.001, 0.005, 0.007, 0.009)
}
if(nobs==1000){
  shape2 <- c(0.009, 0.014, 0.018, 0.020)
  rate2 <- c(0.001, 0.005, 0.006, 0.007)
}


# These lambda values were supplied by Massimo and calcualted before
# simulation study (except for the first one)
if(nobs==100){
  lambda09 <- c(0.1, 0.01641652,0.01605877,0.01320394)
  lambda01 <- c(0.1, 0.8032921,0.4969684, 0.4164075)
}

if(nobs==1000){
  lambda09 <- c(0.1, 0.01318369,0.01087476,0.007295646)
  lambda01 <- c(0.1, 0.5464173, 0.2922722, 0.1808752)
}


# d <- 1
for(d in 1:ndata){

  cat("d = ", d, "\n")

  set.seed(ndata*(datatype-1) + d)
  cat("seed = ", ndata*(datatype-1) + d, "\n")

  print(date())

  if(datatype == 1){
    zi <- sample(1:U, nobs, replace=TRUE)
    yi <- rep(0, nobs)
    muk <- seq(0, 3*(U-1), by=3)
    sigk <- 0.5*rep(1,U)
    for(i in 1:U){
      yi <- yi + rnorm(nobs,muk[i], sigk[i])*(zi==i)
    }
    Utrue <- U
    ygrid <- seq(min(yi)-1, max(yi)+1, length=ndens_y)
    Dtrue <- sapply(1:ndens_y, function(i) sum((1/U)*dnorm(ygrid[i], muk, sigk)))
  }
  if(datatype == 2){
    a1 <- U
    a2 <- 1e-3
    pis <- rdirichlet(1, c(rep(a1,U),rep(a2,K-U)))
    zi <- sample(1:K, nobs, replace=TRUE, prob=pis)
    muk <- rnorm(K, 0, 3)
    sigk <- runif(K, 0, 1)
    yi <- rnorm(nobs, muk[zi], sigk[zi])
    ygrid <- seq(min(yi)-1, max(yi)+1, length=ndens_y)
    # hist(yi, freq=FALSE, breaks=25)
    Utrue <- length(unique(zi))
    Dtrue <- sapply(1:ndens_y, function(i) sum(pis*dnorm(ygrid[i],muk, sigk)))
  }


  Rho_true <- relabel(zi)
  ppt <- matrix(0, ncol=nobs, nrow=nobs)
  for(ii in 1:nobs){
    for(iii in 1:nobs){
      if(zi[ii] == zi[iii]) ppt[ii,iii] <- 1
    }
  }
#  hist(yi, freq=FALSE, breaks=25)

  yi_s <- scale(yi)
  yi_s <- yi
  #ss_results1a, First run

  alpha2_val_fixed <- 1e-5

  m <- mean(yi_s)
  fits1 <- fits2 <- fits3 <- fits4 <- fits5 <- fits6 <- fits7 <- list()

  for(u in 1:length(Uprior)){
    cat("u = ", u, "\n")
    cat("nobs = ", nobs, "\n")

    cat("fit model 1 \n")
    # fits the gamma prior alpha1 unknown, alpha2 fixed
    fits1[[u]] <- informed_mixture(y=yi_s, K=K,
                                   alpha_prior_type='centered',
                                   alpha_prior_dist="gamma",
                                   mu_sigma_prior = 2, # muk ~ N(m0, s20), sigma2k ~ IG(a0, b0)
                                   basemodel=0,
                                   U=Uprior[u],
                                   alpha1_val = 0.5,
                                   alpha2_val = alpha2_val_fixed,
                                   update_alpha1=TRUE,
                                   update_alpha2=FALSE,
                                   a1_gam=10, b1_gam=1/(10*Uprior[u]), # recall that b1_gam is scale
                                   ndens_y=ndens_y,
                                   hierarchy="NO",
                                   m0=mean(yi_s), s20=10^2, a0=3, b0 = 2, # prior mean of sig2 is b0/(a0-1)
                                   niter=niter, nburn=nburn, nthin=nthin)

    kpp1 <- table(fits1[[u]]$kp)/nout
    Kplus_var[d,u] <- sum(((Utrue - as.numeric(names(kpp1)))^2*kpp1))
    Kplus_bias[d,u] <- as.numeric(names(
      which.max(
        table(apply(fits1[[u]]$z, 1, function(x) length(unique(x))))
      )
     )
    ) - Utrue #posterior mode of K+.  Sylvia suggests this
    density_mse[d,u] <- sum((Dtrue - apply(fits1[[u]]$density,2,mean))^2)
    ari[d,u] <- mean(apply(fits1[[u]]$z, 1, function(x) adjustedRandIndex(Rho_true, x)))
    Kplus_salso[d,u] <- length(table(salso(fits1[[u]]$z))) - Utrue
    ccp_diff[d,u] <- mean(((ppt - psm(fits1[[u]]$z))[upper.tri(ppt, diag=FALSE)])^2)
    extra_w[d, u] <- mean(sapply(1:nout, function(x) sum(sort(fits1[[u]]$w[x,], decreasing=TRUE)[-(1:fits1[[u]]$kp[x])])))

    # fits the gamma prior alpha1 fixed, alpha2 unknown
    fits2[[u]] <- informed_mixture(y=yi_s, K=K,
                                   alpha_prior_type='centered',
                                   alpha_prior_dist="gamma",
                                   mu_sigma_prior = 2, # muk ~ N(m0, s20), sigma2k ~ IG(a0, b0)
                                   basemodel=0,
                                   U=Uprior[u],
                                   alpha1_val = 0.5,
                                   alpha2_val = 1.0,
                                   update_alpha1=FALSE,
                                   update_alpha2=TRUE,
                                   a1_gam=10, b1_gam=1/(10*Uprior[u]),
                                   ndens_y=ndens_y,  # recall that b1_gam is scale
                                   hierarchy="NO",
                                   m0=mean(yi_s), s20=10^2, a0=3, b0 = 2, # prior mean of sig2 is b0/(a0-1)
                                   niter=niter, nburn=nburn, nthin=nthin)

    kpp1 <- table(apply(fits2[[u]]$z, 1, function(x) length(unique(x))))/nout
    Kplus_var[d,u+length(Uprior)] <- sum(((Utrue - as.numeric(names(kpp1)))^2*kpp1))
    Kplus_bias[d,u+length(Uprior)] <- as.numeric(names(
      which.max(
        table(apply(fits2[[u]]$z, 1, function(x) length(unique(x))))
      )
     )
    ) - Utrue #posterior mode of K+.  Sylvia suggests this
    density_mse[d,u+length(Uprior)] <- sum((Dtrue - apply(fits2[[u]]$density,2,mean))^2)
    ari[d,u+length(Uprior)] <- mean(apply(fits2[[u]]$z, 1, function(x) adjustedRandIndex(Rho_true, x)))
    Kplus_salso[d,u+length(Uprior)] <- length(table(salso(fits2[[u]]$z))) - Utrue
    ccp_diff[d,u+length(Uprior)] <- mean(((ppt - psm(fits2[[u]]$z))[upper.tri(ppt, diag=FALSE)])^2)
    extra_w[d, u+length(Uprior)] <-mean(sapply(1:nout, function(x) sum(sort(fits2[[u]]$w[x,], decreasing=TRUE)[-(1:fits2[[u]]$kp[x])])))


    # fits the gamma prior alpha1 unknown, alpha2 unknown
    fits3[[u]] <- informed_mixture(y=yi_s, K=K,
                                    alpha_prior_type='centered',
                                    alpha_prior_dist="gamma",
                                    mu_sigma_prior = 2, # muk ~ N(m0, s20), sigma2k ~ IG(a0, b0)
                                    basemodel=0,
                                    U=Uprior[u],
                                    alpha1_val = 0.5,
                                    alpha2_val = alpha2_val_fixed,
                                    update_alpha1=TRUE,
                                    update_alpha2=TRUE,
                                    a1_gam=10, b1_gam=1/(10*Uprior[u]),
                                    a2_gam=1, b2_gam=1/(10*(K-Uprior[u])),
                                    ndens_y=ndens_y,  # recall that b1_gam is scale
                                    hierarchy="NO",
                                    m0=mean(yi_s), s20=10^2, a0=3, b0 = 2, # prior mean of sig2 is b0/(a0-1)
                                    niter=niter, nburn=nburn, nthin=nthin)


    kpp1 <- table(apply(fits3[[u]]$z, 1, function(x) length(unique(x))))/nout
    Kplus_var[d,u+2*length(Uprior)] <- sum(((Utrue - as.numeric(names(kpp1)))^2*kpp1))
    Kplus_bias[d,u+2*length(Uprior)] <- as.numeric(names(
      which.max(
        table(apply(fits3[[u]]$z, 1, function(x) length(unique(x))))
      )
     )
    ) - Utrue #posterior mode of K+.  Sylvia suggests this
    density_mse[d,u+2*length(Uprior)] <- sum((Dtrue - apply(fits3[[u]]$density,2,mean))^2)
    ari[d,u+2*length(Uprior)] <- mean(apply(fits3[[u]]$z, 1, function(x) adjustedRandIndex(Rho_true, x)))
    Kplus_salso[d,u+2*length(Uprior)] <- length(table(salso(fits3[[u]]$z))) - Utrue
    ccp_diff[d,u+2*length(Uprior)] <- mean(((ppt - psm(fits3[[u]]$z))[upper.tri(ppt, diag=FALSE)])^2)
    extra_w[d, u+2*length(Uprior)] <- mean(sapply(1:nout, function(x) sum(sort(fits3[[u]]$w[x,], decreasing=TRUE)[-(1:fits3[[u]]$kp[x])])))


    # fits the PC prior with alpha1 unknown and alpha2 fixed with tail probability 0.1
    fits4[[u]] <- informed_mixture(y=yi_s, K=K,
                                   alpha_prior_type='centered',
                                   alpha_prior_dist="pc",
                                   mu_sigma_prior = 2, # muk ~ N(m0, s20), sigma2k ~ IG(a0, b0)
                                   mylambda = lambda01[u],
                                   U=Uprior[u],
                                   tail.prob=0.1,
                                   alpha1_val = 0.5,
                                   alpha2_val = alpha2_val_fixed,
                                   update_alpha1=TRUE,
                                   update_alpha2=FALSE,
                                   ndens_y=ndens_y,
                                   hierarchy="NO",
                                   m0=mean(yi_s), s20=10^2, a0=3, b0 = 2, # prior mean of sig2 is b0/(a0-1)
                                   niter=niter, nburn=nburn, nthin=nthin)

    kpp1 <- table(apply(fits4[[u]]$z, 1, function(x) length(unique(x))))/nout
    Kplus_var[d,u+2*length(Uprior)+9] <- sum(((Utrue - as.numeric(names(kpp1)))^2*kpp1))
    Kplus_bias[d,u+2*length(Uprior)+9] <- as.numeric(names(
      which.max(
        table(apply(fits4[[u]]$z, 1, function(x) length(unique(x))))
      )
     )
    ) - Utrue #posterior mode of K+.  Sylvia suggests this
    density_mse[d,u+2*length(Uprior)+9] <- sum((Dtrue - apply(fits4[[u]]$density,2,mean))^2)
    ari[d,u+2*length(Uprior)+9] <- mean(apply(fits4[[u]]$z, 1, function(x) adjustedRandIndex(Rho_true, x)))
    if(mean(apply(fits4[[u]]$z,1,sum)) == nobs){
      Kplus_salso[d,u+2*length(Uprior)+9] <- 1 - Utrue
    } else {
      Kplus_salso[d,u+2*length(Uprior)+9] <- length(table(salso(fits4[[u]]$z))) - Utrue
    }
    ccp_diff[d,u+2*length(Uprior)+9] <- mean(((ppt - psm(fits4[[u]]$z))[upper.tri(ppt, diag=FALSE)])^2)
    extra_w[d, u+2*length(Uprior)+9] <- mean(sapply(1:nout, function(x) sum(sort(fits4[[u]]$w[x,], decreasing=TRUE)[-(1:fits4[[u]]$kp[x])])))


    # fits the gamma prior for alpha1 unknown and alpha2.
    # The gamma prior parameters are informed by a probability
    # statement associated with K+ with tail probability 0.1
    # (similar to the PC prior eliciation)
    fits5[[u]] <- informed_mixture(y=yi_s, K=K,
                                   alpha_prior_type='centered',
                                   alpha_prior_dist="gamma",
                                   mu_sigma_prior = 2, # muk ~ N(m0, s20), sigma2k ~ IG(a0, b0)
                                   U=Uprior[u],
                                   alpha1_val = 0.5,
                                   alpha2_val = alpha2_val_fixed,
                                   update_alpha1=TRUE,
                                   update_alpha2=FALSE,
                                   a1_gam=shape[u], b1_gam=1/rate[u],
                                   ndens_y=ndens_y,  # recall that b1_gam is scale
                                   hierarchy="NO",
                                   m0=mean(yi_s), s20=10^2, a0=3, b0 = 2, # prior mean of sig2 is b0/(a0-1)
                                   niter=niter, nburn=nburn, nthin=nthin)

    kpp1 <- table(apply(fits5[[u]]$z, 1, function(x) length(unique(x))))/nout
    Kplus_var[d,u+2*length(Uprior)+13] <- sum(((Utrue - as.numeric(names(kpp1)))^2*kpp1))
    Kplus_bias[d,u+2*length(Uprior)+13] <- as.numeric(names(
      which.max(
        table(apply(fits5[[u]]$z, 1, function(x) length(unique(x))))
      )
    )
    ) - Utrue #posterior mode of K+.  Sylvia suggests this
    density_mse[d,u+2*length(Uprior)+13] <- sum((Dtrue - apply(fits5[[u]]$density,2,mean))^2)
    ari[d,u+2*length(Uprior)+13] <- mean(apply(fits5[[u]]$z, 1, function(x) adjustedRandIndex(Rho_true, x)))
    if(mean(apply(fits5[[u]]$z,1,sum)) == nobs){
      bdif <- 1 - Utrue
    } else {
      bdif <- length(table(salso(fits5[[u]]$z))) - Utrue
    }
    Kplus_salso[d,u+2*length(Uprior)+13] <- length(table(salso(fits5[[u]]$z))) - Utrue
    ccp_diff[d,u+2*length(Uprior)+13] <- mean(((ppt - psm(fits5[[u]]$z))[upper.tri(ppt, diag=FALSE)])^2)
    extra_w[d, u+2*length(Uprior)+13] <- mean(sapply(1:nout, function(x) sum(sort(fits5[[u]]$w[x,], decreasing=TRUE)[-(1:fits5[[u]]$kp[x])])))


    # fits the PC prior with alpha1 unknown and alpha2 fixed and tail probability equal to 0.9
    fits6[[u]] <- informed_mixture(y=yi_s, K=K,
                                   alpha_prior_type='centered',
                                   alpha_prior_dist="pc",
                                   mu_sigma_prior = 2, # muk ~ N(m0, s20), sigma2k ~ IG(a0, b0)
                                   mylambda = lambda09[u],
                                   U=Uprior[u],
                                   tail.prob=0.9,
                                   alpha1_val = 0.5,
                                   alpha2_val = alpha2_val_fixed,
                                   update_alpha1=TRUE,
                                   update_alpha2=FALSE,
                                   ndens_y=ndens_y,  # recall that b1_gam is scale
                                   hierarchy="NO",
                                   m0=mean(yi_s), s20=10^2, a0=3, b0 = 2, # prior mean of sig2 is b0/(a0-1)
                                   niter=niter, nburn=nburn, nthin=nthin)

    kpp1 <- table(apply(fits6[[u]]$z, 1, function(x) length(unique(x))))/nout
    Kplus_var[d,u+2*length(Uprior)+17] <- sum(((Utrue - as.numeric(names(kpp1)))^2*kpp1))
    Kplus_bias[d,u+2*length(Uprior)+17] <- as.numeric(names(
      which.max(
        table(apply(fits6[[u]]$z, 1, function(x) length(unique(x))))
      )
    )
    ) - Utrue #posterior mode of K+.  Sylvia suggests this
    density_mse[d,u+2*length(Uprior)+17] <- sum((Dtrue - apply(fits6[[u]]$density,2,mean))^2)
    ari[d,u+2*length(Uprior)+17] <- mean(apply(fits6[[u]]$z, 1, function(x) adjustedRandIndex(Rho_true, x)))
    if(mean(apply(fits6[[u]]$z,1,sum)) == nobs){
      bdif <- 1 - Utrue
    } else {
      bdif <- length(table(salso(fits6[[u]]$z))) - Utrue
    }
    Kplus_salso[d,u+2*length(Uprior)+17] <- bdif
    ccp_diff[d,u+2*length(Uprior)+17] <- mean(((ppt - psm(fits6[[u]]$z))[upper.tri(ppt, diag=FALSE)])^2)
    extra_w[d, u+2*length(Uprior)+17] <- mean(sapply(1:nout, function(x) sum(sort(fits6[[u]]$w[x,], decreasing=TRUE)[-(1:fits6[[u]]$kp[x])])))


    # fits the gamma prior for alpha1 unknown and alpha2 fixed.
    # The gamma prior parameters are informed by a probability
    # statement associated with K+ with tail probability 0.9
    # Note this is accomplished by using shape2 and rate2
    fits7[[u]] <- informed_mixture(y=yi_s, K=K,
                                   alpha_prior_type='centered',
                                   alpha_prior_dist="gamma",
                                   mu_sigma_prior = 2, # muk ~ N(m0, s20), sigma2k ~ IG(a0, b0)
                                   U=Uprior[u],
                                   alpha1_val = 0.5,
                                   alpha2_val = alpha2_val_fixed,
                                   update_alpha1=TRUE,
                                   update_alpha2=FALSE,
                                   a1_gam=shape2[u], b1_gam=1/rate2[u],
                                   ndens_y=ndens_y,  # recall that b1_gam is scale
                                   hierarchy="NO",
                                   m0=mean(yi_s), s20=10^2, a0=3, b0 = 2, # prior mean of sig2 is b0/(a0-1)
                                   niter=niter, nburn=nburn, nthin=nthin)


    kpp1 <- table(apply(fits7[[u]]$z, 1, function(x) length(unique(x))))/nout
    Kplus_var[d,u+2*length(Uprior)+21] <- sum(((Utrue - as.numeric(names(kpp1)))^2*kpp1))
    Kplus_bias[d,u+2*length(Uprior)+21] <- as.numeric(names(
      which.max(
        table(apply(fits7[[u]]$z, 1, function(x) length(unique(x))))
      )
    )
    ) - Utrue #posterior mode of K+.  Sylvia suggests this
    density_mse[d,u+2*length(Uprior)+21] <- sum((Dtrue - apply(fits7[[u]]$density,2,mean))^2)
    ari[d,u+2*length(Uprior)+21] <- mean(apply(fits7[[u]]$z, 1, function(x) adjustedRandIndex(Rho_true, x)))
    if(mean(apply(fits7[[u]]$z,1,sum)) == nobs){
      bdif <- 1 - Utrue
    } else {
      bdif <- length(table(salso(fits7[[u]]$z))) - Utrue
    }
    Kplus_salso[d,u+2*length(Uprior)+21] <- bdif
    ccp_diff[d,u+2*length(Uprior)+21] <- mean(((ppt - psm(fits7[[u]]$z))[upper.tri(ppt, diag=FALSE)])^2)
    extra_w[d, u+2*length(Uprior)+21] <- mean(sapply(1:nout, function(x) sum(sort(fits7[[u]]$w[x,], decreasing=TRUE)[-(1:fits7[[u]]$kp[x])])))

  }

  # fit Sylvia model i.e., the sparse
  cat("Fit Sylvia's model", "\n")
  pc3 <- informed_mixture(y=yi_s, K=K,
                          alpha_prior_type="sparse",
                          alpha_prior_dist="gamma",
                          mu_sigma_prior = 2, # muk ~ N(m0, s20), sigma2k ~ IG(a0, b0)
                          basemodel=0,
                          m0=mean(yi_s), s20=10^2, a0=3, b0 = 2, # prior mean of sig2 is b0/(a0-1)
                          ndens_y=1000, a_gam=10, b_gam=1/(10*K),
                          niter=niter, nburn=nburn, nthin=nthin)

  kpp1 <- table(apply(pc3$z, 1, function(x) length(unique(x))))/nout
  Kplus_var[d,u+2*length(Uprior)+1] <- sum(((Utrue - as.numeric(names(kpp1)))^2*kpp1))
  Kplus_bias[d,u+2*length(Uprior)+1] <- as.numeric(names(
    which.max(
      table(apply(pc3$z, 1, function(x) length(unique(x))))
    )
   )
  ) - Utrue #posterior mode of K+.  Sylvia suggests this
  density_mse[d,u+2*length(Uprior)+1] <- sum((Dtrue - apply(pc3$density,2,mean))^2)
  ari[d,u+2*length(Uprior)+1] <- mean(apply(pc3$z, 1, function(x) adjustedRandIndex(Rho_true, x)))
  Kplus_salso[d,u+2*length(Uprior)+1] <- length(table(salso(pc3$z))) - Utrue
  ccp_diff[d,u+2*length(Uprior)+1] <- mean(((ppt - psm(pc3$z))[upper.tri(ppt, diag=FALSE)])^2)
  extra_w[d, u+2*length(Uprior)+1] <- mean(sapply(1:nout, function(x) sum(sort(pc3$w[x,], decreasing=TRUE)[-(1:pc3$kp[x])])))


  # Finite mixture using RJMCMC as implemented in mixAK package
  rj <- NMixMCMC(y0 = cbind(yi_s), scale=list(shift=0, scale=1),
                 prior=list(priorK="uniform", Kmax=25),
                 nMCMC=c(nburn, nout, nthin, nout))  # this is the function.
  pdens1 <- NMixPredDensMarg(rj, grid=ygrid)
#  lines(pdens1$x[[1]], pdens1$dens[[1]], col='red')
  Kplus_var[d,u+2*length(Uprior)+2] <- sum(((Utrue - as.numeric(names(rj$propK)))^2*rj$propK))
  Kplus_bias[d,u+2*length(Uprior)+2] <- as.numeric(names(which.max(rj$propK))) - Utrue
  density_mse[d,u+2*length(Uprior)+2] <- sum((Dtrue - pdens1$dens[[1]])^2)


  # Compare fit from mclust?
  mbc <- Mclust(yi_s, G=1:25, verbose=FALSE) #
  dmbc <- densityMclust(yi_s, G=1:25, plot=FALSE, verbose=FALSE) #
  pdmbc <- predict(dmbc, newdata=ygrid)
  Kplus_bias[d,u+2*length(Uprior)+3] <- length(unique(mbc$classification)) - Utrue
  density_mse[d,u+2*length(Uprior)+3] <- sum((Dtrue - pdmbc)^2)
  ari[d,u+2*length(Uprior)+3] <- adjustedRandIndex(Rho_true, mbc$classification)



  # DPM
  dp1 <- dp_mixture(y=yi_s, N=30, m=mean(yi_s), v=10^2, a=3, b=2, alpha=1, ndens_y = 1000,
                    niter=niter,nburn=nburn,nthin=nthin)
  kpp1 <- table(apply(dp1$z, 1, function(x) length(unique(x))))/nout
  Kplus_var[d,u+2*length(Uprior)+4] <- sum(((Utrue - as.numeric(names(kpp1)))^2*kpp1))
  Kplus_bias[d,u+2*length(Uprior)+4] <- as.numeric(names(
    which.max(
      table(apply(dp1$z, 1, function(x) length(unique(x))))
    )
   )
  ) - Utrue #posterior mode of K+.  Sylvia suggests this
  density_mse[d,u+2*length(Uprior)+4] <- sum((Dtrue - apply(dp1$density,2,mean))^2)
  ari[d,u+2*length(Uprior)+4] <- mean(apply(dp1$z, 1, function(x) adjustedRandIndex(Rho_true, x)))
  Kplus_salso[d,u+2*length(Uprior)+4] <- length(table(salso(dp1$z))) - Utrue
  ccp_diff[d,u+2*length(Uprior)+4] <- mean(((ppt - psm(dp1$z))[upper.tri(ppt, diag=FALSE)])^2)
  extra_w[d, u+2*length(Uprior)+4] <- mean(sapply(1:nout, function(x){
                                                            nc <- length(unique(dp1$z[x,]))
                                                            sum(sort(dp1$w[x,], decreasing=TRUE)[-(1:nc)])
                                                           })
                                           )


  # AntMAN package.  This is Raffaele's group.
  mixture_uvn_params = AM_mix_hyperparams_uninorm (m0=mean(yi_s), k0=1, nu0=3, sig02=4)
  components_prior = AM_mix_components_prior_pois (Lambda=U+5)
  weights_prior = AM_mix_weights_prior_gamma(gamma=1/(U+5))
#  components_prior = AM_mix_components_prior_pois (init=3, a=1, b=1)
#  weights_prior = AM_mix_weights_prior_gamma(init=2, a=1, b=1)
  mcmc_params = AM_mcmc_parameters(niter=niter, burnin=nburn, thin=nthin, verbose=0)
  rf <- AM_mcmc_fit(
    y = c(yi_s),
    mix_kernel_hyperparams = mixture_uvn_params,
    mix_components_prior =components_prior,
    mix_weight_prior = weights_prior,
    mcmc_parameters = mcmc_params)

  rf$density <- dp1$density
  for(jj in 1:nout){
    for(ii in 1:length(ygrid)){
      rf$density[jj,ii] <- sum(rf$W[[jj]]*dnorm(ygrid[ii], unlist(rf$mu[[jj]]), sqrt(unlist(rf$sig2[[jj]]))))
    }
  }
  kpp1 <- table(rf$K)/nout
  Kplus_var[d,u+2*length(Uprior)+5] <- sum(((Utrue - as.numeric(names(kpp1)))^2*kpp1))
  Kplus_bias[d,u+2*length(Uprior)+5] <- as.numeric(names(which.max(table(rf$K)))) - Utrue #posterior mode of K+.  Sylvia suggests this
  density_mse[d,u+2*length(Uprior)+5] <- sum((Dtrue - apply(rf$density,2,mean))^2)
  ari[d,u+2*length(Uprior)+5] <- mean(sapply(rf$CI, function(x) adjustedRandIndex(Rho_true, x)))
  rf$zi <- matrix(unlist(rf$CI), nrow=nout, byrow=TRUE)
  Kplus_salso[d,u+2*length(Uprior)+5] <- length(table(salso(rf$zi))) - Utrue
  ccp_diff[d,u+2*length(Uprior)+5] <- mean(((ppt - psm(rf$zi))[upper.tri(ppt, diag=FALSE)])^2)
  extra_w[d, u+2*length(Uprior)+5] <- mean(sapply(1:nout, function(x){
                                                            nc <- length(unique(rf$CI[[x]]))
                                                            sum(sort(rf$W[[x]], decreasing=TRUE)[-(1:nc)])
                                                          })
                                          )


  #dir <- paste0("~/Research/BYU/InformedFiniteMixtures/analysis/ss_results14/Kplus_var_nobs_",nobs,
  dir <- paste0("/home/gpage299/shared/InformedFiniteMixtures/ss_results14/Kplus_var_nobs_",nobs,
                #                "_K_",K,"_vtype_",v_type,"_A_",A,"_datatype_",datatype,".txt")
                "_K_",K,"_datatype_",datatype,"_U_",U,".txt")
  write.table(Kplus_var, file=dir, row.names=FALSE, col.names=TRUE)


  #dir <- paste0("~/Research/BYU/InformedFiniteMixtures/analysis/ss_results14/Kplus_bias_nobs_",nobs,
  dir <- paste0("/home/gpage299/shared/InformedFiniteMixtures/ss_results14/Kplus_bias_nobs_",nobs,
                #                "_K_",K,"_vtype_",v_type,"_A_",A,"_datatype_",datatype,".txt")
                "_K_",K,"_datatype_",datatype,"_U_",U,".txt")
  write.table(Kplus_bias, file=dir, row.names=FALSE, col.names=TRUE)

  #dir <- paste0("~/Research/BYU/InformedFiniteMixtures/analysis/ss_results14/density_mse_nobs_",nobs,
  dir <- paste0("/home/gpage299/shared/InformedFiniteMixtures/ss_results14/density_mse_nobs_",nobs,
                #                "_K_",K,"_vtype_",v_type,"_A_",A,"_datatype_",datatype,".txt")
                "_K_",K,"_datatype_",datatype,"_U_",U,".txt")
  write.table(density_mse, file=dir, row.names=FALSE, col.names=TRUE)

  #dir <- paste0("~/Research/BYU/InformedFiniteMixtures/analysis/ss_results14/ari_mn_nobs_",nobs,
  dir <- paste0("/home/gpage299/shared/InformedFiniteMixtures/ss_results14/ari_mn_nobs_",nobs,
                #                "_K_",K,"_vtype_",v_type,"_A_",A,"_datatype_",datatype,".txt")
                "_K_",K,"_datatype_",datatype,"_U_",U,".txt")
  write.table(ari, file=dir, row.names=FALSE, col.names=TRUE)


  #dir <- paste0("~/Research/BYU/InformedFiniteMixtures/analysis/ss_results14/Kplus_salso_nobs_",nobs,
  dir <- paste0("/home/gpage299/shared/InformedFiniteMixtures/ss_results14/Kplus_salso_nobs_",nobs,
                #                "_K_",K,"_vtype_",v_type,"_A_",A,"_datatype_",datatype,".txt")
                "_K_",K,"_datatype_",datatype,"_U_",U,".txt")
  write.table(Kplus_salso, file=dir, row.names=FALSE, col.names=TRUE)

  #dir <- paste0("~/Research/BYU/InformedFiniteMixtures/analysis/ss_results14/ccp_diff_nobs_",nobs,
  dir <- paste0("/home/gpage299/shared/InformedFiniteMixtures/ss_results14/ccp_diff_nobs_",nobs,
                #                "_K_",K,"_vtype_",v_type,"_A_",A,"_datatype_",datatype,".txt")
                "_K_",K,"_datatype_",datatype,"_U_",U,".txt")
  write.table(ccp_diff, file=dir, row.names=FALSE, col.names=TRUE)

  #dir <- paste0("~/Research/BYU/InformedFiniteMixtures/analysis/ss_results14/extra_w_nobs_",nobs,
  dir <- paste0("/home/gpage299/shared/InformedFiniteMixtures/ss_results14/extra_w_nobs_",nobs,
                #                "_K_",K,"_vtype_",v_type,"_A_",A,"_datatype_",datatype,".txt")
                "_K_",K,"_datatype_",datatype,"_U_",U,".txt")
  write.table(q, file=dir, row.names=FALSE, col.names=TRUE)

  print(Kplus_var[1:d,])
#  print(Kplus_bias[1:d,])
#  print(density_mse[1:d,])

}






if(FALSE){

  #ss_results14, I am now recording two different metrics to get at the idea that there are
  #              "spurious" clusters.  This required that I altered the file that summarizes results

  files <- list.files("~/Research/BYU/InformedFiniteMixtures/analysis/ss_results14")
  dat <- data.frame()
  dat2 <- data.frame()
  for(i in 1:length(files)){

     tmp <- read.table(paste0("~/Research/BYU/InformedFiniteMixtures/analysis/ss_results14/",files[i]), header=TRUE)
     print(files[i])
     print(dim(tmp[apply(tmp, 1, function(x){sum(is.na(x))<5}),]))


    info <- strsplit(files[i], "\\_|\\.txt")[[1]]
#    print(info)
    metric <- paste(info[1], info[2], sep="_")
    nobs <- as.numeric(info[4])
    Kd <- as.numeric(info[6])
    datatype <- as.numeric(info[8])
    U <- as.numeric(info[10])
    ndata <- nrow(tmp)

    dat <- rbind(dat, data.frame(metric=rep(metric,ncol(tmp)),
                                 nobs=rep(nobs, ncol(tmp)),
                                 Kd=rep(Kd, ncol(tmp)),
                                 datatype=rep(datatype,ncol(tmp)),
                                 U = rep(U,ncol(tmp)),
                                 proc=colnames(tmp),
                                 val=apply(tmp,2,mean,na.rm=TRUE)))

    dat2 <- rbind(dat2, data.frame(dataset=rep(1:ndata, ncol(tmp)),
                                   metric=rep(metric,ncol(tmp)*ndata),
                                   nobs=rep(nobs,ncol(tmp)*ndata),
                                   Kd=rep(Kd, ncol(tmp)*ndata),
                                   datatype=rep(datatype,ncol(tmp)*ndata),
                                   U = rep(U, ncol(tmp)*ndata),
                                   proc=rep(colnames(tmp), each=ndata),
                                   val=unlist(tmp)))


  }

  #ggplot
  dat$Klabel <- "K == 15"
  dat$Klabel[dat$Kd==25] <- "K == 25"
  dat$datlabel <- 'U == 1'
  dat$datlabel[dat$U==2] <- 'U == 2'
  dat$datlabel[dat$U==5] <- 'U == 5'
  dat$datlabel[dat$U==10] <- 'U == 10'
  dat$datlabel <- factor(dat$datlabel, levels=c('U == 1', 'U == 2','U == 5',
                                                'U == 10'))
  dat2$Proc <- factor(dat2$proc, levels=c("DPM","FM_Raffaele","FM_RevJump","MBC","Sylvia_Sparse",
                                          "U_1_Gam","U_1_Gam_alpha1F","U_1_Gam_alpha2F","U_1_Gam_ProbKp_01","U_1_Gam_ProbKp_09","U_1_PC_alpha2F_01","U_1_PC_alpha2F_09",
                                          "U_2_Gam","U_2_Gam_alpha1F","U_2_Gam_alpha2F","U_2_Gam_ProbKp_01","U_2_Gam_ProbKp_09","U_2_PC_alpha2F_01","U_2_PC_alpha2F_09",
                                          "U_5_Gam","U_5_Gam_alpha1F","U_5_Gam_alpha2F","U_5_Gam_ProbKp_01","U_5_Gam_ProbKp_09","U_5_PC_alpha2F_01","U_5_PC_alpha2F_09",
                                          "U_10_Gam","U_10_Gam_alpha1F","U_10_Gam_alpha2F","U_10_Gam_ProbKp_01","U_10_Gam_ProbKp_09","U_10_PC_alpha2F_01","U_10_PC_alpha2F_09"),
                                 labels=c("DPM","NormIFFP","FMM","MBC","sFMM",
                                        "AFM_U1_Gam","AFM_U1_Gam_a1F","AFM_U1_Gam_a2F","AFM_U1_QGam01_a2F","AFM_U1_QGam09_a2F","AFM_U1_PC01_a2F","AFM_U1_PC09_a2F",
                                        "AFM_U2_Gam","AFM_U2_Gam_a1F","AFM_U2_Gam_a2F","AFM_U2_QGam01_a2F","AFM_U2_QGam09_a2F","AFM_U2_PC01_a2F","AFM_U2_PC09_a2F",
                                        "AFM_U5_Gam","AFM_U5_Gam_a1F","AFM_U5_Gam_a2F","AFM_U5_QGam01_a2F","AFM_U5_QGam09_a2F","AFM_U5_PC01_a2F","AFM_U5_PC09_a2F",
                                        "AFM_U10_Gam","AFM_U10_Gam_a1F","AFM_U10_Gam_a2F","AFM_U10_QGam01_a2F","AFM_U10_QGam09_a2F","AFM_U10_PC01_a2F","AFM_U10_PC09_a2F"
                                         ))

  dat2$nobs_label <- 'n == 100'
  dat2$nobs_label[dat2$nobs == 1000] <- 'n==1000'
  dat2$datatype_label <- "Data~Type~1"
  dat2$datatype_label[dat$datatype==2] <- "Data~Type~2"
  dat2$Ulabel <- "K[true]^'+' == 1"
  dat2$Ulabel[dat2$U==2] <- "K[true]^'+' == 2"
  dat2$Ulabel[dat2$U==5] <- "K[true]^'+' == 5"
  dat2$Ulabel[dat2$U==10] <- "K[true]^'+' == 10"
  dat2$Ulabel <- factor(dat2$Ulabel, levels=c("K[true]^'+' == 1", "K[true]^'+' == 2",
                                              "K[true]^'+' == 5", "K[true]^'+' == 10"))

  tmp <- strsplit(as.character(dat2$Proc), "\\_|U")
  dat2$Proc2 <- as.character(dat2$Proc)
  dat2$Proc2[sapply(tmp, function(x) x[1]=="AFM")] <- sapply(tmp, function(x) paste0(x[1], "U",x[3]))[sapply(tmp, function(x) x[1]=="AFM")]
  dat2$Proc2 <- as.factor(dat2$Proc2)
  dat2$Proc2 <- factor(dat2$Proc2, levels=c("DPM","NormIFFP","FMM","sFMM","MBC","AFMU1","AFMU2","AFMU5","AFMU10"),
                                   labels=c("DPM","NormIFFP","FMM","sFMM","MBC","aFMMU1","aFMMU2","aFMMU5","aFMMU10"))
  dat2$Ufit <- as.numeric(sapply(tmp,function(x) x[3]))
  dat2$alphaPrior <- as.character(dat2$Proc)
  dat2$alphaPrior[sapply(tmp, function(x) x[1]=="AFM")] <- sapply(tmp, function(x) paste0(x[4],"_",x[5]))[sapply(tmp, function(x) x[1]=="AFM")]
  dat2$alphaPrior <- factor(dat2$alphaPrior, levels=c("DPM","NormIFFP","FMM","sFMM","MBC",
                                                      "Gam_a1F","Gam_a2F","Gam_NA","QGam01_a2F","QGam09_a2F","PC01_a2F","PC09_a2F"),
                                             labels=c("", "", "", "", "",
                                                      "Gam_a1F","Gam_a2F","Gam_NA","QGam01_a2F","QGam09_a2F","PC01_a2F","PC09_a2F"))



  # plots found in the article. Note that I am not including the majority of the model fits
  # as the plots were just to messy

  library(ggplot2)
  library(ggpubr)


  dat_tmp <- dat2[dat2$alphaPrior %in% c("", "Gam_a2F","QGam01_a2F","PC01_a2F", "PC09_a2F") &
                    dat2$Proc %in% c("DPM","NormIFFP","FMM","sFMM",
                                     "AFM_U1_Gam_a2F","AFM_U1_QGam01_a2F","AFM_U1_PC01_a2F","AFM_U1_QGam09_a2F","AFM_U1_PC09_a2F",
                                     "AFM_U2_Gam_a2F","AFM_U2_QGam01_a2F","AFM_U2_PC01_a2F","AFM_U2_QGam09_a2F","AFM_U2_PC09_a2F",
                                     "AFM_U5_Gam_a2F","AFM_U5_QGam01_a2F","AFM_U5_PC01_a2F","AFM_U5_QGam09_a2F","AFM_U5_PC09_a2F",
                                     "AFM_U10_Gam_a2F","AFM_U10_QGam01_a2F","AFM_U10_PC01_a2F","AFM_U10_QGam09_a2F","AFM_U10_PC09_a2F"),]

  U1_proc <- c("U_1_Gam_alpha2F", "U_1_Gam_ProbKp_01", "U_1_PC_alpha2F_01", "U_1_PC_alpha2F_09")
  QPCgam_prior <- c("QGam01_a2F", "PC09_a2F")






  # These include PC tp 0.9, but not QGam prior specification
  ggp1 <- ggplot(data = dat_tmp[dat_tmp$metric=="Kplus_bias" &
                                  dat_tmp$datatype==1 & dat_tmp$U!=1 &
                                  !dat_tmp$proc %in% U1_proc &
                                  dat_tmp$alphaPrior != "QGam01_a2F", ],
                 aes(x=as.factor(alphaPrior), y=val, col=as.factor(Proc2))) +
    geom_boxplot(position = position_dodge(preserve = "total")) +
    facet_grid(Ulabel~nobs_label, label = "label_parsed", scales="free") +
    theme_bw(base_size=15) + geom_hline(yintercept=0, linewidth=1, linetype="dotted") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("", "Gam","PC(0.1)","PC(0.9)")) +
    labs(x="", y ="Bias", col="Procedure")

  ggp2 <- ggplot(data = dat_tmp[dat_tmp$metric=="Kplus_bias" &
                                  dat_tmp$datatype==2 & dat_tmp$U!=1 &
                                  !dat_tmp$proc %in% U1_proc &
                                  dat_tmp$alphaPrior != "QGam01_a2F", ],
                 aes(x=as.factor(alphaPrior), y=val, col=as.factor(Proc2))) +
    geom_boxplot(position = position_dodge(preserve = "total")) +
    facet_grid(Ulabel~nobs_label, label = "label_parsed", scales="free") +
    theme_bw(base_size=15) + geom_hline(yintercept=0, linewidth=1, linetype="dotted") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("", "Gam","PC(0.1)","PC(0.9)")) +
    labs(x="", y ="Bias", col="Procedure")


  ggp3 <- ggplot(data = dat_tmp[dat_tmp$metric=="Kplus_var" &
                                  dat_tmp$datatype==1 & dat_tmp$U!=1 &
                                  !dat_tmp$proc %in% U1_proc &
                                  dat_tmp$alphaPrior != "QGam01_a2F", ],
                 aes(x=as.factor(alphaPrior), y=log(val+1), col=as.factor(Proc2))) +
    geom_boxplot(position = position_dodge(preserve = "total")) +
    facet_grid(Ulabel~nobs_label, label = "label_parsed", scales="free") +
    theme_bw(base_size=15)  +
    scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("", "Gam","PC(0.1)","PC(0.9)")) +
    labs(x="", y=expression(log~~pwss(K^"+")), col="Procedure")

  ggp4 <- ggplot(data = dat_tmp[dat_tmp$metric=="Kplus_var" &
                                  dat_tmp$datatype==2 & dat_tmp$U!=1 &
                                  !dat_tmp$proc %in% U1_proc &
                                  dat_tmp$alphaPrior != "QGam01_a2F", ],
                 aes(x=as.factor(alphaPrior), y=log(val+1), col=as.factor(Proc2))) +
    geom_boxplot(position = position_dodge(preserve = "total")) +
    facet_grid(Ulabel~nobs_label, label = "label_parsed", scales="free") +
    theme_bw(base_size=15)  +
    scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("", "Gam","PC(0.1)","PC(0.9)")) +
    labs(x="", y =expression(log~~pwss(K^"+")), col="Procedure")

  ggp5 <- ggplot(data = dat_tmp[dat_tmp$metric=="ari_mn" &
                                  dat_tmp$datatype==1 & dat_tmp$U!=1 &
                                  !dat_tmp$proc %in% U1_proc &
                                  dat_tmp$alphaPrior != "QGam01_a2F", ],
                 aes(x=as.factor(alphaPrior), y=val, col=as.factor(Proc2))) +
    geom_boxplot(position = position_dodge(preserve = "total")) +
    facet_grid(Ulabel~nobs_label, label = "label_parsed", scales="free") +
    theme_bw(base_size=15)  +
    scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("", "Gam","PC(0.1)","PC(0.9)")) +
    labs(x="", y ="Adjusted Rand Index", col="Procedure")

  ggp6 <- ggplot(data = dat_tmp[dat_tmp$metric=="ari_mn" &
                                  dat_tmp$datatype==2 & dat_tmp$U!=1 &
                                  !dat_tmp$proc %in% U1_proc &
                                  dat_tmp$alphaPrior != "QGam01_a2F", ],
                 aes(x=as.factor(alphaPrior), y=val, col=as.factor(Proc2))) +
    geom_boxplot(position = position_dodge(preserve = "total")) +
    facet_grid(Ulabel~nobs_label, label = "label_parsed", scales="free") +
    theme_bw(base_size=15)  +
    scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("", "Gam","PC(0.1)","PC(0.9)")) +
    labs(x="", y ="Adjusted Rand Index", col="Procedure")

  ggp7 <- ggplot(data = dat_tmp[dat_tmp$metric=="Kplus_salso" &
                                  dat_tmp$datatype==1 & dat_tmp$U!=1 &
                                  !dat_tmp$proc %in% U1_proc &
                                  dat_tmp$alphaPrior != "QGam01_a2F", ],
                 aes(x=as.factor(alphaPrior), y=val, col=as.factor(Proc2))) +
    geom_boxplot(position = position_dodge(preserve = "total")) +
    facet_grid(Ulabel~nobs_label, label = "label_parsed", scales="free") +
    theme_bw(base_size=15)  +
    scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("", "Gam","PC(0.1)","PC(0.9)")) +
    labs(x="", y =expression(K^"+"~~"salso"), col="Procedure")

  ggp8 <- ggplot(data = dat_tmp[dat_tmp$metric=="Kplus_salso" &
                                  dat_tmp$datatype==2 & dat_tmp$U!=1 &
                                  !dat_tmp$proc %in% U1_proc &
                                  dat_tmp$alphaPrior != "QGam01_a2F", ],
                 aes(x=as.factor(alphaPrior), y=val, col=as.factor(Proc2))) +
    geom_boxplot(position = position_dodge(preserve = "total")) +
    facet_grid(Ulabel~nobs_label, label = "label_parsed", scales="free") +
    theme_bw(base_size=15)  +
    scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("", "Gam","PC(0.1)","PC(0.9)")) +
    labs(x="", y =expression(K^"+"~~"salso"), col="Procedure")

  ggp9 <- ggplot(data = dat_tmp[dat_tmp$metric=="ccp_diff" &
                                  dat_tmp$datatype==1 & dat_tmp$U!=1 &
                                  !dat_tmp$proc %in% U1_proc &
                                  dat_tmp$alphaPrior != "QGam01_a2F", ],
                 aes(x=as.factor(alphaPrior), y=val, col=as.factor(Proc2))) +
    geom_boxplot(position = position_dodge(preserve = "total")) +
    facet_grid(Ulabel~nobs_label, label = "label_parsed", scales="free") +
    theme_bw(base_size=15)  +
    scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("", "Gam","PC(0.1)","PC(0.9)")) +
    labs(x="", y ="ccprob_error", col="Procedure")


  ggp10 <- ggplot(data = dat_tmp[dat_tmp$metric=="ccp_diff" &
                                   dat_tmp$datatype==2 & dat_tmp$U!=1 &
                                   !dat_tmp$proc %in% U1_proc &
                                   dat_tmp$alphaPrior != "QGam01_a2F", ],
                  aes(x=as.factor(alphaPrior), y=val, col=as.factor(Proc2))) +
    geom_boxplot(position = position_dodge(preserve = "total")) +
    facet_grid(Ulabel~nobs_label, label = "label_parsed", scales="free") +
    theme_bw(base_size=15)  +
    scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("", "Gam","PC(0.1)","PC(0.9)")) +
    labs(x="", y ="ccprob_error", col="Procedure")



  ggp11 <- ggplot(data = dat_tmp[dat_tmp$metric=="extra_w" &
                                  dat_tmp$datatype==1 & dat_tmp$U!=1 &
                                  !dat_tmp$proc %in% U1_proc &
                                  dat_tmp$alphaPrior != "QGam01_a2F", ],
                 aes(x=as.factor(alphaPrior), y=val, col=as.factor(Proc2))) +
    geom_boxplot(position = position_dodge(preserve = "total")) +
    facet_grid(Ulabel~nobs_label, label = "label_parsed", scales="free") +
    theme_bw(base_size=15)  +
    scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("", "Gam","PC(0.1)","PC(0.9)")) +
    labs(x="", y ="sum extra weights", col="Procedure")


  ggp12 <- ggplot(data = dat_tmp[dat_tmp$metric=="extra_w" &
                                   dat_tmp$datatype==2 & dat_tmp$U!=1 &
                                   !dat_tmp$proc %in% U1_proc &
                                   dat_tmp$alphaPrior != "QGam01_a2F", ],
                  aes(x=as.factor(alphaPrior), y=val, col=as.factor(Proc2))) +
    geom_boxplot(position = position_dodge(preserve = "total")) +
    facet_grid(Ulabel~nobs_label, label = "label_parsed", scales="free") +
    theme_bw(base_size=15)  +
    scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("", "Gam","PC(0.1)","PC(0.9)")) +
    labs(x="", y ="sum extra weights", col="Procedure")


  pdf("~/Research/BYU/InformedFiniteMixtures/latex/figures/asymSimStudyResults_PC01_09_2.pdf",
      width=9, height=9)
  ggp1
  ggp2
  ggp3
  ggp4
  ggp5
  ggp6
  ggp7
  ggp8
  ggp9
  ggp10
  dev.off()

