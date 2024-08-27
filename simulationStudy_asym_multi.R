#!/usr/bin/env Rscript

# Simulation study for curves.
# Generate data similar to the biomechanic data
# Use the pc prior on alpha1 | alpha2.
#
# We will compare methods based on the sparse FMM of Silvia et al
#
# Our comparison will be carried out using five metrics.  They are as follows
#
# 1. true K+ - mode of posterior of K+ (this measures bias)
# 2. sum_{k-1}^K {(U - k)^2*Pr(K+ = k | data)} (this measures accuracy of entire Post of K+ )
# 3. true K+ - K+ estimated from salso partition estimate
# 4. sum_i sum_{l<i} ((I[i~l] - Pr(zi=zl|y))^2
# 5. ARI between estimated and true partition
#
# To match the application, we will generate data using two the vgrf data
#


rm(list=ls())

#install.packages("TruncatedDistributions", repos="http://R-Forge.R-project.org")
library(miscPack) # needed to fit our models version
library(TruncatedDistributions) # needed for the truncated exponential
library(MCMCpack) # needed for the rdirichlet function
library(mixAK)    # needed for NMixMCMC function (RJMCMC finite mixture)
library(BNPmix)   # needed for PYdensity function
library(mclust)   # needed for Mclust and densityMclust functions
library(AntMAN)   # needed for the AM_mcmc_fit function
library(salso)    # needed for salso estimation method
library(mvtnorm)  # needed for rmvnorm function

# install.packages("MCMCpack","mixAK","BNPmix","mclust","AntMAN","salso")

library(splines)
library(spam)
library(ppmSuite)
library(mclust)
library(openxlsx)


sessionInfo()

# set ndata in 100
# set nobs in {100} - 25 per group approximately
# set K in {25}
# set U in {4}
# set datatype in {1}
#       1 - use vgrf data as guide
# set number of knots in {7, 10, 20}
# set Akap in {0.1, 0.5, 1.0}
#
# There are 1 x 1 x 1 x 1 x 3 x 3 = 9  scenarios
#
# nohup ./SimulationStudy_asym.R 25 100 25 1 1 1> out1.txt &

args <- commandArgs(TRUE)
# args <- c(100, 100,  25, 4, 1, 7, 0.5)
print(args)
ndata <- as.numeric(args[1])
nobs <- as.numeric(args[2])
K <- as.numeric(args[3])
Utrue <- as.numeric(args[4])
datatype <- as.numeric(args[5])
ndx <- as.numeric(args[6])
Akap <- as.numeric(args[7])

cat("ndata = ", ndata, "\n")
cat("nobs = ", nobs, "\n")
cat("K = ", K, "\n")
cat("Utrue = ", Utrue, "\n")
cat("datatype = ", datatype, "\n")
cat("ndx = ", ndx, "\n")
cat("Akap = ", Akap, "\n")

ntime <- 100



nsamples <- 1e5 # monte carlo samples to get lambda and alpha density
agrid_length <- 1e5 # denseness of values of alpha

agrid <- exp(seq(log(1e-5), log(1), length.out = agrid_length))


Uprior <- c(2,4,6,8)

Kplus_var <- matrix(NA, nrow=ndata, ncol=2*length(Uprior)+3)
Kplus_bias <- matrix(NA, nrow=ndata, ncol=2*length(Uprior)+3)
ari <- matrix(NA, nrow=ndata, ncol=2*length(Uprior)+3)
Kplus_salso <- matrix(NA, nrow=ndata, ncol=2*length(Uprior)+3)
ccp_diff <- matrix(NA, nrow=ndata, ncol=2*length(Uprior)+3)
extra_w <- matrix(NA, nrow=ndata, ncol=2*length(Uprior)+3)

colnames(Kplus_var) <- colnames(Kplus_bias)  <- colnames(ari) <-
colnames(Kplus_salso) <- colnames(ccp_diff) <- colnames(extra_w) <-
  c(apply(expand.grid(Uprior,c("PC_alpha2F_01")),1,function(x) paste0("U_",as.numeric(x[1]),"_",as.character(x[2]))),
    apply(expand.grid(Uprior,c("PC_alpha2F_09")),1,function(x) paste0("U_",as.numeric(x[1]),"_",as.character(x[2]))),
    "sFMM", "DPM", "FMM")

ndens_y <- 1000

niter <- 50000
nburn <- 40000
nthin <- 10
nout <- (niter - nburn)/nthin



# d <- 1; ndata=5; U <- 3
for(d in 1:ndata){

  cat("d = ", d, "\n")

  set.seed(ndata*(datatype-1) + d)
  cat("seed = ", ndata*(datatype-1) + d, "\n")

  print(date())

  if(datatype == 1){
    # Try using Matt's VGRF data to create curves using B-spline coefficients
#    vgrf <- read.xlsx("~/Research/BYU/BiomechanicCurves/Matt_ACLpatients/data2_25_2023/vgrf.xlsx", startRow=1, colNames=FALSE);
    vgrf <- read.xlsx("~/shared/InformedFiniteMixtures/vgrf.xlsx", startRow=1, colNames=FALSE);
    vgrf <- sapply(vgrf, as.numeric)
    nsub <- (ncol(vgrf)/5)
    ntime <- nrow(vgrf)

    vgrf_mn <- NULL
    subject <- rep(1:nsub, each=5)
    for(i in 1:nsub){
      vgrf_mn <- cbind(vgrf_mn, apply(vgrf[,subject==i], 1, mean, na.rm=TRUE))
    }

    # Fit multivariate mixture to data
    mod1 <- Mclust(t(vgrf_mn), modelNames="EEI", G=4)
    clus_lab <- apply(mod1$z, 1, which.max)

    # compute cluster-means to be used to generate data
    vgrf_clus_mn <- NULL
    for(i in 1:4){
      vgrf_clus_mn <- cbind(vgrf_clus_mn, apply(vgrf_mn[,clus_lab==i], 1, mean, na.rm=TRUE))
    }


    ord <- 3
    beta_dim <- 30
    x <- seq(-3, 3, length=ntime)
    B <- miscPack:::bspline(x=x, xl=range(x)[1], xr=range(x)[2],
                    ndx=beta_dim-3, bdeg=3)

    bbeta <- solve(t(B) %*% B) %*% t(B) %*% vgrf_clus_mn


    zi <- sample(1:4, size=nobs, replace=TRUE)
    yi <- NULL
    for(i in 1:nobs){
      yi <- cbind(yi, rnorm(ntime, B %*% (bbeta[,zi[i]]+rnorm(beta_dim,0,0.03)), 0.01))
    }
  }

  if(datatype == 2){
    # Try using Matt's KFA data to create curves using B-spline coefficients
#    kfa <- read.xlsx("~/Research/BYU/BiomechanicCurves/Matt_ACLpatients/data2_25_2023/kfa.xlsx", startRow=1, colNames=FALSE);
    kfa <- read.xlsx("~/shared/InformedFiniteMixtures/kfa.xlsx", startRow=1, colNames=FALSE);
    kfa <- sapply(kfa, as.numeric)
    nsub <- (ncol(kfa)/5)
    ntime <- nrow(kfa)

    kfa_mn <- NULL
    subject <- rep(1:nsub, each=5)
    for(i in 1:nsub){
      kfa_mn <- cbind(kfa_mn, apply(kfa[,subject==i], 1, mean, na.rm=TRUE))
    }

    # Fit multivariate mixture to data
    mod1 <- Mclust(t(kfa_mn), modelNames="EEI", G=Utrue)
    clus_lab <- apply(mod1$z, 1, which.max)

    # compute cluster-means to be used to generate data
    kfa_clus_mn <- NULL
    for(i in 1:4){
      kfa_clus_mn <- cbind(kfa_clus_mn, apply(kfa_mn[,clus_lab==i], 1, mean, na.rm=TRUE))
    }


    ord <- 3
    beta_dim <- 30
    x <- seq(-3, 3, length=ntime)
    B <- miscPack:::bspline(x=x, xl=range(x)[1], xr=range(x)[2],
                            ndx=beta_dim-3, bdeg=3)

    bbeta <- solve(t(B) %*% B) %*% t(B) %*% kfa_clus_mn


    zi <- sample(1:4, size=nobs, replace=TRUE)
    yi <- NULL
    for(i in 1:nobs){
      yi <- cbind(yi, rnorm(ntime, B %*% (bbeta[,zi[i]]+rnorm(beta_dim,0,2)), 0.05))
    }
  }


  if(datatype == 3){
    # generate data
    ord <- 2
    b_dim <- 10
    sigma2 <- 0.5
    s2 <- 0.5
    sigma2.err <- 0.001 # 0.05

    x <- seq(-3, 3, length=ntime)
    B <- bspline(x=x, xl=range(x)[1], xr=range(x)[2],
               ndx=b_dim-3, bdeg=3)

    D = diff(diag(b_dim), diff=ord)
    R = crossprod(D,D)
    eigendec <- eigen(R)
    ev <- eigendec$values
    V <- eigendec$vectors

    theta.ulc <- matrix(NA, nrow=ncol(B), ncol=Utrue)
    betai <- matrix(NA, nrow=ncol(B), ncol=nobs)
    ind <- c(rep(T, b_dim-2), rep(F,2))
    for (k in 1:Utrue){
      znorm <- rnorm(sum(ind))

      # Generate cluster specific B-splines coefficients
      theta.ulc[,k] <- (V[,ind]%*%diag(sqrt(1/ev[ind]))) %*% znorm  #Cluster B-spline coefficients
    }
    mus <- rnorm(Utrue,0,0.25)
    zi <- sample(1:Utrue, nobs, replace=TRUE)
    yi <- matrix(NA, nrow=ntime, ncol=nobs)
    kappa <- 0.5
    for(i in 1:nobs){
  	  mu <- rnorm(1,0,0.25)
      beta <-0
  	  betai[,i] <- theta.ulc[,zi[i]] + rnorm(ncol(B), 0, kappa)
  	  yi[,i] <- rnorm(ntime, mean=mus[zi[i]] +  sqrt(sigma2)*beta*x + sqrt(s2)*B %*% betai[,i], sd=sqrt(sigma2.err))

    }
  }



  Rho_true <- miscPack:::relabel(zi)
  ppt <- matrix(0, ncol=nobs, nrow=nobs)
  for(ii in 1:nobs){
    for(iii in 1:nobs){
      if(zi[ii] == zi[iii]) ppt[ii,iii] <- 1
    }
  }

  y <- c(yi)
  t <- rep(1:nrow(yi), nobs)
  ids <- rep(1:nobs, each=ntime)

  min_y <- min(y)
  ran_y <- diff(range(y))

  # Rescale so that each curve is on the same scale
  rescale <- function(x) (x - min(x))/diff(range(x))
  y_scaled <- (y-min_y)/ran_y


  A <- 0.001

  alpha2_val_fixed <- 1e-5

  fits1 <- fits2 <- list()

  for(u in 1:length(Uprior)){
    cat("u = ", u, "\n")
    cat("nobs = ", nobs, "\n")

    cat("fit model 1 \n")
    # fits the PC prior with alpha1 unknown and alpha2 fixed with tail probability 0.1
    fits1[[u]] <-pspline_mixture(y=y_scaled, t=t, ids=ids, K=K,
                                 alpha_prior_type="centered",
                                 alpha_prior_dist="pc",
                                 basemodel=0,
                                 U=Uprior[u],
                                 #mylambda = lambda01[u],
                                 alpha1_val = 1,
                                 alpha2_val = 1e-5,
                                 update_alpha1=TRUE,
                                 update_alpha2=FALSE,
                                 tail.prob = 1e-1,
                                 ndx = ndx, q = 3,
                                 ndens_y=1000,
                                 A=A, Akap=Akap,
                                 U_tau= 1/(1*0.31), a_tau=1e-2, U_omega=1/(1*0.31), a_omega=1e-2,
                                 a_gam=1,
                                 niter=niter, nburn=nburn, nthin=nthin)

    kpp1 <- table(fits1[[u]]$kp)/nout
    Kplus_var[d,u] <- sum(((Utrue - as.numeric(names(kpp1)))^2*kpp1))
    Kplus_bias[d,u] <- as.numeric(names(
      which.max(
        table(apply(fits1[[u]]$z, 1, function(x) length(unique(x))))
      )
     )
    ) - Utrue #posterior mode of K+.  Sylvia suggests this
    ari[d,u] <- mean(apply(fits1[[u]]$z, 1, function(x) adjustedRandIndex(Rho_true, x)))
    if(mean(apply(fits1[[u]]$z,1,sum)) == nobs){
      Kplus_salso[d,u] <- 1 - Utrue
    } else {
      Kplus_salso[d,u] <- length(table(salso(fits1[[u]]$z))) - Utrue
    }
    ccp_diff[d,u] <- mean(((ppt - psm(fits1[[u]]$z))[upper.tri(ppt, diag=FALSE)])^2)
    extra_w[d, u] <- mean(sapply(1:nout, function(x) sum(sort(fits1[[u]]$w[x,], decreasing=TRUE)[-(1:fits1[[u]]$kp[x])])))



    # fits the PC prior with alpha1 unknown and alpha2 fixed and tail probability equal to 0.9
    fits2[[u]] <- pspline_mixture(y=y_scaled, t=t, ids=ids, K=K,
                                  alpha_prior_type="centered",
                                  alpha_prior_dist="pc",
                                  basemodel=0,
                                  U=Uprior[u],
                                  #mylambda = lambda09[U],
                                  alpha1_val = 1,
                                  alpha2_val = 1e-5,
                                  update_alpha1=TRUE,
                                  update_alpha2=FALSE,
                                  tail.prob = 0.9,
                                  ndx = ndx, q = 3,
                                  ndens_y=1000,
                                  A=A, Akap=Akap,
                                  U_tau= 1/(1*0.31), a_tau=1e-2, U_omega=1/(1*0.31), a_omega=1e-2,
                                  a_gam=1,
                                  niter=niter, nburn=nburn, nthin=nthin)

    kpp1 <- table(apply(fits2[[u]]$z, 1, function(x) length(unique(x))))/nout
    Kplus_var[d,u+length(Uprior)] <- sum(((Utrue - as.numeric(names(kpp1)))^2*kpp1))
    Kplus_bias[d,u+length(Uprior)] <- as.numeric(names(
      which.max(
        table(apply(fits2[[u]]$z, 1, function(x) length(unique(x))))
      )
    )
    ) - Utrue #posterior mode of K+.  Sylvia suggests this
    ari[d,u+length(Uprior)] <- mean(apply(fits2[[u]]$z, 1, function(x) adjustedRandIndex(Rho_true, x)))
    if(mean(apply(fits2[[u]]$z,1,sum)) == nobs){
      bdif <- 1 - Utrue
    } else {
      bdif <- length(table(salso(fits2[[u]]$z))) - Utrue
    }
    Kplus_salso[d,u+length(Uprior)] <- bdif
    ccp_diff[d,u+length(Uprior)] <- mean(((ppt - psm(fits2[[u]]$z))[upper.tri(ppt, diag=FALSE)])^2)
    extra_w[d, u+length(Uprior)] <- mean(sapply(1:nout, function(x) sum(sort(fits2[[u]]$w[x,], decreasing=TRUE)[-(1:fits2[[u]]$kp[x])])))



  }

  # fit Sylvia model i.e., the sparse
  cat("Fit Sylvia's model", "\n")
  pc3 <- pspline_mixture(y=y_scaled, t=t, ids=ids, K=K,
                         alpha_prior_type="sparse",
                         alpha_prior_dist="gamma",
                         ndx = ndx, q = 3,
                         ndens_y=1000,
                         A=A, Akap=Akap,
                         a_gam=10, b_gam=1/(10*K),
                         niter=niter, nburn=nburn, nthin=nthin)

  kpp1 <- table(apply(pc3$z, 1, function(x) length(unique(x))))/nout
  Kplus_var[d,u+length(Uprior)+1] <- sum(((Utrue - as.numeric(names(kpp1)))^2*kpp1))
  Kplus_bias[d,u+length(Uprior)+1] <- as.numeric(names(
    which.max(
      table(apply(pc3$z, 1, function(x) length(unique(x))))
    )
   )
  ) - Utrue #posterior mode of K+.  Sylvia suggests this
  ari[d,u+length(Uprior)+1] <- mean(apply(pc3$z, 1, function(x) adjustedRandIndex(Rho_true, x)))
  Kplus_salso[d,u+length(Uprior)+1] <- length(table(salso(pc3$z))) - Utrue
  ccp_diff[d,u+length(Uprior)+1] <- mean(((ppt - psm(pc3$z))[upper.tri(ppt, diag=FALSE)])^2)
  extra_w[d, u+length(Uprior)+1] <- mean(sapply(1:nout, function(x) sum(sort(pc3$w[x,], decreasing=TRUE)[-(1:pc3$kp[x])])))


  # DPM
  dp1 <- curve_ppmx(y=cbind(y_scaled), z=t, subject=ids, PPM=TRUE, M=1,
                     nknots=ndx, npredobs=0, Aparm=Akap,
                     modelPriors=c(0.1, 100, 0, 100, 1, 1,1 ,1),
                     consim=1, simParms=c(0,1,1,1,1,1,1),
                     calibrate=0)
  kpp1 <- table(apply(dp1$Si, 1, function(x) length(unique(x))))/nout
  Kplus_var[d,u+length(Uprior)+2] <- sum(((Utrue - as.numeric(names(kpp1)))^2*kpp1))
  Kplus_bias[d,u+length(Uprior)+2] <- as.numeric(names(
    which.max(
      table(apply(dp1$Si, 1, function(x) length(unique(x))))
    )
   )
  ) - Utrue #posterior mode of K+.  Sylvia suggests this
  ari[d,u+length(Uprior)+2] <- mean(apply(dp1$Si, 1, function(x) adjustedRandIndex(Rho_true, x)))
  Kplus_salso[d,u+length(Uprior)+2] <- length(table(salso(dp1$Si))) - Utrue
  ccp_diff[d,u+length(Uprior)+2] <- mean(((ppt - psm(dp1$Si))[upper.tri(ppt, diag=FALSE)])^2)
  extra_w[d, u+length(Uprior)+2] <- mean(sapply(1:nout, function(x){
                                                            nc <- length(unique(dp1$Si[x,]))
                                                            sum(sort(dp1$w[x,], decreasing=TRUE)[-(1:nc)])
                                                           })
                                           )


  # fit static FMM model
  cat("Fit FMM model", "\n")
  pc1 <- pspline_mixture(y=y_scaled, t=t, ids=ids, K=K,
                         alpha_prior_type="centered",
                         alpha_prior_dist="gamma",
                         alpha1_val = 1/K,
                         alpha2_val = 1/K,
                         update_alpha1=FALSE,
                         update_alpha2=FALSE,
                         ndx = ndx, q = 3,
                         ndens_y=1000,
                         A=A, Akap=Akap,
                         a_gam=1,
                         niter=niter, nburn=nburn, nthin=nthin)

  kpp1 <- table(apply(pc1$z, 1, function(x) length(unique(x))))/nout
  Kplus_var[d,u+length(Uprior)+3] <- sum(((Utrue - as.numeric(names(kpp1)))^2*kpp1))
  Kplus_bias[d,u+length(Uprior)+3] <- as.numeric(names(
    which.max(
      table(apply(pc1$z, 1, function(x) length(unique(x))))
    )
  )
  ) - Utrue #posterior mode of K+.  Sylvia suggests this
  ari[d,u+length(Uprior)+3] <- mean(apply(pc1$z, 1, function(x) adjustedRandIndex(Rho_true, x)))
  Kplus_salso[d,u+length(Uprior)+3] <- length(table(salso(pc1$z))) - Utrue
  ccp_diff[d,u+length(Uprior)+3] <- mean(((ppt - psm(pc1$z))[upper.tri(ppt, diag=FALSE)])^2)
  extra_w[d, u+length(Uprior)+3] <- mean(sapply(1:nout, function(x) sum(sort(pc1$w[x,], decreasing=TRUE)[-(1:pc1$kp[x])])))


  #dir <- paste0("~/Research/BYU/InformedFiniteMixtures/analysis/ss_multi_results/Kplus_var_nobs_",nobs,
  dir <- paste0("~/shared/InformedFiniteMixtures/ss_multi_results/Kplus_var_nobs_",nobs,
                "_K_",K,"_datatype_",datatype,"_U_",Utrue,"_ndx_",ndx,"_Akap_",Akap,".txt")
  write.table(Kplus_var, file=dir, row.names=FALSE, col.names=TRUE)


  #dir <- paste0("~/Research/BYU/InformedFiniteMixtures/analysis/ss_multi_results/Kplus_bias_nobs_",nobs,
  dir <- paste0("~/shared/InformedFiniteMixtures/ss_multi_results/Kplus_bias_nobs_",nobs,
                #                "_K_",K,"_vtype_",v_type,"_A_",A,"_datatype_",datatype,".txt")
                "_K_",K,"_datatype_",datatype,"_U_",Utrue,"_ndx_",ndx,"_Akap_",Akap,".txt")
  write.table(Kplus_bias, file=dir, row.names=FALSE, col.names=TRUE)


  #dir <- paste0("~/Research/BYU/InformedFiniteMixtures/analysis/ss_multi_results/ari_mn_nobs_",nobs,
  dir <- paste0("~/shared/InformedFiniteMixtures/ss_multi_results/ari_mn_nobs_",nobs,
                #                "_K_",K,"_vtype_",v_type,"_A_",A,"_datatype_",datatype,".txt")
                "_K_",K,"_datatype_",datatype,"_U_",Utrue,"_ndx_",ndx,"_Akap_",Akap,".txt")
  write.table(ari, file=dir, row.names=FALSE, col.names=TRUE)


  #dir <- paste0("~/Research/BYU/InformedFiniteMixtures/analysis/ss_multi_results/Kplus_salso_nobs_",nobs,
  dir <- paste0("~/shared/InformedFiniteMixtures/ss_multi_results/Kplus_salso_nobs_",nobs,
                #                "_K_",K,"_vtype_",v_type,"_A_",A,"_datatype_",datatype,".txt")
                "_K_",K,"_datatype_",datatype,"_U_",Utrue,"_ndx_",ndx,"_Akap_",Akap,".txt")
  write.table(Kplus_salso, file=dir, row.names=FALSE, col.names=TRUE)

  #dir <- paste0("~/Research/BYU/InformedFiniteMixtures/analysis/ss_multi_results/ccp_diff_nobs_",nobs,
  dir <- paste0("~/shared/InformedFiniteMixtures/ss_multi_results/ccp_diff_nobs_",nobs,
                #                "_K_",K,"_vtype_",v_type,"_A_",A,"_datatype_",datatype,".txt")
                "_K_",K,"_datatype_",datatype,"_U_",Utrue,"_ndx_",ndx,"_Akap_",Akap,".txt")
  write.table(ccp_diff, file=dir, row.names=FALSE, col.names=TRUE)


  print(Kplus_var[1:d,])
#  print(Kplus_bias[1:d,])
#  print(density_mse[1:d,])

}


if(FALSE){

  #ss_multi_results, First simulation study with the "centered" prior.  This simulation
  #             includes U \in {1, 2, 5, 10} and many competitors Silvia, RJ FMM
  #             Raffaele, DPM, and Model based cluster.  Metrics included are
  #             ARI, bias of Kplus, Variance of Kplus, and Density estimate
  #             alpha1 ~ Gamma(1, 0.1*U), alpha2 ~ UN(1e-10, 1e-3)

  files <- list.files("~/Research/BYU/InformedFiniteMixtures/analysis/ss_multi_results")
  dat <- data.frame()
  dat2 <- data.frame()
  for(i in 1:length(files)){

     tmp <- read.table(paste0("~/Research/BYU/InformedFiniteMixtures/analysis/ss_multi_results/",files[i]), header=TRUE)
     print(files[i])
     print(dim(tmp[apply(tmp, 1, function(x){sum(is.na(x))<5}),]))


    info <- strsplit(files[i], "\\_|\\.txt")[[1]]
#    print(info)
    metric <- paste(info[1], info[2], sep="_")
    nobs <- as.numeric(info[4])
    Kd <- as.numeric(info[6])
    datatype <- as.numeric(info[8])
    U <- as.numeric(info[10])
    ndx <- as.numeric(info[12])
    Akap <- as.numeric(info[14])
    ndata <- nrow(tmp)

    dat <- rbind(dat, data.frame(metric=rep(metric,ncol(tmp)),
                                 nobs=rep(nobs, ncol(tmp)),
                                 Kd=rep(Kd, ncol(tmp)),
                                 datatype=rep(datatype,ncol(tmp)),
                                 U = rep(U,ncol(tmp)),
                                 ndx = rep(ndx, ncol(tmp)),
                                 Akap = rep(Akap, ncol(tmp)),
                                 proc=colnames(tmp),
                                 val=apply(tmp,2,mean,na.rm=TRUE)))

    dat2 <- rbind(dat2, data.frame(dataset=rep(1:ndata, ncol(tmp)),
                                   metric=rep(metric,ncol(tmp)*ndata),
                                   nobs=rep(nobs,ncol(tmp)*ndata),
                                   Kd=rep(Kd, ncol(tmp)*ndata),
                                   datatype=rep(datatype,ncol(tmp)*ndata),
                                   U = rep(U, ncol(tmp)*ndata),
                                   ndx = rep(ndx, ncol(tmp)*ndata),
                                   Akap = rep(Akap, ncol(tmp)*ndata),
                                   proc=rep(colnames(tmp), each=ndata),
                                   val=unlist(tmp)))


  }

  #ggplot
  dat$Klabel <- "K == 25"
  dat$datlabel <- 'U == 4'
  dat2$Proc <- factor(dat2$proc, levels=c("DPM","FMM","sFMM",
                                          "U_2_PC_alpha2F_01","U_2_PC_alpha2F_09",
                                          "U_4_PC_alpha2F_01","U_4_PC_alpha2F_09",
                                          "U_6_PC_alpha2F_01","U_6_PC_alpha2F_09",
                                          "U_8_PC_alpha2F_01","U_8_PC_alpha2F_09"),
                                 labels=c("DPM","FMM","sFMM",
                                        "AFM_U2_PC01_a2F","AFM_U2_PC09_a2F",
                                        "AFM_U4_PC01_a2F","AFM_U4_PC09_a2F",
                                        "AFM_U6_PC01_a2F","AFM_U6_PC09_a2F",
                                        "AFM_U8_PC01_a2F","AFM_U8_PC09_a2F"
                                         ))

  dat2$nobs_label <- 'n == 100'
  dat2$nobs_label[dat2$nobs == 1000] <- 'n==1000'
  dat2$datatype_label <- "Data~Type~1"
  dat2$datatype_label[dat$datatype==2] <- "Data~Type~2"
  dat2$datatype_label[dat$datatype==3] <- "Data~Type~3"
  dat2$Ulabel[dat2$U==4] <- "K[true]^'+' == 4"
  dat2$Ulabel <- factor(dat2$Ulabel, levels=c("K[true]^'+' == 4"))
  dat2$Akap[dat2$Akap == 0.999] <- 1.0
  dat2$Akap_label <- "A[0] == 0.1"
  dat2$Akap_label[dat2$Akap == 0.5] <- "A[0] == 0.5"
  dat2$Akap_label[dat2$Akap == 1.0] <- "A[0] == 1.0"
  dat2$ndx_label <- "p == 10"
  dat2$ndx_label[dat2$ndx == 10] <- "p == 13"
  dat2$ndx_label[dat2$ndx == 20] <- "p == 23"



  tmp <- strsplit(as.character(dat2$Proc), "\\_|U")
  dat2$Proc2 <- as.character(dat2$Proc)
  dat2$Proc2[sapply(tmp, function(x) x[1]=="AFM")] <- sapply(tmp, function(x) paste0(x[1], "U",x[3]))[sapply(tmp, function(x) x[1]=="AFM")]
  dat2$Proc2 <- as.factor(dat2$Proc2)
  dat2$Proc2 <- factor(dat2$Proc2, levels=c("DPM","FMM","sFMM","AFMU2","AFMU4","AFMU6","AFMU8"),
                                   labels=c("DPM","FMM","sFMM","aFMMU2","aFMMU4","aFMMU6","aFMMU8"))
  dat2$Ufit <- as.numeric(sapply(tmp,function(x) x[3]))
  dat2$alphaPrior <- as.character(dat2$Proc)
  dat2$alphaPrior[sapply(tmp, function(x) x[1]=="AFM")] <- sapply(tmp, function(x) paste0(x[4],"_",x[5]))[sapply(tmp, function(x) x[1]=="AFM")]
  dat2$alphaPrior <- factor(dat2$alphaPrior, levels=c("DPM","FMM","sFMM","PC01_a2F","PC09_a2F"),
                                             labels=c("", "", "","PC01_a2F","PC09_a2F"))



  # plots found in the article. Note that I am not including the majority of the model fits
  # as the plots were just to messy

  library(ggplot2)
  library(ggpubr)


  dat_tmp <- dat2[dat2$alphaPrior %in% c("", "PC01_a2F", "PC09_a2F") &
                    dat2$Proc %in% c("FMM","sFMM",
                                     "AFM_U2_PC01_a2F","AFM_U2_PC09_a2F",
                                     "AFM_U4_PC01_a2F","AFM_U4_PC09_a2F",
                                     "AFM_U6_PC01_a2F","AFM_U6_PC09_a2F",
                                     "AFM_U8_PC01_a2F","AFM_U8_PC09_a2F"),]

  U1_proc <- c("U_1_Gam_alpha2F", "U_1_Gam_ProbKp_01", "U_1_PC_alpha2F_01", "U_1_PC_alpha2F_09")
  QPCgam_prior <- c("QGam01_a2F", "PC09_a2F")






  # These data type 1 bias
  ggp1 <- ggplot(data = dat_tmp[dat_tmp$metric=="Kplus_bias" &
                                  dat_tmp$datatype==1 & dat_tmp$nobs==100, ],
                 aes(x=as.factor(alphaPrior), y=val, col=as.factor(Proc2))) +
    geom_boxplot(position = position_dodge(preserve = "total")) +
    facet_grid(ndx_label~Akap_label, label = "label_parsed", scales="free") +
    theme_bw(base_size=15) + geom_hline(yintercept=0, linewidth=1, linetype="dotted") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("", "PC(0.1)","PC(0.9)")) +
    labs(x="", y ="Bias", col="Procedure")


  # These data type 2 bias
  ggp2 <- ggplot(data = dat_tmp[dat_tmp$metric=="Kplus_bias" &
                                  dat_tmp$datatype==2 & dat_tmp$nobs==100, ],
                 aes(x=as.factor(alphaPrior), y=val, col=as.factor(Proc2))) +
    geom_boxplot(position = position_dodge(preserve = "total")) +
    facet_grid(ndx_label~Akap_label, label = "label_parsed", scales="free") +
    theme_bw(base_size=15) + geom_hline(yintercept=0, linewidth=1, linetype="dotted") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("", "PC(0.1)","PC(0.9)")) +
    labs(x="", y ="Bias", col="Procedure")

  # These data type 3 bias
  ggp3 <- ggplot(data = dat_tmp[dat_tmp$metric=="Kplus_bias" &
                                  dat_tmp$datatype==3 & dat_tmp$nobs==100, ],
                 aes(x=as.factor(alphaPrior), y=val, col=as.factor(Proc2))) +
    geom_boxplot(position = position_dodge(preserve = "total")) +
    facet_grid(ndx_label~Akap_label, label = "label_parsed", scales="free") +
    theme_bw(base_size=15) + geom_hline(yintercept=0, linewidth=1, linetype="dotted") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("", "PC(0.1)","PC(0.9)")) +
    labs(x="", y ="Bias", col="Procedure")



  pdf("~/Research/BYU/InformedFiniteMixtures/latex/figures/asymSimStudy_multi_Results_bias.pdf",
      width=9, height=9)
  ggp1
  ggp2
  ggp3
  dev.off()


  ggp4 <- ggplot(data = dat_tmp[dat_tmp$metric=="ari_mn" &
                                  dat_tmp$datatype==1 & dat_tmp$nobs==100 , ],
                 aes(x=as.factor(alphaPrior), y=val, col=as.factor(Proc2))) +
    geom_boxplot(position = position_dodge(preserve = "total")) +
    facet_grid(ndx_label~Akap_label, label = "label_parsed", scales="free") +
    theme_bw(base_size=15) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("", "PC(0.1)","PC(0.9)")) +
    labs(x="", y ="ARI", col="Procedure")

  ggp5 <- ggplot(data = dat_tmp[dat_tmp$metric=="ari_mn" &
                                  dat_tmp$datatype==2 & dat_tmp$nobs==100 , ],
                 aes(x=as.factor(alphaPrior), y=val, col=as.factor(Proc2))) +
    geom_boxplot(position = position_dodge(preserve = "total")) +
    facet_grid(ndx_label~Akap_label, label = "label_parsed", scales="free") +
    theme_bw(base_size=15) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("", "PC(0.1)","PC(0.9)")) +
    labs(x="", y ="ARI", col="Procedure")


  ggp6 <- ggplot(data = dat_tmp[dat_tmp$metric=="ari_mn" &
                                  dat_tmp$datatype==3 & dat_tmp$nobs==100 , ],
                 aes(x=as.factor(alphaPrior), y=val, col=as.factor(Proc2))) +
    geom_boxplot(position = position_dodge(preserve = "total")) +
    facet_grid(ndx_label~Akap_label, label = "label_parsed", scales="free") +
    theme_bw(base_size=15) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("", "PC(0.1)","PC(0.9)")) +
    labs(x="", y ="ARI", col="Procedure")

  pdf("~/Research/BYU/InformedFiniteMixtures/latex/figures/asymSimStudy_multi_Results_ari.pdf",
      width=9, height=9)
  ggp4
  ggp5
  ggp6
  dev.off()


  ggp7 <- ggplot(data = dat_tmp[dat_tmp$metric=="Kplus_var" &
                                  dat_tmp$datatype==1 & dat_tmp$nobs==100 , ],
                 aes(x=as.factor(alphaPrior), y=log(val+1), col=as.factor(Proc2))) +
    geom_boxplot(position = position_dodge(preserve = "total")) +
    facet_grid(ndx_label~Akap_label, label = "label_parsed", scales="free") +
    theme_bw(base_size=15) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("", "PC(0.1)","PC(0.9)")) +
    labs(x="", y=expression(log~~pwss(K^"+")), col="Procedure")

  ggp8 <- ggplot(data = dat_tmp[dat_tmp$metric=="Kplus_var" &
                                  dat_tmp$datatype==2 & dat_tmp$nobs==100 , ],
                 aes(x=as.factor(alphaPrior), y=log(val+1), col=as.factor(Proc2))) +
    geom_boxplot(position = position_dodge(preserve = "total")) +
    facet_grid(ndx_label~Akap_label, label = "label_parsed", scales="free") +
    theme_bw(base_size=15) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("", "PC(0.1)","PC(0.9)")) +
    labs(x="", y=expression(log~~pwss(K^"+")), col="Procedure")


  ggp9 <- ggplot(data = dat_tmp[dat_tmp$metric=="Kplus_var" &
                                  dat_tmp$datatype==3 & dat_tmp$nobs==100 , ],
                 aes(x=as.factor(alphaPrior), y=log(val+1), col=as.factor(Proc2))) +
    geom_boxplot(position = position_dodge(preserve = "total")) +
    facet_grid(ndx_label~Akap_label, label = "label_parsed", scales="free") +
    theme_bw(base_size=15) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("", "PC(0.1)","PC(0.9)")) +
    labs(x="", y=expression(log~~pwss(K^"+")), col="Procedure")

  pdf("~/Research/BYU/InformedFiniteMixtures/latex/figures/asymSimStudy_multi_Results_kplusvar.pdf",
      width=9, height=9)
  ggp7
  ggp8
  ggp9
  dev.off()

  ggp10 <- ggplot(data = dat_tmp[dat_tmp$metric=="ccp_diff" &
                                  dat_tmp$datatype==1 & dat_tmp$nobs==100 , ],
                 aes(x=as.factor(alphaPrior), y=val, col=as.factor(Proc2))) +
    geom_boxplot(position = position_dodge(preserve = "total")) +
    facet_grid(ndx_label~Akap_label, label = "label_parsed", scales="free") +
    theme_bw(base_size=15) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("", "PC(0.1)","PC(0.9)")) +
    labs(x="", y ="ccprob_error", col="Procedure")

  ggp11 <- ggplot(data = dat_tmp[dat_tmp$metric=="ccp_diff" &
                                  dat_tmp$datatype==2 & dat_tmp$nobs==100 , ],
                 aes(x=as.factor(alphaPrior), y=val, col=as.factor(Proc2))) +
    geom_boxplot(position = position_dodge(preserve = "total")) +
    facet_grid(ndx_label~Akap_label, label = "label_parsed", scales="free") +
    theme_bw(base_size=15) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("", "PC(0.1)","PC(0.9)")) +
    labs(x="",y ="ccprob_error", col="Procedure")


  ggp12 <- ggplot(data = dat_tmp[dat_tmp$metric=="ccp_diff" &
                                  dat_tmp$datatype==3 & dat_tmp$nobs==100 , ],
                 aes(x=as.factor(alphaPrior), y=val, col=as.factor(Proc2))) +
    geom_boxplot(position = position_dodge(preserve = "total")) +
    facet_grid(ndx_label~Akap_label, label = "label_parsed", scales="free") +
    theme_bw(base_size=15) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2), labels=c("", "PC(0.1)","PC(0.9)")) +
    labs(x="",y ="ccprob_error", col="Procedure")

  pdf("~/Research/BYU/InformedFiniteMixtures/latex/figures/asymSimStudy_multi_Results_ccp_diff.pdf",
      width=9, height=9)
  ggp10
  ggp11
  ggp12
  dev.off()


  # Create example of synthetic dataset employed
  vgrf <- read.xlsx("~/Research/BYU/BiomechanicCurves/Matt_ACLpatients/data2_25_2023/vgrf.xlsx", startRow=1, colNames=FALSE);
  vgrf <- sapply(vgrf, as.numeric)
  nsub <- (ncol(vgrf)/5)
  ntime <- nrow(vgrf)

  vgrf_mn <- NULL
  subject <- rep(1:nsub, each=5)
  for(i in 1:nsub){
    vgrf_mn <- cbind(vgrf_mn, apply(vgrf[,subject==i], 1, mean, na.rm=TRUE))
  }

  # Fit multivariate mixture to data
  mod1 <- Mclust(t(vgrf_mn), modelNames="EEI", G=4)
  clus_lab <- apply(mod1$z, 1, which.max)

  # compute cluster-means to be used to generate data
  vgrf_clus_mn <- NULL
  for(i in 1:4){
    vgrf_clus_mn <- cbind(vgrf_clus_mn, apply(vgrf_mn[,clus_lab==i], 1, mean, na.rm=TRUE))
  }


  ord <- 3
  beta_dim <- 30
  x <- seq(-3, 3, length=ntime)
  B <- miscPack:::bbase(x=x,ndx=beta_dim-3, bdeg=3)

  bbeta <- solve(t(B) %*% B) %*% t(B) %*% vgrf_clus_mn


  zi <- sample(1:4, size=nobs, replace=TRUE)
  yi <- NULL
  for(i in 1:nobs){
    yi <- cbind(yi, rnorm(ntime, B %*% (bbeta[,zi[i]]+rnorm(beta_dim,0,0.03)), 0.01))
  }

  pdf("~/Research/BYU/InformedFiniteMixtures/latex/figures/syntheticData2.pdf")
  matplot(x = seq(0.01,1,by=0.01), y=yi, col=zi, pch=1, type='b', xlab="% Stance Phase",
          ylab="")
  dev.off()
}


