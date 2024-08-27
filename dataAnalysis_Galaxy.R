# I am going to see how our method performs on the galaxy dataset
library(MASS)
library(miscPack)
library(TruncatedDistributions) # pcprior.log_a requires truncated exponential
library(MCMCpack) # needed for the rdirichlet function
library(mixAK) # needed for NMixMCMC function (RJMCMC finite mixture)
library(AntMAN) # needed for the AM_mcmc_fit function
library(salso) # needed for the psm and salso functions
library(fields) # needed for image.plot function
##############
#
# Galaxy data
#
##############
data(galaxies)
yi <- galaxies/1000
summary(yi)
nobs <- length(yi)

K <- 25

niter <- 100000; nburn <- 50000; nthin <- 50; nout <- (niter - nburn)/nthin

# Recall that the mu_sigma_prior has the following arguments
# mu_sigma_prior = 1 - muk ~ N(m0, s20), sigmak ~ UN(0,A)
# mu_sigma_prior = 2 - muk ~ N(m0, s20), sigma2k ~ IG(a0, b0)
# mu_sigma_prior = 3 - (muk, sigma2k) ~ NIG(m0, k0, nu0, s20)

fits <- list()
kp_prior <- list()
psm_list <- list()
U <- c(2:10)
h <- 1
for(i in 1:length(U)){
  for(a2 in 1:1){
    if(a2==1) alpha2_val = 1e-5
    for(tp in 1:2){
      if(tp == 1) tail_prob=0.1
      if(tp == 2) tail_prob=0.5
      set.seed(h)
      fits[[h]] <- informed_mixture(y=yi, K=K,
                                          alpha_prior_type="centered",
                                          alpha_prior_dist="pc",
                                          mu_sigma_prior = 2, # muk ~ N(m0, s20), sigma2k ~ IG(a0, b0)
                                          basemodel=0,
                                          U=U[i],
                                          alpha2_val = alpha2_val,
                                          update_alpha1=TRUE,
                                          update_alpha2=FALSE,
                                          tail.prob = tail_prob, ndens_y=1000,
                                          hierarchy="NO",
                                          m0=mean(yi), s20=10^2, a0=3, b0 = 2, # prior mean of sig2 is b0/(a0-1)
#                                          m0=mean(yi), s20=10^2, a0=3, b0 = 10, # prior mean of sig2 is b0/(a0-1)
                                          niter=niter, nburn=nburn, nthin=nthin)

      print(summary(fits[[h]]$alpha1))
      kp_prior[[h]]<- miscPack:::impliedpcprior.asym.a1.Kplus(lam=fits[[h]]$lambda_val,
                                                              U = U[i],
                                                              K = K,
                                                              alpha1.base=U[i],
                                                              alpha2.base=1e-5,
                                                              alpha2.fixed=1e-5,
                                                              agrid=fits[[h]]$alpha_grid,
                                                              n.obs = length(yi),
                                                              n.samples = 1e4,
                                                              TR=TRUE)$Kplus

      psm_list[[h]] <- psm(fits[[h]]$z)


      h <- h+1

    }
  }
}




alpha1_post_mn <- sapply(fits, function(x) summary(x$alpha1))


pest <- sapply(fits, function(x) salso(x$z, loss="VI"))


ent <- function(x){
  pp<-x[upper.tri(x,diag=FALSE)];
  et1 <- pp*log(pp);
  et2 <- (1-pp)*log(1-pp);
  et1[is.na(et1)] <- 0;
  et2[is.na(et2)] <- 0;
  -sum(et1 + et2)}

kp_entropy <- sapply(psm_list, ent)

purity <- function(fits){
  x <- fits$z
  pest <- salso(x)
  ccp <- psm(x)
  sapply(1:nrow(ccp), function(y) c(mean(ccp[y, pest==pest[y]])*mean(1-ccp[y,pest!=pest[y]])))
}
c_pure <- sapply(fits, purity)




kp_dist <- sapply(fits, function(x) table(x$kp)/nout)

mse_est <- sapply(fits, function(x) sapply(1:nout, function(y) mean((x$mu[y, x$z[y,]] - yi)^2)))

sd_psm <- sapply(psm_list, function(x) apply(x,1,sd))


pdf("~/Research/BYU/InformedFiniteMixtures/latex/figures/postKp_galaxy_2.pdf", height=9, width=6)

  par(mfrow=c(4,2), mar=c( 2.75, 3.5, 1.75, 1.5))# bottome, left, top, right
  UU <- rep(U, each=2)
  tp <- rep(c(0.1, 0.5), times=length(U))
  loc <- rep(c("topright","topleft"), each=length(U))
  for(j in 1:length(UU)){
    if(UU[j] %in% c(3,5,7,10)){
      plot(as.numeric(names( kp_dist[[j]]))+0.1,  kp_dist[[j]], type='h',
           xlim=c(1,max(UU)+2), ylim=c(0,1),
           lwd=2, ylab="", xlab="")
      mtext(text=bquote(U == .(UU[j])*","~ tp == .(tp[j])), side=3, cex=0.75)
      mtext(text=expression(K^'+'), side=1, cex=0.75, line=2)
      mtext(text=expression(Pr(K^'+')), side=2, cex=0.75, line=2)
      axis(1,at=1:12)
      axis(2,at=seq(0,1,by=0.25))
      kpp <- table(kp_prior[[j]])
      lines(as.numeric(names(kpp))-0.1, kpp/sum(kpp), type='h', col='red', lwd=2) # U = 5
      legend(x=loc[j], legend=c("Prior","Posterior"), col=c("red","black"), lty=1, lwd=3)
    }
  }
dev.off()

pdf("~/Research/BYU/InformedFiniteMixtures/latex/figures/postCCP_galaxy_3.pdf", height=9, width=9)
  par(mfrow=c(3,3), mar=c( 1.5, 0.5, 1.25, 4.2))# bottome, left, top, right
  UU <- rep(U, each=2)
  tp <- rep(c(0.1, 0.5), times=length(U))
  mval <- apply(mse_est,2,mean)
  lval <- apply(sd_psm,2,mean)
  for(j in 1:length(UU)){
    if(tp[j] == 0.1){
      print(j)
      image.plot(psm_list[[j]], xaxt="n", yaxt="n")
       mtext(text=bquote(U == .(UU[j])*","~
                         tp == .(tp[j])*"," ~
                         mse == .(round(mval[j]/(K-UU[j]),2))*"," ~
                         sd_ccp == .(round(lval[j],2))),
             side=3, cex=0.75)
    }
  }

dev.off()


pdf("~/Research/BYU/InformedFiniteMixtures/latex/figures/density_est_galaxy.pdf", height=8, width=8)
  par(mfrow=c(1,1))

  clrs <- tim.colors(length(UU))
  hist(yi, breaks=30, freq=FALSE, main="", xlab="galaxy velocity", ylab="density")
  line_func <- function(x) {
    if(tp[x] == 0.1){
      lines(fits[[x]]$ygrid,apply(fits[[x]]$density,2,mean), col=clrs[x], lwd=1.5);
    }
  }
  sapply(1:length(UU), line_func)
  legend(x='topright', legend=paste0("U=",UU[tp==0.1]), col=clrs[tp==0.1], lty=1,
         ncol=3, lwd=1.75)
dev.off()

