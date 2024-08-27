# I am going to see how our method performs on the galaxy dataset
library(MASS)
library(miscPack)
library(MCMCpack) # needed for the rdirichlet function

library(fda)
#library(warptk)
library(salso)
library(parallel)
#install.packages("TruncatedDistributions", repos="http://R-Forge.R-project.org")
library(TruncatedDistributions) # pcprior.log_a requires truncated exponential


##################################################
#
# Exercise Science Movement Data: Curve Clustering
#
##################################################



# Try using Matt's data
library(openxlsx)
vgrf <- read.xlsx("~/Research/BYU/BiomechanicCurves/Matt_ACLpatients/data2_25_2023/vgrf.xlsx", startRow=1, colNames=FALSE);
vgrf <- sapply(vgrf, as.numeric)

nsub <- (ncol(vgrf)/5)
nobs <- nrow(vgrf)

vgrf_mn <- NULL
subject <- rep(1:nsub, each=5)
for(i in 1:nsub){
  vgrf_mn <- cbind(vgrf_mn, apply(vgrf[,subject==i], 1, mean, na.rm=TRUE))
}

# Figure S16
#pdf("~/Research/BYU/InformedFiniteMixtures/latex/figures/vgrf.pdf", height=6, width=8)
#  matplot((1:nobs/nobs), vgrf_mn, type='o',
#          pch=1, lty=1, cex=0.5,
#          ylab="Vertical Ground Reaction Force", xlab="% Stance Phase")
#dev.off()


y <- c(vgrf_mn)
t <- rep(1:nrow(vgrf_mn), nsub)
ids <- rep(1:nsub, each=nobs)

min_y <- min(y)
ran_y <- diff(range(y))

# Rescale so that each curve is on the same scale
rescale <- function(x) (x - min(x))/diff(range(x))
y_scaled <- (y-min_y)/ran_y

matplot(matrix(y_scaled, nrow=nobs, byrow=FALSE), type='o',
        pch=1, col="gray70", lty=1, cex=0.5,
        ylab="Frontal Knee Angle", xlab="Time")




niter <- 60000;nburn <- 10000;nthin <- 50;
nout <- (niter-nburn)/nthin
A <- 0.001
K <- 25

fit_model_given_U <- function(U, tp=1e-1, Akap=0.5, ndx=10, rwd=1){
  set.seed(101)
  pspline_mixture(y=y_scaled, t=t, ids=ids, K=K,
                  #fits1 <- pspline_mixture(y=y_scaled, t=t, ids=ids, K=K,
                  alpha_prior_type="centered",
                  alpha_prior_dist="pc",
                  basemodel=0,
                  U=U,
                  alpha1_val = 1,
                  alpha2_val = 1e-5,
                  update_alpha1=TRUE,
                  update_alpha2=FALSE,
                  tail.prob = tp,
                  ndx = ndx, q = 3, rwd=rwd,
                  ndens_y=1000,
                  A=A, Akap=Akap,
                  U_tau= 1/(1*0.31), a_tau=1e-2, U_omega=1/(1*0.31), a_omega=1e-2,
                  a_gam=1,
                  niter=niter, nburn=nburn, nthin=nthin)
}


plot_curves <- function(fits, U, data){
  pest1 <- salso(fits$z, loss="VI", maxNClusters=13, maxZealousAttempts = 20)
  fit_mat1 <- matrix(apply(fits$line_fit,2,mean)*ran_y + min_y, nrow=nobs, byrow=FALSE)
  nc1 <- length(unique(pest1))
  ord1 <- order(pest1)
  matplot((1:nobs)/nobs, data, pch=1, col=pest1, xlab="", ylab="")
  matplot((1:nobs)/nobs, fit_mat1, type="l", col=pest1, add=TRUE, lty=1,lwd=2)
  mtext(text=paste0("aFMM, U = ", U), side=3, cex=0.75, line=0)
  mtext(text="% Stance Phase", side=1, cex=0.75, line=2)
  mtext(text="Vertical Ground Reaction Force", side=2, cex=0.75, line=2)
}
plot_mean <- function(fits, U, data=vgrf_mn){
  pest1 <- salso(fits$z, loss="VI", maxNClusters=13, maxZealousAttempts = 20)
  fit_mat1 <- matrix(apply(fits$line_fit,2,mean)*ran_y + min_y, nrow=nobs, byrow=FALSE)
  nc1 <- length(unique(pest1))
  ord1 <- order(pest1)
  mncurve1 <- matrix(NA, nrow=nobs, ncol=nc1)
  plot(t/100, c(data), type="n",ylab="",
       xlab="")
  for(jj in 1:nc1){
    mncurve1[,jj] <- apply(fit_mat1[,pest1==jj,drop=FALSE],1,mean)

    lines((1:nobs)/100, mncurve1[,jj], col=jj,lwd=2.5)
  }
  mtext(text=paste0("aFMM, U = ", U), side=3, cex=0.75, line=0)
  mtext(text="% Stance Phase", side=1, cex=0.75, line=2)
  mtext(text="Mean Vertical Ground Reaction Force", side=2, cex=0.75, line=2)
  legend(x="topleft", legend=paste(1:nc1), col=1:nc1, lty=1,ncol=2)

}

date()
fits_vgrf1 <- mclapply(c(2:15), fit_model_given_U, tp=1e-1, Akap=0.50, ndx=10, rwd=2, mc.cores=9)
date()

sapply(fits_vgrf1, function(x) summary(x$kp))

pdf("~/Research/BYU/InformedFiniteMixtures/latex/figures/vgrf_fits_revU3.pdf", height=4, width=8)
  par(mfrow=c(1,2), mar=c( 3.0, 3.5, 1.5, 1.5))# bottom, left, top, right
  plot_curves(fits_vgrf1[[2]], U=3, data=vgrf_mn)
  plot_mean(fits_vgrf1[[2]], U=3)
dev.off()

# Figure 8
pdf("~/Research/BYU/InformedFiniteMixtures/latex/figures/vgrf_fits_revU15.pdf", height=4, width=8)
  par(mfrow=c(1,2), mar=c( 3.0, 3.5, 1.5, 1.5))# bottom, left, top, right
  plot_curves(fits_vgrf1[[14]], U=15, data=vgrf_mn)
  plot_mean(fits_vgrf1[[14]], U=15)
dev.off()

plot(summary(salso(fits_vgrf1[[2]]$z)))
plot(summary(salso(fits_vgrf1[[14]]$z)))

ARI(salso(fits_vgrf1[[2]]$z), salso(fits_vgrf1[[14]]$z))
which(salso(fits_vgrf1[[2]]$z) != salso(fits_vgrf1[[14]]$z))

# Figure 9
library(fields)
pdf("~/Research/BYU/InformedFiniteMixtures/latex/figures/coclustering_rev.pdf", height=3, width=9)
#  par(mfrow=c(1,2), mar=c(2.0, 1.75, 1.25, 1.5))# bottom, left, top, right
  par(mfrow=c(1,3), mai=c(0.35, 0.35, 0.1, 0.35))# bottom, left, top, right

  ord1 <- order(salso(fits_vgrf1[[2]]$z))
  image.plot(psm(fits_vgrf1[[2]]$z)[ord1, ord1],
             axes=FALSE, xlab="", ylab="", col = tim.colors(),
             smallplot= c(.85,0.87,0.075,0.95),
             axis.args=list(cex=0.5))
  axis(1, at = c(0:(nsub-1))/(nsub-1), labels=c(1:nsub)[ord1], tick=TRUE, cex.axis=0.25, las=3)
  axis(2, at = c(0:(nsub-1))/(nsub-1), labels=c(1:nsub)[ord1], tick=TRUE, cex.axis=0.25, las=2)


  ord2 <- order(salso(fits_vgrf1[[14]]$z))
  image.plot(psm(fits_vgrf1[[14]]$z)[ord2, ord2],
             axes=FALSE, xlab="", ylab="", col = tim.colors(),
             smallplot= c(.85,0.87,0.075,0.95))
  axis(1, at = c(0:(nsub-1))/(nsub-1), labels=c(1:nsub)[ord2], tick=TRUE, cex.axis=0.25, las=3)
  axis(2, at = c(0:(nsub-1))/(nsub-1), labels=c(1:nsub)[ord2], tick=TRUE, cex.axis=0.25, las=2)


  coclust_prob <- lapply(fits, function(x) psm(x$z))
  diff_coclust_prob <- coclust_prob[[2]] - coclust_prob[[14]]
  image.plot(diff_coclust_prob[ord2, ord2],
             axes=FALSE, xlab="", ylab="", col = tim.colors())
  axis(1, at = c(0:(nsub-1))/(nsub-1), labels=c(1:nsub)[ord2], tick=TRUE, cex.axis=0.25, las=3)
  axis(2, at = c(0:(nsub-1))/(nsub-1), labels=c(1:nsub)[ord2], tick=TRUE, cex.axis=0.25, las=2)


dev.off()


