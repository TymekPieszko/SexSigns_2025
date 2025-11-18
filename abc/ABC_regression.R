library(abc)
rm(list = ls())
mode <- "all" # all or subset
sim <- read.table("./MP_stats.tsv", header=T)
param <- sim[, "sex"]+3.162e-06
sim <- subset(sim, select = -c(sex, rec)) # Remove sex and rec
statsets <- list(
  "Hi"        = "Hi",
  "Fis"       = "Fis",
  "R2"        = "R2",
  "CL"        = "CL",
  "AR"        = "AR",
  "SX"        = "SX",
  "CL_bias"   = "CL_bias",
  "AR_bias"   = "AR_bias",
  "SX_bias"   = "SX_bias",
  "delta_m"   = "delta_m",
  "dt_ratio"  = "dt_ratio",
  "old"       = c("Hi","Fis","R2"),
  "tcomp"     = c("CL","AR","SX"),
  "tcomp_bias"= c("CL_bias","AR_bias","SX_bias"),
  "treebased" = c("delta_m","dt_ratio"),
  "new"       = c("CL_bias","AR_bias","SX_bias","delta_m","dt_ratio")
)
sim<-sim[names(statsets)[1:11]]

# Select simulated and pseudo-observed data
o <- 10 # Number of pseudo-observed
oidx <- as.vector(sapply(seq(1, nrow(sim), 100), function(x)
  x:(x + o -1))) # First o per param combination
obs <- sim[oidx, ]
sim <- sim[-oidx, ]
obspar <- param[oidx]
param <- param[-oidx]

res<-abc(target=obs[1,],param=log10(param),sumstat=sim,tol=0.01,method="rejection")
as.vector(summary(res))[c(4,2,6)]

coeffs <- c()  
for (statset in statsets) {
  predict.full<-matrix(NA,nrow=nrow(obs),ncol=3)
  for (i in (1:nrow(obs))) {
    res <- abc(target=obs[i,statset],param=log10(param),sumstat=sim[,statset],tol=0.01,method="rejection")
    predict.full[i,]<-as.vector(summary(res))[c(4,2,6)]
  }
  # plot(log10(unique(obspar)),by(predict.full[,1],FUN=mean,obspar))
  # plot(predict.full[,1] ~ log10(obspar))
  m <- lm(predict.full[,1] ~ log10(obspar))
  coeff <- summary(m)$r.squared
  coeffs <- c(coeffs, summary(m)$r.squared)
}
names(coeffs) <- names(statsets)
coeffs <- as.data.frame(t(coeffs))
write.table(coeffs, file=paste0("./out.txt"), sep="\t", row.names=F, col.names=T, quote=F)
coeffs
