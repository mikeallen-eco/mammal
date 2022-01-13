# =========================================================================
#
# 11. Hierarchical models for communities
#
# =========================================================================





# 11.1 Introduction
# ------------------------------------------------------------------------



# 11.2 Simulation of a metacommunity
# ------------------------------------------------------------------------


simComm(type="det/nondet", nsite=30, nrep=3, nspec=100,
        mean.psi=0.25, sig.lpsi=1, mu.beta.lpsi=0, sig.beta.lpsi=0,
        mean.lambda=2, sig.loglam=1, mu.beta.loglam=1, sig.beta.loglam=1,
        mean.p=0.25, sig.lp=1, mu.beta.lp=0, sig.beta.lp=0, show.plot = TRUE)


# Execute function with default arguments
set.seed(1234)
data <- simComm(type="det/nondet", nsite=30, nrep=3, nspec=100,
                mean.psi=0.25, sig.lpsi=1, mu.beta.lpsi=0, sig.beta.lpsi=0,
                mean.lambda=2, sig.loglam=1, mu.beta.loglam=1, sig.beta.loglam=1,
                mean.p=0.25, sig.lp=1, mu.beta.lp=0, sig.beta.lp=0, show.plot = TRUE)
# data <- simComm() # same

str(data)

# Some possibly interesting settings of the function
data <- simComm(nsite = 267, nspec = 190, mean.psi = 0.25, sig.lpsi = 2, mean.p = 0.12, sig.lp = 2) # similar to Swiss MHB
data <- simComm(mean.psi = 1)         # all species occur at every site
data <- simComm(mean.p = 1)           # no measurement error (perfect detection)

# Effect of spatial sample size (nsite) on species richness in sample (Ntotal.fs)
data <- simComm(nsite=50, nspec = 200) # 1-3 are usually missed in sample
data <- simComm(nsite=30, nspec = 200) # 4-6 usually missed
data <- simComm(nsite=10, nspec = 200) # around 30 typically missed

# Check for frequentist characteristics of such statistics
temp <- rep(NA, 100)
for(i in 1:100){
  cat("\nSimrep", i)
  temp[i] <- simComm(nsite=10, nspec = 200, show.plot = F)$Ntotal.fs
}
hist(200-temp, breaks = 30, main = "Number of species in the metacommunity \nthat do not occur in the 10 sampled sites", col = "gold")

# Simulation 1: effects of psi and sd(logit(psi)) on number of species actually occurring in the 50 sampled sites
simrep <- 50                   # Run 50 simulation reps
mpsi <- seq(0.01, 0.25,,10)
slpsi <- seq(0.1, 5,,10)
results1 <- array(NA, dim = c(10, 10, simrep))
for(i in 1:10){      # Loop over levels of factor mean.psi (mpsi)
  for(j in 1:10){    # Loop over levels of factor sig.lpsi (slpsi)
    for(k in 1:simrep){
      cat("\nDim 1:",i, ", Dim 2:", j, ", Simrep", k)
      tmp <-  simComm(nsite=50, nspec = 200, show.plot = F, mean.psi = mpsi[i],
                      sig.lpsi = slpsi[j])
      results1[i,j,k] <- tmp$Ntotal.fs
    }
  }
}


# Simulation 2: effects of p and sd(logit(p)) on the proportion of the species occurring in the 50 sampled sites that are detected at least once
simrep <- 50         # Run 50 simulation reps again
mp <- seq(0.01, 0.25,,10)
slp <- seq(0.1, 5,,10)
results2 <- array(NA, dim = c(10, 10, simrep, 2))
for(i in 1:10){      # Loop over levels of factor mean.p (mp)
  for(j in 1:10){    # Loop over levels of factor sig.lp (slp)
    for(k in 1:simrep){
      cat("\nDim 1:",i, ", Dim 2:", j, ", Simrep", k)
      tmp <-  simComm(nsite=50, nspec = 200, show.plot = F, mean.p = mp[i],
                      sig.lp = slp[j])
      results2[i,j,k,] <- c(tmp$Ntotal.fs, tmp$Ntotal.obs)
    }
  }
}


# Plot these two prediction matrices
par(mfrow = c(1, 2), mar = c(5,5,2,2), cex.lab = 1.5, cex.axis = 1.5)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))

# Plot proportion of species occurring in sampled sites (Fig. 11-4 left)
z1 <- apply(results1/200, c(1,2), mean) # Prop species occurring
image(x=mpsi, y=slpsi, z=z1, col = mapPalette(100), axes = T, xlab = "Occupancy probability (psi)", ylab = "Among-species variability in psi")
contour(x=mpsi, y=slpsi, z=z1, add = T, col = "blue", labcex = 1.5, lwd = 1.5)

# Plot proportion of species detected in sampled sites (Fig. 11-4 right)
z2 <- apply(results2[,,,2] / results2[,,,1], c(1,2), mean)
image(x=mp, y=slp, z=z2, col = mapPalette(100), axes = T, xlab = "Detection probability (p)", ylab = "Among-species variability in p")
contour(x=mp, y=slp, z=z2, add = T, col = "blue", labcex = 1.5, lwd = 1.5)



# 11.3 Metacommunity data from the Swiss breeding bird survey MHB
# ------------------------------------------------------------------------

## Code modified to use the built-in data set MHB2014 instead of the file "MHB_2014.csv"
data(MHB2014)
?MHB2014
str(MHB2014)
# NB some of the data preprocessing on p.644 has already been done.

# Check the detection data in 3D array MHB2014$count: site x rep x species
( nsite <- nrow(MHB2014$sites) )    # number of sites in Swiss MHB
nrep <- 3                           # maximum number of replicate surveys per season
( nspec <- nrow(MHB2014$species) )  # 158 species occur in the 2014 data
dim(MHB2014$count) == c(nsite, nrep, nspec) # check

# Create the detection/nondetection (1/0) array
y <- MHB2014$count ; y[y > 1] <- 1  ## 'Y' replaced with 'y'
str(y)

# Check data for one species, here chaffinch, and pull them out from 3D array
(tmp <- y[, , "Common Chaffinch"])

# Frequency distribution of number of surveys actually carried out per site in 2014
# NB MHB2014$sites$nsurvey gives the number of surveys *planned*.
table(nsurveys <- apply(!is.na(y[,,1]), 1, sum))

# Which site has all NA data in 2014 ?
(NAsites <- which(nsurveys == 0) )

# Observed number of occupied sites
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
# For the 'all NA' site, max returns -Inf with a warning
tmp[tmp == -Inf] <- NA         # Change -Inf to NA
sort(obs.occ <- apply(tmp, 2, sum, na.rm = TRUE))

# Plot species 'occurrence frequency' distribution (not shown)
plot(sort(obs.occ), xlab = "Species number", ylab = "Number of quads with detections")

# Drop data from species that were not observed in 2014
toss.out <- which(obs.occ == 0)
y <- y[,,-toss.out]
obs.occ <- obs.occ[-toss.out]

# Redefine nspec as the number of species observed in 2014: 145
( nspec <- dim(y)[3] )

str(y)

# Get observed number of species per site
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
tmp[tmp == "-Inf"] <- NA
sort(C <- apply(tmp, 1, sum))     # Compute and print sorted species counts

plot(table(C), xlim = c(0, 60), xlab = "Observed number of species", ylab = "Number of quadrats", frame = FALSE)
abline(v = mean(C, na.rm = TRUE), col = "blue", lwd = 3)



# 11.4 Overview of some models for metacommunities
# ------------------------------------------------------------------------



# 11.5 Community models that ignore species identity
# ------------------------------------------------------------------------



# 11.5.1 Simple Poisson regression for the observed community size
# ------------------------------------------------------------------------
# Get covariates (from AHMbook::MHB2014) and standardise them

# Quadrat elevation and forest cover
orig.ele <- MHB2014$sites$elev
(mean.ele <- mean(orig.ele, na.rm = TRUE))
(sd.ele <- sd(orig.ele, na.rm = TRUE))
ele <- (orig.ele - mean.ele) / sd.ele
orig.forest <- MHB2014$sites$forest
(mean.forest <- mean(orig.forest, na.rm = TRUE))
(sd.forest <- sd(orig.forest, na.rm = TRUE))
forest <- (orig.forest - mean.forest) / sd.forest

# Average date and duration of survey
orig.mdate <- apply(MHB2014$date, 1, mean, na.rm = TRUE)
(mean.mdate <- mean(orig.mdate[-NAsites]))   # drop unsurved site
(sd.mdate <- sd(orig.mdate[-NAsites]))
mdate <- (orig.mdate - mean.mdate) / sd.mdate
mdate[NAsites] <- 0                 # impute mean for missing

orig.mdur <- apply(MHB2014$dur, 1, mean, na.rm = TRUE)
(mean.mdur <- mean(orig.mdur[-NAsites]))
(sd.mdur <- sd(orig.mdur[-NAsites]))
mdur <- (orig.mdur - mean.mdur) / sd.mdur
mdur[NAsites] <- 0                  # impute mean for missing

# Bundle data and summarize input data for BUGS
str( win.data <- list(C = C, nsite = length(C), ele = ele, forest = forest,
                      mdate = mdate, mdur = mdur) )

# Specify model in BUGS language
sink("model1.txt")
cat("
model {

# Priors
gamma0 ~ dnorm(0, 0.001)           # Regression intercept
for(v in 1:6){                     # Loop over regression coef's
   gamma[v] ~ dnorm(0, 0.001)
}

# Likelihood for Poisson GLM
for (i in 1:nsite){
   C[i] ~ dpois(lambda[i])
   log(lambda[i]) <- gamma0 + gamma[1] * ele[i] + gamma[2] * pow(ele[i],2) + gamma[3] * forest[i] + gamma[4] * mdate[i] + gamma[5] * pow(mdate[i],2) + gamma[6] * mdur[i]
}
}
",fill = TRUE)
sink()

# Initial values
inits <- function() list(gamma0 = rnorm(1), gamma = rnorm(6))

# Parameters monitored
params <- c("gamma0", "gamma")

# MCMC settings
ni <- 6000   ;   nt <- 4   ;   nb <- 2000   ;   nc <- 3

# Call WinBUGS from R (ART <1 min)
out1 <- bugs(win.data, inits, params, "model1.txt", n.chains = nc,
             n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Call JAGS from R (ART <1 min), check convergence and summarize posteriors
library(jagsUI)
out1J <- jags(win.data, inits, params, "model1.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
traceplot(out1J)    ;    print(out1J, dig = 3)


# Plot posterior distributions for potentially 'ecological' parameters
par(mfrow = c(1,3))
hist(out1J$sims.list$gamma[,1], breaks = 50, col = "grey", main = "", xlab = "Slope of elevation (linear)")
hist(out1J$sims.list$gamma[,2], breaks = 50, col = "grey", main = "", xlab = "Slope of elevation (squared)")
hist(out1J$sims.list$gamma[,3], breaks = 50, col = "grey", main = "", xlab = "Slope of forest")


# Get covariate values for prediction
orig.pred.ele <- seq(250, 2750,, 500)   # 500 vals spread between 250 and 2750
p.ele <- (orig.pred.ele - mean.ele) / sd.ele
orig.pred.forest <- seq(1, 100,, 500)
p.forest <- (orig.pred.forest - mean.forest) / sd.forest

# Compute predictions
nsamp <- out1J$mcmc.info$n.samples
pred.ele <- pred.forest <- array(NA, dim = c(500, nsamp))
for(i in 1:nsamp){
  pred.ele[,i] <- exp(out1J$sims.list$gamma0[i] + out1J$sims.list$gamma[i,1] * p.ele + out1J$sims.list$gamma[i,2]* p.ele^2)
  pred.forest[,i] <- exp(out1J$sims.list$gamma0[i] + out1J$sims.list$gamma[i,3] * p.forest)
}

# Plot posterior mean and a random sample of 100 from posterior of regression
selection <- sample(1:nsamp, 100)
par(mfrow = c(1,3))
matplot(orig.pred.ele, pred.ele[,selection], ylab = "Predicted species count", xlab = "Elevation (m a.s.l.)", type = "l", lty = 1, lwd = 1, col = "grey", ylim = c(0, 50), frame = F)
lines(orig.pred.ele, apply(pred.ele, 1, mean), lwd = 3, col = "blue")
matplot(orig.pred.forest, pred.forest[,selection], ylab = "Predicted species count", xlab = "Forest cover (%)", type = "l", lty = 1, lwd = 1, col = "grey", ylim = c(0, 50), frame = F)
lines(orig.pred.forest, apply(pred.forest, 1, mean), lwd = 3, col = "blue")


# Get observed species richness per site and rep and plot
CC <- apply(y, c(1,2), sum, na.rm = TRUE)
CC[CC == 0] <- NA            # 0 means not surveyed
matplot(t(CC), type = 'l', lty = 1, lwd = 2, xlab = "First to third survey", ylab = "Number of species detected", frame = F)  # Fig. 11?6 right

# Get survey date and survey duration and standardise both
# Survey date (this is Julian date, with day 1 being April 1)
orig.DAT <- MHB2014$date
(mean.date <- mean(orig.DAT, na.rm = TRUE))
(sd.date <- sd(c(orig.DAT), na.rm = TRUE))
DAT <- (orig.DAT - mean.date) / sd.date      # scale
DAT[is.na(DAT)] <- 0                         # impute missings
# Survey duration (in minutes)
orig.DUR <- MHB2014$dur
(mean.dur <- mean(orig.DUR, na.rm = TRUE))
(sd.dur <- sd(c(orig.DUR), na.rm = TRUE))
DUR <- (orig.DUR - mean.dur) / sd.dur        # scale
DUR[is.na(DUR)] <- 0                         # mean impute missings

# Bundle data and summarize
str( win.data <- list(CC = CC, M = nrow(CC), J = ncol(CC), ele = ele, forest = forest, DAT = DAT, DUR = DUR) )


# Specify model in BUGS language
sink("model2.txt")
cat("
model {

# Priors
gamma0 ~ dnorm(0, 0.001)
for(v in 1:6){
  gamma[v] ~ dnorm(0, 0.001)
}

# Likelihood for Poisson GLM
for (i in 1:M){             # Loop over sites
  for(j in 1:J){            # Loop over occasions
    CC[i,j] ~ dpois(lambda[i,j])
    log(lambda[i,j]) <- gamma0 + gamma[1] * ele[i] + gamma[2] * pow(ele[i],2) +
      gamma[3] * forest[i] + gamma[4] * DAT[i,j] + gamma[5] * pow(DAT[i,j],2) +
      gamma[6] * DUR[i,j]
  }
}
}
",fill = TRUE)
sink()

# Initial values
inits <- function() list(gamma0 = rnorm(1), gamma = rnorm(6))

# Parameters monitored
params <- c("gamma0", "gamma")

# MCMC settings
ni <- 6000   ;   nt <- 4   ;   nb <- 2000   ;   nc <- 3

# Call WinBUGS from R (ART 1.7 min)
out2 <- bugs(win.data, inits, params, "model2.txt", n.chains = nc,
             n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Call JAGS from R (ART 0.6 min)
out2J <- jags(win.data, inits, params, "model2.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
traceplot(out2J)   ;   print(out2J, dig = 2)


# 11.5.2 Poisson random effects model for the observed community size
# ------------------------------------------------------------------------
# Bundle and summarize data set
str( win.data <- list(CC = CC, M = nrow(CC), J = ncol(CC), ele = ele, forest = forest, DAT = DAT, DUR = DUR) )

# Specify model in BUGS language
sink("model3.txt")
cat("
model {

# Priors
mugamma0 ~ dnorm(0, 0.001)     # Hyperparameters
taugamma0 <- pow(sd.gamma0,-2)
sd.gamma0 ~ dunif(0, 10)
for(v in 1:6){                 # Parameters
   gamma[v] ~ dnorm(0, 0.001)
}

# Likelihood for Poisson GLMM
for (i in 1:M){                # Loop over sites
  gamma0[i] ~ dnorm(mugamma0, taugamma0)     # site intercepts random now
  for(j in 1:J){               # Loop over repeated measurements
    CC[i,j] ~ dpois(lambda[i,j])
    log(lambda[i,j]) <- gamma0[i] + gamma[1]*ele[i] + gamma[2] * pow(ele[i],2) +
      gamma[3] * forest[i] + gamma[4] * DAT[i,j] + gamma[5] * pow(DAT[i,j],2) +
      gamma[6] * DUR[i,j]
  }
}
}
",fill = TRUE)
sink()

# Initial values
inits <- function() list(gamma0 = rnorm(nrow(CC)), gamma = rnorm(6))

# Parameters monitored
params <- c("mugamma0", "sd.gamma0", "gamma0", "gamma")

# MCMC settings
ni <- 6000   ;   nt <- 4   ;   nb <- 2000   ;   nc <- 3

# Call WinBUGS from R (ART 2.9 min)
out3 <- bugs(win.data, inits, params, "model3.txt", n.chains = nc,
             n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Call JAGS from R (ART 0.7 min)
out3J <- jags(win.data, inits, params, "model3.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
traceplot(out3J, c('mugamma0', 'sd.gamma0', 'gamma'))   ;    print(out3J, dig = 2)


# 11.5.3 N-mixture model for the observed community size
# -----------------------------------------------------
# Bundle and summarize data set
str( win.data <- list(CC = CC, M = nrow(CC), J = ncol(CC), ele = ele, forest = forest, DAT = DAT, DUR = DUR) )

# Specify model in BUGS language
sink("model4.txt")
cat("
model {

# Priors
alpha0 ~ dnorm(0, 0.01)      # Base-line community detection probability
beta0 ~ dnorm(0, 0.01)       # Base-line community size (number of species)
for(v in 1:3){
   alpha[v] ~ dnorm(0, 0.01) # Covariate effects on detection
   beta[v] ~ dnorm(0, 0.01)  # Covariate effects on community size
}

# Likelihood
# Ecological model for true community size
for (i in 1:M){              # Loop over sites
   N[i] ~ dpois(lambda[i])   # Community size
   lambda[i] <- exp(beta0 + beta[1] * ele[i] + beta[2] * pow(ele[i],2) +
      beta[3] * forest[i])

   # Observation model for repeated measurements
   for (j in 1:J){          # Loop over occasions
      CC[i,j] ~ dbin(p[i,j], N[i])
      p[i,j] <- 1 / (1 + exp(-lp[i,j]))
      lp[i,j] <- alpha0 + alpha[1] * DAT[i,j] + alpha[2] * pow(DAT[i,j],2) +
         alpha[3] * DUR[i,j]
   # logit(p) = ... causes undefined real result in WinBUGS (but not JAGS)
   }
}
}
",fill = TRUE)
sink()

# Define function to generate random initial values
Nst <- apply(CC, 1, max, na.rm = TRUE) + 1
Nst[Nst == -Inf] <- max(Nst, na.rm = T)  # Some nonzero val. for unsurv. sites
inits <- function() list(N = Nst, alpha0 = rnorm(1), alpha = rnorm(3), beta0 = rnorm(1), beta = rnorm(3))

# Parameters monitored
params <- c("alpha0", "alpha", "beta0", "beta","N")

# MCMC settings
ni <- 6000   ;   nt <- 4   ;   nb <- 2000   ;   nc <- 3

# Run JAGS from R (ART 1.5 min) in parallel, look at traceplots
#   and summarize posteriors
out4 <- jags(win.data, inits, params, "model4.txt",
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
traceplot(out4, c('alpha0', 'alpha', 'beta0', 'beta'))  ;  print(out4, 3)

print(tmp <- cbind(out3$summary[c(1, 270:272),c(1:3,7)], out4$summary[5:8, c(1:3, 7)]), 3)

plot(orig.ele, out4$summary[9:275, 1], pch = 16, xlab = "Elevation (m)", ylab = "Species richness", ylim = c(0, 60), frame = F)
segments(orig.ele, out4$summary[9:275, 3], orig.ele, out4$summary[9:275, 7])
points(orig.ele+20, C)      # elevation jittered




# 11.6 Community models that retain species identity
# ------------------------------------------------------------------------


# 11.6.1 Simplest community occupancy model: n-fold single species
#        occupancy model with species treated as fixed effects
# ------------------------------------------------------------------------
# Collapse 3D detection/nondetection data to 2D detection frequencies
ysum <- apply(y, c(1,3), sum, na.rm = T) # Collapse to detection frequency
ysum[NAsites,] <- NA                     # Have to NA out sites with NA data

# Bundle and summarize data set
str( win.data <- list(ysum = ysum, M = nrow(ysum), J = MHB2014$sites$nsurvey, nspec = dim(ysum)[2]) )

# Specify model in BUGS language
sink("model5.txt")
cat("
model {

# Priors
for(k in 1:nspec){          # Loop over species
   psi[k] ~ dunif(0, 1)
   p[k] ~ dunif(0, 1)
}

# Ecological model for latent occurrence z (process model)
for(k in 1:nspec){          # Loop over species
   for (i in 1:M) {         # Loop over sites
      z[i,k] ~ dbern(psi[k])
   }
}

# Observation model for observed data y
for(k in 1:nspec){          # Loop over species
   for (i in 1:M) {
      mup[i,k] <- z[i,k] * p[k]
      ysum[i,k] ~ dbin(mup[i,k], J[i])
   }
}

# Derived quantities
for(k in 1:nspec){          # Loop over species
   Nocc.fs[k] <- sum(z[,k]) # Add up number of occupied sites among the 267
}
for (i in 1:M) {            # Loop over sites
   Nsite[i] <- sum(z[i,])   # Add up number of occurring species at each site
}
}
",fill = TRUE)
sink()

# Initial values
zst <- apply(y, c(1,3), max) # Observed occurrence as inits for z
zst[is.na(zst)] <- 1
inits <- function() list(z = zst, psi = rep(0.4, nspec), p = rep(0.4, nspec))

# Parameters monitored
params <- c("psi", "p", "Nsite", "Nocc.fs")

# MCMC settings
ni <- 2500   ;   nt <- 2   ;   nb <- 500   ;   nc <- 3

# Call JAGS from R (ART 2.1 min)
out5 <- jags(win.data, inits, params, "model5.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
par(mfrow = c(4,4))   ;   traceplot(out5)   ;   print(out5, dig = 3)

# Compare observed and estimated site species richness
par(cex = 1.3)
twice <- which(MHB2014$sites$nsurvey == 2)
plot(C[twice], out5$summary[291:557,1][twice], xlab = "Observed number of species", ylab = "Estimated number of species", frame = F, xlim = c(0, 60), ylim = c(0, 70), col = "red", pch = 16)
segments(C[twice], out5$summary[291:557,3][twice], C[twice], out5$summary[291:557,7][twice], col = "red")
points(C[-twice], out5$summary[291:557,1][-twice], col = "blue", pch = 16)
segments(C[-twice], out5$summary[291:557,3][-twice], C[-twice], out5$summary[291:557,7][-twice], col = "blue")
abline(0,1)

# Observed and estimated number of occupied sites for each species
# in a table
cbind(obs.occu = obs.occ, out5$summary[558:702, c(1,3,7)])

# and in a plot
plot(obs.occ, out5$summary[558:702, 1], xlab = "Observed number of occupied sites", ylab = "Estimated version of quantity", ylim = c(0, 267), frame = F, pch = 16)
abline(0,1)
segments(obs.occ, out5$summary[558:702,3], obs.occ, out5$summary[558:702,7], col = "grey", lwd = 2)


# Estimated occupancy and detection probability for each species
plot(out5$summary[1:145,1], out5$summary[146:290,1], xlab = "Occupancy estimate", ylab = "Detection estimate", xlim = c(0,1), ylim = c(0,1), frame = F, pch = 16)
segments(out5$summary[1:145,3], out5$summary[146:290,1], out5$summary[1:145,7], out5$summary[146:290,1], col = "grey", lwd = 2)
segments(out5$summary[1:145,1], out5$summary[146:290,3], out5$summary[1:145,1], out5$summary[146:290,7], col = "grey", lwd = 2)


# 11.6.2 Community occupancy model with bivariate species-specific random effects
# -----------------------------------------------------
# Bundle and summarize data set
str( win.data <- list(ysum = ysum, M = nrow(ysum), J = MHB2014$sites$nsurvey, nspec = dim(ysum)[2], R = matrix(c(5,0,0,1), ncol = 2), df = 3) )


# Specify model in BUGS language
sink("model6.txt")
cat("
model {

# Priors
for(k in 1:nspec){  # Group lpsi and lp together in array eta
   lpsi[k] <- eta[k,1]
   lp[k] <- eta[k,2]
   eta[k, 1:2] ~ dmnorm(mu.eta[], Omega[,])
}
# Hyperpriors
# Priors for mu.lpsi=mu.eta[1] and mu.lp=mu.eta[2]
# probs = community means of occupancy and detection probability
for(v in 1:2){
   mu.eta[v] <- log(probs[v] / (1-probs[v]))
   probs[v] ~ dunif(0,1)
}
# Prior for variance-covariance matrix
Omega[1:2, 1:2] ~ dwish(R[,], df)
Sigma[1:2, 1:2] <- inverse(Omega[,])

# Ecological model for latent occurrence z (process model)
for(k in 1:nspec){
   logit(psi[k]) <- lpsi[k]   # Must take outside of i loop
   for (i in 1:M) {
      z[i,k] ~ dbern(psi[k])
   }
}

# Observation model for observed data y
for(k in 1:nspec){
   logit(p[k]) <- lp[k]       # Needs to be outside of i loop
   for (i in 1:M) {
      mu.p[i,k] <- z[i,k] * p[k]
      ysum[i,k] ~ dbin(mu.p[i,k], J[i])
   }
}

# Derived quantities
rho <- Sigma[1,2] / sqrt(Sigma[1,1] * Sigma[2,2])  # Correlation coefficient
for(k in 1:nspec){
   Nocc.fs[k] <- sum(z[,k])   # Number of occupied sites among the 267
}
for (i in 1:M) {
   Nsite[i] <- sum(z[i,])     # Number of occurring species
}
}
",fill = TRUE)
sink()


# Initial values
zst <- apply(y, c(1,3), max) # Observed occurrence as starting values for z
zst[is.na(zst)] <- 1
inits <- function() list(z = zst, Omega = matrix(c(1,0,0,1), ncol = 2), eta = matrix(0, nrow = nspec, ncol = 2))

# Parameters monitored
params <- c("mu.eta", "probs", "psi", "p", "Nsite", "Nocc.fs", "Sigma", "rho")

# MCMC settings
ni <- 20000   ;   nt <- 15   ;   nb <- 5000   ;   nc <- 3

# Call JAGS from R (ART 12 min), check traceplots and summarize posteriors
out6 <- jags(win.data, inits, params, "model6.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
par(mfrow = c(3,3))   ;   traceplot(out6, c('mu.eta', 'probs', 'Sigma', 'rho'))
print(out6, 3)

# Graphically compare some estimates between fixed- and random-effects model
par(mfrow = c(2,2))      # not shown
# Species-specific occupancy (probability scale)
plot(out5$summary[1:145,1], out6$summary[5:149,1], main = "Species-specific occupancy probability")   ;   abline(0,1)
# Species-specific detection (probability scale)
plot(out5$summary[146:290,1], out6$summary[150:294,1], main = "Species-specific detection probability")   ;   abline(0,1)
# Site-specific species richness
plot(out5$summary[291:557,1], out6$summary[295:561,1], main = "Site-specific species richness (conditional on list of 145 detected)")   ;   abline(0,1)
# Species-specific number of presences
plot(out5$summary[558:702,1], out6$summary[562:706,1], main = "Species-specific number of presences (in 267 sites)")   ;   abline(0,1)

# Estimated occupancy and detection probability for each species
plot(out6$summary[5:149,1], out6$summary[150:294,1], xlab = "Occupancy estimate", ylab = "Detection estimate", xlim = c(0,1), ylim = c(0,1), frame = F, pch = 16)
segments(out6$summary[5:149,3], out6$summary[150:294,1], out6$summary[5:149,7], out6$summary[150:294,1], col = "grey", lwd = 2)
segments(out6$summary[5:149,1], out6$summary[150:294,3], out6$summary[5:149,1], out6$summary[150:294,7], col = "grey", lwd = 2)



# 11.6.3 Modeling species-specific effects in community occupancy models
# ------------------------------------------------------------------------
# Look at distribution of body mass among 145 observed species (see errata)
mass <- MHB2014$species$body.mass[-toss.out] # Get  species mass of observed species
hist(log10(mass), breaks = 40, col = "grey")      # Look at log10
gmass <- as.numeric(log10(mass) %/% 1.3 + 1)      # size groups 1, 2 and 3
gmass[gmass == 4] <- 3                            # Mute swan is group 3, too

# Bundle and summarize data set
str( win.data <- list(ysum = ysum, g = gmass, M = nrow(ysum), J = MHB2014$sites$nsurvey, nspec = dim(ysum)[2]) )

# Specify model in BUGS language
sink("model7.txt")
cat("
model {

# Priors
for(k in 1:nspec){      # loop over species
   lpsi[k] ~ dnorm(mu.lpsi[g[k]], tau.lpsi[g[k]]) # note g-dependence now
   lp[k] ~ dnorm(mu.lp[g[k]], tau.lp[g[k]])
}

# Hyperpriors
for(g in 1:3){          # loop over groups (g)
   mu.lpsi[g] <- logit(mu.psi[g])      # everything is indexed g now
   mu.lp[g] <- logit(mu.p[g])
   mu.psi[g] ~ dunif(0,1)
   mu.p[g] ~ dunif(0,1)
   tau.lpsi[g] <- pow(sd.lpsi[g], -2)
   sd.lpsi[g] ~ dunif(0,5)
   tau.lp[g] <- pow(sd.lp[g], -2)
   sd.lp[g] ~ dunif(0,5)
}

# Ecological model for latent occurrence z (process model)
for(k in 1:nspec){      # no change at all down here in model
   logit(psi[k]) <- lpsi[k]
   for (i in 1:M) {
      z[i,k] ~ dbern(psi[k])
   }
}

# Observation model for observed data ysum
for(k in 1:nspec){      # Loop over species
   logit(p[k]) <- lp[k]
   for (i in 1:M) {
      mu.px[i,k] <- z[i,k] * p[k]  # call mu.px to avoid conflict with above
      ysum[i,k] ~ dbin(mu.px[i,k], J[i])
   }
}

# Derived quantities
for(k in 1:nspec){          # Loop over species
   Nocc.fs[k] <- sum(z[,k]) # Number of occupied sites among the 267
}
for (i in 1:M) {            # Loop over sites
   Nsite[i] <- sum(z[i,])   # Number of occurring species at each site
}
}
",fill = TRUE)
sink()

# Initial values
zst <- apply(y, c(1,3), max)
zst[is.na(zst)] <- 1
inits <- function() list(z = zst)

# Parameters monitored
params <- c("mu.psi", "mu.lpsi", "sd.lpsi", "mu.p", "mu.lp", "sd.lp")

# MCMC settings
ni <- 6000   ;   nt <- 2   ;   nb <- 2000   ;   nc <- 3

# Call JAGS from R (ART 6 min), look at convergence and summarize posteriors
out7 <- jags(win.data, inits, params, "model7.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
par(mfrow = c(3,3))   ;   traceplot(out7)
print(out7, dig = 3)


# Bundle and summarize data set
logmass <- as.numeric(log10(mass))         # Take log10 of body mass
str( win.data <- list(ysum = ysum, logmass = logmass, M = nrow(ysum), nsite = nrow(ysum),
                      J = MHB2014$sites$nsurvey, nspec = dim(ysum)[2]) )


# Specify model in BUGS language
sink("model8.txt")
cat("
model {

# Priors
for(k in 1:nspec){              # loop over species
   lpsi[k] ~ dnorm(mu.lpsi[k], tau.lpsi[k]) # now all indexed by k, not g
   tau.lpsi[k] <- 1/var.lpsi[k]
   lp[k] ~ dnorm(mu.lp[k], tau.lp[k])
   tau.lp[k] <- 1/var.lp[k]
   mu.lpsi[k] <- delta0.lpsi + delta1.lpsi * logmass[k]
   mu.lp[k] <- delta0.lp + delta1.lp * logmass[k]
   log(var.lpsi[k]) <- phi0.lpsi + phi1.lpsi * logmass[k]
   log(var.lp[k]) <- phi0.lp + phi1.lp * logmass[k]
}
# Priors for regression params for means
delta0.lpsi ~ dnorm(0, 0.01)
delta1.lpsi ~ dnorm(0, 0.01)
delta0.lp ~ dnorm(0, 0.01)
delta1.lp ~ dnorm(0, 0.01)
# Priors for regression params for variances
phi0.lpsi ~ dnorm(0, 0.01)
phi1.lpsi ~ dnorm(0, 0.01)
phi0.lp ~ dnorm(0, 0.01)
phi1.lp ~ dnorm(0, 0.01)

# Ecological model for latent occurrence z (process model)
for(k in 1:nspec){
   logit(psi[k]) <- lpsi[k]
   for (i in 1:M) {
      z[i,k] ~ dbern(psi[k])
   }
}

# Observation model for observed data ysum
for(k in 1:nspec){              # Loop over species
   logit(p[k]) <- lp[k]
   for (i in 1:M) {
      mu.p[i,k] <- z[i,k] * p[k]
      ysum[i,k] ~ dbin(mu.p[i,k], J[i])
   }
}

# Derived quantities
for(k in 1:nspec){          # Loop over species
   Nocc.fs[k] <- sum(z[,k]) # Number of occupied sites among the 267
}
for (i in 1:M) {            # Loop over sites ## see errata
   Nsite[i] <- sum(z[i,])   # Number of occurring species at each site
}
}
",fill = TRUE)
sink()

# Initial values
zst <- apply(y, c(1,3), max)
zst[is.na(zst)] <- 1
inits <- function() list(z = zst, delta0.lpsi = rnorm(1), delta1.lpsi = rnorm(1),
                         delta0.lp = rnorm(1), delta1.lp = rnorm(1), phi0.lpsi = rnorm(1),
                         phi1.lpsi = rnorm(1), phi0.lp = rnorm(1), phi1.lp = rnorm(1))

# Parameters monitored
params <- c("delta0.lpsi", "delta1.lpsi", "delta0.lp", "delta1.lp", "phi0.lpsi",
            "phi1.lpsi", "phi0.lp", "phi1.lp", "psi", "p", "Nocc.fs", "Nsite")

# MCMC settings
ni <- 12000   ;   nt <- 2   ;   nb <- 2000   ;   nc <- 3

# Call JAGS from R (ART 12 min), look at convergence and summarize posteriors
out8 <- jags(win.data, inits, params, "model8.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
par(mfrow = c(3,3))
traceplot(out8, c('delta0.lpsi', 'delta1.lpsi', 'delta0.lp', 'delta1.lp', 'phi0.lpsi', 'phi1.lpsi', 'phi0.lp', 'phi1.lp'))
print(out8, dig = 3)


# Get covariate values for prediction
predm <- seq(10, 10000,,500)        # Predict for mass of 10g to 10 kg
pred.logm <- log10(predm)

# Compute predictions (all in one array)
tmp <- out8$sims.list               # Grab simulation list
nsamp <- out8$mcmc.info$n.samples   # Number of MCMC samples
pred <- array(NA, dim = c(500, nsamp, 4)) # Array for predictions
# [,,1] mu.psi, [,,2] mu.p, [,,3] var.lpsi, [,,4] var.lp
for(i in 1:nsamp){                  # Fill array
  pred[,i,1] <- plogis(tmp$delta0.lpsi[i] + tmp$delta1.lpsi[i] * pred.logm)
  pred[,i,2] <- plogis(tmp$delta0.lp[i] + tmp$delta1.lp[i] * pred.logm)
  pred[,i,3] <- exp(tmp$phi0.lpsi[i] + tmp$phi1.lpsi[i] * pred.logm)
  pred[,i,4] <- exp(tmp$phi0.lp[i] + tmp$phi1.lp[i] * pred.logm)
}

# Plot posterior mean and a random sample of 100 from posterior of regression
selection <- sample(1:nsamp, 100)   # Choose random sample of MCMC output
par(mfrow = c(2,2), mar = c(5,5,2,2))
matplot(predm, pred[,selection,1], ylab = "Occupancy mean", xlab = "Body mass (g)", type = "l", lty = 1, lwd = 1, col = "grey", ylim = c(0, 0.4), frame = F)
lines(predm, apply(pred[,,1], 1, mean), lwd = 3, col = "blue")
matplot(predm, pred[,selection,2], ylab = "Detection mean", xlab = "Body mass (g)", type = "l", lty = 1, lwd = 1, col = "grey", ylim = c(0, 0.8), frame = F)
lines(predm, apply(pred[,,2], 1, mean), lwd = 3, col = "blue")
matplot(predm, pred[,selection,3], ylab = "Occupancy variance", xlab = "Body mass (g)", type = "l", lty = 1, lwd = 1, col = "grey", ylim = c(0, 8), frame = F)
lines(predm, apply(pred[,,3], 1, mean), lwd = 3, col = "blue")
matplot(predm, pred[,selection,4], ylab = "Detection variance", xlab = "Body mass (g)", type = "l", lty = 1, lwd = 1, col = "grey", ylim = c(0, 8), frame = F)
lines(predm, apply(pred[,,4], 1, mean), lwd = 3, col = "blue")



# 11.6.4 Modeling species richness in a two-step analysis
# ------------------------------------------------------------------------
# Extract estimates of N from model 5
N.pm <- out5$summary[291:557, 1]       # Posterior means of Nsite
N.psd <- out5$summary[291:557, 2]      # ... posterior sd's of Nsite
N.cri <- out5$summary[291:557, c(3,7)] # ... CRL's of Nsite

# Plot estimates as a function of elevation
elev <- MHB2014$sites$elev
plot(elev, N.pm, xlab = "Altitude (m a.s.l.)", ylab = "Estimated avian species richness", ylim = c(0, 70), frame = F)
segments(elev, N.cri[,1], elev, N.cri[,2], col = "grey")
lines(smooth.spline(N.pm ~ elev, w = 1 / N.psd), col = "grey", lwd = 3)


# Bundle and summarize data set
pred.ele <- (seq(200, 2750,5) - mean.ele) / sd.ele # elevation standardised
str(win.data <- list(ele = ele, N = N.pm, psd = N.psd, n = length(N.pm), pred.ele = pred.ele, npred = length(pred.ele)))

# Define model in BUGS language
sink("meta.analysis.txt")
cat("
model{

# Priors
for(v in 1:4){         # Priors for intercept and polynomial coefficients
   beta[v] ~ dnorm(0, 0.0001)
}
tau.site <- pow(sd.site, -2)
sd.site ~ dunif(0,10)

# Likelihood
for (i in 1:n){
   N[i] ~ dnorm(muN[i], tau.psd[i]) # Measurement error model for estimated N
   tau.psd[i] <- pow(psd[i], -2)    # 'Known' part of residual: meas. error
   muN[i] <- beta[1] + beta[2] * ele[i] + beta[3] * pow(ele[i],2) +
   beta[4] * pow(ele[i],3) + eps.site[i] # add another source of uncertainty
   eps.site[i] ~ dnorm(0, tau.site) # this is the usual 'residual'
}
# Get predictions for plot
for(i in 1:npred){
   Npred[i] <- beta[1] + beta[2] * pred.ele[i] + beta[3] * pow(pred.ele[i],2) + beta[4] * pow(pred.ele[i],3)
}
} # end model
",fill=TRUE)
sink()

# Initial values, params monitored, and MCMC settings
inits <- function() list(beta = rnorm(4))
params <- c("beta", "sd.site", "Npred")
ni <- 12000   ;   nt <- 10   ;   nb <- 2000   ;   nc <- 3

# Call JAGS and summarize posterior
out <- jags(win.data, inits, params, "meta.analysis.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
print(out, 3)

lines(seq(200, 2750,5), out$mean$Npred, col = "blue", lwd = 3)
matlines(seq(200,2750,5), out$summary[6:516,c(3, 7)], col = "blue", lwd = 2, lty= "dashed")


par(mfrow = c(3, 3), mar = c(5,4,3,2))
for(i in 1:267){
  for(j in 1:267){
    plot(jitter(out5$sims.list$Nsite[,i]), jitter(out5$sims.list$Nsite[,j]),
         main = paste("Joint posterior sites", i, "and", j))
    #   browser()
  }
}



# 11.7 The Dorazio-Royle (DR) community occupancy model with data augmentation
# ----------------------------------------------------------------------------


# 11.7.1 The simplest DR community model with data augmentation
# ------------------------------------------------------------------------
# Augment data set (DA part)
nz <- 150                # Number of potential species in superpopulation
M <- nspec + nz          # Size of augmented data set ('superpopulation')
yaug <- cbind(ysum, array(0, dim=c(nsite, nz))) # Add all zero histories

# Bundle and summarize data set
str( win.data <- list(yaug = yaug, nsite = nrow(ysum), nrep = MHB2014$sites$nsurvey, M = M, nspec = nspec, nz = nz) )

# Specify model in BUGS language
sink("model9.txt")
cat("
model {

# Priors to describe heterogeneity among species in community
for(k in 1:M){                  # Loop over all species in augmented list
  lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
  lp[k] ~ dnorm(mu.lp, tau.lp)
}

# Hyperpriors to describe full community
omega ~ dunif(0,1)              # Data augmentation or 'occupancy' parameter
mu.lpsi ~ dnorm(0,0.001)        # Community mean of occupancy (logit)
mu.lp ~ dnorm(0,0.001)          # Community mean of detection (logit)
tau.lpsi <- pow(sd.lpsi, -2)
sd.lpsi ~ dunif(0,5)            # Species heterogeneity in logit(psi)
tau.lp <- pow(sd.lp, -2)
sd.lp ~ dunif(0,5)              # Species heterogeneity in logit(p)

# Superpopulation process:this is the 'paramater expansion' part of PX-DA
for(k in 1:M){
  w[k] ~ dbern(omega)           # Metacommunity membership indicator
}                               # (or data augmentation variable)

# Ecological model for latent occurrence z (process model)
for(k in 1:M){
  mu.psi[k] <- w[k] * psi[k]    # species not part of community zeroed out for z
  logit(psi[k]) <- lpsi[k]
  for (i in 1:nsite) {
    z[i,k] ~ dbern(mu.psi[k])
  }
}

# Observation model for observed detection frequencies
for(k in 1:M){
  logit(p[k]) <- lp[k]
  for (i in 1:nsite) {
    mu.p[i,k] <- z[i,k] * p[k]  # non-occurring species are zeroed out for p
    yaug[i,k] ~ dbin(mu.p[i,k], nrep[i])
  }
}

# Derived quantities
for(k in 1:M){
   Nocc.fs[k] <- sum(z[,k])     # Number of occupied sites among the 267
}
for (i in 1:nsite) {
   Nsite[i] <- sum(z[i,])       # Number of occurring species at each site
}
n0 <- sum(w[(nspec+1):(nspec+nz)]) # Number of unseen species in metacommunity
Ntotal <- sum(w[])              # Total metacommunity size (= nspec + n0)
}
",fill = TRUE)
sink()

# Initial values
wst <- rep(1, nspec+nz)                   # Simply set everybody at 'occurring'
zst <- array(1, dim = c(nsite, nspec+nz)) # ditto for z
inits <- function() list(z = zst, w = wst, lpsi = rnorm(n = nspec+nz), lp = rnorm(n = nspec+nz))

# Parameters monitored
params <- c("mu.lpsi", "sd.lpsi", "mu.lp", "sd.lp", "psi", "p", "Nsite", "Ntotal", "omega", "n0")

# MCMC settings
ni <- 22000   ;   nt <- 2   ;   nb <- 2000   ;   nc <- 3

# Call JAGS from R (ART 62 min), check convergence and summarize posteriors
out9 <- jags(win.data, inits, params, "model9.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
par(mfrow = c(2,2)) ; traceplot(out9, c('mu.lpsi', 'sd.lpsi', 'mu.lp', 'sd.lp'))
print(out9, dig = 3)


# Plot posterior distribution of site-specific species richness (Nsite)
par(mfrow = c(3,3), mar = c(5,4,3,2))
for(i in 1:267){
  plot(table(out9$sims.list$Nsite[,i]), main = paste("Quadrat", i),
       xlab = "Local species richness", ylab = "", frame = F,
       xlim = c((min(C[i], out9$sims.list$Nsite[,i], na.rm = T)-2),
                max(out9$sims.list$Nsite[,i]) ))
  abline(v = C[i], col = "grey", lwd = 4)
  browser()
}

# Plot it only for a selection of sites
par(mfrow = c(3,3), mar = c(5,4,3,2))
for(i in c(9, 32, 162, 12, 27, 30, 118, 159, 250)){
  plot(table(out9$sims.list$Nsite[,i]), main = paste("Quadrat", i),
       xlab = "Local species richness", ylab = "", frame = F,
       xlim = c((min(C[i], out9$sims.list$Nsite[,i], na.rm = T)-2),
                max(out9$sims.list$Nsite[,i]) ))
  abline(v = C[i], col = "grey", lwd = 4)
}


# Plot posterior distribution of total species richness (Ntotal)
plot(table(out9$sims.list$Ntotal), main = "", ylab = "", xlab = "Avian metacommunity size in Swiss MHB survey (267 1km2 quadrats)", frame = F, xlim = c(144, 245))
abline(v = nspec, col = "grey", lwd = 4)



# 11.7.2 Dorazio-Royle community model with covariates
# ------------------------------------------------------------------------
# Augment data set: choose one of two different priors on Ntotal
nz <- 250                 # Use for vague prior on Ntotal: M = 395
nz <- 215 - nspec         # Use for informative prior on Ntotal: M = 215
yaug <- array(0, dim=c(nsite, nrep, nspec+nz)) # array with only zeroes
yaug[,,1:nspec] <- y      # copy into it the observed data

# Create same NA pattern in augmented species as in the observed species
missings <- is.na(yaug[,,1]) # e.g., third survey in high-elevation quads
for(k in (nspec+1):(nspec+nz)){
  yaug[,,k][missings] <- NA
}

# Bundle and summarize data
str(win.data <- list(y = yaug, nsite = dim(y)[1], nrep = dim(y)[2], nspec = dim(y)[3], nz = nz, M = nspec + nz, ele = ele, forest = forest, DAT = DAT, DUR = DUR) )


# Specify model in BUGS language
sink("model10.txt")
cat("
model {

# Priors
omega ~ dunif(0,1)
# Priors for species-specific effects in occupancy and detection
for(k in 1:M){
  lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)    # Hyperparams describe community
  betalpsi1[k] ~ dnorm(mu.betalpsi1, tau.betalpsi1)
  betalpsi2[k] ~ dnorm(mu.betalpsi2, tau.betalpsi2)
  betalpsi3[k] ~ dnorm(mu.betalpsi3, tau.betalpsi3)
  lp[k] ~ dnorm(mu.lp, tau.lp)
  betalp1[k] ~ dnorm(mu.betalp1, tau.betalp1)
  betalp2[k] ~ dnorm(mu.betalp2, tau.betalp2)
  betalp3[k] ~ dnorm(mu.betalp3, tau.betalp3)
}

# Hyperpriors
# For the model of occupancy
mu.lpsi ~ dnorm(0,0.01)
tau.lpsi <- pow(sd.lpsi, -2)
sd.lpsi ~ dunif(0,8)   # as always, bounds of uniform chosen by trial and error
mu.betalpsi1 ~ dnorm(0,0.1)
tau.betalpsi1 <- pow(sd.betalpsi1, -2)
sd.betalpsi1 ~ dunif(0, 4)
mu.betalpsi2 ~ dnorm(0,0.1)
tau.betalpsi2 <- pow(sd.betalpsi2, -2)
sd.betalpsi2 ~ dunif(0,2)
mu.betalpsi3 ~ dnorm(0,0.1)
tau.betalpsi3 <- pow(sd.betalpsi3, -2)
sd.betalpsi3 ~ dunif(0,2)

# For the model of detection
mu.lp ~ dnorm(0,0.1)
tau.lp <- pow(sd.lp, -2)
sd.lp ~ dunif(0, 2)
mu.betalp1 ~ dnorm(0,0.1)
tau.betalp1 <- pow(sd.betalp1, -2)
sd.betalp1 ~ dunif(0,1)
mu.betalp2 ~ dnorm(0,0.1)
tau.betalp2 <- pow(sd.betalp2, -2)
sd.betalp2 ~ dunif(0,1)
mu.betalp3 ~ dnorm(0,0.1)
tau.betalp3 <- pow(sd.betalp3, -2)
sd.betalp3 ~ dunif(0,1)

# Superpopulation process: Ntotal species sampled out of M available
for(k in 1:M){
   w[k] ~ dbern(omega)
}

# Ecological model for true occurrence (process model)
for(k in 1:M){
  for (i in 1:nsite) {
    logit(psi[i,k]) <- lpsi[k] + betalpsi1[k] * ele[i] +
      betalpsi2[k] * pow(ele[i],2) + betalpsi3[k] * forest[i]
    mu.psi[i,k] <- w[k] * psi[i,k]
    z[i,k] ~ dbern(mu.psi[i,k])
  }
}

# Observation model for replicated detection/nondetection observations
for(k in 1:M){
  for (i in 1:nsite){
    for(j in 1:nrep){
      logit(p[i,j,k]) <- lp[k] + betalp1[k] * DAT[i,j] +
        betalp2[k] * pow(DAT[i,j],2) + betalp3[k] * DUR[i,j]
      mu.p[i,j,k] <- z[i,k] * p[i,j,k]
      y[i,j,k] ~ dbern(mu.p[i,j,k])
    }
  }
}

# Derived quantities
#for(k in 1:M){
#   Nocc.fs[k] <- sum(z[,k])       # Number of occupied sites among the 267
#}
for (i in 1:nsite){
   Nsite[i] <- sum(z[i,])          # Number of occurring species at each site
}
n0 <- sum(w[(nspec+1):(nspec+nz)]) # Number of unseen species
Ntotal <- sum(w[])                 # Total metacommunity size

# Vectors to save (S for ?save?; discard posterior samples for
# all minus 1 of the potential species to save disk space)
# we do this for nz = 250 (i.e., M = 395)
lpsiS[1:(nspec+1)] <- lpsi[1:(nspec+1)]
betalpsi1S[1:(nspec+1)] <- betalpsi1[1:(nspec+1)]
betalpsi2S[1:(nspec+1)] <- betalpsi2[1:(nspec+1)]
betalpsi3S[1:(nspec+1)] <- betalpsi3[1:(nspec+1)]
lpS[1:(nspec+1)] <- lp[1:(nspec+1)]
betalp1S[1:(nspec+1)] <- betalp1[1:(nspec+1)]
betalp2S[1:(nspec+1)] <- betalp2[1:(nspec+1)]
betalp3S[1:(nspec+1)] <- betalp3[1:(nspec+1)]
}
",fill = TRUE)
sink()


# Initial values
wst <- rep(1, nspec+nz)                   # Simply set everybody at occurring
zst <- array(1, dim = c(nsite, nspec+nz)) # ditto
inits <- function() list(z = zst, w = wst, lpsi = rnorm(n = nspec+nz), betalpsi1 = rnorm(n = nspec+nz), betalpsi2 = rnorm(n = nspec+nz), betalpsi3 = rnorm(n = nspec+nz), lp = rnorm(n = nspec+nz), betalp1 = rnorm(n = nspec+nz), betalp2 = rnorm(n = nspec+nz), betalp3 = rnorm(n = nspec+nz))

# Set 1
params1 <- c("omega", "mu.lpsi", "sd.lpsi", "mu.betalpsi1", "sd.betalpsi1", "mu.betalpsi2", "sd.betalpsi2", "mu.betalpsi3", "sd.betalpsi3", "mu.lp", "sd.lp", "mu.betalp1", "sd.betalp1", "mu.betalp2", "sd.betalp2", "mu.betalp3", "sd.betalp3", "Ntotal", "Nsite")

# MCMC settings
ni <- 15000   ;   nt <- 10   ;   nb <- 5000   ;   nc <- 3

# Run JAGS, check convergence and summarize posteriors
out101 <- jags(win.data, inits, params1, "model10.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
par(mfrow = c(2, 2))
traceplot(out101, c(c("omega", "mu.lpsi", "sd.lpsi", "mu.betalpsi1", "sd.betalpsi1", "mu.betalpsi2", "sd.betalpsi2", "mu.betalpsi3", "sd.betalpsi3", "mu.lp", "sd.lp", "mu.betalp1", "sd.betalp1", "mu.betalp2", "sd.betalp2", "mu.betalp3", "sd.betalp3", "Ntotal")) )

# Set 2
params2 <- c("mu.lpsi", "sd.lpsi", "mu.betalpsi1", "sd.betalpsi1", "mu.betalpsi2", "sd.betalpsi2", "mu.betalpsi3", "sd.betalpsi3", "lpsi", "betalpsi1", "betalpsi2", "betalpsi3", "lp", "betalp1", "betalp2", "betalp3", "z", "w")
ni <- 12000   ;   nt <- 20   ;   nb <- 2000   ;   nc <- 3
out102 <- jags.basic(win.data, inits, params2, "model10.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
library(coda)
all10 <- as.matrix(out102) # Put output from 3 chains into a matrix
summary(out102)            # May take a loooong time
gelman.diag(out102)        # ditto


# Comparison of main hyperparameters when M = 215 and with M = 395
# (not all code to produce this output is shown)
print(cbind(out10.215$summary[1:17,c(1:3, 7)], out10.395$summary[1:17, c(1:3, 7)]), 2)


out10 <- out101


par(mfrow = c(1,2))       # Fig. 11-16
psi.sample <- plogis(rnorm(10^6, mean = out10$mean$mu.lpsi, sd = out10$mean$sd.lpsi))
p.sample <- plogis(rnorm(10^6, mean = out10$mean$mu.lp, sd = out10$mean$sd.lp))
hist(psi.sample, freq = F, breaks = 50, col = "grey", xlab = "Species occupancy probability", ylab = "Density", main = "")
hist(p.sample, freq = F, breaks = 50, col = "grey", xlab = "Species detection probability", ylab = "Density", main = "")
summary(psi.sample)   ;   summary(p.sample)

par(mfrow = c(2,4))  # Among-species variability in parameters (not shown)
hist(out10$sims.list$sd.lpsi, breaks = 100, col = "grey", xlim = c(0,6), main = "Occupancy: intercept")
abline(v = mean(out10$sims.list$sd.lpsi), col = "blue", lwd = 3)
hist(out10$sims.list$sd.betalpsi1, breaks = 100, col = "grey", xlim = c(0,3), main = "Occupancy: linear effect of elevation")
abline(v = mean(out10$sims.list$sd.betalpsi1), col = "blue", lwd = 3)
hist(out10$sims.list$sd.betalpsi2, breaks = 100, col = "grey", xlim = c(0,3), main = "Occupancy: quadratic effect of elevation")
abline(v = mean(out10$sims.list$sd.betalpsi2), col = "blue", lwd = 3)
hist(out10$sims.list$sd.betalpsi3, breaks = 100, col = "grey", xlim = c(0,3), main = "Occupancy: linear effect of forest cover")
abline(v = mean(out10$sims.list$sd.betalpsi3), col = "blue", lwd = 3)
hist(out10$sims.list$sd.lp, breaks = 100, col = "grey", xlim = c(0,2), main = "Detection: intercept")
abline(v = mean(out10$sims.list$sd.lp), col = "blue", lwd = 3)
hist(out10$sims.list$sd.betalp1, breaks = 100, col = "grey", xlim = c(0,1), main = "Detection: linear effect of survey date")
abline(v = mean(out10$sims.list$sd.betalp1), col = "blue", lwd = 3)
hist(out10$sims.list$sd.betalp2, breaks = 100, col = "grey", xlim = c(0,1), main = "Detection: quadratic linear effect of survey date")
abline(v = mean(out10$sims.list$sd.betalp2), col = "blue", lwd = 3)
hist(out10$sims.list$sd.betalp3, breaks = 100, col = "grey", xlim = c(0,1), main = "Detection: linear effect of survey duration")
abline(v = mean(out10$sims.list$sd.betalp3), col = "blue", lwd = 3)


# Visualize covariate mean relationships for the average species
o.ele <- seq(200, 2500,,500)               # Get covariate values for prediction
o.for <- seq(0, 100,,500)
o.dat <- seq(15, 120,,500)
o.dur <- seq(100, 420,,500)
ele.pred <- (o.ele - mean.ele) / sd.ele
for.pred <- (o.for - mean.forest) / sd.forest
dat.pred <- (o.dat - mean.date) / sd.date
dur.pred <- (o.dur - mean.dur) / sd.dur

# Predict occupancy for elevation and forest and detection for date and duration
# Put all fourpredictions into a single
str( tmp <- out10$sims.list )              # grab MCMC samples
nsamp <- length(tmp[[1]])    # number of mcmc samples
predC <- array(NA, dim = c(500, nsamp, 4)) # "C" for 'community mean'
for(i in 1:nsamp){
  predC[,i,1] <- plogis(tmp$mu.lpsi[i] + tmp$mu.betalpsi1[i] * ele.pred +
                          tmp$mu.betalpsi2[i] * ele.pred^2 )
  predC[,i,2] <- plogis(tmp$mu.lpsi[i] + tmp$mu.betalpsi3[i] * for.pred)
  predC[,i,3] <- plogis(tmp$mu.lp[i] + tmp$mu.betalp1[i] * dat.pred +
                          tmp$mu.betalp2[i] * dat.pred^2 )
  predC[,i,4] <- plogis(tmp$mu.lp[i] + tmp$mu.betalp3[i] * dur.pred)
}

# Get posterior means and 95% CRIs and plot (Fig. 11?17)
pmC <- apply(predC, c(1,3), mean)
criC <- apply(predC, c(1,3), function(x) quantile(x, prob = c(0.025, 0.975)))

par(mfrow = c(2, 2))
plot(o.ele, pmC[,1], col = "blue", lwd = 3, type = 'l', lty = 1, frame = F, ylim = c(0, 0.05), xlab = "Elevation (m a.s.l)", ylab = "Community mean occupancy")
matlines(o.ele, t(criC[,,1]), col = "grey", lty = 1)
plot(o.for, pmC[,2], col = "blue", lwd = 3, type = 'l', lty = 1, frame = F, ylim = c(0, 0.05), xlab = "Forest cover", ylab = "Community mean occupancy")
matlines(o.for, t(criC[,,2]), col = "grey", lty = 1)
plot(o.dat, pmC[,3], col = "blue", lwd = 3, type = 'l', lty = 1, frame = F, ylim = c(0.2, 0.8), xlab = "Survey date", ylab = "Community mean detection")
matlines(o.dat, t(criC[,,3]), col = "grey", lty = 1)
plot(o.dur, pmC[,4], col = "blue", lwd = 3, type = 'l', lty = 1, frame = F, ylim = c(0.2, 0.8), xlab = "Survey duration", ylab = "Community mean detection")
matlines(o.dur, t(criC[,,4]), col = "grey", lty = 1)


# Plot posterior distribution of site-specific species richness (Nsite)
par(mfrow = c(3,3), mar = c(5,4,3,2))
for(i in 1:267){
  plot(table(out10$sims.list$Nsite[,i]), main = paste("Quadrat", i),
       xlab = "Local species richness", ylab = "", frame = F,
       xlim = c((min(C[i], out10$sims.list$Nsite[,i], na.rm = T)-2),
                max(out10$sims.list$Nsite[,i]) ))
  abline(v = C[i], col = "grey", lwd = 4)
  browser()
}

# Plot it only for a selection of sites (Fig. 11-18)
par(mfrow = c(3,3), mar = c(5,4,3,2))
for(i in c(9, 32, 162, 12, 27, 30, 118, 159, 250)){
  plot(table(out10$sims.list$Nsite[,i]), main = paste("Quadrat", i),
       xlab = "Local species richness", ylab = "", frame = F,
       xlim = c((min(C[i], out10$sims.list$Nsite[,i], na.rm = T)-2),
                max(out10$sims.list$Nsite[,i]) ))
  abline(v = C[i], col = "grey", lwd = 4)
}

# Plot Nsite estimates under models 9 & 10 vs. elevation (Fig. 11-19)
offset <- 30    # Set off elevation for better visibility
plot(elev, out9$mean$Nsite, xlab = "Elevation (metres)", ylab = "Community size estimate (Nsite)", frame = F, ylim = c(0,60), pch = 16) # black: model 9
lines(smooth.spline(out9$mean$Nsite ~ elev), lwd = 3)
points(elev+offset, out10$mean$Nsite, pch = 16, col = "blue") # red: model 10
lines(smooth.spline(out10$mean$Nsite ~ elev), lwd = 3, col = "blue")


str(all10)                    # look at the MCMC output
pm <- apply(all10, 2, mean)    # Get posterior means and 95% CRIs
cri <- apply(all10, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs


# Effects of date (linear and quadratic) and of duration on detection
#par(mfrow = c(1,3), cex.lab = 1.3, cex.axis = 1.3) # Can put all three in one
par(mfrow = c(1,2), cex.lab = 1.3, cex.axis = 1.3)
# Date linear (Fig. 11 ? 20 left)
plot(pm[1:145], 1:145, xlim = c(-1.5, 1.5), xlab = "Parameter estimate", ylab = "Species number", main = "Effect of date (linear) on detection", pch = 16)
abline(v = 0, lwd = 2, col = "black")
segments(cri[1, 1:145], 1:145, cri[2, 1:145], 1:145, col = "grey", lwd = 1)
sig1 <- (cri[1, 1:145] * cri[2, 1:145]) > 0
segments(cri[1, 1:145][sig1 == 1], (1:145)[sig1 == 1], cri[2, 1:145][sig1 == 1], (1:145)[sig1 == 1], col = "blue", lwd = 2)
abline(v = out101$summary[11,1], lwd = 3, col = "red")
abline(v = out101$summary[11,c(3,7)], lwd = 2, col = "red", lty = 2)



# Date quadratic (not shown)
plot(pm[216:360], 1:145, xlim = c(-1.5, 1.5), xlab = "Parameter estimate", ylab = "Species number", main = "Effect of date (quadratic) on detection", pch = 16)
abline(v = 0, lwd = 2, col = "black")
segments(cri[1, 216:360], 1:145, cri[2, 216:360], 1:145, col = "grey", lwd = 1)
sig2 <- (cri[1, 216:360] * cri[2, 216:360]) > 0
segments(cri[1, 216:360][sig2 == 1], (1:145)[sig2 == 1], cri[2, 216:360][sig2 == 1], (1:145)[sig2 == 1], col = "blue", lwd = 2)
abline(v = out101$summary[13,1], lwd = 3, col = "red")
abline(v = out101$summary[13, c(3,7)], lwd = 3, col = "red", lty = 2)


# Survey duration (Fig. 11-20 right)
plot(pm[431:575], 1:145, xlim = c(-0.5, 1), xlab = "Parameter estimate", ylab = "Species number", main = "Effect of survey duration on detection", pch = 16)
abline(v = 0, lwd = 2, col = "black")
segments(cri[1, 431:575], 1:145, cri[2, 431:575], 1:145, col = "grey", lwd = 1)
sig3 <- (cri[1, 431:575] * cri[2, 431:575]) > 0
segments(cri[1, 431:575][sig3 == 1], (1:145)[sig3 == 1], cri[2, 431:575][sig3 == 1], (1:145)[sig3 == 1], col = "blue", lwd = 2)
abline(v = out101$summary[15,1], lwd = 3, col = "red")
abline(v = out101$summary[15, c(3,7)], lwd = 3, col = "red", lty = 2)


# Effects of elevation (linear and quadratic) and of forest on occupancy
# par(mfrow = c(1,3), cex.lab = 1.3, cex.axis = 1.3) # can do all in one
# Effect of elevation (linear) on occupancy probability (Fig. 11-21)
plot(pm[646:790], 1:145, xlim = c(-8, 8), xlab = "Parameter estimate", ylab = "Species number", main = "Effect of elevation (linear) on occupancy", pch = 16)
abline(v = 0, lwd = 2, col = "black")
segments(cri[1, 646:790], 1:145, cri[2, 646:790], 1:145, col = "grey", lwd = 1)
sig4 <- (cri[1, 646:790] * cri[2, 646:790]) > 0
segments(cri[1, 646:790][sig4 == 1], (1:145)[sig4 == 1], cri[2, 646:790][sig4 == 1], (1:145)[sig4 == 1], col = "blue", lwd = 2)
abline(v = out101$summary[3,1], lwd = 3, col = "red")
abline(v = out101$summary[3,c(3,7)], lwd = 3, col = "red", lty = 2)


# Effect of elevation (quadratic) on occupancy probability (Fig. 11-22)
plot(pm[861:1005], 1:145, xlim = c(-4, 2), xlab = "Parameter estimate", ylab = "Species number", main = "Effect of elevation (quadratic) on occupancy", pch = 16)
abline(v = 0, lwd = 2, col = "black")
segments(cri[1, 861:1005], 1:145, cri[2, 861:1005], 1:145, col = "grey", lwd=1)
sig5 <- (cri[1, 861:1005] * cri[2, 861:1005]) > 0
segments(cri[1, 861:1005][sig5 == 1], (1:145)[sig5 == 1], cri[2, 861:1005][sig5 == 1], (1:145)[sig5 == 1], col = "blue", lwd = 2)
abline(v = out101$summary[5,1], lwd = 3, col = "red")
abline(v = out101$summary[5,c(3,7)], lwd = 3, col = "red", lty = 2)


# Effect of forest (linear) on occupancy probability (Fig. 11-23)
plot(pm[1076:1220], 1:145, xlim = c(-3, 4), xlab = "Parameter estimate", ylab = "Species number", main = "Effect of forest cover on occupancy", pch = 16)
abline(v = 0, lwd = 2, col = "black")
segments(cri[1, 1076:1220], 1:145, cri[2, 1076:1220],1:145, col = "grey", lwd=1)
sig6 <- (cri[1, 1076:1220] * cri[2, 1076:1220]) > 0
segments(cri[1, 1076:1220][sig6 == 1], (1:145)[sig6 == 1], cri[2, 1076:1220][sig6 == 1], (1:145)[sig6 == 1], col = "blue", lwd = 2)
abline(v = out101$summary[7,1], lwd = 3, col = "red")
abline(v = out101$summary[7,c(3,7)], lwd = 3, col = "red", lty = 2)
negsig6 <- (cri[1, 1076:1220] < 0 & cri[2, 1076:1220] < 0) == 1 # sig negative
possig6 <- (cri[1, 1076:1220] > 0 & cri[2, 1076:1220] > 0) == 1 # sig positive


# Predict detection for date and duration and occupancy for elevation and forest
# for each of the 145 observed species
predS <- array(NA, dim = c(500, nspec, 4))   # covariate value x species x response, "S" for 'species'
p.coef <- cbind(lp=pm[1292:1436], betalp1 = pm[1:145], betalp2 = pm[216:360], betalp3 = pm[431:575])
psi.coef <- cbind(lpsi=pm[1507:1651], betalpsi1 = pm[646:790], betalpsi2 = pm[861:1005], betalpsi3 = pm[1076:1220])

for(i in 1:nspec){          # Loop over 145 observed species
  predS[,i,1] <- plogis(p.coef[i,1] + p.coef[i,2] * dat.pred +
                          p.coef[i,3] * dat.pred^2 )     # p ~ date
  predS[,i,2] <- plogis(p.coef[i,1] + p.coef[i,4] * dur.pred) # p ~ duration
  predS[,i,3] <- plogis(psi.coef[i,1] + psi.coef[i,2] * ele.pred +
                          psi.coef[i,3] * ele.pred^2 )     # psi ~ elevation
  predS[,i,4] <- plogis(psi.coef[i,1] + psi.coef[i,4] * for.pred) # psi ~ forest
}

# Plots for detection probability and survey date and duration (Fig. 11-24)
par(mfrow = c(1,2), cex.lab = 1.3, cex.axis = 1.3)
plot(o.dat, predS[,1,1], lwd = 3, type = 'l', lty = 1, frame = F,
     ylim = c(0, 1), xlab = "Survey date (1 = 1 April)",
     ylab = "Detection probability")
for(i in 2:145){
  lines(o.dat, predS[,i,1], col = i, lwd = 3)
}

plot(o.dur, predS[,1,2], lwd = 3, type = 'l', lty = 1, frame = F,
     ylim = c(0, 1), xlab = "Survey duration (min)",
     ylab = "Detection probability")
for(i in 2:145){
  lines(o.dur, predS[,i,2], col = i, lwd = 3)
}


# Plots for occupancy probability and elevation and forest cover (Fig. 11-25)
par(mfrow = c(1,2), cex.lab = 1.3, cex.axis = 1.3)
plot(o.ele, predS[,1,3], lwd = 3, type = 'l', lty = 1, frame = F,
     ylim = c(0, 1), xlab = "Elevation (m a.s.l.)",
     ylab = "Occupancy probability")
for(i in 2:145){
  lines(o.ele, predS[,i,3], col = i, lwd = 3)
}

plot(o.for, predS[,1,4], lwd = 3, type = 'l', lty = 1, frame = F,
     ylim = c(0, 1), xlab = "Forest cover (%)", ylab = "Occupancy probability")
for(i in 2:145){
  lines(o.for, predS[,i,4], col = i, lwd = 3)
}



# 11.8 Inferences based on the estimated Z matrix: similarity among sites and species
# -----------------------------------------------------------------------------------


# Plug MCMC samples for full z matrix into 3D array
str(all10)
nsite <- 267
nspec <- 215
nsamp <- dim(all10)[1]        # 1200 MCMC samples
z <- array(NA, dim = c(nsite, nspec, nsamp))
Jacc <- array(NA, dim = c(nsite, nspec, nsamp))
for(j in 1:nsamp){    # Fill z matrix by column (default)
  cat(paste("\nMCMC sample", j, "\n"))
  z[,,j] <- all10[j, 1937:59341]
}

# Restrict computations to observed species
zobs <- z[,1:145,]      # Species 1 to 145

# Compute Jaccard index for sites and for species
Jsite <- array(NA, dim = c(nsite, nsamp))
Jspec <- array(NA, dim = c(145, nsamp))


# Choose reference site and species for Jaccard indices
ref.site <- 1         # Just choose first site
ref.species <- 13     # European Sparrowhawk (check object 'obs.occ')

# Get posterior distributions for Jsite and Jspec (for references)
for(k in 1:nsamp){
  for(i in 1:nsite){ # Jaccard index for sites (in terms of shared species)
    Jsite[i,k] <- sum(zobs[ref.site,,k] * zobs[i,,k]) /
      (sum(zobs[ref.site,,k]) + sum(zobs[i,,k]) -
         sum(zobs[ref.site,,k] * zobs[i,,k]))
  }
  for(i in 1:(nspec-nz)){ # Jacc. index for species (in terms of shared sites)
    Jspec[i,k] <- sum(zobs[,ref.species,k] * zobs[,i,k]) /
      (sum(zobs[,ref.species,k]) + sum(zobs[,i,k]) -
         sum(zobs[,ref.species,k] * zobs[,i,k]))
  }
}
# NA's arise when a site has no species or a species no sites

# Get posterior means, standard deviations and 95% CRI
# Jaccard index for sites, compared to reference site 1
pm <- apply(Jsite, 1, mean, na.rm = TRUE)  # Post. mean of Jsite wrt. site 1
psd <- apply(Jsite, 1, sd, na.rm = TRUE)   # Post. mean of Jsite wrt. site 1
cri <- apply(Jsite, 1, function(x) quantile(x, prob = c(0.025, 0.975), na.rm = TRUE)) # CRI
cbind('post. mean' = pm, 'post. sd' = psd, '2.5%' = cri[1,], '97.5%' = cri[2,])


# Make a map of Jaccard site indices (Fig. 11-26)
x <- 3        # poportional size of plotting symbol
plot(MHB2014$sites$coordx, MHB2014$sites$coordy, xlab = "x coordinate", ylab = "y coordinate", cex = x*pm, asp = 1, pch = 16)
points(MHB2014$sites$coordx[which(pm == 1)], MHB2014$sites$coordy[which(pm == 1)], cex = x*pm, col = "red", pch = 16)


# Jaccard index for species, compared to reference species
# (species 13, European Sparrowhawk)
pm <- apply(Jspec, 1, mean, na.rm = TRUE)  # Post. mean of Jspec wrt. species 1
psd <- apply(Jspec, 1, sd, na.rm = TRUE)   # Post. mean of Jspec wrt. species 1
cri <- apply(Jspec, 1, function(x) quantile(x, prob = c(0.025, 0.975), na.rm = TRUE)) # CRI
tmp <- cbind('post. mean' = pm, 'post. sd' = psd, '2.5%' = cri[1,], '97.5%' = cri[2,])
rownames(tmp) <- names(obs.occ)
print(tmp])          # print in systematic order
print(tmp[rev(order(tmp[,1])),]) # print in order of decreasing Jacc. values
plot(1:145, tmp[rev(order(tmp[,1])),1])   # can also plot




# 11.9 Species richness maps and species accumulation curves
# ------------------------------------------------------------------------


# Get Swiss landscape data and standardise covariates as for model 10
library(unmarked)
data(Switzerland)
ch <- Switzerland
ELE <- (ch$elevation - mean.ele) / sd.ele
FOREST <- (ch$forest - mean.forest) / sd.forest

nsamp <- nrow(all10)            # 1200   ..... far too many
nkm2 <- length(ch[[1]])         # 42275, that's a LOT!
select.samp <- sort(sample(1:nsamp, 50)) # Chose random sample of 50
nsamp <- length(select.samp)    # new sample size 50

# Create posterior predictive distribution for Z for Swiss landscape
str( zCH <- array(NA, dim = c(nkm2, 215, nsamp)) ) # BIG array !
W <- all10[,1722:1936]          # Grab MCMC samples from w
LPSI <- all10[,1507:1721]       # Grab MCMC samples from logit(psi)
BETALPSI1 <- all10[,646:860]    # Grab MCMC samples from betalpsi1
BETALPSI2 <- all10[,861:1075]   # Grab MCMC samples from betalpsi2
BETALPSI3 <- all10[,1076:1290]  # Grab MCMC samples from betalpsi3
for(i in 1:nkm2){               # takes about 5 mins !
  cat(paste("\nQuadrat", i, "\n"))
  for(u in 1:length(select.samp)){
    psi <- W[select.samp[u],] * plogis(LPSI[select.samp[u],] +
                                         BETALPSI1[select.samp[u],] * ELE[i] +
                                         BETALPSI2[select.samp[u],] * ELE[i]^2 +
                                         BETALPSI3[select.samp[u],] * FOREST[i] )
    zCH[i,,u] <- rbinom(215, 1, psi)
  }
}

# Compute posterior distribution of species richness by collapsing z array
SR <- apply(zCH, c(1,3), sum)   # posterior distribution
pmSR <- apply(SR, 1, mean)      # posterior mean
sdSR <- apply(SR, 1, sd)        # posterior standard deviation


library(raster)
library(rgdal)
par(mfrow = c(1,2), mar = c(2,2,3,5))
# Posterior mean map
r1 <- rasterFromXYZ(data.frame(x = ch$x, y = ch$y, z = pmSR))
elev <- rasterFromXYZ(cbind(ch$x, ch$y,ch$elevation))
elev[elev > 2250] <- NA         # Mask areas > 2250 m a.s.l.
r1 <- mask(r1, elev)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
plot(r1, col = mapPalette(100), axes = F, box = FALSE, main ="")
lakes <- readOGR(".", "lakes")
rivers <- readOGR(".", "rivers")
border <- readOGR(".", "border")
plot(rivers, col = "dodgerblue", add = TRUE)
plot(border, col = "transparent", lwd = 1.5, add = TRUE)
plot(lakes, col = "skyblue", border = "royalblue", add = TRUE)

# Posterior standard deviation map
r1 <- rasterFromXYZ(data.frame(x = ch$x, y = ch$y, z = sdSR))
elev <- rasterFromXYZ(cbind(ch$x, ch$y,ch$elevation))
elev[elev > 2250] <- NA         # Mask areas > 2250 m a.s.l.
r1 <- mask(r1, elev)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
plot(r1, col = mapPalette(100), axes = F, box = FALSE, main ="")
lakes <- readOGR(".", "lakes")
rivers <- readOGR(".", "rivers")
border <- readOGR(".", "border")
plot(rivers, col = "dodgerblue", add = TRUE)
plot(border, col = "transparent", lwd = 1.5, add = TRUE)
plot(lakes, col = "skyblue", border = "royalblue", add = TRUE)


# Get 3,000 posterior samples of omega, and the mean and sd hyperparameters
omega <- out101$sims.list$omega
mu.lpsi <- out101$sims.list$mu.lpsi
str( sd.lpsi <- out101$sims.list$sd.lpsi )    # Confirms we have 3,000 draws

# compute posterior predictions of species occurrence probabilities
nsites <- 100
ndraws <- length(omega)
Nmax <- 215
psi <- matrix(NA, nrow=ndraws, ncol=Nmax)
for (i in 1:ndraws) {
  w <- rbinom(215, 1, omega[i])
  psi[i,] <- w * plogis(rnorm(Nmax, mean = mu.lpsi[i], sd=sd.lpsi[i]))
}

# compute posterior predictions of species presence at each site
z <- array(NA, dim=c(ndraws, Nmax, nsites))
for (i in 1:ndraws) {
  for (j in 1:Nmax) {
    z[i,j, ] <- rbinom(nsites, size=1, prob=psi[i,j])
  }
}

# compute posterior predictions of cumulative number of species present
Ntot <- matrix(NA, nrow=ndraws, ncol=nsites)
for (i in 1:ndraws) {
  for (j in 1:nsites) {
    zsum <- rep(NA, Nmax)
    if (j>1) {
      zsum <- apply(z[i, , 1:j], 1, sum)
    }
    else {
      zsum <- z[i, , 1]
    }
    Ntot[i,j] <- sum(zsum>0)
  }
}                        # takes about 4 min

# compute summary stats of species accumulation curve
nSpeciesPresent <- matrix(NA, nrow=3, ncol=nsites)
for (j in 1:nsites) {
  x <- Ntot[,j]
  nSpeciesPresent[1, j] <- mean(x)
  nSpeciesPresent[2:3, j] <- quantile(x, probs=c(0.05, 0.95))
}

# Plot species accumulation curve
ylim = c(min(nSpeciesPresent[2,]), max(nSpeciesPresent[3,]))
plot(1:nsites, nSpeciesPresent[1,], pch=16, ylim=ylim, type="b",
     xlab="Number of sample locations", ylab="Number of occurring species",
     las=1, cex.axis=1.2, cex.lab=1.5, cex=1.2, frame = F)
segments(1:nsites, nSpeciesPresent[2,], 1:nsites, nSpeciesPresent[3,])



# 11.10 Community N-mixture models
# ------------------------------------------------------------------------

# Organize counts in 3D array: site x rep x species
Yc <- MHB2014$counts
str(Yc)

# Observed maximum and mean maximum count per species
tmp <- apply(Yc, c(1,3), max, na.rm = TRUE)
tmp[tmp == -Inf] <- NA         # 1 quadrat with NA data in 2014
sort(round(meanmax <- apply(tmp, 2, mean, na.rm = TRUE), 3)) # mean of max
sort(obs.max.C <- apply(tmp, 2, max, na.rm = TRUE))          # max

# Plot observed species abundance distribution
plot(sort(meanmax), xlab = "Species number", ylab = "Mean maximum count")

# Spatio-temporal patterns in counts (mean over sites)
tmp <- apply(Yc, c(2,3), mean, na.rm = TRUE)
matplot(log10(tmp+0.1), type = "l", lty = 1, lwd = 3, xlab = "MHB survey 1 - 3", ylab = "log10 of mean count over sites", frame = F, cex.lab = 1.3, cex.axis = 1.3)

# Drop data from 13 species not observed in 2014
toss.out <- which(obs.max.C == 0)   # list of species not seen
Yc <- Yc[,,-toss.out]               # toss them out
obs.max.C <- obs.max.C[-toss.out]
( nspec <- dim(Yc)[3] )             # Redefine nspec as 145

# So here are our data
str(Yc)
plot(table(Yc))   # Extremely skewed distribution of observed counts


# Bundle and summarize data set
str(win.data <- list(Yc = Yc, nsite = dim(Yc)[1], nrep = dim(Yc)[2],
                     nspec = dim(Yc)[3], ele = ele, forest = forest, DAT = DAT, DUR = DUR))


# Specify model in BUGS language
sink("model11.txt")
cat("
model {

# Community priors (with hyperparameters) for species-specific parameters
for(k in 1:nspec){
  phi[k] ~ dunif(0,1)                              # Zero-inflation
  alpha0[k] ~ dnorm(mu.alpha0, tau.alpha0)         # Detection intercepts
  beta0[k] ~ dnorm(mu.beta0, tau.beta0)            # Abundance intercepts
  for(v in 1:3){
    alpha[k, v] ~ dnorm(mu.alpha[v], tau.alpha[v]) # Slopes detection
    beta[k, v] ~ dnorm(mu.beta[v], tau.beta[v])    # Slopes abundance
  }
}

# Hyperpriors for community hyperparameters
# abundance model
mu.beta0 ~ dunif(-1, 2)
tau.beta0 <- pow(sd.beta0, -2)
sd.beta0 ~ dunif(0, 3)
for(v in 1:3){
  mu.beta[v] ~ dunif(-1.5, 1)
  tau.beta[v] <- pow(sd.beta[v], -2)
}
sd.beta[1] ~ dunif(0, 3)
sd.beta[2] ~ dunif(0, 1.5)
sd.beta[3] ~ dunif(0, 1)

# detection model
mu.alpha0 ~ dunif(-2, 0)
tau.alpha0 <- pow(sd.alpha0, -2)
sd.alpha0 ~ dunif(0, 2)
for(v in 1:3){
  mu.alpha[v] ~ dunif(-0.5, 0.5)
  tau.alpha[v] <- pow(sd.alpha[v], -2)
}
sd.alpha[1] ~ dunif(0, 0.8)
sd.alpha[2] ~ dunif(0, 0.5)
sd.alpha[3] ~ dunif(0, 0.3)

# Ecological model for true abundance (process model)
for(k in 1:nspec){
  for (i in 1:nsite){
    a[i,k] ~ dbern(phi[k])   # zero-inflation
    N[i,k] ~ dpois(a[i,k] * lambda[i,k])
    log(lambda[i,k]) <- beta0[k] + beta[k,1] * ele[i] +
      beta[k,2] * pow(ele[i],2) + beta[k,3] * forest[i]
    # Compute presence/absence matrix z (for N > 0) from latent abundance
    z[i,k] <- step(N[i,k]-1)  # returns TRUE if N >= 0
  }
}

# Observation model for replicated counts
for(k in 1:nspec){
  for (i in 1:nsite){
    for (j in 1:nrep){
      Yc[i,j,k] ~ dbin(p[i,j,k], N[i,k])
      logit(p[i,j,k]) <- alpha0[k] + alpha[k,1] * DAT[i,j] +
        alpha[k,2] * pow(DAT[i,j],2) + alpha[k,3] * DUR[i,j]
    }
  }
}

# Other derived quantities
for(k in 1:nspec){
  mlambda[k] <- phi[k] * exp(beta0[k]) # Expected abundance on natural scale
  logit(mp[k]) <- alpha0[k]     # Mean detection on natural scale
  Nocc.fs[k] <- sum(z[,k])      # Number of occupied sites among the 267
}
for (i in 1:nsite) {
  Nsite[i] <- sum(z[i,])        # Number of occurring species at each site
}
}
",fill = TRUE)
sink()

# Initial values
ast <- matrix(rep(1, nspec*nsite), nrow = nsite)
some.more <- 5          # May have to play with this until JAGS is happy
Nst <- apply(Yc, c(1,3), max, na.rm = T) + some.more
Nst[Nst == '-Inf'] <- 20          # May have to play with this, too
Nst <- Nst
inits <- function()list(a = ast, N = Nst)

# OR: use inits at earlier solutions (greatly speeds up convergence)
pm <- out11$mean     # Pull out posterior means from earlier run
inits <- function() list(a = ast, N = Nst, alpha0 = rnorm(nspec), beta0 = rnorm(nspec), alpha = matrix(rnorm(n = nspec*3), ncol = 3), beta = matrix(rnorm(n = nspec*3), ncol = 3), mu.beta0 = pm$mu.beta0, sd.beta0 = pm$sd.beta0, mu.beta = pm$mu.beta, sd.beta = pm$sd.beta, mu.alpha0 = pm$mu.alpha0, sd.alpha0 = pm$sd.alpha0, mu.alpha = pm$mu.alpha, sd.alpha = pm$sd.alpha )

# Parameters monitored
params <- c("phi", "mp", "mlambda", "alpha0", "beta0", "alpha", "beta", "mu.beta0", "sd.beta0", "mu.beta", "sd.beta", "mu.alpha0", "sd.alpha0", "mu.alpha", "sd.alpha", "Nsite")

# MCMC settings
ni <- 60000   ;   nt <- 30   ;   nb <- 30000   ;   nc <- 3

# Call JAGS from R (BRT XXX min), check convergence and summarize posteriors
out11 <- jags(win.data, inits, params, "model11.txt", n.chains = nc,
              n.thin = nt, n.iter = ni, n.burnin = nb, parallel = FALSE)
par(mfrow = c(3,3))   ;    traceplot(out11, c("mu.beta0", "sd.beta0", "mu.beta", "sd.beta", "mu.alpha0", "sd.alpha0", "mu.alpha", "sd.alpha") )
print(out11, 2)


summary(p.sample <- plogis(rnorm(10^6, mean = -1.170, sd = 0.980)) )
hist(p.sample, breaks = 50, col = "grey", xlab = "Per-individual detection probability", freq = FALSE)


# Predict detection for date and duration and occupancy for elevation and forest
# for each of the 145 observed species
predI <- array(NA, dim = c(500, nspec, 4))   # covariate value x species x response, "I" for 'individual' (as opposed to 'species' in model 10)
pm <- out11$mean            # Grab posterior means from model 11
for(i in 1:nspec){          # Loop over 145 observed species
  predI[,i,1] <- plogis(pm$alpha0[i] + pm$alpha[i,1] * dat.pred +
                          pm$alpha[i,2] * dat.pred^2 )     # p ~ date
  predI[,i,2] <- plogis(pm$alpha0[i] + pm$alpha[i,3] * dur.pred) # p ~ duration
  predI[,i,3] <- pm$phi[i] * exp(pm$beta0[i] + pm$beta[i,1] * ele.pred +
                                   pm$beta[i,2] * ele.pred^2 )     # psi ~ elevation
  predI[,i,4] <- pm$phi[i] * exp(pm$beta0[i] + pm$beta[i,3] * for.pred) # psi ~ forest
}

# Plots for detection probability and survey date and duration (Fig. 11-29)
par(mfrow = c(1,2), cex.lab = 1.3, cex.axis = 1.3)
plot(o.dat, predI[,1,1], lwd = 3, type = 'l', lty = 1, frame = F,
     ylim = c(0, 1), xlab = "Survey date (1 = 1 April)",
     ylab = "Per-individual detection probability")
for(i in 2:145){
  lines(o.dat, predI[,i,1], col = i, lwd = 3)
}

plot(o.dur, predI[,1,2], lwd = 3, type = 'l', lty = 1, frame = F,
     ylim = c(0, 1), xlab = "Survey duration (min)",
     ylab = "Per-individual detection probability")
for(i in 2:145){
  lines(o.dur, predI[,i,2], col = i, lwd = 3)
}


# Plots for expected abundance and elevation and forest cover (Fig. 11-30)
par(mfrow = c(1,2), cex.lab = 1.3, cex.axis = 1.3)
plot(o.ele, predI[,1,3], lwd = 3, type = 'l', lty = 1, frame = F,
     ylim = c(0, 60), xlab = "Elevation (m a.s.l.)",ylab = "Expected abundance")
for(i in 2:145){
  lines(o.ele, predI[,i,3], col = i, lwd = 3)
}

plot(o.for, predI[,1,4], lwd = 3, type = 'l', lty = 1, frame = F,
     ylim = c(0, 60), xlab = "Forest cover (%)", ylab = "Expected abundance")
for(i in 2:145){
  lines(o.for, predI[,i,4], col = i, lwd = 3)
}