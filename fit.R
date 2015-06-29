### code to fit the statistical model using NIMBLE

# note - to create an MCMC for a second model within the same R session you may need to give all the related model and MCMC objects new names as overwriting old names may cause problems

source("config")

source("model.R")

setwd(dataDir)
setwd('lyford')
load('data.Rda')

tsplot <- function(x) plot(seq_along(x), x, type = 'l')

library(nimble)
# library(nimble,lib.loc='/tmp/nimble/devel')

# notes on steps for using NIMBLE:

#1 - create the model with nimbleModel() - takes a couple minutes
#2 - configure the MCMC samplers with configureMCMC() - takes about 15 sec.
#3 - create the MCMC algorithm to run in R with buildMCMC() - takes a couple minutes
#4 - compile the model to C++ with compileNimble() - takes 2-3 minutes
#5 - compile the MCMC to run in C++ with compileNimble() - takes
#6 - run the C++ version of the MCMC - 5k iterations of the more complicated models takes 10-15 minutes

# the time for step #5 is longer than it should be - I'm waiting for Perry to look into this; also it uses a bunch of memory that it shouldn't - 5-6 Gb so do on a machine with a bunch of RAM

# initial model
m <- nimbleModel(body(ringModel),
                                   constants = list(n = n, nT = nT, nDBH = nDBH,
                                       nWidths = nWidths, dbh_tree_id = dbh_tree_id, 
                                       dbh_year_id = dbh_year_id, incr_tree_id = incr_tree_id,
                                      incr_year_id = incr_year_id, last_time = last_time ))


m$setData(list(logDobs = logDobs, logXobs = logXobs))
m$setInits(list(X = matrix(1, n, nT), D0 = rep(5, n),
                     D = matrix(.1, n, nT), beta0 = .1, beta = rep(.1, n),
                     beta_t = rep(.1, nT), sig_x_obs = .1,
                     sig_d_obs = .1,
                     beta_sd = .1, beta_t_sd = .1, sig_x = .1))

spec <- configureMCMC(m, thin = 50)
#spec$samplerSpecs[[206]]$control$scale=.01
spec$addMonitors(c('D', 'X', 'beta_t', 'beta_t_sd', 'beta_sd'))

Rmcmc <- buildMCMC(spec)
cm <- compileNimble(m)
Cmcmc <- compileNimble(Rmcmc, project = m)

# some timing from NIMBLE 0.3-1
# 505 sec for Rmodel (8.5 mins)
# 300 for configureMCMC
# 140 for buildMCMC
# 270 to compile model
# 325 to compile MCMC

Cmcmc$run(25000)

out <- as.matrix(Cmcmc$mvSamples)[101:500, ]

# t model

m <- nimbleModel(body(ringModelt),
                                   constants = list(n = n, nT = nT, nDBH = nDBH,
                                       nWidths = nWidths, dbh_tree_id = dbh_tree_id,
                                       dbh_year_id = dbh_year_id, incr_tree_id = incr_tree_id,
                                      incr_year_id = incr_year_id, last_time = last_time ))


m$setData(list(logDobs = logDobs, logXobs = logXobs))
m$setInits(list(X = matrix(1, n, nT), D0 = rep(5, n),
                     D = matrix(.1, n, nT), beta0 = .1, beta = rep(.1, n),
                     beta_t = rep(.1, nT), sig_x_obs = .1,
                     sig_d_obs = .1,
                     beta_sd = .1, beta_t_sd = .1, sig_x = .1))


spec <- configureMCMC(m, thin = 50)
#spec$samplerSpecs[[206]]$control$scale=.01
spec$addMonitors(c('D', 'X', 'beta_t', 'beta_t_sd', 'beta_sd'))

Rmcmc <- buildMCMC(spec)
cm <- compileNimble(m)
Cmcmc <- compileNimble(Rmcmc, project = m)

Cmcmc$run(25000)))

out <- as.matrix(Cmcmc$mvSamples)

spec$addSampler('RW_block', list(targetNodes = c('sig_x', 'sig_x_obs'), adaptInterval = 100))
Rmcmc2 <- buildMCMC(spec)
Cmcmc2 <- compileNimble(Rmcmc2, project = mt, resetFunctions = TRUE)

Cmcmc2$run(25000)))
out2 <- as.matrix(Cmcmc2$mvSamples)

spec$removeSamplers(9797)
spec$addSampler('slice', list(targetNode = c('sig_x_obs'), adaptInterval = 100))
spec$addSampler('slice', list(targetNode = c('sig_x'), adaptInterval = 100))
Rmcmc3 <- buildMCMC(spec)
Cmcmc3 <- compileNimble(Rmcmc3, project = m, resetFunctions = TRUE)

Cmcmc3$run(25000)))
out3 <- as.matrix(Cmcmc3$mvSamples)

#save(outt, outt2, outt3, file = '~/research/presentations/umn15/tmodel_output.Rda')
if(FALSE) {
    cor(out[ , c('beta0','beta_sd','beta_t_sd','sig_x','sig_x_obs','sig_d_obs')])
    
    spec$addSampler('RW_block', list(targetNodes = c('sig_x_obs', 'sig_x'), adaptInterval = 100))
    Rmcmc2 <- buildMCMC(spec)
                                        # need to reset the nimbleFunctions in order to add the new MCMC
    Cmcmc2 <- compileNimble(Rmcmc2, project = m, resetFunctions = TRUE)
    
    cm$setData(list(logDobs = logDobs, logXobs = logXobs))
    cm$setInits(list(X = matrix(1, n, nT), D0 = rep(5, n),
                     D = matrix(.1, n, nT), beta0 = .1, beta = rep(.1, n),
                     beta_t = rep(.1, nT), sig_x_obs = .1, sig_d_obs = .1,
                     beta_sd = .1, beta_t_sd = .1, sig_x = .1))
    
    Cmcmc2$run(50000)))
}

# t model with season effects and 1987/1991/1992

m <- nimbleModel(body(ringModeltDate),
                                   constants = list(n = n, nT = nT, nDBH = nDBH,
                                       nWidths = nWidths, dbh_tree_id = dbh_tree_id, dbh_day_id = dbh_day_id,
                                       dbh_year_id = dbh_year_id, incr_tree_id = incr_tree_id,
                                      incr_year_id = incr_year_id, last_time = last_time ))))


m$setData(list(logDobs = logDobs, logXobs = logXobs))
m$setInits(list(X = matrix(1, n, nT), D0 = rep(5, n),
                     D = matrix(.1, n, nT), beta0 = .1, beta = rep(.1, n),
                     beta_t = rep(.1, nT), sig_x_obs = .1,
                     sig_d_obs = .1, b0 = .5, b1 = 10,
                     beta_sd = .1, beta_t_sd = .1, sig_x = .3))
# if start sig_x lower, sig_x_obs can get trapped high

spec <- configureMCMC(m, thin = 50)))
#spec$samplerSpecs[[206]]$control$scale=.01
spec$addMonitors(c('D', 'X', 'beta_t', 'beta_t_sd', 'beta_sd'))

Rmcmc <- buildMCMC(spec)))
cm <- compileNimble(m)))
Cmcmc <- compileNimble(Rmcmc, project = m)))

Cmcmc$run(5000)))

out <- as.matrix(Cmcmc$mvSamples)

# t model with season effects and 1987/1991/1992 and sapling size constraints

m <- nimbleModel(body(ringModeltDateSapl),
                                   constants = list(n = n, nT = nT, nDBH = nDBH, nSaplings = nSaplings,
                                       nWidths = nWidths, dbh_tree_id = dbh_tree_id,
                                       dbh_day_id = dbh_day_id, sapling_tree_id = sapling_tree_id,
                                       sapling_year_id = sapling_year_id, max_size = max_size,
                                       dbh_year_id = dbh_year_id, incr_tree_id = incr_tree_id,
                                      incr_year_id = incr_year_id, last_time = last_time ))))


m$setData(list(logDobs = logDobs, logXobs = logXobs, sapling_constraint = rep(1, nSaplings)))
m$setInits(list(X = matrix(1, n, nT), D0 = rep(0, n),
                     beta0 = .1, beta = rep(.1, n),
                     beta_t = rep(.1, nT), sig_x_obs = .1,
                     sig_d_obs = .1, b0 = .5, b1 = 10,
                     beta_sd = .1, beta_t_sd = .1, sig_x = .3))
# if start sig_x lower, sig_x_obs can get trapped high

spec <- configureMCMC(m, thin = 50)))
#spec$samplerSpecs[[206]]$control$scale=.01
spec$addMonitors(c('D', 'X', 'beta_t', 'beta_t_sd', 'beta_sd'))

Rmcmc <- buildMCMC(spec)))
cm <- compileNimble(m)))
Cmcmc <- compileNimble(Rmcmc, project = ms)))

Cmcmc$run(5000)))

out <- as.matrix(Cmcmc$mvSamples)

### model with taxa/size increment effects

m <- nimbleModel(body(ringModeltDateSaplTaxon),
                                   constants = list(n = n, nT = nT, nDBH = nDBH, nSaplings = nSaplings, nTaxa = nTaxa, taxon = taxon, openDBH = 25,
                                       nWidths = nWidths, dbh_tree_id = dbh_tree_id,
                                       dbh_day_id = dbh_day_id, sapling_tree_id = sapling_tree_id,
                                       sapling_year_id = sapling_year_id, max_size = max_size,
                                       dbh_year_id = dbh_year_id, incr_tree_id = incr_tree_id,
                                      incr_year_id = incr_year_id, last_time = last_time ))))


m$setData(list(logDobs = logDobs, logXobs = logXobs, sapling_constraint = rep(1, nSaplings)))
m$setInits(list(X = matrix(1, n, nT), D0 = rep(0, n),
                     beta0 = .1, beta = rep(.1, n),
                     beta_t = rep(.1, nT), sig_x_obs = .1,
                     sig_d_obs = .1, b0 = .5, b1 = 10,
                     beta_sd = .1, beta_t_sd = .1, sig_x = .3,
                 beta_spp = rep(.1, nTaxa), beta_slope = 0.01, beta_spp_sd = 0.1))
# if start sig_x lower, sig_x_obs can get trapped high

spec <- configureMCMC(m, thin = 50)))
#spec$samplerSpecs[[206]]$control$scale=.01
spec$addMonitors(c('D', 'X', 'beta_t', 'beta_t_sd', 'beta_sd', 'beta_spp', 'beta_slope', 'beta_spp_sd', 'b0', 'b1' ))

Rmcmc <- buildMCMC(spec)))
cm <- compileNimble(m)))
Cmcmc <- compileNimble(Rmcmc, project = m)))

Cmcmc$run(5000)))
#Cmcmc$run(25000)))

out <- as.matrix(Cmcmc$mvSamples)

### model with autocorr increments

m <- nimbleModel(body(ringModeltDateSaplTaxonAR),
                                   constants = list(n = n, nT = nT, nDBH = nDBH, nSaplings = nSaplings,
                                       nTaxa = nTaxa, taxon = taxon, openDBH = 25,
                                       nWidths = nWidths, dbh_tree_id = dbh_tree_id,
                                       dbh_day_id = dbh_day_id, sapling_tree_id = sapling_tree_id,
                                       sapling_year_id = sapling_year_id, max_size = max_size,
                                       dbh_year_id = dbh_year_id, incr_tree_id = incr_tree_id,
                                      incr_year_id = incr_year_id, last_time = last_time ))))


m$setData(list(logDobs = logDobs, logXobs = logXobs, sapling_constraint = rep(1, nSaplings)))
m$setInits(list(logX = log(matrix(1, n, nT)), D0 = rep(0, n),
                     beta0 = .1, beta = rep(.1, n),
                     beta_t = rep(.1, nT), sig_x_obs = .1,
                     sig_d_obs = .1, b0 = .5, b1 = 10,
                     beta_sd = .1, beta_t_sd = .1, sig_x = .3,
                                 beta_sd = .1, beta_t_sd = .1, sig_x = .3,
                 beta_spp = rep(.1, nTaxa), beta_slope = 0.01, beta_spp_sd = 0.1,    
                     nu = 5, rho = 0.9))
# if start sig_x lower, sig_x_obs can get trapped high

spec <- configureMCMC(m, thin = 50)))
#spec$samplerSpecs[[206]]$control$scale=.01
spec$addMonitors(c('D', 'logX', 'beta_t', 'beta_t_sd', 'beta_sd', 'beta_spp', 'beta_slope', 'beta_spp_sd', 'b0', 'b1', 'rho', 'nu' ))

Rmcmc <- buildMCMC(spec)))
cm <- compileNimble(m)))
Cmcmc <- compileNimble(Rmcmc, project = m)))

set.seed(0)
Cmcmc$run(25000)))

out <- as.matrix(Cmcmc$mvSamples)

### model with autocorr increments, testing dinterval

m <- nimbleModel(body(ringModeltDateSaplTaxonARcens),
                                   constants = list(n = n, nT = nT, nDBH = nDBH, nSaplings = nSaplings,
                                       nTaxa = nTaxa, taxon = taxon, openDBH = 25,
                                       nWidths = nWidths, dbh_tree_id = dbh_tree_id,
                                       dbh_day_id = dbh_day_id, sapling_tree_id = sapling_tree_id,
                                       sapling_year_id = sapling_year_id, max_size = max_size,
                                       dbh_year_id = dbh_year_id, incr_tree_id = incr_tree_id,
                                      incr_year_id = incr_year_id, last_time = last_time ))))


m$setData(list(logDobs = logDobs, logXobs = logXobs, not_cens = rep(0, nSaplings)))
m$setInits(list(logX = log(matrix(1, n, nT)), D0 = rep(0, n),
                     beta0 = .1, beta = rep(.1, n),
                     beta_t = rep(.1, nT), sig_x_obs = .1,
                     sig_d_obs = .1, b0 = .5, b1 = 10,
                     beta_sd = .1, beta_t_sd = .1, sig_x = .3,
                                 beta_sd = .1, beta_t_sd = .1, sig_x = .3,
                 beta_spp = rep(.1, nTaxa), beta_slope = 0.01, beta_spp_sd = 0.1,    
                     nu = 5, rho = 0.9))
# if start sig_x lower, sig_x_obs can get trapped high

spec <- configureMCMC(m, thin = 1)))
#spec$samplerSpecs[[206]]$control$scale=.01
spec$addMonitors(c('D', 'logX', 'beta_t', 'beta_t_sd', 'beta_sd', 'beta_spp', 'beta_slope', 'beta_spp_sd', 'b0', 'b1', 'rho', 'nu' ))

Rmcmc <- buildMCMC(spec)))
cm <- compileNimble(m)))
Cmcmc <- compileNimble(Rmcmc, project = m)))

Cmcmc$run(25000)))

out <- as.matrix(Cmcmc$mvSamples)

stop()

### OLD CODE


N <- nrow(data)
nT <- length(ringCols)
logDobs <- log(cbind(matrix(NA, N, nT-2), data$DBH11, data$DBH12))
logXobs <- log(as.matrix(data[ , ringCols]))

code <- body(ringModel) # for Nimble

# run MCMC

require(R2jags, quietly = TRUE) 
out <- jags(data = list(logDobs = logDobs, logXobs = logXobs, N = N, nT = nT),
  parameters.to.save = c('D', 'X', 'beta0', 'beta', 'beta_t',
    'w_sd', 'v_sd', 'tau2_sd', 'vt_sd', 'sigma2_sd'), n.chains = 1,
  n.iter = 5000, n.burnin = 1000, model.file = ringModel, DIC = FALSE)

# MCMC diagnostics
tsplot <- function(x) plot(seq_along(x), x, type = 'l')

out.mcmc <- as.mcmc(out)[[1]]

tsplot(out.mcmc[ , 'beta0'])

beta <- out.mcmc[ , paste('beta[', 1:N, ']', sep = '')]

par(mfrow = c(4, 4))
for(i in 1:N)
  tsplot(beta[ , i])

beta_t <- out.mcmc[ , paste('beta_t[', 1:nT, ']', sep = '')]

par(mfrow = c(5, 6))
for(i in 1:nT)
  tsplot(beta_t[ , i])

par(mfrow = c(2,3))
tsplot(out.mcmc[ , 'w_sd'])
tsplot(out.mcmc[ , 'v_sd'])
tsplot(out.mcmc[ , 'tau2_sd'])
tsplot(out.mcmc[ , 'vt_sd'])
tsplot(out.mcmc[ , 'sigma2_sd'])

id <- 1
Xvals <- out.mcmc[ , paste('X[', id, ',', 1:nT, ']', sep = '')]
plot(1:nT, exp(logXobs[id, ]), type = 'l')
lines(1:nT, colMeans(Xvals), col = 'red')
lines(1:nT, apply(Xvals, 2, quantile, .025), col = 'blue')
lines(1:nT, apply(Xvals, 2, quantile, .975), col = 'blue')


