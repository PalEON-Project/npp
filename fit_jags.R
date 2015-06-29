source("config")

source("model_jags.R")

setwd(dataDir)
setwd('lyford')
load('data.Rda')

library(R2jags)

setwd(projectDir)

Xinit <- matrix(1, n, nT)
for(i in 1:n) {
    if(last_time[i] < nT) 
        Xinit[i,last_time[i]:nT] <- NA
}
system.time({
out <- jags(data = list(n = n, nT = nT, nDBH = nDBH, nWidths = nWidths,
                         dbh_tree_id = dbh_tree_id, dbh_year_id = dbh_year_id, incr_tree_id = incr_tree_id,
                                      incr_year_id = incr_year_id, last_time = last_time,logDobs = logDobs, logXobs = logXobs), inits = list(list(X = Xinit, D0 = rep(5, n), beta0 = .1, beta = rep(.1, n),
                     beta_t = rep(.1, nT), sig_x_obs = .1, sig_d_obs = .1,
                     beta_sd = .1, beta_t_sd = .1, sig_x = .1)), parameters.to.save = c("X", "D", "D0", "sig_d_obs", "sig_x_obs", "sig_x","beta0","beta_t"), n.chains = 1, n.iter = 50000, n.thin = 50, DIC = FALSE, model = "model_jags.bug")})

# init model 10:07:20 - 10:35:40: 28 minutes
# fit 10:35:40-11:07: 33 minutes
# postprocess: ?
# total: 61 minutes


# try with taxon model that is slow in NIMBLE

source("config")

setwd(dataDir)
setwd('lyford')
load('data.Rda')

library(R2jags)

setwd(projectDir)
Xinit <- matrix(1, n, nT)
for(i in 1:n) {
    if(last_time[i] < nT) 
        Xinit[i,last_time[i]:nT] <- NA
}

inits = list(X = Xinit, D0 = rep(0, n),
                     beta0 = .1, beta = rep(.1, n),
                     beta_t = rep(.1, nT), sig_x_obs = .1,
                     sig_d_obs = .1, b0 = .5, b1 = 10,
                     beta_sd = .1, beta_t_sd = .1, sig_x = .3,
                 beta_spp = rep(.1, nTaxa), beta_slope = 0.01, beta_spp_sd = 0.1)
system.time({
out <- jags(data = list(n = n, nT = nT, nDBH = nDBH, nWidths = nWidths,nTaxa = nTaxa, taxon = taxon, openDBH = 25,
                         dbh_day_id = dbh_day_id,dbh_tree_id = dbh_tree_id, dbh_year_id = dbh_year_id, incr_tree_id = incr_tree_id,
                                      incr_year_id = incr_year_id, last_time = last_time,logDobs = logDobs, logXobs = logXobs), inits = list(inits), parameters.to.save = c('D', 'X', 'beta_t', 'beta_t_sd', 'beta_sd', 'beta_spp', 'beta_slope', 'beta_spp_sd', 'b0', 'b1'), n.chains = 1, n.iter = 50, n.thin = 1, DIC = FALSE, model = "model_jags2.bug")
})

