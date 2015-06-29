ringModel <- function() {

    # dbh likelihood
    for(j in 1:nDBH) {
        logDobs[j] ~ dnorm(log(D[dbh_tree_id[j], dbh_year_id[j]]), sd = sig_d_obs)
    }

    # increment likelihood for trees with rings
    for(i in 1:nWidths) {
        logXobs[i] ~ dnorm(log(X[incr_tree_id[i], incr_year_id[i]]), sd = sig_x_obs)
    }

    # process evolution
    for(i in 1:n) {
        D0[i] ~ dunif(-200, 300)
        D[i, 1] <- D0[i] + X[i, 1]/10
        
        for(j in 1:last_time[i]) {
            X[i, j] ~ dlnorm(beta[i] + beta_t[j], sd = sig_x)
        }
        for(j in (last_time[i]+1):nT) {
            # BUGS (and NIMBLE by default) should ignore if
            # indexing is decreasing
            X[i, j] <- 0
        }
        
        for(j in 2:nT) {
            D[i, j] <- D[i, j-1] + X[i, j]/10
        }
        beta[i] ~ dnorm(beta0, sd = beta_sd)
         
    }

    for(j in 1:nT) {
        beta_t[j] ~ dnorm(0, sd = beta_t_sd)
    }

    
    beta0 ~ dnorm(0, .00001)
    sig_d_obs ~ dunif(0, 1000)
    sig_x_obs ~ dunif(0, 1000)

    sig_x ~ dunif(0, 1000)
    beta_sd ~ dunif(0, 1000)
    beta_t_sd ~ dunif(0, 1000)
}


ringModelt <- function() {

    # dbh likelihood
    for(j in 1:nDBH) {
        logDobs[j] ~ dt(log(D[dbh_tree_id[j], dbh_year_id[j]]), tau_d_obs, 3)
    }

    # increment likelihood for trees with rings
    for(i in 1:nWidths) {
        logXobs[i] ~ dnorm(log(X[incr_tree_id[i], incr_year_id[i]]), sd = sig_x_obs)
    }

    # process evolution
    for(i in 1:n) {
        D0[i] ~ dunif(-200, 300)
        D[i, 1] <- D0[i] + X[i, 1]/10
        
        for(j in 1:last_time[i]) {
            X[i, j] ~ dlnorm(beta[i] + beta_t[j], sd = sig_x)
        }
        for(j in (last_time[i]+1):nT) {
            # BUGS (and NIMBLE by default) should ignore if
            # indexing is decreasing
            X[i, j] <- 0
        }
        
        for(j in 2:nT) {
            D[i, j] <- D[i, j-1] + X[i, j]/10
        }
        beta[i] ~ dnorm(beta0, sd = beta_sd)
         
    }

    for(j in 1:nT) {
        beta_t[j] ~ dnorm(0, sd = beta_t_sd)
    }
    tau_d_obs <- 1/(sig_d_obs * sig_d_obs)
    
    beta0 ~ dnorm(0, .00001)
    sig_d_obs ~ dunif(0, 1000)
    sig_x_obs ~ dunif(0, 1000)

    sig_x ~ dunif(0, 1000)
    beta_sd ~ dunif(0, 1000)
    beta_t_sd ~ dunif(0, 1000)
}


ringModeltDate <- function() {

    b0 ~ dnorm(0, 1)
    b1 ~ dunif(5, 25)
        
    
    # dbh likelihood
    for(j in 1:nDBH) {
        # -1 is for begin of grow season
        logDobs[j] ~ dt(log(D[dbh_tree_id[j], dbh_year_id[j]-1] + 
                                expit(b0 + b1*(dbh_day_id[j] - 0.5)) *
                                    X[dbh_tree_id[j], dbh_year_id[j]]/10),
                            tau_d_obs, 3)
    }

    # increment likelihood for trees with rings
    for(i in 1:nWidths) {
        logXobs[i] ~ dnorm(log(X[incr_tree_id[i], incr_year_id[i]]), sd = sig_x_obs)
    }

    # process evolution
    for(i in 1:n) {
        D0[i] ~ dunif(-200, 300)
        D[i, 1] <- D0[i] + X[i, 1]/10
        
        for(j in 1:last_time[i]) {
            X[i, j] ~ dlnorm(beta[i] + beta_t[j], sd = sig_x)
        }
        for(j in (last_time[i]+1):nT) {
            # BUGS (and NIMBLE by default) should ignore if
            # indexing is decreasing
            X[i, j] <- 0
        }
        
        for(j in 2:nT) {
            D[i, j] <- D[i, j-1] + X[i, j]/10
        }
        beta[i] ~ dnorm(beta0, sd = beta_sd)
         
    }

    for(j in 1:nT) {
        beta_t[j] ~ dnorm(0, sd = beta_t_sd)
    }
    tau_d_obs <- 1/(sig_d_obs * sig_d_obs)
    
    beta0 ~ dnorm(0, .00001)
    sig_d_obs ~ dunif(0, 1000)
    sig_x_obs ~ dunif(0, 1000)

    sig_x ~ dunif(0, 1000)
    beta_sd ~ dunif(0, 1000)
    beta_t_sd ~ dunif(0, 1000)
}

ringModeltDateSapl <- function() {

    b0 ~ dnorm(0, 1)
    b1 ~ dunif(5, 25)
        
    
    # dbh likelihood
    for(j in 1:nDBH) {
        # -1 is for begin of grow season
        logDobs[j] ~ dt(log(D[dbh_tree_id[j], dbh_year_id[j]-1] + 
                                expit(b0 + b1*(dbh_day_id[j] - 0.5)) *
                                    X[dbh_tree_id[j], dbh_year_id[j]]/10),
                            tau_d_obs, 3)
    }

    # increment likelihood for trees with rings
    for(i in 1:nWidths) {
        logXobs[i] ~ dnorm(log(X[incr_tree_id[i], incr_year_id[i]]), sd = sig_x_obs)
    }

    # constraint on trees not yet appearing in a census
    for(i in 1:nSaplings) {
        sapling_ind[i] <- D[sapling_tree_id[i], sapling_year_id[i]] < max_size[i]
        sapling_constraint[i] ~ dconstraint(sapling_ind[i])
    }

    # process evolution
    for(i in 1:n) {
        D0[i] ~ dunif(-200, 300)
        D[i, 1] <- D0[i] + X[i, 1]/10
        
        for(j in 1:last_time[i]) {
            X[i, j] ~ dlnorm(beta[i] + beta_t[j], sd = sig_x)
        }
        for(j in (last_time[i]+1):nT) {
            # BUGS (and NIMBLE by default) should ignore if
            # indexing is decreasing
            X[i, j] <- 0
        }
        
        for(j in 2:nT) {
            D[i, j] <- D[i, j-1] + X[i, j]/10
        }
        beta[i] ~ dnorm(beta0, sd = beta_sd)
         
    }

    for(j in 1:nT) {
        beta_t[j] ~ dnorm(0, sd = beta_t_sd)
    }
    tau_d_obs <- 1/(sig_d_obs * sig_d_obs)
    
    beta0 ~ dnorm(0, .00001)
    sig_d_obs ~ dunif(0, 1000)
    sig_x_obs ~ dunif(0, 1000)

    sig_x ~ dunif(0, 1000)
    beta_sd ~ dunif(0, 1000)
    beta_t_sd ~ dunif(0, 1000)
}

ringModeltDateSaplTaxon <- function() {

    b0 ~ dnorm(0, 1)
    b1 ~ dunif(5, 25)
        
    
    # dbh likelihood
    for(j in 1:nDBH) {
        # -1 is for begin of grow season
        logDobs[j] ~ dt(log(D[dbh_tree_id[j], dbh_year_id[j]-1] + 
                                expit(b0 + b1*(dbh_day_id[j] - 0.5)) *
                                    X[dbh_tree_id[j], dbh_year_id[j]]/10),
                            tau_d_obs, 3)
    }

    # increment likelihood for trees with rings
    for(i in 1:nWidths) {
        logXobs[i] ~ dnorm(log(X[incr_tree_id[i], incr_year_id[i]]), sd = sig_x_obs)
    }

    # constraint on trees not yet appearing in a census
    for(i in 1:nSaplings) {
        sapling_ind[i] <- D[sapling_tree_id[i], sapling_year_id[i]] < max_size[i]
        sapling_constraint[i] ~ dconstraint(sapling_ind[i])
    }

    # process evolution
    for(i in 1:n) {
        D0[i] ~ dunif(-200, 300)
        D[i, 1] <- D0[i] + X[i, 1]/10

        mnX[i, 1] <- beta[i] + beta_slope * (D0[i] - openDBH) * (D0[i] > openDBH) + beta_t[1]
        for(j in 2:last_time[i]) {
            mnX[i, j] <- beta[i] + beta_slope * (D[i, j-1] - openDBH) * (D[i, j-1] > openDBH) + beta_t[j]
        }
        for(j in 1:last_time[i]) {
            X[i, j] ~ dlnorm(mnX[i, j], sd = sig_x)
        }
        for(j in (last_time[i]+1):nT) {
            # BUGS (and NIMBLE by default) should ignore if
            # indexing is decreasing
            X[i, j] <- 0
        }
        
        for(j in 2:nT) {
            D[i, j] <- D[i, j-1] + X[i, j]/10
        }
        beta[i] ~ dnorm(beta_spp[taxon[i]], sd = beta_sd)
         
    }
    for(j in 1:nTaxa) {
        beta_spp[j] ~ dnorm(beta0, sd = beta_spp_sd)
    }

    for(j in 1:nT) {
        beta_t[j] ~ dnorm(0, sd = beta_t_sd)
    }
    tau_d_obs <- 1/(sig_d_obs * sig_d_obs)
    
    beta0 ~ dnorm(0, .00001)
    beta_slope ~ dunif(0, .3)
    sig_d_obs ~ dunif(0, 1000)
    sig_x_obs ~ dunif(0, 1000)

    sig_x ~ dunif(0, 1000)
    beta_sd ~ dunif(0, 1000)
    beta_t_sd ~ dunif(0, 1000)
    beta_spp_sd ~ dunif(0, 1000)
}


ringModeltDateSaplTaxonAR <- function() {

    b0 ~ dnorm(0, 1)
    b1 ~ dunif(5, 25)
        
    
    # dbh likelihood
    for(j in 1:nDBH) {
        # -1 is for begin of grow season
        logDobs[j] ~ dt(log(D[dbh_tree_id[j], dbh_year_id[j]-1] + 
                                expit(b0 + b1*(dbh_day_id[j] - 0.5)) *
                                    exp(logX[dbh_tree_id[j], dbh_year_id[j]])/10),
                            sigma = sig_d_obs, df = 3)
    }

    # increment likelihood for trees with rings
    for(i in 1:nWidths) {
        logXobs[i] ~ dnorm(logX[incr_tree_id[i], incr_year_id[i]], sd = sig_x_obs)
    }

    # constraint on trees not yet appearing in a census
    for(i in 1:nSaplings) {
        sapling_ind[i] <- D[sapling_tree_id[i], sapling_year_id[i]] < max_size[i]
        sapling_constraint[i] ~ dconstraint(sapling_ind[i])
    }

    # process evolution
    for(i in 1:n) {
        D0[i] ~ dunif(-200, 300)
        D[i, 1] <- D0[i] + exp(logX[i, 1])/10

        mnX[i, 1] <- beta[i] + beta_slope * (D0[i] - openDBH) * (D0[i] > openDBH) + beta_t[1]
        for(j in 2:last_time[i]) {
            mnX[i, j] <- beta[i] + beta_slope * (D[i, j-1] - openDBH) * (D[i, j-1] > openDBH) + beta_t[j]
        }
        # 2.23 is largest observed log increment
        logX[i, 1] ~ dt(mnX[i, 1], sigma = sig_x, df = nu)
        for(j in 2:last_time[i]) {
            # solely to avoid the 60 char issue in deparse
            # mnlogX[i, j] <- mnX[i, j] + rho*(logX[i, j-1] - mnX[i, j])
            # logX[i, j] ~ dt(mnlogX[i,j], tau_x, df = nu)
            logX[i, j] ~ dt(mnX[i, j] + rho*(logX[i, j-1] - mnX[i, j]), sigma = sig_x, df = nu)
        }
        for(j in (last_time[i]+1):nT) {
            # BUGS (and NIMBLE by default) should ignore if
            # indexing is decreasing
            logX[i, j] <- -100
        }
        
        for(j in 2:nT) {
            D[i, j] <- D[i, j-1] + exp(logX[i, j])/10
        }
        beta[i] ~ dnorm(beta_spp[taxon[i]], sd = beta_sd)
         
    }
    for(j in 1:nTaxa) {
        beta_spp[j] ~ dnorm(beta0, sd = beta_spp_sd)
    }

    for(j in 1:nT) {
        beta_t[j] ~ dnorm(0, sd = beta_t_sd)
    }
   # tau_d_obs <- 1/(sig_d_obs * sig_d_obs)
   # tau_x <- 1/(sig_x * sig_x)
    
    rho ~ dunif(0, 1)
    nu ~ dunif(1, 20)
    
    beta0 ~ dnorm(0, .00001)
    beta_slope ~ dunif(0, 0.3)
    sig_d_obs ~ dunif(0, 1000)
    sig_x_obs ~ dunif(0, 1000)

    sig_x ~ dunif(0, 1000)
    beta_sd ~ dunif(0, 1000)
    beta_t_sd ~ dunif(0, 1000)
    beta_spp_sd ~ dunif(0, 1000)
}

# using dinterval not dconstraint -- seems to work fine
ringModeltDateSaplTaxonARcens <- function() {

    b0 ~ dnorm(0, 1)
    b1 ~ dunif(5, 25)
        
    
    # dbh likelihood
    for(j in 1:nDBH) {
        # -1 is for begin of grow season
        logDobs[j] ~ dt(log(D[dbh_tree_id[j], dbh_year_id[j]-1] + 
                                expit(b0 + b1*(dbh_day_id[j] - 0.5)) *
                                    exp(logX[dbh_tree_id[j], dbh_year_id[j]])/10),
                            tau_d_obs, 3)
    }

    # increment likelihood for trees with rings
    for(i in 1:nWidths) {
        logXobs[i] ~ dnorm(logX[incr_tree_id[i], incr_year_id[i]], sd = sig_x_obs)
    }

    # constraint on trees not yet appearing in a census
    for(i in 1:nSaplings) {
        not_cens[i] ~ dinterval(D[sapling_tree_id[i], sapling_year_id[i]], max_size[i])
    }

    # process evolution
    for(i in 1:n) {
        D0[i] ~ dunif(-200, 300)
        D[i, 1] <- D0[i] + exp(logX[i, 1])/10

        mnX[i, 1] <- beta[i] + beta_slope * (D0[i] - openDBH) * (D0[i] > openDBH) + beta_t[1]
        for(j in 2:last_time[i]) {
            mnX[i, j] <- beta[i] + beta_slope * (D[i, j-1] - openDBH) * (D[i, j-1] > openDBH) + beta_t[j]
        }
        logX[i, 1] ~ dt(mnX[i, 1], tau_x, df = nu)
        for(j in 2:last_time[i]) {
            # solely to avoid the 60 char issue in deparse
            mnlogX[i, j] <- mnX[i, j] + rho*(logX[i, j-1] - mnX[i, j])
            logX[i, j] ~ dt(mnlogX[i, j], tau_x, df = nu)
        }
        for(j in (last_time[i]+1):nT) {
            # BUGS (and NIMBLE by default) should ignore if
            # indexing is decreasing
            logX[i, j] <- -100
        }
        
        for(j in 2:nT) {
            D[i, j] <- D[i, j-1] + exp(logX[i, j])/10
        }
        beta[i] ~ dnorm(beta_spp[taxon[i]], sd = beta_sd)
         
    }
    for(j in 1:nTaxa) {
        beta_spp[j] ~ dnorm(beta0, sd = beta_spp_sd)
    }

    for(j in 1:nT) {
        beta_t[j] ~ dnorm(0, sd = beta_t_sd)
    }
    tau_d_obs <- 1/(sig_d_obs * sig_d_obs)
    tau_x <- 1/(sig_x * sig_x)
    
    rho ~ dunif(0, 1)
    nu ~ dunif(1, 20)
    
    beta0 ~ dnorm(0, .00001)
    beta_slope ~ dunif(0, 0.3)
    sig_d_obs ~ dunif(0, 1000)
    sig_x_obs ~ dunif(0, 1000)

    sig_x ~ dunif(0, 1000)
    beta_sd ~ dunif(0, 1000)
    beta_t_sd ~ dunif(0, 1000)
    beta_spp_sd ~ dunif(0, 1000)
}
