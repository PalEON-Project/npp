# note this functions uses various objects from data.Rda as global variables
treePlot <-  function(focal, pm, logX = FALSE) {
    focal_neil <- census$id[census$stat_id == focal]
    focal_audrey <- census$census_id[census$stat_id == focal]
    
    xdat <- exp(logXobs[incr_tree_id == focal])
    xdatt <- incr_year_id[incr_tree_id == focal]
    
    ddat <- exp(logDobs[dbh_tree_id == focal])
    ddatt <- dbh_year_id[dbh_tree_id == focal]
    
    nm <- dimnames(out)[[2]]
    dinds <- grep(paste0("D[", focal, ", "), nm, fixed = TRUE)
    if(logX) {
            xinds <- grep(paste0("logX[", focal, ", "), nm, fixed = TRUE)
        } else  xinds <- grep(paste0("X[", focal, ", "), nm, fixed = TRUE)
    
    yrs <- 1960:2012
    
    mx = max(c(ddat, pm[dinds]), na.rm = TRUE)
    par(mfrow = c(2,1))
    plot(yrs, rep(0, nT), type = 'n', ylim = c(0, mx+5))
    points(yrs[ddatt], ddat, col = 'red')
    lines(yrs, pm[dinds])
    title(c(focal_neil, focal_audrey))
    if(!is.na(focal_neil))
        points(2012, treeMeta$dbh[treeMeta$id == focal_neil], pch = 'x', col = 'red')
    points(censusYrs, rep(5, length(censusYrs)), pch='x')
    
    plot(yrs, rep(0, nT), type = 'n', ylim = c(0, 10))
    if(!is.na(focal_neil))
        points(yrs[xdatt], xdat, col = 'red')
    if(logX)    lines(yrs, exp(pm[xinds])) else   lines(yrs, pm[xinds])

    if(!is.na(focal_neil)) {
        print(c(census$dbh11[focal], treeMeta$dbh[treeMeta$id == focal_neil]))
        print(c(xdat[xdatt==51], xdat[xdatt==52], xdat[xdatt==53]))
    }
}


# example usage: loop over the 198 trees by their 'stat_id' (see 'census' object)
if(FALSE) {
    burnin <- 100
    samples <- as.matrix(Cmcmc$mvSamples)
    pm <- colMeans(samples[(burnin+1):nrow(samples)])
    pdf('treePlots.pdf')
    for(focal in 1:nrow(census)) {
        treePlot(focal, pm, logX = TRUE)
    }
    dev.off()
}
