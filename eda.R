par(mfrow = c(5,5))

tmp = treeMeta
tmp = tmp[order(tmp$dbh), ]
ids = tmp$id

for(i in 1:25) {
    sub= incrData[incrData$id == ids[i],]
    wh = which(treeMeta$id == ids[i])
    d = treeMeta[wh, 'dbh']
    sub = sub[order(sub$year), ]
    sub$incr[sub$incr > 6] = 6
    plot(sub$year, sub$incr, xlim = c(1960, 2012), ylim = c(0,6), type = 'l')
}


merg = merge(incrData, treeMeta)

tmp = merg[merg$year == 2012,]
plot(tmp$dbh, tmp$incr)
tmp = merg[merg$year == 2011,]
plot(tmp$dbh, tmp$incr)
tmp = merg[merg$year == 2010,]
plot(tmp$dbh, tmp$incr)
tmp = merg[merg$year == 2009,]
plot(tmp$dbh, tmp$incr)

par(mfrow=c(2,3))
yr = 2009
mx = max(merg$incr[merg$year == yr])
tmp = merg[merg$year == yr & merg$species == 'ACRU',]
plot(tmp$dbh, tmp$incr, xlim = c(10,70), ylim = c(0,mx))
tmp = merg[merg$year == yr & merg$species == 'QURU',]
plot(tmp$dbh, tmp$incr, xlim = c(10,70), ylim = c(0,mx))
tmp = merg[merg$year == yr & merg$species == 'BEAL',]
plot(tmp$dbh, tmp$incr, xlim = c(10,70), ylim = c(0,mx))
tmp = merg[merg$year == yr & merg$species == 'FAGR',]
plot(tmp$dbh, tmp$incr, xlim = c(10,70), ylim = c(0,mx))
tmp = merg[merg$year == yr & merg$species == 'TSCA',]
plot(tmp$dbh, tmp$incr, xlim = c(10,70), ylim = c(0,mx))

dtmp = census
dtmp$dbh69[dtmp$cond69 != "L"] = NA
dtmp$dbh75[dtmp$cond75 != "L"] = NA
dtmp$dbh91[dtmp$cond91 != "L"] = NA
dtmp$dbh01[dtmp$cond01 != "L"] = NA
dtmp$dbh11[dtmp$cond11 != "L"] = NA
dtmp <- dtmp[,c('census_id', 'dbh69','dbh75','dbh91','dbh01','dbh11')]
dtmp = merge(dtmp, treeMeta[ , c('census_id', 'dbh')], by.x = 'census_id', by.y = 'census_id', all.x= TRUE, all.y = FALSE)

par(mfrow=c(5,5), mai = c(.4,.4,.3,.1))
for(i in 1:25) {
 #   mx = max(dtmp[i, 2:7], na.rm =TRUE)
    print(i)
    if(sum(!is.na(dtmp[i,2:7]))) {
        plot(c(1969,1975,1991,2001,2011,2012), dtmp[i, 2:7], xlab = 'yr', ylab='',
             xlim = c(1968,2013))
        title(dtmp$census_id[i])
    } else plot(0,0)
}

getID <- function(census_id, neil_id, useCensus = TRUE) {
    if(useCensus) {
        tmp <- census[census$census_id == census_id, , drop = FALSE]
    } else {
        tmp <- census[census$id == neil_id, , drop = FALSE]
    }
    return(tmp[['stat_id']])
}


out = outs
pm <- colMeans(out)

makePlot <-  function(focal, pm, logX = FALSE) {
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

censusYrs <- as.numeric(strsplit(censusYears, ',')[[1]])

logX = TRUE
focal=0
focal = focal+1
makePlot(focal, pm, logX)


cor(outt[,'sig_x_obs'], outt)


# look at beta[i] and beta_t[j] values

# look at whether empirical dbh11-dbh01 suggests dbh was measured at begin or end of period

# look at sig_d_obs - how big for different sizes of d?
# pretty big

# plot growth from full census data as fxn of taxa and size
censusFull <- read.csv(censusFile)
names(censusFull)[names(censusFull) == "treeid"] <- "census_id"

# 1969
cf = censusFull
cf <- cf[!is.na(cf$cond01) & !is.na(cf$cond11) & !is.na(cf$dbh11) & !is.na(cf$dbh01)
        & cf$cond01 == "L" & cf$cond11 == "L", ]

par(mfrow = c(3, 3))
plot(cf$dbh01, 10*(cf$dbh11 - cf$dbh01) / 10, main = 'all')
for(sp in c('ACRU','BEAL','BELE','FAGR','HAVI','PIST','QURU','TSCA')) 
plot(cf$dbh01[cf$species == sp], 10*(cf$dbh11[cf$species == sp] - cf$dbh01[cf$species == sp]) / 10, main = sp, xlim = c(5,75), ylim = c(-5,6.5))

# log scale
par(mfrow = c(3, 3))
plot(cf$dbh01, log(10*(cf$dbh11 - cf$dbh01) / 10), main = 'all')
for(sp in c('ACRU','BEAL','BELE','FAGR','HAVI','PIST','QURU','TSCA')) 
plot(cf$dbh01[cf$species == sp], log(10*(cf$dbh11[cf$species == sp] - cf$dbh01[cf$species == sp]) / 10), main = sp, xlim = c(5,75), ylim = c(-2,log(6.5)))


# incr from ring data as fxn of taxa and size
par(mfrow = c(3, 3))
tmp <- merge(incrData, treeMeta[,c('id','species','dbh')], all.x = TRUE, all.y = FALSE)
plot(tmp$dbh[tmp$year > 2007], tmp$incr[tmp$year > 2007], main = 'all')
for(sp in c('ACRU','BEAL','BEAL','FAGR','HAVI','PIST','QURU','TSCA')) 
plot(tmp$dbh[tmp$species == sp &tmp$year > 2007], tmp$incr[tmp$species==sp&tmp$year > 2007], main = sp, xlim = c(5,75), ylim = c(-5,6.5))
#plot(tmp$dbh[tmp$species == sp &tmp$year > 2007], log(tmp$incr[tmp$species==sp&tmp$year > 2007]), main = sp, xlim = c(5,75))



tsplot <- function(x) plot(1:length(x),x,type='l')

par(mfrow=c(3,4))
tsplot(tmp[, 'beta0'])
tsplot(tmp[, 'sig_x'])
tsplot(tmp[, 'sig_x_obs'])
tsplot(tmp[, 'sig_d_obs'])
tsplot(tmp[, 'D0[5]'])
tsplot(tmp[, 'D0[10]'])
tsplot(tmp[, 'D0[15]'])
set.seed(0); ids <- sample(1:n, 5, replace = T); tt <- sample(1:nT, 5, replace = T)
for(i in 1:5)
    tsplot(tmp[,paste0("X[", ids[i], ", ", tt[i], "]")])

# get census from process.R but without screening for aliveness (thrugh line 62)

cfBig <- census[(!is.na(census$dbh69) & census$dbh69 > 20) | (!is.na(census$dbh75) & census$dbh75 > 20) | (!is.na(census$dbh91) & census$dbh91 > 20) | (!is.na(census$dbh01) & census$dbh01 > 20) | (!is.na(census$dbh11) & census$dbh11 > 20), ]
tmp <- merge(cfBig, treeMeta[!is.na(treeMeta$census_id),], by.x = 'census_id', by.y = 'census_id', all.x = TRUE, all.y = FALSE) 


# plot increments by size

load('/var/tmp/paciorek/paleon/npp/data/lyford/data.Rda')

dat = data.frame(incr)

dat <- cbind(dimnames(dat)[[1]], rep(0, nrow(dat)), dat)

for(i in 1:nrow(dat))
    dat[i,2] <- treeMeta$dbh[treeMeta$id == dat[i,1]]

names(dat)[1:2] = c('id','dbh')
dat = dat[order(dat$dbh), ]

tmp1 = dat[dat$dbh < 20, ]
tmp2 = dat[dat$dbh >= 20 & dat$dbh < 30,]
tmp3 = dat[dat$dbh >= 30 & dat$dbh < 45,]
tmp4 = dat[dat$dbh >= 45,]

par(mfrow=c(2,2))

for(i in 1:4) {
    dat <- get(paste0("tmp",i))
    plot(1960:2012, dat[1, 3:ncol(dat)], type = 'l', ylim = c(0, 5))
    for(j in 2:5) # nrow(dat))
        lines(1960:2012, dat[j, 3:ncol(dat)], col = j)
}


par(mfrow=c(4,4))
start = 1
end = 5
for(i in 1:16) {
    plot(1960:2012, dat[start, 3:ncol(dat)], type = 'l', ylim = c(0, 5))
    for(j in (start+1):end) # nrow(dat))
        lines(1960:2012, dat[j, 3:ncol(dat)], col = j)
    start = start + 5
    end = end + 5
}

start = 0
par(mfrow=c(1,1))
start = start + 1; plot(1960:2012, log(dat[start, 3:ncol(dat)]), type = 'l', ylim = c(-3,log(5))); points(1960:2012, log(dat[start,3:ncol(dat)]))
