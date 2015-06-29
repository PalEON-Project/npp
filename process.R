#!/usr/bin/Rscript

source("config")

setwd(projectDir)
centroids <- read.csv('centroids.csv')
nPlots <- nrow(centroids)

ftPerMeter <- 3.2808399

lastYear <- 2012
firstYear <- 1960
years <- firstYear:lastYear
nT <- length(years)

growDays <- 30+31+30+31+31+30+31

library(dplR)
library(fields)
library(reshape2)

setwd(dataDir)
setwd('lyford')

# encapsulate firstYear in getTimeIndex closure
getTimeIndex_gen <- function(firstYear) {
    function(year) {
        year - firstYear + 1
    }
}
getTimeIndex <- getTimeIndex_gen(firstYear)


## read data in

censusFull <- read.csv(censusFile, stringsAsFactors = FALSE)
ringMeta <- read.csv(ringMetaFile, stringsAsFactors = FALSE)
treeMeta <- read.csv(treeMetaFile, skip = 2, stringsAsFactors = FALSE)
censusDates <- read.csv(censusSampleDatesFile, stringsAsFactors = FALSE)
    
# remove second tree that has Audrey id of 1863; this one is not consistent with Audrey's data
treeMeta <- treeMeta[!(treeMeta$Site == "LF2" & treeMeta$Tree.Number == 27), ]

rwFiles <- list.files("RW")
rwFiles <- rwFiles[grep(".rw$", rwFiles)]
rwData <- list()
for(fn in rwFiles) {
    id <- gsub(".rw", "", fn)
    rwData[[id]] <- t(read.tucson(file.path("RW", fn)))  # rows are tree, cols are times
}

## initial processing of input datasets
names(censusFull)[names(censusFull) == "treeid"] <- "census_id"

ringMeta$SITE <- as.numeric(gsub("LF", "", ringMeta$SITE))
ringMeta$X <-  NULL
names(ringMeta) <- tolower(names(ringMeta))

names(treeMeta) <- tolower(names(treeMeta))
treeMeta$site <- as.numeric(gsub("LF", "", treeMeta$site))
wh <- which(names(treeMeta) == "tag")
names(treeMeta)[wh] <- "census_id" # to match census file
treeMeta$id <- treeMeta$site*100 + treeMeta$tree.number

## initial processing of census sample dates
threshold <- function(vals, minDate, lower = TRUE) {
    if(!is(minDate, "Date")) minDate <- as.Date(minDate)
    if(length(minDate) == 1) minDate <- rep(minDate, length(vals))
    if(lower) {
        vals[julian(vals) < julian(minDate)] <- minDate[julian(vals) < julian(minDate)]
    } else vals[julian(vals) > julian(minDate)] <- minDate[julian(vals) > julian(minDate)]
    return(vals)
}

inds <- which(censusDates$X1987.1992_start == "")
censusDates$X1987.1992_start[inds] <-  censusDates$X1987.1992_end[inds]
censusDates$yr1991 <- as.numeric(gsub(".*(19[89][0-9])$", "\\1", censusDates$X1987.1992_start))

start <- as.Date("1969-3-31")
end <- as.Date("1969-11-1")
tmp1 <- as.Date(censusDates$X1969_start, format = "%m/%d/%Y")
tmp1 <- threshold(tmp1, start)
tmp1 <- threshold(tmp1, end, lower = FALSE)
tmp2 <- as.Date(censusDates$X1969_end, format = "%m/%d/%Y")
tmp2 <- threshold(tmp2, start)
tmp2 <- threshold(tmp2, end, lower = FALSE)
censusDates$day69 <- (julian(tmp1) + julian(tmp2)) / 2 - julian(start)

start <- as.Date("1975-3-31")
end <- as.Date("1975-11-1")
censusDates$X1975[is.na(censusDates$X1975)] <- "5/15/1975" # rough average of other blocks
tmp1 <- as.Date(censusDates$X1975, format = "%m/%d/%Y")
tmp1 <- threshold(tmp1, start)
tmp1 <- threshold(tmp1, end, lower = FALSE)
censusDates$day75 <- julian(tmp1) - julian(start)

start <- as.Date(paste0(censusDates$yr1991, "-3-31"))
end <- as.Date(paste0(censusDates$yr1991, "-11-1"))
tmp1 <- as.Date(censusDates$X1987.1992_start, format = "%m/%d/%Y")
tmp1 <- threshold(tmp1, start)
tmp1 <- threshold(tmp1, end, lower = FALSE)
tmp2 <- as.Date(censusDates$X1987.1992_end, format = "%m/%d/%Y")
tmp2 <- threshold(tmp2, start)
tmp2 <- threshold(tmp2, end, lower = FALSE)
censusDates$day91 <- (julian(tmp1) + julian(tmp2)) / 2 - julian(start)

start <- as.Date("2001-3-31")
end <- as.Date("2001-11-1")
tmp1 <- as.Date(censusDates$X2001_start, format = "%m/%d/%Y")
tmp1 <- threshold(tmp1, start)
tmp1 <- threshold(tmp1, end, lower = FALSE)
tmp2 <- as.Date(censusDates$X2001_end, format = "%m/%d/%Y")
tmp2 <- threshold(tmp2, start)
tmp2 <- threshold(tmp2, end, lower = FALSE)
censusDates$day01 <- (julian(tmp1) + julian(tmp2)) / 2 - julian(start)

start <- as.Date("2011-3-31")
end <- as.Date("2011-11-1")
tmp1 <- as.Date(censusDates$X2011_start, format = "%m/%d/%Y")
tmp1 <- threshold(tmp1, start)
tmp1 <- threshold(tmp1, end, lower = FALSE)
tmp2 <- as.Date(censusDates$X2011_end, format = "%m/%d/%Y")
tmp2 <- threshold(tmp2, start)
tmp2 <- threshold(tmp2, end, lower = FALSE)
censusDates$day11 <- (julian(tmp1) + julian(tmp2)) / 2 - julian(start)


## get census trees only in plots
census <- list(); length(census) <- nPlots
for(i in seq_len(nPlots)) {
    dist <- rdist(centroids[i, 2:3], censusFull[ , c('xsite', 'ysite')])
    census[[i]] <- censusFull[dist < plotRadius*ftPerMeter, ]
    census[[i]]$site <- i
    census[[i]]$dist_census <- dist[dist < plotRadius*ftPerMeter]
}

census <- do.call(rbind, census)
dbhCols <- grep("dbh", names(census))
condCols <- grep("cond", names(census))

existFun <- function(x) {
    dbhs <- x[ , dbhCols]
    conds <- x[ , condCols]
    numAlive <- sum(conds == "L", na.rm = TRUE)
    numAlive[is.na(numAlive)] <- 0
    numAlive[]
    return(numAlive > 0)
}

numAlive <- apply(census, 1, function(x) sum(x[condCols] == "L", na.rm = TRUE))
noDbh <- apply(census[ , dbhCols], 1, function(x) sum(is.na(x)) == length(dbhCols))
numAlive[noDbh] <- 0
census <- census[numAlive > 0, ]
# for now throw out trees never alive in a census - assume these were never big so limited impact on biomass increment

censusYears <- strsplit(censusYears, ",")[[1]]
censusYears2digit <- substring(censusYears, 3, 4)
nCensus <- length(censusYears)


census <- merge(census, treeMeta, all.x = TRUE, all.y = FALSE)

census <- merge(census, censusDates[ , c('Block', 'yr1991', paste0('day', censusYears2digit))], by.x = 'block', by.y = 'Block', all.x = TRUE, all.y = FALSE)

census <- census[order(census$site, census$id), ]
census$stat_id <- 1:nrow(census)


## extract dbh values for living trees

dbhCols <- paste0('dbh', censusYears2digit)
condCols <- paste0('cond', censusYears2digit)
dayCols <- paste0('day', censusYears2digit)

cols <- c("census_id", "id", "stat_id", "species", dbhCols, condCols)
tmp <- census[ , cols]
dbh <- melt(tmp, id.vars = c("species","census_id", "stat_id", "id", condCols))
dbh$yr=substring(dbh$variable,4,5)

cols <- c("census_id", "id", "stat_id", "species", dayCols, condCols)
tmp <- census[ , cols]
day <- melt(tmp, id.vars = c("species","census_id", "stat_id", "id", condCols))
day$yr=substring(day$variable,4,5)

names(day)[names(day) == 'value'] <- 'day'

dbh <- merge(dbh, day[ , c('stat_id', 'day', 'yr')], all.x = TRUE, all.y = FALSE)
# get in terms of 4-digit year, not 2-digit
match <- data.frame(yr = censusYears2digit, year = as.numeric(censusYears)) 
dbh <- merge(dbh, match, all.x = TRUE, all.y = FALSE)

colMatch = data.frame(censusYears2digit, which(names(dbh) %in% condCols))
names(colMatch) = c('yr','col')
dbh2 = merge(dbh, colMatch, by.x= 'yr', by.y = 'yr')
dbh2$status <- dbh2[cbind(1:nrow(dbh2), dbh2$col)]
dbh2 <- merge(dbh2, census[ , c('stat_id', 'yr1991')], all.x = TRUE, all.y = FALSE)
dbh2$year[dbh2$year == 1991] <- dbh2$yr1991[dbh2$year == 1991]

dbh <- subset(dbh2, status == "L", c('yr', 'year', 'species', 'day', 'census_id', 'id', 'stat_id','value'))
dbh <- subset(dbh, !is.na(value))
## dbh is a core data object going into the stat model

saplings <- subset(dbh2, is.na(status))
## saplings is a core data object going into the stat model

# sum over multiple cores to get tree-level increment
# insert 0 for NA when match non-NA in another core

combineCores <- function(rw) {
    zeroGrowthFlag <- -1
    id <- as.numeric(substring(dimnames(rw)[[1]], 3, 5))
    rw <- rw[!is.na(id), ]
    id <- id[!is.na(id)]
    orient <- substring(dimnames(rw)[[1]], 6, 7)
    trees <- unique(id)
    plot <- unique(as.numeric(substring(dimnames(rw)[[1]], 3, 3)))
    if(length(plot) > 1) stop("multiple plots in object")
    vals <- matrix(zeroGrowthFlag, length(trees), ncol(rw)) # -1 is code for no ring 
    dimnames(vals)[[2]] <- dimnames(rw)[[2]]
    for(i in seq_along(trees)) {
        tree_rw <- rw[id == trees[i], , drop = FALSE]
        cat("processing ", plot, trees[i], substring(dimnames(tree_rw)[[1]], 6, 6), "\n")
        anyNonNA <- apply(tree_rw, 2, function(x) sum(!is.na(x)) > 0)
        vals[i, anyNonNA] <- (colSums(tree_rw, na.rm = TRUE) / (0.5*nrow(tree_rw)))[anyNonNA]
        if(!anyNonNA[1]) {
            lastZero <- max(which(cumsum(vals[i,]) == -(1:ncol(rw))))
            vals[i, 1:lastZero] <- NA  # NA is for before tree existed
        }
    }
    dimnames(vals)[[1]] <- trees
    return(vals)
}

treeRW <- lapply(rwData, combineCores)

# combine rw data to get increment dataset
ringYears <- sapply(treeRW, function(x) as.numeric(dimnames(x)[[2]]))
nTreesWithCores <- sum(sapply(treeRW, function(x) nrow(x)))
incr <- matrix(NA, nTreesWithCores, nT)
dimnames(incr)[[2]] <- years

start <- 1
for(i in seq_along(treeRW)) {
    end <-  start -1 + nrow(treeRW[[i]])
    include <- which(ringYears[[i]] <= lastYear & ringYears[[i]] >= firstYear)
    incr[start:end, as.character(ringYears[[i]][include])] <- treeRW[[i]][ , include]
 
    start <- end + 1
}
dimnames(incr)[[1]] <- unlist(sapply(treeRW, function(x) dimnames(x)[[1]]))

# need to figure out lastDate and put in census
census$lastYear <- 2012
census$lastYear[census$cond11 != "L" & !is.na(census$cond11)] <- 2010
census$lastYear[census$cond01 != "L" & !is.na(census$cond01)] <- 2000
census$lastYear[census$cond91 != "L" & !is.na(census$cond91)] <- census$yr1991[census$cond91 != "L" & !is.na(census$cond91)] - 1
census$lastYear[census$cond75 != "L" & !is.na(census$cond75)] <- 1974
census$lastYear[census$cond69 != "L" & !is.na(census$cond69)] <- 1968

# now walk through trees with rings
firstNoGrowth <- function(x) {
    tmp <- which(x == -1)
    if(length(tmp)) return(as.numeric(names(x)[min(tmp)])) else return(NA)
}
zeroRing <- apply(incr, 1, firstNoGrowth)
zeroRing[is.na(zeroRing)] <- lastYear + 1 # kludge so that 2012 is given as last year in code below for trees with ring in 2012

censusMatches <-  which(census$id %in% names(zeroRing))

census$lastYear[censusMatches] <- zeroRing[as.character(census$id[censusMatches])] - 1


# temporary restrict to one site to quicken calcs
if(F){
census <- subset(census, site == 1)
dbh <- dbh[substring(dbh$id, 1,1)=="1", ]
wh <- substring(dimnames(incr)[[1]], 1 ,1) == "1"
incr <- incr[wh,]
}

tbl <- table(census$species)
taxaMatch <- data.frame(tbl); names(taxaMatch) <- c('species', 'count')
taxaMatch$taxon <- 1:nrow(taxaMatch)

nDBH <- nrow(dbh)
logDobs <- log(dbh$value)
dbh_tree_id <- dbh$stat_id
dbh_day_id <- dbh$day / growDays
dbh_year_id <- getTimeIndex(dbh$year)

nSaplings <- nrow(saplings)
sapling_tree_id <- saplings$stat_id
sapling_year_id <- getTimeIndex(saplings$year)
max_size <- ifelse(saplings$year == 1969, 4.5, 5)

incrMelted <- melt(incr)
names(incrMelted) <- c('id', 'year', 'incr')
incrData <- merge(incrMelted, census[ , c('id', 'stat_id')])
incrData <- subset(incrData, !is.na(incr))
incrData <- subset(incrData, incr != -1)

nWidths <- nrow(incrData)
incr_tree_id <- incrData$stat_id
incr_year_id <- getTimeIndex(incrData$year)
logXobs <- log(incrData$incr)
# fudge 0 incr for now
logXobs[logXobs < log(.02)] <- log(.02)

tmp <- merge(census[ , c('stat_id', 'species')], taxaMatch[ , c('species', 'taxon')], all.x = TRUE, all.y = FALSE)
tmp <- tmp[order(tmp$stat_id), ]
taxon <- tmp$taxon
nTaxa <- nrow(taxaMatch)

n <- nrow(census)
dead <- census$lastYear < lastYear

last_time <- getTimeIndex(census$lastYear)

save(n, nT, nDBH, nWidths, dbh_tree_id, dbh_day_id, dbh_year_id, incr_tree_id, incr_year_id,
     logDobs, logXobs, last_time, census, dbh, incrData, incr, treeMeta, ringMeta,
     nSaplings, sapling_tree_id, sapling_year_id, max_size, taxon, nTaxa,
     file = 'data.Rda')

####################
# old code
####################

if(F){
treeIDs <- lapply(rwData, function(x) dimnames(x)[[1]])

rwYears <- min(sapply(yrs, function(x) x[1])):lastYear
nYears <-  length(rwYears)

trees <- unique(treeMeta$Tag)
# plot these on A's map


for(i in nPlots) {
    diam[[i]] <- cbind(census[[i]]$treeid, matrix(0, nrow = nrow(census[[i]]), ncol = nCensus))
    dimnames(diam[[i]])[[2]] <- c('treeid', censusYears)
    colnames <- paste0('dbh', substring(censusYears, 3, 4))   
    diam[[i]][ , 2:ncol(diam[[i]])] <- as.matrix(census[[i]][ , colnames])
}






if(F) {
    plot(censusFull$xsite, censusFull$ysite)
    for(i in 1:3){
        points(centroids[i, 2:3], col = 'blue', pch = 'X')
        points(census[[i]]$xsite, census[[i]]$ysite, col='red')
    }

    tmp = merge(treeMeta, rbind(census[[1]],census[[2]],census[[3]])[, c('treeid', 'dbh11')], all.x=T,all.y=F)
}

}

