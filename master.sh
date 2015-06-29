#!/bin/bash
# controller script for fitting models to tree ring data
# to estimate biomass increment

# currently focuses on Lyford plot at Harvard Forest

# this is not intended to be run as a full script as later components depend on earlier ones having finished; rather it is intended to allow one to run all of the steps of the model fitting/analysis

# this code is being run under R 3.2.0 and with package versioning controlled by packrat
# restore any packages that are not installed on the system
Rscript -e "require(packrat); packrat::restore()"

source config

export OMP_NUM_THREADS=1

########################################################################
# create directories    ------------------------------------------------
########################################################################

if [ ! -e $outputDir ]; then
    mkdir -p $outputDir
fi
if [ ! -e $dataDir ]; then
    mkdir -p $dataDir
fi
if [ ! -e $tmpDir ]; then
    mkdir -p $tmpDir
fi

########################################################################
# setup for down/uploading from Wiki   ---------------------------------
########################################################################

# get cookie with Wiki authentication info
wget --post-data="u=${WIKI_USERNAME}&p=${WIKI_PASSWORD}&sectok=b7649cb05e87dc0f45d9d91851a7b097&id=start&r=1&do=login" --save-cookies=${dataDir}/my-cookies.txt --keep-session-cookies https://paleon.geography.wisc.edu/doku.php/dw__login

export cookieArgs="--load-cookies=${dataDir}/my-cookies.txt --save-cookies=${dataDir}/my-cookies.txt --keep-session-cookies"

########################################################################
# download Lyford data ------------------------------------------------
########################################################################

cd $dataDir
if [ ! =e lyford ]; then
    mkdir lyford
fi
cd lyford

if [ ! -e lyford.zip ]; then
    wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3B${lyfordVersion}" -O lyford.zip
fi

unzip lyford.zip

if [ ! -e lyford_with_audrey.zip ]; then
    wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3B${audreyCensusVersion}" -O lyford_with_audrey.zip
fi

unzip lyford_with_audrey.zip


if [ ! -e lyford_sample_dates.csv ]; then
    wget $cookieArgs "https://paleon.geography.wisc.edu/lib/exe/fetch.php/data_and_products%3Blyford_sample_dates.csv" -O lyford_sample_dates.csv
fi


ln -s Lyford_Data_13m/LyfordAllPlots.csv
ln -s Lyford_Data_13m/RW

# TMP: for tree 325
mv RW/LF3_ACRU.rw{,.bak}
cp  /tmp/LF3_ACRU.rw RW/.

# this won't run via Rscript for some reason (need to do manually in R)
cd $projectDir
./exportExcelData.R

./process.R

./fit.R

# plot.R contains a function that can be used to make plots of the dbh and increment evolution over time for each tree - plotting both obs and model fit
