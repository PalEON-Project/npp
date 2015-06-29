#!/usr/bin/Rscript
library(XLConnect)

source("config")

setwd(file.path(dataDir, "lyford", "Lyford_Data_13m"))
wb = loadWorkbook("LyfordPlots_LogsheetAll_FINAL.xls")
incr = readWorksheet(wb, sheet = "Logsheet", header = TRUE)

setwd(file.path('..', 'PalEON_LyfordBiomassData', 'yBarkerPlotkinDataFileMap'))
wb = loadWorkbook("Lyford_Paleon.xlsx")
census = readWorksheet(wb, sheet = "Sheet1", header = TRUE)

setwd(file.path("..", ".."))

write.csv(incr, file = "lyford_13m_increment_meta.csv")
write.csv(census, file = "lyford_census.csv")
