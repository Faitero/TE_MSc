#### Create full data table of species lifted alus
## Script will detect the variation trend of each alu exon and if it consistence on 3SS strenght or lowering
library(ggplot2)
require(ggplot2)
require(reshape)
library(system)
debugonce(devtools::install)

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
if (!require("gridGraphics")) {
  install.packages("gridGraphics", dependencies = TRUE)
  library(gridGraphics)
}
if (!require("GMD")) {
  install.packages("GMD", dependencies = TRUE)
  library(GMD)
}

if (!require("ColorPalette")) {
  install.packages("ColorPalette", dependencies = TRUE)
  library(ColorPalette)
}
if (!require("scales")) {
  install.packages("scales", dependencies = TRUE)
  library(scales)
}
if (!require("dplyr")) {
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}



if (!require("devtools")) {
  install.packages("devtools", dependencies = TRUE)
  library(devtools)
}



median.quartile <- function(x){
  out <- quantile(x, probs = c(0.25,0.5,0.75))
  names(out) <- c("ymin","y","ymax")
  return(out) 
}

args<-commandArgs(TRUE)


setwd('/Users/Igor/Dropbox (UCL-MN Team)/hnRNPC.NMD (1).manuscript/Alus_primates_evolution/Data')
setwd('/home/igor/Dropbox (UCL-MN Team)/hnRNPC.NMD (1).manuscript/Alus_primates_evolution/Data')  ## Curro


hg19 <- read.table(("All_Aluexons_3SS_hg19_Distance_hg19_hg19.tab"), sep="\t",  header = TRUE)
hg19$region <- "hg19"

hg38 <- read.table(("All_Aluexons_3SS_hg19_Distance_hg19_hg38.tab"), sep="\t",  header = TRUE)
hg38$region <- "hg38"

panPan1<- read.table("All_Aluexons_3SS_hg19_Distance_hg38_panPan1.tab", sep="\t", header = TRUE)
panPan1$region <- "panPan1"

panTro4<- read.table("All_Aluexons_3SS_hg19_Distance_hg38_panTro4.tab", sep="\t", header = TRUE)
panTro4$region <- "panTro4"

rheMac3<- read.table("All_Aluexons_3SS_hg19_Distance_hg38_rheMac3.tab", sep="\t",  header = TRUE)
rheMac3$region <- "rheMac3"

tarSyr2<- read.table("All_Aluexons_3SS_hg19_Distance_hg38_tarSyr2.tab", sep="\t",  header = TRUE)
tarSyr2$region <- "tarSyr2"

otoGar1<- read.table("All_Aluexons_3SS_hg19_Distance_hg19_otoGar1.tab", sep="\t",  header = TRUE)
otoGar1$region <- "otoGar1"

calJac3<- read.table("All_Aluexons_3SS_hg19_Distance_hg19_calJac3.tab", sep="\t",  header = TRUE)
calJac3$region <- "calJac3"

micMur1<- read.table("All_Aluexons_3SS_hg19_Distance_hg19_micMur1.tab", sep="\t",  header = TRUE)
micMur1$region <- "micMur1"

nomLeu1<- read.table("All_Aluexons_3SS_hg19_Distance_hg19_nomLeu1.tab", sep="\t",  header = TRUE)
nomLeu1$region <- "nomLeu1"

papHam1<- read.table("All_Aluexons_3SS_hg19_Distance_hg19_papHam1.tab", sep="\t",  header = TRUE)
papHam1$region <- "papHam1"

tupBel1<- read.table("All_Aluexons_3SS_hg19_Distance_hg19_tupBel1.tab", sep="\t",  header = TRUE)
tupBel1$region <- "tupBel1"

#genome_list <-list(hg19, hg38, panTro4, panPan1, nomLeu1, papHam1, rheMac3, calJac3, tarSyr2, micMur1, otoGar1, tupBel1) 
#whole <- rbind(
whole <- rbind(hg19, panPan1, panTro4, rheMac3, tarSyr2, otoGar1, calJac3, micMur1, nomLeu1, papHam1, tupBel1)

colMax <- function(whole) sapply(whole, max, na.rm = TRUE) ## get biger value in a column
colMin <- function(whole) sapply(whole, min, na.rm = TRUE) ## get min value in a column
colMax(whole$X3SSS)

delete <- subset(whole, select=c(X3SSS, WU, U1, U2))
colMax(delete)    ## max 3sss 14.69
colMin(delete)

#write.table(whole_final, file = "whole_final.bed", sep = "\t", na = "NA", row.names = FALSE, col.names = TRUE)


names(hg19)<- c("chr_hg19",  "start_hg19",  "end_hg19",  "aluexon",  "position_hg19",	"strand_hg19",	"distance_hg19",	"3SSS_hg19",	"WU_hg19",	"U1_hg19",	"U2_hg19", "region_hg19")

names(hg38)<- c("chr_hg38",  "start_hg38",	"end_hg38",	"aluexon",	"position_hg38",	"strand_hg38",	"distance_hg38",	"3SSS_hg38",	"WU_hg38",	"U1_hg38",	"U2_hg38", "region_hg38")

names(panPan1)<- c("chr_panPan1",	"start_panPan1",	"end_panPan1",	"aluexon",	"position_panPan1",	"strand_panPan1",	"distance_panPan1",	"3SSS_panPan1",	"WU_panPan1",	"U1_panPan1",	"U2_panPan1", "region_panPan1")

names(panTro4)<- c("chr_panTro4",	"start_panTro4",	"end_panTro4",	"aluexon",	"position_panTro4",	"strand_panTro4",	"distance_panTro4",	"3SSS_panTro4",	"WU_panTro4",	"U1_panTro4",	"U2_panTro4", "region_panTro4")

names(rheMac3)<- c("chr_rheMac3",	"start_rheMac3",	"end_rheMac3",	"aluexon",	"position_rheMac3",	"strand_rheMac3",	"distance_rheMac3",	"3SSS_rheMac3",	"WU_rheMac3",	"U1_rheMac3",	"U2_rheMac3", "region_rheMac3")

names(tarSyr2)<- c("chr_tarSyr2",	"start_tarSyr2",	"end_tarSyr2",	"aluexon",	"position_tarSyr2",	"strand_tarSyr2",	"distance_tarSyr2",	"3SSS_tarSyr2",	"WU_tarSyr2",	"U1_tarSyr2",	"U2_tarSyr2", "region_tarSyr2")

names(otoGar1)<- c("chr_otoGar1",	"start_otoGar1",	"end_otoGar1",	"aluexon",	"position_otoGar1",	"strand_otoGar1",	"distance_otoGar1",	"3SSS_otoGar1",	"WU_otoGar1",	"U1_otoGar1",	"U2_otoGar1", "region_otoGar1")

names(calJac3)<- c("chr_calJac3",	"start_calJac3",	"end_calJac3",	"aluexon",	"position_calJac3",	"strand_calJac3",	"distance_calJac3",	"3SSS_calJac3",	"WU_calJac3",	"U1_calJac3",	"U2_calJac3", "region_calJac3")

names(micMur1)<- c("chr_micMur1",	"start_micMur1",	"end_micMur1",	"aluexon",	"position_micMur1",	"strand_micMur1",	"distance_micMur1",	"3SSS_micMur1",	"WU_micMur1",	"U1_micMur1",	"U2_micMur1", "region_micMur1")

names(nomLeu1)<- c("chr_nomLeu1",	"start_nomLeu1",	"end_nomLeu1",	"aluexon",	"position_nomLeu1",	"strand_nomLeu1",	"distance_nomLeu1",	"3SSS_nomLeu1",	"WU_nomLeu1",	"U1_nomLeu1",	"U2_nomLeu1", "region_nomLeu1")

names(papHam1)<- c("chr_papHam1",	"start_papHam1",	"end_papHam1",	"aluexon",	"position_papHam1",	"strand_papHam1",	"distance_papHam1",	"3SSS_papHam1",	"WU_papHam1",	"U1_papHam1",	"U2_papHam1", "region_papHam1")

names(tupBel1)<- c("chr_tupBel1",	"start_tupBel1",	"end_tupBel1",	"aluexon",	"position_tupBel1",	"strand_tupBel1",	"distance_tupBel1",	"3SSS_tupBel1",	"WU_tupBel1",	"U1_tupBel1",	"U2_tupBel1", "region_tupBel1")


#("hg19"="1", "hg38"="2", "panTro4"="3", "panPan1"="4", "nomLeu1"="5", "papHam1"="6", "rheMac3"="7", "calJac3"="8", "tarSyr2"="9", "micMur1"="10", "otoGar1"="11", "tupBel1"="12")  

genome_list <-list(hg19, hg38, panTro4, panPan1, nomLeu1, papHam1, rheMac3, calJac3, tarSyr2, micMur1, otoGar1, tupBel1) 

total <- read.table("All_Aluexons_3SS_hg19_DQ.bed", sep="\t")


colnames(total)<- c("chr_hg19_t",  "start_hg19_t",	"end_hg19_t",	"aluexon",	"position_hg19_t",	"strand_hg19_t",	"distance_hg19_t",	"chr_alu_hg19_t",	"start_alu_hg19_t", "end_alu_hg19_t", "alu_element_hg19_t", "alu_substitutions_hg19_t",	"alu_strand_hg19_t",	
                    "number")
total$number<-NULL

#forJan <- read.table("All_Aluexons_3SS_hg19_DQ.bed", sep="\t") ## My initial table
# initial_table <- read.table("All_Aluexons.bed")                ## Jan Initial table
# forJan1 <- merge.data.frame(forJan, initial_table, all.x=TRUE, by.x=4, by.y=4)
# write.table(forJan1, file = "ForJan.tab", sep = "\t", na = "NA", row.names = FALSE, col.names = TRUE)


library(plyr)
for (genome in genome_list) {
  
  total <- merge.data.frame(total, genome, by="aluexon", all=TRUE) 
  
  
}
#write.table(total, file = "total.tab", sep = "\t", na = "NA", row.names = FALSE, col.names = TRUE)

rnames <- total[,1]
temp <- data.frame(total[,seq(24, ncol(total), 11)])
rownames(temp) <- rnames


#temp <- temp[1:10,]
temp1 <- data.frame(total[,1])
#temp1 <- temp1[1:10,]
for (i in 1:nrow(temp)) {
  lista <- list()  
  lista <- unname(unlist(temp[i,]))
  lista_n <- lista[!is.na(lista)]
  temp1$furthest_all[i] <- as.character(list(lista))
  temp1$furthest[i] <- as.character(tail(lista_n, n=1))
  
}
head(temp1)
class(temp1)

names(temp1) <- c("aluexon", "furthest_all", "furthest")
total_furthest <- merge.data.frame(total, temp1, by="aluexon", all=TRUE) 
total_furthest <- cbind(total, temp1)

unique(total_furthest$furthest)
queries_furthest <- total_furthest
# Substitute hg38 by hg19. And delete hg38 columns
queries_furthest$furthest[queries_furthest$furthest == "hg38"] <- "hg19"
queries_furthest <- queries_furthest[, - c(25:35)]
unique(queries_furthest$furthest)

###################################
###### Check if the 3SSS is increasing the strenght or decresing
#################################
#################### 
#rnames <- total[,1]
#queries <- data.frame(total)
#rownames(queries) <- rnames
#head(queries)
## Same six groups like in One query but in this case we take in acount if the 3SS is higher in human or in caljac3


queries2 <- data.frame(queries_furthest)
#queries2 <- total

## Change colnames to introduce X in front of 3SS score...
names(queries2)[names(queries2) == '3SSS_hg19'] <- 'X3SSS_hg19'
names(queries2)[names(queries2) == '3SSS_panTro4'] <- 'X3SSS_panTro4'
names(queries2)[names(queries2) == '3SSS_panPan1'] <- 'X3SSS_panPan1'
names(queries2)[names(queries2) == '3SSS_nomLeu1'] <- 'X3SSS_nomLeu1'
names(queries2)[names(queries2) == '3SSS_papHam1'] <- 'X3SSS_papHam1'
names(queries2)[names(queries2) == '3SSS_rheMac3'] <- 'X3SSS_rheMac3'
names(queries2)[names(queries2) == '3SSS_calJac3'] <- 'X3SSS_calJac3'
names(queries2)[names(queries2) == '3SSS_tarSyr2'] <- 'X3SSS_tarSyr2'
names(queries2)[names(queries2) == '3SSS_micMur1'] <- 'X3SSS_micMur1'
names(queries2)[names(queries2) == '3SSS_otoGar1'] <- 'X3SSS_otoGar1'
names(queries2)[names(queries2) == '3SSS_tupBel1'] <- 'X3SSS_tupBel1'
colnames(queries2)
# Substitute NAs in X3SSS_sp with zeros
#queries2[,c(20,31,42,53,64,75,86,97,108,119,130)][is.na(queries2[c(20,31,42,53,64,75,86,97,108,119,130)])] <- as.numeric("0")  


test <- queries2[,c(20,31,42,53,64,75,86,97,108,119,130)]
head(test)

                                                                 


#queries2[,c(20,31,42,53,64,75,86,97,108,119,130)][queries2[,c(20,31,42,53,64,75,86,97,108,119,130)] == "NA"   )] <- 0





## Transfor 
head(queries2)
#d <- transform(d, new = min / count2.freq) # example

## 1. via `[` and character indexes
d[, "new"] <- d[, "min"] / d[, "count2.freq"]

## 2. via `[` with numeric indices
d[, 3] <- d[, 1] / d[, 2]

## 3. via `$`
d$new <- d$min / d$count2.freq



setwd('/media/igor/DATA/UCL/Evolution_Alus/New3SS/Output_R' )


test$cinco3SS_f <- test$cinco3SS
test$accumulate <- as.double(0)
test$cuatro3SS_f <- apply(test[,c('cuatro3SS', 'cinco3SS', 'accumulate')], 1, function(x) {TREND(x[1],x[2], x[3]) })
test$accumulate <- apply(test[,c('cuatro3SS_f', 'accumulate')], 1, function(x) TREND_accum(x[1],x[2]))
test$tres3SS_f <- apply(test[,c('tres3SS', 'cuatro3SS', 'accumulate')], 1, function(x) TREND(x[1],x[2], x[3]))
test$accumulate <- apply(test[,c('tres3SS_f', 'accumulate')], 1, function(x) TREND_accum(x[1],x[2]))

test$dos3SS_f <- apply(test[,c('dos3SS', 'tres3SS', 'accumulate')], 1, function(x) TREND(x[1],x[2], x[3]))
test$accumulate <- apply(test[,c('dos3SS_f', 'accumulate')], 1, function(x) TREND_accum(x[1],x[2]))


test


TREND_accum(NULL, NULL)
TREND_accum(1, NULL)
TREND_accum(NULL, 1)
TREND_accum(NA, -1)
TREND_accum(NA, NA)
TREND_accum(1, -1)
TREND_accum(1, 0.00)
TREND_accum(0.00, -1)
TREND_accum(0.00, 0.000)
TREND_accum(1, 3)
TREND_accum(3, -1)

debugonce(devtools::install)
options(error = browser())

TREND_accum <- function(actual, accumulate){
  
  #print(paste0(actual))
  #print(paste0(accumulate))
  actual <- as.numeric(actual)
  accumulate <- as.numeric(accumulate)

  if(is.null(actual) & is.null(accumulate)) {return(0)}
  if(is.null(actual) & accumulate == 0) {return(0)}
  if(actual=="NULL" & accumulate == 0) {return(0)}
  if(!is.null(actual) & !is.null(accumulate)) {return(as.numeric(actual) + accumulate)}
  #if((is.na(actual) | is.null(actual)) & !is.null(accumulate)) {return(accumulate)}
  
  if(is.null(actual) & !is.null(accumulate)) {return(as.numeric(accumulate))}    
  if(!is.null(actual) & is.null(accumulate)) {return(as.numeric(actual))}

  
}


TREND <- function(actual, old, accumulate){
  
#   if(is.null(actual) | is.null(old) | is.null(accumulate)) { print(paste0("NULLLLL   Actual: ", actual, "OLD: ", old, "ACUMULATE", accumulate))
#     return(NULL)}
#   if(!is.double(actual) | !is.double(old) | !is.double(accumulate)) { print(paste0("DOUBLEEEE   Actual: ", actual, "  OLD: ", old, "  ACUMULATE", accumulate))
#     return(NULL)}
#   if(!is.numeric(actual) | !is.numeric(old) | !is.numeric(accumulate)) { print(paste0("NOT NUMERIC  Actual: ", actual, "  OLD: ", old, "  ACUMULATE", accumulate))
#     return(NULL)}
  if(is.null(actual)) {return(NULL)}
  if(actual == 0) {return(NULL)}
  if(is.null(actual) & is.null(old) & is.null(accumulate)) {return(NULL)}
  ### tests
  
  if(!is.null(actual) & is.null(old) & !is.null(accumulate)) {return(actual - accumulate)}
  if(!is.null(actual) & is.null(old) & is.null(accumulate)) {return(actual)}
  
#   if(actual > 0 & old > 0 & !is.null(accumulate)) {return((actual-old) + accumulate)}
#   if(actual > 0 & old < 0 & !is.null(accumulate)) {return((actual-old) + accumulate)}
#   if(actual < 0 & old < 0 & !is.null(accumulate)) {return((actual-old) + accumulate)}
#   if(actual < 0 & old > 0 & !is.null(accumulate)) {return((actual-old) + accumulate)}
#   
  
  else {
    return((actual-old) + accumulate)
  }
  
}


if(actual > 0 & old > 0 & !is.null(accumulate)) {return((actual-old) + accumulate)}
if(actual > 0 & old < 0 & !is.null(accumulate)) {return((actual-old) + accumulate)}
if(actual < 0 & old < 0 & !is.null(accumulate)) {return((actual-old) + accumulate)}
if(actual < 0 & old > 0 & !is.null(accumulate)) {return((actual-old) + accumulate)}



queries <- queries2

queries2<- queries[, c('aluexon', 'X3SSS_hg19', 'X3SSS_panTro4', 'X3SSS_panPan1', 'X3SSS_nomLeu1', 'X3SSS_papHam1', 'X3SSS_rheMac3', 'X3SSS_calJac3', 'X3SSS_tarSyr2', 'X3SSS_micMur1', 'X3SSS_otoGar1', 'X3SSS_tupBel1')]
head(queries2)

rownames(queries2) <- queries2$aluexon

queries2 <- as.matrix(queries2) 


## 
##as.numeric(unlist(queries2))

queries2$X3SSS_tupBel1_f <- queries2$X3SSS_tupBel1
queries2$accumulate <- (0.00)
queries2$X3SSS_otoGar1_f <- apply(queries2[ , c('X3SSS_otoGar1','X3SSS_tupBel1', 'accumulate')], 1, function(x) TREND(x[1],x[2],x[3]))
queries2$accumulate  <- apply(queries2[ , c('X3SSS_otoGar1_f', 'accumulate')], 1, function(x) TREND_accum(x[1],x[2]))
queries2$X3SSS_otoGar1_f <- as.character(queries2$X3SSS_otoGar1_f)
queries2$accumulate<- as.numeric(queries2$accumulate)

queries2$X3SSS_micMur1_f <- apply(queries2[ , c('X3SSS_micMur1', 'X3SSS_otoGar1', 'accumulate')], 1, function(x) as.numeric(TREND(x[1],x[2],x[3])))
queries2$accumulate  <- apply(queries2[ , c('X3SSS_micMur1_f', 'accumulate')], 1, function(x) as.numeric(TREND_accum(x[1],x[2])))
#queries2$X3SSS_micMur1_f <- as.character(queries2$X3SSS_micMur1_f)
#queries2$accumulate<- as.character(queries2$accumulate)

queries2$X3SSS_tarSyr2_f <- apply(queries2[ , c('X3SSS_tarSyr2', 'X3SSS_micMur1', 'accumulate')], 1, function(x) as.numeric(TREND(x[1],x[2],x[3])))
queries2$accumulate   <- apply(queries2[ , c('X3SSS_tarSyr2_f', 'accumulate')], 1, function(x) as.numeric(TREND_accum(x[1],x[2])))
#queries2$X3SSS_tarSyr2_f<- as.character(queries2$X3SSS_tarSyr2_f)
#queries2$accumulate<- as.character(queries2$accumulate)

queries2$X3SSS_calJac3_f <- apply(queries2[ , c('X3SSS_calJac3', 'X3SSS_tarSyr2', 'accumulate')], 1, function(x) as.numeric(TREND(x[1],x[2],x[3])))
queries2$accumulate   <- apply(queries2[ , c('X3SSS_calJac3_f', 'accumulate')], 1, function(x) as.numeric(TREND_accum(x[1],x[2])))
#queries2$X3SSS_calJac3_f <- as.character(queries2$X3SSS_calJac3_f)
#queries2$accumulate<- as.character(queries2$accumulate)

queries2$X3SSS_rheMac3_f <- apply(queries2[ , c('X3SSS_rheMac3', 'X3SSS_calJac3', 'accumulate')], 1, function(x) as.numeric(TREND(x[1],x[2],x[3])))
queries2$accumulate   <- apply(queries2[ , c('X3SSS_rheMac3_f', 'accumulate')], 1, function(x) as.numeric(TREND_accum(x[1],x[2])))
#queries2$X3SSS_rheMac3_f <- as.character(queries2$X3SSS_rheMac3_f)
#queries2$accumulate<- as.character(queries2$accumulate)

queries2$X3SSS_papHam1_f <- apply(queries2[ , c('X3SSS_papHam1', 'X3SSS_rheMac3', 'accumulate')], 1, function(x) as.numeric(TREND(x[1],x[2],x[3])))
queries2$accumulate   <- apply(queries2[ , c('X3SSS_papHam1_f', 'accumulate')], 1, function(x) as.numeric(TREND_accum(x[1],x[2])))
#queries2$X3SSS_papHam1_f <- as.character(queries2$X3SSS_papHam1_f)
#queries2$accumulate<- as.character(queries2$accumulate)

queries2$X3SSS_nomLeu1_f <- apply(queries2[ , c('X3SSS_nomLeu1', 'X3SSS_papHam1', 'accumulate')], 1, function(x) as.numeric(TREND(x[1],x[2],x[3])))
queries2$accumulate   <- apply(queries2[ , c('X3SSS_nomLeu1_f', 'accumulate')], 1, function(x) as.numeric(TREND_accum(x[1],x[2])))
#queries2$X3SSS_nomLeu1_f <- as.character(queries2$X3SSS_nomLeu1_f)
#queries2$accumulate<- as.character(queries2$accumulate)

queries2$X3SSS_panPan1_f <- apply(queries2[ , c('X3SSS_panPan1', 'X3SSS_nomLeu1', 'accumulate')], 1, function(x) as.numeric(TREND(x[1],x[2],x[3])))
queries2$accumulate   <- apply(queries2[ , c('X3SSS_panPan1_f', 'accumulate')], 1, function(x) as.numeric(TREND_accum(x[1],x[2])))
#queries2$X3SSS_panPan1_f <- as.character(queries2$X3SSS_panPan1_f)
#queries2$accumulate<- as.character(queries2$accumulate)

queries2$X3SSS_panTro4_f <- apply(queries2[ , c('X3SSS_panTro4', 'X3SSS_panPan1', 'accumulate')], 1, function(x) as.numeric(TREND(x[1],x[2],x[3])))
queries2$accumulate   <- apply(queries2[ , c('X3SSS_panTro4_f', 'accumulate')], 1, function(x) as.numeric(TREND_accum(x[1],x[2])))
#queries2$X3SSS_panTro4_f <- as.character(queries2$X3SSS_panTro4_f)
#queries2$accumulate<- as.character(queries2$accumulate)

queries2$X3SSS_hg19_f <- apply(queries2[ , c('X3SSS_hg19', 'X3SSS_panTro4', 'accumulate')], 1, function(x) as.numeric(TREND(x[1],x[2],x[3])))
queries2$accumulate   <- apply(queries2[ , c('X3SSS_hg19_f', 'accumulate')], 1, function(x) as.numeric(TREND_accum(x[1],x[2])))
#queries2$X3SSS_hg19_f <- as.character(queries2$X3SSS_hg19_f)
#queries2$accumulate<- as.character(queries2$accumulate)

sapply(queries2, class)
sapply(queries2, mode)
attributes(queries2)
class(queries2)

variations <- queries2[, c('aluexon', 'X3SSS_hg19_f', 'X3SSS_hg19', 'X3SSS_panTro4_f', 'X3SSS_panTro4', 'X3SSS_panPan1_f', 'X3SSS_panPan1', 'X3SSS_nomLeu1_f', 'X3SSS_nomLeu1', 'X3SSS_papHam1_f', 'X3SSS_papHam1', 'X3SSS_rheMac3_f', 'X3SSS_rheMac3', 'X3SSS_calJac3_f', 'X3SSS_calJac3', 'X3SSS_tarSyr2_f', 'X3SSS_tarSyr2', 'X3SSS_micMur1_f', 'X3SSS_micMur1' ,'X3SSS_otoGar1_f', 'X3SSS_otoGar1', 'X3SSS_tupBel1_f',  'X3SSS_tupBel1', 'accumulate')]
head(variations)


variations <- queries2[, c('aluexon', 'X3SSS_hg19_f', 'X3SSS_panTro4_f', 'X3SSS_panPan1_f', 'X3SSS_nomLeu1_f', 'X3SSS_papHam1_f', 'X3SSS_rheMac3_f', 'X3SSS_calJac3_f', 'X3SSS_tarSyr2_f', 'X3SSS_micMur1_f', 'X3SSS_otoGar1_f', 'X3SSS_tupBel1_f')]
head(variations)


data_frame <- queries2
data_frame <- variations
ncol(data_frame)
head(data_frame[,138])
head(hg19)

data_frame <- as.matrix(queries2)
data_frame <- unlist(data_frame[,'X3SSS_hg19_f'])

data_frame <- data.frame(X3SSS_hg19_f = unlist(X3SSS_hg19_f))

is.null(data_frame[,149])

### New function without hg38
wide_to_long <- function(data_frame, new_column){
  
  # Function transform the wide table of lift ovet to long table in where each row has a factor that represent its procedence and used for later plotting
  #new_column <- "cal_high_3ss"
  #data_frame<-queries
  
  colnames(data_frame) <- NULL
  hg19_q <- data.frame(data_frame[,1], data_frame[,14:24], data_frame[,149])
  hg19_q <- cbind(hg19, data_frame$X3SSS_hg19_f)
  pantro4_q <- data.frame(data_frame[,1], data_frame[,25:35])
  panPan1_q <-  data.frame(data_frame[,1], data_frame[,36:46])
  nomLeu1_q <- data.frame(data_frame[,1], data_frame[,47:57])
  papHam1_q <- data.frame(data_frame[,1], data_frame[,58:68])
  rheMac3_q <- data.frame(data_frame[,1], data_frame[,69:79])
  calJac3_q <- data.frame(data_frame[,1], data_frame[,80:90])
  tarSyr2_q <- data.frame(data_frame[,1], data_frame[,91:101])
  micMur1_q <- data.frame(data_frame[,1], data_frame[,102:112])
  otoGar1_q <- data.frame(data_frame[,1], data_frame[,113:123])
  tupBel1_q <- data.frame(data_frame[,1], data_frame[,124:134])
  
  out_data_frame <- rbind(hg19_q, pantro4_q, panPan1_q, nomLeu1_q, papHam1_q, rheMac3_q, calJac3_q, tarSyr2_q, micMur1_q, otoGar1_q, tupBel1_q)
  
  colnames(out_data_frame) <- c("aluexon", "chr", "start", "end", "position", "strand", "distance_to_alu", "X3SSS", "WU", "U1", "U2", "region")# , "group")  
  
  out_data_frame$group <- as.factor(new_column)
  out_data_frame <- out_data_frame[complete.cases(out_data_frame),]
  
  return(out_data_frame)
  
}


#test <- queries2[,c(20,31,42,53,64,75,86,97,108,119,140,130,138,139)]
head(test)



library(dplyr)
test1 <- filter(test, abs(test$uno3SS) <= abs(test$dos3ss*abs(1.2)))
head(test1)



X3ss_higher_hum <- subset(queries, queries$X3SSS_hg19 > queries$X3SSS_calJac3) # 2123
X3ss_higher_calJac <- subset(queries, queries$X3SSS_hg19 <= queries$X3SSS_calJac3) # 1031


#and we perform the same clustering:
## CalJack 3ss higher than 3
X3ss_higher_hum_cal_high_3ss <- subset(X3ss_higher_hum, X3ss_higher_hum$"X3SSS_calJac3" >=3) # 825
X3ss_higher_calJac_cal_high_3ss <- subset(X3ss_higher_calJac, X3ss_higher_calJac$"X3SSS_calJac3" >=3) # 581
## CalJack lower than 3
X3ss_higher_hum_cal_low_3ss <- subset(X3ss_higher_hum, X3ss_higher_hum$"X3SSS_calJac3" <3) # 1298
X3ss_higher_calJac_cal_low_3ss <- subset(X3ss_higher_calJac, X3ss_higher_calJac$"X3SSS_calJac3" <3) # 450

## Check that the number of row is the same in all of them
nrow(X3ss_higher_hum) == nrow(X3ss_higher_hum_cal_high_3ss) + nrow(X3ss_higher_hum_cal_low_3ss) 
nrow(X3ss_higher_calJac) == nrow(X3ss_higher_calJac_cal_high_3ss) + nrow(X3ss_higher_calJac_cal_low_3ss)


X3ss_higher_hum_cal_high_3ss <- wide_to_long(X3ss_higher_hum_cal_high_3ss, "X3ss_higher_hum_cal_high_3ss")
X3ss_higher_calJac_cal_high_3ss <- wide_to_long(X3ss_higher_calJac_cal_high_3ss, "X3ss_higher_calJac_cal_high_3ss")
X3ss_higher_hum_cal_low_3ss <- wide_to_long(X3ss_higher_hum_cal_low_3ss, "X3ss_higher_hum_cal_low_3ss")
X3ss_higher_calJac_cal_low_3ss <- wide_to_long(X3ss_higher_calJac_cal_low_3ss, "X3ss_higher_calJac_cal_low_3ss")



final_data_frame <- rbind(X3ss_higher_hum_cal_high_3ss, X3ss_higher_calJac_cal_high_3ss, X3ss_higher_hum_cal_low_3ss, X3ss_higher_calJac_cal_low_3ss)
## Check that the number of row is the same in all of them
nrow(final_data_frame) == nrow(X3ss_higher_hum_cal_high_3ss) + nrow(X3ss_higher_calJac_cal_high_3ss) + nrow(X3ss_higher_hum_cal_low_3ss) + nrow(X3ss_higher_calJac_cal_low_3ss)


final_data_frame[] <- lapply(final_data_frame, as.character) ## Transfor all to character

final_data_frame[, 7:11] <- sapply(final_data_frame[, c(7:11)], as.numeric)  ## rename some columns to number
final_data_frame$group <- as.factor(final_data_frame$group)
final_data_frame$region <- as.factor(final_data_frame$region)

setwd('/media/igor/DATA/UCL/Evolution_Alus/New3SS/Output_R' )


p_meds <- ddply(final_data_frame, .(region, group), summarise, med = mean(X3SSS))
p_meds


a<-  ggplot(final_data_frame, aes(factor(region), X3SSS), alpha = 1, colour = "black") + 
  geom_boxplot() + 
  #scale_y_log10() + 
  theme_bw() +
  ggtitle("median quartile 25%, 50%, 75%") + 
  xlab("") + 
  ylab("length (nt)") +
  theme(text=element_text(size=12),axis.text=element_text(size=12), axis.title=element_text(size=12,face="plain")) +
  scale_x_discrete(limits=c("hg19", "hg38", "panTro4", "panPan1", "nomLeu1", "papHam1", "rheMac3", "calJac3", "tarSyr2", "micMur1", "otoGar1", "tupBel1"))

a 



count_table <- data.frame(table(final_data_frame$region, final_data_frame$group))
colnames(count_table) <- c("region", "group", "freq")
count_table


library(plyr)
### quartile <- median.quartile(clusters$length)     Nejc what its clusters

dodge <- position_dodge(width = 0.9)


pdf("3SSS_clustering2.pdf", width=20, height=13)
ALL<- ggplot(final_data_frame, aes(factor(region), X3SSS, colour = group), alpha = 1, width = 0.5) + 
  #geom_violin(position = dodge) + 
  geom_boxplot(width=1.1, position = dodge, alpha = 1, outlier.shape = NA ) +    ### uncoment if you need the boxplot inside
  theme_bw() +
  scale_y_continuous(limits = c(-15, 16)) +
  #scale_y_log10() +
  ggtitle("3SSS_OWN_clustering2.pdf") + 
  xlab("") + 
  ylab("3ss score") +
  scale_colour_manual(values=c("#00441b", "#006d2c", "#0570b0", "#3690c0", "#810f7c", "#88419d", "#8c6bb1"))+  #("#00441b", "#006d2c", "#238b45", "#0570b0", "#3690c0", "#74a9cf", "#810f7c", "#88419d", "#8c6bb1")
  #scale_colour_manual(values=c("#242B38","#2165E8","#2CA3B5")) +
  theme(text=element_text(size=12),axis.text=element_text(size=12), axis.title=element_text(size=12,face="plain")) +
  geom_text(data = p_meds, aes(x = region, y = med, label = format(med, digits=2)), size = 2.5, vjust = -1.5, position = dodge) +
  geom_text(data = count_table, aes(x = region, y = -15, label = freq), size = 3, vjust = 0, position = dodge, angle=45) +
  scale_x_discrete(limits=c("hg19", "hg38", "panTro4", "panPan1", "nomLeu1", "papHam1", "rheMac3", "calJac3", "tarSyr2", "micMur1", "otoGar1", "tupBel1"))

ALL
dev.off()

