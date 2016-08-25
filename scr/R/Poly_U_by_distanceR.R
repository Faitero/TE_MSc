#### Graphs for Alu exons evolution
library(ggplot2)
require(ggplot2)

## Created on 23 april 2016
## @author: Igor Ruiz de los Mozos


median.quartile <- function(x){
  out <- quantile(x, probs = c(0.25,0.5,0.75))
  names(out) <- c("ymin","y","ymax")
  return(out) 
}

args<-commandArgs(TRUE)

##################################################
## Lowest
##################################################


##### Whole sequence 320

hg38 <- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg19_hg38_valid-320.tab", sep="\t")
#TDP43 <- read.table(args[1], sep="\t")
hg38$length <- hg38$V4 
hg38$region <- "hg38"

panPan1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg38_panPan1_valid-320.tab", sep="\t")
panPan1$length <- panPan1$V4
panPan1$region <- "panPan1"

panTro4<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg38_panTro4_valid-320.tab", sep="\t")
panTro4$length <- panTro4$V4 
panTro4$region <- "panTro4"


rheMac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg38_rheMac3_valid-320.tab", sep="\t")
rheMac3$length <- rheMac3$V4 
rheMac3$region <- "rheMac3"


tarSyr2<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg38_tarSyr2_valid-320.tab", sep="\t")
tarSyr2$length <- tarSyr2$V4 
tarSyr2$region <- "tarSyr2"


otoGar1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg19_otoGar1_valid-320.tab", sep="\t")
otoGar1$length <- otoGar1$V4 
otoGar1$region <- "otoGar1"


calJac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg19_calJac3_valid-320.tab", sep="\t")
calJac3$length <- calJac3$V4 
calJac3$region <- "calJac3"


micMur1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg19_micMur1_valid-320.tab", sep="\t")
micMur1$length <- micMur1$V4 
micMur1$region <- "micMur1"

nomLeu1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg19_nomLeu1_valid-320.tab", sep="\t")
nomLeu1$length <- nomLeu1$V4 
nomLeu1$region <- "nomLeu1"


papHam1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg19_papHam1_valid-320.tab", sep="\t")
papHam1$length <- papHam1$V4 
papHam1$region <- "papHam1"


tupBel1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg19_tupBel1_valid-320.tab", sep="\t")
tupBel1$length <- tupBel1$V4 
tupBel1$region <- "tupBel1"


W_lowest <- rbind(hg38, panPan1, panTro4, rheMac3, tarSyr2, otoGar1, calJac3, micMur1, nomLeu1, papHam1, tupBel1)
W_lowest$alu_exon <- "Whole_lowest"



#### First 70 nt = U1  ####


hg38 <- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg19_hg38_valid-70.tab", sep="\t")
#TDP43 <- read.table(args[1], sep="\t")
hg38$length <- hg38$V4 
hg38$region <- "hg38"

panPan1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg38_panPan1_valid-70.tab", sep="\t")
panPan1$length <- panPan1$V4
panPan1$region <- "panPan1"

panTro4<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg38_panTro4_valid-70.tab", sep="\t")
panTro4$length <- panTro4$V4 
panTro4$region <- "panTro4"


rheMac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg38_rheMac3_valid-70.tab", sep="\t")
rheMac3$length <- rheMac3$V4 
rheMac3$region <- "rheMac3"


tarSyr2<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg38_tarSyr2_valid-70.tab", sep="\t")
tarSyr2$length <- tarSyr2$V4 
tarSyr2$region <- "tarSyr2"


otoGar1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg19_otoGar1_valid-70.tab", sep="\t")
otoGar1$length <- otoGar1$V4 
otoGar1$region <- "otoGar1"


calJac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg19_calJac3_valid-70.tab", sep="\t")
calJac3$length <- calJac3$V4 
calJac3$region <- "calJac3"


micMur1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg19_micMur1_valid-70.tab", sep="\t")
micMur1$length <- micMur1$V4 
micMur1$region <- "micMur1"

nomLeu1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg19_nomLeu1_valid-70.tab", sep="\t")
nomLeu1$length <- nomLeu1$V4 
nomLeu1$region <- "nomLeu1"


papHam1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg19_papHam1_valid-70.tab", sep="\t")
papHam1$length <- papHam1$V4 
papHam1$region <- "papHam1"


tupBel1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg19_tupBel1_valid-70.tab", sep="\t")
tupBel1$length <- tupBel1$V4 
tupBel1$region <- "tupBel1"


U1_lowest <- rbind(hg38, panPan1, panTro4, rheMac3, tarSyr2, otoGar1, calJac3, micMur1, nomLeu1, papHam1, tupBel1)
U1_lowest$alu_exon <- "U1_lowest"


######## last 250 nt U2 #########


hg38 <- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg19_hg38_valid-250.tab", sep="\t")
#TDP43 <- read.table(args[1], sep="\t")
hg38$length <- hg38$V4 
hg38$region <- "hg38"

panPan1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg38_panPan1_valid-250.tab", sep="\t")
panPan1$length <- panPan1$V4
panPan1$region <- "panPan1"

panTro4<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg38_panTro4_valid-250.tab", sep="\t")
panTro4$length <- panTro4$V4 
panTro4$region <- "panTro4"


rheMac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg38_rheMac3_valid-250.tab", sep="\t")
rheMac3$length <- rheMac3$V4 
rheMac3$region <- "rheMac3"


tarSyr2<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg38_tarSyr2_valid-250.tab", sep="\t")
tarSyr2$length <- tarSyr2$V4 
tarSyr2$region <- "tarSyr2"


otoGar1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg19_otoGar1_valid-250.tab", sep="\t")
otoGar1$length <- otoGar1$V4 
otoGar1$region <- "otoGar1"


calJac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg19_calJac3_valid-250.tab", sep="\t")
calJac3$length <- calJac3$V4 
calJac3$region <- "calJac3"


micMur1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg19_micMur1_valid-250.tab", sep="\t")
micMur1$length <- micMur1$V4 
micMur1$region <- "micMur1"

nomLeu1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg19_nomLeu1_valid-250.tab", sep="\t")
nomLeu1$length <- nomLeu1$V4 
nomLeu1$region <- "nomLeu1"


papHam1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg19_papHam1_valid-250.tab", sep="\t")
papHam1$length <- papHam1$V4 
papHam1$region <- "papHam1"


tupBel1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/lowest/All_Aluexons_3SS_Q_lowest_hg19_tupBel1_valid-250.tab", sep="\t")
tupBel1$length <- tupBel1$V4 
tupBel1$region <- "tupBel1"




U2_lowest <- rbind(hg38, panPan1, panTro4, rheMac3, tarSyr2, otoGar1, calJac3, micMur1, nomLeu1, papHam1, tupBel1)
U2_lowest$alu_exon <- "U2_lowest"




##################################################
## low
##################################################


##### Whole sequence 320

hg38 <- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg19_hg38_valid-320.tab", sep="\t")
#TDP43 <- read.table(args[1], sep="\t")
hg38$length <- hg38$V4 
hg38$region <- "hg38"

panPan1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg38_panPan1_valid-320.tab", sep="\t")
panPan1$length <- panPan1$V4
panPan1$region <- "panPan1"

panTro4<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg38_panTro4_valid-320.tab", sep="\t")
panTro4$length <- panTro4$V4 
panTro4$region <- "panTro4"


rheMac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg38_rheMac3_valid-320.tab", sep="\t")
rheMac3$length <- rheMac3$V4 
rheMac3$region <- "rheMac3"


tarSyr2<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg38_tarSyr2_valid-320.tab", sep="\t")
tarSyr2$length <- tarSyr2$V4 
tarSyr2$region <- "tarSyr2"


otoGar1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg19_otoGar1_valid-320.tab", sep="\t")
otoGar1$length <- otoGar1$V4 
otoGar1$region <- "otoGar1"


calJac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg19_calJac3_valid-320.tab", sep="\t")
calJac3$length <- calJac3$V4 
calJac3$region <- "calJac3"


micMur1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg19_micMur1_valid-320.tab", sep="\t")
micMur1$length <- micMur1$V4 
micMur1$region <- "micMur1"

nomLeu1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg19_nomLeu1_valid-320.tab", sep="\t")
nomLeu1$length <- nomLeu1$V4 
nomLeu1$region <- "nomLeu1"


papHam1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg19_papHam1_valid-320.tab", sep="\t")
papHam1$length <- papHam1$V4 
papHam1$region <- "papHam1"


tupBel1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg19_tupBel1_valid-320.tab", sep="\t")
tupBel1$length <- tupBel1$V4 
tupBel1$region <- "tupBel1"


W_low <- rbind(hg38, panPan1, panTro4, rheMac3, tarSyr2, otoGar1, calJac3, micMur1, nomLeu1, papHam1, tupBel1)
W_low$alu_exon <- "Whole_low"



#### First 70 nt = U1  ####


hg38 <- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg19_hg38_valid-70.tab", sep="\t")
#TDP43 <- read.table(args[1], sep="\t")
hg38$length <- hg38$V4 
hg38$region <- "hg38"

panPan1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg38_panPan1_valid-70.tab", sep="\t")
panPan1$length <- panPan1$V4
panPan1$region <- "panPan1"

panTro4<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg38_panTro4_valid-70.tab", sep="\t")
panTro4$length <- panTro4$V4 
panTro4$region <- "panTro4"


rheMac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg38_rheMac3_valid-70.tab", sep="\t")
rheMac3$length <- rheMac3$V4 
rheMac3$region <- "rheMac3"


tarSyr2<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg38_tarSyr2_valid-70.tab", sep="\t")
tarSyr2$length <- tarSyr2$V4 
tarSyr2$region <- "tarSyr2"


otoGar1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg19_otoGar1_valid-70.tab", sep="\t")
otoGar1$length <- otoGar1$V4 
otoGar1$region <- "otoGar1"


calJac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg19_calJac3_valid-70.tab", sep="\t")
calJac3$length <- calJac3$V4 
calJac3$region <- "calJac3"


micMur1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg19_micMur1_valid-70.tab", sep="\t")
micMur1$length <- micMur1$V4 
micMur1$region <- "micMur1"

nomLeu1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg19_nomLeu1_valid-70.tab", sep="\t")
nomLeu1$length <- nomLeu1$V4 
nomLeu1$region <- "nomLeu1"


papHam1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg19_papHam1_valid-70.tab", sep="\t")
papHam1$length <- papHam1$V4 
papHam1$region <- "papHam1"


tupBel1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg19_tupBel1_valid-70.tab", sep="\t")
tupBel1$length <- tupBel1$V4 
tupBel1$region <- "tupBel1"


U1_low <- rbind(hg38, panPan1, panTro4, rheMac3, tarSyr2, otoGar1, calJac3, micMur1, nomLeu1, papHam1, tupBel1)
U1_low$alu_exon <- "U1_low"


######## last 250 nt U2 #########


hg38 <- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg19_hg38_valid-250.tab", sep="\t")
#TDP43 <- read.table(args[1], sep="\t")
hg38$length <- hg38$V4 
hg38$region <- "hg38"

panPan1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg38_panPan1_valid-250.tab", sep="\t")
panPan1$length <- panPan1$V4
panPan1$region <- "panPan1"

panTro4<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg38_panTro4_valid-250.tab", sep="\t")
panTro4$length <- panTro4$V4 
panTro4$region <- "panTro4"


rheMac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg38_rheMac3_valid-250.tab", sep="\t")
rheMac3$length <- rheMac3$V4 
rheMac3$region <- "rheMac3"


tarSyr2<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg38_tarSyr2_valid-250.tab", sep="\t")
tarSyr2$length <- tarSyr2$V4 
tarSyr2$region <- "tarSyr2"


otoGar1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg19_otoGar1_valid-250.tab", sep="\t")
otoGar1$length <- otoGar1$V4 
otoGar1$region <- "otoGar1"


calJac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg19_calJac3_valid-250.tab", sep="\t")
calJac3$length <- calJac3$V4 
calJac3$region <- "calJac3"


micMur1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg19_micMur1_valid-250.tab", sep="\t")
micMur1$length <- micMur1$V4 
micMur1$region <- "micMur1"

nomLeu1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg19_nomLeu1_valid-250.tab", sep="\t")
nomLeu1$length <- nomLeu1$V4 
nomLeu1$region <- "nomLeu1"


papHam1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg19_papHam1_valid-250.tab", sep="\t")
papHam1$length <- papHam1$V4 
papHam1$region <- "papHam1"


tupBel1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/low/All_Aluexons_3SS_Q_low_hg19_tupBel1_valid-250.tab", sep="\t")
tupBel1$length <- tupBel1$V4 
tupBel1$region <- "tupBel1"




U2_low <- rbind(hg38, panPan1, panTro4, rheMac3, tarSyr2, otoGar1, calJac3, micMur1, nomLeu1, papHam1, tupBel1)
U2_low$alu_exon <- "U2_low"



##################################################
## moderate
##################################################


##### Whole sequence 320

hg38 <- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg19_hg38_valid-320.tab", sep="\t")
#TDP43 <- read.table(args[1], sep="\t")
hg38$length <- hg38$V4 
hg38$region <- "hg38"

panPan1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg38_panPan1_valid-320.tab", sep="\t")
panPan1$length <- panPan1$V4
panPan1$region <- "panPan1"

panTro4<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg38_panTro4_valid-320.tab", sep="\t")
panTro4$length <- panTro4$V4 
panTro4$region <- "panTro4"


rheMac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg38_rheMac3_valid-320.tab", sep="\t")
rheMac3$length <- rheMac3$V4 
rheMac3$region <- "rheMac3"


tarSyr2<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg38_tarSyr2_valid-320.tab", sep="\t")
tarSyr2$length <- tarSyr2$V4 
tarSyr2$region <- "tarSyr2"


otoGar1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg19_otoGar1_valid-320.tab", sep="\t")
otoGar1$length <- otoGar1$V4 
otoGar1$region <- "otoGar1"


calJac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg19_calJac3_valid-320.tab", sep="\t")
calJac3$length <- calJac3$V4 
calJac3$region <- "calJac3"


micMur1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg19_micMur1_valid-320.tab", sep="\t")
micMur1$length <- micMur1$V4 
micMur1$region <- "micMur1"

nomLeu1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg19_nomLeu1_valid-320.tab", sep="\t")
nomLeu1$length <- nomLeu1$V4 
nomLeu1$region <- "nomLeu1"


papHam1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg19_papHam1_valid-320.tab", sep="\t")
papHam1$length <- papHam1$V4 
papHam1$region <- "papHam1"


tupBel1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg19_tupBel1_valid-320.tab", sep="\t")
tupBel1$length <- tupBel1$V4 
tupBel1$region <- "tupBel1"


W_moderate <- rbind(hg38, panPan1, panTro4, rheMac3, tarSyr2, otoGar1, calJac3, micMur1, nomLeu1, papHam1, tupBel1)
W_moderate$alu_exon <- "Whole_moderate"



#### First 70 nt = U1  ####


hg38 <- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg19_hg38_valid-70.tab", sep="\t")
#TDP43 <- read.table(args[1], sep="\t")
hg38$length <- hg38$V4 
hg38$region <- "hg38"

panPan1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg38_panPan1_valid-70.tab", sep="\t")
panPan1$length <- panPan1$V4
panPan1$region <- "panPan1"

panTro4<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg38_panTro4_valid-70.tab", sep="\t")
panTro4$length <- panTro4$V4 
panTro4$region <- "panTro4"


rheMac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg38_rheMac3_valid-70.tab", sep="\t")
rheMac3$length <- rheMac3$V4 
rheMac3$region <- "rheMac3"


tarSyr2<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg38_tarSyr2_valid-70.tab", sep="\t")
tarSyr2$length <- tarSyr2$V4 
tarSyr2$region <- "tarSyr2"


otoGar1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg19_otoGar1_valid-70.tab", sep="\t")
otoGar1$length <- otoGar1$V4 
otoGar1$region <- "otoGar1"


calJac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg19_calJac3_valid-70.tab", sep="\t")
calJac3$length <- calJac3$V4 
calJac3$region <- "calJac3"


micMur1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg19_micMur1_valid-70.tab", sep="\t")
micMur1$length <- micMur1$V4 
micMur1$region <- "micMur1"

nomLeu1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg19_nomLeu1_valid-70.tab", sep="\t")
nomLeu1$length <- nomLeu1$V4 
nomLeu1$region <- "nomLeu1"


papHam1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg19_papHam1_valid-70.tab", sep="\t")
papHam1$length <- papHam1$V4 
papHam1$region <- "papHam1"


tupBel1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg19_tupBel1_valid-70.tab", sep="\t")
tupBel1$length <- tupBel1$V4 
tupBel1$region <- "tupBel1"


U1_moderate <- rbind(hg38, panPan1, panTro4, rheMac3, tarSyr2, otoGar1, calJac3, micMur1, nomLeu1, papHam1, tupBel1)
U1_moderate$alu_exon <- "U1_moderate"


######## last 250 nt U2 #########


hg38 <- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg19_hg38_valid-250.tab", sep="\t")
#TDP43 <- read.table(args[1], sep="\t")
hg38$length <- hg38$V4 
hg38$region <- "hg38"

panPan1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg38_panPan1_valid-250.tab", sep="\t")
panPan1$length <- panPan1$V4
panPan1$region <- "panPan1"

panTro4<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg38_panTro4_valid-250.tab", sep="\t")
panTro4$length <- panTro4$V4 
panTro4$region <- "panTro4"


rheMac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg38_rheMac3_valid-250.tab", sep="\t")
rheMac3$length <- rheMac3$V4 
rheMac3$region <- "rheMac3"


tarSyr2<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg38_tarSyr2_valid-250.tab", sep="\t")
tarSyr2$length <- tarSyr2$V4 
tarSyr2$region <- "tarSyr2"


otoGar1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg19_otoGar1_valid-250.tab", sep="\t")
otoGar1$length <- otoGar1$V4 
otoGar1$region <- "otoGar1"


calJac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg19_calJac3_valid-250.tab", sep="\t")
calJac3$length <- calJac3$V4 
calJac3$region <- "calJac3"


micMur1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg19_micMur1_valid-250.tab", sep="\t")
micMur1$length <- micMur1$V4 
micMur1$region <- "micMur1"

nomLeu1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg19_nomLeu1_valid-250.tab", sep="\t")
nomLeu1$length <- nomLeu1$V4 
nomLeu1$region <- "nomLeu1"


papHam1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg19_papHam1_valid-250.tab", sep="\t")
papHam1$length <- papHam1$V4 
papHam1$region <- "papHam1"


tupBel1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/moderate/All_Aluexons_3SS_Q_moderate_hg19_tupBel1_valid-250.tab", sep="\t")
tupBel1$length <- tupBel1$V4 
tupBel1$region <- "tupBel1"




U2_moderate <- rbind(hg38, panPan1, panTro4, rheMac3, tarSyr2, otoGar1, calJac3, micMur1, nomLeu1, papHam1, tupBel1)
U2_moderate$alu_exon <- "U2_moderate"




##################################################
## high
##################################################


##### Whole sequence 320

hg38 <- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg19_hg38_valid-320.tab", sep="\t")
#TDP43 <- read.table(args[1], sep="\t")
hg38$length <- hg38$V4 
hg38$region <- "hg38"

panPan1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg38_panPan1_valid-320.tab", sep="\t")
panPan1$length <- panPan1$V4
panPan1$region <- "panPan1"

panTro4<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg38_panTro4_valid-320.tab", sep="\t")
panTro4$length <- panTro4$V4 
panTro4$region <- "panTro4"


rheMac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg38_rheMac3_valid-320.tab", sep="\t")
rheMac3$length <- rheMac3$V4 
rheMac3$region <- "rheMac3"


tarSyr2<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg38_tarSyr2_valid-320.tab", sep="\t")
tarSyr2$length <- tarSyr2$V4 
tarSyr2$region <- "tarSyr2"


otoGar1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg19_otoGar1_valid-320.tab", sep="\t")
otoGar1$length <- otoGar1$V4 
otoGar1$region <- "otoGar1"


calJac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg19_calJac3_valid-320.tab", sep="\t")
calJac3$length <- calJac3$V4 
calJac3$region <- "calJac3"


micMur1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg19_micMur1_valid-320.tab", sep="\t")
micMur1$length <- micMur1$V4 
micMur1$region <- "micMur1"

nomLeu1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg19_nomLeu1_valid-320.tab", sep="\t")
nomLeu1$length <- nomLeu1$V4 
nomLeu1$region <- "nomLeu1"


papHam1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg19_papHam1_valid-320.tab", sep="\t")
papHam1$length <- papHam1$V4 
papHam1$region <- "papHam1"


tupBel1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg19_tupBel1_valid-320.tab", sep="\t")
tupBel1$length <- tupBel1$V4 
tupBel1$region <- "tupBel1"


W_high <- rbind(hg38, panPan1, panTro4, rheMac3, tarSyr2, otoGar1, calJac3, micMur1, nomLeu1, papHam1, tupBel1)
W_high$alu_exon <- "Whole_high"



#### First 70 nt = U1  ####


hg38 <- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg19_hg38_valid-70.tab", sep="\t")
#TDP43 <- read.table(args[1], sep="\t")
hg38$length <- hg38$V4 
hg38$region <- "hg38"

panPan1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg38_panPan1_valid-70.tab", sep="\t")
panPan1$length <- panPan1$V4
panPan1$region <- "panPan1"

panTro4<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg38_panTro4_valid-70.tab", sep="\t")
panTro4$length <- panTro4$V4 
panTro4$region <- "panTro4"


rheMac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg38_rheMac3_valid-70.tab", sep="\t")
rheMac3$length <- rheMac3$V4 
rheMac3$region <- "rheMac3"


tarSyr2<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg38_tarSyr2_valid-70.tab", sep="\t")
tarSyr2$length <- tarSyr2$V4 
tarSyr2$region <- "tarSyr2"


otoGar1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg19_otoGar1_valid-70.tab", sep="\t")
otoGar1$length <- otoGar1$V4 
otoGar1$region <- "otoGar1"


calJac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg19_calJac3_valid-70.tab", sep="\t")
calJac3$length <- calJac3$V4 
calJac3$region <- "calJac3"


micMur1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg19_micMur1_valid-70.tab", sep="\t")
micMur1$length <- micMur1$V4 
micMur1$region <- "micMur1"

nomLeu1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg19_nomLeu1_valid-70.tab", sep="\t")
nomLeu1$length <- nomLeu1$V4 
nomLeu1$region <- "nomLeu1"


papHam1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg19_papHam1_valid-70.tab", sep="\t")
papHam1$length <- papHam1$V4 
papHam1$region <- "papHam1"


tupBel1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg19_tupBel1_valid-70.tab", sep="\t")
tupBel1$length <- tupBel1$V4 
tupBel1$region <- "tupBel1"


U1_high <- rbind(hg38, panPan1, panTro4, rheMac3, tarSyr2, otoGar1, calJac3, micMur1, nomLeu1, papHam1, tupBel1)
U1_high$alu_exon <- "U1_high"


######## last 250 nt U2 #########


hg38 <- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg19_hg38_valid-250.tab", sep="\t")
#TDP43 <- read.table(args[1], sep="\t")
hg38$length <- hg38$V4 
hg38$region <- "hg38"

panPan1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg38_panPan1_valid-250.tab", sep="\t")
panPan1$length <- panPan1$V4
panPan1$region <- "panPan1"

panTro4<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg38_panTro4_valid-250.tab", sep="\t")
panTro4$length <- panTro4$V4 
panTro4$region <- "panTro4"


rheMac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg38_rheMac3_valid-250.tab", sep="\t")
rheMac3$length <- rheMac3$V4 
rheMac3$region <- "rheMac3"


tarSyr2<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg38_tarSyr2_valid-250.tab", sep="\t")
tarSyr2$length <- tarSyr2$V4 
tarSyr2$region <- "tarSyr2"


otoGar1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg19_otoGar1_valid-250.tab", sep="\t")
otoGar1$length <- otoGar1$V4 
otoGar1$region <- "otoGar1"


calJac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg19_calJac3_valid-250.tab", sep="\t")
calJac3$length <- calJac3$V4 
calJac3$region <- "calJac3"


micMur1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg19_micMur1_valid-250.tab", sep="\t")
micMur1$length <- micMur1$V4 
micMur1$region <- "micMur1"

nomLeu1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg19_nomLeu1_valid-250.tab", sep="\t")
nomLeu1$length <- nomLeu1$V4 
nomLeu1$region <- "nomLeu1"


papHam1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg19_papHam1_valid-250.tab", sep="\t")
papHam1$length <- papHam1$V4 
papHam1$region <- "papHam1"


tupBel1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/high/All_Aluexons_3SS_Q_high_hg19_tupBel1_valid-250.tab", sep="\t")
tupBel1$length <- tupBel1$V4 
tupBel1$region <- "tupBel1"




U2_high <- rbind(hg38, panPan1, panTro4, rheMac3, tarSyr2, otoGar1, calJac3, micMur1, nomLeu1, papHam1, tupBel1)
U2_high$alu_exon <- "U2_high"




##################################################
## highest
##################################################


##### Whole sequence 320

hg38 <- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg19_hg38_valid-320.tab", sep="\t")
#TDP43 <- read.table(args[1], sep="\t")
hg38$length <- hg38$V4 
hg38$region <- "hg38"

panPan1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg38_panPan1_valid-320.tab", sep="\t")
panPan1$length <- panPan1$V4
panPan1$region <- "panPan1"

panTro4<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg38_panTro4_valid-320.tab", sep="\t")
panTro4$length <- panTro4$V4 
panTro4$region <- "panTro4"


rheMac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg38_rheMac3_valid-320.tab", sep="\t")
rheMac3$length <- rheMac3$V4 
rheMac3$region <- "rheMac3"


tarSyr2<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg38_tarSyr2_valid-320.tab", sep="\t")
tarSyr2$length <- tarSyr2$V4 
tarSyr2$region <- "tarSyr2"


otoGar1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg19_otoGar1_valid-320.tab", sep="\t")
otoGar1$length <- otoGar1$V4 
otoGar1$region <- "otoGar1"


calJac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg19_calJac3_valid-320.tab", sep="\t")
calJac3$length <- calJac3$V4 
calJac3$region <- "calJac3"


micMur1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg19_micMur1_valid-320.tab", sep="\t")
micMur1$length <- micMur1$V4 
micMur1$region <- "micMur1"

nomLeu1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg19_nomLeu1_valid-320.tab", sep="\t")
nomLeu1$length <- nomLeu1$V4 
nomLeu1$region <- "nomLeu1"


papHam1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg19_papHam1_valid-320.tab", sep="\t")
papHam1$length <- papHam1$V4 
papHam1$region <- "papHam1"


tupBel1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg19_tupBel1_valid-320.tab", sep="\t")
tupBel1$length <- tupBel1$V4 
tupBel1$region <- "tupBel1"


W_highest <- rbind(hg38, panPan1, panTro4, rheMac3, tarSyr2, otoGar1, calJac3, micMur1, nomLeu1, papHam1, tupBel1)
W_highest$alu_exon <- "Whole_highest"



#### First 70 nt = U1  ####


hg38 <- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg19_hg38_valid-70.tab", sep="\t")
#TDP43 <- read.table(args[1], sep="\t")
hg38$length <- hg38$V4 
hg38$region <- "hg38"

panPan1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg38_panPan1_valid-70.tab", sep="\t")
panPan1$length <- panPan1$V4
panPan1$region <- "panPan1"

panTro4<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg38_panTro4_valid-70.tab", sep="\t")
panTro4$length <- panTro4$V4 
panTro4$region <- "panTro4"


rheMac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg38_rheMac3_valid-70.tab", sep="\t")
rheMac3$length <- rheMac3$V4 
rheMac3$region <- "rheMac3"


tarSyr2<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg38_tarSyr2_valid-70.tab", sep="\t")
tarSyr2$length <- tarSyr2$V4 
tarSyr2$region <- "tarSyr2"


otoGar1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg19_otoGar1_valid-70.tab", sep="\t")
otoGar1$length <- otoGar1$V4 
otoGar1$region <- "otoGar1"


calJac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg19_calJac3_valid-70.tab", sep="\t")
calJac3$length <- calJac3$V4 
calJac3$region <- "calJac3"


micMur1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg19_micMur1_valid-70.tab", sep="\t")
micMur1$length <- micMur1$V4 
micMur1$region <- "micMur1"

nomLeu1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg19_nomLeu1_valid-70.tab", sep="\t")
nomLeu1$length <- nomLeu1$V4 
nomLeu1$region <- "nomLeu1"


papHam1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg19_papHam1_valid-70.tab", sep="\t")
papHam1$length <- papHam1$V4 
papHam1$region <- "papHam1"


tupBel1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg19_tupBel1_valid-70.tab", sep="\t")
tupBel1$length <- tupBel1$V4 
tupBel1$region <- "tupBel1"


U1_highest <- rbind(hg38, panPan1, panTro4, rheMac3, tarSyr2, otoGar1, calJac3, micMur1, nomLeu1, papHam1, tupBel1)
U1_highest$alu_exon <- "U1_highest"


######## last 250 nt U2 #########


hg38 <- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg19_hg38_valid-250.tab", sep="\t")
#TDP43 <- read.table(args[1], sep="\t")
hg38$length <- hg38$V4 
hg38$region <- "hg38"

panPan1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg38_panPan1_valid-250.tab", sep="\t")
panPan1$length <- panPan1$V4
panPan1$region <- "panPan1"

panTro4<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg38_panTro4_valid-250.tab", sep="\t")
panTro4$length <- panTro4$V4 
panTro4$region <- "panTro4"


rheMac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg38_rheMac3_valid-250.tab", sep="\t")
rheMac3$length <- rheMac3$V4 
rheMac3$region <- "rheMac3"


tarSyr2<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg38_tarSyr2_valid-250.tab", sep="\t")
tarSyr2$length <- tarSyr2$V4 
tarSyr2$region <- "tarSyr2"


otoGar1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg19_otoGar1_valid-250.tab", sep="\t")
otoGar1$length <- otoGar1$V4 
otoGar1$region <- "otoGar1"


calJac3<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg19_calJac3_valid-250.tab", sep="\t")
calJac3$length <- calJac3$V4 
calJac3$region <- "calJac3"


micMur1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg19_micMur1_valid-250.tab", sep="\t")
micMur1$length <- micMur1$V4 
micMur1$region <- "micMur1"

nomLeu1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg19_nomLeu1_valid-250.tab", sep="\t")
nomLeu1$length <- nomLeu1$V4 
nomLeu1$region <- "nomLeu1"


papHam1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg19_papHam1_valid-250.tab", sep="\t")
papHam1$length <- papHam1$V4 
papHam1$region <- "papHam1"


tupBel1<- read.table("/media/igor/DATA/UCL/Evolution_Alus/New3SS/quantiles/highest/All_Aluexons_3SS_Q_highest_hg19_tupBel1_valid-250.tab", sep="\t")
tupBel1$length <- tupBel1$V4 
tupBel1$region <- "tupBel1"




U2_highest <- rbind(hg38, panPan1, panTro4, rheMac3, tarSyr2, otoGar1, calJac3, micMur1, nomLeu1, papHam1, tupBel1)
U2_highest$alu_exon <- "U2_highest"





##############################
####### final plots ##########
##############################

whole_table <- rbind(W_highest, U1_highest, U2_highest, W_high, U1_high, U2_high, W_moderate, U1_moderate, U2_moderate, W_low, U1_low, U2_low, W_lowest, U1_lowest, U2_lowest)
### quartile <- median.quartile(clusters$length)     Nejc what its clusters
quartile <- median.quartile(whole_table$length)


#pdf(paste(args[6],".pdf", sep=""))

a<-  ggplot(whole_table, aes(factor(region), length), alpha = 0.3, colour = "black") + 
      geom_boxplot() + 
      #scale_y_log10() + 
      theme_bw() +
      ggtitle("median quartile 25%, 50%, 75%") + 
      xlab("") + 
      ylab("length (nt)") +
      theme(text=element_text(size=12),axis.text=element_text(size=12), axis.title=element_text(size=12,face="plain")) +
      scale_x_discrete(limits=c("hg38", "panTro4", "panPan1", "nomLeu1", "papHam1", "rheMac3", "calJac3", "tarSyr2", "micMur1", "otoGar1", "tupBel1"))

a 



all_table <- data.frame(table(whole_table$region))
all_table 

library(plyr)
p_meds <- ddply(whole_table, .(region, alu_exon), summarise, med = mean(length))


dodge <- position_dodge(width = 0.9)

ALL<- ggplot(whole_table, aes(factor(region), length, colour = alu_exon), alpha = 1, width = 0.5) + 
  #geom_violin(position = dodge) + 
  geom_boxplot(width=1.1, position = dodge, alpha = 1, outlier.shape = NA ) +    ### uncoment if you need the boxplot inside
  theme_bw() +
  scale_y_continuous(limits = c(3, 33)) +
  #scale_y_log10() +
  ggtitle("U stretches") + 
  xlab("") + 
  ylab("U stretches count") +
  scale_colour_manual(values=c("#00441b", "#006d2c", "#238b45", "#41ae76", "#66c2a4", "#023858", "#045a8d", "#0570b0", "#3690c0", "#74a9cf", "#67001f", "#980043", "#ce1256", "#e7298a", "#df65b0" ))+
  theme(text=element_text(size=12),axis.text=element_text(size=12), axis.title=element_text(size=12,face="plain")) +
  geom_text(data = p_meds, aes(x = region, y = (med), label = format(med, digits=1)), size = 2.5, vjust = -1.5, position = dodge) +
  scale_x_discrete(limits=c("hg38", "panTro4", "panPan1", "nomLeu1", "papHam1", "rheMac3", "calJac3", "tarSyr2", "micMur1", "otoGar1", "tupBel1")) #+
  #scale_x_discrete(limits=c("W_highest", "U1_highest", "U2_highest"))
  

ALL

###################################
## Without whole
###################################


whole_table <- rbind( U1_highest, U2_highest, U1_high, U2_high, U1_moderate, U2_moderate, U1_low, U2_low, U1_lowest, U2_lowest)
### quartile <- median.quartile(clusters$length)     Nejc what its clusters
quartile <- median.quartile(whole_table$length)


#whole_table <- factor(c("U1_highest", "U2_highest", "U1_high", "U2_high", "U1_moderate", "U2_moderate", "U1_low", "U2_low", "U1_lowest", "U2_lowest"), levels=c("U1_highest", "U2_highest", "U1_high", "U2_high", "U1_moderate", "U2_moderate", "U1_low", "U2_low", "U1_lowest", "U2_lowest"))

#pdf(paste(args[6],".pdf", sep=""))

a<-  ggplot(whole_table, aes(factor(region), length), alpha = 0.3, colour = "black") + 
  geom_boxplot() + 
  #scale_y_log10() + 
  theme_bw() +
  ggtitle("median quartile 25%, 50%, 75%") + 
  xlab("") + 
  ylab("length (nt)") +
  theme(text=element_text(size=12),axis.text=element_text(size=12), axis.title=element_text(size=12,face="plain")) +
  scale_x_discrete(limits=c("hg38", "panTro4", "panPan1", "nomLeu1", "papHam1", "rheMac3", "calJac3", "tarSyr2", "micMur1", "otoGar1", "tupBel1"))

a 



all_table <- data.frame(table(whole_table$region))
all_table 

library(plyr)
p_meds <- ddply(whole_table, .(region, alu_exon), summarise, med = median(length))


dodge <- position_dodge(width = 0.9)


ALL<- ggplot(whole_table, aes(factor(region), length, colour = alu_exon), alpha = 1, width = 0.5) + 
  #geom_violin(position = dodge) + 
  geom_boxplot(width=1.1, position = dodge, alpha = 1, outlier.shape = NA ) +    ### uncoment if you need the boxplot inside
  theme_bw() +
  scale_y_continuous(limits = c(3, 33)) +
  #scale_y_log10() +
  ggtitle("U stretches") + 
  xlab("") + 
  ylab("U stretches count") +
  #scale_colour_manual(values=c("#00441b", "#0570b0", "#810f7c"))+
  scale_colour_manual(values=c("#00441b", "#006d2c", "#238b45", "#41ae76", "#66c2a4", "#023858", "#045a8d", "#0570b0", "#3690c0", "#74a9cf"))+
  #scale_colour_manual(values=c("#242B38","#2165E8","#2CA3B5")) +
  theme(text=element_text(size=12),axis.text=element_text(size=12), axis.title=element_text(size=12,face="plain")) +
  geom_text(data = p_meds, aes(x = region, y = (med), label = format(med, digits=1)), size = 2.5, vjust = -1.5, position = dodge) +
  scale_x_discrete(limits=c("hg38", "panTro4", "panPan1", "nomLeu1", "papHam1", "rheMac3", "calJac3", "tarSyr2", "micMur1", "otoGar1", "tupBel1"))

ALL



####### Data count #####

all_table <- data.frame(table(whole_table$region, whole_table$alu_exon))
all_table 
#all_table[1,2]

