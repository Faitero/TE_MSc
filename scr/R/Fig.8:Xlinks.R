### First use /media/igor/DATA/UCL/Evolution_Alus/New3SS/CLIPs_on_evolution_groups/get_tables_for CLIP.sh to create the tables of xlinks on aluexons

library(ggplot2)
require(ggplot2)
require(reshape)
library(system)
library(boot)
library(parallel)
#library(doParallel)
#registerDoParallel(cores=detectCores(all.tests=TRUE))

if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}
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
if (!require("GGally")) {
  install.packages("GGally", dependencies = TRUE)
  library(GGally)
}


### Function

PolyU_clasification <- function(data_frame) {

  count <- nrow(data_frame)
  less7 <- subset(data_frame, data_frame$"WU_hg19" < 7)
  less7$Utrack <- "less7"
  betw7.8 <- subset(data_frame, data_frame$"WU_hg19" == 7 | data_frame$"WU_hg19" == 8)
  betw7.8$Utrack <- "betw7.8"
  betw9.10 <- subset(data_frame, data_frame$"WU_hg19" == 9 | data_frame$"WU_hg19" == 10)
  betw9.10$Utrack <- "betw9.10"
  more10 <- subset(data_frame, data_frame$"WU_hg19" > 10)
  more10$Utrack <- "more10"

  result <- rbind(less7, betw7.8, betw9.10, more10)
  count2 <- nrow(result)
  count == count2
  #names(result) <- c("alu_element_hg19_t", "alu_substitutions_hg19_t", "X3SSS_hg19", "WU_hg19", "Exon_Originated", "CLIPS", "Utrack")


  return(data.frame(result))
}




setwd('/Users/Igor/Dropbox (UCL-MN Team)/hnRNPC.NMD (1).manuscript/Alus_primates_evolution/Data')
setwd('/home/igor/Dropbox (UCL-MN Team)/hnRNPC.NMD (1).manuscript/Alus_primates_evolution/Data')  ## Curro

whole_final_5sp <- read.table("whole_final_20160620_5species.csv", sep="\t", header=TRUE)
whole_final_5sp$aluexon <- as.character(whole_final_5sp[,1])

setwd('/media/igor/DATA/UCL/Evolution_Alus/New3SS/CLIPs_on_evolution_groups/New_Proteins')

total_xlinks <- read.table("total_xlinks_aluexons.tab", sep="\t", header=FALSE)

# paste ../last_tables/total_hnRNPC4964_G.bed ../last_tables/total_PTB4174_G.bed ../last_tables/total_U2AF65_PTB_kd_G.bed ../last_tables/total_U2AF65_ctrl_G.bed ../last_tables/total_U2AF65_wt_G.bed ../last_tables/total_U2AF65_hnRNPC_kd_G.bed total_HUR_4361_G.bed total_TIA1_3959_G.bed total_TIAL1_3936_G.bed > total_xlinks_aluexons.tab



### Merge CLIP table with previous whole table
temp <- total_xlinks[,c(4,7,37, 47,57,67)] #17 ,77,87)] Remove PTB, U2AF65_ctrl, TIA1 and TIAL1
names(temp) <- c("aluexon", "hnRNPC",  "U2AF65_ctrl", "U2AF65_wt", "U2AF65_hnRNPC_kd", "HUR")
temp$aluexon <- as.character(temp[,1])

whole_final_5sp_CLIPs <- merge.data.frame(whole_final_5sp, temp, by="aluexon", all.x=TRUE)


plots_xlinks_mini <- whole_final_5sp_CLIPs[,c(11,12,20,21,70,73:78)]
head(plots_xlinks_mini)
nrow(plots_xlinks_mini)
rownames(plots_xlinks_mini) <- whole_final_5sp_CLIPs[,1]


## Clasify in base of U-track lenght

plots_xlinks_hnRNPC<- plots_xlinks_mini[,c(1:6,7)]
head(plots_xlinks_hnRNPC)
plots_xlinks_hnRNPC <- PolyU_clasification(plots_xlinks_hnRNPC)
names(plots_xlinks_hnRNPC) <- c("alu_element_hg19_t", "alu_substitutions_hg19_t", "X3SSS_hg19", "WU_hg19", "Furthest", "Exon_Originated", "CLIPS", "Utrack")
plots_xlinks_hnRNPC$protein <- "hnRNPC"
head(plots_xlinks_hnRNPC)

## We don't use U2AF65_ctrl

plots_xlinks_U2AF65_wt<- plots_xlinks_mini[,c(1:6,9)]
head(plots_xlinks_U2AF65_wt)
plots_xlinks_U2AF65_wt <- PolyU_clasification(plots_xlinks_U2AF65_wt)
names(plots_xlinks_U2AF65_wt) <- c("alu_element_hg19_t", "alu_substitutions_hg19_t", "X3SSS_hg19", "WU_hg19", "Furthest", "Exon_Originated", "CLIPS", "Utrack")
plots_xlinks_U2AF65_wt$protein <- "U2AF65_wt"
head(plots_xlinks_U2AF65_wt)

plots_xlinks_U2AF65_hnRNPC_kd<- plots_xlinks_mini[,c(1:6,10)]
head(plots_xlinks_U2AF65_hnRNPC_kd)
plots_xlinks_U2AF65_hnRNPC_kd <- PolyU_clasification(plots_xlinks_U2AF65_hnRNPC_kd)
names(plots_xlinks_U2AF65_hnRNPC_kd) <- c("alu_element_hg19_t", "alu_substitutions_hg19_t", "X3SSS_hg19", "WU_hg19", "Furthest", "Exon_Originated", "CLIPS", "Utrack")
plots_xlinks_U2AF65_hnRNPC_kd$protein <- "U2AF65_hnRNPC_kd"
head(plots_xlinks_U2AF65_hnRNPC_kd)

plots_xlinks_HUR<- plots_xlinks_mini[,c(1:6,11)]
head(plots_xlinks_HUR)
plots_xlinks_HUR <- PolyU_clasification(plots_xlinks_HUR)
names(plots_xlinks_HUR) <- c("alu_element_hg19_t", "alu_substitutions_hg19_t", "X3SSS_hg19", "WU_hg19", "Furthest", "Exon_Originated", "CLIPS", "Utrack")
plots_xlinks_HUR$protein <- "HUR"
head(plots_xlinks_HUR)


final_xlinks <- rbind(plots_xlinks_hnRNPC, plots_xlinks_U2AF65_wt, plots_xlinks_U2AF65_hnRNPC_kd, plots_xlinks_HUR)

levels(final_xlinks$Utrack) <- as.factor(c("less7", "bet7.8", "bet9.10", "more10"))

factor(final_xlinks$Utrack)
nrow(final_xlinks)

final_xlinks$Exon_Originated <- as.factor(final_xlinks$Exon_Originated)
final_xlinks$alu_element_hg19_t <- as.factor(final_xlinks$alu_element_hg19_t)
final_xlinks$alu_substitutions_hg19_t <- as.factor(final_xlinks$alu_substitutions_hg19_t)
final_xlinks$Furthest <- as.character(final_xlinks$Furthest)
final_xlinks$CLIPS <- as.numeric(final_xlinks$CLIPS)
head(final_xlinks)
nrow(final_xlinks)

### final_xlinks final long data frame to plot ... now filtering



final_xlinks_filter <- subset(final_xlinks, !is.na(final_xlinks$Exon_Originated) & final_xlinks$CLIPS > 3 & Furthest == "calJac3") #
nrow(final_xlinks_filter)
head(final_xlinks_filter)



## row bind and remove duplicates
# whole_final_5sp_CLIPs_Top <- rbind(whole_final_5sp_CLIPs_Top_hnRNPC, whole_final_5sp_CLIPs_Top_U2AF65_ctrl, whole_final_5sp_CLIPs_Top_U2AF65_wt, whole_final_5sp_CLIPs_Top_U2AF65_hnRNPC_kd, whole_final_5sp_CLIPs_Top_HUR) # , whole_final_5sp_CLIPs_Top_TIA1, whole_final_5sp_CLIPs_Top_TIAL1)
# duplicated(whole_final_5sp_CLIPs_Top)
# whole_final_5sp_CLIPs_Top <- whole_final_5sp_CLIPs_Top[!duplicated(whole_final_5sp_CLIPs_Top),]
#


### Now ussing 3ss >5 not the top 30% !!!!!!! take care
whole_final_5sp_CLIPs_Top <- whole_final_5sp_CLIPshigh5
##!!!!!!!
plots_xlinks <- whole_final_5sp_CLIPs_Top

## subset 3ss >3 and MArmoset originated
whole_final_5sp_CLIPs_Top_marmoset_high3sshuman <- subset(whole_final_5sp_CLIPs_Top, whole_final_5sp_CLIPs_Top$X3SSS_hg19 >=3 & whole_final_5sp_CLIPs_Top$furthest == "calJac3") # 1190 cool
plots_xlinks <- whole_final_5sp_CLIPs_Top_marmoset_high3sshuman
nrow(plots_xlinks)

# ## OR subset 3ss >3 all species
# whole_final_5sp_CLIPs_Top_high3sshuman <- subset(whole_final_5sp_CLIPs_Top, whole_final_5sp_CLIPs_Top$X3SSS_hg19 >=3) # 2295 cool
# plots_xlinks <- whole_final_5sp_CLIPs_Top_high3sshuman
# nrow(whole_final_5sp_CLIPs_Top_high3sshuman)




#Mean
#p_meds <- ddply(final_xlinks_filter, .(protein, Utrack), summarise, med = mean(CLIPS, na.rm = TRUE))
#
# getmode <- function(v) {
#   uniqv <- unique(v)
#   uniqv[which.max(tabulate(match(v, uniqv)))]
# }

#p_meds <- ddply(final_xlinks_filter, .(protein, Exon_Originated), summarise, med = getmode(CLIPS)
#Median
p_meds <- ddply(final_xlinks_filter, .(protein, Exon_Originated), summarise, med = median(CLIPS, na.rm = TRUE))
p_meds


count_table <- data.frame(table(final_xlinks_filter$protein, final_xlinks_filter$Exon_Originated))
colnames(count_table) <- c("protein", "Exon_Originated", "freq")
count_table


library(plyr)

dodge <- position_dodge(width = 0.9)

plotXlinks<- ggplot(final_xlinks_filter, aes(factor(protein), CLIPS, fill = Exon_Originated), alpha = 1, width = 0.5) +
  geom_violin(position = dodge, alpha = 0.6) +
  geom_boxplot(width=0.3, position = dodge, alpha = 1, outlier.shape = NA) +    ### uncoment if you need the boxplot inside
  theme_bw() +
  #theme(panel.background = element_rect(fill='white', colour='black')) +
  scale_y_continuous(limits = c(-1 , 250)) +
  #scale_y_log10() +
  ggtitle("Xlinks Utrack classified") +
  xlab("") +
  ylab("Xlinks") +
  #scale_colour_manual(values=c("#00441b", "#006d2c", "#238b45", "#41ae76", "#66c2a4", "#023858", "#045a8d", "#0570b0", "#3690c0", "#74a9cf", "#67001f", "#980043", "#ce1256", "#e7298a", "#df65b0")) +                                             #"#00441b", "#006d2c", "#0570b0", "#3690c0", "#810f7c", "#88419d", "#8c6bb1"))+  #("#00441b", "#006d2c", "#238b45", "#0570b0", "#3690c0", "#74a9cf", "#810f7c", "#88419d", "#8c6bb1")
  #scale_fill_manual(values=c("#00441b", "#006d2c", "#238b45", "#41ae76", "#66c2a4", "#023858", "#045a8d", "#0570b0", "#3690c0", "#74a9cf", "#67001f", "#980043", "#ce1256", "#e7298a", "#df65b0")) +
  scale_fill_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
  theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold")) +
  geom_text(data = p_meds, aes(x = protein, y = med, label = format((med+1), digits=2)), size = 3, vjust = -1.5, position = dodge, color="white") +
  geom_text(data = count_table, aes(x = protein, y = -1 , label = freq), size = 3, vjust = 0, position = dodge, angle=45) +
  scale_x_discrete(limits=c("hnRNPC", "U2AF65_wt", "U2AF65_hnRNPC_kd", "HUR")) # , labels=c("hg19"="Human", "panTro4"="Chimp", "nomLeu1"="Gibbon", "rheMac3"="Rhesus", "calJac3"="Marmoset"))


plotXlinks

ggsave("CLIPS_aluexons_Originated_high3clips-high3hum_marmoset.pdf", width=20, height=13)

getwd()



plots_xlinks_long <- rbind(plots_xlinks_hnRNPC, plots_xlinks_U2AF65_wt, plots_xlinks_U2AF65_hnRNPC_kd, plots_xlinks_HUR, plots_xlinks_TIA1, plots_xlinks_TIAL1)


dodge <- position_dodge(width = 0.9)
#Mean
#p_meds <- ddply(final_data_frame, .(region, group), summarise, med = mean(X3SSS, na.rm = TRUE))
#Mode
#p_meds <- ddply(final_data_frame, .(region, group), summarise, med = getmode(X3SSS, na.rm = TRUE))
#Median
p_meds <- ddply(plots_xlinks_long, .(Exon_Originated, protein), summarise, med = median(xlinks, na.rm = TRUE))
p_meds

a<-  ggplot(plots_xlinks_long, aes(factor(protein), xlinks), alpha = 1, colour = "black") +
  geom_boxplot() +
  #scale_y_log10() +
  theme_bw() +
  ggtitle("median quartile 25%, 50%, 75%") +
  xlab("") +
  ylab("length (nt)") +
  theme(text=element_text(size=12),axis.text=element_text(size=12), axis.title=element_text(size=12,face="plain")) +
  scale_x_discrete(limits=c("hnRNPC", "U2AF65_wt", "U2AF65_hnRNPC_kd", "HUR", "TIA1", "TIAL1"))

a

count_table <- data.frame(table(plots_xlinks_long$protein, plots_xlinks_long$Exon_Originated))
colnames(count_table) <- c("protein", "Exon_Originated", "freq")
count_table


library(plyr)


dodge <- position_dodge(width =1)

pdf_name <- paste0("3SS_" , experiment_name, ".pdf")
#pdf(pdf_name , width=20, height=13)
ggtitle_name <- paste0("3SS ", experiment_name)

plotxlinks<- ggplot(plots_xlinks_long, aes(factor(protein), xlinks, fill = Exon_Originated), alpha = 1, width = 0.5) +
  geom_violin(position = dodge, alpha = 0.6) +
  geom_boxplot(width=0.5, position = dodge, alpha = 1, outlier.shape = NA) +    ### uncoment if you need the boxplot inside
  theme_bw() +
  #theme(panel.background = element_rect(fill='white', colour='black')) +
  scale_y_continuous(limits = c(-2, 100)) +
  #scale_y_log10() +
  #ggtitle(ggtitle_name) +
  xlab("") +
  ylab("XLinks Counts") +
  #scale_colour_manual(values=c("#00441b", "#006d2c", "#238b45", "#41ae76", "#66c2a4", "#023858", "#045a8d", "#0570b0", "#3690c0", "#74a9cf", "#67001f", "#980043", "#ce1256", "#e7298a", "#df65b0")) +                                             #"#00441b", "#006d2c", "#0570b0", "#3690c0", "#810f7c", "#88419d", "#8c6bb1"))+  #("#00441b", "#006d2c", "#238b45", "#0570b0", "#3690c0", "#74a9cf", "#810f7c", "#88419d", "#8c6bb1")
  scale_fill_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
  theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold")) +
  geom_text(data = p_meds, aes(x = Exon_Originated, y = med, label = format((med+1), digits=2)), size = 3, vjust = -1.5, position = dodge, color="white") +
  geom_text(data = count_table, aes(x = protein, y = -1.5 , label = freq), size = 3, vjust = 0, position = dodge, angle=45) +
  scale_x_discrete(limits=c("hnRNPC", "U2AF65_wt", "U2AF65_hnRNPC_kd", "HUR", "TIA1", "TIAL1"))


plotxlinks

# ggsave("CLIPS_aluexons_more30_xlinks.pdf", width=20, height=13)



p_meds <- ddply(plots_xlinks_long_fract, .(Exon_Originated, protein), summarise, med = median(xlinks, na.rm = TRUE))
p_meds


plotxlinks<- ggplot(plots_xlinks_long_fract, aes(factor(protein), xlinks, fill = Exon_Originated), alpha = 1, width = 0.5) +
  geom_violin(position = dodge, alpha = 0.6) +
  geom_boxplot(width=0.3, position = dodge, alpha = 1, outlier.shape = NA) +    ### uncoment if you need the boxplot inside
  theme_bw() +
  #theme(panel.background = element_rect(fill='white', colour='black')) +
  scale_y_continuous(limits = c(-1, 15)) +
  #scale_y_log10() +
  ggtitle("Fold Change") +
  xlab("") +
  ylab("Fold Change") +
  #scale_colour_manual(values=c("#00441b", "#006d2c", "#238b45", "#41ae76", "#66c2a4", "#023858", "#045a8d", "#0570b0", "#3690c0", "#74a9cf", "#67001f", "#980043", "#ce1256", "#e7298a", "#df65b0")) +                                             #"#00441b", "#006d2c", "#0570b0", "#3690c0", "#810f7c", "#88419d", "#8c6bb1"))+  #("#00441b", "#006d2c", "#238b45", "#0570b0", "#3690c0", "#74a9cf", "#810f7c", "#88419d", "#8c6bb1")
  scale_fill_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
  theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold")) +
  #geom_text(data = p_meds, aes(x = Exon_Originated, y = med, label = format((med+1), digits=2)), size = 3, vjust = -1.5, position = dodge, color="white") +
  #geom_text(data = count_table, aes(x = region, y = 15 , label = freq), size = 3, vjust = 0, position = dodge, angle=45) +
  scale_x_discrete(limits=c("hnRNPC_frac_U2AF65", "U2AF65_frac_U2AF65_hnRNPCKD"), labels=c("hnRNPC_frac_U2AF65"="hnRNPC/U2AF65", "U2AF65_frac_U2AF65_hnRNPCKD"="U2AF65/U2AF65_hnRNPCKD"))
  # scale_x_discrete(limits=c("hg19", "panTro4", "panPan1", "nomLeu1", "papHam1", "rheMac3", "calJac3"), labels=c("hg19"="Human", "panTro4"="Chimp", "panPan1"="Bonobo", "nomLeu1"="Gibbon", "papHam1"="Baboon", "rheMac3"="Rhesus", "calJac3"="Marmoset"))


plotxlinks

ggsave("FC_alus_no_filtering.pdf", width=20, height=13)




ggsave(pdf_name, width=20, height=13)




rownames(plots_xlinks_mini) <- plots_xlinks_mini[,1]
plots_xlinks_mini$Exon_Originated <- as.factor(plots_xlinks_mini$Exon_Originated)
head(plots_xlinks_mini)

# Boot straping on xlinks
plots_xlinks_long <- rbind(plots_xlinks_hnRNPC, plots_xlinks_U2AF65_wt, plots_xlinks_U2AF65_hnRNPC_kd)
# 2721 obs!!! OK
tot <- data.frame(stringsAsFactors=FALSE)
#names(tot) <- c("Median", "Normal_low", "Normal_hight", "Basic_low", "Basic_hight", "Region", "Group")

tot <- data.frame(Median=as.numeric(),
                  Normal_low=as.numeric(),
                  Normal_hight=as.numeric(),
                  Basic_low=as.numeric(),
                  Basic_hight=as.numeric(),
                  Region=as.character(),
                  Group=as.character(),
                  stringsAsFactors=FALSE)

lista <- c("hnRNPC", "hnRNPC", "U2AF65_wt", "U2AF65_hnRNPC_kd")

head(plots_xlinks_long)
unique(plots_xlinks_long$Exon_Originated)
unique(plots_xlinks_long$protein)

########## For 3SS
for (variable in lista){

  print(variable)
  toboot <- plots_xlinks_long[plots_xlinks_long$protein == variable,]

  aggregate(toboot$xlinks, by=list(var=toboot$Exon_Originated), median)
  boot.statistics = function(data, indices){
    d <- data[ indices, ]
    d <- aggregate(d$xlinks, by=list(var=d$Exon_Originated), median)
    return(d$x)
  }

  #toboot$group


  bootobject1 <- boot(data=toboot, statistic = boot.statistics, R=20000,  parallel = "multicore", ncpus=4)

  plot(bootobject1, index=3)
  print(bootobject1)
  plot(bootobject1)

  Constant <- boot.ci(bootobject1, conf=0.95, type = c("norm", "basic"), index=1) ## Constant Exon
  Emerging <- boot.ci(bootobject1, conf=0.95, type = c("norm", "basic"), index=2) ## Emerging Exon
  Evolving <- boot.ci(bootobject1, conf=0.95, type = c("norm", "basic"), index=3) ## Evolving Exon

  line_Constant  <- c(bootobject1$t0[1], Constant$normal[2], Constant$normal[3], Constant$basic[4], Constant$basic[5], variable, "Constant")
  line_Emerging <- c(bootobject1$t0[2], Emerging$normal[2], Emerging$normal[3], Emerging$basic[4], Emerging$basic[5], variable, "Emerging")
  line_Evolving<- c(bootobject1$t0[3], Evolving$normal[2], Evolving$normal[3], Evolving$basic[4], Evolving$basic[5], variable, "Evolving")



  tot[,1:5] <- sapply(tot[,1:5], as.numeric)
  tot[,6:7] <- sapply(tot[,6:7], as.character)
  tot <- rbind(tot, line_Constant, line_Emerging, line_Evolving)
  names(tot) <- c("Median", "Normal_low", "Normal_hight", "Basic_low", "Basic_hight", "Region", "Group")

}

tot <- tot[-c(1:3),]

setwd('/home/igor/Dropbox (UCL-MN Team)/hnRNPC.NMD (1).manuscript/Alus_primates_evolution/Data')  ## Curro
write.table(tot, file='Bootstraping_xlinks_evolution_alus.tab', sep="\t", row.names = FALSE, col.names = TRUE)



# Multiplot pairwaise Xlinks comparision
library(ggplot2)
#data(iris)
#ggpairs(plots_xlinks_mini, colour='Exon_Originated', alpha=0.4)
tiff("RBP_aluexons.tiff" , width = 1000, height = 1000, units = "px")
GGally::ggpairs(plots_xlinks_mini, aes(colour = Exon_Originated, alpha=0.4))
dev.off()


#  Multiplot custom pairwaise Xlinks comparision
library(gridExtra)

A <- ggplot(plots_xlinks_mini, aes(plots_xlinks_mini$hnRNPC, plots_xlinks$U2AF65_ctrl)) +
  geom_point(aes(color=Exon_Originated, alpha=1/100)) +
  theme_bw() +
  scale_y_continuous(limits = c(0 , 1000)) +
  scale_x_continuous(limits = c(0 , 1000)) +
  #ggtitle("Aluexons binding") +
  xlab("hnRNP C") +
  ylab("U2AF65_ctrl") +
  scale_colour_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
  theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold"))


  B <- ggplot(plots_xlinks_mini, aes(plots_xlinks_mini$hnRNPC, plots_xlinks$U2AF65_wt)) +
    geom_point(aes(color=Exon_Originated, alpha=1/100)) +
    theme_bw() +
    scale_y_continuous(limits = c(0 , 1000)) +
    scale_x_continuous(limits = c(0 , 1000)) +
    #ggtitle("Aluexons binding") +
    xlab("hnRNPC") +
    ylab("U2AF65_wt") +
    scale_colour_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold"))




  C <- ggplot(plots_xlinks_mini, aes(plots_xlinks_mini$hnRNPC, plots_xlinks$U2AF65_hnRNPC_kd)) +
    geom_point(aes(color=Exon_Originated, alpha=1/100)) +
    theme_bw() +
    scale_y_continuous(limits = c(0 , 1000)) +
    scale_x_continuous(limits = c(0 , 1000)) +
    #ggtitle("Aluexons binding") +
    xlab("hnRNPC") +
    ylab("U2AF65_hnRNPC_kd") +
    scale_colour_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold"))



  D <- ggplot(plots_xlinks_mini, aes(plots_xlinks_mini$hnRNPC, plots_xlinks$U2AF65_PTB_kd)) +
    geom_point(aes(color=Exon_Originated, alpha=1/100)) +
    theme_bw() +
    scale_y_continuous(limits = c(0 , 1000)) +
    scale_x_continuous(limits = c(0 , 1000)) +
    #ggtitle("Aluexons binding") +
    xlab("hnRNPC") +
    ylab("U2AF65_PTB_kd") +
    scale_colour_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold"))



  E <- ggplot(plots_xlinks_mini, aes(plots_xlinks_mini$hnRNPC, plots_xlinks$PTB)) +
    geom_point(aes(color=Exon_Originated, alpha=1/100)) +
    theme_bw() +
    scale_y_continuous(limits = c(0 , 1000)) +
    scale_x_continuous(limits = c(0 , 1000)) +
    #ggtitle("Aluexons binding") +
    xlab("hnRNPC") +
    ylab("PTB") +
    scale_colour_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold"))




  F <- ggplot(plots_xlinks_mini, aes(plots_xlinks_mini$U2AF65_ctrl, plots_xlinks$U2AF65_wt)) +
    geom_point(aes(color=Exon_Originated, alpha=1/100)) +
    theme_bw() +
    scale_y_continuous(limits = c(0 , 1000)) +
    scale_x_continuous(limits = c(0 , 1000)) +
    #ggtitle("Aluexons binding") +
    xlab("U2AF65_ctrl") +
    ylab("U2AF65_wt") +
    scale_colour_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold"))



  G <- ggplot(plots_xlinks_mini, aes(plots_xlinks_mini$U2AF65_ctrl, plots_xlinks$U2AF65_hnRNPC_kd)) +
    geom_point(aes(color=Exon_Originated, alpha=1/100)) +
    theme_bw() +
    scale_y_continuous(limits = c(0 , 1000)) +
    scale_x_continuous(limits = c(0 , 1000)) +
    #ggtitle("Aluexons binding") +
    xlab("U2AF65_ctrl") +
    ylab("U2AF65_hnRNPC_kd") +
    scale_colour_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold"))



  H <- ggplot(plots_xlinks_mini, aes(plots_xlinks_mini$U2AF65_ctrl, plots_xlinks$U2AF65_PTB_kd)) +
    geom_point(aes(color=Exon_Originated, alpha=1/100)) +
    theme_bw() +
    scale_y_continuous(limits = c(0 , 1000)) +
    scale_x_continuous(limits = c(0 , 1000)) +
    #ggtitle("Aluexons binding") +
    xlab("U2AF65_ctrl") +
    ylab("U2AF65_PTB_kd") +
    scale_colour_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold"))



  I <- ggplot(plots_xlinks_mini, aes(plots_xlinks_mini$U2AF65_ctrl, plots_xlinks$PTB)) +
    geom_point(aes(color=Exon_Originated, alpha=1/100)) +
    theme_bw() +
    scale_y_continuous(limits = c(0 , 1000)) +
    scale_x_continuous(limits = c(0 , 1000)) +
    #ggtitle("Aluexons binding") +
    xlab("U2AF65_ctrl") +
    ylab("PTB") +
    scale_colour_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold"))



  J <- ggplot(plots_xlinks_mini, aes(plots_xlinks_mini$U2AF65_wt, plots_xlinks$U2AF65_hnRNPC_kd)) +
    geom_point(aes(color=Exon_Originated, alpha=1/100)) +
    theme_bw() +
    scale_y_continuous(limits = c(0 , 1000)) +
    scale_x_continuous(limits = c(0 , 1000)) +
    #ggtitle("Aluexons binding") +
    xlab("U2AF65_wt") +
    ylab("U2AF65_hnRNPC_kd") +
    scale_colour_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold"))



  K <- ggplot(plots_xlinks_mini, aes(plots_xlinks_mini$U2AF65_wt, plots_xlinks$U2AF65_PTB_kd)) +
    geom_point(aes(color=Exon_Originated, alpha=1/100)) +
    theme_bw() +
    scale_y_continuous(limits = c(0 , 1000)) +
    scale_x_continuous(limits = c(0 , 1000)) +
    #ggtitle("Aluexons binding") +
    xlab("U2AF65_wt") +
    ylab("U2AF65_PTB_kd") +
    scale_colour_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold"))



  L <- ggplot(plots_xlinks_mini, aes(plots_xlinks_mini$U2AF65_wt, plots_xlinks$PTB)) +
    geom_point(aes(color=Exon_Originated, alpha=1/100)) +
    theme_bw() +
    scale_y_continuous(limits = c(0 , 1000)) +
    scale_x_continuous(limits = c(0 , 1000)) +
    #ggtitle("Aluexons binding") +
    xlab("U2AF65_wt") +
    ylab("PTB") +
    scale_colour_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold"))




  M <- ggplot(plots_xlinks_mini, aes(plots_xlinks_mini$U2AF65_hnRNPC_kd, plots_xlinks$U2AF65_PTB_kd)) +
    geom_point(aes(color=Exon_Originated, alpha=1/100)) +
    theme_bw() +
    scale_y_continuous(limits = c(0 , 1000)) +
    scale_x_continuous(limits = c(0 , 1000)) +
    #ggtitle("Aluexons binding") +
    xlab("U2AF65_hnRNPC_kd") +
    ylab("U2AF65_PTB_kd") +
    scale_colour_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold"))



  N <- ggplot(plots_xlinks_mini, aes(plots_xlinks_mini$U2AF65_hnRNPC_kd, plots_xlinks$PTB)) +
    geom_point(aes(color=Exon_Originated, alpha=1/100)) +
    theme_bw() +
    scale_y_continuous(limits = c(0 , 1000)) +
    scale_x_continuous(limits = c(0 , 1000)) +
    #ggtitle("Aluexons binding") +
    xlab("U2AF65_hnRNPC_kd") +
    ylab("PTB") +
    scale_colour_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold"))



  O <- ggplot(plots_xlinks_mini, aes(plots_xlinks_mini$U2AF65_PTB_kd, plots_xlinks$PTB)) +
    geom_point(aes(color=Exon_Originated, alpha=1/100)) +
    theme_bw() +
    scale_y_continuous(limits = c(0 , 1000)) +
    scale_x_continuous(limits = c(0 , 1000)) +
    #ggtitle("Aluexons binding") +
    xlab("U2AF65_PTB_kd") +
    ylab("PTB") +
    scale_colour_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold"))


  P <- ggplot(plots_xlinks_mini, aes(plots_xlinks_mini$hnRNPC, plots_xlinks$HUR)) +
    geom_point(aes(color=Exon_Originated, alpha=1/100)) +
    theme_bw() +
    scale_y_continuous(limits = c(0 , 1000)) +
    scale_x_continuous(limits = c(0 , 1000)) +
    #ggtitle("Aluexons binding") +
    xlab("hnRNPC") +
    ylab("HUR") +
    scale_colour_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold"))
  P


  Q<-  ggplot(plots_xlinks_mini, aes(plots_xlinks_mini$U2AF65_ctrl, plots_xlinks$HUR)) +
    geom_point(aes(color=Exon_Originated, alpha=1/100)) +
    theme_bw() +
    scale_y_continuous(limits = c(0 , 1000)) +
    scale_x_continuous(limits = c(0 , 1000)) +
    #ggtitle("Aluexons binding") +
    xlab("U2AF65_ctrl") +
    ylab("HUR") +
    scale_colour_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold"))

  R<-ggplot(plots_xlinks_mini, aes(plots_xlinks_mini$U2AF65_wt, plots_xlinks$HUR)) +
    geom_point(aes(color=Exon_Originated, alpha=1/100)) +
    theme_bw() +
    scale_y_continuous(limits = c(0 , 1000)) +
    scale_x_continuous(limits = c(0 , 1000)) +
    #ggtitle("Aluexons binding") +
    xlab("U2AF65_wt") +
    ylab("HUR") +
    scale_colour_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold"))




  S<-ggplot(plots_xlinks_mini, aes(plots_xlinks_mini$U2AF65_hnRNPC_kd, plots_xlinks$HUR)) +
    geom_point(aes(color=Exon_Originated, alpha=1/100)) +
    theme_bw() +
    scale_y_continuous(limits = c(0 , 1000)) +
    scale_x_continuous(limits = c(0 , 1000)) +
    #ggtitle("Aluexons binding") +
    xlab("U2AF65_hnRNPC_kd") +
    ylab("HUR") +
    scale_colour_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold"))



  T<-ggplot(plots_xlinks_mini, aes(plots_xlinks_mini$U2AF65_PTB_kd, plots_xlinks$HUR)) +
    geom_point(aes(color=Exon_Originated, alpha=1/100)) +
    theme_bw() +
    scale_y_continuous(limits = c(0 , 1000)) +
    scale_x_continuous(limits = c(0 , 1000)) +
    #ggtitle("Aluexons binding") +
    xlab("U2AF65_PTB_kd") +
    ylab("HUR") +
    scale_colour_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold"))

  U<-ggplot(plots_xlinks_mini, aes(plots_xlinks_mini$PTB, plots_xlinks$HUR)) +
    geom_point(aes(color=Exon_Originated, alpha=1/100)) +
    theme_bw() +
    scale_y_continuous(limits = c(0 , 1000)) +
    scale_x_continuous(limits = c(0 , 1000)) +
    #ggtitle("Aluexons binding") +
    xlab("PTB") +
    ylab("HUR") +
    scale_colour_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold"))

  U


  ng <- nullGrob()


#tiff("RBP_correlation.tiff" , width = 1000, height = 1000, units = "px")
## export as PDF 40X40

grid.arrange(A, ng, ng, ng, ng, ng,
             B, F, ng, ng, ng, ng,
             C, G, J, ng, ng, ng,
             D, H, K, M, ng, ng,
             E, I, L, N, O, ng,
             P, Q, R, S, T, U,
             ncol=6, nrow =6)


#dev.off()

### Last plots with HUR..hnRNPC

P <- ggplot(plots_xlinks_mini, aes(plots_xlinks_mini$hnRNPC, plots_xlinks$HUR)) +
  geom_point(aes(color=Exon_Originated, size=WU_hg19)) +
  theme_bw() +
  scale_y_continuous(limits = c(0 , 250)) +
  scale_x_continuous(limits = c(0 , 250)) +
  #ggtitle("Aluexons binding") +
  xlab("hnRNPC") +
  ylab("HUR") +
  scale_colour_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
  theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold"))
P

head(plots_xlinks_mini)

P <- ggplot(plots_xlinks_mini, aes(plots_xlinks_mini$hnRNPC, plots_xlinks$HUR, shape=alu_element_hg19_t), alpha=1/1000)+  #, shape=alu_element_hg19_t
  geom_point(aes(color=Exon_Originated, size=WU_hg19)) +
  #scale_shape(factor=(alu_element_hg19_t)) +
  scale_shape_manual(values=1:nlevels(plots_xlinks_mini$alu_element_hg19_t)) +
  theme_bw() +
  scale_y_continuous(limits = c(0 , 250)) +
  scale_x_continuous(limits = c(0 , 250)) +
  #ggtitle("Aluexons binding") +
  xlab("hnRNPC") +
  ylab("HUR") +
  scale_colour_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
  theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold"))
P


### Violin plots


head(plots_xlinks_mini)

plots_xlinks_hnRNPC<- plots_xlinks_mini[,c(5,6)]
plots_xlinks_hnRNPC$protein <- "hnRNPC"
names(plots_xlinks_hnRNPC) <- c( "Exon_Originated", "xlinks", "protein")
plots_xlinks_U2AF65_wt<- plots_xlinks_fraction[,c(1,73,78)]
plots_xlinks_U2AF65_wt$protein <- "U2AF65_wt"
names(plots_xlinks_U2AF65_wt) <- c("aluexon", "Exon_Originated", "xlinks", "protein")
plots_xlinks_U2AF65_hnRNPC_kd<- plots_xlinks_fraction[,c(1,73,79)]
plots_xlinks_U2AF65_hnRNPC_kd$protein <- "U2AF65_hnRNPC_kd"
names(plots_xlinks_U2AF65_hnRNPC_kd) <- c("aluexon", "Exon_Originated", "xlinks", "protein")

plots_xlinks_hnRNPC_frac_U2AF65 <- plots_xlinks_fraction[,c(1,73,80)]
plots_xlinks_hnRNPC_frac_U2AF65$protein <- as.factor("hnRNPC_frac_U2AF65")
names(plots_xlinks_hnRNPC_frac_U2AF65) <- c("aluexon", "Exon_Originated", "xlinks", "protein")



plots_xlinks_long <- rbind(plots_xlinks_hnRNPC, plots_xlinks_U2AF65_wt, plots_xlinks_U2AF65_hnRNPC_kd)
plots_xlinks_long_fract <- rbind(plots_xlinks_hnRNPC_frac_U2AF65, plots_xlinks_U2AF65_frac_U2AF65_hnRNPCKD)




dodge <- position_dodge(width = 0.9)
#Mean
#p_meds <- ddply(final_data_frame, .(region, group), summarise, med = mean(X3SSS, na.rm = TRUE))
#Mode
#p_meds <- ddply(final_data_frame, .(region, group), summarise, med = getmode(X3SSS, na.rm = TRUE))
#Median
p_meds <- ddply(plots_xlinks_long, .(Exon_Originated, protein), summarise, med = median(xlinks, na.rm = TRUE))
p_meds

a<-  ggplot(plots_xlinks_long_fract, aes(factor(protein), xlinks), alpha = 1, colour = "black") +
  geom_boxplot() +
  #scale_y_log10() +
  theme_bw() +
  ggtitle("median quartile 25%, 50%, 75%") +
  xlab("") +
  ylab("length (nt)") +
  theme(text=element_text(size=12),axis.text=element_text(size=12), axis.title=element_text(size=12,face="plain")) +
  scale_x_discrete(limits=c("U2AF65_frac_hnRNPC", "U2AF65_frac_U2AF65_hnRNPCKD"))

a

# count_table <- data.frame(table(final_data_frame$region, final_data_frame$group))
# colnames(count_table) <- c("region", "group", "freq")
# count_table


library(plyr)
### quartile <- median.quartile(clusters$length)     Nejc what its clusters

dodge <- position_dodge(width =1)

pdf_name <- paste0("3SS_" , experiment_name, ".pdf")
#pdf(pdf_name , width=20, height=13)
ggtitle_name <- paste0("3SS ", experiment_name)

plotxlinks<- ggplot(plots_xlinks_long, aes(factor(protein), xlinks, fill = Exon_Originated), alpha = 1, width = 0.5) +
  geom_violin(position = dodge, alpha = 0.6) +
  geom_boxplot(width=0.3, position = dodge, alpha = 1, outlier.shape = NA) +    ### uncoment if you need the boxplot inside
  theme_bw() +
  #theme(panel.background = element_rect(fill='white', colour='black')) +
  scale_y_continuous(limits = c(0, 600)) +
  #scale_y_log10() +
  #ggtitle(ggtitle_name) +
  xlab("") +
  ylab("XLinks Counts") +
  #scale_colour_manual(values=c("#00441b", "#006d2c", "#238b45", "#41ae76", "#66c2a4", "#023858", "#045a8d", "#0570b0", "#3690c0", "#74a9cf", "#67001f", "#980043", "#ce1256", "#e7298a", "#df65b0")) +                                             #"#00441b", "#006d2c", "#0570b0", "#3690c0", "#810f7c", "#88419d", "#8c6bb1"))+  #("#00441b", "#006d2c", "#238b45", "#0570b0", "#3690c0", "#74a9cf", "#810f7c", "#88419d", "#8c6bb1")
  scale_fill_manual(values=c("#2C4259", "#F42837","#36B1BF")) +
  theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold")) +
  #geom_text(data = p_meds, aes(x = Exon_Originated, y = med, label = format((med+1), digits=2)), size = 3, vjust = -1.5, position = dodge, color="white") +
  #geom_text(data = count_table, aes(x = region, y = 15 , label = freq), size = 3, vjust = 0, position = dodge, angle=45) +
  scale_x_discrete(limits=c("hnRNPC", "U2AF65_wt", "U2AF65_hnRNPC_kd"))


plotxlinks

ggsave("CLIPS_aluexons_more30_xlinks.pdf", width=20, height=13)






