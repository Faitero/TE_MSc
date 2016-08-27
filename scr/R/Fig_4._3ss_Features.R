


##################################
########    Figure 4A       ######
##################################

#### histogram Alu exons evolution
library(ggplot2)
require(ggplot2)



alus <- read.table('/media/igor/DATA/UCL/Evolution_Alus/LiftOver_bedPositions/3SS_Alus/All_Aluexons_3SS_C_distance_to_alu.bed', sep="\t")

distance_end_alu <- alus$V5 < 0


Fig4A<- qplot(alus$V5, geom="histogram",
      binwidth = 1,
      main = "Density of 3SS in alu element",
      xlab = "Alu element position",
      fill=I("grey"),
      col=I("grey")) +
      #alpha=I(.2)) +
      theme_bw() +
      theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold"))
     #,
      #xlim=c(0,320))

Fig4A

## Change theme
theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold")),

## Save filgure
setwd('./Results')
ggsave("3SSS _distance_withind_Alu.pdf", width=20, height=13)




##################################
########    Figure 4B       ######
##################################



## Load rich WIDE table
setwd('./Results')  ## Curro

whole_final_5sp <- read.table("whole_final_5sp.WIDE.RICH.tab.txt", sep="\t", header=TRUE)

queries <- whole_final_5sp
colnames(queries)


## Contingency table by the number of regions
#whole <- as.factor(whole$region)
#class(whole)
count_table <- data.frame(table(whole$region))
colnames(count_table) <- c("region", "freq")
count_table

### line inside box plot could be the mean, mode or median 
#p_meds <- ddply(whole, .(region), summarise, med = mean(X3SSS, na.rm = TRUE))
#Mode
#p_meds <- ddply(whole, .(region), summarise, med = getmode(X3SSS, na.rm = TRUE))
#Median
p_meds <- ddply(whole, .(region), summarise, med = median(X3SSS, na.rm = TRUE))
p_meds


a<-  ggplot(whole, aes(factor(region), X3SSS), alpha = 1, colour = "black") +
  geom_boxplot() + 
  #scale_y_log10() + 
  theme_bw() +
  ggtitle("median quartile 25%, 50%, 75%") + 
  xlab("") + 
  ylab("length (nt)") +
  theme(text=element_text(size=12),axis.text=element_text(size=12), axis.title=element_text(size=12,face="plain")) +
  scale_x_discrete(limits=c("hg19", "panTro4", "nomLeu1", "rheMac3", "calJac3"))                                                # Stablish legend and ortder of genomes

a 

## Overlap several ggplots in one with region as factor
dodge <- position_dodge(width = 0.9)

## Violin plot thta represent the 3´ss score on each specie 
fig4B<- ggplot(whole, aes(factor(region), X3SSS), alpha = 1) +
  geom_violin( position = dodge, alpha = 1, , width = 0.9) + 
  geom_boxplot( position = dodge, alpha = 1, outlier.shape = NA, width = 0.3) + 
  
  theme_bw() +
  scale_y_continuous(limits = c(-40, 20)) +
  #scale_y_log10() +
  ggtitle("3SS score Aluexons") + 
  xlab("") + 
  ylab("3SS score") +
  theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold")) +
  geom_text(data = p_meds, aes(x = region, y = med, label = format(med, digits=2)), size = 4, vjust = -1.5, position = dodge) +
  geom_text(data = count_table, aes(x = region, y = 18 , label = freq), size = 5, vjust = 0, position = dodge, angle=45) +
  scale_x_discrete(limits=c("calJac3", "rheMac3", "nomLeu1","panTro4" ,"hg19"), labels=c("hg19"="Human", "panTro4"="Chimp", "nomLeu1"="Gibbon", "rheMac3"="Rhesus", "calJac3"="Marmoset"))

fig4B

setwd('./Results')
ggsave("3SSS Species.pdf", width=20, height=13)



#################################
#### Figure 4C   ### Heat maps Furthest clustered
#################################

total_furthest <- queries
unique(total_furthest$furthest)
#("hg19"="1", "panTro4"="2", "panPan1"="3", "nomLeu1"="4", "papHam1"="5", "rheMac3"="6", "calJac3"="7")  
col2 <- as.integer(revalue(as.character(total_furthest$furthest), c("hg19"="1", "panTro4"="2", "nomLeu1"="3", "rheMac3"="4", "calJac3"="5")))
furthest_clustered <- cbind(total_furthest, col2)

### Order by furthest
# and WU reversed
furthest_clustered <- furthest_clustered[order(furthest_clustered$col2, -furthest_clustered$WU_hg19),]
# and 3SS reversed hight to less
furthest_clustered <- furthest_clustered[order(furthest_clustered$col2, -furthest_clustered$"X3SSS_hg19"),]

furthest_clustered$col2
furthest_clustered$furthest

furthest_clustered  <- subset(furthest_clustered, furthest_clustered$"X3SSS_hg19">3)     ### Use only those that the 3ss is higher than 3
# 3SS
#furthest_clustered_3 <- data.frame(furthest_clustered$"X3SSS_hg19", furthest_clustered$WU_hg19, furthest_clustered[,seq(20, ncol(furthest_clustered), 11)])
furthest_clustered_3 <- data.frame( furthest_clustered[,seq(20, ncol(furthest_clustered), 11)])
head(furthest_clustered_3)
matdata_furthest_clustered_3 <- data.matrix(furthest_clustered_3)
head(matdata_furthest_clustered_3)
# PolyU
#Divergence_clustered_U <- data.frame(Divergence_clustered$aluexon, Divergence_clustered$alu_substitutions_hg19_t, Divergence_clustered[,seq(21, ncol(total_heat), 11)])
#furthest_clustered_U <- data.frame(furthest_clustered$"X3SSS_hg19", furthest_clustered$WU_hg19, furthest_clustered[,seq(21, ncol(furthest_clustered), 11)])
furthest_clustered_U <- data.frame(furthest_clustered[,seq(21, ncol(furthest_clustered), 11)])
head(furthest_clustered_U)
matdata_furthest_clustered_U <- data.matrix(furthest_clustered_U)
#matdata_furthest_clustered_U <- data.matrix(furthest_clustered_U[3:13])
head(matdata_furthest_clustered_U)



my_palette_3 <- colorRampPalette(c("#FFF8B4", "#FFF88A", "#FFF662",  "#4A75D6", "#4202B2"))(n = 499)

# my_palette_3 <- colorRampPalette(c("#FFF8B4", "#FFF662",  "#4A75D6"))(n = 499)
col_breaks_3 = c(seq(-15, -10,length=100), # for Yeloow
                 seq(-9.90,-2,length=100), # for light YElllow
                 seq(-1.90,3,length=100), # for Grey
                 seq(3.10,10,length=100), # for Bluis
                 seq(10.10,15,length=100)) # for Dark Blue



mat_data_WU <- matdata_furthest_clustered_U
mat_data_3ss <- matdata_furthest_clustered_3 
experiment <- "Clustered_by_Furthest_specie_"



library( pheatmap)
pdf("YellowBlue2_with_PDF2.pdf", width=40, height=26)
tiff("YellowBlue2.tiff" , width = 1000, height = 1000, units = "px")
Fig4C<-  heatmap.2(mat_data_3ss,
                         #cellnote = mat_data,  # same data set for cell labels
                         main = "3SS clustered by furthest and U track lenght", # heat map title
                         #notecol="black",      # change font color of cell labels to black
                         key = TRUE,
                         density.info="density",  #,"density","none"),
                         densadj = 0.5,
                         #density.info="none",  # turns off density plot inside color legend
                         trace="none",         # turns off trace lines inside the heat map
                         margins =c(7, 7),     # widens margins around plot
                         denscol="black",
                         keysize=1.5, 
                         col=my_palette_3,       # use on color palette defined earlier
                         breaks=col_breaks_3,    # enable color transition at specified limits
                         #dendrogram=NULL,     # only draw a row dendrogram
                         RowSideColors=as.character(as.numeric(furthest_clustered$col2)),
                         #colRow= as.character(as.numeric(furthest_clustered$col2)),        ## color categorize rows
                         Rowv=NA,
                         Colv=NA)      # turn off column clustering
dev.off()
ggsave("Fig4C_Heatmap.pdf", width=20, height=13)



