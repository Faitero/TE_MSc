
## In this scipt I have created two funcions
#  1. WIDE to LONG allow to extract any column from WIDE table, factorise the column/elemnt studied and return a LONG table ready to plot
#  2. Function that is feed with dataframe in LONG format of AluID 3SSS UTRack_whole UtrackLeft UtraclRight
# Call Funcion passing the data_frame to plot
#
# Usage:
#         plot_custom_bar_plots(final_data_frame, experiment_name)


setwd('/Results' )

## Set name of the experiment to label and save the plots
experiment_name <- "Exons Marmoset Originated"

## Call Funcion passing the data_frame to plot
plot_custom_bar_plots(final_data_frame, experiment_name)

#####################################   Functions   #############################################################

############################################
## Function creates a LONG table (need to plot) from the WIDE
############################################

wide_to_long <- function(data_frame, new_column){
  
  # Function transform the wide table of lift ovet to long table in where each row has a factor that represent its procedence and used for later plotting
  #new_column <- "cal_high_3ss"
  #data_frame<-queries
  
  colnames(data_frame) <- NULL
  hg19_q <- data.frame(data_frame[,1], data_frame[,14:24])
  pantro4_q <- data.frame(data_frame[,1], data_frame[,25:35])
  nomLeu1_q<-  data.frame(data_frame[,1], data_frame[,36:46])
  rheMac3_q <- data.frame(data_frame[,1], data_frame[,47:57])
  calJac3_q <- data.frame(data_frame[,1], data_frame[,58:68])
  
  # <- data.frame(data_frame[,1], data_frame[,69:79])
  # <- data.frame(data_frame[,1], data_frame[,80:90])
  #rest <- data.frame(data_frame[,1], data_frame[,91:101])
  
  
  
  out_data_frame <- rbind(hg19_q, pantro4_q, nomLeu1_q, rheMac3_q, calJac3_q)
  
  colnames(out_data_frame) <- c("aluexon", "chr", "start", "end", "position", "strand", "distance_to_alu", "X3SSS", "WU", "U1", "U2", "region")# , "group")  
  
  out_data_frame$group <- as.factor(new_column)
  out_data_frame <- out_data_frame[complete.cases(out_data_frame),]
  
  return(out_data_frame)
  
}

############################################
## Plots for X3SSS WU U1 and U2
############################################

## Function to plot 3SSS WU U1 and U2 at the same time ## Now ussing median  
plot_custom_bar_plots <- function(final_data_frame, experiment_name) {
  
  #### 3SSS## Plots
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  #Mean
  #p_meds <- ddply(final_data_frame, .(region, group), summarise, med = mean(X3SSS, na.rm = TRUE))
  #Mode
  #p_meds <- ddply(final_data_frame, .(region, group), summarise, med = getmode(X3SSS, na.rm = TRUE))
  #Median
  p_meds <- ddply(final_data_frame, .(region, group), summarise, med = median(X3SSS, na.rm = TRUE))
  p_meds
  
  a<-  ggplot(final_data_frame, aes(factor(region), X3SSS), alpha = 1, colour = "black") + 
    geom_boxplot() + 
    #scale_y_log10() + 
    theme_bw() +
    ggtitle("median quartile 25%, 50%, 75%") + 
    xlab("") + 
    ylab("length (nt)") +
    theme(text=element_text(size=12),axis.text=element_text(size=12), axis.title=element_text(size=12,face="plain")) +
    scale_x_discrete(limits=c("hg19", "panTro4", "nomLeu1", "rheMac3", "calJac3"))
  
  a 
  
  count_table <- data.frame(table(final_data_frame$region, final_data_frame$group))
  colnames(count_table) <- c("region", "group", "freq")
  count_table
  
  

  ### quartile <- median.quartile(clusters$length)     Nejc what its clusters
  
  dodge <- position_dodge(width = 0.9)
  
  pdf_name <- paste0("3SS_" , experiment_name, ".pdf")
  #pdf(pdf_name , width=20, height=13)
  ggtitle_name <- paste0("3SS ", experiment_name)
  
  plot3SS<- ggplot(final_data_frame, aes(factor(region), X3SSS, fill = group), alpha = 1, width = 0.5) + 
    geom_violin(position = dodge, alpha = 0.6) + 
    geom_boxplot(width=0.3, position = dodge, alpha = 1, outlier.shape = NA) +    ### uncoment if you need the boxplot inside
    theme_bw() +
    #theme(panel.background = element_rect(fill='white', colour='black')) +
    scale_y_continuous(limits = c(-10 , 16)) +
    #scale_y_log10() +
    ggtitle(ggtitle_name) + 
    xlab("") + 
    ylab("3SS score") +
    #scale_colour_manual(values=c("#00441b", "#006d2c", "#238b45", "#41ae76", "#66c2a4", "#023858", "#045a8d", "#0570b0", "#3690c0", "#74a9cf", "#67001f", "#980043", "#ce1256", "#e7298a", "#df65b0")) +                                             #"#00441b", "#006d2c", "#0570b0", "#3690c0", "#810f7c", "#88419d", "#8c6bb1"))+  #("#00441b", "#006d2c", "#238b45", "#0570b0", "#3690c0", "#74a9cf", "#810f7c", "#88419d", "#8c6bb1")
    scale_fill_manual(values=c("#F42837","#36B1BF","#2C4259")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold")) +
    geom_text(data = p_meds, aes(x = region, y = med, label = format((med+1), digits=2)), size = 3, vjust = -1.5, position = dodge, color="white") +
    geom_text(data = count_table, aes(x = region, y = 15 , label = freq), size = 3, vjust = 0, position = dodge, angle=45) +
    scale_x_discrete(limits=c("hg19", "panTro4", "nomLeu1", "rheMac3", "calJac3"), labels=c("hg19"="Human", "panTro4"="Chimp", "nomLeu1"="Gibbon", "rheMac3"="Rhesus", "calJac3"="Marmoset"))
  
  
  plot3SS
  
  ggsave(pdf_name, width=20, height=13)
  
  #### U track Whole aluexon WU ##
  
  #Mean
  #p_meds <- ddply(final_data_frame, .(region, group), summarise, med = mean(WU, na.rm = TRUE))
  #Mode
  #p_meds <- ddply(final_data_frame, .(region, group), summarise, med = getmode(WU, na.rm = TRUE))
  #Median
  p_meds <- ddply(final_data_frame, .(region, group), summarise, med = median(WU, na.rm = TRUE))
  p_meds
  
  a<-  ggplot(final_data_frame, aes(factor(region), WU), alpha = 1, colour = "black") + 
    geom_boxplot() + 
    #scale_y_log10() + 
    theme_bw() +
    ggtitle("median quartile 25%, 50%, 75%") + 
    xlab("") + 
    ylab("length (nt)") +
    theme(text=element_text(size=12),axis.text=element_text(size=12), axis.title=element_text(size=12,face="plain")) +
    scale_x_discrete(limits=c("hg19", "panTro4", "nomLeu1", "rheMac3", "calJac3"))
  
  a 
  
  count_table <- data.frame(table(final_data_frame$region, final_data_frame$group))
  colnames(count_table) <- c("region", "group", "freq")
  count_table
  
  
  library(plyr)
  ### quartile <- median.quartile(clusters$length)     Nejc what its clusters
  
  dodge <- position_dodge(width = 0.9)
  
  pdf_name <- paste0("WU_" , experiment_name, ".pdf")
  ggtitle_name <- paste0("WU ", experiment_name)
  
  plotWU<- ggplot(final_data_frame, aes(factor(region), WU, fill = group), alpha = 1, width = 0.5) + 
    geom_violin(position = dodge, alpha = 0.6) + 
    geom_boxplot(width=0.3, position = dodge, alpha = 0.6, outlier.shape = NA) +    ### uncoment if you need the boxplot inside
    theme_bw() +
    #theme(panel.background = element_rect(fill='white', colour='black')) +
    scale_y_continuous(limits = c(0 , 30)) +
    #scale_y_log10() +
    ggtitle(ggtitle_name) + 
    xlab("") + 
    ylab("U track lenght") +
    #scale_colour_manual(values=c("#00441b", "#006d2c", "#238b45", "#41ae76", "#66c2a4", "#023858", "#045a8d", "#0570b0", "#3690c0", "#74a9cf", "#67001f", "#980043", "#ce1256", "#e7298a", "#df65b0")) +                                             #"#00441b", "#006d2c", "#0570b0", "#3690c0", "#810f7c", "#88419d", "#8c6bb1"))+  #("#00441b", "#006d2c", "#238b45", "#0570b0", "#3690c0", "#74a9cf", "#810f7c", "#88419d", "#8c6bb1")
    scale_fill_manual(values=c("#F42837","#36B1BF","#2C4259")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold")) +
    geom_text(data = p_meds, aes(x = region, y = med, label = format((med+1), digits=2)), size = 3, vjust = -1.5, position = dodge, color="white") +
    geom_text(data = count_table, aes(x = region, y = 0 , label = freq), size = 3, vjust = 0, position = dodge, angle=45) +
    scale_x_discrete(limits=c("hg19", "panTro4", "nomLeu1", "rheMac3", "calJac3"), labels=c("hg19"="Human", "panTro4"="Chimp", "nomLeu1"="Gibbon", "rheMac3"="Rhesus", "calJac3"="Marmoset"))
  
  
  plotWU
  
  ggsave(pdf_name, width=20, height=13)
  
  #### U1 track aluexon  ##
  
  #Mean
  #p_meds <- ddply(final_data_frame, .(region, group), summarise, med = mean(U1, na.rm = TRUE))
  #Mode
  #p_meds <- ddply(final_data_frame, .(region, group), summarise, med = getmode(U1, na.rm = TRUE))
  #Median
  p_meds <- ddply(final_data_frame, .(region, group), summarise, med = median(U1, na.rm = TRUE))
  p_meds
  
  a<-  ggplot(final_data_frame, aes(factor(region), U1), alpha = 1, colour = "black") + 
    geom_boxplot() + 
    #scale_y_log10() + 
    theme_bw() +
    ggtitle("median quartile 25%, 50%, 75%") + 
    xlab("") + 
    ylab("length (nt)") +
    theme(text=element_text(size=12),axis.text=element_text(size=12), axis.title=element_text(size=12,face="plain")) +
    scale_x_discrete(limits=c("hg19", "panTro4", "nomLeu1", "rheMac3", "calJac3"))
  
  a 
  
  count_table <- data.frame(table(final_data_frame$region, final_data_frame$group))
  colnames(count_table) <- c("region", "group", "freq")
  count_table
  
  
  library(plyr)
  ### quartile <- median.quartile(clusters$length)     Nejc what its clusters
  
  dodge <- position_dodge(width = 0.9)
  
  pdf_name <- paste0("U1_" , experiment_name, ".pdf")
  ggtitle_name <- paste0("U1 ", experiment_name)
  
  plotU1<- ggplot(final_data_frame, aes(factor(region), U1, fill = group), alpha = 1, width = 0.5) + 
    geom_violin(position = dodge, alpha = 0.6) + 
    geom_boxplot(width=0.3, position = dodge, alpha = 1, outlier.shape = NA) +    ### uncoment if you need the boxplot inside
    theme_bw() +
    #theme(panel.background = element_rect(fill='white', colour='black')) +
    scale_y_continuous(limits = c(0 , 30)) +
    #scale_y_log10() +
    ggtitle(ggtitle_name) + 
    xlab("") + 
    ylab("U track lenght") +
    #scale_colour_manual(values=c("#00441b", "#006d2c", "#238b45", "#41ae76", "#66c2a4", "#023858", "#045a8d", "#0570b0", "#3690c0", "#74a9cf", "#67001f", "#980043", "#ce1256", "#e7298a", "#df65b0")) +                                             #"#00441b", "#006d2c", "#0570b0", "#3690c0", "#810f7c", "#88419d", "#8c6bb1"))+  #("#00441b", "#006d2c", "#238b45", "#0570b0", "#3690c0", "#74a9cf", "#810f7c", "#88419d", "#8c6bb1")
    scale_fill_manual(values=c("#F42837","#36B1BF","#2C4259")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold")) +
    geom_text(data = p_meds, aes(x = region, y = med, label = format((med+1), digits=2)), size = 2, vjust = -1.5, position = dodge, color="white") +
    geom_text(data = count_table, aes(x = region, y = 0 , label = freq), size = 3, vjust = 0, position = dodge, angle=45) +
    scale_x_discrete(limits=c("hg19", "panTro4", "nomLeu1", "rheMac3", "calJac3"), labels=c("hg19"="Human", "panTro4"="Chimp", "nomLeu1"="Gibbon", "rheMac3"="Rhesus", "calJac3"="Marmoset"))
  
  
  plotU1
  
  ggsave(pdf_name, width=20, height=13)
  
  
  #### U2 track aluexon  ##
  
  #Mean
  #p_meds <- ddply(final_data_frame, .(region, group), summarise, med = mean(U2, na.rm = TRUE))
  #Mode
  #p_meds <- ddply(final_data_frame, .(region, group), summarise, med = getmode(U2, na.rm = TRUE))
  #Median
  p_meds <- ddply(final_data_frame, .(region, group), summarise, med = median(U2, na.rm = TRUE))
  p_meds
  
  
  a<-  ggplot(final_data_frame, aes(factor(region), U2), alpha = 1, colour = "black") + 
    geom_boxplot() + 
    #scale_y_log10() + 
    theme_bw() +
    ggtitle("median quartile 25%, 50%, 75%") + 
    xlab("") + 
    ylab("length (nt)") +
    theme(text=element_text(size=12),axis.text=element_text(size=12), axis.title=element_text(size=12,face="plain")) +
    scale_x_discrete(limits=c("hg19", "panTro4", "nomLeu1", "rheMac3", "calJac3"))
  
  a 
  
  count_table <- data.frame(table(final_data_frame$region, final_data_frame$group))
  colnames(count_table) <- c("region", "group", "freq")
  count_table
  
  
  library(plyr)
  ### quartile <- median.quartile(clusters$length)     Nejc what its clusters
  
  dodge <- position_dodge(width = 0.9)
  
  pdf_name <- paste0("U2_" , experiment_name, ".pdf")
  ggtitle_name <- paste0("U2 ", experiment_name)
  
  plotU2<- ggplot(final_data_frame, aes(factor(region), U2, fill = group), alpha = 1, width = 0.5) + 
    geom_violin(position = dodge, alpha = 1, alpha = 0.6) + 
    geom_boxplot(width=0.3, position = dodge, alpha = 1, outlier.shape = NA) +    ### uncoment if you need the boxplot inside
    theme_bw() +
    #theme(panel.background = element_rect(fill='white', colour='black')) +
    scale_y_continuous(limits = c(0 , 30)) +
    #scale_y_log10() +
    ggtitle(ggtitle_name) + 
    xlab("") + 
    ylab("U track lenght") +
    #scale_colour_manual(values=c("#00441b", "#006d2c", "#238b45", "#41ae76", "#66c2a4", "#023858", "#045a8d", "#0570b0", "#3690c0", "#74a9cf", "#67001f", "#980043", "#ce1256", "#e7298a", "#df65b0")) +                                             #"#00441b", "#006d2c", "#0570b0", "#3690c0", "#810f7c", "#88419d", "#8c6bb1"))+  #("#00441b", "#006d2c", "#238b45", "#0570b0", "#3690c0", "#74a9cf", "#810f7c", "#88419d", "#8c6bb1")
    scale_fill_manual(values=c("#F42837","#36B1BF","#2C4259")) +
    theme(text=element_text(size=12),axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, vjust=-.5, face="bold"), plot.title=element_text(vjust=1, size=14, face="bold")) +
    geom_text(data = p_meds, aes(x = region, y = med, label = format((med+1), digits=2)), size = 2, vjust = -1.5, position = dodge, color="white") +
    geom_text(data = count_table, aes(x = region, y = 0 , label = freq), size = 3, vjust = 0, position = dodge, angle=45) +
    scale_x_discrete(limits=c("hg19", "panTro4", "nomLeu1", "rheMac3", "calJac3"), labels=c("hg19"="Human", "panTro4"="Chimp", "nomLeu1"="Gibbon", "rheMac3"="Rhesus", "calJac3"="Marmoset"))
  
  
  plotU2
  
  ggsave(pdf_name, width=20, height=13)
  
  
  
  dev.off()
} ## End function

