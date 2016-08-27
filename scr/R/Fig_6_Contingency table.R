
#################################################
####### histogram of factor frecuencies
#################################################




## Load rich WIDE table
setwd('./Results')  ## Curro

whole_final_5sp <- read.table("whole_final_5sp.WIDE.RICH.tab.txt", sep="\t", header=TRUE)



## Factorize Substitutions
whole_final_5sp$alu_substitutions_hg19_t <- as.factor(whole_final_5sp$alu_substitutions_hg19_t)
whole_final_5sp$alu_substitutions_hg19_t <- factor(whole_final_5sp$alu_substitutions_hg19_t, levels= c("highest", "high", "moderate", "low" ,"lowest", "-1"))
levels(whole_final_5sp$alu_substitutions_hg19_t)

## Aplly the same subset that in the rest of the study
whole_final_5sp <- subset(whole_final_5sp, whole_final_5sp$furthest == "calJac3" & whole_final_5sp$X3SSS_hg19 > 3)
whole_final_5sp <- subset(whole_final_5sp,  whole_final_5sp$X3SSS_hg19 > 3)




## Diferent contingency tables created and studied one by one
par(mfrow=c(3,1))

evolving <- subset(whole_final_5sp, whole_final_5sp$Exon_Originated == "Evolving Exon")
evolving <- evolving$alu_element_hg19_t
evolving <- evolving$dist
evolving <- evolving$UCSCtype
evolving <- evolving$alu_substitutions_hg19_t
evolving <- evolving$chr_hg19_t
evolving <- evolving$WU_hg19
length(evolving)
barplot(table(evolving), main = "Evolving Alu-exons", las=2)

stable <- subset(whole_final_5sp, whole_final_5sp$Exon_Originated == "Constant Exon")
stable <- stable$alu_element_hg19_t
stable <- stable$dist
stable <- stable$UCSCtype
stable <- stable$alu_substitutions_hg19_t
stable <- stable$chr_hg19_t
stable <- stable$WU_hg19
length(stable)
barplot(table(stable), main = "Stable Alu-exons", las=2)

emerging <- subset(whole_final_5sp, whole_final_5sp$Exon_Originated == "Emerging Exon")
emerging <- emerging$alu_element_hg19_t
emerging <- emerging$dist
emerging <- emerging$UCSCtype
emerging <- emerging$alu_substitutions_hg19_t
emerging <- emerging$chr_hg19_t
emerging <- emerging$WU_hg19
length(emerging)
barplot(table(emerging), main = "Emerging Alu-exons", las=2)







#################################################
####### Contingency tables
#################################################

library(ca)
library(vcd)
library(gmodels)


######### Categorical data over U track


## Select only U track longet that 5
whole_final_5sp <- subset(whole_final_5sp, whole_final_5sp$WU_hg19 >= 5)


## Exon originated
#count_table <- data.frame(table(whole_final_5sp$WU_hg19, whole_final_5sp$Exon_Originated))
count_table <- table(whole_final_5sp$WU_hg19, whole_final_5sp$Exon_Originated)
count_table <- structable(whole_final_5sp$Exon_Originated ~ whole_final_5sp$WU_hg19)


## Alu family
count_table <- table(whole_final_5sp$WU_hg19, whole_final_5sp$alu_element_hg19_t)
count_table <- structable(whole_final_5sp$alu_element_hg19_t ~ whole_final_5sp$WU_hg19)

## Substitution
whole_final_5sp <- subset(whole_final_5sp, whole_final_5sp$alu_substitutions_hg19_t !=-1 & whole_final_5sp$WU_hg19 >= 4)
whole_final_5sp$alu_substitutions_hg19_t <- factor(whole_final_5sp$alu_substitutions_hg19_t, levels= c("highest", "high", "moderate", "low" ,"lowest"))
count_table <- structable(whole_final_5sp$alu_substitutions_hg19_t ~ whole_final_5sp$WU_hg19)

## Furthest
count_table <- table(whole_final_5sp$WU_hg19, whole_final_5sp$furthest)


#### Triple contingency table
count_table <- xtabs(~whole_final_5sp$alu_substitutions_hg19_t+whole_final_5sp$WU_hg19+whole_final_5sp$Exon_Originated, data=whole_final_5sp)



## Some summary statistics
#ftable(count_table)
print(count_table)
summary(count_table)

## Factor asociation
assocstats(count_table)

## Table of proportions
prop.table(count_table, 1)

## Marginal Pearson and Chi
margin.table(count_table, 1)

## Plot contingency table
mosaic(count_table, shade = TRUE, legend = TRUE , rot_labels=c(0,0,0,0))   #  main=      las=
assoc(count_table, shade = TRUE, legend=TRUE, rot_labels=c(0,90,0,0))


## CrossTable
CrossTable(whole_final_5sp$WU_hg19, whole_final_5sp$Exon_Originated)

## Chi test
chisq.test(count_table)




######### Categorical data over 3´ss


## Now with the 3'ss in human
#whole_final_5sp$cat_X3SSS_hg19 <- cut(whole_final_5sp$X3SSS_hg19, seq(-40,15,1))
#whole_final_5sp$cat_X3SSS_hg19 <- cut(whole_final_5sp$X3SSS_hg19, c(-40,-30,-20,-10,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))

######### Categorical data over U track
whole_final_5sp$cat_X3SSS_hg19 <- cut(whole_final_5sp$X3SSS_hg19, c(3,4,5,6,7,8,9,10,11,12,13,14,15))


## Exon originated
#count_table <- data.frame(table(whole_final_5sp$WU_hg19, whole_final_5sp$Exon_Originated))
count_table <- table(whole_final_5sp$WU_hg19, whole_final_5sp$Exon_Originated)
count_table <- structable(whole_final_5sp$Exon_Originated ~ whole_final_5sp$cat_X3SSS_hg19)


## Alu family
count_table <- table(whole_final_5sp$WU_hg19, whole_final_5sp$alu_element_hg19_t)
count_table <- structable(whole_final_5sp$alu_element_hg19_t ~ whole_final_5sp$WU_hg19)

## Substitution
whole_final_5sp <- subset(whole_final_5sp, whole_final_5sp$alu_substitutions_hg19_t !=-1 & whole_final_5sp$WU_hg19 >= 4)
whole_final_5sp$alu_substitutions_hg19_t <- factor(whole_final_5sp$alu_substitutions_hg19_t, levels= c("highest", "high", "moderate", "low" ,"lowest"))
count_table <- structable(whole_final_5sp$alu_substitutions_hg19_t ~ whole_final_5sp$cat_X3SSS_hg19)

## Furthest
count_table <- table(whole_final_5sp$WU_hg19, whole_final_5sp$furthest)


#### Triple contingency table
count_table <- xtabs(~whole_final_5sp$alu_substitutions_hg19_t+whole_final_5sp$cat_X3SSS_hg19+whole_final_5sp$Exon_Originated, data=whole_final_5sp)


## Some summary stats
#ftable(count_table)
print(count_table)
summary(count_table)


## Create asociations
assocstats(count_table)

## Proportions table
prop.table(count_table, 1)

## Pearson porb. over residuals
margin.table(count_table, 1)

## Plot Contingeny results
mosaic(count_table, shade = TRUE, legend = TRUE , rot_labels=c(0,0,0,0))   #  main=      las=
assoc(count_table, shade = TRUE, legend=TRUE, rot_labels=c(0,90,0,0))

