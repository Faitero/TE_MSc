#!/bin/bash -l
cd '/media/igor/DATA/UCL/Evolution_Alus/Join_final_tables' 

####################################################
## get the final table to publication
####################################################

awk '{print $2"\t"$3"\t"$4"\t"$1"\t""\t"$6}' '/media/igor/DATA/UCL/Evolution_Alus/Join_final_tables/whole_final_20160620_5species.csv' > whole_final_20160620_5species_3ss.csv

# Remove first line of heathers
tail -n +2 whole_final_20160620_5species_3ss.csv > whole_final_20160620_5species_3ss_temp.csv && mv whole_final_20160620_5species_3ss_temp.csv whole_final_20160620_5species_3ss.csv
sed 's/\"//g' whole_final_20160620_5species_3ss.csv > whole_final_20160620_5species_3ss_temp.csv && mv whole_final_20160620_5species_3ss_temp.csv whole_final_20160620_5species_3ss.bed

### Check they overlap on antisense -S and get the cordinates from the alu elements
bedtools intersect -S -wo -a whole_final_20160620_5species_3ss.bed -b '/media/igor/DATA/UCL/Evolution_Alus/Raw_Alus/rmsk_hg19_full_family_Alu_elements.bed' > whole_final_20160620_5species_3ss_aluelement.bed

awk '{print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$4"\t"$11}' whole_final_20160620_5species_3ss_aluelement.bed > whole_final_20160620_5species_3ss_aluelement_temp.bed && mv whole_final_20160620_5species_3ss_aluelement_temp.bed whole_final_20160620_5species_3ss_aluelement.bed
## Remove weird \tSINE_Alu
sed 's/\tSINE_Alu//g' whole_final_20160620_5species_3ss_aluelement.bed > whole_final_20160620_5species_3ss_aluelement_temp.bed && mv whole_final_20160620_5species_3ss_aluelement_temp.bed whole_final_20160620_5species_3ss_aluelement.bed

## get the substitutions rate in alu elements
bedtools intersect -s -wo -a whole_final_20160620_5species_3ss.bed -b  '/media/igor/DATA/UCL/Evolution_Alus/Join_final_tables/hg19_allAlus.bed'  > substitutions.bed
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$10"\t"$5}' substitutions.bed > substitutions_temp.bed  && mv substitutions_temp.bed substitutions.bed

## Rest -2nt on positive strand on 3ss and 5ss
python '/home/igor/Dropbox (UCL-MN Team)/AptanaWorkSpaceIgor/2016.02.02_Alus_Jan_Evolution_Hrnpc/get_aluexon_corrected-2nt_positive strand.py' '/media/igor/DATA/UCL/Evolution_Alus/New3SS/3SS_Alus/All_Aluexons.bed' '/media/igor/DATA/UCL/Evolution_Alus/Join_final_tables/All_Aluexons_positive_corrected.bed'


## Merge in R  whole_final_20160620_5species_3ss_aluelement.bed  and All_Aluexons_positive_corrected.bed with '/media/igor/DATA/UCL/Evolution_Alus/Join_final_tables/whole_final_20160620_5species.csv'
## merge_final_tables.R








