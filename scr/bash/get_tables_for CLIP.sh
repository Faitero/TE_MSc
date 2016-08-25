#!/bin/bash -l

## Created on 13 May 2016
## @author: Igor Ruiz de los Mozos

## Script will download iCLIP data from iCount and compute coverage on each Alu evolitionary path
## Return tables of xlinks sites per bed file. Two column per protein
## It will have G files (single genome mapper) and B files (multimapers)

## most of its work is to call xlinks_to_coverage.sh script



#################################
### hnRNPC#######################
#################################

#### Arrange data Evolutionnary path on column 78 ($78)
awk '{print $2"\t"$3"\t"$4"\t"$78"\t"$1"\t"$6"\t"$7}' '/home/igor/Dropbox (UCL-MN Team)/hnRNPC.NMD (1).manuscript/Alus_primates_evolution/Data/whole_final_20160620_5species.csv' > './Evolution_3SS_Aluexons.bed'

# Remove first line of headers and double quotes
tail -n +2 './Evolution_3SS_Aluexons.bed' > './Evolution_Aluexons_temp.bed' && mv './Evolution_Aluexons_temp.bed' './Evolution_3SS_Aluexons.bed'
sed 's/\"//g' './Evolution_3SS_Aluexons.bed' > './Evolution_3SS_Aluexons_temp.bed' && mv './Evolution_3SS_Aluexons_temp.bed' './Evolution_3SS_Aluexons.bed'
## Bs file (multimapers)
bash xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/4964_Ule_hnRNPC_all_bestof/bedGraph_cDNA/group_4964_Ule-hnRNPC-all-bestof_sum_B_hg19--ensembl59_from_1164-1165-1166-1854-1855-...601-602_bedGraph-cDNA.bed.gz Evolution_3SS_Aluexons.bed total_by_Evolution_hnRNPC_B.bed

cat ./total_by_Evolution_hnRNPC_B.bed | grep 'Constant' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat ./total_by_Evolution_hnRNPC_B.bed | grep 'Emerging' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat ./total_by_Evolution_hnRNPC_B.bed | grep 'Evolving' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat ./total_by_Evolution_hnRNPC_B.bed | grep 'NA' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'

## Gs file (single mapper on the genemoe)
bash xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/4964_Ule_hnRNPC_all_bestof/bedGraph_cDNA/group_4964_Ule-hnRNPC-all-bestof_sum_G_hg19--ensembl59_from_1164-1165-1166-1854-1855-...601-602_bedGraph-cDNA-hits-in-genome.bed.gz  Evolution_3SS_Aluexons.bed total_by_Evolution_hnRNPC_G.bed

cat ./total_by_Evolution_hnRNPC_G.bed | grep 'Constant' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat ./total_by_Evolution_hnRNPC_G.bed | grep 'Emerging' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat ./total_by_Evolution_hnRNPC_G.bed | grep 'Evolving' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat ./total_by_Evolution_hnRNPC_G.bed | grep 'NA' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'



###################
## By Substitutions colulmn ($12)

awk '{print $2"\t"$3"\t"$4"\t"$12"\t"$1"\t"$6"\t"$7}' '/home/igor/Dropbox (UCL-MN Team)/hnRNPC.NMD (1).manuscript/Alus_primates_evolution/Data/whole_final_20160620_5species.csv' > './Substitutions_3SS_Aluexons.bed'

tail -n +2 './Substitutions_3SS_Aluexons.bed' > './Substitutions_Aluexons_temp.bed' && mv './Substitutions_Aluexons_temp.bed' './Substitutions_3SS_Aluexons.bed'
sed 's/\"//g' './Substitutions_3SS_Aluexons.bed' > './Substitutions_3SS_Aluexons_temp.bed' && mv './Substitutions_3SS_Aluexons_temp.bed' './Substitutions_3SS_Aluexons.bed'
## Bs
bash xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/4964_Ule_hnRNPC_all_bestof/bedGraph_cDNA/group_4964_Ule-hnRNPC-all-bestof_sum_B_hg19--ensembl59_from_1164-1165-1166-1854-1855-...601-602_bedGraph-cDNA.bed.gz Substitutions_3SS_Aluexons.bed total_by_substitutions_hnRNPC_B.bed

### Calculate average of xlinks on selected exons
cat total_by_substitutions_hnRNPC_B.bed | grep 'highest' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat total_by_substitutions_hnRNPC_B.bed | grep 'high' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat total_by_substitutions_hnRNPC_B.bed | grep 'moderate' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat total_by_substitutions_hnRNPC_B.bed | grep 'low' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat total_by_substitutions_hnRNPC_B.bed | grep 'lowest' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'


## Gs
bash xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/4964_Ule_hnRNPC_all_bestof/bedGraph_cDNA/group_4964_Ule-hnRNPC-all-bestof_sum_G_hg19--ensembl59_from_1164-1165-1166-1854-1855-...601-602_bedGraph-cDNA-hits-in-genome.bed.gz  Substitutions_3SS_Aluexons.bed total_by_substitutions_hnRNPC_G.bed

cat total_by_substitutions_hnRNPC_G.bed | grep 'highest' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat total_by_substitutions_hnRNPC_G.bed | grep 'high' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat total_by_substitutions_hnRNPC_G.bed | grep 'moderate' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat total_by_substitutions_hnRNPC_G.bed | grep 'low' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat total_by_substitutions_hnRNPC_G.bed | grep 'lowest' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'



#################################
### PTB   ####################### 3286
#################################



#### Evolution ($78)
awk '{print $2"\t"$3"\t"$4"\t"$78"\t"$1"\t"$6"\t"$7}' '/home/igor/Dropbox (UCL-MN Team)/hnRNPC.NMD (1).manuscript/Alus_primates_evolution/Data/whole_final_20160620_5species.csv' > './Evolution_3SS_Aluexons.bed'

# Remove first line of headers and double quotes
tail -n +2 './Evolution_3SS_Aluexons.bed' > './Evolution_Aluexons_temp.bed' && mv './Evolution_Aluexons_temp.bed' './Evolution_3SS_Aluexons.bed'
sed 's/\"//g' './Evolution_3SS_Aluexons.bed' > './Evolution_3SS_Aluexons_temp.bed' && mv './Evolution_3SS_Aluexons_temp.bed' './Evolution_3SS_Aluexons.bed'
## Bs
bash xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/3286_PTB-CLIP/bedGraph_cDNA/group_3286_PTB-CLIP_sum_B_hg19--ensembl59_from_2325-2326-2327-2328_bedGraph-cDNA.bed.gz Evolution_3SS_Aluexons.bed total_by_Evolution_PTB_B.bed

cat ./total_by_Evolution_PTB_B.bed | grep 'Constant' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat ./total_by_Evolution_PTB_B.bed | grep 'Emerging' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat ./total_by_Evolution_PTB_B.bed | grep 'Evolving' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat ./total_by_Evolution_PTB_B.bed | grep 'NA' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'

## Gs
bash xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/3286_PTB-CLIP/bedGraph_cDNA/group_3286_PTB-CLIP_sum_G_hg19--ensembl59_from_2325-2326-2327-2328_bedGraph-cDNA-hits-in-genome.bed.gz  Evolution_3SS_Aluexons.bed total_by_Evolution_PTB_G.bed

cat ./total_by_Evolution_PTB_G.bed | grep 'Constant' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat ./total_by_Evolution_PTB_G.bed | grep 'Emerging' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat ./total_by_Evolution_PTB_G.bed | grep 'Evolving' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat ./total_by_Evolution_PTB_G.bed | grep 'NA' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'



###################
## By Substitutions ($12)

awk '{print $2"\t"$3"\t"$4"\t"$12"\t"$1"\t"$6"\t"$7}' '/home/igor/Dropbox (UCL-MN Team)/hnRNPC.NMD (1).manuscript/Alus_primates_evolution/Data/whole_final_20160620_5species.csv' > './Substitutions_3SS_Aluexons.bed'

tail -n +2 './Substitutions_3SS_Aluexons.bed' > './Substitutions_Aluexons_temp.bed' && mv './Substitutions_Aluexons_temp.bed' './Substitutions_3SS_Aluexons.bed'
sed 's/\"//g' './Substitutions_3SS_Aluexons.bed' > './Substitutions_3SS_Aluexons_temp.bed' && mv './Substitutions_3SS_Aluexons_temp.bed' './Substitutions_3SS_Aluexons.bed'
## Bs
bash xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/3286_PTB-CLIP/bedGraph_cDNA/group_3286_PTB-CLIP_sum_B_hg19--ensembl59_from_2325-2326-2327-2328_bedGraph-cDNA.bed.gz Substitutions_3SS_Aluexons.bed total_by_substitutions_PTB_B.bed

### Calculate average of xlinks on selected exons
cat total_by_substitutions_PTB_B.bed | grep 'highest' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat total_by_substitutions_PTB_B.bed | grep 'high' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat total_by_substitutions_PTB_B.bed | grep 'moderate' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat total_by_substitutions_PTB_B.bed | grep 'low' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat total_by_substitutions_PTB_B.bed | grep 'lowest' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'


## Gs
bash xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/3286_PTB-CLIP/bedGraph_cDNA/group_3286_PTB-CLIP_sum_G_hg19--ensembl59_from_2325-2326-2327-2328_bedGraph-cDNA-hits-in-genome.bed.gz  Substitutions_3SS_Aluexons.bed total_by_Evolution_PTB_G.bed

cat total_by_Evolution_PTB_G.bed | grep 'highest' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat total_by_Evolution_PTB_G.bed | grep 'high' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat total_by_Evolution_PTB_G.bed | grep 'moderate' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat total_by_Evolution_PTB_G.bed | grep 'low' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat total_by_Evolution_PTB_G.bed | grep 'lowest' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'



#################################
### PTB   ####################### 4174
#################################



#### Evolution ($78)
awk '{print $2"\t"$3"\t"$4"\t"$78"\t"$1"\t"$6"\t"$7}' '/home/igor/Dropbox (UCL-MN Team)/hnRNPC.NMD (1).manuscript/Alus_primates_evolution/Data/whole_final_20160620_5species.csv' > './Evolution_3SS_Aluexons.bed'

# Remove first line of headers and double quotes
tail -n +2 './Evolution_3SS_Aluexons.bed' > './Evolution_Aluexons_temp.bed' && mv './Evolution_Aluexons_temp.bed' './Evolution_3SS_Aluexons.bed'
sed 's/\"//g' './Evolution_3SS_Aluexons.bed' > './Evolution_3SS_Aluexons_temp.bed' && mv './Evolution_3SS_Aluexons_temp.bed' './Evolution_3SS_Aluexons.bed'
## Bs
bash xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/4174_PTB_LUjh27a_all/bedGraph_cDNA/group_4174_PTB-LUjh27a-all_sum_B_hg19--ensembl59_from_3416-3417-3418-3419_bedGraph-cDNA.bed.gz Evolution_3SS_Aluexons.bed total_by_Evolution_PTB_B.bed

cat ./total_by_Evolution_PTB_B.bed | grep 'Constant' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat ./total_by_Evolution_PTB_B.bed | grep 'Emerging' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat ./total_by_Evolution_PTB_B.bed | grep 'Evolving' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat ./total_by_Evolution_PTB_B.bed | grep 'NA' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'

## Gs
bash xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/4174_PTB_LUjh27a_all/bedGraph_cDNA/group_4174_PTB-LUjh27a-all_sum_G_hg19--ensembl59_from_3416-3417-3418-3419_bedGraph-cDNA-hits-in-genome.bed.gz Evolution_3SS_Aluexons.bed total_by_Evolution_PTB_G.bed

cat ./total_by_Evolution_PTB_G.bed | grep 'Constant' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat ./total_by_Evolution_PTB_G.bed | grep 'Emerging' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat ./total_by_Evolution_PTB_G.bed | grep 'Evolving' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat ./total_by_Evolution_PTB_G.bed | grep 'NA' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'



###################
## By Substitutions ($12)

awk '{print $2"\t"$3"\t"$4"\t"$12"\t"$1"\t"$6"\t"$7}' '/home/igor/Dropbox (UCL-MN Team)/hnRNPC.NMD (1).manuscript/Alus_primates_evolution/Data/whole_final_20160620_5species.csv' > './Substitutions_3SS_Aluexons.bed'

tail -n +2 './Substitutions_3SS_Aluexons.bed' > './Substitutions_Aluexons_temp.bed' && mv './Substitutions_Aluexons_temp.bed' './Substitutions_3SS_Aluexons.bed'
sed 's/\"//g' './Substitutions_3SS_Aluexons.bed' > './Substitutions_3SS_Aluexons_temp.bed' && mv './Substitutions_3SS_Aluexons_temp.bed' './Substitutions_3SS_Aluexons.bed'
## Bs
bash xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/4174_PTB_LUjh27a_all/bedGraph_cDNA/group_4174_PTB-LUjh27a-all_sum_B_hg19--ensembl59_from_3416-3417-3418-3419_bedGraph-cDNA.bed.gz Substitutions_3SS_Aluexons.bed total_by_substitutions_PTB_B.bed

### Calculate average of xlinks on selected exons
cat total_by_substitutions_PTB_B.bed | grep 'highest' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat total_by_substitutions_PTB_B.bed | grep 'high' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat total_by_substitutions_PTB_B.bed | grep 'moderate' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat total_by_substitutions_PTB_B.bed | grep 'low' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat total_by_substitutions_PTB_B.bed | grep 'lowest' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'


## Gs
bash xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/4174_PTB_LUjh27a_all/bedGraph_cDNA/group_4174_PTB-LUjh27a-all_sum_G_hg19--ensembl59_from_3416-3417-3418-3419_bedGraph-cDNA-hits-in-genome.bed.gz Substitutions_3SS_Aluexons.bed total_by_Evolution_PTB_G.bed

cat total_by_Evolution_PTB_G.bed | grep 'highest' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat total_by_Evolution_PTB_G.bed | grep 'high' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat total_by_Evolution_PTB_G.bed | grep 'moderate' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat total_by_Evolution_PTB_G.bed | grep 'low' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat total_by_Evolution_PTB_G.bed | grep 'lowest' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'

##############################
#### Merge tables   ##########
##############################


## Evolution
awk '{print $2"\t"$3"\t"$4"\t"$78"\t"$1"\t"$6"\t"$7}' '/home/igor/Dropbox (UCL-MN Team)/hnRNPC.NMD (1).manuscript/Alus_primates_evolution/Data/whole_final_20160620_5species.csv' > './Evolution_3SS_Aluexons.bed'

# Remove first line of headers and double quotes
tail -n +2 './Evolution_3SS_Aluexons.bed' > './Evolution_Aluexons_temp.bed' && mv './Evolution_Aluexons_temp.bed' './Evolution_3SS_Aluexons.bed'
sed 's/\"//g' './Evolution_3SS_Aluexons.bed' > './Evolution_3SS_Aluexons_temp.bed' && mv './Evolution_3SS_Aluexons_temp.bed' './Evolution_3SS_Aluexons.bed'

## Bs
bash xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/4964_Ule_hnRNPC_all_bestof/bedGraph_cDNA/group_4964_Ule-hnRNPC-all-bestof_sum_B_hg19--ensembl59_from_1164-1165-1166-1854-1855-...601-602_bedGraph-cDNA.bed.gz Evolution_3SS_Aluexons.bed total_by_Evolution_hnRNPC4964_B.bed
## Bs
bash xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/4174_PTB_LUjh27a_all/bedGraph_cDNA/group_4174_PTB-LUjh27a-all_sum_B_hg19--ensembl59_from_3416-3417-3418-3419_bedGraph-cDNA.bed.gz Evolution_3SS_Aluexons.bed total_by_Evolution_PTB4174_B.bed

paste total_by_Evolution_hnRNPC4964_B.bed total_by_Evolution_PTB4174_B.bed | awk '{print $5"\t"$15"\t"$4"\t"$14"\t"$7"\t"$17"\t"($7+1)/($17+1)}' > Evolution_RBP_ratio_B.tab

cat ./Evolution_RBP_ratio_B.tab | grep 'Constant' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat ./Evolution_RBP_ratio_B.tab | grep 'Emerging' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat ./Evolution_RBP_ratio_B.tab | grep 'Evolving' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat ./Evolution_RBP_ratio_B.tab | grep 'NA' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'

## Gs
bash xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/4964_Ule_hnRNPC_all_bestof/bedGraph_cDNA/group_4964_Ule-hnRNPC-all-bestof_sum_G_hg19--ensembl59_from_1164-1165-1166-1854-1855-...601-602_bedGraph-cDNA-hits-in-genome.bed.gz  Evolution_3SS_Aluexons.bed total_by_Evolution_hnRNPC4964_G.bed
## Gs
bash xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/4174_PTB_LUjh27a_all/bedGraph_cDNA/group_4174_PTB-LUjh27a-all_sum_G_hg19--ensembl59_from_3416-3417-3418-3419_bedGraph-cDNA-hits-in-genome.bed.gz Evolution_3SS_Aluexons.bed total_by_Evolution_PTB4174_G.bed


paste total_by_Evolution_hnRNPC4964_G.bed total_by_Evolution_PTB4174_G.bed | awk '{print $5"\t"$15"\t"$4"\t"$14"\t"$7"\t"$17"\t"($7+1)/($17+1)}' > Evolution_RBP_ratio_G.tab

cat ./Evolution_RBP_ratio_G.tab | grep 'Constant' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat ./Evolution_RBP_ratio_G.tab | grep 'Emerging' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat ./Evolution_RBP_ratio_G.tab | grep 'Evolving' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat ./Evolution_RBP_ratio_G.tab | grep 'NA' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'

## Substitutions

awk '{print $2"\t"$3"\t"$4"\t"$12"\t"$1"\t"$6"\t"$7}' '/home/igor/Dropbox (UCL-MN Team)/hnRNPC.NMD (1).manuscript/Alus_primates_evolution/Data/whole_final_20160620_5species.csv' > './Substitutions_3SS_Aluexons.bed'

# Remove first line of headers and double quotes
tail -n +2 './Substitutions_3SS_Aluexons.bed' > './Substitutions_Aluexons_temp.bed' && mv './Substitutions_Aluexons_temp.bed' './Substitutions_3SS_Aluexons.bed'
sed 's/\"//g' './Substitutions_3SS_Aluexons.bed' > './Substitutions_3SS_Aluexons_temp.bed' && mv './Substitutions_3SS_Aluexons_temp.bed' './Substitutions_3SS_Aluexons.bed'

## Bs
bash xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/4964_Ule_hnRNPC_all_bestof/bedGraph_cDNA/group_4964_Ule-hnRNPC-all-bestof_sum_B_hg19--ensembl59_from_1164-1165-1166-1854-1855-...601-602_bedGraph-cDNA.bed.gz Substitutions_3SS_Aluexons.bed total_by_Substitutions_hnRNPC4964_B.bed
## Bs
bash xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/4174_PTB_LUjh27a_all/bedGraph_cDNA/group_4174_PTB-LUjh27a-all_sum_B_hg19--ensembl59_from_3416-3417-3418-3419_bedGraph-cDNA.bed.gz Substitutions_3SS_Aluexons.bed total_by_Substitutions_PTB4174_B.bed

paste total_by_Substitutions_hnRNPC4964_B.bed total_by_Substitutions_PTB4174_B.bed | awk '{print $5"\t"$15"\t"$4"\t"$14"\t"$7"\t"$17"\t"($7+1)/($17+1)}' > Substitutions_RBP_ratio_B.tab

cat Substitutions_RBP_ratio_B.tab | grep 'highest' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat Substitutions_RBP_ratio_B.tab | grep 'high' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat Substitutions_RBP_ratio_B.tab | grep 'moderate' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat Substitutions_RBP_ratio_B.tab | grep 'low' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat Substitutions_RBP_ratio_B.tab | grep 'lowest' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'


## Gs
bash xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/4964_Ule_hnRNPC_all_bestof/bedGraph_cDNA/group_4964_Ule-hnRNPC-all-bestof_sum_G_hg19--ensembl59_from_1164-1165-1166-1854-1855-...601-602_bedGraph-cDNA-hits-in-genome.bed.gz  Substitutions_3SS_Aluexons.bed total_by_Substitutions_hnRNPC4964_G.bed
## Gs
bash xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/4174_PTB_LUjh27a_all/bedGraph_cDNA/group_4174_PTB-LUjh27a-all_sum_G_hg19--ensembl59_from_3416-3417-3418-3419_bedGraph-cDNA-hits-in-genome.bed.gz Substitutions_3SS_Aluexons.bed total_by_Substitutions_PTB4174_G.bed


paste total_by_Substitutions_hnRNPC4964_G.bed total_by_Substitutions_PTB4174_G.bed | awk '{print $5"\t"$15"\t"$4"\t"$14"\t"$7"\t"$17"\t"($7+1)/($17+1)}' > Substitutions_RBP_ratio_G.tab

cat Substitutions_RBP_ratio_G.tab | grep 'highest' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat Substitutions_RBP_ratio_G.tab | grep 'high' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat Substitutions_RBP_ratio_G.tab | grep 'moderate' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat Substitutions_RBP_ratio_G.tab | grep 'low' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'
cat Substitutions_RBP_ratio_G.tab | grep 'lowest' |  awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }'


###############################
####### All CLIPS to table for R   ##### This is the last ONEEE
###############################
awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$1"\t"$6"\t"$7}' '/home/igor/Dropbox (UCL-MN Team)/hnRNPC.NMD (1).manuscript/Alus_primates_evolution/Data/whole_final_20160620_5species.csv' > './Total_evolution_3SS_Aluexons.bed'

# Remove first line of headers and double quotes
tail -n +2 './Total_evolution_3SS_Aluexons.bed' > './Total_evolution_3SS_Aluexons_temp.bed' && mv './Total_evolution_3SS_Aluexons_temp.bed' './Total_evolution_3SS_Aluexons.bed'
sed 's/\"//g' './Total_evolution_3SS_Aluexons.bed' > './Total_evolution_3SS_Aluexons_temp.bed' && mv './Total_evolution_3SS_Aluexons_temp.bed' './Total_evolution_3SS_Aluexons.bed'

## hnRNPC Group 4964 Ule_hnRNPC_all_bestof 
bash ../xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/4964_Ule_hnRNPC_all_bestof/bedGraph_cDNA/group_4964_Ule-hnRNPC-all-bestof_sum_G_hg19--ensembl59_from_1164-1165-1166-1854-1855-...601-602_bedGraph-cDNA-hits-in-genome.bed.gz  Total_evolution_3SS_Aluexons.bed total_hnRNPC4964_G.bed

## PTB Group 4174  PTB_LUjh27a_all
bash ../xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/4174_PTB_LUjh27a_all/bedGraph_cDNA/group_4174_PTB-LUjh27a-all_sum_G_hg19--ensembl59_from_3416-3417-3418-3419_bedGraph-cDNA-hits-in-genome.bed.gz Total_evolution_3SS_Aluexons.bed total_PTB4174_G.bed

## U2AF65_PTB_kd   ## 20130917_LUJH33_1and4_U2AF65_PTB_kd
bash ../xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/4106_20130917_LUJH33_1and4_U2AF65_PTB_kd/bedGraph_cDNA/group_4106_20130917-LUJH33-1and4-U2AF65-PTB-kd_sum_G_hg19--ensembl59_from_3786-3789_bedGraph-cDNA-hits-in-genome.bed.gz Total_evolution_3SS_Aluexons.bed total_U2AF65_PTB_kd_G.bed

## U2AF65_ctrl     ## 20130917_LUJH32_1and4_U2AF65_ctrl
bash ../xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/4105_20130917_LUJH32_1and4_U2AF65_ctrl/bedGraph_cDNA/group_4105_20130917-LUJH32-1and4-U2AF65-ctrl_sum_G_hg19--ensembl59_from_3790-3793_bedGraph-cDNA-hits-in-genome.bed.gz Total_evolution_3SS_Aluexons.bed total_U2AF65_ctrl_G.bed

## U2AF65_wt       ## all_U2AF65_Hela_wt1_hg19 
bash ../xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/4014_all_U2AF65_Hela_wt1_hg19/bedGraph_cDNA/group_4014_all-U2AF65-Hela-wt1-hg19_sum_G_hg19--ensembl59_from_1685-1686-3413-3414_bedGraph-cDNA-hits-in-genome.bed.gz Total_evolution_3SS_Aluexons.bed total_U2AF65_wt_G.bed

## U2AF65_hnRNPC_kd   ## all_U2AF65_Hela_hnRNPCkd_hg19 
bash ../xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/3999_all_U2AF65_Hela_hnRNPCkd_hg19/bedGraph_cDNA/group_3999_all-U2AF65-Hela-hnRNPCkd-hg19_sum_G_hg19--ensembl59_from_1681-1682-1683-3411-3412_bedGraph-cDNA-hits-in-genome.bed.gz Total_evolution_3SS_Aluexons.bed total_U2AF65_hnRNPC_kd_G.bed  

paste total_hnRNPC4964_G.bed total_PTB4174_G.bed total_U2AF65_PTB_kd_G.bed total_U2AF65_ctrl_G.bed total_U2AF65_wt_G.bed total_U2AF65_hnRNPC_kd_G.bed > total_xlinks_aluexons.tab

#### For Hur, TIA and TIAL    cd /media/igor/DATA/UCL/Evolution_Alus/New3SS/CLIPs_on_evolution_groups/New_Proteins

bash ../xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/4361_all_HuR_FlpIn293_hg19/bedGraph_cDNA/group_4361_all-HuR-FlpIn293-hg19_sum_G_hg19--ensembl59_from_4349-4351-4352-4353-4354-...57-4358_bedGraph-cDNA-hits-in-genome.bed.gz Total_evolution_3SS_Aluexons.bed total_HUR_4361_G.bed
bash ../xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/3959_all_TIA1_Hela_hg19/bedGraph_cDNA/group_3959_all-TIA1-Hela-hg19_sum_G_hg19--ensembl59_from_1-21-2293-2294-2295-464_bedGraph-cDNA-hits-in-genome.bed.gz Total_evolution_3SS_Aluexons.bed total_TIA1_3959_G.bed
bash ../xlinks_to_coverage.sh http://icount.fri.uni-lj.si/results/groups/3936_all_TIAL1_Hela_hg19/bedGraph_cDNA/group_3936_all-TIAL1-Hela-hg19_sum_G_hg19--ensembl59_from_2329-2330-2331-248-26-3-4-466-467_bedGraph-cDNA-hits-in-genome.bed.gz Total_evolution_3SS_Aluexons.bed total_TIAL1_3936_G.bed

paste ../last_tables/total_hnRNPC4964_G.bed ../last_tables/total_PTB4174_G.bed ../last_tables/total_U2AF65_PTB_kd_G.bed ../last_tables/total_U2AF65_ctrl_G.bed ../last_tables/total_U2AF65_wt_G.bed ../last_tables/total_U2AF65_hnRNPC_kd_G.bed total_HUR_4361_G.bed total_TIA1_3959_G.bed total_TIAL1_3936_G.bed > total_xlinks_aluexons.tab



