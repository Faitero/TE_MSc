#!/bin/bash -l

## Created on 15 May 2016
## @author: Igor Ruiz de los Mozos

## Feed with a Full Alu Evolutionary table in tab format
#  WebLogo alignemt on Alu evolutionary paths From -60 nt of the 3´s to 15 nt upstream
#  Dreme motives on Alu evolutionary paths uspreaam and downstream the 3´s. From -60 nt of the 3´s to 15 nt upstream
#

# Script produces all the plots on Figure 7


## Prepare file

## Arange columns
awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$78"\t"$6"\t"$11"\t"$12"\t"$74"\t"$75"\t"$76}' './Data/whole_final_3sshigh3_marmoset.tab' > './Data/Clasified_3SS_Aluexons.bed'

# Remove first line of heathers
tail -n +2 './Data/Clasified_3SS_Aluexons.bed' > './Clasified_Aluexons_temp.bed' && mv './Data/Clasified_Aluexons_temp.bed' './Data/Clasified_3SS_Aluexons.bed'

# Remove quotes
sed 's/\"//g' './Data/Clasified_3SS_Aluexons.bed' > './Data/Clasified_3SS_Aluexons_temp.bed' && mv './Data/Clasified_3SS_Aluexons_temp.bed' './Data/Clasified_3SS_Aluexons.bed'



#################################################
############## WebLogo Alignment    3ss >3 and Marmoset Originated
############## Now from -60 to +15 web logo
#################################################

## Get -l 60 -r 15 flanking 3ss ## Using bedtools to get all the rest of the columns
bedtools slop -i './Data/Clasified_3SS_Aluexons.bed' -g '/media/igor/DATA/UCL/Genomes/Human/hg19.chrom.sizes.txt' -l 60 -r 15 -s > './Data/Clasified_3SS_Aluexons60_15.bed'
##
python './Data/scr/python/get_fasta_species.py' './Data/Clasified_3SS_Aluexons60_15.bed' "hg19" './Data/Clasified_3SS_Aluexons60_15.fasta'

## Logo Constitutive

cat './Data/Clasified_3SS_Aluexons.bed' | grep "Constant" > './Data/Clasified_3SS_Aluexons_constant.bed'
bedtools slop -i './Data/Clasified_3SS_Aluexons_constant.bed' -g '/media/igor/DATA/UCL/Genomes/Human/hg19.chrom.sizes.txt' -l 60 -r 15 -s > './Data/Clasified_3SS_Aluexons60_15_constant.bed'
python './Data/scr/python/get_fasta_species.py' './Data/Clasified_3SS_Aluexons60_15_constant.bed' "hg19" './Data/Clasified_3SS_Aluexons60_15_constant.fasta'

## Logo Emerging

cat './Data/Clasified_3SS_Aluexons.bed' | grep "Emerging" > './Data/Clasified_3SS_Aluexons_Emerging.bed'
bedtools slop -i './Data/Clasified_3SS_Aluexons_Emerging.bed' -g '/media/igor/DATA/UCL/Genomes/Human/hg19.chrom.sizes.txt' -l 60 -r 15 -s > './Data/Clasified_3SS_Aluexons60_15_Emerging.bed'
python './Data/scr/python/get_fasta_species.py' './Data/Clasified_3SS_Aluexons60_15_Emerging.bed' "hg19" './Data/Clasified_3SS_Aluexons60_15_Emerging.fasta'


## Logo Evolving

cat './Data/Clasified_3SS_Aluexons.bed' | grep "Evolving" > './Data/Clasified_3SS_Aluexons_Evolving.bed'
bedtools slop -i './Data/Clasified_3SS_Aluexons_Evolving.bed' -g '/media/igor/DATA/UCL/Genomes/Human/hg19.chrom.sizes.txt' -l 60 -r 15 -s > './Data/Clasified_3SS_Aluexons60_15_Evolving.bed'
python './Data/scr/python/get_fasta_species.py' './Data/Clasified_3SS_Aluexons60_15_Evolving.bed' "hg19" './Data/Clasified_3SS_Aluexons60_15_Evolving.fasta'




#################################################
############## DREME motif 45 nt Downstream 3ss or 45 nt Upstream 3ss        -norc    strand specific
################################################# cd /media/igor/DATA/UCL/Evolution_Alus/New3SS/CLIPs_on_evolution_groups/Downstream_Upstream



awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$78"\t"$6"\t"$11"\t"$12"\t"$74"\t"$75"\t"$76}' '/home/igor/Dropbox (UCL-MN Team)/hnRNPC.NMD (1).manuscript/Alus_primates_evolution/Data/whole_final_3sshigh3_marmoset.tab' > './Data/Clasified_3SS_Aluexons.bed'

# Remove first line of heathers
tail -n +2 './Data/Clasified_3SS_Aluexons.bed' > './Clasified_Aluexons_temp.bed' && mv './Clasified_Aluexons_temp.bed' './Data/Clasified_3SS_Aluexons.bed'
sed 's/\"//g' './Data/Clasified_3SS_Aluexons.bed' > './Data/Clasified_3SS_Aluexons_temp.bed' && mv './Data/Clasified_3SS_Aluexons_temp.bed' './Data/Clasified_3SS_Aluexons.bed'

#################################################
####### Upstream
#################################################

bedtools slop -i './Data/Clasified_3SS_Aluexons.bed' -g '/media/igor/DATA/UCL/Genomes/Human/hg19.chrom.sizes.txt' -l 45 -r 0 -s > './Data/Clasified_3SS_Aluexons45_0.bed'
python './Data/scr/python/get_fasta_species.py' './Data/Clasified_3SS_Aluexons45_0.bed' "hg19" './Data/Clasified_3SS_Aluexons45_0.fasta'

## Logo Constitutive

cat './Data/Clasified_3SS_Aluexons.bed' | grep "Constant" > './Data/Clasified_3SS_Aluexons_constant.bed'
bedtools slop -i './Data/Clasified_3SS_Aluexons_constant.bed' -g '/media/igor/DATA/UCL/Genomes/Human/hg19.chrom.sizes.txt' -l 45 -r 0 -s > './Data/Clasified_3SS_Aluexons45_0_constant.bed'
python './Data/scr/python/get_fasta_species.py' './Data/Clasified_3SS_Aluexons45_0_constant.bed' "hg19" './Data/Clasified_3SS_Aluexons45_0_constant.fasta'

## Logo Emerging

cat './Data/Clasified_3SS_Aluexons.bed' | grep "Emerging" > './Data/Clasified_3SS_Aluexons_Emerging.bed'
bedtools slop -i './Data/Clasified_3SS_Aluexons_Emerging.bed' -g '/media/igor/DATA/UCL/Genomes/Human/hg19.chrom.sizes.txt' -l 45 -r 0 -s > './Data/Clasified_3SS_Aluexons45_0_Emerging.bed'
python './Data/scr/python/get_fasta_species.py' './Data/Clasified_3SS_Aluexons45_0_Emerging.bed' "hg19" './Data/Clasified_3SS_Aluexons45_0_Emerging.fasta'


## Logo Evolving

cat './Data/Clasified_3SS_Aluexons.bed' | grep "Evolving" > './Data/Clasified_3SS_Aluexons_Evolving.bed'
bedtools slop -i './Data/Clasified_3SS_Aluexons_Evolving.bed' -g '/media/igor/DATA/UCL/Genomes/Human/hg19.chrom.sizes.txt' -l 45 -r 0 -s > './Data/Clasified_3SS_Aluexons45_0_Evolving.bed'
python './Data/scr/python/get_fasta_species.py' './Data/Clasified_3SS_Aluexons45_0_Evolving.bed' "hg19" './Data/Clasified_3SS_Aluexons45_0_Evolving.fasta'


#################################################
### Now with the constitutive as negative control


## Evolving
/home/igor/Programs/meme_4.10.0/bin/dreme -mink 3 -maxk 6 -e 1000 -m 5 -o Evolving_Constant_Negative_Upstream -p './Data/Clasified_3SS_Aluexons45_0_Evolving.fasta' -n './Data/Clasified_3SS_Aluexons45_0_constant.fasta' -norc

## Emerging

/home/igor/Programs/meme_4.10.0/bin/dreme -mink 3 -maxk 6 -e 1000 -m 5 -o Emerging_Constant_Negative_Upstream -p './Data/Clasified_3SS_Aluexons45_0_Emerging.fasta' -n './Data/Clasified_3SS_Aluexons45_0_constant.fasta' -norc


#################################################
### Now with the Evolving as negative control

## Emerging

/home/igor/Programs/meme_4.10.0/bin/dreme -mink 3 -maxk 6 -e 1000 -m 5 -o Emerging_Evolving_Negative_Upstream -p './Data/Clasified_3SS_Aluexons45_0_Emerging.fasta' -n './Data/Clasified_3SS_Aluexons45_0_Evolving.fasta' -norc

## Constant

/home/igor/Programs/meme_4.10.0/bin/dreme -mink 3 -maxk 6 -e 1000 -m 5 -o Contant_Evolving_Negative_Upstream -p './Data/Clasified_3SS_Aluexons45_0_constant.fasta' -n './Data/Clasified_3SS_Aluexons45_0_Evolving.fasta' -norc



#################################################
####### Downstream
#################################################

bedtools slop -i './Data/Clasified_3SS_Aluexons.bed' -g '/media/igor/DATA/UCL/Genomes/Human/hg19.chrom.sizes.txt' -l 0 -r 45 -s > './Data/Clasified_3SS_Aluexons0_45.bed'
python './Data/scr/python/get_fasta_species.py' './Data/Clasified_3SS_Aluexons0_45.bed' "hg19" './Data/Clasified_3SS_Aluexons0_45.fasta'

## Logo Constitutive

cat './Data/Clasified_3SS_Aluexons.bed' | grep "Constant" > './Data/Clasified_3SS_Aluexons_constant.bed'
bedtools slop -i './Data/Clasified_3SS_Aluexons_constant.bed' -g '/media/igor/DATA/UCL/Genomes/Human/hg19.chrom.sizes.txt' -l 0 -r 45 -s > './Data/Clasified_3SS_Aluexons0_45_constant.bed'
python './Data/scr/python/get_fasta_species.py' './Data/Clasified_3SS_Aluexons0_45_constant.bed' "hg19" './Data/Clasified_3SS_Aluexons0_45_constant.fasta'

## Logo Emerging

cat './Data/Clasified_3SS_Aluexons.bed' | grep "Emerging" > './Data/Clasified_3SS_Aluexons_Emerging.bed'
bedtools slop -i './Data/Clasified_3SS_Aluexons_Emerging.bed' -g '/media/igor/DATA/UCL/Genomes/Human/hg19.chrom.sizes.txt' -l 0 -r 45 -s > './Data/Clasified_3SS_Aluexons0_45_Emerging.bed'
python './Data/scr/python/get_fasta_species.py' './Data/Clasified_3SS_Aluexons0_45_Emerging.bed' "hg19" './Data/Clasified_3SS_Aluexons0_45_Emerging.fasta'


## Logo Evolving

cat './Data/Clasified_3SS_Aluexons.bed' | grep "Evolving" > './Data/Clasified_3SS_Aluexons_Evolving.bed'
bedtools slop -i './Data/Clasified_3SS_Aluexons_Evolving.bed' -g '/media/igor/DATA/UCL/Genomes/Human/hg19.chrom.sizes.txt' -l 0 -r 45 -s > './Data/Clasified_3SS_Aluexons0_45_Evolving.bed'
python './Data/scr/python/get_fasta_species.py' './Data/Clasified_3SS_Aluexons0_45_Evolving.bed' "hg19" './Data/Clasified_3SS_Aluexons0_45_Evolving.fasta'


#################################################
### Now with the constitutive as negative control


## Evolving
/home/igor/Programs/meme_4.10.0/bin/dreme -mink 3 -maxk 6 -e 1000 -m 5 -o Evolving_Constant_Negative_Downstream -p './Data/Clasified_3SS_Aluexons0_45_Evolving.fasta' -n './Data/Clasified_3SS_Aluexons0_45_constant.fasta' -norc

## Emerging

/home/igor/Programs/meme_4.10.0/bin/dreme -mink 3 -maxk 6 -e 1000 -m 5 -o Emerging_Constant_Negative_Downstream -p './Data/Clasified_3SS_Aluexons0_45_Emerging.fasta' -n './Data/Clasified_3SS_Aluexons0_45_constant.fasta' -norc


#################################################
### Now with the Evolving as negative control

## Emerging

/home/igor/Programs/meme_4.10.0/bin/dreme -mink 3 -maxk 6 -e 1000 -m 5 -o Emerging_Evolving_Negative_Downstream -p './Data/Clasified_3SS_Aluexons0_45_Emerging.fasta' -n './Data/Clasified_3SS_Aluexons0_45_Evolving.fasta' -norc

## Constant

/home/igor/Programs/meme_4.10.0/bin/dreme -mink 3 -maxk 6 -e 1000 -m 5 -o Contant_Evolving_Negative_Downstream -p './Data/Clasified_3SS_Aluexons0_45_constant.fasta' -n './Data/Clasified_3SS_Aluexons0_45_Evolving.fasta' -norc








