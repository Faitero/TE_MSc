#!/bin/bash -l

## Created on 27 May 2016
## @author: Igor Ruiz de los Mozos

#################################################
############## Main Script of this research
#################################################
##
##  Script is feed by a URL containing a bedgraph file of RBP xlinks per nucleotide on the human genome,
##  and bed file of Alu exon positions
##  Return bed file with number of xlinks sites on thta fall on each Alu exon
##
##  Usage:
##          bash xlinks_to_coverage.sh URL_bedgraph.bed  InputFile.bed OutputFile.bed
##

## Arguments
URLxlink=$1
bed_file=$2
out_bed_file=$3

CLIP_name=$(basename "$URLxlink")
CLIP_name="${CLIP_name%.*}"
echo "URLxlinks" ${URLxlink}
echo "$CLIP_name"

filename=$(basename "$bed_file")
extension="${filename##*.}"
filename="${filename%.*}"
echo "bed file: " ${bed_file}
echo "$filename"
echo "$extension"

### This script gets the distance of 3ss (7th column) from alu and create a 320 imaginary alu element. Info that will pass to output must be in column 4
python '/home/igor/Dropbox (UCL-MN Team)/AptanaWorkSpaceIgor/2016.02.02_Alus_Jan_Evolution_Hrnpc/get_aluexon_from_distance_from_alu2.py' "$bed_file" "$filename"_320.bed

## Get RBP CLIP file from URL. Download a bedgraph file with coount per nucleotide )xlinks)
wget "$URLxlink" -O "$CLIP_name".bed.gz
gunzip "$CLIP_name".bed.gz

## Convert to bed
python '/home/igor/Dropbox (UCL-MN Team)/AptanaWorkSpaceIgor/Converters/BEDgraph2BED.py' "$CLIP_name".bed "$CLIP_name"_BED.bed
rm "$CLIP_name".bed

## Split counts and positive and negative lines
python '/home/igor/Dropbox (UCL-MN Team)/AptanaWorkSpaceIgor/Converters/BED2BED_no_counts_G&B.py' "$CLIP_name"_BED.bed "$CLIP_name"_BED_NC.bed
rm "$CLIP_name"_BED.bed

## Sort optimization for posterior usage
sort -k1,1 -k2,2n "$CLIP_name"_BED_NC.bed > "$CLIP_name"_BED_NC_sort.bed
rm "$CLIP_name"_BED_NC.bed

sort -k1,1 -k2,2n "$filename"_320.bed > "$filename"_sort.bed
rm "$filename"_320.bed

## bedtools to get the coverage of xlilnks per bed file (alu exon)
bedtools coverage -s -sorted -a "$filename"_sort.bed -b "$CLIP_name"_BED_NC_sort.bed > "$out_bed_file"
rm -r "$filename"_sort.bed "$CLIP_name"_BED_NC_sort.bed 

