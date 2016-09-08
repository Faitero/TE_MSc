# Alu Evolutionary Analysis


###
### 1.  Main Pipeline

This pipeline will produce all the analysis and the plots from this thesis. Each script contain detailed information on usage and actions are commented line by line.
@@ NOTE -- harsh coding or vital script are noted with **

All the programs here described are piped by `bash run.sh` script  [ ** run.sh](run.sh) It start with some data reordering, filtration and  3´ss testing. Select only Alus that are on antisense, sort and get unique. Those unique were intersected with Alu TE from repeat masker. Then we measure the distance from the 3´ss to the Alu TE end.

Same process were taken for random Alus.

Subsequently it calls all bellow scripts until plots.



#### Usage:
`bash run.sh`

#### Scripts called
|Link|Description|
|---|-----------|
|[ get_start_position_from_bed_corrected.py ](scr/python/get_start_position_from_bed_corrected.py )| get the start position from a bed file. Used to get the 3'ss of Alu exonsCheck the best 3SS score in 2nt upstream and 2nt downstream     |
|[ ** get_best_3SS.py  ](scr/python/get_best_3SS.py )|   Check the best 3SS score in 2nt upstream and 2nt downstream   |
|[ add_ID_to_bed.py ](scr/python/add_ID_to_bed.py )|   ADD ID to column 4th. Id is formed by the joining of Bed positions   |
|[ ** 3SS_distance_from_alu.py ](scr/python/3SS_distance_from_alu.py )|  The script will get the distance between 3ss and Aluexon start     |
|[ get3SS_from_random_alu.py ]( scr/python/get3SS_from_random_alu.py )|   get the 3'SS 20 nt inside alu element    |


###
### 2.  Lift Over and Process full table

Main script [ ** lift_and_procces.sh  ](scr/bash/lift_and_procces.sh )

Script is feed by a bed file with Alu exons 3´ss position separated by tab

    chr      start      end     Alu_exon_ID      Alu_class      strand      3´ss_distance_to_Alu

#
1. Lift over the 3´ss to different genomes
2. Split the bed file on individual files. Each bed line to a different bed file
3. Get MaxEntSplice site score.  Predicted max entropy used
4. Get fasta sequence, check that is correct and measure the longest U streech  - Whole alu
5. Get fasta sequence, check that is correct and measure the longest U streech  - right arm
6. Get fasta sequence, check that is correct and measure the longest U streech  - left arm
7. Return a tabular table with all of those results ordered in columns:
#

        chr     start   end      aluexon    position    strand      distance_to_alu     X3SSS   LongestUTrack   UTrack_Left     UTrack_right

#
#### Usage:
` bash ./src/bash/lift_and_procces.sh Aluexons_3SS_hg19_Distance.bed OutDIR  `

#### Scripts called
|Link|Description|
|---|-----------|
|[ ** lift_over_specie.py]( scr/python/lift_over_specie.py)|  Lift over bed file to a new bed file specifying the specie conversion. Optional flag to get fasta from those liftovers Seq|
|[ split_bed_record.py ]( scr/python/split_bed_record.py)|   Return and split each line in separate file with file name equal to string representing the bed position on the genome   |
|[ flankBEDpositionsStrandSpecific.py ]( scr/python/flankBEDpositionsStrandSpecific.py )|  The script will flank the region in both directions in a new bed file.   |
|[ ** get_fasta_species.py ]( scr/python/get_fasta_species.py )|  Get fasta sequence from a specified genome  |
|[ check_test_sequence.py ]( scr/python/check_test_sequence.py )|  Check that the fasta sequence is appropriate for downstream analysis    |
|[ ** findLongestStrech.py ]( scr/python/findLongestStrech.py)|  Find the longest stretch of a given letter in a string (case insensitive)    |
|[ ** get_aluexon_from_distance_from_alu2.py ]( scr/python/get_aluexon_from_distance_from_alu2.py)|  The script will get the distance between 3ss and Alu start. Then it will create a synthetic Alu element covering all the predicted Alu seq (at least 320 nt).    |


### 3.  3´ss features

Main script [  ** Fig_4._3ss_Features.R  ]( /scr/R/Fig_4._3ss_Features.R )

Exploration of 3´splice site characteristics.

Fig 4A .- Histogram of the distance from the 3´splice site to the start of Alu exon.

Fig 4B .- Violin plots of the density of 3´splice site MaxEntScan score on the furthest Specie that we could find Alu homologous sequences.

Fig 4C .- Heat map of the 3´splice site MaxEntScan score on each specie.

#### Usage:
` Rscript ./scr/R/Fig_4._3ss_Features.R `

### 4.  3´ss coupled with U track lenght

Main script [ ** Fig_5.R ]( scr/R/Fig_5.R)


#### Clasification of Alu exon by the evolutive path.


1. First we get the most distant specie in where we could find homologous Alu sequences inserted (FURTHERST)

2. Then we classified the Alu elements based on the evolutionary dynamics of their 3’ss.

   - ‘Emerging’ Alu exons have a 3’ss with a score less than 3 in the most distant species.
   - ‘Stable’ Alu-exons have a 3’ss higher than 3 in the species most distant to human, and its strength increased towards human by less than 1.
   - ‘Evolving’ Alu-exons have a 3’ss higher than 3 in the species most distant to human, and its strength in human is more than 1+(score in the distant species). For example, if the score in marmoset is 2.5, then the Alu exon is considered as ‘emerging’, if it is 4 in marmoset and 4.5 in human, then it’s considered as ‘stable’, and if it’s 4 in marmoset and 6 in human, then it’s considered as ‘evolving’.

3. Plot 3 splice site strength, Longest U track on Alu, longest U track on left arm, longest U track on right arm



#### Usage:
` Rscipt .src/R/Fig_5.R `


### 5.  Contingency tables

Main script [ Fig_6_Contingency table.R  ]( scr/R/Fig_6_Contingency%20table.R )

To classify Alu elements by the divergence of their sequence in human genome compared to the Alu consensus, I used the nucleotide difference / 1000nt, which is provided for each Alu element in the human genome by the RepeatMasker table ((Smit et al. 2010), hg19, Repeat Library 20090604).). I split them in five equal quartiles ([23-99], [100-121], [122-146], [147-168] and [169-269] substitutions per 1,000 nucleotides).

Then I apply a categorical Pearson test into a two-way contingency table, comparing 3´ splice site strength and U track length against number of nucleotide substitutions.
To transform continuous data to categorical I have just get integer values from the 3´ splice site strength and from the U track length (already discrete data type).

This contingency table were tested against Pearson residual probability using mosaic R package (Pruim, 2016).


#### Usage:
` Rscript ./src/R/Fig_6_Contingency%20table.R  `

### 6.  Alu 3´ss alignment and motif discovery

Main script [ Alu_motives.sh ]( scr/bash/Alu_motives.sh )

Grab needed columns from Wide.tab table

Feed with a Full Alu Evolutionary table in tab format.
WebLogo alignment on Alu evolutionary paths From -60 nt of the 3´s to 15 nt upstream.

Dreme motives on Alu evolutionary paths upstream and downstream the 3´s. From -45 nt of the 3´s to 45 nt upstream
#

Script produces all the plots on Figure 7
#### Usage:
` bash  ./scr/bash/Alu_motives.sh `

#### Scripts called
|Link|Description|
|---|-----------|
|[ get_fasta_species.py ]( scr/python/get_fasta_species.py )|  Get fasta sequence from a specified genome  |
|[ meme athgorithm ](http://meme-suite.org/ )|   Motif discovery tool HMM based   |
|[ Tomtom ](http://meme-suite.org/tools/tomtom )|   Search on JASPAR RNA motif database  |
|[Figure 7](scr/bash/Fig_8.Xlinks.sh)|       3´ss Alignments  - Motifs        |


### 7.  iCLIP RBP binding data

Main script [ get_tables_for_CLIP.sh  ]( scr/bash/get_tables_for%20CLIP.sh )

Script download Xlink sites from iCount iCLIP web server.

Use function [ ** xlinks_to_coverage.sh ]( scr/bash/xlinks_to_coverage.sh ) to assign Xlink count to each Alu exon on a bed file.

Finally R script [ Fig.8:Xlinks.R ]( scr/R/Fig.8:Xlinks.R ) filter the data and plot Violins densities of each class.

#### Usage:
` bash  ./scr/bash/xlinks_to_coverage.sh  `

#### Scripts called
|Link|Description|
|---|-----------|
|[ ** xlinks_to_coverage.sh ]( scr/bash/xlinks_to_coverage.sh )|   Function to assign xlinks to BED file   |
|[ ** get_tables_for_CLIP.sh  ]( scr/bash/get_tables_for%20CLIP.sh )|   Pipeline of iCLIP    |
|[ Fig.8:Xlinks.R ]( scr/R/Fig.8:Xlinks.R )|   Plot iCLIP data   |

### 8.  Plots
  Most of the plots were obtained on R using RStudio IDE

|Link|Description|
|---|-----------|
|[Fig_4._3ss_Features.R](scr/R/Fig_4._3ss_Features.R)| 3´ss position - 3´ss density - 3´ss heatmap     |
|[Figure 5](scr/R/Fig_5.R)| 3´ss strengthening coupled with U lengthening            |
|[Fig_6_Contingency%20table.R](scr/R/Fig_6_Contingency%20table.R)|      3´ss and U track Contingency table          |
|[Fig_7.motifs.sh](scr/bash/Fig_7.motifs.sh)|       3´ss Alignments  - Motifs        |
|[Fig.8:Xlinks.R](scr/bash/Fig.8:Xlinks.R)|      RBP Xlinks on Alu evolutionary paths and U track lengths      |


###Source Code Overview
![module diagram](Data/Structure.png "Source Code Overview")
