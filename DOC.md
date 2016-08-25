# Alu Evolutionary Analysis

# Plots
## Plots
### Plots
#### Plots
##### Plots
@@ Plots
@ Plots


% porcentaje
$ dolar
& and

### 1.  Main Pipeline

This pipeline will produce all the analysis and the plots from this thesis. Each script contain detailed information on usage and action are comented line by line.

All the programes here described are piped by `run.sh` script  [run.sh](run.sh) It start with some data reordering, filtration and  3´ss testing. Select only Alus that are on antisense, sort and get uniques. Those unique were intersected with Alu TE from repeat masker. Then we measure the distance from the 3´ss to the Alu TE end.

Same procccess were taken for random Alus.

Posteriously call all bellow scripts.

#### Usage:
`bash run.sh`

#### Scripts called
|Link|Description|
|---|-----------|
|[ get_start_position_from_bed_corrected.py ](scr/python/get_start_position_from_bed_corrected.py )| get the start position from a bed file. Used to get the 3'ss of alu exonsCheck the best 3SS score in 2nt upstream and 2nt downstream     |
|[ get_best_3SS.py  ](scr/python/get_best_3SS.py )|   Check the best 3SS score in 2nt upstream and 2nt downstream   |
|[ add_ID_to_bed.py ](scr/python/add_ID_to_bed.py )|   ADD ID to column 4th. Id is formed by the joining of Bed positions   |
|[ 3SS_distance_from_alu.py ](scr/python/3SS_distance_from_alu.py )|  The script will get the distance between 3ss and Aluexon start     |
|[ get3SS_from_random_alu.py ]( scr/python/get3SS_from_random_alu.py )|   get the 3'SS 20 nt inside alu element    |



### 5.  Lift Over and Proccess full table

Main script [ lift_and_procces.sh  ](scr/bash/lift_and_procces.sh )

##  Script is feed by a bed file with Alu exons 3´ss position separated by tab
chr      start   end     Alu_exon_ID    Alu_class    strand  3´ss_distance_to_Alu


1) Lift over the 3´ss to diferent genomes
2) Split the bed file on individual files. Each bed line to a different bed file
3) Get MaxEntSplice site score
4) Get fasta sequence, check that is correct and measure the longest U streech  - Whole alu
5) Get fasta sequence, check that is correct and measure the longest U streech  - right arm
6) Get fasta sequence, check that is correct and measure the longest U streech  - left arm
7) Return a tabular table with all of those results ordered in columns:
        chr start end  aluexon  position  strand  distance_to_alu  X3SSS  LongestUTrack  UTrack_Left UTrack_right


#### Usage:
` bash ./src/bash/lift_and_procces.sh Aluexons_3SS_hg19_Distance.bed OutDIR  `

#### Scripts called
|Link|Description|
|---|-----------|
|[ lift_over_specie.py]( scr/python/lift_over_specie.py)|  lift over bed file to a new bed file specifying the specie conversion. Optional flag to get fasta from those liftovers Seq|
|[ split_bed_record.py ]( scr/python/split_bed_record.py)|      |
|[ flankBEDpositionsStrandSpecific.py ]( scr/python/flankBEDpositionsStrandSpecific.py )|  The script will flank the region in both directions in a new bed file.   |
|[ get_fasta_species.py ]( scr/python/get_fasta_species.py )|  Get fasta sequence from a specified genome  |
|[ check_test_sequence.py ]( scr/python/check_test_sequence.py )|  Check that the fasta sequence is apropiate for downstream analysis    |
|[findLongestStrech.py ]( scr/python/findLongestStrech.py)|  Find the longest strech of a given letter in a string (case insensitive)    |
|[ get_aluexon_from_distance_from_alu2.py ]( scr/python/get_aluexon_from_distance_from_alu2.py)|  The script will get the distance between 3ss and Alu start    |








#################################################
############## Main Script of this research
#################################################
##
##  Script is feed by a bed file with Alu exons 3´ss position. chr      start   end     Alu_exon_ID    Alu_class    strand  3´ss_distance_to_Alu
##
##  It will do:
##                  1) Lift over the 3´ss to diferent genomes
##                  2) Split the bed file on individual files. Each bed line to a different bed file
##                  3) Get MaxEntSplice site score
##                  4) Get fasta sequence, check that is correct and measure the longest U streech  - Whole alu
##                  5) Get fasta sequence, check that is correct and measure the longest U streech  - right arm
##                  6) Get fasta sequence, check that is correct and measure the longest U streech  - left arm
##                  7) Return a tabular table with all of those results
##
##  Output file will be named as the input file but end on .tab
##
##
## Usage:
##          ./lift_and_procces.sh Aluexons_3SS_hg19_Distance.bed OutDIR










### 5.  iCLIP data

Main script [   ](  )
#### Usage:
`   `

#### Scripts called
|Link|Description|
|---|-----------|
|[ ]( )|      |
|[  ]( )|      |
|[  ](  )|      |
|[  ](  )|      |
|[  ](  )|      |



### 5.  iCLIP data

Main script [   ](  )
#### Usage:
`   `

#### Scripts called
|Link|Description|
|---|-----------|
|[ ]( )|      |
|[  ]( )|      |
|[  ](  )|      |
|[  ](  )|      |
|[  ](  )|      |


### 5.  iCLIP data

Main script [   ](  )
#### Usage:
`   `

#### Scripts called
|Link|Description|
|---|-----------|
|[ ]( )|      |
|[  ]( )|      |
|[  ](  )|      |
|[  ](  )|      |
|[  ](  )|      |



### 5.  iCLIP data

Main script [   ](  )
#### Usage:
`   `

#### Scripts called
|Link|Description|
|---|-----------|
|[ ]( )|      |
|[  ]( )|      |
|[  ](  )|      |
|[  ](  )|      |
|[  ](  )|      |



### 5.  iCLIP data

Main script [   ](  )
#### Usage:
`   `

#### Scripts called
|Link|Description|
|---|-----------|
|[ ]( )|      |
|[  ]( )|      |
|[  ](  )|      |
|[  ](  )|      |
|[  ](  )|      |



### 6.  Plots
  Most of the plots were obtained on R using RStudio IDE

|Link|Description|
|---|-----------|
|[Figure 4](scr/R/Fig_4.R)| 3´ss position - 3´ss density - 3´ss heatmap     |
|[Figure 5](scr/R/Fig_5.R)| 3´ss strengthening coupled with U lengthening            |
|[Figure 6](scr/R/Fig_6.R)|      3´ss and U track Contigency table          |
|[Figure 7](scr/bash/Fig_8.Xlinks.sh)|       3´ss Alignments  - Motifs        |
|[Figure 8](scr/bash/Fig_7.motifs.sh)|      RBP Xlinks on Alu evolutionary paths and U track lengths      |


###Source Code Overview
![module diagram](Structure.png "Source Code Overview")


## example usage
Just run `bash run.sh` to go over all the preograms



