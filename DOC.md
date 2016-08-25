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

All the programes here described are piped by `run.sh` script  [run.sh](run.sh) It start with some data reordering, filtration and  3´ss testing. Select only Alus that are on antisense, sort and get uniques. Those unique were intersected with Alu TE from repeat masker. Then we measure the distance from the 3´ss to the Alu TE end.

Same procccess were taken for 

#### Usage:
`bash run.sh`

#### Scripts called
|Link|Description|
|---|-----------|
|[ get_start_position_from_bed_corrected.py ](scr/python/get_start_position_from_bed_corrected.py )| get the start position from a bed file. Used to get the 3'ss of alu exonsCheck the best 3SS score in 2nt upstream and 2nt downstream     |
|[ get_best_3SS.py  ](scr/python/get_best_3SS.py )|   Check the best 3SS score in 2nt upstream and 2nt downstream   |
|[ add_ID_to_bed.py ](scr/python/add_ID_to_bed.py )|   ADD ID to column 4th. Id is formed by the joining of Bed positions   |
|[ 3SS_distance_from_alu.py ](scr/python/3SS_distance_from_alu.py )|  The script will get the distance between 3ss and Aluexon start     |
|[  ](  )|      |




### 5.  iCLIP data

Main script [get_tables_for%20CLIP.sh](scr/bash/get_tables_for%20CLIP.sh) calculates number of xlinks on a bed file.
Script is feed by a URL containing a bedgraph file of RBP xlinks and the overlaping Alu bed file 
Returns bed file with number of xlinks sites that fall on each Alu exon

#### Usage:
`bash xlinks_to_coverage.sh URL_bedgraph.bed  InputFile.bed OutputFile.bed`

#### Scripts called
|Link|Description|
|---|-----------|
|[ xlinks_to_coverage.sh ](scr/bash/xlinks_to_coverage.sh )| Return bed file with number of xlinks sites on thta fall on each Alu exon     |
|[  ](scr/bash/Fig_8.Xlinks.sh )|      |
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



