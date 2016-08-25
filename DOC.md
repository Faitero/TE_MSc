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



