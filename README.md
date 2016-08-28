# TE_MSc

### Example usage

Just run `bash run.sh` to go over all the programs and plot all the figures

## Software version

Pipeline was develop and tested on:

GNU bash, version 4.3.11(1)-release (x86_64-pc-linux-gnu)

Python 2.7.6 [GCC 4.8.2] on linux2

Pycharm IDE PyCharm 2016.1.2
Build #PY-145.844, built on April 8, 2016
JRE: 1.8.0_40-release-b132 x86_64
JVM: OpenJDK 64-Bit Server VM by JetBrains s.r.o

R version 3.3.1 (2016-06-21) -- "Bug in Your Hair"
Platform: x86_64-pc-linux-gnu (64-bit)

Rstudio IDE Version 0.99.489 – © 2009-2015 RStudio, Inc.
Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/538.1 (KHTML, like Gecko) RStudio Safari/538.1 Qt/5.4.0

## Dependencies

The folowing programs must be on the path to be called by the pipeline:

1.Bedtools [   ]( )

2.[ MaxEntScan  ]( https://github.com/Congenica/maxentscan )
NOTE!! splicemodels folder  is needed to run the program. Refer to [ MaxEntScan  ]( https://github.com/Congenica/maxentscan ) for instalation and usage.


3.R Library must be availables for R

    library(vcd)
    library(gmodels)
    library(plyr)
    library(pheatmap)
    library(ggplot2)
    require(gplots)
    require(reshape)
    library(system)
    library(parallel)
    library(GMD)
    library(RColorBrewer)
    library(gridGraphics)
    library(ColorPalette)
    library(scales)
    library(GGally)
    library(sys)
    library(Bio) import SeqIO



4.[ meme_4.10.0 ](http://meme-suite.org/ )   Motif discovery tool HMM based. Check intalation manual.

5.[ Tomtom ](http://meme-suite.org/tools/tomtom )  Search on JASPAR RNA motif database.

6.[ Weblogo3.0  ](http://weblogo.threeplusone.com/manual.html ) sequence alignment tool. Check Manual and Readme on instalation and dependencies needed.

7.Access to iCLIP data at  [ iCOUNT  ]( http://icount.biolab.si/) web server




