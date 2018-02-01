# START App 

START app allows users to upload their gene (or protein, miRNA, etc) expression data and interactively explore and visualize the data and analysis results. It provides simple analysis pipelines as well as the ability to visualize already analyzed data.

The app is hosted on Shinyapps.io here:
<https://kcvi.shinyapps.io/START/>
(please read terms & conditions before uploading sensitive and private data)

Please cite us!

[Nelson, JW, Sklenar J, Barnes AP, Minnier J. (2016) "The START App: A Web-Based RNAseq Analysis and Visualization Resource." Bioinformatics.  doi: 10.1093/bioinformatics/btw624.](http://bioinformatics.oxfordjournals.org/content/early/2016/09/27/bioinformatics.btw624.abstract)

## Run locally

To run this app locally on your machine, download R or RStudio and run the following commands once to set up the environment:
```

install.packages(c("reshape2","ggplot2","ggthemes","gplots","ggvis","dplyr","tidyr","DT",
                   "RColorBrewer","pheatmap","shinyBS","plotly",
                   "markdown","NMF","scales","heatmaply"))
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite(c("limma","edgeR"))

```

You may now run the shiny app with just one command in R:

```
shiny::runGitHub("STARTapp", "jminnier")
```

# Authors

Jonathan Nelson, Jiri Sklenar, Anthony Barnes, Jessica Minnier (code author and maintainer).
*The Knight Cardiovascular Institute and OHSU-PSU School of Public Health, Oregon Health & Science University, Portland, OR 97239-3098, USA.*

We would appreciate reports of any issues with the app via the issues option of 
[Github](https://github.com/jminnier/STARTapp) or by emailing start.app.help-at-gmail.com.

# Instructions

Instructions can be found here: <https://github.com/jminnier/STARTapp/blob/master/instructions/Instructions.md> 

# Licensing

This shiny code is licensed under the GPLv3. Please see the file LICENSE.txt for
information.

    START (Shiny Transcriptome Analysis Resource Tool) App
    Shiny App for analysis and visualization of transcriptome data.
    Copyright (C) 2016  Jessica Minnier

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    You may contact the author of this code, Jessica Minnier, at <minnier@ohsu.edu>
    
# Code adapted for use in app:

- Linear regression p-value extraction code from <https://github.com/ohsu-computational-biology/R-utils>
- `voom` method to allow linear analysis of RNA-seq read counts (Law et al 2014)
- DESeq2 vignette for improving PCA plots (Love et al 2014)


Law, CW, Chen, Y, Shi, W, and Smyth, GK (2014). Voom: precision weights unlock
linear model analysis tools for RNA-seq read counts. Genome Biology 15, R29.

Love MI, Huber W, and Anders S (2014). Moderated
  estimation of fold change and dispersion for RNA-Seq data with DESeq2.
  Genome Biology, 15, 550.

# DOI

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.160626.svg)](https://doi.org/10.5281/zenodo.160626)
