---
output: 
  html_document: 
    theme: united
    toc: yes
---

**The START app allows users to visualize RNA-seq data starting with count data.**

- *Explore* the app's features with the example data set pre-loaded by clicking on the tabs above.
- *Upload* your data in the "Input Data" tab.

# Instructions

The app is hosted on the website: https://kcvi.shinyapps.io/START/ 

Code can be found on github: https://github.com/jminnier/STARTapp

To run this app locally on your machine, download R or RStudio and run the following commands once to set up the environment:

```
install.packages(c("reshape2","ggplot2","ggthemes","gplots","ggvis","dplyr","tidyr","DT",
                   "RColorBrewer","pheatmap","shinyBS","plotly","markdown","NMF","scales","heatmaply"))
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite(c("limma","edgeR"))

```
You may now run the shiny app with just one command in R:

```
shiny::runGitHub("STARTapp", "jminnier")
```

<a name="inputdata"></a> 
## Input Data 

You may use this app by

1. Exploring the pre-loaded example data set. This is a pre-loaded mouse RNA-seq example for exploring the app's features.
2. Upload your own data that is either
    i. Count data (or log2-expression data)
    ii. Analyzed data = expression data + p-values and fold changes.
3. Uploading an .RData file containing your data that was previously downloaded from a START app session.

<a name="dataformat"></a> 
### Data Format 

- Must be a .CSV *comma-separated-value* file (you may export from Excel).
- File must have a header row.
- First/Left-hand column(s) must be gene identifiers.
- Format expression column names as `GROUPNAME_REPLICATE#`, e.g. `Treat_1, Treat_2, Treat_3, Control_1, Control_2, High_1, High_2`


#### Count or Expression Data
- Each row denotes a gene, each column denotes a sample.

![](examplecounts.png)

Count data contains read counts for each gene for each sample, along with gene identifiers.

Analysis: When raw counts are uploaded, the data is then analyzed by the app. The app uses the voom method from the ‘limma’ Bioconductor package to transform the raw counts into logged and normalized intensity values. These values are then analyzed via linear regression where gene intensity is regressed on the group factor. P-values from all pairwise regression tests for group effect are computed and Benjamini-Hochberg false discovery rate adjusted p-values are computed for each pairwise comparison. The “log2cpm” values are the log2-counts-per-million values. The “log2cpm_voom” values are the normalized logcpm values from the voom method.  Both methods use an offset of 0.5, which means 0.5 is added to all count values before normalizing (in the case of voom) and log transforming so that 0 counts have non infinite values.

Example file: https://github.com/jminnier/STARTapp/blob/master/data/examplecounts_short.csv

#### Analyzed Data
- Each row denotes a gene, each column denotes a sample.
- Additional columns provide Fold Changes and P-values

![](exampleanalysisdata.png)


Analyzed data must contain some kind of expression measure for each sample (i.e. counts, normalized intensities, CPMs), and a set of p-values with corresponding fold changes for those p-values. For instance, if you have a p-value for the comparison of group1 vs group2, you can upload the observed fold change or log2(fold change) between group1 vs group2. If you have a more complex design and do not have fold changes readily available, you may upload the test statistics or other similar measures of effect size as placeholders. The fold changes are mainly used in the volcano plots. We recommend uploading p-values that are adjusted for multiple comparisons (such as q-values from the qvalue package, or adjusted p-values from p.adjust() function in R).

Example file: <https://github.com/jminnier/STARTapp/blob/master/data/exampleanalysisres_short.csv>


<a name="rdata"></a> 

#### *TIP*: Save Data for Future Upload

After submitting a raw data or analyzed file, you may download the .csv file with the analysis results for your own use (or to upload as an “analyzed data”) or more conveniently click the button “Save Results as RData File for Future Upload” so that you may easily and quickly upload your data to the START app in the future under the “RData from previous START upload” option with one click. 


After uploading your data to START, click red button
![](ex_click_rdata.png)

to download an .RData file to upload your data to START with one click.

Next time use the "Input Data" tab --> "START RData file" option.


<a name="vis"></a> 
# Visualizations

## Group Plots

<a name="pcaplots"></a>
### PCA Plot

This plot uses Principal Component Analysis (PCA) to calculate the principal components of the expression data using data from all genes. Euclidean distances between expression values are used. Samples are projected on the first two principal components (PCs) and the percent variance explained by those PCs are displayed along the x and y axes. Ideally your samples will cluster by group identifier.

 
### Sample Distance Heatmap

This plot displays unsupervised clustering of the Euclidean distances between samples using data from all genes. Again your data should cluster by group.
 

<a name="analysisplots"></a>
## Analysis Plots

These plots use the p-values and fold changes to visualize your data.

<a name="volcano"></a>
### Volcano Plot

This is a scatter plot log fold changes vs –log10(p-values) so that genes with the largest fold changes and smallest p-values are shown on the extreme top left and top right of the plot. Hover over points to see which gene is represented by each point.
 (<https://en.wikipedia.org/wiki/Volcano_plot_(statistics)>)

 

<a name="scatterplots"></a>
### Scatter Plot

This is a scatter plot of average gene expression in one group against another group. This allows the viewer to observe which genes have the largest differences between two groups. The smallest distances will be along the diagonal line, and points far away from the diagonal show the most differences. Hover over points to see which gene is represented by each point.

 
<a name="boxplots"></a>
## Gene Expression Boxplot

Use the search bar to look up genes in your data set. For selected gene(s) the stripchart (dotplot) and boxplots of the expression values are presented for each group. You may plot one or multiple genes along side each other. Hover over points for more information about the data.

<a name="heatmaps"></a>
## Heatmap

A heatmap of expression values are shown, with genes and samples arranged by unsupervised clustering. You may filter on test results as well as P-value cutoffs. By default the top 100 genes (with lowest P-values) are shown.

