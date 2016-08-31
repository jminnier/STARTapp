---
output: html_document
---

This app allows users to visualize RNA-seq data starting with count data.

- Begin exploring the app with the example data set pre-loaded by clicking on the tabs above.
- Upload your data in the "Input Data" tab

#### <a name="features"></a> Features

Explore your data with visualizations:

- clustering (PCA plots, heatmaps)
- group comparisons (scatterplots, volcano plots)
- gene-level boxplots of expression values

<img src="explot_pca.png" alt="PCA Plot" style="height: 200px;"/>
<img src="explot_boxplot.png" alt="Box Plot" style="height: 200px;"/>
<img src="explot_heatmap.png" alt="Heatmap" style="height: 200px;"/>
<img src="explot_volcano.png" alt="Volcano Plot" style="height: 200px;"/>

#### <a name="dataformats"></a> Data Format

- Must be a .CSV *comma-separated-value* file (you may export from Excel).
- File must have a header row.
- First/Left-hand column(s) must be gene identifiers.
- Format expression column names as `GROUPNAME_REPLICATE#`: `Group1_1, Group1_2, Group2_1, Group2_2`


**Count or Expression Data**
- Each row denotes a gene, each column denotes a sample.

![](examplecounts.png)

**Analyzed Data**
- Each row denotes a gene, each column denotes a sample.
- Additional columns provide Fold Changes and P-values.

![](exampleanalysisdata.png)

#### <a name="savedata"></a> TIP: Save Data for Future Upload

After uploading your data to START, click red button "Save Results as START RData File for Future Upload" and use
this .RData file to upload your data to START with one click.


#### <a name="help"></a>Help

Additional help information is provided under the "Help" tab.

#### App Info
This app has been developed by Jessica Minnier, Jiri Sklenar, Anthony Paul Barnes, and Jonathan Nelson
of Oregon Health & Science University, Knight Cardiovascular Institute and School of Public Health.

The source code of START is available on [Github](https://github.com/jminnier/STARTapp).

We would appreciate reports of any issues with the app via the issues option of 
[Github](https://github.com/jminnier/STARTapp) or by emailing minnier-at-ohsu-dot-edu.


