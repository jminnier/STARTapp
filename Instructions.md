# Instructions for the use of the START App


## Input Data


### Example: RNA-Seq gene counts

This is a pre-loaded mouse RNA-seq example for exploring the app's features.

### RData from previous START upload

You may upload an .RData file that you previously downloaded from the START app.

### Upload Data

You may upload your data in two formats, raw counts, or analyzed data.

You must include at least one gene identifier or name.

Your file must have a header row.

The names of your expression/counts must be in the format:
group1_1 where "group1" is the name if your group and "1" is the replicate id number. These must be separated by an underscore. The program determines the names of your groups from these columns.

#### Raw Data: Gene Counts

Raw counts contain read counts for each gene for each sample, along with gene identifiers.

#### Analyzed Data

Analyzedcontains some kind of expression measure for each sample (i.e. counts, normalized intensities, CPMs), and a set of p-values with corresponding fold changes for those p-values. For instance, if you have a p-value for the comparison of group1 vs group2, you can upload the observed fold change or log2(fold change) between group1 vs group2. If you have a more complex design and do not have fold changes readily available, you may upload the test statistics or other similar measures of effect size as placeholders. The fold changes are mainly used in the volcano plots.

