source("helpers.R")
source("fun-input-analyze-data.R")


x <- load_existing_rdata("data/mousecounts_example.RData")
str(x)

# Normal example data
alldata  <- read_csv("data/mousecounts_example.csv")
tmp = extract_count_data(alldata)
analysis_method="edgeR"
numgeneids = 2
tmp2 <- analyze_expression_data(alldata[1:200,], analysis_method="edgeR", numgeneids = 2)
glimpse(tmp2$data_long)
tmp2 <- analyze_expression_data(alldata[1:200,], analysis_method="voom", numgeneids = 2)
glimpse(tmp2$data_long)
tmp2 <- analyze_expression_data(alldata[1:200,], analysis_method="linear_model", numgeneids = 2)
glimpse(tmp2$data_long)

# One group
testdata  <- read_csv("data/testdata_counts_onegroup.csv")
tmp3 <- analyze_expression_data(testdata, analysis_method="edgeR", numgeneids = 2)
glimpse(tmp3$data_long)
glimpse(tmp3$results)
tmp3 <- analyze_expression_data(testdata, analysis_method="voom", numgeneids = 2)
glimpse(tmp3$data_long)
glimpse(tmp3$results)
tmp3 <- analyze_expression_data(testdata, analysis_method="linear_model", numgeneids = 2)
glimpse(tmp3$data_long)
glimpse(tmp3$results)

# non-counts
testdata  <- read_csv("data/testdata_noncounts.csv")
tmp4 <- analyze_expression_data(testdata, analysis_method="edgeR", numgeneids = 2)
head(tmp4$results)
dotplot_fun(data_long = tmp4$data_long,
            geneids = tmp4$geneids,
            genelabel="unique_id",
            sel_group=tmp4$group_names,
            sel_gene=tmp4$geneids$unique_id[1:2],
            ytype="log2cpm")
gene_pcaplot(data_long= tmp4$data_long,
             valuename= "log2cpm",
             sampleid= tmp4$sampledata$sampleid,
             groupdat= tmp4$sampledata[,"group",drop=FALSE],
             pcnum = 1:2,
             colorfactor="group")


# One replication
testdata  <- read_csv("data/testdata_counts_onerep.csv")
tmp3 <- analyze_expression_data(testdata, analysis_method="edgeR", numgeneids = 2)
# need to add option for this somehow

# alldata <- read_csv("~/Downloads/joe515.csv")
# tmp = extract_count_data(alldata)
# tmp3 <- analyze_expression_data(alldata, analysis_method="edgeR")

# analyzed data not working
# sample heatmaps not working


testdata <- read_csv("data/testdata_analyzed_onecomparison.csv")
data_analyzed = load_analyzed_data(testdata, tmpgenecols = 1:2, tmpexprcols = 3:12,
                         tmpfccols = 13, tmppvalcols = 14, tmpqvalcols = 15, isfclogged = TRUE)
tmpdatlong = data_analyzed$data_long
(tmpynames = tmpdatlong%>%select(-unique_id,-sampleid,-group,-one_of("rep"))%>%colnames())
(tmpgroups = data_analyzed$group_names)
(tmpsamples = as.character(data_analyzed$sampledata$sampleid))
tmpdat = data_analyzed$results
(tmptests = unique(as.character(tmpdat$test)))

# Uploaded data
testdata <- load("data/testdata_counts_prot_uploaded.RData") # start_list
testdata <- load_existing_rdata("data/testdata_counts_prot_uploaded.RData")

rdata_filepath <- "data/mousecounts_example.RData"
testdata <- load_existing_rdata(rdata_filepath)
