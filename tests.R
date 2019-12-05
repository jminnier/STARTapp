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




testdata  <- read_csv("data/testdata_counts_onerep.csv")
tmp3 <- analyze_expression_data(testdata, analysis_method="edgeR", numgeneids = 2)
# need to add option for this somehow

alldata <- read_csv("~/Downloads/joe515.csv")
tmp = extract_count_data(alldata)
tmp3 <- analyze_expression_data(alldata, analysis_method="edgeR")

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
