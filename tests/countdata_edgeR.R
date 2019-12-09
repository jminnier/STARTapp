app <- ShinyDriver$new("../")
app$snapshotInit("countdata_edgeR")

app$setInputs(data_file_type = "upload")
app$uploadFile(datafile = "../data/examplecounts_short.csv") # <-- This should be the path to the file, relative to the app's tests/ directory
app$setInputs(upload_data = "click")
app$setInputs(input_collapse_panel = "analysis_panel")
app$setInputs(sel_gene = "ENSMUSG00000008028_1700008O03Rik")
app$setInputs(sel_gene = c("ENSMUSG00000008028_1700008O03Rik", "ENSMUSG00000036151_Tm6sf2"))
app$setInputs(action_heatmaps = "click")
#app$snapshot()
app$setInputs(datafilter_groups = c("group1", "group2"))
app$setInputs(datafilter_groups = c("group1", "group2", "group3"))
