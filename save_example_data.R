## ==================================================================================== ##
# START Shiny App for analysis and visualization of transcriptome data.
# Copyright (C) 2016  Jessica Minnier
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
# You may contact the author of this code, Jessica Minnier, at <minnier@ohsu.edu>
## ==================================================================================== ##
## 
source("helpers.R")
source("fun-input-analyze-data.R")

alldata  <- read_csv("data/mousecounts_example.csv")

analyzed_data <- analyze_expression_data(alldata, analysis_method="edgeR", numgeneids = 2)

write.csv(analyzed_data$data_results_table,file="data/mousecounts_example_analyzed.csv",quote=FALSE,row.names=FALSE)
write.csv(analyzed_data$data_results_table[1:100,],file="data/exampleanalysisres_short.csv",quote=FALSE,row.names=FALSE)
write.csv(alldata[1:100,],"data/examplecounts_short.csv",row.names = FALSE)

list2env(analyzed_data, envir=.GlobalEnv)

# LOADED DATA FOR EXAMPLE
save(countdata,group_names,sampledata,results,data_long,geneids,data_results_table,
     file="data/mousecounts_example.RData")



