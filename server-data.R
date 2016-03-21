## ==================================================================================== ##
# START Shiny App for analysis and visualization of transcriptome data.
# Copyright (C) 2014-2016  Jessica Minnier
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
## ================================================================================== ##
## ALL DATA variables
## ================================================================================== ##		

#why is this reactive?
AllRNAdatReactive <- reactive({
  outdat=mousedata[,c("gene.name","gene.id","sample","count","log2cpm","cpm")]
  colnames(outdat) = c("gene.name","gene.id","sample.id","count","log2cpm.edgeR.adjusted","cpm.bowtie.raw")
  return(outdat)		
})

output$downloadData <- downloadHandler(filename = c('all_data.csv'),
                                       content = function(file) {write.csv(AllRNAdatReactive(), file, row.names=FALSE)})



output$outdat <- DT::renderDataTable({
  tmpdat = AllRNAdatReactive()
  DT::datatable(tmpdat)
})
