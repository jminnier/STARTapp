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

seqdata <- read.csv("../old_data/mousecounts_example0.csv",stringsAsFactors = FALSE)
numgeneids <- 2

#scramble real data and replace NAs with random counts
seqdatascr = seqdata
set.seed(100)
tmpind = sample(1:nrow(seqdata))
seqdatascr[,1:numgeneids] = seqdata[tmpind,1:numgeneids]
counts = seqdatascr[,-(1:numgeneids)]
tmpind = which(is.na(counts),arr.ind=T)
if(length(tmpind)>0) {
  set.seed(100)
  tmprepl = round(runif(nrow(tmpind),min=0,max=10))
  counts[tmpind] <- tmprepl
}
seqdatascr = cbind(seqdatascr[,1:numgeneids],counts)
write.csv(seqdatascr,"data/mousecounts_example.csv",row.names = FALSE)

seqdata  <- read.csv("data/mousecounts_example.csv",stringsAsFactors = FALSE)
numgeneids <- 2

sampleid <- colnames(seqdata[,-(1:numgeneids)])
tmpnames <- do.call(rbind,strsplit(sampleid,"_",fixed=TRUE))
group_names <- unique(tmpnames[,1])
group <- tmpnames[,1]
rep_id <- tmpnames[,2]
sampledata = data.frame(sampleid,group,rep_id)

print("analyzing data: input data")

countdata <- seqdata
counts <- countdata[,-(1:numgeneids),drop=FALSE]
geneids <- countdata[,1:numgeneids,drop=FALSE]
tmpkeep = which(apply(is.na(geneids),1,mean)<1) #remove rows with no gene identifiers
print(paste0("Num genes kept: ",length(tmpkeep)," of ", nrow(geneids)))

counts = counts[tmpkeep,,drop=FALSE]
geneids = geneids[tmpkeep,,drop=FALSE]
countdata = countdata[tmpkeep,,drop=FALSE]

geneids = geneids%>%unite_("unique_id",colnames(geneids),remove = FALSE)
#add in catch if geneids not unique?

#handle NAs 
counts[which(is.na(counts),arr.ind=T)] <- 0 #allow choice of this or removal
rownames(counts) = geneids$unique_id

design <- model.matrix(~0+sampledata$group) # allow selection of reference group
colnames(design) = levels(as.factor(sampledata$group))

log2cpm <- cpm(counts, prior.count=0.5, log=TRUE)

#voom+limma
dge <- DGEList(counts=counts) #TMM normalization first
dge <- calcNormFactors(dge)
# v <- voom(dge,design,plot=FALSE)
v <- voom(dge,design,plot=FALSE,normalize.method = "cyclicloess")
# v <- voom(counts,design,plot=TRUE,normalize="quantile") #use this to allow different normalization

#fit <- lmFit(v,design)
#fit <- eBayes(fit)

expr_data = v$E

tmpgroup = sampledata$group
#contrasts(tmpgroup)

lmobj_res = list()
for(ii in 1:length(group_names)) {
  grp = relevel(tmpgroup,ref= group_names[ii] )
  lm.obj = lm(t(expr_data) ~ grp)
  beta.lm = t(lm.obj$coefficients)
  pval.lm = t(lm.pval(lm.obj)$pval)
  pval.adj.lm = apply(pval.lm,2,p.adjust,method="BH")
  
  colnames(beta.lm) = colnames(pval.lm) = colnames(pval.adj.lm) = gsub("grp","",colnames(beta.lm))
  
  tmpout = cbind(melt(beta.lm[,-1,drop=FALSE]),melt(pval.lm[,-1,drop=FALSE])$value,melt(pval.adj.lm[,-1,drop=FALSE])$value)
  colnames(tmpout) = c("unique_id","numer_group","logFC","P.Value","adj.P.Val")
  tmpout$denom_group = group_names[ii]
  tmpout$test = with(tmpout, paste(numer_group,denom_group,sep="/"))
  tmpout = tmpout[,c("unique_id","test","denom_group","numer_group","logFC","P.Value","adj.P.Val")]
  
  
  lmobj_res[[ii]] = tmpout
}
lmobj_res = do.call(rbind,lmobj_res)

pvals = lmobj_res%>%select(unique_id,test,adj.P.Val)%>%spread(test,adj.P.Val)
logfcs = lmobj_res%>%select(unique_id,test,logFC)%>%spread(test,logFC)

colnames(pvals)[-1] = paste0("padj_",colnames(pvals)[-1])
colnames(logfcs)[-1] = paste0("logFC_",colnames(logfcs)[-1])
tmpdat = cbind(geneids,log2cpm)
tmpdat = left_join(tmpdat,logfcs)
tmpdat = left_join(tmpdat,pvals)

data_results_table = tmpdat%>%select(-unique_id) #save this into csv

tmpexprdata = data.frame("unique_id" =geneids$unique_id,expr_data)
tmpcountdata = data.frame("unique_id"=geneids$unique_id,countdata)

tmplog2cpm = data.frame("unique_id"=geneids$unique_id,log2cpm)
log2cpm_long = melt(tmplog2cpm,variable.name = "sampleid",value.name="log2cpm")

countdata_long = melt(tmpcountdata,variable.name = "sampleid",value.name="count")
#countdata_long$log2count = log2(countdata_long$count+.25)

exprdata_long = melt(tmpexprdata,variable.name = "sampleid",value.name="log2cpm_voom")
data_long = left_join(countdata_long,log2cpm_long)
data_long = left_join(data_long,exprdata_long)
data_long$group = do.call(rbind,strsplit(as.character(data_long$sampleid),"_",fixed=TRUE))[,1]

tmpgeneidnames = colnames(geneids%>%select(-unique_id))
data_long = data_long%>%select(-one_of(tmpgeneidnames))


print('analyze data: done')



# tmpcounts = data.frame("unique_id"=geneids$unique_id,counts)
# countdata_long = melt(tmpcounts,variable.name = "sampleid",value.name="count")
# tmplog2cpm = data.frame("unique_id"=geneids$unique_id,log2cpm)
# log2cpm_long = melt(tmplog2cpm,variable.name = "sampleid",value.name="log2cpm")
# data_long = left_join(countdata_long,log2cpm_long)
# 
# data_long$group = do.call(rbind,strsplit(as.character(data_long$sampleid),"_",fixed=TRUE))[,1]
# #expr_data = tmplog2cpm[,-1]
# 
# print('analyze data: done')


data = countdata
results = lmobj_res

#after running analysis pipeline, export this code to another example construction file
save(countdata,group_names,sampledata,results,data_long,geneids,expr_data,
     file="data/mousecounts_example_analysis_results.RData")



# use for example when including p-values and fold changes
load('data/mousecounts_example_analysis_results.RData')

pvals = results%>%select(unique_id,test,adj.P.Val)%>%spread(test,adj.P.Val)
logfcs = results%>%select(unique_id,test,logFC)%>%spread(test,logFC)

colnames(pvals)[-1] = paste0("padj_",colnames(pvals)[-1])
colnames(logfcs)[-1] = paste0("logFC_",colnames(logfcs)[-1])
tmpdat = cbind(geneids,log2cpm)
tmpdat = left_join(tmpdat,logfcs)
tmpdat = left_join(tmpdat,pvals)

# For example uploads into app

example_data_results = tmpdat%>%select(-unique_id)
save(example_data_results,file="data/mousecounts_example_analyzed.RData")
write.csv(example_data_results,file="data/mousecounts_example_analyzed.csv",quote=FALSE,row.names=FALSE)

write.csv(example_data_results[1:100,],file="data/exampleanalysisres_short.csv",quote=FALSE,row.names=FALSE)
write.csv(seqdata[1:100,],"data/examplecounts_short.csv",row.names = FALSE)
