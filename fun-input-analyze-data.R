
required_data_names <- c("group_names","sampledata","results","data_long","geneids","data_results_table")

load_existing_rdata <- function(rdata_filepath) {
  
  loaded_datanames <- load(rdata_filepath,envir=environment())
  missing_datanames <- setdiff(required_data_names,loaded_datanames)
  validate(
    need(length(missing_datanames)==0, 
         paste("The data file does not contain all the required data objects for this version of the START app.
         Please reload your data using counts/analyzed data and re-save the .RData file.\nData objects missing:", 
               paste0(missing_datanames,collapse=", ")))
  )
  
  
}

# rdata_filepath <- "data/mousecounts_example_analysis_results.RData"
# load_existing_rdata(rdata_filepath)


extract_count_data <- function(alldata, tmpexprcols, tmpgenecols) {
  
  #split expression names into groups
  sampleid <- colnames(alldata[,tmpexprcols])
  tmpnames <- do.call(rbind,strsplit(sampleid,"_",fixed=TRUE))
  group_names <- unique(tmpnames[,1])
  group <- tmpnames[,1]
  rep_id <- tmpnames[,2]
  sampledata = data.frame(sampleid,group,rep_id)
  
  countdata <- alldata[,tmpexprcols,drop=FALSE]
  geneids <- alldata[,tmpgenecols,drop=FALSE]
  
  tmpkeep = which(apply(is.na(geneids),1,mean)<1) #remove rows with no gene identifiers
  print(paste0("Num genes kept after removing empty geneids: ",
               length(tmpkeep)," of ", nrow(geneids)))
  
  validate(need(length(tmpkeep)>0,
                message = "Your data is empty. Please check file format is .csv. 
                                  You may need a non-empty gene identifier column."))
  
  geneids = geneids[tmpkeep,,drop=FALSE]
  countdata = countdata[tmpkeep,,drop=FALSE]
  alldata = countdata[tmpkeep,,drop=FALSE]
  
  # Create unique identifier
  geneids = geneids%>%unite_("unique_id",colnames(geneids),remove = FALSE)
  
  #if geneids not unique
  if(length(unique(geneids$unique_id))<nrow(geneids)) {
    geneids = geneids%>%group_by(unique_id)%>%
      mutate(rn = row_number(unique_id),
             new = ifelse(rn==1,unique_id,paste(unique_id,rn,sep="_")))%>%
      ungroup()%>%mutate(unique_id=new)%>%select(-rn,-new)
  }
  
  rownames(countdata) = geneids$unique_id

  return(list(countdata=countdata, 
              geneids=geneids, 
              group_names=group_names, 
              sampledata=sampledata,
              alldata=alldata))
  
}



analyze_expression_data <- function(alldata, analysis_method, numgeneids = 0) {
  # catch incorrect gene id error, only works if geneids are 1:numbeneids and no other columns are characters
  numgeneids <- max(numgeneids, max(which(sapply(alldata,class)=="character")))
  validate(need(numgeneids>0,
                message = "You have no columns with characters, check that you have a least one column of gene ids as the first column in your file."))
  tmpgenecols = 1:numgeneids
  tmpexprcols = setdiff(1:ncol(alldata),tmpgenecols)
  
  validate(need(length(tmpexprcols)>0,
                message = "Your last column has characters. Check that your count data is numeric and that your gene ids are in the
                first (left) columns only."))
  
  tmpcount <- extract_count_data(alldata, tmpexprcols, tmpgenecols)
  list2env(tmpcount)
  
  #add filter for max # counts
  
  #handle NAs, update this later
  countdata[which(is.na(countdata),arr.ind=T)] <- 0 #allow choice of this or removal
  
  #analyze data
  
  # Check if data appears to be integer counts. If not, skip voom.
  datacounts <- function(input) {
    remainder = sum(apply(input,2,function(k) sum(k%%1,na.rm=T)),na.rm=T)
    if (remainder ==0) {
      TRUE
    } else {
      FALSE
    }
  }
  
  #do not perform voom/edgeR on non-counts and assume log2 uploaded intensities
  dovoom= datacounts(countdata)
  
  # if(not_counts(countdata)){print("Warning: You are uploading data that does not appear to be counts, the analysis pipeline will not be valid!")}
  # validate(
  #   not_counts(countdata)
  # )
  
  print("analyze data: counts")
  
  # Only one group
  if(nlevels(sampledata$group)<2) {
    design <- matrix(1,nrow=nrow(sampledata),ncol=1)
    colnames(design) = "(Intercept)"
  }else{
    design <- model.matrix(~0+sampledata$group) # allow selection of reference group
    colnames(design) = levels(as.factor(sampledata$group))
  }
  
  if(dovoom){
    #voom+limma
    dge <- DGEList(counts=countdata) #TMM normalization first
    dge <- calcNormFactors(dge)
    log2cpm <- cpm(dge, prior.count=0.5, log=TRUE)
    if(max(colSums(design)==1)) {
      # if only one replicate for each group
      v <- voom(dge,normalize.method = "cyclicloess")
    }else{
      v <- voom(dge,design,plot=FALSE,normalize.method = "cyclicloess")
    }
    
    # v <- voom(countdata,design,plot=TRUE,normalize="quantile") #use this to allow different normalization
    #fit <- lmFit(v,design)
    #fit <- eBayes(fit)
    
    expr_data = v$E
  }else{
    print("already normalized")
    countdata2 = countdata
    # crude check for logged data, unlikely to have a logged value >1000
    if(max(countdata)>1000) countdata2 = log2(countdata+0.5)
    log2cpm = countdata2
    expr_data = countdata2
  }
  
  
  
  tmpgroup = sampledata$group
  #contrasts(tmpgroup)
  if(length(group_names)==1) { #If only one group no tests
    lmobj_res = data.frame(matrix(NA,nrow=nrow(expr_data),ncol=6))
    colnames(lmobj_res) = c("test","dneom_group","numer_group","logFC","P.Value","adj.P.Val")
    lmobj_res = cbind("unique_id"=geneids$unique_id,lmobj_res)
    lmobj_res$numer_group = group_names[1]
    lmobj_res$test = "None"
  }else{
    lmobj_res = list()
    for(ii in 1:length(group_names)) {
      grp = relevel(tmpgroup,ref= group_names[ii] )
      design <- model.matrix(~grp)
      
      if(analysis_method=="edgeR") {
        dge <- estimateDisp(dge,design)
        fit <- glmQLFit(dge,design)
        beta <- fit$coefficients[,-1,drop=FALSE]
        pval <- sapply(2:(ncol(design)),
                       function(k) {glmQLFTest(fit,k)$table[,"PValue"]})
      }else{
        lm.obj = lm(t(expr_data) ~ grp)
        beta = t(lm.obj$coefficients)[,-1,drop=FALSE]
        pval = t(lm.pval(lm.obj)$pval)[,-1,drop=FALSE]
      }
      pval.adj = apply(pval,2,p.adjust,method="BH")
      
      colnames(beta) = colnames(pval) = colnames(pval.adj) = 
        gsub("grp","",colnames(beta))
      
      tmpout = cbind(melt(beta),
                     melt(pval)$value,
                     melt(pval.adj)$value)
      colnames(tmpout) = c("unique_id","numer_group","logFC","P.Value","adj.P.Val")
      tmpout$denom_group = group_names[ii]
      tmpout$test = with(tmpout, paste(numer_group,denom_group,sep="/"))
      tmpout = tmpout[,c("unique_id","test","denom_group","numer_group",
                         "logFC","P.Value","adj.P.Val")]
      
      lmobj_res[[ii]] = tmpout
    }
    lmobj_res = do.call(rbind,lmobj_res)
  }
  
  # matrix of pvalues with each column a type of test, same for logfc
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
  if(length(tmpgeneidnames)>0) {
    data_long = data_long%>%select(-one_of(tmpgeneidnames))
  }
  
  print('analyze data: done')
  
  
  return(list('group_names'=group_names,
              'sampledata'=sampledata,
              "results"=lmobj_res,
              "data_long"=data_long, 
              "geneids"=geneids, 
              "data_results_table"=data_results_table))
  
  
  
}


load_analyzed_data <- function(alldata, tmpgenecols, tmpexprcols, tmpfccols, tmppvalcols, tmpqvalcols, isfclogged) {

  tmpcount <- extract_count_data(alldata, tmpexprcols, tmpgenecols)
  list2env(tmpcount)
  
  tmpfc = alldata[,tmpfccols,drop=F]
  if(isfclogged=="No (Log my data please)") {log2(tmpfc)}
  
  fcdata = cbind("unique_id"=geneids$unique_id,tmpfc)
  pvaldata = cbind("unique_id"=geneids$unique_id,alldata[,tmppvalcols,drop=F])
  qvaldata = cbind("unique_id"=geneids$unique_id,alldata[,tmpqvalcols,drop=F])
  
  tmpnames = paste(colnames(fcdata),colnames(qvaldata),sep=":")[-1]
  colnames(fcdata)[-1] = tmpnames
  colnames(pvaldata)[-1] = tmpnames
  colnames(qvaldata)[-1] = tmpnames
  
  fcdatalong = fcdata%>%gather(key = "test",value = "logFC",-1)
  pvaldatalong = pvaldata%>%gather(key = "test",value = "P.Value",-1)
  qvaldatalong = qvaldata%>%gather(key = "test",value = "adj.P.Val",-1)
  tmpres = full_join(fcdatalong,pvaldatalong)
  tmpres = full_join(tmpres,qvaldatalong)
  
  tmpdat = cbind("unique_id"=geneids$unique_id,expr_data)
  tmpdatlong = tmpdat%>%gather(key="sampleid",value="expr",-1)
  data_long = left_join(tmpdatlong,sampledata%>%select(sampleid,group))
  # add summized means by group/unique id for scatterplot
  
  tmpres$test = as.character(tmpres$test)
  
  return(list("group_names"=group_names,
              "sampledata"=sampledata,
              "results"=tmpres,
              "data_long"=data_long, 
              "geneids"=geneids,
              "data_results_table"=alldata))
}