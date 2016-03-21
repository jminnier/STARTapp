
#Julja's code
lm.pval = function(lm.obj){
  # T-test p-value calculation for each parameter's estimated coefficient 
  #  for each response in linear regression output from lm()
  # Code adapted from summary.lm 
  # test: is lm.obj an output from lm()?
  if (is.null(lm.obj$terms) || is.null(lm.obj$qr)) {
    stop("Invalid 'lm' object:  no 'terms' or 'qr' component") }
  if (is.na(lm.obj$df.residual) || 
      (nrow(lm.obj$qr$qr) - lm.obj$rank) != lm.obj$df.residual) {
    warning("Residual degrees of freedom in object suggest this is not an \"lm\" fit") }
  if (lm.obj$rank==0) {
    stop("Regression rank zero: no significance to calculate") }
  # test: one response or many?
  m = !is.null(dim(lm.obj$coefficients)) 
  # extract statistics from lm() output
  p = lm.obj$rank
  Qr = lm.obj$qr
  f = lm.obj$fitted.values
  r = lm.obj$residuals
  rdf = lm.obj$df.residual
  w = lm.obj$weights
  # calculate proportion of variation explained by model
  if(m){ # multiple responses
    if( is.null(w) ){ # non-weighted regression
      mss = if( attr(lm.obj$terms, "intercept")){ 
        colSums( (f - matrix(nrow=nrow(f),ncol=ncol(f),data=colMeans(f),byrow=T))^2 ) 
      } else { colSums( f^2 ) }
      rss =  colSums(r^2) 
    } else { # weighted regression
      mss = if( attr(lm.obj$terms, "intercept")){ 
        colSums(w * (f - matrix(nrow=nrow(f),ncol=ncol(f),data=colSums(w * f/colSums(w)),byrow=T)  )^2 )
      } else {
        colSums(w * f^2)}
      rss = colSums(w * r^2)
      r = sqrt(w) * r
    }
  } else { # single response
    if( is.null(w) ){ # non-weighted regression
      mss = if( attr(lm.obj$terms, "intercept")){ sum( (f - mean(f))^2 ) 
      } else { sum( f^2 ) }
      rss =  sum(r^2) 
    } else { # weighted regression
      mss = if( attr(lm.obj$terms, "intercept")){ 
        sum(w * (f - sum(w * f/sum(w)))^2 )
      } else {
        sum(w * f^2)}
      rss = sum(w * r^2)
      r = sqrt(w) * r
    }
  }
  resvar = rss/rdf
  # calculate squared error per parameter
  p1 = 1:p
  R = chol2inv(Qr$qr[p1, p1, drop = FALSE])
  se = if(m){
    # matrix resvar for proper output dimensions
    # in case of square data, default elementwise vector * matrix is bycol
    #  as required here
    sqrt( diag(R) * matrix(nrow=p,ncol=length(rss),data=resvar,byrow=T) )
    
  }else{ sqrt( diag(R) * resvar ) }
  # calculate T-statistics and p-values from coefficient estimates and SE
  est = if(m){ lm.obj$coefficients[Qr$pivot[p1],] }
  else{ lm.obj$coefficients[Qr$pivot[p1]] }
  tval = est/se
  pval =  2 * pt(abs(tval), rdf, lower.tail = FALSE)
  # calculate overall model performance
  if( p == attr(lm.obj$terms, "intercept") ){ #only intercept parameter
    r.squared = adj.r.squared = 0 
  }else{
    df.int = if( attr(lm.obj$terms,"intercept")){ 1 } else { 0 }
    r.squared = mss / (mss + rss)
    adj.r.squared = 1 - (1-r.squared)*( (nrow(Qr$qr)-df.int)/rdf )
    f.statistic = (mss/(p - df.int))/resvar 
    f.df = p - df.int; f.dendf = rdf 
    p.F = pf(f.statistic, f.df, f.dendf, lower.tail=F)
  }
  if(m){ #removed I() since AsIs complicates things
    ans = list(pval=pval, r.squared=r.squared, adj.r.squared=adj.r.squared, p.F=p.F, f.statistic=f.statistic, f.df=f.df, f.dendf=f.dendf)
  } else {
    ans = list(pval=pval, r.squared=r.squared, adj.r.squared=adj.r.squared, p.F=p.F, f.statistic=f.statistic, f.df=f.df, f.dendf=f.dendf)
  }
  return(ans)
}


