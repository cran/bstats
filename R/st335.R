## to compute the confidence intervals of the regression parameters

model.check <- function(lmobj){
  MNAME <- lmobj$call
  stopifnot(class(lmobj)=='lm')
  p = length(lmobj$coef)
  n = length(lmobj$resi)
  h <- hat(model.matrix(lmobj)) # (1) hat matrix: leverage values
  SSE = sum((lmobj$resi)^2)
  RMS = SSE/lmobj$df.residual
  Radj = summary(lmobj)$adj.r.squared
  Cp = SSE/(summary(lmobj)$sigma)^2+2*p-n
  AIC = n*log2(SSE/n)+2*p
  BIC = n * log2(SSE/n)+p*log2(n)
  AICc = AIC + 2.*(p+2.)*(p+3.)/(n-p-3.)
  PRESS = sum( (lmobj$resi/(1-h))^2 )
  GOF = c(RMS=RMS,Radj=Radj,Cp=Cp,PRESS=PRESS,AIC=AIC,BIC=BIC, AICc=AICc)
  nam <- names(lmobj$model)
  if(k <- match("(weights)", nam, nomatch = F)){
    warning("Assumption checks are not supported for weighted regressions!")
    RVAL <- list(data.name = MNAME, Criteria = GOF)
  }else{
    ## normality test
    nout = shapiro.test(rstandard(lmobj))
    dname = "W"; pv = nout$p.value; tmethod = nout$method; stat = nout$statistic
    ## test of autocorrelation
    nout = dw.test(lmobj)
    dname = c(dname,"DW"); pv = c(pv,nout$p.value);
    tmethod = c(tmethod,nout$method);
    stat = c(stat,nout$statistic)
    ## bp test for homoscedasticity
    nout = bptest(lmobj)
    dname = c(dname,"BP"); pv = c(pv,nout$p.value);
    tmethod = c(tmethod,nout$method);
    stat = c(stat,nout$statistic)
    ## white test for homoscedasticity
    nout = white.test(lmobj)
    dname = c(dname,"LM"); pv = c(pv,nout$p.value);
    tmethod = c(tmethod,nout$method);
    stat = c(stat,nout$statistic)
    
    out = data.frame(statistic = stat, p.value=pv, Description=tmethod)
    row.names(out) = dname
    RVAL <- list(data.name = MNAME,  Criteria = GOF, HTests = out,VIF=vif(lmobj))
  }
  
  return(RVAL)
}


lm.ci <- function(lmobj,level=0.95)
  {
    if(class(lmobj)!='lm')
      stop("Invalid class type!")
    df0 = lmobj$df
    tmp = summary(lmobj)$coe
    tnames = row.names(tmp)
    out = as.matrix(tmp)
    t0 = abs(qt(0.5*(1.0 + level), df0))
    ll = out[,1] - t0*out[,2]
    ul = out[,1] + t0*out[,2]
    out = data.frame(Estimate=out[,1],lower.limits=ll,upper.limits=ul)
    row.names(out) = tnames
    out
  }

model.test <- function(fmobj,rmobj,alpha=0.05)
  {
    if(class(fmobj)!='lm' || class(rmobj)!='lm')
      stop("Invalid class type(s)!")
    tm1 = names(fmobj$model)[-1]#attr(fmobj$terms,"term.labels")
    tm2 = names(rmobj$model)[-1]#attr(rmobj$terms,"term.labels")
    if(length(tm1)>length(tm2)){
      sele = match(tm2,tm1)
    }else{
      sele = match(tm1,tm2)
      tmp = fmobj; fmobj = rmobj; rmobj = tmp;
    }
    FM = fmobj$call$formula
    RM = rmobj$call$formula
    tm1 = names(fmobj$model)[-1]#attr(fmobj$terms,"term.labels")
    tm2 = names(rmobj$model)[-1]#attr(rmobj$terms,"term.labels")
    if(any(is.na(sele)))
      stop("The two models are not 'Full model' and 'Reduced model'.")
    fmout = anova(fmobj)
    rmout = anova(rmobj)
    fmdf = fmout[nrow(fmout),1]
    rmdf = rmout[nrow(rmout),1]
    fmsse = fmout[nrow(fmout),2]
    rmsse = rmout[nrow(rmout),2]
    F = (rmsse-fmsse)/(rmdf-fmdf)/fmsse*fmdf
    pv = 1-pf(F,rmdf-fmdf,fmdf)
    cv = qf(1-alpha,rmdf-fmdf,fmdf)
    c("Test Stat."=F, "p-value"=pv, "Critical value"=cv)
    out = structure(list(F=F, df1 = rmdf-fmdf, df2 = fmdf,
      pvalue=pv, CriticalValue=cv))
    cat("\n\nNull hypothesis (H0):\t\t the reduced model (RM) is adequate.")
    cat("\nAlternative hypothesis (H1):\t the full model (FM) is adequate.")
    cat("\n\nF-statistic: ", format(F,3),
        " on ", round(rmdf-fmdf,0), " and ",
        round(fmdf,0), "DF, p-value: ", format(pv,10))
    cat("\n\nRM: ", paste(RM[2],RM[1],RM[3]) )
    cat("\nFM: ", paste(FM[2],FM[1],FM[3]),"\n\n")
    invisible(out)
  }


influential.plot <- function(lmobj,type='hadi',ID=FALSE,col=1){
  if(class(lmobj)!='lm')
    stop("Invalid class type!")
  type = match.arg(tolower(type),
    c("hadi","potential-residual","dfits","leverage","cook","hat"))
  stdres = rstandard(lmobj)
  par(mfrow=c(1,1))
  h <- hat(model.matrix(lmobj)) # (1) hat matrix: leverage values
  p = lmobj$rank-1
  n = length(h)
  tmp = anova(lmobj)
  k = nrow(tmp)
  di2 = (lmobj$res)^2/tmp[k,2]
  Hi = h/(1-h) + (p+1)/(1-h)*di2/(1-di2) # (2) Hadi's 
  DFIT = dffits(lmobj) # Welsch and Kuh Measure
  CookD = cooks.distance(lmobj) # Cook's distance
  tmp0 = influence.measures(lmobj)
  infpts = (apply(signif(tmp0$is.inf),1,sum) > 0)
  infpts = c(1:n)[infpts]
  ##  if(is.null(col)) {
    ##    col = rep(1,n); col[infpts] = 2;
  ##  }
  out = cbind(Ri = stdres, Leverage=h,Hadi = Hi, DFIT=DFIT,CookD=CookD)
  if(type=="hadi"){
    Index = 1:length(di2)
    plot(Hi~Index,col=col)
    if(ID)    identify(Hi~Index)
  }else if(type=='potential-residual'){
    di2 = (lmobj$res)^2/tmp[k,2]
    Y = h/(1-h)
    X =  (p+1)/(1-h)*di2/(1-di2)
    plot(Y~X, ylab="Potential",xlab="Residual",col=col)
    if(ID)    identify(Y~X)
  }else if(type == "leverage" || type == "hat"){
    plot(h,type='h',col=col)
    tmp = 2*(p+1)/n
    abline(h=tmp, col='gray')
  }else if(type == "dfits"){
    tmp = 2*sqrt((p+1)/(n-p-1)) 
    plot(DFIT,type='h',col=col)
    abline(h=c(-tmp,tmp), col='gray')
  }else if(type == "cook"){
    plot(CookD,type='h',col=col)
    abline(h=1, col='gray')
  }
  res = structure(list(measures = out, influence.points=infpts))
  invisible(res)
}



residual.plot <- function(lmobj,type='fitted',col=1){
  if(class(lmobj)!='lm')
    stop("Invalid class type!")
  type = match.arg(tolower(type),c("fitted","index","predictor","qqplot"))
  stdres = rstandard(lmobj)
  if(type=="fitted"){
    fittedvalue = lmobj$fitted
    par(mfrow=c(1,1))
    plot(stdres~fittedvalue,ylab="Standardized Residuals",xlab="Fitted values",col=col)
    abline(h=c(-2,2),col="gray")
  }else if(type=="index"){
    Index = 1:nrow(lmobj$mode)
    par(mfrow=c(1,1))
    plot(stdres~Index,ylab="Standardized Residuals",xlab="Index",col=col)
    abline(h=c(-2,2),col="gray")
  }else if(type=='predictor'){
    k = ncol(lmobj$model)-1
    if(!is.null(lmobj$weights)) k = k - 1
    xnames = names(lmobj$model)
    i = ceiling(sqrt(k))
    j = ceiling(k/i)
    par(mfrow=c(j,i))
    for(l in 2:(k+1)){
      x = lmobj$mode[,l]
      plot(stdres~x,ylab="Standardized Residuals",xlab=xnames[l],col=col)
      abline(h=c(-2,2),col="gray")
    }
  }else if(type=='qqplot'){
    par(mfrow=c(1,1))
    qqnorm(stdres, ylab="Standardized Residuals",
           xlab="Normal Scores", main=lmobj$call)
    qqline(stdres)
  }else stop("'type' not supported!")
}

predictor.plot <- function(lmobj,type='av',ID=FALSE,col=1){
  if(class(lmobj)!='lm')
    stop("Invalid class type!")
  type = match.arg(tolower(type),c("av","rc"))
  k = ncol(lmobj$model)-1
  if(!is.null(lmobj$weights)) k = k - 1
  if(k<2) stop("This is a simple linear regression model ('p=1')")
  xnames = names(lmobj$model)
  i = ceiling(sqrt(k))
  j = ceiling(k/i)
  par(mfrow=c(j,i))

  if(type=="av"){
    for(l in 2:(k+1)){
      xname2 = xnames[-l]
      model1 = paste(xname2[1],"~")
      model2 = paste(xnames[l],"~")
      if(k==2){
        model1 = paste(model1,xname2[2])
        model2 = paste(model2,xname2[2])
      }else{
        for(i in 2:(k-1)){
          model1 = paste(model1,xname2[k],"+")
          model2 = paste(model2,xname2[k],"+")
        }
        model1 = paste(model1,xname2[k])
        model2 = paste(model2,xname2[k])
      }
      yres = lm(model1,data=lmobj$model)$res
      xres = lm(model2,data=lmobj$model)$res
      
      plot(yres~xres,,col=col,
           ylab=paste(xnames[1],"-Residuals",sep=""),
           xlab=paste(xnames[l],"-Residuals",sep=""))
      abline(lm(yres~xres)$coef,col='gray')
      if(ID)    identify(yres~xres)
    }
  }else if(type=='rc'){
    betas = lmobj$coef[-1]
    xs = lmobj$model[,-1]
    xnames = names(lmobj$model)[-1]
    for(l in 1:length(betas)){
      xres = betas[l]*xs[,l]
      yres = lmobj$res + xres
        
      plot(yres~xres,,col=col,
           ylab="Residuals+Component",
           xlab=xnames[l])
      abline(lm(yres~xres)$coef,col='gray')
      if(ID)    identify(yres~xres)
    }
  }
}

##  the following code computes sum(e^2)/(n-1)
.myvar <- function(x) sum(x^2)/(length(x)-1)

wls <- function(lmobj,group){
  if(class(lmobj)!='lm')
    stop("Invalid class type!")
  if(is.factor(group)) group = as.numeric(group)
  res = lmobj$res
  sigs = tapply(res,group,.myvar)
  w = mean(res^2)/sigs[group]
  update(lmobj,weights=w)
}

ac <- function(lmobj,type='cochrane', ...)
UseMethod("ac")

ac.default <- function(lmobj,type='cochrane', ...)
stop("No default method for ac. Sorry.")

ac.lm <- function(lmobj,type='cochrane',...){
  if(class(lmobj)!='lm')
    stop("Invalid class type!")
  type = match.arg(tolower(type),c("cochrane","iterative"))
  res = lmobj$res; n=length(res)
  rhohat = sum(res[-1]*res[-n])/sum(res^2)
  mydata = lmobj$model; k =ncol(mydata); l = nrow(mydata)
  if(!is.null(lmobj$weights))
    mydata = mydata[,-k]
  data1 = mydata[-l,]; data2 = mydata[-1,]
  mydata2 = data2 - rhohat*data1
  out = update(lmobj, data=mydata2)
  S0 = sum(out$res^2)
  beta0 = out$coe[1]/(1-rhohat)
  coef = out$coef; coef[1] = beta0
  dwout = dw.test(out)
  if(type=="iterative"){
    a = max(0.0,rhohat-.2); b=min(1.0, rhohat+.2)
    rhos = seq(a,b,length=100)
    for(rhat in rhos){
      mydata2 = data2 - rhat*data1
      out2 = update(lmobj, data=mydata2)
      S1 = sum(out2$res^2)
      if(S1 < S0){
        beta0 = out2$coe[1]/(1-rhat)
        coef = out2$coef; coef[1] = beta0
        dwout = dw.test(out2)
        rhohat = rhat
        out = out2
        S0 = S1
      }
    }

  }
  out = list(coefficients = coef, rhohat=rhohat, dwtest=dwout,lmobj=out)
  ##  invisible(out)
  out
}

vif <- function(object, ...)
UseMethod("vif")

vif.default <- function(object, ...)
stop("No default method for vif. Sorry.")

vif.lm <- function(object, ...) {
  nam <- names(object$model)
  if(k <- match("(weights)", nam, nomatch = F))
    stop("weighted regressions are not supported")
  
  V <- summary(object)$cov.unscaled
  Vi <- crossprod(model.matrix(object))
  nam <- names(coef(object))
  if(k <- match("(Intercept)", nam, nomatch = F)) {
    v1 <- diag(V)[-k]
    v2 <- (diag(Vi)[-k] - Vi[k, -k]^2/Vi[k,k])
    nam <- nam[-k]
  } else {
    v1 <- diag(V)
    v2 <- diag(Vi)
    warning("No intercept term detected. Results may surprise.")
  }
  ##  VIF = structure(v1*v2, names = nam)
  ##  VIF greater than 10 needs to be investigated
  Vi = object$model[,-1]; 
  V = eigen(cor(Vi))$values
  kappa = V[1]/V[length(V)]
  ## eigenvalues suppose to be not small if no presence of
  ##  collinarity.  Zero if perfect multicollinearity. CIndex =
  ##  structure(V,names=nam)
  stat = data.frame(VIF = v1*v2, eigenvalue=V)
  SumEigen = sum(1/V)
#  cat("\n\n (Multi-)Collineary presents if any of VIF > 10;")
#  cat("\n Or if The Condition number is larger than 15;")
#  cat("\n Or if sum(1/eigenvalue) is larger 5p (=",
#      5*length(V), ")\n\n")
  rnames = c("Condition Number","Sum(1/eigenvalue)")
  stat2 = c(kappa, SumEigen)
  Comments = c("Collinear if CN > 15",paste("Collinear if Sum > ", 5*length(V)))
  stats = data.frame(Stat = stat2, Comments = Comments)
  row.names(stats) = rnames
  list(VIF.Eigenvalue = stat, Stat=stats )
} 
