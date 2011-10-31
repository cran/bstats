## to compute the confidence intervals of the regression parameters

lm.ci <- function(lmobj,level=0.95)
  {
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
    fmout = anova(fmobj)
    rmout = anova(rmobj)
    fmdf = fmout[nrow(fmout),1]
    rmdf = rmout[nrow(rmout),1]
    fmsse = fmout[nrow(fmout),2]
    rmsse = rmout[nrow(rmout),2]
    if(rmsse<fmsse){
      tmp = rmsse; rmsse=fmsse; fmsse=tmp;
      tmp = rmdf; rmdf=fmdf; fmdf=tmp;
    }
    F = (rmsse-fmsse)/(rmdf-fmdf)/fmsse*fmdf
    pv = 1-pf(F,rmdf-fmdf,fmdf)
    cv = qf(1-alpha,rmdf-fmdf,fmdf)
    c("Test Stat."=F, "p-value"=pv, "Critical value"=cv)
  }


influential.plot <- function(lmobj,type='hadi',ID=FALSE){
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
  out = cbind(Ri = stdres, Leverage=h,Hadi = Hi, DFIT=DFIT,CookD=CookD)
  if(type=="hadi"){
    Index = 1:length(di2)
    plot(Hi~Index)
    if(ID)    identify(Hi~Index)
  }else if(type=='potential-residual'){
    di2 = (lmobj$res)^2/tmp[k,2]
    Y = h/(1-h)
    X =  (p+1)/(1-h)*di2/(1-di2)
    plot(Y~X, ylab="Potential",xlab="Residual")
    if(ID)    identify(Y~X)
  }else if(type == "leverage" || type == "hat"){
    plot(h,type='h')
    tmp = 2*(p+1)/n
    abline(h=tmp, col='gray')
  }else if(type == "dfits"){
    tmp = 2*sqrt((p+1)/(n-p-1)) 
    plot(DFIT,type='h')
    abline(h=c(-tmp,tmp), col='gray')
  }else if(type == "cook"){
    plot(CookD,type='h')
    abline(h=1, col='gray')
  }
  invisible(round(out,3))
}



residual.plot <- function(lmobj,type='fitted'){
  type = match.arg(tolower(type),c("fitted","index","predictor","qqplot"))
  stdres = rstandard(lmobj)
  if(type=="fitted"){
    fittedvalue = lmobj$fitted
    par(mfrow=c(1,1))
    plot(stdres~fittedvalue,ylab="Standardized Residuals",xlab="Fitted values")
    abline(h=c(-2,2),col="gray")
  }else if(type=="index"){
    Index = 1:nrow(lmobj$mode)
    par(mfrow=c(1,1))
    plot(stdres~Index,ylab="Standardized Residuals",xlab="Index")
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
      plot(stdres~x,ylab="Standardized Residuals",xlab=xnames[l])
      abline(h=c(-2,2),col="gray")
    }
  }else if(type=='qqplot'){
    par(mfrow=c(1,1))
    qqnorm(stdres, ylab="Standardized Residuals",
           xlab="Normal Scores", main=lmobj$call)
    qqline(stdres)
  }else stop("'type' not supported!")
}

predictor.plot <- function(lmobj,type='av',ID=FALSE){
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
      
      plot(yres~xres,
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
        
      plot(yres~xres,
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
  structure(v1*v2, names = nam)
} 
