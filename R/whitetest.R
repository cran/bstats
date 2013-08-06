##  Testing constant variance  (http://en.wikipedia.org/wiki/White_test)

## To test for constant variance one undertakes an auxiliary
## regression analysis: this regresses the squared residuals from the
## original regression model onto a set of regressors that contain the
## original regressors, the cross-products of the regressors and the
## squared regressors. One then inspects the R2. The LM test statistic
## is the product of the R2 value and sample size:

##  LM = n*R^2

##  This follows a chi-square distribution, with degrees of freedom
##  equal to the number of estimated parameters (in the auxiliary
##  regression) minus one.

##1. White, H. (1980). A Heteroskedasticity-Consistent Covariance
## Matrix Estimator and a Direct Test for
## Heteroskedasticity. Econometrica 48 (4): 817–838. JSTOR
## 1912934. MR575027.

## 2. Kim, E.H.; Morse, A.; Zingales, L. (2006). What Has Mattered to
## Economics since 1970. Journal of Economic Perspectives 20 (4):
## 189–202. doi:10.1257/jep.20.4.189.


white.test <- function(lmobj)
{
  stopifnot(class(lmobj)=='lm')
  mydata = lmobj$model
  mydata[,1] = lmobj$residual^2
  fml = lmobj$call$formula
  formula1 = paste(fml[2],fml[1],fml[3])
  pvs = attr(lmobj$terms,"term.labels")
  k = length(pvs)
  for(i in 1:k){
    tmp <- NULL
    for(j in 1:nchar(pvs[i])){
      tmp1 <- substr(pvs[i],j,j)
      if(tmp1 == ":")
        tmp <- paste(tmp, "*", sep='')
      else
        tmp <- paste(tmp, tmp1, sep='')
    }
    pvs[i] <- tmp
  }
  formula2 <- paste(fml[2],fml[1])
  for(i in 1:k){
    if(i>1)
      formula2 <- paste(formula2, "+", sep='')
    formula2 <- paste(formula2, "I(", pvs[i],")",sep='')

    for(j in i:k)
      formula2 = paste(formula2,"+I(",pvs[i],"*",pvs[j],")", sep='')
  }
  ##  print(formula1)
  out = lm(as.formula(formula2),data=mydata)
  n = length(lmobj$fit)
  LM = summary(out)$r.squared * n  
  names(LM) <- "White"
  df <- out$rank - 1
  names(df) <- "df";
  RVAL <- list(statistic = LM,
      parameter = df,
      method = "White test for constant variance",
      p.value= pchisq(LM,df,lower.tail=FALSE),
      data.name=NULL)

  class(RVAL) <- "htest"
  return(RVAL)
}

