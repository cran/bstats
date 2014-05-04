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

## 2014/05/04: add an option to white.test per Oleg's suggestion (Oleg
## is a Professor teaching econometrics and forecasting at Kiev
## Shevchenko university oleg_komashko@ukr.net.  THANKS.)


white.test <- function(lmobj, squares.only=FALSE)
{
    stopifnot(class(lmobj)=='lm')
    mydata <- lmobj$model
    mydata[,1] <- lmobj$residual^2
    fml <- lmobj$call$formula
    formula1 <- paste(fml[2],fml[1],fml[3])
    pvs <- attr(lmobj$terms,"term.labels")
    k <- length(pvs);
    n <- length(lmobj$fit)
    
    for(i in 1:k){
        tmp <- NULL;
        if(substr(pvs[i],1,2)=="I("){
            tmp2 <- substr(pvs[i],3, nchar(pvs[i])-1);
        }else{
            tmp2 <- pvs[i];
        }
        for(j in 1:nchar(tmp2)){
            tmp1 <- substr(tmp2,j,j)
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
        if(squares.only){
            formula2 <- paste(formula2, "+I(", pvs[i], 
                              "*", pvs[i], ")", sep = "")
        }else{
            for(j in i:k)
                formula2 <- paste(formula2,"+I(",pvs[i],
                                  "*",pvs[j],")", sep='')
        }
    }

    method <- ifelse(squares.only,
                     "White test for constant variance, squares only",
                     "White test for constant variance")

    out <- lm(as.formula(formula2),data=mydata)
    if(summary(out)$r.squared == 1.0){
        RVAL <- NULL;
        warning("Test failed.  Possible reasons:\n\t (1) collinearity, or (2) sample size is not big enough for the White's test.");
    }else{
        LM = summary(out)$r.squared * n  
        names(LM) <- "White"
        df <- out$rank - 1
        names(df) <- "df";
        RVAL <- list(statistic = LM,
                     parameter = df,
                     method = method,
                     p.value= pchisq(LM,df,lower.tail=FALSE),
                     data.name=NULL)
        class(RVAL) <- "htest"
    }
    return(RVAL)
}

