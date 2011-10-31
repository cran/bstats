
scb <- function(x,alpha=0.05){
  name <- deparse(substitute(x))
  if(!is.object(x))stop("'x' must be an R object.")
  if(alpha>1||alpha<0)stop("Invalid confidence/significance level.")
  if(alpha>.5) alpha=1-alpha
  clsname = class(x)
  if(clsname=='edf'){
    epsn = sqrt(.5/x$n*log(2./alpha))
    LFn = x$y-epsn; LFn[LFn<0] = 0
    UFn = x$y+epsn; UFn[UFn>1]=1
  }
  return(structure(list(y=x$y, l=LFn,u=UFn,x=x$x,n = x$n, 
                        data.name = name
                        ), class = "scb"))
}

print.scb <- function (x, digits = NULL, ...) 
{
  cat("\nData: ", x$data.name, 
      " (", x$n, " obs.);", "\n", sep = "")
  print(summary(as.data.frame(x[c("x", "y")])), digits = digits, 
        ...)
  invisible(x)
}

plot.scb  <- 
function (x, main = NULL, xlab = NULL, ylab = "EDF", type = "l", lwd=1,
          col=1,zero.line = TRUE, ...) 
{
  if (is.null(xlab)) 
    xlab <- paste("N =", x$n)
  if (is.null(main)) 
    main <- deparse(x$data.name)
  plot.default(x, type = type, lwd=lwd,
               main = main, xlab = xlab, ylab = ylab, ...)
  mycol = c("#11223311","#33332211","#44223322","#33124255",
    "#77445511","#99663322","#77665533")
  cord.x = c(x$x,rev(x$x));
  cord.y = c(x$l,rev(x$u))
  polygon(cord.x,cord.y,col=mycol[1],border=col,lty=2)#'aliceblue'

  lines(x$x,x$l,col='gray')
  lines(x$x,x$u,col='gray')
  if (zero.line) 
    abline(h = 0, lwd = 0.1, col = "gray")
  invisible(NULL)
}

lines.scb  <- function (x,lwd=1,col=1,lty=1, ...) 
{
  lines(x$x,x$y,col=col,lwd=lwd,lty=lty)
  mycol = c("#11223311","#33332211","#44223322","#33124255",
    "#77445511","#99663322","#77665533")
  cord.x = c(x$x,rev(x$x));
  cord.y = c(x$l,rev(x$u))
  polygon(cord.x,cord.y,col=mycol[1],border=col,lty=2)#'aliceblue'
  lines(x$x,x$l,col='gray')
  lines(x$x,x$u,col='gray')
  invisible(NULL)
}

