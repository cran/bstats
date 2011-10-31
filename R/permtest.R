### The following functions are to compute the empirical distribution
### function and construct a confidence band.
.edfperm <- function(m,x,x0){
  n = length(x)
  n0 = sample(n)
  x1 = x[n0[1:m]];
  x2 = x[n0[(m+1):n]];
  Fx1 = edf(x1,x0); Fx2 = edf(x2,x0);
  max(abs(Fx1$y-Fx2$y))
}

perm.test <- function(x, INDEX,iter=10000){
  if(any(is.na(x)))stop("'x' contains missing value(s).");
  name <- deparse(substitute(x))
  n = length(x); 
  z = as.factor(INDEX)
  if(length(levels(z))!=2)stop("'INDEX' should have two levels.")
  if(length(INDEX)!= n)stop("'x' and 'INDEX' have different length!")
  xmin = min(x); xmax = max(x)
  sele = INDEX == levels(z)[1]
  m = sum(sele);
  x0 = seq(xmin,xmax,length=256)
  Fx1 = edf(x[sele],x0); Fx2 = edf(x[!sele],x0);
  D = max(abs(Fx1$y-Fx2$y))
  z = apply(as.matrix(rep(m,iter),ncol=1),1,.edfperm,x=x,x0=x0)
  pv = mean(z>D)
  c(D=D,p.value=pv,iter=iter)
}


