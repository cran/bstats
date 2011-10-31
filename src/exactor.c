#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <stdio.h>
#include <math.h>
#include "R_ext/Applic.h"

// ISNA(x): true for R's NA only
// ISNAN(x): true for R's NA and IEEE NaN
// R_FINITE(x): false for Inf,-Inf,NA,NaN
// R_IsNaN(x): true for NaN but not NA

// Rprintf:  printing from a C routine compiled into R

//recursive version to compute factorial of n := n!
int factorial(int n){
  return n>=1 ? n * factorial(n-1) : 1;
}

double g1(double p, int m1, int n11, double a[], double alpha){
  int i;
  double p1=0.0, p2=0.0;
  for(i=0;i<n11;i++)
    p1 += a[i] * pow(p, i);
  for(i=n11;i<m1+1;i++){
    p1 += a[i] * pow(p, i);
    p2 += a[i] * pow(p, i);
  }
  return p2/p1 - 0.5 * alpha;
}

double g2(double p, int m1, int n11, double a[], double alpha){
  int i;
  double p1=0.0, p2=0.0;
  for(i=0;i<n11+1;i++){
    p1 += a[i] * pow(p, i);
    p2 += a[i] * pow(p, i);
  }
  for(i=n11+1;i<m1+1;i++){
    p1 += a[i] * pow(p, i);
  }
  return p2/p1 - 0.5 * alpha;
}

double dg1(double p, int m1, int n11, double a[]){
  int i;
  double p1, p2=0.0, dp1=0.0, dp2=0.0;

  p1 = a[0];
  for(i=1;i<n11;i++){
    p1 += a[i] * pow(p, i);
    dp1 += i * a[i] * pow(p, i-1);
  }
  for(i=n11;i<m1+1;i++){
    p1 += a[i] * pow(p, i);
    dp1 += i * a[i] * pow(p, i-1);
    p2 += a[i] * pow(p, i);
    dp2 += i * a[i] * pow(p, i-1);
  }
  return (dp2 * p1 - p2 * dp1)/(p1 * p1);
}

double dg2(double p, int m1, int n11, double a[]){
  int i;
  double p1, p2=0.0, dp1=0.0, dp2=0.0;

  p1 = a[0];
  p2 = a[0];
  for(i=1;i<n11+1;i++){
    p1 += a[i] * pow(p, i);
    dp1 += i * a[i] * pow(p, i-1);
    dp2 += i * a[i] * pow(p, i-1);
  }
  for(i=n11+1;i<m1+1;i++){
    p1 += a[i] * pow(p, i);
    dp1 += i * a[i] * pow(p, i-1);
  }
  return (dp2 * p1 - p2 * dp1)/(p1 * p1);
}

/* use out to pass the initial value and the final estimate

   The Newton method was used, fast but not stable.  We change to the
   bisection method instead.  We take small value a = 1e-16 and b=1e6.
   If not in the range, simply use 0 or 100000+
 */

void orexactl(int *counts, double *alpha, double *out)
{
  int i,n11,n12,n21,n22, n1, n2, m1;
  n11 = counts[0];
  n12 = counts[1];
  n21 = counts[2];
  n22 = counts[3];
  n1 = n11+n12;
  n2 = n21+n22;
  m1 = n11+n21;
  double delta = 1.0, p0 = out[0],f0, fa, fb,pa,pb;
  double a[m1 + 1]; // to store the coefficients
  

  for(i=0;i < m1+1;i++)
    a[i] = choose(n1, i) * choose(n2, m1 - i);
  
  i = 0;
  pa = 1.e-16; pb = 1.e6;
  f0 = g1(p0, m1, n11, a, alpha[0]);
  if(f0>0.) pb = p0; else pa = p0;
  fa = g1(pa, m1, n11, a, alpha[0]);
  fb = g1(pb, m1, n11, a, alpha[0]);   

  while(i<10000 && fabs(delta) >0.00001){
    if(fa >=0.){
      p0 = pa; break; //exit
    }else if(fb<=0.){
      p0 = pb; break;
    }else{
      p0 = 0.5 * (pa + pb);      
      f0 = g1(p0, m1, n11, a, alpha[0]);
      if(f0>0.){
	pb = p0; 
	fb = g1(pb, m1, n11, a, alpha[0]);   
      }else{
	pa = p0;
	fa = g1(pa, m1, n11, a, alpha[0]);
      }
      i++;
    }
  }

  out[0] = p0;
}

void orexactu(int *counts, double *alpha, double *out)
{
  int i,n11,n12,n21,n22, n1, n2, m1;
  n11 = counts[0];
  n12 = counts[1];
  n21 = counts[2];
  n22 = counts[3];
  n1 = n11+n12;
  n2 = n21+n22;
  m1 = n11+n21;
  double delta = 1.0, p0 = out[0],f0, fa, fb,pa,pb;
  double a[m1 + 1]; // to store the coefficients
  

  for(i=0;i < m1+1;i++)
    a[i] = choose(n1, i) * choose(n2, m1 - i);
  
  i = 0;
  pa = 1.e-16; pb = 1.e6;
  f0 = g2(p0, m1, n11, a, alpha[0]);
  if(f0 < 0.) pb = p0; else pa = p0;
  fa = g2(pa, m1, n11, a, alpha[0]);
  fb = g2(pb, m1, n11, a, alpha[0]);   

  while(i<10000 && fabs(delta) >0.00001){
    if(fa <=0.){
      p0 = pa; break; //exit
    }else if(fb>=0.){
      p0 = pb; break;
    }else{
      p0 = 0.5 * (pa + pb);      
      f0 = g2(p0, m1, n11, a, alpha[0]);
      if(f0<0.){
	pb = p0; 
	fb = g2(pb, m1, n11, a, alpha[0]);   
      }else{
	pa = p0;
	fa = g2(pa, m1, n11, a, alpha[0]);
      }
      i++;
    }
  }
  //  p0 = g2(10., m1, n11, a, alpha[0]);

  out[0] = p0;
}
