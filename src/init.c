/*
 *  Part of R package Genome
 *  Copyright (C) 2009-2010  B. Wang
 *
 *  Unlimited use and distribution (see LICENCE).
 */

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void F77_SUB(pan)(double *A, int *M, double *C, int *N, double *RESULT);
void F77_SUB(nrlogit)(double *x0, double *betas, double *ps, int *n);

static const R_FortranMethodDef FortEntries[] = {
  {"pan", (DL_FUNC) &F77_SUB(pan), 5},
  {"nrlogit", (DL_FUNC) &F77_SUB(nrlogit),  4},
  {NULL, NULL, 0}
};


void R_init_bstats(DllInfo *dll)
{
  //    R_registerRoutines(dll, NULL, NULL, callMethods, NULL);
  R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
