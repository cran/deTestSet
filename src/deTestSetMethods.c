#include<R_ext/Rdynload.h>
#ifndef R_R_H
# include <R.h>
#endif

#ifndef R_EXT_DYNLOAD_H_
# include <R_ext/Rdynload.h>
#endif

#include "de.h"

#include <Rinternals.h>
#include <stdlib.h> // for NULL

/*
    .C	R_CMethodDef
    .Call	R_CallMethodDef
    .Fortran	R_FortranMethodDef
    .External	R_ExternalMethodDef
	
	void (*)(double *, double *, double *, double *, double *,
                  double *, int *),              jac 
		     void (*)(double *, double *, double *, double *, double *, int *,
                  double *, int *),               func  
     FCN(N,X,Y,F,RPAR,IPAR)	
     jac(m,t,y,ml,mu,jac,ldjac,rpar,ipar) in bmd	 
*/



/* .Call calls */
extern SEXP call_mebdfi(SEXP, SEXP, SEXP, SEXP, SEXP,
		SEXP, SEXP, SEXP, SEXP, SEXP, 
    SEXP, SEXP, SEXP, SEXP, SEXP, 
    SEXP, SEXP, SEXP, SEXP, SEXP,
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP call_gambim(SEXP, SEXP, SEXP, SEXP, SEXP,
		SEXP, SEXP, SEXP, SEXP, SEXP,
    SEXP, SEXP, SEXP, SEXP, SEXP,
    SEXP, SEXP, SEXP, SEXP, SEXP, 
    SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP call_dop(SEXP, SEXP, SEXP, SEXP, SEXP,
		SEXP, SEXP, SEXP, SEXP, SEXP, 
    SEXP, SEXP, SEXP, SEXP, SEXP, 
    SEXP, SEXP, SEXP);

static const R_CallMethodDef callMethods[]  = {
  {"call_mebdfi", (DL_FUNC) &call_mebdfi, 27},
  {"call_gambim", (DL_FUNC) &call_gambim, 25},
  {"call_dop",    (DL_FUNC) &call_dop   , 18},
  {NULL, NULL, 0}
};


/* .Fortran calls */
void F77_NAME(fekinit)(int *, double *, double *, double *);
void F77_NAME(twobinit)(int *, double *, double *, double *);
void F77_NAME(pleiasoln)(int *, double *);    
void F77_NAME(beamsoln)(int *, double *);    
void F77_NAME(andsoln)(int *, double *);    
void F77_NAME(cranksoln)(int *, double *);    
void F77_NAME(feksoln)(int *, double *);    
void F77_NAME(emepsoln)(int *, double *);    
void F77_NAME(nandsoln)(int *, double *);    
void F77_NAME(polsoln)(int *, double *);    
void F77_NAME(ringsoln)(int *, double *);    
void F77_NAME(tubesoln)(int *, double *);    
void F77_NAME(twobsoln)(int *, double *);    
 
R_FortranMethodDef FEntries[] = {
    {"fekinit",      (DL_FUNC) &F77_SUB(fekinit),      4},
    {"twobinit",     (DL_FUNC) &F77_SUB(twobinit),     4},
    {"pleiasoln",    (DL_FUNC) &F77_SUB(pleiasoln),    2},
    {"beamsoln",     (DL_FUNC) &F77_SUB(beamsoln),     2},
    {"andsoln",      (DL_FUNC) &F77_SUB(andsoln),      2},
    {"cranksoln",    (DL_FUNC) &F77_SUB(cranksoln),    2},
    {"feksoln",      (DL_FUNC) &F77_SUB(feksoln),      2},
    {"emepsoln",     (DL_FUNC) &F77_SUB(emepsoln),     2},
    {"nandsoln",     (DL_FUNC) &F77_SUB(nandsoln),     2},
    {"polsoln",      (DL_FUNC) &F77_SUB(polsoln),      2},
    {"ringsoln",     (DL_FUNC) &F77_SUB(ringsoln),     2},
    {"tubesoln",     (DL_FUNC) &F77_SUB(tubesoln),     2},
    {"twobsoln",     (DL_FUNC) &F77_SUB(twobsoln),     2},
    {NULL, NULL, 0}
};

/* Initialization ---------------------------------------------------------- */
void R_init_deTestSet(DllInfo *dll) {

  // register entry points
  R_registerRoutines(dll, NULL, callMethods, FEntries, NULL);

  R_useDynamicSymbols(dll, TRUE); // enable dynamic searching
} 
