/* Define some global variables and functions that operate on some of them */

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#include "de.h"

/*==================================================
some functions for keeping track of how many SEXPs
are PROTECTed, and UNPROTECTing them in the
case of a FORTRAN stop.
==================================================*/
 
long int N_Protected;

void init_N_Protect(void) { N_Protected = 0; }

void incr_N_Protect(void) { N_Protected++; }

void unprotect_all(void) { UNPROTECT((int) N_Protected); }

void my_unprotect(int n) {
    UNPROTECT(n);
    N_Protected -= n;
}

/* Globals :*/
SEXP R_deriv_func;
SEXP R_jac_func;
SEXP R_event_func;

SEXP R_envir;

SEXP R_res_func;
SEXP R_mas_func;
SEXP R_daejac_func;

SEXP de_gparms;

/*======================================================
SEXP initialisation functions
=======================================================*/

void initglobals(int nt) {

  PROTECT(Time = NEW_NUMERIC(1));                  incr_N_Protect();
  PROTECT(Y = allocVector(REALSXP,(n_eq)));        incr_N_Protect();
  PROTECT(YOUT = allocMatrix(REALSXP,ntot+1,nt));  incr_N_Protect();
}

void initdaeglobals(int nt) {
  PROTECT(Time = NEW_NUMERIC(1));                    incr_N_Protect();
  PROTECT(Rin  = NEW_NUMERIC(2));                    incr_N_Protect();
  PROTECT(Y = allocVector(REALSXP,n_eq));            incr_N_Protect();
  PROTECT(YPRIME = allocVector(REALSXP,n_eq));       incr_N_Protect();
  PROTECT(YOUT = allocMatrix(REALSXP,ntot+1,nt));    incr_N_Protect();
}

/*======================================================
Parameter initialisation functions
note: forcing initialisation function is in forcings.c
=======================================================*/

void initParms(SEXP Initfunc, SEXP Parms) {
  if (Initfunc == NA_STRING) return;
  if (inherits(Initfunc, "NativeSymbol")) {
    init_func_type *initializer;
    PROTECT(de_gparms = Parms); incr_N_Protect();
    initializer = (init_func_type *) R_ExternalPtrAddr(Initfunc);
    initializer(Initdeparms);
  }
}


void Initdeparms(int *N, double *parms) {
  int i, Nparms;
  Nparms = LENGTH(de_gparms);
  if ((*N) != Nparms) {
    warning("Number of parameters passed to solver, %i; number in DLL, %i\n",
      Nparms, *N);
    PROBLEM "Confusion over the length of parms"
    ERROR;
  } else {
    for (i = 0; i < *N; i++) parms[i] = REAL(de_gparms)[i];
  }
}
  
SEXP get_de_gparms(void) {
  return de_gparms;
}

/*==================================================
 Termination 
===================================================*/

/* an error occurred - save output in YOUT2 */
void returnearly (int Print) {
  int j, k;
  if (Print) 
    warning("Returning early. Results are accurate, as far as they go\n");
  PROTECT(YOUT2 = allocMatrix(REALSXP,ntot+1,it)); incr_N_Protect();
  for (k = 0; k < it; k++)
    for (j = 0; j < ntot+1; j++)
      REAL(YOUT2)[k*(ntot+1) + j] = REAL(YOUT)[k*(ntot+1) + j];
}   

/* add ISTATE and RSTATE */
void terminate(int istate, int ilen, int ioffset, int rlen, int roffset) {

  int k;
  
  PROTECT(ISTATE = allocVector(INTSXP, ilen)); incr_N_Protect();
  for (k = 0; k < ilen-1; k++) INTEGER(ISTATE)[k+1] = iwork[k +ioffset];
  INTEGER(ISTATE)[0] = istate;  
        
  PROTECT(RWORK = allocVector(REALSXP, rlen)); incr_N_Protect();
  for (k = 0; k < rlen; k++) REAL(RWORK)[k] = rwork[k+roffset];
  if (istate > 0) {
    setAttrib(YOUT, install("istate"), ISTATE);
    setAttrib(YOUT, install("rstate"), RWORK);
  }
  else  {
    setAttrib(YOUT2, install("istate"), ISTATE);
    setAttrib(YOUT2, install("rstate"), RWORK);
  }
}

/*==================================================
 extracting elements from a list
===================================================*/

SEXP getListElement(SEXP list, const char *str) {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;
  for (i = 0; i < length(list); i++)
    if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}

/*==================================================
 output initialisation function

 out and ipar are used to pass output variables
 (number set by nout) followed by other input
 by R-arguments rpar, ipar
 ipar[0]: number of output variables,
 ipar[1]: length of rpar,
 ipar[2]: length of ipar
===================================================*/

/* Initialise output - output variables calculated in C-code ... */

void initOutC(int isDll, int neq, SEXP nOut, SEXP Rpar, SEXP Ipar) {
  int j;
  /* initialise output when a dae ... */   
  /*  output always done here in C-code (<-> lsode, vode)... */

  nout  = INTEGER(nOut)[0];
  ntot  = n_eq+nout;
  
  if (isDll == 1) {                /* function is a dll */
    lrpar = nout + LENGTH(Rpar);   /* length of rpar */
    lipar = 3    + LENGTH(Ipar);   /* length of ipar */
  } else {                         /* function is not a dll */
    lipar = 3;
    lrpar = nout;
  }
  out   = (double*) R_alloc(lrpar, sizeof(double));
  ipar  = (int*)    R_alloc(lipar, sizeof(int));

  if (isDll == 1) {
    ipar[0] = nout;                /* first 3 elements of ipar are special */
    ipar[1] = lrpar;
    ipar[2] = lipar;

    /* other elements of ipar are set in R-function lsodx via argument *ipar* */
    for (j = 0; j < LENGTH(Ipar); j++) ipar[j+3] = INTEGER(Ipar)[j];

    /* first nout elements of rpar reserved for output variables
       other elements are set in R-function lsodx via argument *rpar* */
    for (j = 0; j < nout;         j++) out[j] = 0.;
    for (j = 0; j < LENGTH(Rpar); j++) out[nout+j] = REAL(Rpar)[j];
   }
}
