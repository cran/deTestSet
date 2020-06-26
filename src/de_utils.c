/* Define some global variables and functions that operate on some of them */

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#include "de.h"
#include "externalptr.h"

/*==================================================
some functions for keeping track of how many SEXPs
are PROTECTed, and UNPROTECTing them in the
case of a FORTRAN stop. REDUNDANT NOW
==================================================
 
long int N_Protected;
void init_N_Protect(void) { N_Protected = 0; }
void incr_N_Protect(void) { N_Protected++; }
void unprotect_all(void) { UNPROTECT((int) N_Protected); }
void my_unprotect(int n) {
    UNPROTECT(n);
    N_Protected -= n;
}
*/


/* Globals :*/
//SEXP R_deriv_func;
//SEXP R_jac_func;
//SEXP R_event_func;

//SEXP R_envir;

//SEXP R_res_func;
//SEXP R_mas_func;
//SEXP R_daejac_func;

//SEXP de_gparms;

/*======================================================
SEXP initialisation functions - removed - done in solver functions
 Used to be in initglobals and initdaeglobals
=======================================================*/

/*======================================================
Parameter initialisation is done in the solver function
 (used to be in initParms)
=======================================================*/

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
 Termination - Redundant now done in solvers
 used to be returnearly and terminate
===================================================*/

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

