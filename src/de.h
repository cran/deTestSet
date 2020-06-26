#ifndef R_R_H
#  include <R.h>
#endif

#ifndef R_DEFINES_H
#  include <Rdefines.h>
#endif

#ifndef R_INTERNALS_H_
#  include <Rinternals.h>
#endif


/*============================================================================
  global R variables 
============================================================================*/

#ifndef EXTERN
# define EXTERN extern
#endif

/*============================================================================
  global R variables 
============================================================================*/
EXTERN SEXP YOUT, YOUT2, ISTATE, RWORK, IROOT;    /* returned to R */
EXTERN SEXP Y, YPRIME , Rin;

/*============================================================================
  global C variables 
============================================================================*/
EXTERN int    it, n_eq; 
EXTERN int    *iwork;   
EXTERN double *rwork;
EXTERN double *xytmp, *xdytmp, tin, tout;

EXTERN long int mu;
EXTERN long int ml;
EXTERN long int nrowpd;

/* output in DLL globals */
EXTERN int nout, ntot, isOut, lrpar, lipar, *ipar;
EXTERN double *out;

/* forcings  */
EXTERN long int nforc;  /* the number of forcings */
EXTERN double *tvec;
EXTERN double *fvec;
EXTERN int    *ivec;
EXTERN int    fmethod;
EXTERN int    finit;

EXTERN int    *findex;
EXTERN double *intpol;
EXTERN int    *maxindex;

EXTERN double *forcings;

/* events - not yet used */
EXTERN double tEvent;
EXTERN int iEvent, nEvent, typeevent, rootevent;

EXTERN double *timeevent, *valueevent;
EXTERN int *svarevent, *methodevent;

/*============================================================================
 type definitions for C functions
============================================================================*/
typedef void C_deriv_func_type(int*, double*, double*, double*, double*, int*);
EXTERN C_deriv_func_type *DLL_deriv_func;

EXTERN  C_deriv_func_type *deriv_func;  

typedef void C_res_func_type(double*, double*, double*, double*, double*,
                             int*, double*, int*);
EXTERN C_res_func_type* DLL_res_func;


/* this is in compiled code */
typedef void init_func_type (void (*)(int*, double*));
/* C_deriv_func_type *deriv_func;*/ 

EXTERN double             *tt, *ytmp;
EXTERN int                isDll;

EXTERN  int maxt;

/*============================================================================
  solver R- global functions 
============================================================================*/

/* DAE globals */
EXTERN  SEXP R_res_func;
EXTERN  SEXP R_daejac_func;
EXTERN  SEXP R_deriv_func;
EXTERN  SEXP R_jac_func;
EXTERN  SEXP R_mas_func;
EXTERN  SEXP R_envir;
EXTERN  SEXP R_event_func;

EXTERN  SEXP de_gparms;
EXTERN SEXP getListElement(SEXP list, const char* str);

/*============================================================================ 
  C- utilities, functions 
============================================================================
void init_N_Protect(void);
void incr_N_Protect(void);
void unprotect_all(void);
void my_unprotect(int);
void returnearly (int);
void terminate(int, int, int, int, int);
void initParms(SEXP Initfunc, SEXP Parms);*/

/* declarations for initialisations */
void Initdeparms(int*, double*);
void Initdeforc(int*, double*);
void initOutC(int isDll, int neq, SEXP nOut, SEXP Rpar, SEXP Ipar);

/* not yet used in mebdfi... */
void initglobals(int);
void initdaeglobals(int);

/* the forcings and event functions - latter not yet implemented */
void updatedeforc(double*);
int initForcings(SEXP list);
int initEvents(SEXP list, SEXP);
void updateevent(double*, double*, int*);



/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                         DECLARATIONS for time lags
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*==========================================
  R-functions
==========================================*/
EXTERN int noff;

EXTERN SEXP getPastValue   (SEXP T, SEXP nr);
EXTERN SEXP getPastGradient(SEXP T, SEXP nr);

/*==========================================
  C- utilities, functions
==========================================*/
/* Hermitian interpolation */
double Hermite (double t0, double t1, double y0, double y1, double dy0,
                double dy1, double t);

double dHermite(double t0, double t1, double y0, double y1, double dy0,
                double dy1, double t);

int initLags(SEXP elag, int solver, int nroot);

/* history vectors  */
void inithist(int max, int maxlags, int solver, int nroot);

void updatehistini(double t, double *y, double *dY, double *rwork, int *iwork);
void updatehist(double t, double *y, double *dy, double *rwork, int *iwork);

int nexthist(int i);
double interpolate(int i, int k, double t0, double t1, double t,
  double *Yh, int nq);


/*==========================================
  Global variables for history arrays
==========================================*/
/* time delays */
EXTERN int interpolMethod;  /* for time-delays : 1 = hermite; 2=dense */

EXTERN int indexhist, indexlag, endreached, starthist;
EXTERN double *histvar, *histdvar, *histtime, *histhh, *histsave;
EXTERN int    *histord;
EXTERN int    histsize, offset;
EXTERN int    initialisehist, lyh, lhh, lo;

#undef EXTERN
