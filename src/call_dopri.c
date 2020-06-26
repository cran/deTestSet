#include <time.h>
#include <string.h>
#include "de.h"
#include "externalptr.h"

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Explicit runge-Kutta of order 8(5,3)
   due to Dormand and Prince, with stepsize control and dense output
   
   The C-wrappers that provide the interface between FORTRAN codes and R-code 
   are: C_deriv_func_dop: interface with R-code "func", passes derivatives  
        C_deriv_out_dop : interface with R-code "func", passes derivatives + output variables  
  
   C_deriv_func_forc_dop provides the interface between the function specified in
   a DLL and the integrator, in case there are forcing functions.
   
 karline soetaert
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*  some globals */
 int type ;     /* 1 = dopri 8, 2 = dopri5, 3 = cashkarp */
 int lrc;

/* definition of the calls to the FORTRAN subroutines in file cash.f, dopri853 and 
   dopri5.f */
		     
void F77_NAME(cashkarp)( int *,
         void (*)(int *, double *, double *, double *,
                              double *, int *),         // func
		     double *, double *, double *, double *,
		     double *, int *,  
 	       void (*)(int *, double *, double *, double *, int *, double *,  
			            int *, int *, double *, int *, int *, double *),   // solout
		     int *, double *, int *, int *, int*, double *, int*, int*, double*);

void F77_NAME(dopri6)( int *,
         void (*)(int *, double *, double *, double *,
                              double *, int *),         // func
		     double *, double *, double *, double *,
		     double *, int *,  
 	       void (*)(int *, double *, double *, double *, int *, double *,  
			            int *, int *, double *, int *, int *, double *),   // solout
		     int *, double *, int *, int *, int*, double *, int*, int*, double*);

void F77_NAME(dop853)( int *,
         void (*)(int *, double *, double *, double *,
                              double *, int *),         // func
		     double *, double *, double *, double *,
		     double *, int *,  
 	       void (*)(int *, double *, double *, double *, int *, double *,  
			            int *, int *, double *, int *, int *, double *),   // solout
		     int *, double *, int *, int *, int*, double *, int*, int*);

/* continuous output formula */
void F77_NAME(contd5)(int *, double *, double *,int *, int *, double *);
void F77_NAME(contd5ck)(int *, double *, double *,int *, int *, double *);
void F77_NAME(contd8)(int *, double *, double *,int *, int *, double *);

/* wrapper above the derivate function in a dll that first estimates the
values of the forcing functions */

static void C_deriv_func_forc_dop (int *neq, double *t, double *y,
                         double *ydot, double *yout, int *iout)
{
  updatedeforc(t);
  DLL_deriv_func(neq, t, y, ydot, yout, iout);
}

/* interface between FORTRAN function call and R function
   Fortran code calls C_deriv_func_dop(N, t, y, ydot, yout, iout) 
   R code called as R_deriv_func(time, y) and returns ydot */

static void C_deriv_func_dop (int *neq, double *t, double *y,
                          double *ydot, double *yout, int *iout)
{
  int i;
  SEXP R_fcall, Time, ans;

  for (i = 0; i < *neq; i++)  REAL(Y)[i] = y[i];

  PROTECT(Time = ScalarReal(*t));                  
  PROTECT(R_fcall = lang3(R_deriv_func,Time,Y));   
  PROTECT(ans = eval(R_fcall, R_envir));           

  for (i = 0; i < *neq; i++)   ydot[i] = REAL(ans)[i];

  UNPROTECT(3);
}

/* deriv output function - for ordinary output variables */

static void C_deriv_out_dop (int *nOut, double *t, double *y, 
                       double *ydot, double *yout)
{
  int i;
  SEXP R_fcall, Time, ans;
  
  for (i = 0; i < n_eq; i++)  
      REAL(Y)[i] = y[i];
     
  PROTECT(Time = ScalarReal(*t));                   
  PROTECT(R_fcall = lang3(R_deriv_func,Time, Y));   
  PROTECT(ans = eval(R_fcall, R_envir));            

  for (i = 0; i < *nOut; i++) yout[i] = REAL(ans)[i + n_eq];

  UNPROTECT(3);                                  
}      

/* save output in R-variables */

static void saveOut (double t, double *y) {
  int j;
  
    REAL(YOUT)[(it)*(ntot+1)] = t;
	  for (j = 0; j < n_eq; j++)
	    REAL(YOUT)[(it)*(ntot + 1) + j + 1] = y[j];

    /* if ordinary output variables: call function again */
    if (nout>0)   {
      if (isDll == 1)   /* output function in DLL */
        deriv_func (&n_eq, &t, y, xdytmp, out, ipar) ;
      else
        C_deriv_out_dop(&nout, &t, y, xdytmp, out);  
      for (j = 0; j < nout; j++) 
        REAL(YOUT)[(it)*(ntot + 1) + j + n_eq + 1] = out[j];
    }                
}

/* function called by Fortran to check for output */
static void C_solout(int * nr, double * told, double *t, double * y, int * neq, 
  double * con, int *icomp, int * nd, double * rpar, int * ipar, int * irtrn, double *xout)
{
  if (*told == *t) return;

  while (*told <= tt[it] && tt[it] < *t) {
    if (type == 1)      F77_CALL(contd8)   (neq, &tt[it], con, icomp, nd, ytmp);
    else if (type == 2) F77_CALL(contd5)   (neq, &tt[it], con, icomp, nd, ytmp);
    else if (type == 3) F77_CALL(contd5ck) (neq, &tt[it], con, icomp, nd, ytmp);
    saveOut(tt[it], ytmp);	     
    it++;
  }
}

/* give name to data types */
typedef void C_solout_type (int *, double *, double *, double *,
  int *, double *, int *, int *, double *, int *, int *, double *) ;


/* MAIN C-FUNCTION, CALLED FROM R-code */

SEXP call_dop(SEXP y, SEXP times, SEXP derivfunc, SEXP parms, SEXP rtol,
		SEXP atol, SEXP rho, SEXP initfunc, 
    SEXP verbose, SEXP rWork, SEXP iWork, 
    SEXP nOut, SEXP lRw, SEXP lIw, 
    SEXP Rpar, SEXP Ipar, SEXP flist, SEXP Type)

{
/******************************************************************************/
/******                   DECLARATION SECTION                            ******/
/******************************************************************************/

  int  j, k, nt, latol, lrtol, lrw, liw;
  int  isForcing;
  double *Atol, *Rtol, *ww;
  int itol, iout, idid, mflag;

  /* pointers to functions passed to FORTRAN */
  C_solout_type         *solout = NULL;
  
/******************************************************************************/
/******                         STATEMENTS                               ******/
/******************************************************************************/

/*                      #### initialisation ####                              */    
  int nprot = 0;

  n_eq = LENGTH(y);             /* number of equations */ 
  nt   = LENGTH(times);         /* number of output times */
  
  tt = (double *) R_alloc(nt, sizeof(double));
  for (j = 0; j < nt; j++) tt[j] = REAL(times)[j];
  
  mflag = INTEGER(verbose)[0];
  type  = INTEGER(Type)[0];     /* 1 = dopri 8, 2 = dopri5, 3 = cashkarp */
  lrc = 4;
  if (type == 1) lrc = 8;
  
  /* is function a dll ?*/
  isDll = inherits(derivfunc, "NativeSymbol");

  /* initialise output ... */
  initOutC(isDll, n_eq, nOut, Rpar, Ipar);

  /* copies of variables that will be changed in the FORTRAN subroutine */
  xytmp = (double *) R_alloc(n_eq, sizeof(double));
  for (j = 0; j < n_eq; j++) xytmp[j] = REAL(y)[j];

  ytmp = (double *) R_alloc(n_eq, sizeof(double));
 
  latol = LENGTH(atol);
  Atol = (double *) R_alloc((int) latol, sizeof(double));
  for (j = 0; j < latol; j++) Atol[j] = REAL(atol)[j];

  lrtol = LENGTH(rtol);
  Rtol = (double *) R_alloc((int) lrtol, sizeof(double));
  for (j = 0; j < lrtol; j++) Rtol[j] = REAL(rtol)[j];

  /* tolerance specifications */
  if (latol == 1 ) itol = 0;
  else             itol = 1;
  
  /* work vectors */
  ww = (double *) R_alloc(n_eq, sizeof(double));
  for (j = 0; j < n_eq; j++) ww[j] = 0.;

  liw = INTEGER (lIw)[0];
  iwork = (int *) R_alloc(liw, sizeof(int));
  for (j=0; j<LENGTH(iWork); j++) iwork[j] = INTEGER(iWork)[j];
  for (j=LENGTH(iWork); j<liw; j++) iwork[j] = 0;
  iwork[2] = mflag;

  lrw = INTEGER(lRw)[0];
  rwork = (double *) R_alloc(lrw, sizeof(double));
  for (j=0; j<length(rWork); j++) rwork[j] = REAL(rWork)[j];
  for (j=length(rWork); j<lrw; j++) rwork[j] = 0.;

  /* initialise global R-variables... initglobals (nt); */
  
  PROTECT(Y = allocVector(REALSXP,(n_eq)));                nprot++;
  PROTECT(YOUT = allocMatrix(REALSXP,ntot+1,nt));          nprot++;
  
  /* Initialization of Parameters, Forcings (DLL) initParms (initfunc, parms);*/
  if (initfunc != NA_STRING) {
    if (inherits(initfunc, "NativeSymbol")) {
      init_func_type *initializer;
      PROTECT(de_gparms = parms);                          nprot++;
      initializer = (init_func_type *) R_ExternalPtrAddrFn_(initfunc);
      initializer(Initdeparms);
    }
  }

  isForcing = initForcings(flist);

  if (nout > 0 ) {
     xdytmp= (double *) R_alloc(n_eq, sizeof(double));
     for (j = 0; j < n_eq; j++) xdytmp[j] = 0.; 
  }
  
 /* pointers to functions deriv_func, jac_func, jac_vec, root_func, passed to FORTRAN */
  if (isDll)  { /* DLL address passed to FORTRAN */
      deriv_func = (C_deriv_func_type *) R_ExternalPtrAddrFn_(derivfunc);  
	  
 	   /* overruling deriv_func if forcing */
      if (isForcing) {
        DLL_deriv_func = deriv_func;
        deriv_func = (C_deriv_func_type *) C_deriv_func_forc_dop;
      }
  } else {
      /* interface function between FORTRAN and C/R passed to FORTRAN */
      deriv_func = (C_deriv_func_type *) C_deriv_func_dop; 
      /* needed to communicate with R */
      R_deriv_func = derivfunc;
      R_envir = rho;
  }
  solout = C_solout;              
  iout = 2;                           /* solout called after each step */
  idid = 0;

/*                   ####      integration     ####                           */    
  it   = 0;
  tin  = REAL(times)[0];
  tout = REAL(times)[nt-1];

  saveOut (tin, xytmp);               /* save initial condition */ 
//  it++;

  if (type == 1)
    F77_CALL(dop853) ( &n_eq, deriv_func, &tin, xytmp, &tout,  
		     Rtol, Atol, &itol, solout, &iout,
		     rwork, &lrw, iwork, &liw, out, ipar, &idid);
  else if (type == 2)
    F77_CALL(dopri6) ( &n_eq, deriv_func, &tin, xytmp, &tout,  
		     Rtol, Atol, &itol, solout, &iout,
		     rwork, &lrw, iwork, &liw, out, ipar, &idid, ww);
  else if (type == 3)
    F77_CALL(cashkarp) ( &n_eq, deriv_func, &tin, xytmp, &tout,  
		     Rtol, Atol, &itol, solout, &iout,
		     rwork, &lrw, iwork, &liw, out, ipar, &idid, ww);
  if (idid == -1)  
     warning("input is not consistent");
  else if (idid == -2)   
     warning("larger maxsteps needed");
  else if (idid == -3)   
     warning("step size becomes too small");
  else if (idid == -4)   
     warning("problem is probably stiff - interrupted");
  else if (idid == -5)
     warning("stopped based on error estimation");
  else if (idid == -6)
     warning("stopped based on conditioning estimation; the problem is stiff (sigma > 100)");
  else if (idid == -7)
     warning("stopped based on conditioning estimation; the stepsize is restricted only by stability reason");
  else if (idid == -8)
     warning("stopped based on conditioning estimation; the stepsize is restricted only by stability reason and kappa > 1");
  else if (idid == -9)
     warning("stopped based on conditioning estimation; kappa >   1e20");

/*                   ####  an error occurred   ####                           */    
  if (idid < 0 ){
//    returnearly (1);
    warning("Returning early. Results are accurate, as far as they go\n");
    PROTECT(YOUT2 = allocMatrix(REALSXP,ntot+1,it));              nprot++;
    for (k = 0; k < it; k++)
      for (j = 0; j < ntot+1; j++)
        REAL(YOUT2)[k*(ntot+1) + j] = REAL(YOUT)[k*(ntot+1) + j];
  }  
  
  saveOut (tin, xytmp);              /* save final condition */

/*                   ####   returning output   ####                           */    
  rwork[0] = rwork[6];
  rwork[1] = rwork[6];
  rwork[2] = tin ;
  // terminate(idid,5,16,5,0);   
  
  PROTECT(ISTATE = allocVector(INTSXP, 5));                nprot++;
  for (k = 0; k < 4; k++) INTEGER(ISTATE)[k+1] = iwork[k +16];
  INTEGER(ISTATE)[0] = idid;  
    
  PROTECT(RWORK = allocVector(REALSXP, 5));                nprot++;
  for (k = 0; k < 5; k++) REAL(RWORK)[k] = rwork[k];
  if (idid > 0) {
      setAttrib(YOUT, install("istate"), ISTATE);
      setAttrib(YOUT, install("rstate"), RWORK);
  }
  else  {
      setAttrib(YOUT2, install("istate"), ISTATE);
      setAttrib(YOUT2, install("rstate"), RWORK);
  }

  
/*                   ####     termination      ####                           */    
  UNPROTECT(nprot);
  if (idid > 0)
    return(YOUT);
  else
    return(YOUT2);
}

