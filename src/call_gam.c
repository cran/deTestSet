#include <time.h>
#include <string.h>
#include "de.h"
#include "externalptr.h"

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   differential algebraic equation solvers gamd and bimd.

   The C-wrappers that provide the interface between FORTRAN codes and R-code
   are: C_deriv_func_gb: interface with R-code "func", passes derivatives
        C_deriv_out_gb : interface with R-code "func", passes derivatives + output variables
        C_jac_func_gb  : interface with R-code "jacfunc", passes jacobian

   C_deriv_func_forc_gb provides the interface between the function specified in
   a DLL and the integrator, in case there are forcing functions.

   karline soetaert
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int maxt;
  C_deriv_func_type     *deriv_func;  

/* definition of the calls to the FORTRAN subroutines in file gamd.f**/

void F77_NAME(gamd)( int *,
         void (*)(int *, double *, double *, double *,
                              double *, int *),         // func
		     double *, double *, double *, double *,
		     double *, double *, int *,
		     void (*)(int *, double *, double *, int *, int*, double *,  // jac
			            int *, double*, int*),
		     int *, int *, int *,
 	       void (*)(int *, double *, int *, double *, int *),              // mas
		     int *, int *, int *,
 	       void (*)(int *, double *, double *, double *,  // solout
			            int *, int *, int *, double *, int *, int* ),
		     int *, double *, int *, int *, int*, double *, int*, int*);

void F77_NAME(contout)(int *, double *, double *, double *, int *, int *, double *);

/* definition of the calls to the FORTRAN subroutines in file bim.f**/

void F77_NAME(bimd)( int *,
         void (*)(int *, double *, double *, double *,
                              double *, int *),         // func
		     double *, double *, double *, double *,
		     double *, double *, int *,
		     void (*)(int *, double *, double *, int *, int*, double *,  // jac
			            int *, double*, int*),
		     int *, int *, int *,
 	       void (*)(int *, double *, int *, double *, int *),  // mas
		     int *, int *, int *,
 	       void (*)(int *, int *, int *, double *,  double *,  // soloutbim
			            double *, double *, double *, double *, int *, int *),
		     int *, double *, int *, int *, int*, double *, int*, int*);

void F77_NAME(contsolall)( double *, int *, int *, double *, double *, double *, double *);


/* Each succesful timestep, the solver enters here.
   Check if needs interpolating to output time */



/* wrapper above the derivate function in a dll that first estimates the
values of the forcing functions */

static void C_deriv_func_forc_gb (int *neq, double *t, double *y,
                         double *ydot, double *yout, int *iout)
{
  updatedeforc(t);
  DLL_deriv_func(neq, t, y, ydot, yout, iout);
}

/* interface between FORTRAN function call and R function
   Fortran code calls C_deriv_func_gb(N, t, y, ydot, yout, iout)
   R code called as R_deriv_func(time, y) and returns ydot */

static void C_deriv_func_gb (int *neq, double *t, double *y,
                          double *ydot, double *yout, int *iout)
{
  int i;
  SEXP R_fcall, Time, ans;

  for (i = 0; i < *neq; i++)  REAL(Y)[i] = y[i];

  PROTECT(Time = ScalarReal(*t));                  incr_N_Protect();
  PROTECT(R_fcall = lang3(R_deriv_func,Time,Y));   incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));           incr_N_Protect();

  for (i = 0; i < *neq; i++)   ydot[i] = REAL(ans)[i];

  my_unprotect(3);
}

/* deriv output function - for ordinary output variables */

static void C_deriv_out_gb (int *nOut, double *t, double *y,
                       double *ydot, double *yout)
{
  int i;
  SEXP R_fcall, Time, ans;

  for (i = 0; i < n_eq; i++)
      REAL(Y)[i] = y[i];

  PROTECT(Time = ScalarReal(*t));                   incr_N_Protect();
  PROTECT(R_fcall = lang3(R_deriv_func,Time, Y));   incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));            incr_N_Protect();

  for (i = 0; i < *nOut; i++) yout[i] = REAL(ans)[i + n_eq];

  my_unprotect(3);
}
/* the mass matrix function */
static void C_mas_func (int *neq, double *am, int *lmas,
                             double *yout, int *iout)
{
  int i;
  SEXP NEQ, LM, R_fcall, ans;

  PROTECT(NEQ = NEW_INTEGER(1));                  incr_N_Protect();
  PROTECT(LM = NEW_INTEGER(1));                   incr_N_Protect();

                              INTEGER(NEQ)[0] = *neq;
                              INTEGER(LM) [0] = *lmas;
  PROTECT(R_fcall = lang3(R_mas_func,NEQ,LM));   incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));         incr_N_Protect();

  for (i = 0; i <*lmas * *neq; i++)   am[i] = REAL(ans)[i];

  my_unprotect(4);
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
        C_deriv_out_gb(&nout, &t, y, xdytmp, out);
      for (j = 0; j < nout; j++)
        REAL(YOUT)[(it)*(ntot + 1) + j + n_eq + 1] = out[j];
    }
}

/* function called by Fortran to check for output */

static void C_solout_gam (int * neq, double * tp, double * yp, double * ff,
  int *nt, int * dblk, int * ord, double * rpar, int * ipar, int * irtrn)
{
  *irtrn = 1;
  while ((tp[0] <= tt[it]) && (tt[it] < tp[*dblk])) {
 	  F77_CALL(contout) (neq, &tt[it], tp, ff, dblk, nt, ytmp);
    saveOut(tt[it], ytmp);
	  it++;
	  if (it >= maxt) break;
  }
}

/* function called by Fortran to check for output */
static void C_solout_bim (int * m, int *k, int * ord,
   double * t0, double * tstep, double * y, double * f,
   double *dd, double * rpar, int * ipar, int * irtrn)
{
  *irtrn = 1;
  while ((*t0 <= tt[it]) && (tt[it] < tstep[*k-1])) {
 	  F77_CALL(contsolall) (&tt[it], m, k, t0, tstep, dd, ytmp);
    saveOut(tt[it], ytmp);
	  it++;
	  if (it >= maxt) break;

  }
}

/* interface to jacobian function */

static void C_jac_func_gb (int *neq, double *t, double *y, int *ml,
	    int *mu, double *pd,  int *nrowpd, double *yout, int *iout)
{
  int i;
  SEXP R_fcall, Time, ans;

  for (i = 0; i < *neq; i++) REAL(Y)[i] = y[i];

  PROTECT(Time = ScalarReal(*t));                 incr_N_Protect();
  PROTECT(R_fcall = lang3(R_jac_func,Time,Y));    incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));          incr_N_Protect();

  for (i = 0; i < *neq * *nrowpd; i++)  pd[i] = REAL(ans)[i];

  my_unprotect(3);
}

/* give name to data types */
typedef void C_jac_func_type_gb  (int *, double *, double *, int *, int *, double *,
		                    int *, double *, int *);

typedef void C_solout_type (int *, double *, double *, double *,
  int *, int *, int *, double *, int *, int *) ;

typedef void C_solout_type_bim (int *, int *, int *, double *, double *,
                double *, double *, double *, double *, int *, int *) ;

typedef void C_mas_type (int *, double *, int *, double *, int *);

/* MAIN C-FUNCTION, CALLED FROM R-code */

SEXP call_gambim(SEXP y, SEXP times, SEXP derivfunc, SEXP parms, SEXP rtol,
		SEXP atol, SEXP rho, SEXP Tcrit, SEXP jacfunc, SEXP initfunc,
    SEXP verbose, SEXP LRW, SEXP rWork, SEXP iWork, SEXP jT,
    SEXP nOut, SEXP Nrmas, SEXP masfunc, SEXP ML, SEXP MU, SEXP Hini,
    SEXP Rpar, SEXP Ipar, SEXP flist, SEXP Type)

{
/******************************************************************************/
/******                   DECLARATION SECTION                            ******/
/******************************************************************************/

  int  j, nt, latol, lrtol, imas, mlmas, mumas, type;
  int  isForcing, runOK;
  double *Atol, *Rtol, hini;
  int itol, ijac, ml, mu, iout, idid, liw, lrw, sum;

  /* pointers to functions passed to FORTRAN */
  C_jac_func_type_gb    *jac_func_gb = NULL;
  C_solout_type         *solout = NULL;
  C_solout_type_bim     *solout_bim = NULL;
  C_mas_type            *mas_func = NULL;

/******************************************************************************/
/******                         STATEMENTS                               ******/
/******************************************************************************/

/*                      #### initialisation ####                              */
  init_N_Protect();

  type  = INTEGER(Type)[0];     /* jacobian type */
  ijac  = INTEGER(jT)[0];       /* jacobian type */
  n_eq = LENGTH(y);             /* number of equations */
  nt   = LENGTH(times);         /* number of output times */
  maxt = nt;

  tt = (double *) R_alloc(nt, sizeof(double));
  for (j = 0; j < nt; j++) tt[j] = REAL(times)[j];

//  mflag = INTEGER(verbose)[0];

  /* is function a dll ?*/
  isDll = inherits(derivfunc, "NativeSymbol");

  /* initialise output ... */
  initOutC(isDll, n_eq, nOut, Rpar, Ipar);

  /* copies of variables that will be changed in the FORTRAN subroutine */
  xytmp = (double *) R_alloc(n_eq, sizeof(double));
  for (j = 0; j < n_eq; j++) xytmp[j] = REAL(y)[j];

  ytmp = (double *) R_alloc(n_eq, sizeof(double));
  for (j = 0; j < n_eq; j++) ytmp[j] = 0.;

  latol = LENGTH(atol);
  Atol = (double *) R_alloc((int) latol, sizeof(double));
  for (j = 0; j < latol; j++) Atol[j] = REAL(atol)[j];

  lrtol = LENGTH(rtol);
  Rtol = (double *) R_alloc((int) lrtol, sizeof(double));
  for (j = 0; j < lrtol; j++) Rtol[j] = REAL(rtol)[j];

  /* tolerance specifications */
  if (latol == 1 ) itol = 0;
  else             itol = 1;

  /* mass matrix */
  imas  = INTEGER(Nrmas)[0];
  mlmas = INTEGER(Nrmas)[1];
  mumas = INTEGER(Nrmas)[2];

  /* work vectors */
  if (type == 1) {
    liw = 27;
    lrw = 21;
  } else  {       //  if (type == 2)
    liw = n_eq + 40;
    lrw = INTEGER(LRW)[0];
 }

  iwork = (int *) R_alloc(liw, sizeof(int));
  for (j=0; j<LENGTH(iWork); j++) iwork[j] = INTEGER(iWork)[j];
  for (j=LENGTH(iWork); j<liw; j++) iwork[j] = 0;

  rwork = (double *) R_alloc(lrw, sizeof(double));
  for (j=0; j<length(rWork); j++) rwork[j] = REAL(rWork)[j];
  for (j=length(rWork); j<lrw; j++) rwork[j] = 0.;

  ml = INTEGER(ML)[0];
  mu = INTEGER(MU)[0];
  hini = REAL(Hini)[0];

  /* initialise global R-variables...  */
  initglobals (nt);

  /* Initialization of Parameters, Forcings (DLL) */
  initParms(initfunc, parms);
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
        deriv_func = (C_deriv_func_type *) C_deriv_func_forc_gb;
      }
  } else {
      /* interface function between FORTRAN and C/R passed to FORTRAN */
      deriv_func = (C_deriv_func_type *) C_deriv_func_gb;
      /* needed to communicate with R */
      R_deriv_func = derivfunc;
      R_envir = rho;
  }

  if (!isNull(jacfunc))   {
      if (isDll)
	      jac_func_gb = (C_jac_func_type_gb *) R_ExternalPtrAddrFn_(jacfunc);
	    else  {
	      R_jac_func = jacfunc;
	      jac_func_gb= C_jac_func_gb;
	    }
    }

  if (!isNull(masfunc))   {
	      R_mas_func = masfunc;
	      mas_func= C_mas_func;
  }

  solout = C_solout_gam;
  solout_bim = C_solout_bim;
  iout = 1;                           /* solout used */
  idid = 0;

/*                   ####      integration     ####                           */
  it   = 0;
  tin  = REAL(times)[0];
  tout = REAL(times)[nt-1];

  saveOut (tin, xytmp);               /* save initial condition */
  it = it +1;

  if (type == 1)
    F77_CALL(gamd) ( &n_eq, deriv_func, &tin, xytmp, &tout, &hini,
		     Rtol, Atol, &itol, jac_func_gb, &ijac, &ml, &mu,
         mas_func, &imas, &mlmas, &mumas, solout, &iout,
		     rwork, &lrw, iwork, &liw, out, ipar, &idid);
  else if (type == 2)
    F77_CALL(bimd) ( &n_eq, deriv_func, &tin, &tout, xytmp, &hini,
		     Rtol, Atol, &itol, jac_func_gb, &ijac, &ml, &mu,
         mas_func, &imas, &mlmas, &mumas, solout_bim, &iout,
		     rwork, &lrw, iwork, &liw, out, ipar, &idid);


  if (idid == -1)
     warning("input is not consistent");
  else if (idid == -2)
     warning("larger maxsteps needed");
  else if (idid == -3)
     warning("step size becomes too small");
  else if (idid == -4)
     warning("matrix is repeatedly singular");

/*                   ####  an error occurred   ####                           */
  if(it <= nt-1) saveOut (tin, xytmp);              /* save final condition */
  if (idid < 0 ) {
    it = it-1;
    returnearly (1);
  }

/*                   ####   returning output   ####                           */

/* feval */

  PROTECT(RWORK = allocVector(REALSXP, 3)); incr_N_Protect();
  REAL(RWORK)[0] = hini;
  REAL(RWORK)[1] = hini;

  REAL(RWORK)[2] = tin;

  PROTECT(ISTATE = allocVector(INTSXP, 6)); incr_N_Protect();
  INTEGER(ISTATE)[0] = idid;

/* nsteps */
  sum = 0;
  runOK = 0;
  if (type == 1)  {
    for (j = 11; j < 23; j++) sum = sum + iwork[j];
    INTEGER(ISTATE)[1] = sum;

/* feval */
    INTEGER(ISTATE)[2] = iwork[9];

/* jacobian eval */
    INTEGER(ISTATE)[3] = iwork[10];

/* LU decomp */
    INTEGER(ISTATE)[4] = iwork[23];

/* number rejected steps */
    sum = 0;
    for (j = 11; j < 15; j++) sum = sum + iwork[j];
    INTEGER(ISTATE)[5] = INTEGER(ISTATE)[1]- sum;
	if(idid > 0) runOK = 1;
  } else {
    for (j = 20; j < 25; j++) sum = sum + iwork[j];
    INTEGER(ISTATE)[1] = sum;

/* feval */
    INTEGER(ISTATE)[2] = iwork[11];

/* jacobian eval */
    INTEGER(ISTATE)[3] = iwork[12];

/* LU decomp */
    INTEGER(ISTATE)[4] = iwork[13];

/* number rejected steps */
    sum = 0;
    for (j = 25; j < 29; j++) sum = sum + iwork[j];
    INTEGER(ISTATE)[5] = INTEGER(ISTATE)[1]- sum;
    if(idid >= 0)  runOK = 1;

}


  if (runOK) {
    setAttrib(YOUT, install("istate"), ISTATE);
    setAttrib(YOUT, install("rstate"), RWORK);
  }
  else  {
    setAttrib(YOUT2, install("istate"), ISTATE);
    setAttrib(YOUT2, install("rstate"), RWORK);
  }

/*                   ####     termination      ####                           */
  unprotect_all();
  if (runOK)
    return(YOUT);
  else
    return(YOUT2);
}





