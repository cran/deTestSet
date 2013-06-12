\name{dopri5}
\alias{dopri5}
\alias{cashkarp}

\title{Dormand-Prince or CashCarp Runge-Kutta of Order (4)5 }

\description{
  Solves the initial value problem for systems of
  ordinary differential equations (ODE) in the form: \deqn{dy/dt = 
  f(t,y)} 
  
  The \R function \code{dopri5} provides an interface to the Fortran ODE
  solver DOPRI5, written by E. Hairer and G. Wanner.
  
  It implements the explicit Runge-Kutta method of order 4(5)  
      due to Dormand & Prince with stepsize contral and dense output
  
  The \R function \code{cashkarp} provides an interface to the Fortran ODE
  solver CASHCARP, written by J. Cash and F. Mazzia.
  
  It implements the explicit Runge-Kutta method of order 4(5)  
      due to Cash-Carp, with stepsize contral and dense output
  
  The system of ODE's is written as an \R function or can be defined in
  compiled code that has been dynamically loaded.
  
}
\usage{
dopri5   (y, times, func, parms, rtol = 1e-6, atol = 1e-6,
  verbose = FALSE, hmax = NULL, hini = hmax, ynames = TRUE, 
  maxsteps = 10000, dllname = NULL, initfunc = dllname, 
  initpar=parms, rpar = NULL, ipar = NULL, nout = 0, 
  outnames = NULL, forcings = NULL, initforc = NULL, fcontrol = NULL, ...)

cashkarp (y, times, func, parms, rtol = 1e-6, atol = 1e-6,
  verbose = FALSE, hmax = NULL, hini = hmax, ynames = TRUE, 
  maxsteps = 10000, dllname = NULL, initfunc = dllname, initpar = parms, 
  rpar = NULL, ipar = NULL, nout = 0, outnames = NULL, forcings = NULL,
  initforc = NULL, fcontrol = NULL, stiffness = 2, ...)

}

\arguments{
  \item{y }{the initial (state) values for the ODE system. If \code{y}
    has a name attribute, the names will be used to label the output
    matrix.
  }
  \item{times }{time sequence for which output is wanted; the first
    value of \code{times} must be the initial time; if only one step is
    to be taken; set \code{times} = \code{NULL}.
  }
  \item{func }{either an \R-function that computes the values of the
    derivatives in the ODE system (the \emph{model definition}) at time
    t, or a character string giving the name of a compiled function in a
    dynamically loaded shared library.
  
    If \code{func} is an \R-function, it must be defined as:
     \code{func <- function(t, y, parms,...)}.  \code{t} is the current time
    point in the integration, \code{y} is the current estimate of the
    variables in the ODE system.  If the initial values \code{y} has a
    \code{names} attribute, the names will be available inside \code{func}.
    \code{parms} is a vector or list of parameters; ... (optional) are
    any other arguments passed to the function.
  
    The return value of \code{func} should be a list, whose first
    element is a vector containing the derivatives of \code{y} with
    respect to \code{time}, and whose next elements are global values
    that are required at each point in \code{times}.  The derivatives
    should be specified in the same order as the state variables \code{y}.

    If \code{func} is
    a string, then \code{dllname} must give the name of the shared
    library (without extension) which must be loaded before
    \code{lsode()} is called. See package vignette \code{"compiledCode"}
    for more details.
  }
  \item{parms }{vector or list of parameters used in \code{func} or
    \code{jacfunc}.
  }
  \item{rtol }{relative error tolerance, either a
    scalar or an array as long as \code{y}. See details.
  }
  \item{atol }{absolute error tolerance, either a scalar or an array as
    long as \code{y}. See details.
  }
  \item{verbose }{if \code{TRUE}: full output to the screen, e.g. will
    print the \code{diagnostiscs} of the integration - if the method becomes
    stiff it will rpint a message.
  }
  \item{hmax }{an optional maximum value of the integration stepsize. If
    not specified, \code{hmax} is set to the largest difference in
    \code{times}. 
  }
  \item{hini }{initial step size to be attempted.
  }
  \item{ynames }{logical, if \code{FALSE} names of state variables are not
    passed to function \code{func}; this may speed up the simulation especially
    for multi-D models.
  }
  \item{maxsteps }{maximal number of steps taken by the solver, 
    \emph{for the entire integration}. This is different from the settings
    of this argument in the solvers from package deSolve!
  }
  \item{dllname }{a string giving the name of the shared library
    (without extension) that contains all the compiled function or
    subroutine definitions refered to in \code{func} and
    \code{jacfunc}. See vignette \code{"compiledCode"}
    from package \code{deSolve}.
  }
  \item{initfunc }{if not \code{NULL}, the name of the initialisation function
    (which initialises values of parameters), as provided in
    \file{dllname}. See vignette \code{"compiledCode"}
    from package \code{deSolve}.
  }
  \item{initpar }{only when \file{dllname} is specified and an
    initialisation function \code{initfunc} is in the dll: the
    parameters passed to the initialiser, to initialise the common
    blocks (FORTRAN) or global variables (C, C++).
  }
  \item{rpar }{only when \file{dllname} is specified: a vector with
    double precision values passed to the dll-functions whose names are
    specified by \code{func} and \code{jacfunc}.
  }
  \item{ipar }{only when \file{dllname} is specified: a vector with
    integer values passed to the dll-functions whose names are specified
    by \code{func} and \code{jacfunc}.
  }
  \item{nout }{only used if \code{dllname} is specified and the model is
    defined in compiled code: the number of output variables calculated
    in the compiled function \code{func}, present in the shared
    library. Note: it is not automatically checked whether this is
    indeed the number of output variables calculed in the dll - you have
    to perform this check in the code - See vignette \code{"compiledCode"}
    from package \code{deSolve}.
  }
  \item{outnames }{only used if \file{dllname} is specified and
    \code{nout} > 0: the names of output variables calculated in the
    compiled function \code{func}, present in the shared library.
    These names will be used to label the output matrix.
  }
  \item{forcings }{only used if \file{dllname} is specified: a list with
    the forcing function data sets, each present as a two-columned matrix,
    with (time,value); interpolation outside the interval
    [min(\code{times}), max(\code{times})] is done by taking the value at
    the closest data extreme.

    See \link[deSolve]{forcings} or package vignette \code{"compiledCode"}.
  }
  \item{initforc }{if not \code{NULL}, the name of the forcing function
    initialisation function, as provided in
    \file{dllname}. It MUST be present if \code{forcings} has been given a
    value.
    See \link[deSolve]{forcings} or package vignette \code{"compiledCode"}.
  }
  \item{fcontrol }{A list of control parameters for the forcing functions.
    See \link[deSolve]{forcings} or vignette \code{compiledCode}.
  }
  \item{stiffness }{How the stiffness of the solution should be estimated.
   Default = stiffness based on eigenvalue approximation; 
   when = \code{stiffness = 0}: no stiffness estimate; when = \code{stiffness = 1} 
   or \code{-1}: all stiffness estimates calculated ; when = \code{stiffness = 2} or \code{-2}: stiffness based on
   eigenvalue approximation;  when = \code{stiffness = 3} or \code{-3}: stiffness
   based on error estimate; when = \code{stiffness = 4} or \code{-4}: stiffness
   based on conditioning. Positive values of \code{stiffness} will cause the
   integration to stop; negative values will continue anyway.
  }
  \item{... }{additional arguments passed to \code{func} and
    \code{jacfunc} allowing this to be a generic function.
  }
}
\value{
  A matrix of class \code{deSolve} with up to as many rows as elements
  in \code{times} and as many columns as elements in \code{y} plus the number of "global"
  values returned in the next elements of the return from \code{func},
  plus and additional column for the time value.  There will be a row
  for each element in \code{times} unless the FORTRAN routine `lsoda'
  returns with an unrecoverable error. If \code{y} has a names
  attribute, it will be used to label the columns of the output value.

}
\author{Karline Soetaert <karline.soetaert@nioz.nl> 
}
\examples{

## =======================================================================
## Example :
##   The Arenstorff orbit model
## =======================================================================

Arenstorff <- function(t, y, parms) {

  D1 <- ((y[1]+mu)^2+y[2]^2)^(3/2)
  D2 <- ((y[1]-(1-mu))^2+y[2]^2)^(3/2)

  dy1 <- y[3]
  dy2 <- y[4]
  dy3 <- y[1] + 2*y[4]-(1-mu)*(y[1]+mu)/D1 -mu*(y[1]-(1-mu))/D2
  dy4 <- y[2] - 2*y[3]-(1-mu)*y[2]/D1 - mu*y[2]/D2

  list(c(dy1,dy2,dy3,dy4))
}

#-----------------------------
# parameters, initial values and times
#-----------------------------
mu    <- 0.012277471
yini  <- c(x = 0.994, y = 0, dx = 0, 
  dy = -2.00158510637908252240537862224)
times <- seq(0, 18, 0.01)

#-----------------------------
# solve the model
#-----------------------------

#out <- dopri5 (times=times, y=yini, func = Arenstorff, parms=NULL )
out  <- cashkarp (times = times, y = yini, func = Arenstorff, parms = NULL )
plot(out[,c("x", "y")], type = "l", lwd = 2, main = "Arenstorff")

#-----------------------------
# First and last value should be the same
#-----------------------------

times <- c(0, 17.0652165601579625588917206249)

Test  <- dopri5 (times = times, y = yini, func = Arenstorff, parms = NULL)

diagnostics(Test)

}
\references{
E. Hairer, S.P. Norsett AND G. Wanner, Solving Ordinary
      Differential Equations I. Nonstiff Problems. 2nd Edition. 
      Springer Series In Computational Mathematics, 
       SPRINGER-VERLAG (1993)


}
\details{
  The work is done by the FORTRAN subroutine \code{dop853}, whose
  documentation should be consulted for details. The implementation
  is based on the Fortran 77 version fromOctober 11, 2009.
    
  The input parameters \code{rtol}, and \code{atol} determine the
  \bold{error control} performed by the solver, which roughly keeps the 
  local error of y(i) below rtol(i)*abs(y(i))+atol(i).  
   
  The diagnostics of the integration can be printed to screen
  by calling \code{\link[deSolve]{diagnostics}}. If \code{verbose} = \code{TRUE},
  the diagnostics will written to the screen at the end of the integration.

  See vignette("deSolve") from the \code{deSolve} package for an 
  explanation of each element in the vectors
  containing the diagnostic properties and how to directly access them.

  \bold{Models} may be defined in compiled C or FORTRAN code, as well as
  in an R-function. See package vignette \code{"compiledCode"} from package 
  \code{deSolve} for details.

  Information about linking forcing functions to compiled code is in 
  \link[deSolve]{forcings} (from package \code{deSolve}).

}
\seealso{
  \itemize{
    \item \code{\link[deSolve]{ode}} for a general interface to most of the ODE solvers
      from package \code{deSolve},
    \item \code{\link[deSolve]{ode.1D}} for integrating 1-D models,
    \item \code{\link[deSolve]{ode.2D}} for integrating 2-D models,
    \item \code{\link[deSolve]{ode.3D}} for integrating 3-D models,
    \item \code{\link{mebdfi}} for integrating DAE models,
    \item \code{\link{bimd}} for blended implicit methods,
    \item \code{\link{gamd}} for the generalised adams method
  }

  \code{\link[deSolve]{diagnostics}} to print diagnostic messages.
}
\keyword{math}