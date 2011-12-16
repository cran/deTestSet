\name{dae}
\alias{dae}

\title{General Solver for Differential Algebraic Equations}

\description{Solves a system of differential algebraic equations; a wrapper
around the implemented DAE solvers}

\usage{dae(y, times, parms, dy, res = NULL, func = NULL,
    method = c("mebdfi", "daspk", "radau", "gamd", "bimd"), ...)
}
\arguments{
  \item{y }{the initial (state) values for the DAE system, a vector. If
    \code{y} has a name attribute, the names will be used to label the
    output matrix.
  }
  \item{times }{time sequence for which output is wanted; the first
    value of \code{times} must be the initial time.
  }
  \item{parms }{vector or list of parameters used in \code{res}
  }
  \item{dy }{the initial derivatives of the state variables of the DAE
    system.  
  }
  \item{func }{to be used if the model is an ODE, or a DAE written in linearly
    implicit form (M y' = f(t, y)). 
    \code{func} should be an \R-function that computes the
    values of the derivatives in the ODE system (the \emph{model
    definition}) at time t.
  
    \code{func} must be defined as:  \code{func <- function(t, y, parms,...)}.
    \cr \code{t} is the current time point in the
    integration, \code{y} is the current estimate of the variables in
    the ODE system.  If the initial values \code{y} has a \code{names}
    attribute, the names will be available inside \code{func}, unless
    \code{ynames} is FALSE.  \code{parms} is a vector or list of
    parameters.  \code{...} (optional) are any other arguments passed to
    the function.
    
    The return value of \code{func} should be a list,
    whose first element is a vector containing the derivatives of
    \code{y} with respect to \code{time}, and whose next elements are
    global values that are required at each point in \code{times}.
    The derivatives should be specified in the same order as the specification
    of the state variables \code{y}.

    } 
  \item{res }{either an \R-function that computes the
    residual function F(t,y,y') of the DAE system (the model
    defininition) at time \code{t}, or a character string giving the
    name of a compiled function in a dynamically loaded shared library.
      
    If \code{res} is a user-supplied R-function, it must be defined as:
    \code{res <- function(t, y, dy, parms, ...)}.
    
    Here \code{t} is the current time point in the integration, \code{y}
    is the current estimate of the variables in the DAE system,
    \code{dy} are the corresponding derivatives.  If the initial
    \code{y} or \code{dy} have a \code{names} attribute, the names will be
    available inside \code{res}, unless \code{ynames} is \code{FALSE}.
    \code{parms} is a vector of parameters.
    
    The return value of \code{res} should be a list, whose first element
    is a vector containing the residuals of the DAE system,
    i.e. delta = F(t,y,y'), and whose next elements contain output
    variables that are required at each point in \code{times}.
  
    If \code{res} is a string, then \code{dllname} must give the name of
    the shared library (without extension) which must be loaded before
    \code{dae()} is called (see deSolve package vignette \code{"compiledCode"}
    for more information).
  }
  \item{method }{the solver to use, either a string (\code{"mebdfi"},
    \code{"daspk"}), \code{"radau"}, \code{"gamd"} or a function that performs integration.
  }
  \item{... }{additional arguments passed to the solvers.
  }
}

\value{
  A matrix of class \code{deSolve} with up to as many rows as elements in
  \code{times} and as many
  columns as elements in \code{y} plus the number of "global" values
  returned in the second element of the return from \code{res}, plus an
  additional column (the first) for the time value.  There will be one
  row for each element in \code{times} unless the integrator returns
  with an unrecoverable error.  If \code{y} has a names attribute, it
  will be used to label the columns of the output value.

}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}

\details{
  This is simply a wrapper around the various dae solvers.
  
  See package vignette for information about specifying the model in
  compiled code.
  
  See the selected integrator for the additional options.
}

\seealso{
  \itemize{
    \item \code{\link[deSolve]{ode}} for a wrapper around the ode solvers,
    \item \code{\link[deSolve]{ode.band}} for solving models with a banded
      Jacobian,
    \item \code{\link[deSolve]{ode.1D}} for integrating 1-D models,
    \item \code{\link[deSolve]{ode.2D}} for integrating 2-D models,
    \item \code{\link[deSolve]{ode.3D}} for integrating 3-D models,
    \item \code{\link{mebdfi}}, \code{\link[deSolve]{daspk}},\code{\link[deSolve]{radau}},\code{\link{gamd}},\code{\link{bimd}}, for the  dae solvers
  }
  \code{\link[deSolve]{diagnostics}} to print diagnostic messages.
}
\keyword{math}
\examples{

## =======================================================================
## Chemical problem
## =======================================================================

daefun <- function(t, y, dy, parms) {
  with (as.list(c(y, dy, parms)), {
    res1 <- dA + dAB + lambda * A
    res2 <- dAB + dB
    alg <- B * A - K * AB
  list(c(res1, res2, alg), sumA = A + AB)
  })
}

parms <- c(lambda = 0.1, K = 1e-4)

yini  <- with(as.list(parms),
     c(A = 1, AB = 1, B = K))
dyini <- c(dA = 0, dAB = 0, dB = 0)

times <- 0:100

print(system.time(
out <-dae (y=yini, dy = dyini, times = times, res = daefun,
           parms = parms, method = "daspk")
))

plot(out, ylab = "conc.", xlab = "time", lwd = 2)
mtext("IVP DAE", side = 3, outer = TRUE, line = -1)

}

