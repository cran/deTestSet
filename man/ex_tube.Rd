\name{tube}
\alias{tube}
\title{Water Tube System, Mechanics problem, DAE of Index 2}
\description{
The tube problem describes the water flow through a tube system, taking into 
account turbulence and the roughness of the tube walls.

It is an index 2 system of 49 non-linear Differential-Algebraic Equations.

}
\usage{
tube (times = seq(0, 17.0*3600, by = 100), yini = NULL, dyini = NULL, 
      parms = list(), printmescd = TRUE, method = radau,   
      atol = 1e-6, rtol = 1e-6, maxsteps = 1e+05, ...) 
}

\arguments{
  \item{yini }{the initial (state) values for the DE system. If \code{y}
    has a name attribute, the names will be used to label the output
    matrix.
  }
  \item{dyini }{the initial derivatives of the state variables of the DE
    system.  
  }
  \item{times }{time sequence for which output is wanted; the first
    value of \code{times} must be the initial time.
  }
  \item{parms }{list of parameters that overrule the default parameter values
  }
  \item{method }{the solver to use; only \code{mebdfi} available for now
  }
  \item{maxsteps }{maximal number of steps per output interval taken by 
    the solver
  }
 \item{atol }{absolute error tolerance, either a scalar or a vector,
    one value for each y.
  }
  \item{rtol }{relative error tolerance, either a scalar or a vector,
    one value for each y,
  }
  \item{printmescd }{if TRUE  the mixed error significant digits computed using
    the  reference solution at time 1e13  are printed
    } 
  \item{... }{additional arguments passed to the solver .
  }
}
\value{
  A matrix of class \code{deSolve} with up to as many rows as elements in
  \code{times} and as many
  columns as elements in \code{yini}, plus an additional column (the first) 
  for the time value.

  There will be one row for each element in \code{times} unless the
  solver returns with an unrecoverable error.  If
  \code{yini} has a names attribute, it will be used to label the columns
  of the output value.
      
}
\author{Karline Soetaert <karline.soetaert@nioz.nl>

Francesca Mazzia
}
\examples{
out <- tube()
plot(out, lwd = 2, which = 1:9)
plot(out, which = "phi3.4", lwd = 2, xlim = c(10000, 60000),
 ylim = c(0.000145, 0.000185))

# compare with reference solution
max(abs(out[nrow(out),-1]- reference("tube")))    
}
\references{
\url{http://www.dm.uniba.it/~testset}


}
\details{

    parameter <- c(nu = 1.31e-6, g = 9.8, rho = 1.0e3, rcrit = 2.3e3,
            length= 1.0e3, k = 2.0e-4, d= 1.0e0, b = 2.0e2)
  
}
\keyword{utilities}

