\name{E5}
\alias{E5}
\title{E5 Model for Chemical Pyrolysis, ODE}
\description{

It is an ODE, 4 equations

}
\usage{
E5 (times = c(0, 10^(seq(-5, 13, by = 0.1))), yini = NULL, 
    parms = list(), printmescd = TRUE, 
    atol = 1.11e-24, rtol = 1e-6,  maxsteps = 1e5, ...)
}

\arguments{
  \item{yini }{the initial (state) values for the DE system. If \code{y}
    has a name attribute, the names will be used to label the output
    matrix.
  }
  \item{times }{time sequence for which output is wanted; the first
    value of \code{times} must be the initial time.
  }
  \item{parms }{list of parameters that overrule the default parameter values
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
  \item{maxsteps }{maximal number of steps per output interval taken by 
    the solver
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

Francesca Mazzia  <mazzia@dm.uniba.it>
}
\examples{
out <- E5()
plot(out, lwd = 2, log = "xy")
# compare with reference solution
out[nrow(out),-1] - reference("E5")

}
\references{
\url{http://archimede.dm.uniba.it/~testset}

}
\note{
This model is implemented in R
}

\details{

  The default parameters are:
    A = 7.89e-10, B = 1.1e7, C = 1.13e3, M = 1e6
}
\keyword{utilities}

