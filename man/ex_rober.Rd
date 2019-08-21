\name{rober}
\alias{rober}
\title{Autocatalytic Chemical Reaction of Robertson, ODE}
\description{
Describes the kinetics of an autocatalytic reation.

It is an ODE, 3 equations

}
\usage{
rober  (times = 10^(seq(-5, 11, by = 0.1)), yini = NULL, 
        parms = list(), printmescd = TRUE,   
        atol = 1e-14, rtol = 1e-10, maxsteps = 1e5, ...)  
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
  \item{maxsteps }{maximal number of steps per output interval taken by 
    the solver
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
out <- rober()
plot(out, lwd = 2, log = "x")

# compare to reference solution
out[nrow(out),-1] - reference("rober")

}
\references{
\url{http://archimede.dm.uniba.it/~testset}

}
\note{
This model is implemented in R
}

\details{

  The default parameters are:
       k1 = 0.04, k2 = 3e7, k3 = 1e4
}
\keyword{utilities}

