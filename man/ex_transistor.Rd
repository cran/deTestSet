\name{transistor}
\alias{transistor}
\title{The Transistor Amplifier, Index 1 DAE}
\description{

Electrical circuit model for the transistor amplifier. 

It is an index 1 DAE, 8 equations

}
\usage{
transistor (times = seq(0, 0.2, 0.001), yini = NULL, dyini = NULL,
            parms = list(), printmescd = TRUE, method = mebdfi,
            atol = 1e-6, rtol = 1e-6, maxsteps = 1e5, ...) 
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
  \item{method }{the solver to use
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
out <- transistor()
plot(out, lwd = 2)

out[nrow(out),-1]-reference("transistor")
}
\references{
\url{http://archimede.dm.uniba.it/~testset}

}
\note{
This model is implemented in FORTRAN
}

\details{

  The default parameters are:
       ub=6, uf=0.026, alpha=0.99, beta=1e-6,
                r0=1000, r1=9000, r2=9000, r3=9000, r4 = 9000,
                r5=9000, r6=9000, r7=9000, r8 = 9000, r9 = 9000,
                c1=1e-6, c2=2e-6, c3=3e-6, c4=4e-6, c5=5e-6
}
\keyword{utilities}

