\name{twobit}
\alias{twobit}
\title{The Two Bit Adding Unit, Index 1 DAE}
\description{
Computes the sum of two base-2 numbers, each two digits long, and a carry bit.
These numbers are fed into the circuit in the form of input signals.

Index 1 DAE of dimension 350

}
\usage{
twobit (times = seq(0, 320, by = 0.5), yini = NULL, dyini = NULL, 
  method = "radau",  atol = 1e-4, rtol = 1e-4, 
  maxsteps = 1e5, hmax = 0.1, ...) 
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
  \item{method }{the solver to use
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
  \item{hmax }{maximal size of step; if too large: will fail.
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
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>

Francesca Mazzia
}
\examples{

out <- twobit()
plot(out, lwd = 2, which = c("x49", "x130", "x148"), mfrow = c(3, 1))

# compare with reference solution
max(abs(out[nrow(out),-1] - reference("twobit")))    
}
\references{
\url{http://www.dm.uniba.it/~testset}


}
\note{
This model is implemented in FORTRAN
}

\details{

  This model has no parameters
}
\keyword{utilities}

