\name{caraxis}
\alias{caraxis}
\title{The Car Axis Mechanical Problem, Index 3 DAE}
\description{
A rather simple multibody system, describing the behavior of a car
axis on a bumpy road.

It is a differential algebraic equation of index 3

}
\usage{
caraxis (times = seq(0, 3, by = 0.01), yini = NULL, dyini = NULL,
   method = "mebdfi", parms = list(), ...)  
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

Francesca Mazzia <mazzia@dm.uniba.it>
}
\examples{
out <- caraxis()
plot(out, lwd = 2, mfrow = c(3, 4))

# compare with reference solution
out[nrow(out),2:11]-reference("caraxis")
}
\references{
\url{http://www.dm.uniba.it/~testset}


}
\details{

  The default parameters are:
    eps = 1e-2, M = 10, L = 1, L0 = 0.5,
          r   = 0.1,  w = 10, g = 1
  
}
\keyword{utilities}

