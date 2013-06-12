\name{pleiades}
\alias{pleiades}
\title{Motion of Inextensible Elastic Beam, ODE}
\description{
The pleiades problem is a problem from celestial mechanics, describing the
motion of seven stars in the plane of coordinates xi,
yi and masses mi = i (i = 1, ... 7).

It is a set of nonstiff ordinary differential equations of dimension 28.

}
\usage{
pleiades (times = seq(0, 3.0, by = 0.01), yini = NULL, 
          printmescd = TRUE, method = lsoda,
          atol = 1e-6, rtol = 1e-6, ...)
}

\arguments{
  \item{yini }{the initial (state) values for the DE system. If \code{y}
    has a name attribute, the names will be used to label the output
    matrix.
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
out <- pleiades()
par(mfrow = c(3,3))
for (i in 1:7) plot(out[,i+1], out[,i+8], type = "l", main = paste("body ",i),
     xlab = "x", ylab = "y")


plot(0, 0 , type = "n", main = "ALL",
     xlab = "x", ylab = "y", xlim = c(-3, 4), ylim = c(-4, 5))
for (i in 1:7) lines(out[,i+1], out[,i+8], col = i, lwd = 2)

# compare with reference solution
max(abs(out[nrow(out),-1]- reference("pleiades")))

}
\references{
\url{http://www.dm.uniba.it/~testset}

}
\details{
There are no parameters  
}
\keyword{utilities}

