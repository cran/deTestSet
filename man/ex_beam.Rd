\name{beam}
\alias{beam}
\title{Motion of Inextensible Elastic Beam, ODE}
\description{
The beam modulator problem is a problem from mechanics, describing the 
motion of an elastic beam, supposed inextensible, of length 1 and thin.

It is an ordinary differential equation of dimension 80.

}
\usage{
beam (times=seq(0, 5, by = 0.05), yini = NULL, 
      printmescd = TRUE, method = gamd, 
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
  \item{printmescd }{if TRUE  the mixed error significant digits computed using
    the  reference solution at time 5 are printed
    }
  \item{atol }{absolute error tolerance, either a scalar or a vector,
    one value for each y.
  }
  \item{rtol }{relative error tolerance, either a scalar or a vector,
    one value for each y,
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

Francesca Mazzia <mazzia@dm.uniba.it>
}
\examples{
out <- beam()
plot(out, col = "darkblue", lwd = 2, which = 1:16)
mtext(side = 3, line = -1.5, "beam", cex = 1.25, outer = TRUE)

image(out[,-1])

# compare with reference solution
refsol <- reference("beam")
max(abs(out[nrow(out),-1] - refsol)/refsol)  

}
\references{
\url{http://www.dm.uniba.it/~testset}

}
\details{
There are no parameters  
}
\keyword{utilities}

