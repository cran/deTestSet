\name{nand}
\alias{nand}
\title{Nand Gate, Index 1 IDE}
\description{

It is an index 1 IDE, 14 equations

}
\usage{
nand  (times = 0:80, yini =NULL, dyini = NULL,
  parms = list(), method = "mebdfi", maxsteps = 1e5, ...) 
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
out <- nand(method = "daspk")
plot(out, lwd = 2, which = 1:9)

# compare with reference solution
max(abs(out[nrow(out),-1] - reference("nand")))    
}
\references{
\url{http://www.dm.uniba.it/~testset}


}
\note{
This model is implemented in FORTRAN
}

\details{

  The default parameters are:
       RGS = 4, RGD = 4, RBS = 10, RBD = 10,
       CGS = 0.6e-4, CGD = 0.6e-4, CBD = 2.4e-5, CBS = 2.4e-5,
       C9 = 0.5e-4, DELTA = 0.2e-1, CURIS = 1.e-14, VTH = 25.85,
       VDD = 5., VBB = -2.5 
}
\keyword{utilities}

