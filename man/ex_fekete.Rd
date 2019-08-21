\name{fekete}
\alias{fekete}
\title{Elliptic Fekete Points, Mechanical Problem, Index 2 DAE}
\description{
The fekete problem computes the elliptic Fekete points.  

 
}
\usage{
fekete (times = seq(0, 1e3, by = 10 ), yini = NULL, dyini = NULL, 
       parms=list(), printmescd = TRUE, method = mebdfi, 
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
  \item{method }{the solver to use
  }
   \item{rtol }{relative error tolerance, either a scalar or a vector,
    one value for each y,
  }
  \item{atol }{absolute error tolerance, either a scalar or a vector,
    one value for each y.
  }
   \item{printmescd }{if TRUE  the mixed error significant digits computed using
    the  reference solution at time 0.1  are printed
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

Francesca Mazzia <mazzia@dm.uniba.it>
}
\examples{
out <- fekete()
plot(out, lwd = 2, which = 1:20)

# reference run compared with output at end of interval for first 7 components
out1 <- fekete(times = c(0, 1000))
max(abs(out1[nrow(out1),-1] - reference("fekete")))    
}
\references{
\url{http://archimede.dm.uniba.it/~testset}
}
\details{

There are no parameters
}
\keyword{utilities}

