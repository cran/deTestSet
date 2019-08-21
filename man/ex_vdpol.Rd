\name{vdpol}
\alias{vdpol}
\title{van der Pol Equation, Nonlinear Vacuum Tube Circuit, ODE}
\description{

Problem originating from electronics, describing the behavior of 
nonlinear vacuum tube circuots.

It is an ODE, 2 equations. 
}
\usage{
vdpol  (times = 0:2000, yini = NULL,  
        parms = list(), printmescd = TRUE,
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
  \item{parms }{list of parameters that overrule the default parameter values
  }
  \item{atol }{absolute error tolerance, either a scalar or a vector,
    one value for each y.
  }
  \item{rtol }{relative error tolerance, either a scalar or a vector,
    one value for each y,
  }
  \item{printmescd }{if TRUE  the mixed error significant digits computed using
    the  reference solution at time 5 are printed
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
out <- vdpol()
plot(out, lwd = 2, which = 1)

# compare to reference solution
out[nrow(out),-1] - reference("vdpol")

}
\references{
\url{http://archimede.dm.uniba.it/~testset}

}
\note{
This model is implemented in R
}

\details{

  The default parameters are:
       mu=1000
       
  The default initial conditions are: y1 = 2, y2 = 0     
}
\keyword{utilities}

