\name{orego}
\alias{orego}
\title{The Oregonator Chemistry Model, ODE}
\description{
Chemical model implementing the Belousov-Zhabotinskii reaction.                                  

It is an ODE, 3 equations

}
\usage{
orego (times = 0:360, yini = NULL, parms = list(),  ...) 
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
out <- orego()
plot(out, lwd = 2, log = "y")

# compare with exact solution
out[nrow(out),-1] - reference("orego")

}
\references{
\url{http://www.dm.uniba.it/~testset}


}
\note{
This model is implemented in R
}

\details{

  The default parameters are:
     k1 = 77.27, k2 = 8.375e-6, k3 = 77.27, k4 = 0.161
}
\keyword{utilities}

