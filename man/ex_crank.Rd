\name{crank}
\alias{crank}
\title{Slider Crank Mechanical Problem, Index 2 DAE}
\description{
The crank problem is a constrained mechanical system including both rigid
and elastic bodies

It is a differential algebraic equation of index 2, 24 equations.

}
\usage{
crank (times=seq(0, 0.1, by = 0.001), yini = NULL, dyini = NULL, 
  parms = list(), method = "mebdfi", maxsteps = 1e+06, 
  options = list(), ...) 
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
  \item{options }{a list which specifies the initial conditions used \code{ini},
    whether the problem is stiff, and the damping. The default is 
    \code{list(ini=1,stiff=0,damp=0)} 
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

Francesca Mazzia  <mazzia@dm.uniba.it>
}
\examples{
out <- crank()
plot(out, lwd = 2, which = 1:9)

# compare  with reference solution (only  the first seven component)
refsol <- reference("crank")
max(abs(out[nrow(out),2:8] - refsol[1:7])/refsol[1:7])     
}
\references{
\url{http://www.dm.uniba.it/~testset}

Simeon, B.: Modelling a flexible slider crank
          mechanism by a mixed system of DAEs and PDEs,
          Math. Modelling of Systems 2, 1-18 (1996);

}
\details{

  The default parameters are:
   M1 = 0.36, M2 = 0.151104, M3 = 0.075552, 
   L1 = 0.15, L2 = 0.30,  J1 = 0.002727, J2 = 0.0045339259,
   EE = 0.20e12, NUE = 0.30, BB = 0.0080,   HH = 0.0080,
   RHO = 7870.0,  GRAV = 0.0, OMEGA = 150.0 
  
  There are two default initial conditions - set with options(ini=x)
  
}
\keyword{utilities}

