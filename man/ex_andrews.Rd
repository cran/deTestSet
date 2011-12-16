\name{andrews}
\alias{andrews}
\title{Andrews Squeezing Mechanism, Index 3 DAE}
\description{
The andrews problem describes the motion of 7 rigid bodies connected 
by joints without friction

It is a non-stiff second order differential algebraic equation of index
3, consisting of 14 differential and 13 algebraic equations

}
\usage{
andrews (times = seq(0, 0.03, by = 0.001), yini = NULL, 
         dyini = NULL, parms = list(), method = "mebdfi", 
         atol = 1e-7, rtol = 1e-7, maxsteps = 1e+05, ...) 
}

\arguments{
  \item{yini }{the initial (state) values for the DE system. If \code{y}
    has a name attribute, the names will be used to label the output
    matrix.
  }
  \item{dyini }{the initial derivatives of the state variables of the DE
    system.  
  }
  \item{atol }{absolute error tolerance, either a scalar or a vector,
    one value for each y.
  }
  \item{rtol }{relative error tolerance, either a scalar or a vector,
    one value for each y,
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
out <- andrews()
plot(out, lwd = 2, which = 1:9)

refsol <- reference("andrews")
max(abs(out[nrow(out),-1] - refsol)/refsol)    
}
\references{
\url{http://www.dm.uniba.it/~testset}
}
\details{

The default parameters are:
parameter <- c(m1 = .04325, m2 = .00365, m3 = .02373, m4 = .00706 ,
            m5 = .07050, m6 = .00706, m7 = .05498 ,
            xa = -.06934, ya = -.00227 ,
            xb = -0.03635, yb = .03273 ,
            xc = .014 , yc = .072, c0 = 4530 ,
            i1 = 2.194e-6, i2 = 4.410e-7, i3 = 5.255e-6, i4 = 5.667e-7,
            i5 = 1.169e-5, i6 = 5.667e-7, i7 = 1.912e-5,
            d = 28e-3, da = 115e-4,e=2e-2, ea = 1421e-5,
            rr = 7e-3, ra = 92e-5, l0 = 7785e-5,
            ss = 35e-3, sa = 1874e-5, sb = 1043e-5, sc = 18e-3, sd = 2e-2,
            ta = 2308e-5, tb = 916e-5, u = 4e-2, ua = 1228e-5, ub = 449e-5,
            zf = 2e-2, zt = 4e-2,fa=1421e-5, mom = 33e-3)
  
  
}
\keyword{utilities}

