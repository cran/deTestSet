\name{pollution}
\alias{pollution}
\title{Pollution Problem, from Chemistry, ODE}
\description{

This IVP is a stiff system of 20 non-linear Ordinary Differential Equations. 

It is the chemical reaction part of the air pollution model developed at 
The Dutch National Institute of Public Health and Environmental
Protection (RIVM) and it is described by Verwer in [Ver94]. 

The parallel-IVP-algorithm group of CWI contributed this problem to the test set.
The software part of the problem is in the file pollu.f available at [MM08].
}
\usage{
pollution (times = seq(0, 60, 1), yini = NULL, 
           parms = list(), printmescd = TRUE, method = mebdfi, 
           atol = 1e-6, rtol = 1e-6,  ...)   
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
out <- pollution()
plot(out, lwd = 2, which = 1:9)

out1 <- pollution(times = 0:60)

# compare with reference solution
max(abs(out1[nrow(out1),-1] -  reference("pollution")))    
}
\references{
\url{http://archimede.dm.uniba.it/~testset}

[MM08] F. Mazzia and C. Magherini. Test Set for Initial Value Problem Solvers, 
release 2.4. Department of Mathematics, University of Bari and INdAM, 
Research Unit of Bari, February 2008.

[Ver94] J.G. Verwer. Gauss-Seidel iteration for stiff ODEs from chemical 
kinetics. SIAM J. Sci.bComput., 15(5):1243 -- 1259,
}
\note{
This model is implemented in R
}

\details{

  The default parameters are:
      k1 = .35, k2 = .266e2, k3 = .123e5, k4 = .86e-3, k5 = .82e-3,
      k6 = .15e5, k7 = .13e-3, k8 = .24e5,k9 = .165e5,
      k10 = .9e4, k11 = .22e-1, k12 = .12e5, k13 = .188e1,
      k14 = .163e5, k15 = .48e7, k16 = .35e-3, k17 = .175e-1,
      k18 = .1e9, k19 = .444e12, k20 = .124e4, k21 = .21e1,
      k22 = .578e1, k23 = .474e-1, k24 = .178e4, k25 = .312e1
}
\keyword{utilities}

