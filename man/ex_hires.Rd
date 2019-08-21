\name{hires}
\alias{hires}
\title{High Irradiance Response model, from Plant Physiology, ODE}
\description{

This IVP is a stiff system of 8 non-linear Ordinary Differential Equations. 

It was proposed by Schafer in 1975 [Sch75]. 

The name HIRES was given by Hairer & Wanner [HW96]. It refers to 'High Irradiance
RESponse', which is described by this ODE. 

The parallel-IVP-algorithm group of CWI contributed this problem to the test 
set. The software part of the problem is in the file hires.f available at [MM08].

}
\usage{
hires (yini = NULL, times = seq(0, 321.8122, by = 321.8122/500), 
       parms = list(), printmescd = TRUE, method = mebdfi, 
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
   \item{atol }{absolute error tolerance, either a scalar or a vector,
    one value for each y.
  }
  \item{rtol }{relative error tolerance, either a scalar or a vector,
    one value for each y,
  }
  \item{method }{the solver to use
  }
  \item{parms }{list of parameters that overrule the default parameter values
  }
  \item{printmescd }{if TRUE  the mixed error significant digits computed using
    the  reference solution at time 321.8122 are printed}
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
out <- hires()
plot(out, lwd = 2)

# compare with reference solution
out1 <- hires(times = c(0, 321.8122))

max(abs(out1[nrow(out1),-1] - reference("hires")))
}
\references{
\url{http://archimede.dm.uniba.it/~testset}

[Got77] B.A. Gottwald. MISS - ein einfaches Simulations-System fur biologische und chemische
Prozesse. EDV in Medizin und Biologie, 3:85-90, 1977.

[HW96] E. Hairer and G. Wanner. Solving Ordinary Differential Equations II: 
Stiff and Differential-
algebraic Problems. Springer-Verlag, second revised edition, 1996.

[MM08] F. Mazzia and C. Magherini. Test Set for Initial Value Problem Solvers, 
release 2.4. Department of Mathematics, University of Bari and INdAM, 
Research Unit of Bari, February
2008. Available at http://www.dm.uniba.it/testset.

[Sch75] E. Schafer. A new approach to explain the 'high irradiance responses' 
of photomorphogenesis on the basis of phytochrome. 
J. of Math. Biology, 2:41 - 56, 1975.

[SL98] J.J.B. de Swart and W.M. Lioen. Collecting real-life problems to test 
solvers for implicit differential equations. CWI Quarterly, 11(1):83 - 100, 1998.
}
\note{
This model is implemented in R
}

\details{

  The default parameters are:
     k1 = 1.71, k2 = 0.43, k3 = 8.32, k4 = 0.69, k5 = 0.035,
       k6 = 8.32, k7 = 280, k8 = 0.69, k9 = 0.69, Oks = 0.0007
}
\keyword{utilities}

