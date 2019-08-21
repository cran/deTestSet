\name{ring}
\alias{ring}
\title{The Ring Modulator Problem, Electrical Circuit Analysis, ODE}
\description{
The probelm describes the behavior of the ring modulator, an electrical
circuit.

The type of the problem depends on the parameter Cs. If Cs is not equal 0, then it is a
 stiff system of 15
non-linear ordinary differential equations. 

For Cs = 0 we have a DAE of index 2, consisting of 11 differential equations 
and 4 algebraic equations. The numerical results presented here refer to
Cs = 2 e-12. The problem has been taken from [KRS92], where the approach of 
Horneber [Hor76] is followed. 
The parallel-IVP-algorithm group of CWI contributed this problem to the test set.
The software part of the problem is in the file ringmod.f available at [MM08].
}
\usage{
ring (times = seq(0, 0.001, by = 5e-06), yini = NULL, dyini = NULL, 
      parms = list(), printmescd = TRUE, method = mebdfi,
      atol = 1e-8, rtol = 1e-8, maxsteps = 1e+06, ...) 
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
  \item{atol }{absolute error tolerance, either a scalar or a vector,
    one value for each y.
  }
  \item{rtol }{relative error tolerance, either a scalar or a vector,
    one value for each y,
  }
  \item{maxsteps }{maximal number of steps per output interval taken by 
    the solver
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

Francesca Mazzia <mazzia@dm.uniba.it>
}
\examples{
out <- ring()
plot(out, col = "darkblue", lwd = 2)
mtext(side = 3, line = -1.5, "RING modulator",cex = 1.25, outer = TRUE)

# compare with reference solution
max(abs(out[nrow(out),-1]- reference("ring")))    
}
\references{
\url{http://archimede.dm.uniba.it/~testset}

[Hor76] E.H. Horneber. Analyse nichtlinearer RLCU-Netzwerke mit Hilfe der 
gemischten Potentialfunktion
mit einer systematischen Darstellung der Analyse nichtlinearer dynamischer Netzwerke.
PhD thesis, Universitat Kaiserslautern, 1976.

[KRS92] W. Kampowski, P. Rentrop, and W. Schmidt. Classiffcation and 
numerical simulation of electric circuits. Surveys on Mathematics for 
Industry, 2(1):23--65, 1992.

[MM08] F. Mazzia and C. Magherini. Test Set for Initial Value Problem Solvers, 
release 2.4. Department
of Mathematics, University of Bari and INdAM, Research Unit of Bari, February
2008
}
\details{

  The default parameters are:
   M1 = 0.36, M2 = 0.151104, M3 = 0.075552, 
       L1 = 0.15, L2 = 0.30,  J1 = 0.002727, J2 = 0.0045339259,
       EE = 0.20e12, NUE= 0.30, BB = 0.0080,   HH = 0.0080,
       RHO= 7870.0,  GRAV= 0.0, OMEGA = 150.0 
  
  There are two default initial conditions - set with options(ini=x)
  
}
\keyword{utilities}

