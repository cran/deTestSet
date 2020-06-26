c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     This file is part of the Test Set for IVP solvers
c     http://www.dm.uniba.it/~testset/
c
c        Transistor Amplifier
c        index 1 DAE of dimension 8
c
c     DISCLAIMER: see
c     http://www.dm.uniba.it/~testset/disclaimer.php
c
c     The most recent version of this source file can be found at
c     http://www.dm.uniba.it/~testset/src/problems/transamp.f
c
c     This is revision
c     $Id: transamp.F,v 1.3 2006/10/25 08:21:22 testset Exp $
c
c-----------------------------------------------------------------------

      SUBROUTINE transpar(daeparms)

      EXTERNAL daeparms
      double precision parms(19)
      common /transcom/parms

      call daeparms(19, parms)
      return
      end


c----------------------------------------------------------------------
c     residual function
c----------------------------------------------------------------------

      SUBROUTINE transres (T,Y,YPRIME,cj,DELTA,IERR,RPAR,IPAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PARAMETER (NEQN = 8)
      DIMENSION Y(NEQN),DELTA(NEQN),YPRIME(NEQN),IPAR(*),RPAR(*)
      double precision ub,uf,alpha,beta,r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,
     +                 pi, c1,c2,c3,c4,c5
      parameter (pi=3.1415926535897931086244d0)

      common /transcom/ub, uf, alpha, beta, r0, r1, r2, r3, r4,
     +          r5, r6, r7, r8, r9, c1, c2, c3, c4, c5

C
      CALL transFunc(NEQN,T,Y,DELTA,RPAR,IPAR)
C
      DELTA(1) = -C1*YPRIME(1)+C1*YPRIME(2) -DELTA(1)
      DELTA(2) =  C1*YPRIME(1)-C1*YPRIME(2) -DELTA(2)
      DELTA(3) = -C2*YPRIME(3)              -DELTA(3)
      DELTA(4) = -C3*YPRIME(4)+C3*YPRIME(5) -DELTA(4)
      DELTA(5) =  C3*YPRIME(4)-C3*YPRIME(5) -DELTA(5)
      DELTA(6) = -C4*YPRIME(6)              -DELTA(6)
      DELTA(7) = -C5*YPRIME(7)+C5*YPRIME(8) -DELTA(7)
      DELTA(8) =  C5*YPRIME(7)-C5*YPRIME(8) -DELTA(8)

      RETURN
      END

c----------------------------------------------------------------------------

      subroutine transfunc(neqn,t,y,dy,rpar,ipar)
      IMPLICIT double precision (A-H,O-Z)
      double precision t,y(neqn),dy(neqn), rpar(*)
        integer ipar(*)

      double precision ub,uf,alpha,beta,r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,
     +                 pi,uet,fac1,fac2,c1,c2,c3,c4,c5
      parameter (pi=3.1415926535897931086244d0)

      common /transcom/ub, uf, alpha, beta, r0, r1, r2, r3, r4,
     +          r5, r6, r7, r8, r9, c1, c2, c3, c4, c5

      uet   = 0.1d0*sin(200d0*pi*t)
      fac1  = beta*(exp((y(2)-y(3))/uf)-1d0)
      fac2  = beta*(exp((y(5)-y(6))/uf)-1d0)

      dy(1) = (y(1)-uet)/r0
      dy(2) = y(2)/r1+(y(2)-ub)/r2+(1d0-alpha)*fac1
      dy(3) = y(3)/r3-fac1
      dy(4) = (y(4)-ub)/r4+alpha*fac1
      dy(5) = y(5)/r5+(y(5)-ub)/r6+(1d0-alpha)*fac2
      dy(6) = y(6)/r7-fac2
      dy(7) = (y(7)-ub)/r8+alpha*fac2
      dy(8) = y(8)/r9

      return
      end
c----------------------------------------------------------------------
c     jacobian function
c----------------------------------------------------------------------


      subroutine transjac(neqn,t,y,ml,mu,dfdy,ldim,rpar,ipar)
      IMPLICIT double precision (A-H,O-Z)
      integer n,i,ipar(*)
        parameter(n=8, mebnd=8)
      double precision t,y(n),dfdy(mebnd,n),rpar(*)

      double precision ub,uf,alpha,beta,r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,
     +                 pi,c1,c2,c3,c4,c5
      parameter (pi=3.1415926535897931086244d0)

      common /transcom/ub, uf, alpha, beta, r0, r1, r2, r3, r4,
     +          r5, r6, r7, r8, r9, c1, c2, c3, c4, c5
c
      fac1p = beta*dexp((y(2)-y(3))/uf)/uf
      fac2p = beta*dexp((y(5)-y(6))/uf)/uf
      do 10 i=1,8
         dfdy(1,i) = 0d0
         dfdy(3,i) = 0d0
         dfdy(4,i) = 0d0
   10 continue

      dfdy(1,3) = -(1d0-alpha)*fac1p
      dfdy(1,6) = -(1d0-alpha)*fac2p
      dfdy(2,1) = 1d0/r0
      dfdy(2,2) = 1d0/r1+1d0/r2+(1d0-alpha)*fac1p
      dfdy(2,3) = 1d0/r3+fac1p
      dfdy(2,4) = 1d0/r4
      dfdy(2,5) = 1d0/r5+1d0/r6+(1d0-alpha)*fac2p
      dfdy(2,6) = 1d0/r7+fac2p
      dfdy(2,7) = 1d0/r8
      dfdy(2,8) = 1d0/r9
      dfdy(3,2) = -fac1p
      dfdy(3,3) = -alpha*fac1p
      dfdy(3,5) = -fac2p
      dfdy(3,6) = -alpha*fac2p
      dfdy(4,2) = alpha*fac1p
      dfdy(4,5) = alpha*fac2p
      return
      end

c-----------------------------------------------------------------------
      subroutine transsoln(neqn,y)
      integer neqn
      double precision y(neqn)
c
c computed on Cray C90 using Cray double precision
c Solving Transistor Amplifier using PSIDE
c
c User input:
c
c give relative error tolerance: 1d-14
c give absolute error tolerance: 1d-14
c
c
c Integration characteristics:
c
c    number of integration steps       16061
c    number of accepted steps          15824
c    number of f evaluations          401944
c    number of Jacobian evaluations      458
c    number of LU decompositions        4884
c
c CPU-time used:                         182.44 sec

      y(  1) = -0.5562145012262709d-002
      y(  2) =  0.3006522471903042d+001
      y(  3) =  0.2849958788608128d+001
      y(  4) =  0.2926422536206241d+001
      y(  5) =  0.2704617865010554d+001
      y(  6) =  0.2761837778393145d+001
      y(  7) =  0.4770927631616772d+001
      y(  8) =  0.1236995868091548d+001

      return
      end

