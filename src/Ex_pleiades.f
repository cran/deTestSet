c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     This file is part of the Test Set for IVP solvers
c     http://www.dm.uniba.it/~testset/
c
c        PLEI problem
c        ODE of dimension 28
c
c     DISCLAIMER: see
c     http://www.dm.uniba.it/~testset/disclaimer.php
c
c     The most recent version of this source file can be found at
c     http://www.dm.uniba.it/~testset/src/problems/plei.f
c
c     This is revision
c     $Id: plei.F,v 1.2 2006/10/02 10:29:14 testset Exp $
c
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     derivative function
c-----------------------------------------------------------------------
      subroutine pleiafunc(neqn,t,y,f,rpar,ipar)
      integer neqn,ierr,ipar(*)
      double precision t,y(neqn),yprime(neqn),f(neqn),rpar(*)

      integer mj,i,j
      double precision sumx,sumy,rij,rij32

      do 20 i=1,7
         sumx = 0d0
         sumy = 0d0
         do 10 j=1,7
            mj = j
            rij = (y(i)-y(j))**2+(y(i+7)-y(j+7))**2
            rij32 = rij**(3d0/2d0)
            if (j.ne.i) then
               sumx = sumx+mj*(y(j)-y(i))/rij32
               sumy = sumy+mj*(y(j+7)-y(i+7))/rij32
            endif
   10    continue
         f(i+14) = sumx
         f(i+21) = sumy
   20 continue
      do 30 i=1,14
         f(i) = y(i+14)
   30 continue

      return
      end
c-----------------------------------------------------------------------
c     jacobian function
c-----------------------------------------------------------------------
      subroutine pleiajac(neqn,t,y,ml,mu,dfdy,ldim,rpar,ipar)
      integer ldim,neqn,ierr,ipar(*),ml,mu
      double precision t,y(neqn),dfdy(ldim,neqn),rpar(*)

      integer i,j,mi,mj
      double precision rij,rij32,rij52,fjh,sumxx,sumxy,sumyy

      do 20 i=1,neqn
         do 10 j=1,neqn
            dfdy(i,j)=0d0
   10    continue
   20 continue
      do 30 i=1,14
         dfdy(i,14+i)=1d0
   30 continue
      do 50 i=2,7
         mi=i
         do 40 j=1,i-1
            mj=j
            rij=(y(i)-y(j))**2+(y(i+7)-y(j+7))**2
            rij32=rij**1.5d0
            rij52=rij**2.5d0
            fjh=(1d0-3*(y(j)-y(i))**2/rij)/rij32
            dfdy(i+14,j)=mj*fjh
            dfdy(j+14,i)=mi*fjh
            fjh=(1d0-3*(y(j+7)-y(i+7))**2/rij)/rij32
            dfdy(i+21,j+7)=mj*fjh
            dfdy(j+21,i+7)=mi*fjh
            fjh=-3*(y(j)-y(i))*(y(j+7)-y(i+7))/rij52
            dfdy(i+14,j+7)=mj*fjh
            dfdy(j+14,i+7)=mi*fjh
            dfdy(i+21,j)=mj*fjh
            dfdy(j+21,i)=mi*fjh
   40    continue
   50 continue
      do 70 i=1,7
         sumxx=0d0
         sumxy=0d0
         sumyy=0d0
         do 60 j=1,7
            if (j.ne.i) then
               sumxx=sumxx+dfdy(i+14,j)
               sumxy=sumxy+dfdy(i+14,j+7)
               sumyy=sumyy+dfdy(i+21,j+7)
            endif
   60    continue
         dfdy(i+14,i)=-sumxx
         dfdy(i+14,i+7)=-sumxy
         dfdy(i+21,i)=-sumxy
         dfdy(i+21,i+7)=-sumyy
   70 continue
      return
      end

c-----------------------------------------------------------------------
      subroutine pleiasoln(neqn, y)
      integer neqn
      double precision y(28)
c
c computed at Cray C90 using Cray double precision
c Solving Pleiades problem using PSIDE
c
c User input:
c
c give relative error tolerance: 1d-16
c give absolute error tolerance: 1d-16
c
c
c Integration characteristics:
c
c    number of integration steps       12401
c    number of accepted steps          12397
c    number of f evaluations          111573
c    number of Jacobian evaluations        1
c    number of LU decompositions         712
c
c CPU-time used:                         343.29 sec

      y(  1) =  0.3706139143970502d+000
      y(  2) =  0.3237284092057233d+001
      y(  3) = -0.3222559032418324d+001
      y(  4) =  0.6597091455775310d+000
      y(  5) =  0.3425581707156584d+000
      y(  6) =  0.1562172101400631d+001
      y(  7) = -0.7003092922212495d+000
      y(  8) = -0.3943437585517392d+001
      y(  9) = -0.3271380973972550d+001
      y( 10) =  0.5225081843456543d+001
      y( 11) = -0.2590612434977470d+001
      y( 12) =  0.1198213693392275d+001
      y( 13) = -0.2429682344935824d+000
      y( 14) =  0.1091449240428980d+001
      y( 15) =  0.3417003806314313d+001
      y( 16) =  0.1354584501625501d+001
      y( 17) = -0.2590065597810775d+001
      y( 18) =  0.2025053734714242d+001
      y( 19) = -0.1155815100160448d+001
      y( 20) = -0.8072988170223021d+000
      y( 21) =  0.5952396354208710d+000
      y( 22) = -0.3741244961234010d+001
      y( 23) =  0.3773459685750630d+000
      y( 24) =  0.9386858869551073d+000
      y( 25) =  0.3667922227200571d+000
      y( 26) = -0.3474046353808490d+000
      y( 27) =  0.2344915448180937d+001
      y( 28) = -0.1947020434263292d+001

      return
      end

