c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     This file is part of the Test Set for IVP solvers
c     http://www.dm.uniba.it/~testset/
c
c        Ring Modulator (ODE case)
c        ODE of dimension 15
c
c     DISCLAIMER: see
c     http://www.dm.uniba.it/~testset/disclaimer.php
c
c     The most recent version of this source file can be found at
c     http://www.dm.uniba.it/~testset/src/problems/ringmod.f
c
c     This is revision
c     $Id: ringmod.F,v 1.2 2006/10/02 10:29:14 testset Exp $
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c----------------------------------------------------------------------
c     parameter common block initialisation
c----------------------------------------------------------------------

      SUBROUTINE ringpar(deparms)

      EXTERNAL deparms
      double precision parms(16)
      common /ringcom/parms

      call deparms(16, parms)
      return
      end

c----------------------------------------------------------------------
c     residual function
c----------------------------------------------------------------------

      subroutine ringres(t,y,yprime,cj,f,ires,rpar,ipar)
      integer neqn,ires,ipar(*)
      parameter  (neqn=15)
      double precision t,y(neqn),yprime(neqn),f(neqn),rpar(*),cj

      call ringfuncierr(neqn,t,y,f,ires, rpar,ipar)
      do i = 1,15
        f(i) = yprime(i) - f(i)
      enddo
      end subroutine

c----------------------------------------------------------------------
c     derivative function
c----------------------------------------------------------------------

      subroutine ringfunc(neqn,t,y,f,rpar,ipar)
      integer neqn,ierr,ipar(*)
      double precision t,y(neqn),f(neqn),rpar(*)

      double precision c,cs,cp,r,rp,lh,ls1,ls2,ls3,rg1,rg2,rg3,ri,rc,
     +                 gamma,delta,pi,
     +                 uin1,uin2,ud1,ud2,ud3,ud4,qud1,qud2,qud3,qud4
      common /ringcom/c, cs, cp, r, rp, lh, ls1, ls2, ls3,
     +     rg1, rg2, rg3, ri, rc, gamma, delta
      parameter (pi=3.141592653589793238462643383d0)

      uin1   = 0.5d0*sin(2d3*pi*t)
      uin2   = 2d0*sin(2d4*pi*t)
      ud1    = +y(3)-y(5)-y(7)-uin2
      ud2    = -y(4)+y(6)-y(7)-uin2
      ud3    = +y(4)+y(5)+y(7)+uin2
      ud4    = -y(3)-y(6)+y(7)+uin2

c     prevent overflow
c     (native NEC SX double precision .le. 1d75)
c      if (delta*max(ud1,ud2,ud3,ud4).gt.172d0) then
c         ierr = -1
c         return
c      endif
c
c      (double precisione ieee .le. 1d304)
c      if (delta*max(ud1,ud2,ud3,ud4).gt.708d0) then
c         ierr = -1
c         WRITE(MSG, *)"AN ERROR OCCURRED in RING, at time", T
c         call rexit(MSG)
c         return
c      endif

      qud1   = gamma*(exp(delta*ud1)-1d0)
      qud2   = gamma*(exp(delta*ud2)-1d0)
      qud3   = gamma*(exp(delta*ud3)-1d0)
      qud4   = gamma*(exp(delta*ud4)-1d0)

      f(1)  = (y(8)-0.5d0*y(10)+0.5d0*y(11)+y(14)-y(1)/r)/c
      f(2)  = (y(9)-0.5d0*y(12)+0.5d0*y(13)+y(15)-y(2)/r)/c
      f(3)  = (y(10)-qud1+qud4)/cs
      f(4)  = (-y(11)+qud2-qud3)/cs
      f(5)  = (y(12)+qud1-qud3)/cs
      f(6)  = (-y(13)-qud2+qud4)/cs
      f(7)  = (-y(7)/rp+qud1+qud2-qud3-qud4)/cp
      f(8)  = -y(1)/lh
      f(9)  = -y(2)/lh
      f(10) = (0.5d0*y(1)-y(3)-rg2*y(10))/ls2
      f(11) = (-0.5d0*y(1)+y(4)-rg3*y(11))/ls3
      f(12) = (0.5d0*y(2)-y(5)-rg2*y(12))/ls2
      f(13) = (-0.5d0*y(2)+y(6)-rg3*y(13))/ls3
      f(14) = (-y(1)+uin1-(ri+rg1)*y(14))/ls1
      f(15) = (-y(2)-(rc+rg1)*y(15))/ls1

      return
      end


      subroutine ringfuncierr(neqn,t,y,f,ierr,rpar,ipar)
      integer neqn,ierr,ipar(*)
      double precision t,y(neqn),f(neqn),rpar(*)

      double precision c,cs,cp,r,rp,lh,ls1,ls2,ls3,rg1,rg2,rg3,ri,rc,
     +                 gamma,delta,pi,
     +                 uin1,uin2,ud1,ud2,ud3,ud4,qud1,qud2,qud3,qud4
      common /ringcom/c, cs, cp, r, rp, lh, ls1, ls2, ls3,
     +     rg1, rg2, rg3, ri, rc, gamma, delta
      parameter (pi=3.141592653589793238462643383d0)

      uin1   = 0.5d0*sin(2d3*pi*t)
      uin2   = 2d0*sin(2d4*pi*t)
      ud1    = +y(3)-y(5)-y(7)-uin2
      ud2    = -y(4)+y(6)-y(7)-uin2
      ud3    = +y(4)+y(5)+y(7)+uin2
      ud4    = -y(3)-y(6)+y(7)+uin2

c     prevent overflow
c     (native NEC SX double precision .le. 1d75)
c      if (delta*max(ud1,ud2,ud3,ud4).gt.172d0) then
c         ierr = -1
c         return
c      endif
c
c      (double precisione ieee .le. 1d304)
      if (delta*max(ud1,ud2,ud3,ud4).gt. 304d0) then
         ierr = -1
c         WRITE(MSG, *)"AN ERROR OCCURRED in RING, at time", T
c         call rexit(MSG)
         return
      endif

      qud1   = gamma*(exp(delta*ud1)-1d0)
      qud2   = gamma*(exp(delta*ud2)-1d0)
      qud3   = gamma*(exp(delta*ud3)-1d0)
      qud4   = gamma*(exp(delta*ud4)-1d0)

      f(1)  = (y(8)-0.5d0*y(10)+0.5d0*y(11)+y(14)-y(1)/r)/c
      f(2)  = (y(9)-0.5d0*y(12)+0.5d0*y(13)+y(15)-y(2)/r)/c
      f(3)  = (y(10)-qud1+qud4)/cs
      f(4)  = (-y(11)+qud2-qud3)/cs
      f(5)  = (y(12)+qud1-qud3)/cs
      f(6)  = (-y(13)-qud2+qud4)/cs
      f(7)  = (-y(7)/rp+qud1+qud2-qud3-qud4)/cp
      f(8)  = -y(1)/lh
      f(9)  = -y(2)/lh
      f(10) = (0.5d0*y(1)-y(3)-rg2*y(10))/ls2
      f(11) = (-0.5d0*y(1)+y(4)-rg3*y(11))/ls3
      f(12) = (0.5d0*y(2)-y(5)-rg2*y(12))/ls2
      f(13) = (-0.5d0*y(2)+y(6)-rg3*y(13))/ls3
      f(14) = (-y(1)+uin1-(ri+rg1)*y(14))/ls1
      f(15) = (-y(2)-(rc+rg1)*y(15))/ls1

      return
      end
c----------------------------------------------------------------------
c     solution at default settings
c----------------------------------------------------------------------

      subroutine ringsoln(neqn,y)
      integer neqn
      double precision y(neqn)
c
c determined with RADAU
c
c User input:
c
c give relative error tolerance: 1d-12
c give absolute error tolerance: 1d-12
c give initial stepsize: 1d-14
c
c Integration characteristics:
c
c    number of integration steps       27195
c    number of accepted steps          23541
c    number of f evaluations          632631
c    number of Jacobian evaluations     9468
c    number of LU decompositions       22386

      y(  1) = -0.2339057358486745d-001
      y(  2) = -0.7367485485540825d-002
      y(  3) =  0.2582956709291169d+000
      y(  4) = -0.4064465721283450d+000
      y(  5) = -0.4039455665149794d+000
      y(  6) =  0.2607966765422943d+000
      y(  7) =  0.1106761861269975d+000
      y(  8) =  0.2939904342435596d-006
      y(  9) = -0.2840029933642329d-007
      y( 10) =  0.7267198267264553d-003
      y( 11) =  0.7929487196960840d-003
      y( 12) = -0.7255283495698965d-003
      y( 13) = -0.7941401968526521d-003
      y( 14) =  0.7088495416976114d-004
      y( 15) =  0.2390059075236570d-004

      return
      end

