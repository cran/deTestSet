c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     This file is part of the Test Set for IVP solvers
c     http://www.dm.uniba.it/~testset/
c
c        Problem HIRES
c        ODE of dimension 8
c
c     DISCLAIMER: see
c     http://www.dm.uniba.it/~testset/disclaimer.php
c
c     The most recent version of this source file can be found at
c     http://www.dm.uniba.it/~testset/src/problems/hires.f
c
c     This is revision
c     $Id: hires.F,v 1.2 2006/10/02 10:29:14 testset Exp $
c
c-----------------------------------------------------------------------


c----------------------------------------------------------------------
c     parameter common block initialisation
c----------------------------------------------------------------------

      SUBROUTINE hirespar(deparms)
      EXTERNAL daeparms
      
      double precision parms(10)
      common / hirescom/ parms

        CALL deparms(10,parms)
        
      END SUBROUTINE hirespar

c-----------------------------------------------------------------------
c     derivative function
c-----------------------------------------------------------------------

      subroutine hiresfun(neqn,t,y,f,rpar,ipar)
      integer neqn,ierr,ipar(*)
      double precision t,y(neqn),yprime(neqn),f(neqn),rpar(*)
      double precision k1, k2, k3, k4, k5, k6, k7, k8, k9, oks 
      common/hirescom/ k1, k2, k3, k4, k5, k6, k7, k8, k9, oks 

      f(1) = -k1*y(1) + k2*y(2) + k6*y(3) + oks
      f(2) =  k1*y(1) - (k2 + k3)*y(2)
      f(3) = -(k6+k1)*y(3) + k2*y(4) + k5*y(5)
      f(4) = k3*y(2) + k1*y(3) - (k4+k2)*y(4)
      f(5) = -(k5+k1)*y(5) + k2*(y(6)+y(7))
      f(6) = -k7*y(6)*y(8) + k8*y(4) + k1*y(5) - k2*y(6) + k8*y(7)
      f(7) = k7*y(6)*y(8) - (k2+k8+k9)*y(7)
      f(8) = -f(7)

      return
      end
c-----------------------------------------------------------------------
c     jacobian function
c-----------------------------------------------------------------------
      subroutine hiresjac(neqn,t,y,dfdy,ml,mu,ldim,rpar,ipar)
      integer ldim,neqn,ierr,ipar(*),mu,ml
      double precision t,y(neqn),yprime(neqn),dfdy(ldim,neqn),rpar(*)
      double precision k1, k2, k3, k4, k5, k6, k7, k8, k9, oks 
      common/hirescom/ k1, k2, k3, k4, k5, k6, k7, k8, k9, oks 

      integer i,j

      do 20 j=1,neqn
         do 10 i=1,neqn
            dfdy(i,j)=0d0
   10    continue
   20 continue

      dfdy(1,1) = -k1
      dfdy(1,2) = k2
      dfdy(1,3) = k6
      dfdy(2,1) = k1
      dfdy(2,2) = -(k2+k3)
      dfdy(3,3) = -(k6+k1)
      dfdy(3,4) = k2
      dfdy(3,5) = k5
      dfdy(4,2) = k3
      dfdy(4,3) = k1
      dfdy(4,4) = -(k4+k2)
      dfdy(5,5) = -(k5+k1)
      dfdy(5,6) = k2
      dfdy(5,7) = k2
      dfdy(6,4) = k8
      dfdy(6,5) = k1
      dfdy(6,6) = -k7*y(8)-k2
      dfdy(6,7) = k8
      dfdy(6,8) = -k7*y(6)
      dfdy(7,6) = k7*y(8)
      dfdy(7,7) =- (k2+k8+k9)
      dfdy(7,8) = k7*y(6)
      dfdy(8,6) = -k7*y(8)
      dfdy(8,7) = (k2+k8+k9)
      dfdy(8,8) = -k7*y(6)

      return
      end
c-----------------------------------------------------------------------
      subroutine hiressol(neqn,t,y)
      integer neqn
      double precision t,y(neqn)
c
c computed using true double precision RADAU5 on Cray C90
c          uround = work(1) = 1.01d-19
c          rtol = atol = h0 = 1.1d-18
c
      y(1) = 0.7371312573325668d-3
      y(2) = 0.1442485726316185d-3
      y(3) = 0.5888729740967575d-4
      y(4) = 0.1175651343283149d-2
      y(5) = 0.2386356198831331d-2
      y(6) = 0.6238968252742796d-2
      y(7) = 0.2849998395185769d-2
      y(8) = 0.2850001604814231d-2

      return
      end

