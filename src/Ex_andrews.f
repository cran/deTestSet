c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     This file is derived from the Test Set for IVP solvers
c     http://www.dm.uniba.it/~testset/
c
c        Andrews'' squeezing mechanism (in index 3 formulation)
c        index 3 DAE of dimension 27
c
c     DISCLAIMER: see
c     http://www.dm.uniba.it/~testset/disclaimer.php
c
c     The most recent version of this source file can be found at
c     http://www.dm.uniba.it/~testset/src/problems/andrews.f
c
c     This is revision
c     $Id: andrews.F,v 1.2 2006/10/02 10:29:13 testset Exp $
c
c-----------------------------------------------------------------------

c----------------------------------------------------------------------
c     parameter common block initialisation
c----------------------------------------------------------------------

      SUBROUTINE andpar(daeparms)
      EXTERNAL daeparms
      
      double precision parms(42)
      common / andcom/ parms

        CALL daeparms(42,parms)
        
      END SUBROUTINE andpar

c-----------------------------------------------------------------------
c     residual function
c-----------------------------------------------------------------------

      SUBROUTINE andres(T,Y,YPRIME,CJ,DELTA,IERR,RPAR,IPAR)
      
      implicit none
      DOUBLE PRECISION T, Y(27), DELTA(27), YPRIME(27),RPAR(*), CJ
      INTEGER I, J, IERR, N, IPAR(*)
C
      N = 27
      CALL andfunc(N, T, Y, DELTA, rpar, ipar)
C
      DO J = 1,14
         DELTA(J) = YPRIME(J)-DELTA(J)
      ENDDO
      DO I=15,N
         DELTA(I) = -DELTA(I)
      ENDDO
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE andfunc(NEQN, T, Y, DY, Rpar, Ipar)
      implicit none

      DOUBLE PRECISION T, Y(27), DY(27), rpar(*)
      integer i,j,neqn,Ipar(*)
      DOUBLE PRECISION m1,m2,m3,m4,m5,m6,m7,xa,ya,xb,yb,xc,yc,c0,
     +     i1,i2,i3,i4,i5,i6,i7,d,da,e,ea,rr,ra,l0,
     +     ss,sa,sb,sc,sd,ta,tb,u,ua,ub,zf,zt,fa,mom
      common/andcom/ m1,m2,m3,m4,m5,m6,m7,xa,ya,xb,yb,xc,yc,c0,
     +     i1,i2,i3,i4,i5,i6,i7,d,da,e,ea,rr,ra,l0,
     +     ss,sa,sb,sc,sd,ta,tb,u,ua,ub,zf,zt,fa,mom

      DOUBLE PRECISION sibe,sith,siga,siph,side,siom,siep,
     +     cobe,coth,coga,coph,code,coom,coep,
     +     sibeth,siphde,siomep,cobeth,cophde,coomep,
     +     bep,thp,php,dep,omp,epp,
     +     m(7,7),ff(7),gp(6,7),g(6),xd,yd,lang,force,fx,fy

      sibe = dsin(y(1))
      sith = dsin(y(2))
      siga = dsin(y(3))
      siph = dsin(y(4))
      side = dsin(y(5))
      siom = dsin(y(6))
      siep = dsin(y(7))
c
      cobe = dcos(y(1))
      coth = dcos(y(2))
      coga = dcos(y(3))
      coph = dcos(y(4))
      code = dcos(y(5))
      coom = dcos(y(6))
      coep = dcos(y(7))
c
      sibeth = dsin(y(1)+y(2))
      siphde = dsin(y(4)+y(5))
      siomep = dsin(y(6)+y(7))
c
      cobeth = dcos(y(1)+y(2))
      cophde = dcos(y(4)+y(5))
      coomep = dcos(y(6)+y(7))
c
      bep = y(8)
      thp = y(9)
      php = y(11)
      dep = y(12)
      omp = y(13)
      epp = y(14)
c
      do 20 j = 1,7
         do 10 i = 1,7
            m(i,j) = 0d0
 10      continue
 20   continue
c
      m(1,1) = m1*ra**2 + m2*(rr**2-2*da*rr*coth+da**2) + i1 + i2
      m(2,1) = m2*(da**2-da*rr*coth) + i2
      m(2,2) = m2*da**2 + i2
      m(3,3) = m3*(sa**2+sb**2) + i3
      m(4,4) = m4*(e-ea)**2 + i4
      m(5,4) = m4*((e-ea)**2+zt*(e-ea)*siph) + i4
      m(5,5) = m4*(zt**2+2*zt*(e-ea)*siph+(e-ea)**2) + m5*(ta**2+tb**2)
     &     + i4 + i5
      m(6,6) = m6*(zf-fa)**2 + i6
      m(7,6) = m6*((zf-fa)**2-u*(zf-fa)*siom) + i6
      m(7,7) = m6*((zf-fa)**2-2*u*(zf-fa)*siom+u**2) + m7*(ua**2+ub**2)
     &     + i6 + i7

      do 40 j=2,7
         do 30 i=1,j-1
            m(i,j) = m(j,i)
 30      continue
 40   continue
c
      xd = sd*coga + sc*siga + xb
      yd = sd*siga - sc*coga + yb
      lang  = dsqrt ((xd-xc)**2 + (yd-yc)**2)
      force = - c0 * (lang - l0)/lang
      fx = force * (xd-xc)
      fy = force * (yd-yc)
      ff(1) = mom - m2*da*rr*thp*(thp+2*bep)*sith
      ff(2) = m2*da*rr*bep**2*sith
      ff(3) = fx*(sc*coga - sd*siga) + fy*(sd*coga + sc*siga)
      ff(4) = m4*zt*(e-ea)*dep**2*coph
      ff(5) = - m4*zt*(e-ea)*php*(php+2*dep)*coph
      ff(6) = - m6*u*(zf-fa)*epp**2*coom
      ff(7) = m6*u*(zf-fa)*omp*(omp+2*epp)*coom
c
      do 60 j=1,7
         do 50 i=1,6
            gp(i,j) = 0d0
 50      continue
 60   continue
c
      gp(1,1) = - rr*sibe + d*sibeth
      gp(1,2) = d*sibeth
      gp(1,3) = - ss*coga
      gp(2,1) = rr*cobe - d*cobeth
      gp(2,2) = - d*cobeth
      gp(2,3) = - ss*siga
      gp(3,1) = - rr*sibe + d*sibeth
      gp(3,2) = d*sibeth
      gp(3,4) = - e*cophde
      gp(3,5) = - e*cophde + zt*side
      gp(4,1) = rr*cobe - d*cobeth
      gp(4,2) = - d*cobeth
      gp(4,4) = - e*siphde
      gp(4,5) = - e*siphde - zt*code
      gp(5,1) = - rr*sibe + d*sibeth
      gp(5,2) = d*sibeth
      gp(5,6) = zf*siomep
      gp(5,7) = zf*siomep - u*coep
      gp(6,1) = rr*cobe - d*cobeth
      gp(6,2) = - d*cobeth
      gp(6,6) = - zf*coomep
      gp(6,7) = - zf*coomep - u*siep
c
      g(1) = rr*cobe - d*cobeth - ss*siga - xb
      g(2) = rr*sibe - d*sibeth + ss*coga - yb
      g(3) = rr*cobe - d*cobeth - e*siphde - zt*code - xa
      g(4) = rr*sibe - d*sibeth + e*cophde - zt*side - ya
      g(5) = rr*cobe - d*cobeth - zf*coomep - u*siep - xa
      g(6) = rr*sibe - d*sibeth - zf*siomep + u*coep - ya
c
      do 70 i=1,14
         dy(i) = y(i+7)
   70 continue

      do 100 i=15,21
         dy(i) = -ff(i-14)
         do 80 j=1,7
            dy(i) = dy(i)+m(i-14,j)*y(j+14)
   80    continue
         do 90 j=1,6
            dy(i) = dy(i)+gp(j,i-14)*y(j+21)
   90    continue
  100 continue
      do 110 i=22,27
         dy(i) = g(i-21)
  110 continue

      return
      end

c----------------------------------------------------------------------
c     jacobian function
c----------------------------------------------------------------------

      SUBROUTINE andjac(T,Y,YPRIME,DFDY,CON,RPAR,IPAR)
      INTEGER NEQN,MN
      PARAMETER (NEQN=27, MN = 27)
      DOUBLE PRECISION T, Y(NEQN), DFDY(MN,NEQN),CON


c-----------------------------------------------------------------------
c     the Jacobian computed here is an approximation, see p. 540 of
c     Hairer & Wanner `solving ordinary differential equations II'
c-----------------------------------------------------------------------
      DOUBLE PRECISION m1,m2,m3,m4,m5,m6,m7,xa,ya,xb,yb,xc,yc,c0,
     +     i1,i2,i3,i4,i5,i6,i7,d,da,e,ea,rr,ra,l0,
     +     ss,sa,sb,sc,sd,ta,tb,u,ua,ub,zf,zt,fa,mom
      common/andcom/ m1,m2,m3,m4,m5,m6,m7,xa,ya,xb,yb,xc,yc,c0,
     +     i1,i2,i3,i4,i5,i6,i7,d,da,e,ea,rr,ra,l0,
     +     ss,sa,sb,sc,sd,ta,tb,u,ua,ub,zf,zt,fa,mom

      DOUBLE PRECISION sibe,siga,siph,side,siom,siep,
     +     cobe,coth,coga,code,coep,
     +     sibeth,siphde,siomep,cobeth,cophde,coomep,
     +     m(7,7),gp(6,7)


      sibe = dsin(y(1))
      siga = dsin(y(3))
      siph = dsin(y(4))
      side = dsin(y(5))
      siom = dsin(y(6))
      siep = dsin(y(7))
c
      cobe = dcos(y(1))
      coth = dcos(y(2))
      coga = dcos(y(3))
      code = dcos(y(5))
      coep = dcos(y(7))
c
      sibeth = dsin(y(1)+y(2))
      siphde = dsin(y(4)+y(5))
      siomep = dsin(y(6)+y(7))
c
      cobeth = dcos(y(1)+y(2))
      cophde = dcos(y(4)+y(5))
      coomep = dcos(y(6)+y(7))
c
      do 51 j = 1,7
         do 52 i = 1,7
            m(i,j) = 0d0
 52      continue
 51      continue
c
      m(1,1) = m1*ra**2 + m2*(rr**2-2*da*rr*coth+da**2) + i1 +i2
      m(2,1) = m2*(da**2-da*rr*coth) + i2
      m(2,2) = m2*da**2 + i2
      m(3,3) = m3*(sa**2+sb**2) + i3
      m(4,4) = m4*(e-ea)**2 + i4
      m(5,4) = m4*((e-ea)**2+zt*(e-ea)*siph) + i4
      m(5,5) = m4*(zt**2+2*zt*(e-ea)*siph+(e-ea)**2) + m5*(ta**2+tb**2)
     +         + i4 + i5
      m(6,6) = m6*(zf-fa)**2 + i6
      m(7,6) = m6*((zf-fa)**2-u*(zf-fa)*siom) + i6
      m(7,7) = m6*((zf-fa)**2-2*u*(zf-fa)*siom+u**2) + m7*(ua**2+ub**2)
     +         + i6 + i7

      do 40 j=2,7
         do 30 i=1,j-1
            m(i,j) = m(j,i)
   30    continue
   40 continue
c
      do 60 j=1,7
         do 50 i=1,6
            gp(i,j) = 0d0
   50    continue
   60 continue
c
      gp(1,1) = - rr*sibe + d*sibeth
      gp(1,2) = d*sibeth
      gp(1,3) = - ss*coga
      gp(2,1) = rr*cobe - d*cobeth
      gp(2,2) = - d*cobeth
      gp(2,3) = - ss*siga
      gp(3,1) = - rr*sibe + d*sibeth
      gp(3,2) = d*sibeth
      gp(3,4) = - e*cophde
      gp(3,5) = - e*cophde + zt*side
      gp(4,1) = rr*cobe - d*cobeth
      gp(4,2) = - d*cobeth
      gp(4,4) = - e*siphde
      gp(4,5) = - e*siphde - zt*code
      gp(5,1) = - rr*sibe + d*sibeth
      gp(5,2) = d*sibeth
      gp(5,6) = zf*siomep
      gp(5,7) = zf*siomep - u*coep
      gp(6,1) = rr*cobe - d*cobeth
      gp(6,2) = - d*cobeth
      gp(6,6) = - zf*coomep
      gp(6,7) = - zf*coomep - u*siep
c
      do 80 j=1,neqn
         do 70 i=1,neqn
            dfdy(i,j) = 0d0
 70      continue
 80   continue
      do 90 i=1,14
         dfdy(i,i+7) = 1d0
 90   continue
      do 110 i=1,7
         do 100 j=1,7
            dfdy(14+j,14+i) = m(j,i)
 100     continue
 110  continue
      do 130 i=1,6
         do 120 j=1,7
            dfdy(14+j,21+i) = gp(i,j)
 120     continue
 130  continue
      do 150 i=1,7
         do 140 j=1,6
            dfdy(21+j,i) = gp(j,i)
 140     continue
 150  continue
C
      do i=1,neqn
         do  j=1,neqn
            dfdy(i,j) = -dfdy(i,j)
         enddo
      enddo
c compute pd = -df/dy + con*M
      do j=1,14
         dfdy(j,j) = 1.0d0/con+dfdy(j,j)
      enddo
c
      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine andsoln(neqn,y)
      integer neqn
      double precision y(neqn)
c
c computed at Cray C90, using Cray double precision:
c Solving Andrews` squeezing mechanism using PSIDE
c  User input:
c
c  give relative error tolerance: 1d-14
c  give absolute error tolerance: 1d-14
c
c  Integration characteristics:
c
c     number of integration steps        1112
c     number of accepted steps           1022
c     number of f evaluations           26000
c     number of Jacobian evaluations      173
c     number of LU decompositions        2056
c
c  CPU-time used:                          56.74 sec


      y(  1) =  0.1581077119629904d+2
      y(  2) = -0.1575637105984298d+2
      y(  3) =  0.4082224013073101d-1
      y(  4) = -0.5347301163226948d+0
      y(  5) =  0.5244099658805304d+0
      y(  6) =  0.5347301163226948d+0
      y(  7) =  0.1048080741042263d+1
      y(  8) =  0.1139920302151208d+4
      y(  9) = -0.1424379294994111d+4
      y( 10) =  0.1103291221937134d+2
      y( 11) =  0.1929337464421385d+2
      y( 12) =  0.5735699284790808d+0
      y( 13) = -0.1929337464421385d+2
      y( 14) =  0.3231791658026955d+0
      y( 15) = -0.2463176316945196d+5
      y( 16) =  0.5185037701610329d+5
      y( 17) =  0.3241025686413781d+6
      y( 18) =  0.5667493645176213d+6
      y( 19) =  0.1674362929479361d+5
      y( 20) = -0.5667493645176222d+6
      y( 21) =  0.9826520791458422d+4
      y( 22) =  0.1991753333731910d+3
      y( 23) = -0.2975531228015052d+2
      y( 24) =  0.2306654119098399d+2
      y( 25) =  0.3145271365475927d+2
      y( 26) =  0.2264249232082739d+2
      y( 27) =  0.1161740700019673d+2

      return
      end
