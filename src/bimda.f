C Karline: changed call to sollu: passes only first element of ijob
C NOW:   call sollu(m,theta,ldlu,Z(1,1),mljac,mujac,ipvt,ijob(1))
C WAS:   call sollu(m,theta,ldlu,Z(1,1),mljac,mujac,ipvt,ijob)

C -----------------------------------------------------------------------------------
C     THIS MODULE IS PART OF THE CODE BIMD.
C     THE CODE BIMD NUMERICALLY SOLVES (STIFF) DIFFERENTIAL ODE 
C     PROBLEMS OR LINEARLY IMPLICIT DAE PROBLEMS OF INDEX UP TO 3 
C     WITH CONSTANT MASS MATRIX
C
C     Copyright (C)2005-2007   
C
C     Authors: CECILIA MAGHERINI (cecilia.magherini@ing.unipi.it)
C              LUIGI   BRUGNANO  (brugnano@math.unifi.it) 
C
C
C     This program is free software; you can redistribute it and/or
C     modify it under the terms of the GNU General Public License
C     as published by the Free Software Foundation; either version 2
C     of the License, or (at your option) any later version.
C
C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     Licensed under The GNU General Public License, Version 2 or later.
C       http://www.gnu.org/licenses/info/GPLv2orLater.html
C
C     You should have received a copy of the GNU General Public License
C     along with this program; if not, write to the Free Software
C     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
C     USA.
C -----------------------------------------------------------------------------------



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                         c
c     In case the ISNAN function is not supported.        c
c                                                         c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function isnan(A)
      implicit none
      double precision A, B
      logical X
      B = A
      X = ( A .EQ. B )
      isnan = (.NOT.X)
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     BLENDSTEP
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine blendstep4(m,y0,f0,Y,F,h,theta,ipvt,Z,gamma,
     &                      ldlu,mljac,mujac,ijob,imas,
     &                      ldmas,mlmas,mumas,M0,MZ)
c
c     Blended iteration for the 4th order method
c

      implicit none

      integer k
      parameter (k=3)

c     Input parameters
      integer m,ipvt(m),ldlu,mljac,mujac,ijob(2),imas,ldmas,mlmas,mumas
      double precision  y0(m),f0(m),theta(ldlu,m),gamma,h,M0(ldmas,m)

c     I/O parameters
      double precision Y(m,k),F(m,k),MZ(m,k)

c     Output parameters
      double precision Z(m,k)

c     Local variables
      integer i,j

C ---------------------------------------------------------------------------------------
C 4th order BIM
C ---------------------------------------------------------------------------------------
      double precision DA4_1_1,DA4_1_2,DA4_1_3,DA4_1_4,
     &     DA4_2_1,DA4_2_2,DA4_2_3,DA4_2_4,
     &     DA4_3_1,DA4_3_2,DA4_3_3,DA4_3_4,
     &     DB4_1_1,DB4_1_2,DB4_1_3,DB4_1_4,
     &     DB4_2_1,DB4_2_2,DB4_2_3,DB4_2_4,
     &     DB4_3_1,DB4_3_2,DB4_3_3,DB4_3_4,
     &     A24_1_1,A24_1_2,A24_1_3,A24_1_4,
     &     A24_2_1,A24_2_2,A24_2_3,A24_2_4,
     &     A24_3_1,A24_3_2,A24_3_3,A24_3_4,
     &     B24_1_1,B24_2_1,B24_3_1

      parameter(
     & DA4_1_1 = -102133D0/405D3,
     & DA4_1_2 =   98743D0/18D4,
     & DA4_1_3 =   -7387D0/225D2,
     & DA4_1_4 =  +51709D0/162D4,
     & DA4_2_1 = -950353D0/81D4,
     & DA4_2_2 =   +7387D0/9D3,
     & DA4_2_3 =   10613D0/18D3,
     & DA4_2_4 =  -96031D0/405D3,
     & DA4_3_1 =  -22613D0/3D4,
     & DA4_3_2 =  -22161D0/2D4,
     & DA4_3_3 =  +22161D0/1D4,
     & DA4_3_4 =  -21257D0/6D4)

      parameter(
     & A24_1_1 = -302867D0/405D3,
     & A24_1_2 =  81257D0/18D4,
     & A24_1_3 =  7387D0/225D2,
     & A24_1_4 = -51709D0/162D4,
     & A24_2_1 =  140353D0/81D4,
     & A24_2_2 = -7387D0/9D3,
     & A24_2_3 =  7387D0/18D3,
     & A24_2_4 =  96031D0/405D3,
     & A24_3_1 = -7387D0/3D4,
     & A24_3_2 =  22161D0/2D4,
     & A24_3_3 = -22161D0/1D4,
     & A24_3_4 =  81257D0/6D4)
 

      parameter(
     & DB4_1_1 =  919D0/135D2, 
     & DB4_1_2 =  4589D0/3D4,
     & DB4_1_3 = -37D0/12D1, 
     & DB4_1_4 =  3D0/4D1,
     & DB4_2_1 =  115387D0/27D4,
     & DB4_2_2 =  17D0/15D0,
     & DB4_2_3 = -6161D0/3D4, 
     & DB4_2_4 = -1D0/15D0,
     & DB4_3_1 =  3D0/8D0,
     & DB4_3_2 =  9D0/8D0,
     & DB4_3_3 =  9D0/8D0,
     & DB4_3_4 = -3637D0/1D4,
     & B24_1_1 =  7387D0/27D3,
     & B24_2_1 = -7387D0/27D4, 
     & B24_3_1 =  0D0)
C ---------------------------------------------------------------------------------------

      goto(10,20) imas+1

10    continue
c     ODE CASE

c       Z=[y0 Y]*(DA)'
      do i=1,m
            Z(i,1)=y0(i)* DA4_1_1+Y(i,1)*DA4_1_2+Y(i,2)*DA4_1_3+
     &           Y(i,3)*DA4_1_4
            Z(i,2)=y0(i)* DA4_2_1+Y(i,1)*DA4_2_2+Y(i,2)*DA4_2_3+
     &           Y(i,3)*DA4_2_4
            Z(i,3)=y0(i)* DA4_3_1+Y(i,1)*DA4_3_2+Y(i,2)*DA4_3_3+
     &           Y(i,3)*DA4_3_4
      end do

c     Z=Z-h*[f0 F]*(dB)'

      do i=1,m
           Z(i,1) = Z(i,1)-
     &           h*(f0(i)* DB4_1_1+F(i,1)*DB4_1_2+F(i,2)*DB4_1_3+
     &              F(i,3)*DB4_1_4)
           Z(i,2) = Z(i,2)-
     &           h*(f0(i)* DB4_2_1+F(i,1)*DB4_2_2+F(i,2)*DB4_2_3+
     &              F(i,3)*DB4_2_4)
           Z(i,3) = Z(i,3)-
     &           h*(f0(i)* DB4_3_1+F(i,1)*DB4_3_2+F(i,2)*DB4_3_3+
     &              F(i,3)*DB4_3_4)
      end do

c     theta*Z=Z
      call sollu(m,theta,ldlu,Z(1,1),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,2),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,3),mljac,mujac,ipvt,ijob(1))

c      Z=Z+[y0 Y]*(A2)'
      do i=1,m
            Z(i,1)=Z(i,1)+
     &           y0(i)* A24_1_1+Y(i,1)*A24_1_2+Y(i,2)*A24_1_3+
     &          Y(i,3)* A24_1_4
            Z(i,2)=Z(i,2)+
     &           y0(i)* A24_2_1+Y(i,1)*A24_2_2+Y(i,2)*A24_2_3+
     &          Y(i,3)* A24_2_4
            Z(i,3)=Z(i,3)+
     &           y0(i)* A24_3_1+Y(i,1)*A24_3_2+Y(i,2)*A24_3_3+
     &          Y(i,3)* A24_3_4
      end do

c     Z=Z-h*[f0 F]*(B2)'
      do i=1,m
          Z(i,1)=Z(i,1)-h*(f0(i)*B24_1_1 + F(i,1)*gamma)
          Z(i,2)=Z(i,2)-h*(f0(i)*B24_2_1 + F(i,2)*gamma)
          Z(i,3)=Z(i,3)-h*(f0(i)*B24_3_1 + F(i,3)*gamma)
      end do

c     theta*Z=Z
      call sollu(m,theta,ldlu,Z(1,1),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,2),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,3),mljac,mujac,ipvt,ijob(1))

      goto 30

20    continue
c     DAE CASE

c      Z=[y0 Y]*(DA)'
      do i=1,m
            Z(i,1)=y0(i)* DA4_1_1+Y(i,1)*DA4_1_2+Y(i,2)*DA4_1_3+
     &           Y(i,3)*DA4_1_4
            Z(i,2)=y0(i)* DA4_2_1+Y(i,1)*DA4_2_2+Y(i,2)*DA4_2_3+
     &           Y(i,3)*DA4_2_4
            Z(i,3)=y0(i)* DA4_3_1+Y(i,1)*DA4_3_2+Y(i,2)*DA4_3_3+
     &           Y(i,3)*DA4_3_4
      end do

c     MZ_i=M0*Z_i, i=1,2,3

      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,1),MZ(1,1),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,2),MZ(1,2),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,3),MZ(1,3),ijob(2))

c     MZ=MZ-h*[f0 F]*(dB)'
      do i=1,m
           MZ(i,1) = MZ(i,1)-
     &              h*(f0(i)* DB4_1_1+F(i,1)*DB4_1_2+F(i,2)*DB4_1_3+
     &                 F(i,3)*DB4_1_4)
           MZ(i,2) = MZ(i,2)-
     &              h*(f0(i)* DB4_2_1+F(i,1)*DB4_2_2+F(i,2)*DB4_2_3+
     &                 F(i,3)*DB4_2_4)
           MZ(i,3) = MZ(i,3)-
     &              h*(f0(i)* DB4_3_1+F(i,1)*DB4_3_2+F(i,2)*DB4_3_3+
     &                 F(i,3)*DB4_3_4)
      end do

c     theta*MZ=MZ
      call sollu(m,theta,ldlu,MZ(1,1),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,2),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,3),mljac,mujac,ipvt,ijob(1))

c      MZ=MZ+[y0 Y]*(A2)'
      do i=1,m
            MZ(i,1)=MZ(i,1)+
     &           y0(i)* A24_1_1+Y(i,1)*A24_1_2+Y(i,2)*A24_1_3+
     &          Y(i,3)* A24_1_4
            MZ(i,2)=MZ(i,2)+
     &           y0(i)* A24_2_1+Y(i,1)*A24_2_2+Y(i,2)*A24_2_3+
     &          Y(i,3)* A24_2_4
            MZ(i,3)=MZ(i,3)+
     &           y0(i)* A24_3_1+Y(i,1)*A24_3_2+Y(i,2)*A24_3_3+
     &          Y(i,3)* A24_3_4
      end do

c     Z_i=M0*MZ_i, i=1,2,3

      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,1),Z(1,1),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,2),Z(1,2),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,3),Z(1,3),ijob(2))

c     Z=Z-h*[f0 F]*(B2)'
      do i=1,m
         Z(i,1)=Z(i,1)-h*(f0(i)*B24_1_1 + F(i,1)*gamma)
         Z(i,2)=Z(i,2)-h*(f0(i)*B24_2_1 + F(i,2)*gamma)
         Z(i,3)=Z(i,3)-h*(f0(i)*B24_3_1 + F(i,3)*gamma)
      end do

c     theta*Z=Z
      call sollu(m,theta,ldlu,Z(1,1),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,2),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,3),mljac,mujac,ipvt,ijob(1))

30    continue

c
c     Y = Y-Z
c
      do j=1,k
         do i=1,m
           Y(i,j) = Y(i,j) - Z(i,j)
         end do
      end do

      return
      end

C -------------------------------------------------------------------

      subroutine blendstep6(m,y0,f0,Y,F,h,theta,ipvt,Z,gamma,
     &                      ldlu,mljac,mujac,ijob,imas,
     &                                ldmas,mlmas,mumas,M0,MZ)
c
c     Blended iteration for the 6th order method
c

      implicit none

      integer k
      parameter (k=4)

c     Input parameters
      integer m,ipvt(m),ldlu,mljac,mujac,ijob(2),imas,ldmas,mlmas,mumas
      double precision  y0(m),f0(m),theta(ldlu,m),gamma,h,M0(ldmas,m)

c     I/O parameters
      double precision Y(m,k),F(m,k),MZ(m,k)

c     Output parameters
      double precision Z(m,k)

c     Local variables
      integer i,j
C ---------------------------------------------------------------------------------------
C 6th order BIM
C ---------------------------------------------------------------------------------------

      double precision DA6_1_1,DA6_1_2,DA6_1_3,DA6_1_4,DA6_1_5,
     &     DA6_2_1,DA6_2_2,DA6_2_3,DA6_2_4,DA6_2_5,
     &     DA6_3_1,DA6_3_2,DA6_3_3,DA6_3_4,DA6_3_5,
     &     DA6_4_1,DA6_4_2,DA6_4_3,DA6_4_4,DA6_4_5,
     &     DB6_1_1,DB6_1_2,DB6_1_3,DB6_1_4,DB6_1_5,
     &     DB6_2_1,DB6_2_2,DB6_2_3,DB6_2_4,DB6_2_5,
     &     DB6_3_1,DB6_3_2,DB6_3_3,DB6_3_4,DB6_3_5,
     &     DB6_4_1,DB6_4_2,DB6_4_3,DB6_4_4,DB6_4_5,
     &     A26_1_1,A26_1_2,A26_1_3,A26_1_4,A26_1_5,
     &     A26_2_1,A26_2_2,A26_2_3,A26_2_4,A26_2_5,
     &     A26_3_1,A26_3_2,A26_3_3,A26_3_4,A26_3_5,
     &     A26_4_1,A26_4_2,A26_4_3,A26_4_4,A26_4_5,
     &     B26_1_1,B26_2_1,B26_3_1,B26_4_1


      parameter(
     &  DA6_1_1 = -1171629D0/512D4,
     &  DA6_1_2 =  607997D0/96D4,
     &  DA6_1_3 = -597981D0/128D4,
     &  DA6_1_4 = +4241D0/64D3,
     &  DA6_1_5 = -55133D0/1536D4,
     &  DA6_2_1 = -307277D0/32D4,
     &  DA6_2_2 = +4241D0/12D3,
     &  DA6_2_3 = +92723D0/8D4,
     &  DA6_2_4 = -12723D0/2D4,
     &  DA6_2_5 = +80579D0/96D4,
     &  DA6_3_1 = -5853693D0/512D4,
     &  DA6_3_2 = -4241D0/32D4,
     &  DA6_3_3 = +1234131D0/128D4,
     &  DA6_3_4 =  137637D0/32D4,
     &  DA6_3_5 = -1217167D0/512D4,
     &  DA6_4_1 = -24241D0/2D4,
     &  DA6_4_2 = +4241D0/375D1,
     &  DA6_4_3 = -12723D0/5D3,
     &  DA6_4_4 = +4241D0/125D1,
     &  DA6_4_5 = -1841D0/24D2)


      parameter(
     &  A26_1_1 = -3948371D0/512D4,
     &  A26_1_2 = 352003D0/96D4,
     &  A26_1_3 = 597981D0/128D4,
     &  A26_1_4 = -4241D0/64D3,
     &  A26_1_5 = 55133D0/1536D4,
     &  A26_2_1 = -12723D0/32D4,
     &  A26_2_2 = -4241D0/12D3,
     &  A26_2_3 = -12723D0/8D4,
     &  A26_2_4 = 12723D0/2D4,
     &  A26_2_5 = -80579D0/96D4,
     &  A26_3_1 = 733693D0/512D4,
     &  A26_3_2 = 4241D0/32D4,
     &  A26_3_3 = -1234131D0/128D4,
     &  A26_3_4 =  182363D0/32D4,
     &  A26_3_5 =  1217167D0/512D4,
     &  A26_4_1 =  4241D0/2D4,
     &  A26_4_2 = -4241D0/375D1,
     &  A26_4_3 =  12723D0/5D3,
     &  A26_4_4 = -4241D0/125D1,
     &  A26_4_5 =  4241D0/24D2)


      parameter(
     &  DB6_1_1 =  60311D0/1152D4,
     &  DB6_1_2 =  7853D0/225D2,
     &  DB6_1_3 = -49D0/6D1,
     &  DB6_1_4 =  161D0/36D1,
     &  DB6_1_5 = -73D0/72D1,
     &  DB6_2_1 =  257831D0/72D4,
     &  DB6_2_2 =  46D0/45D0,
     &  DB6_2_3 = -241D0/5D3,
     &  DB6_2_4 = -14D0/45D0,
     &  DB6_2_5 =  7D0/9D1,
     &  DB6_3_1 =  1850413D0/384D4,
     &  DB6_3_2 =  133D0/12D1,
     &  DB6_3_3 =  23D0/2D1,
     &  DB6_3_4 = -1837D0/375D1,
     &  DB6_3_5 =  1D0/24D1,
     &  DB6_4_1 =  14D0/45D0,
     &  DB6_4_2 =  64D0/45D0,
     &  DB6_4_3 =  8D0/15D0,
     &  DB6_4_4 =  64D0/45D0,
     &  DB6_4_5 = -24169D0/45D3,
     &  B26_1_1 =  343521D0/128D4,
     &  B26_2_1 =  4241D0/8D4,
     &  B26_3_1 = -131471D0/128D4,
     &  B26_4_1 =  0D0)

C ---------------------------------------------------------------------------------------

      goto(10,20) imas+1

10    continue
c     ODE case

c      Z=[y0 Y]*(DA)'
      do i=1,m
            Z(i,1)=y0(i)* DA6_1_1+Y(i,1)*DA6_1_2+Y(i,2)*DA6_1_3+
     &           Y(i,3)*DA6_1_4+Y(i,4)*DA6_1_5
            Z(i,2)=y0(i)* DA6_2_1+Y(i,1)*DA6_2_2+Y(i,2)*DA6_2_3+
     &           Y(i,3)*DA6_2_4+Y(i,4)*DA6_2_5
            Z(i,3)=y0(i)* DA6_3_1+Y(i,1)*DA6_3_2+Y(i,2)*DA6_3_3+
     &           Y(i,3)*DA6_3_4+Y(i,4)*DA6_3_5
            Z(i,4)=y0(i)* DA6_4_1+Y(i,1)*DA6_4_2+Y(i,2)*DA6_4_3+
     &           Y(i,3)*DA6_4_4+Y(i,4)*DA6_4_5
      end do

c     Z=Z-h*[f0 F]*(dB)'
      do i=1,m
          Z(i,1)=Z(i,1)-
     &           h*(f0(i)* DB6_1_1+F(i,1)*DB6_1_2+F(i,2)*DB6_1_3+
     &             F(i,3)*DB6_1_4+F(i,4)*DB6_1_5)
          Z(i,2)=Z(i,2)-
     &           h*(f0(i)* DB6_2_1+F(i,1)*DB6_2_2+F(i,2)*DB6_2_3+
     &              F(i,3)*DB6_2_4+F(i,4)*DB6_2_5)
          Z(i,3)=Z(i,3)-
     &           h*(f0(i)* DB6_3_1+F(i,1)*DB6_3_2+F(i,2)*DB6_3_3+
     &              F(i,3)*DB6_3_4+F(i,4)*DB6_3_5)
          Z(i,4)=Z(i,4)-
     &           h*(f0(i)* DB6_4_1+F(i,1)*DB6_4_2+F(i,2)*DB6_4_3+
     &              F(i,3)*DB6_4_4+F(i,4)*DB6_4_5)
      end do

c     theta*Z=Z
      call sollu(m,theta,ldlu,Z(1,1),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,2),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,3),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,4),mljac,mujac,ipvt,ijob(1))

c      Z=Z+[y0 Y]*(A2)'
      do i=1,m
            Z(i,1)=Z(i,1)+
     &           y0(i)* A26_1_1+Y(i,1)*A26_1_2+Y(i,2)*A26_1_3+
     &          Y(i,3)* A26_1_4+Y(i,4)*A26_1_5
            Z(i,2)=Z(i,2)+
     &           y0(i)* A26_2_1+Y(i,1)*A26_2_2+Y(i,2)*A26_2_3+
     &          Y(i,3)* A26_2_4+Y(i,4)*A26_2_5
            Z(i,3)=Z(i,3)+
     &           y0(i)* A26_3_1+Y(i,1)*A26_3_2+Y(i,2)*A26_3_3+
     &          Y(i,3)* A26_3_4+Y(i,4)*A26_3_5
            Z(i,4)=Z(i,4)+
     &           y0(i)* A26_4_1+Y(i,1)*A26_4_2+Y(i,2)*A26_4_3+
     &          Y(i,3)* A26_4_4+Y(i,4)*A26_4_5
      end do

c     Z=Z-h*[f0 F]*(B2)'
      do i=1,m
         Z(i,1)=Z(i,1)-h*(f0(i)*B26_1_1 + F(i,1)*gamma)
         Z(i,2)=Z(i,2)-h*(f0(i)*B26_2_1 + F(i,2)*gamma)
         Z(i,3)=Z(i,3)-h*(f0(i)*B26_3_1 + F(i,3)*gamma)
         Z(i,4)=Z(i,4)-h*(f0(i)*B26_4_1 + F(i,4)*gamma)
      end do

c     theta*Z=Z
      call sollu(m,theta,ldlu,Z(1,1),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,2),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,3),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,4),mljac,mujac,ipvt,ijob(1))

      goto 30

20    continue
c     DAE case

c      Z=[y0 Y]*(DA)'
      do i=1,m
            Z(i,1)=y0(i)* DA6_1_1+Y(i,1)*DA6_1_2+Y(i,2)*DA6_1_3+
     &           Y(i,3)*DA6_1_4+Y(i,4)*DA6_1_5
            Z(i,2)=y0(i)* DA6_2_1+Y(i,1)*DA6_2_2+Y(i,2)*DA6_2_3+
     &           Y(i,3)*DA6_2_4+Y(i,4)*DA6_2_5
            Z(i,3)=y0(i)* DA6_3_1+Y(i,1)*DA6_3_2+Y(i,2)*DA6_3_3+
     &           Y(i,3)*DA6_3_4+Y(i,4)*DA6_3_5
            Z(i,4)=y0(i)* DA6_4_1+Y(i,1)*DA6_4_2+Y(i,2)*DA6_4_3+
     &           Y(i,3)*DA6_4_4+Y(i,4)*DA6_4_5
      end do

c     MZ_i=M0*Z_i, i=1,2,3,4

      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,1),MZ(1,1),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,2),MZ(1,2),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,3),MZ(1,3),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,4),MZ(1,4),ijob(2))

c     MZ=MZ-h*[f0 F]*(dB)'
      do i=1,m
         MZ(i,1)=MZ(i,1)-
     &          h*(f0(i)* DB6_1_1+F(i,1)*DB6_1_2+F(i,2)*DB6_1_3+
     &             F(i,3)*DB6_1_4+F(i,4)*DB6_1_5)
         MZ(i,2)=MZ(i,2)-
     &          h*(f0(i)* DB6_2_1+F(i,1)*DB6_2_2+F(i,2)*DB6_2_3+
     &             F(i,3)*DB6_2_4+F(i,4)*DB6_2_5)
         MZ(i,3)=MZ(i,3)-
     &          h*(f0(i)* DB6_3_1+F(i,1)*DB6_3_2+F(i,2)*DB6_3_3+
     &             F(i,3)*DB6_3_4+F(i,4)*DB6_3_5)
         MZ(i,4)=MZ(i,4)-
     &          h*(f0(i)* DB6_4_1+F(i,1)*DB6_4_2+F(i,2)*DB6_4_3+
     &             F(i,3)*DB6_4_4+F(i,4)*DB6_4_5)
      end do

c     theta*MZ=MZ
      call sollu(m,theta,ldlu,MZ(1,1),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,2),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,3),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,4),mljac,mujac,ipvt,ijob(1))

c      MZ=MZ+[y0 Y]*(A2)'
      do i=1,m
            MZ(i,1)=MZ(i,1)+
     &           y0(i)* A26_1_1+Y(i,1)*A26_1_2+Y(i,2)*A26_1_3+
     &          Y(i,3)* A26_1_4+Y(i,4)*A26_1_5
            MZ(i,2)=MZ(i,2)+
     &           y0(i)* A26_2_1+Y(i,1)*A26_2_2+Y(i,2)*A26_2_3+
     &          Y(i,3)* A26_2_4+Y(i,4)*A26_2_5
            MZ(i,3)=MZ(i,3)+
     &           y0(i)* A26_3_1+Y(i,1)*A26_3_2+Y(i,2)*A26_3_3+
     &          Y(i,3)* A26_3_4+Y(i,4)*A26_3_5
            MZ(i,4)=MZ(i,4)+
     &           y0(i)* A26_4_1+Y(i,1)*A26_4_2+Y(i,2)*A26_4_3+
     &          Y(i,3)* A26_4_4+Y(i,4)*A26_4_5
      end do

c     Z_i=M0*MZ_i, i=1,2,3,4

      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,1),Z(1,1),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,2),Z(1,2),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,3),Z(1,3),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,4),Z(1,4),ijob(2))

c     Z=Z-h*[f0 F]*(B2)'
      do i=1,m
           Z(i,1)=Z(i,1)-h*(f0(i)*B26_1_1 + F(i,1)*gamma)
           Z(i,2)=Z(i,2)-h*(f0(i)*B26_2_1 + F(i,2)*gamma)
           Z(i,3)=Z(i,3)-h*(f0(i)*B26_3_1 + F(i,3)*gamma)
           Z(i,4)=Z(i,4)-h*(f0(i)*B26_4_1 + F(i,4)*gamma)
      end do

c     theta*Z=Z
      call sollu(m,theta,ldlu,Z(1,1),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,2),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,3),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,4),mljac,mujac,ipvt,ijob(1))

30    continue

c
c     Y = Y-Z
c
      do j=1,k
         do i=1,m
           Y(i,j) = Y(i,j) - Z(i,j)
         end do
      end do

      return
      end

C -------------------------------------------------------------------

      subroutine blendstep8(m,y0,f0,Y,F,h,theta,ipvt,Z,gamma,
     &                            ldlu,mljac,mujac,ijob,imas,
     &                               ldmas,mlmas,mumas,M0,MZ)
c
c     Blended iteration for the 8th order method
c

      implicit none

      integer k
      parameter (k=6)

c     Input parameters
      integer m,ipvt(m),ldlu,mljac,mujac,ijob(2),imas,ldmas,mlmas,mumas
      double precision  y0(m),f0(m),theta(ldlu,m),gamma,h,M0(ldmas,m)
c     I/O parameters
      double precision Y(m,k),F(m,k),MZ(m,k)

c     Output parameters
      double precision Z(m,k)

c     Local variables
      integer i,j
C ---------------------------------------------------------------------------------------
C 8th order BIM
C ---------------------------------------------------------------------------------------

      double precision
     &     DA8_1_1,DA8_1_2,DA8_1_3,DA8_1_4,DA8_1_5,DA8_1_6,DA8_1_7,
     &     DA8_2_1,DA8_2_2,DA8_2_3,DA8_2_4,DA8_2_5,DA8_2_6,DA8_2_7,
     &     DA8_3_1,DA8_3_2,DA8_3_3,DA8_3_4,DA8_3_5,DA8_3_6,DA8_3_7,
     &     DA8_4_1,DA8_4_2,DA8_4_3,DA8_4_4,DA8_4_5,DA8_4_6,DA8_4_7,
     &     DA8_5_1,DA8_5_2,DA8_5_3,DA8_5_4,DA8_5_5,DA8_5_6,DA8_5_7,
     &     DA8_6_1,DA8_6_2,DA8_6_3,DA8_6_4,DA8_6_5,DA8_6_6,DA8_6_7,
     &     DB8_1_1,DB8_1_2,DB8_1_3,DB8_1_4,DB8_1_5,DB8_1_6,DB8_1_7,
     &     DB8_2_1,DB8_2_2,DB8_2_3,DB8_2_4,DB8_2_5,DB8_2_6,DB8_2_7,
     &     DB8_3_1,DB8_3_2,DB8_3_3,DB8_3_4,DB8_3_5,DB8_3_6,DB8_3_7,
     &     DB8_4_1,DB8_4_2,DB8_4_3,DB8_4_4,DB8_4_5,DB8_4_6,DB8_4_7,
     &     DB8_5_1,DB8_5_2,DB8_5_3,DB8_5_4,DB8_5_5,DB8_5_6,DB8_5_7,
     &     DB8_6_1,DB8_6_2,DB8_6_3,DB8_6_4,DB8_6_5,DB8_6_6,DB8_6_7,
     &     A28_1_1,A28_1_2,A28_1_3,A28_1_4,A28_1_5,A28_1_6,A28_1_7,
     &     A28_2_1,A28_2_2,A28_2_3,A28_2_4,A28_2_5,A28_2_6,A28_2_7,
     &     A28_3_1,A28_3_2,A28_3_3,A28_3_4,A28_3_5,A28_3_6,A28_3_7,
     &     A28_4_1,A28_4_2,A28_4_3,A28_4_4,A28_4_5,A28_4_6,A28_4_7,
     &     A28_5_1,A28_5_2,A28_5_3,A28_5_4,A28_5_5,A28_5_6,A28_5_7,
     &     A28_6_1,A28_6_2,A28_6_3,A28_6_4,A28_6_5,A28_6_6,A28_6_7,
     &     B28_1_1,B28_2_1,B28_3_1,B28_4_1,B28_5_1,B28_6_1


      parameter(
     &  DA8_1_1  =   -24312887D0/62208D3,
     &  DA8_1_2  =     9595787D0/1296D4,
     &  DA8_1_3  =   - 680419D0/20736D2,
     &  DA8_1_4  =   - 263717D0/23328D2,
     &  DA8_1_5  =   + 578429D0/41472D2,
     &  DA8_1_6  =   - 147157D0/2592D3,
     &  DA8_1_7  =   + 4150993D0/46656D4,
     &  DA8_2_1  =   -1366847D0/972D3,   
     &  DA8_2_2  =   +496837D0/405D3,
     &  DA8_2_3  =    165733D0/648D3,
     &  DA8_2_4  =   +24769D0/3645D2,
     &  DA8_2_5  =   -71393D0/324D3,
     &  DA8_2_6  =   +1457D0/162D2,
     &  DA8_2_7  =   -403589D0/2916D4,
     &  DA8_3_1  =   -801511D0/768D3,
     &  DA8_3_2  =   +4371D0/16D4,
     &  DA8_3_3  =   +48081D0/128D3,
     &  DA8_3_4  =    11057D0/96D2,                
     &  DA8_3_5  =   -161727D0/256D3,
     &  DA8_3_6  =   +4371D0/32D3,
     &  DA8_3_7  =   -10199D0/64D4,
     &  DA8_4_1  =   -388381D0/486D3,
     &  DA8_4_2  =   -85963D0/2025D2,
     &  DA8_4_3  =   +2914D0/10125D0,
     &  DA8_4_4  =   +71393D0/18225D1,
     &  DA8_4_5  =    145973D0/162D3,
     &  DA8_4_6  =   -16027D0/405D2,
     &  DA8_4_7  =   +141329D0/3645D3,
     &  DA8_5_1  =   -64281311D0/62208D3,
     &  DA8_5_2  =   -106361D0/2592D3,
     &  DA8_5_3  =   +893141D0/20736D2,
     &  DA8_5_4  =   -2466701D0/23328D2,
     &  DA8_5_5  =   +7187381D0/41472D2,
     &  DA8_5_6  =    241859D0/2592D3,
     &  DA8_5_7  =   -11695339D0/93312D3,
     &  DA8_6_1  =   -13457D0/12D3,         
     &  DA8_6_2  =   +4371D0/5D3,
     &  DA8_6_3  =   -4371D0/16D2,
     &  DA8_6_4  =   +1457D0/3D2,
     &  DA8_6_5  =   -4371D0/8D2,
     &  DA8_6_6  =   +4371D0/1D3,
     &  DA8_6_7  =   -31393D0/4D4)



      parameter(
     &   A28_1_1  =   -37895113D0/62208D3,
     &   A28_1_2  =    3364213D0/1296D4,
     &   A28_1_3  =    680419D0/20736D2,
     &   A28_1_4  =    263717D0/23328D2,
     &   A28_1_5  =    -578429D0/41472D2,
     &   A28_1_6  =    147157D0/2592D3,
     &   A28_1_7  =    -4150993D0/46656D4,
     &   A28_2_1  =    394847D0/972D3,
     &   A28_2_2  =    -496837D0/405D3,
     &   A28_2_3  =    482267D0/648D3,
     &   A28_2_4  =    -24769D0/3645D2,
     &   A28_2_5  =    71393D0/324D3,
     &   A28_2_6  =    -1457D0/162D2,
     &   A28_2_7  =    403589D0/2916D4,
     &   A28_3_1  =    33511D0/768D3,
     &   A28_3_2  =    -4371D0/16D4,
     &   A28_3_3  =    -48081D0/128D3,
     &   A28_3_4  =    -1457D0/96D2,
     &   A28_3_5  =    161727D0/256D3,
     &   A28_3_6  =    -4371D0/32D3,
     &   A28_3_7  =    10199D0/64D4,
     &   A28_4_1  =    -97619D0/486D3,
     &   A28_4_2  =    85963D0/2025D2,
     &   A28_4_3  =    -2914D0/10125D0,
     &   A28_4_4  =    -71393D0/18225D1,
     &   A28_4_5  =    16027D0/162D3,
     &   A28_4_6  =    16027D0/405D2,
     &   A28_4_7  =    -141329D0/3645D3,
     &   A28_5_1  =    2073311D0/62208D3,
     &   A28_5_2  =    106361D0/2592D3,
     &   A28_5_3  =    -893141D0/20736D2,
     &   A28_5_4  =    2466701D0/23328D2,
     &   A28_5_5  =    -7187381D0/41472D2,
     &   A28_5_6  =    2350141D0/2592D3,
     &   A28_5_7  =    11695339D0/93312D3,
     &   A28_6_1  =    1457D0/12D3,
     &   A28_6_2  =    -4371D0/5D3,
     &   A28_6_3  =    4371D0/16D2,
     &   A28_6_4  =    -1457D0/3D2,
     &   A28_6_5  =    4371D0/8D2,
     &   A28_6_6  =    -4371D0/1D3,
     &   A28_6_7  =    71393D0/4D4)



      parameter(
     &   DB8_1_1  =    87907D0/62208D1,
     &   DB8_1_2  =    3587D0/18D3,
     &   DB8_1_3  =    -1141D0/288D1,
     &   DB8_1_4  =    67D0/54D1,
     &   DB8_1_5  =    109D0/288D1,
     &   DB8_1_6  =    -2D0/45D0,
     &   DB8_1_7  =    91D0/864D1,
     &   DB8_2_1  =    203071D0/42525D1,
     &   DB8_2_2  =    2158D0/1575D0,
     &   DB8_2_3  =    -52291D0/126D3,
     &   DB8_2_4  =    -52D0/945D0,
     &   DB8_2_5  =    23D0/252D0,
     &   DB8_2_6  =    -82D0/1575D0,
     &   DB8_2_7  =    199D0/189D2,
     &   DB8_3_1  =    19177D0/64D3,
     &   DB8_3_2  =    81D0/5D1,
     &   DB8_3_3  =    27D0/32D1,
     &   DB8_3_4  =    1643D0/2D3,
     &   DB8_3_5  =    -243D0/32D1,
     &   DB8_3_6  =    27D0/1D2,
     &   DB8_3_7  =    -67D0/16D2,
     &   DB8_4_1  =    716549D0/3402D3,
     &   DB8_4_2  =    2368D0/1575D0,
     &   DB8_4_3  =    104D0/315D0,
     &   DB8_4_4  =    320D0/189D0,
     &   DB8_4_5  =    -78191D0/126D3,
     &   DB8_4_6  =    128D0/1575D0,
     &   DB8_4_7  =    -64D0/4725D0,
     &   DB8_5_1  =    37053349D0/108864D3,
     &   DB8_5_2  =    1739D0/126D1,
     &   DB8_5_3  =    2713D0/4032D0,
     &   DB8_5_4  =    853D0/756D0,
     &   DB8_5_5  =    4463D0/4032D0,
     &   DB8_5_6  =    -40391D0/126D3,
     &   DB8_5_7  =    -787D0/6048D1,
     &   DB8_6_1  =    41D0/14D1,
     &   DB8_6_2  =    54D0/35D0,
     &   DB8_6_3  =    27D0/14D1,
     &   DB8_6_4  =    68D0/35D0,
     &   DB8_6_5  =    27D0/14D1,
     &   DB8_6_6  =    54D0/35D0,
     &   DB8_6_7  =    -6099D0/14D3,
     &   B28_1_1  =    24769D0/124416D0,
     &   B28_2_1  =    -18941D0/1215D2,
     &   B28_3_1  =    -1457D0/64D3,
     &   B28_4_1  =    42253D0/486D3,
     &   B28_5_1  =    -365707D0/15552D3,
     &   B28_6_1  =    0D0)




C ---------------------------------------------------------------------------------------

      goto(10,20) imas+1

10    continue
c     ODE case

c      Z=[y0 Y]*(DA)'
      do i=1,m
            Z(i,1)=y0(i)* DA8_1_1+Y(i,1)*DA8_1_2+Y(i,2)*DA8_1_3+
     &           Y(i,3)*DA8_1_4+Y(i,4)*DA8_1_5+Y(i,5)*DA8_1_6+
     &           Y(i,6)*DA8_1_7
            Z(i,2)=y0(i)* DA8_2_1+Y(i,1)*DA8_2_2+Y(i,2)*DA8_2_3+
     &           Y(i,3)*DA8_2_4+Y(i,4)*DA8_2_5+Y(i,5)*DA8_2_6+
     &           Y(i,6)*DA8_2_7
            Z(i,3)=y0(i)* DA8_3_1+Y(i,1)*DA8_3_2+Y(i,2)*DA8_3_3+
     &           Y(i,3)*DA8_3_4+Y(i,4)*DA8_3_5+Y(i,5)*DA8_3_6+
     &           Y(i,6)*DA8_3_7
            Z(i,4)=y0(i)* DA8_4_1+Y(i,1)*DA8_4_2+Y(i,2)*DA8_4_3+
     &           Y(i,3)*DA8_4_4+Y(i,4)*DA8_4_5+Y(i,5)*DA8_4_6+
     &           Y(i,6)*DA8_4_7
            Z(i,5)=y0(i)* DA8_5_1+Y(i,1)*DA8_5_2+Y(i,2)*DA8_5_3+
     &           Y(i,3)*DA8_5_4+Y(i,4)*DA8_5_5+Y(i,5)*DA8_5_6+
     &           Y(i,6)*DA8_5_7
            Z(i,6)=y0(i)* DA8_6_1+Y(i,1)*DA8_6_2+Y(i,2)*DA8_6_3+
     &           Y(i,3)*DA8_6_4+Y(i,4)*DA8_6_5+Y(i,5)*DA8_6_6+
     &           Y(i,6)*DA8_6_7
      end do

c     Z=Z-h*[f0 F]*(dB)'
      do i=1,m
          Z(i,1)=Z(i,1)-
     &      h*(f0(i)* DB8_1_1+F(i,1)*DB8_1_2+F(i,2)*DB8_1_3+
     &         F(i,3)*DB8_1_4+F(i,4)*DB8_1_5+F(i,5)*DB8_1_6+
     &         F(i,6)*DB8_1_7)
          Z(i,2)=Z(i,2)-
     &      h*(f0(i)* DB8_2_1+F(i,1)*DB8_2_2+F(i,2)*DB8_2_3+
     &         F(i,3)*DB8_2_4+F(i,4)*DB8_2_5+F(i,5)*DB8_2_6+
     &         F(i,6)*DB8_2_7)
          Z(i,3)=Z(i,3)-
     &      h*(f0(i)* DB8_3_1+F(i,1)*DB8_3_2+F(i,2)*DB8_3_3+
     &         F(i,3)*DB8_3_4+F(i,4)*DB8_3_5+F(i,5)*DB8_3_6+
     &         F(i,6)*DB8_3_7)
          Z(i,4)=Z(i,4)-
     &      h*(f0(i)* DB8_4_1+F(i,1)*DB8_4_2+F(i,2)*DB8_4_3+
     &         F(i,3)*DB8_4_4+F(i,4)*DB8_4_5+F(i,5)*DB8_4_6+
     &         F(i,6)*DB8_4_7)
          Z(i,5)=Z(i,5)-
     &      h*(f0(i)* DB8_5_1+F(i,1)*DB8_5_2+F(i,2)*DB8_5_3+
     &         F(i,3)*DB8_5_4+F(i,4)*DB8_5_5+F(i,5)*DB8_5_6+
     &         F(i,6)*DB8_5_7)
          Z(i,6)=Z(i,6)-
     &      h*(f0(i)* DB8_6_1+F(i,1)*DB8_6_2+F(i,2)*DB8_6_3+
     &         F(i,3)*DB8_6_4+F(i,4)*DB8_6_5+F(i,5)*DB8_6_6+
     &         F(i,6)*DB8_6_7)
      end do

c     theta*Z=Z
      call sollu(m,theta,ldlu,Z(1,1),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,2),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,3),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,4),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,5),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,6),mljac,mujac,ipvt,ijob(1))

c      Z=Z+[y0 Y]*(A2)'
      do i=1,m
            Z(i,1)=Z(i,1)+
     &        y0(i)* A28_1_1+Y(i,1)*A28_1_2+Y(i,2)*A28_1_3+
     &       Y(i,3)* A28_1_4+Y(i,4)*A28_1_5+Y(i,5)*A28_1_6+
     &       Y(i,6)* A28_1_7
            Z(i,2)=Z(i,2)+
     &        y0(i)* A28_2_1+Y(i,1)*A28_2_2+Y(i,2)*A28_2_3+
     &       Y(i,3)* A28_2_4+Y(i,4)*A28_2_5+Y(i,5)*A28_2_6+
     &       Y(i,6)* A28_2_7
            Z(i,3)=Z(i,3)+
     &        y0(i)* A28_3_1+Y(i,1)*A28_3_2+Y(i,2)*A28_3_3+
     &       Y(i,3)* A28_3_4+Y(i,4)*A28_3_5+Y(i,5)*A28_3_6+
     &       Y(i,6)* A28_3_7
            Z(i,4)=Z(i,4)+
     &        y0(i)* A28_4_1+Y(i,1)*A28_4_2+Y(i,2)*A28_4_3+
     &       Y(i,3)* A28_4_4+Y(i,4)*A28_4_5+Y(i,5)*A28_4_6+
     &       Y(i,6)* A28_4_7
            Z(i,5)=Z(i,5)+
     &        y0(i)* A28_5_1+Y(i,1)*A28_5_2+Y(i,2)*A28_5_3+
     &       Y(i,3)* A28_5_4+Y(i,4)*A28_5_5+Y(i,5)*A28_5_6+
     &       Y(i,6)* A28_5_7
            Z(i,6)=Z(i,6)+
     &        y0(i)* A28_6_1+Y(i,1)*A28_6_2+Y(i,2)*A28_6_3+
     &       Y(i,3)* A28_6_4+Y(i,4)*A28_6_5+Y(i,5)*A28_6_6+
     &       Y(i,6)* A28_6_7
      end do

c     Z=Z-h*[f0 F]*(B2)'
      do i=1,m
          Z(i,1)=Z(i,1)-h*(f0(i)*B28_1_1 + F(i,1)*gamma)
          Z(i,2)=Z(i,2)-h*(f0(i)*B28_2_1 + F(i,2)*gamma)
          Z(i,3)=Z(i,3)-h*(f0(i)*B28_3_1 + F(i,3)*gamma)
          Z(i,4)=Z(i,4)-h*(f0(i)*B28_4_1 + F(i,4)*gamma)
          Z(i,5)=Z(i,5)-h*(f0(i)*B28_5_1 + F(i,5)*gamma)
          Z(i,6)=Z(i,6)-h*(f0(i)*B28_6_1 + F(i,6)*gamma)
      end do

c     theta*Z=Z
      call sollu(m,theta,ldlu,Z(1,1),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,2),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,3),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,4),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,5),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,6),mljac,mujac,ipvt,ijob(1))

      goto 30

20    continue
c     DAE case

c      Z=[y0 Y]*(DA)'
      do i=1,m
            Z(i,1)=y0(i)* DA8_1_1+Y(i,1)*DA8_1_2+Y(i,2)*DA8_1_3+
     &           Y(i,3)*DA8_1_4+Y(i,4)*DA8_1_5+Y(i,5)*DA8_1_6+
     &           Y(i,6)*DA8_1_7
            Z(i,2)=y0(i)* DA8_2_1+Y(i,1)*DA8_2_2+Y(i,2)*DA8_2_3+
     &           Y(i,3)*DA8_2_4+Y(i,4)*DA8_2_5+Y(i,5)*DA8_2_6+
     &           Y(i,6)*DA8_2_7
            Z(i,3)=y0(i)* DA8_3_1+Y(i,1)*DA8_3_2+Y(i,2)*DA8_3_3+
     &           Y(i,3)*DA8_3_4+Y(i,4)*DA8_3_5+Y(i,5)*DA8_3_6+
     &           Y(i,6)*DA8_3_7
            Z(i,4)=y0(i)* DA8_4_1+Y(i,1)*DA8_4_2+Y(i,2)*DA8_4_3+
     &           Y(i,3)*DA8_4_4+Y(i,4)*DA8_4_5+Y(i,5)*DA8_4_6+
     &           Y(i,6)*DA8_4_7
            Z(i,5)=y0(i)* DA8_5_1+Y(i,1)*DA8_5_2+Y(i,2)*DA8_5_3+
     &           Y(i,3)*DA8_5_4+Y(i,4)*DA8_5_5+Y(i,5)*DA8_5_6+
     &           Y(i,6)*DA8_5_7
            Z(i,6)=y0(i)* DA8_6_1+Y(i,1)*DA8_6_2+Y(i,2)*DA8_6_3+
     &           Y(i,3)*DA8_6_4+Y(i,4)*DA8_6_5+Y(i,5)*DA8_6_6+
     &           Y(i,6)*DA8_6_7
      end do

c     MZ_i=M0*Z_i, i=1,2,...,6

      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,1),MZ(1,1),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,2),MZ(1,2),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,3),MZ(1,3),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,4),MZ(1,4),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,5),MZ(1,5),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,6),MZ(1,6),ijob(2))

c     MZ=MZ-h*[f0 F]*(dB)'
      do i=1,m
          MZ(i,1)=MZ(i,1)-
     &      h*(f0(i)* DB8_1_1+F(i,1)*DB8_1_2+F(i,2)*DB8_1_3+
     &         F(i,3)*DB8_1_4+F(i,4)*DB8_1_5+F(i,5)*DB8_1_6+
     &         F(i,6)*DB8_1_7)
          MZ(i,2)=MZ(i,2)-
     &      h*(f0(i)* DB8_2_1+F(i,1)*DB8_2_2+F(i,2)*DB8_2_3+
     &         F(i,3)*DB8_2_4+F(i,4)*DB8_2_5+F(i,5)*DB8_2_6+
     &         F(i,6)*DB8_2_7)
          MZ(i,3)=MZ(i,3)-
     &      h*(f0(i)* DB8_3_1+F(i,1)*DB8_3_2+F(i,2)*DB8_3_3+
     &         F(i,3)*DB8_3_4+F(i,4)*DB8_3_5+F(i,5)*DB8_3_6+
     &         F(i,6)*DB8_3_7)
          MZ(i,4)=MZ(i,4)-
     &      h*(f0(i)* DB8_4_1+F(i,1)*DB8_4_2+F(i,2)*DB8_4_3+
     &         F(i,3)*DB8_4_4+F(i,4)*DB8_4_5+F(i,5)*DB8_4_6+
     &         F(i,6)*DB8_4_7)
          MZ(i,5)=MZ(i,5)-
     &      h*(f0(i)* DB8_5_1+F(i,1)*DB8_5_2+F(i,2)*DB8_5_3+
     &         F(i,3)*DB8_5_4+F(i,4)*DB8_5_5+F(i,5)*DB8_5_6+
     &         F(i,6)*DB8_5_7)
          MZ(i,6)=MZ(i,6)-
     &      h*(f0(i)* DB8_6_1+F(i,1)*DB8_6_2+F(i,2)*DB8_6_3+
     &         F(i,3)*DB8_6_4+F(i,4)*DB8_6_5+F(i,5)*DB8_6_6+
     &         F(i,6)*DB8_6_7)
      end do

c     theta*MZ=MZ
      call sollu(m,theta,ldlu,MZ(1,1),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,2),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,3),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,4),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,5),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,6),mljac,mujac,ipvt,ijob(1))

c      MZ=MZ+[y0 Y]*(A2)'
      do i=1,m
            MZ(i,1)=MZ(i,1)+
     &        y0(i)* A28_1_1+Y(i,1)*A28_1_2+Y(i,2)*A28_1_3+
     &       Y(i,3)* A28_1_4+Y(i,4)*A28_1_5+Y(i,5)*A28_1_6+
     &       Y(i,6)* A28_1_7
            MZ(i,2)=MZ(i,2)+
     &        y0(i)* A28_2_1+Y(i,1)*A28_2_2+Y(i,2)*A28_2_3+
     &       Y(i,3)* A28_2_4+Y(i,4)*A28_2_5+Y(i,5)*A28_2_6+
     &       Y(i,6)* A28_2_7
            MZ(i,3)=MZ(i,3)+
     &        y0(i)* A28_3_1+Y(i,1)*A28_3_2+Y(i,2)*A28_3_3+
     &       Y(i,3)* A28_3_4+Y(i,4)*A28_3_5+Y(i,5)*A28_3_6+
     &       Y(i,6)* A28_3_7
            MZ(i,4)=MZ(i,4)+
     &        y0(i)* A28_4_1+Y(i,1)*A28_4_2+Y(i,2)*A28_4_3+
     &       Y(i,3)* A28_4_4+Y(i,4)*A28_4_5+Y(i,5)*A28_4_6+
     &       Y(i,6)* A28_4_7
            MZ(i,5)=MZ(i,5)+
     &        y0(i)* A28_5_1+Y(i,1)*A28_5_2+Y(i,2)*A28_5_3+
     &       Y(i,3)* A28_5_4+Y(i,4)*A28_5_5+Y(i,5)*A28_5_6+
     &       Y(i,6)* A28_5_7
            MZ(i,6)=MZ(i,6)+
     &        y0(i)* A28_6_1+Y(i,1)*A28_6_2+Y(i,2)*A28_6_3+
     &       Y(i,3)* A28_6_4+Y(i,4)*A28_6_5+Y(i,5)*A28_6_6+
     &       Y(i,6)* A28_6_7
      end do

c     Z_i=M0*MZ_i, i=1,2,...,6

      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,1),Z(1,1),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,2),Z(1,2),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,3),Z(1,3),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,4),Z(1,4),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,5),Z(1,5),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,6),Z(1,6),ijob(2))

c     Z=Z-h*[f0 F]*(B2)'
      do i=1,m
          Z(i,1)=Z(i,1)-h*(f0(i)*B28_1_1 + F(i,1)*gamma)
          Z(i,2)=Z(i,2)-h*(f0(i)*B28_2_1 + F(i,2)*gamma)
          Z(i,3)=Z(i,3)-h*(f0(i)*B28_3_1 + F(i,3)*gamma)
          Z(i,4)=Z(i,4)-h*(f0(i)*B28_4_1 + F(i,4)*gamma)
          Z(i,5)=Z(i,5)-h*(f0(i)*B28_5_1 + F(i,5)*gamma)
          Z(i,6)=Z(i,6)-h*(f0(i)*B28_6_1 + F(i,6)*gamma)
      end do

c     theta*Z=Z
      call sollu(m,theta,ldlu,Z(1,1),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,2),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,3),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,4),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,5),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,6),mljac,mujac,ipvt,ijob(1))

30    continue

c
c     Y = Y-Z
c
      do j=1,k
         do i=1,m
           Y(i,j) = Y(i,j) - Z(i,j)
         end do
      end do

      return
      end

C -------------------------------------------------------------------

      subroutine blendstep10(m,y0,f0,Y,F,h,theta,ipvt,Z,gamma,
     &                             ldlu,mljac,mujac,ijob,imas,
     &                                ldmas,mlmas,mumas,M0,MZ)
c
c     Blended iteration for the 10th order method
c

      implicit none

      integer k
      parameter (k=8)

c     Input parameters
      integer m,ipvt(m),ldlu,mljac,mujac,ijob(2),imas,ldmas,mlmas,mumas
      double precision  y0(m),f0(m),theta(ldlu,m),gamma,h,M0(ldmas,m)

c     I/O parameters
      double precision Y(m,k),F(m,k),MZ(m,k)

c     Output parameters
      double precision Z(m,k)

c     Local variables
      integer i,j

C ---------------------------------------------------------------------------------------
C 10th order BIM
C ---------------------------------------------------------------------------------------

      double precision
     &   DA10_1_1,DA10_1_2,DA10_1_3,DA10_1_4,DA10_1_5,DA10_1_6,DA10_1_7,
     &   DA10_1_8,DA10_1_9,
     &   DA10_2_1,DA10_2_2,DA10_2_3,DA10_2_4,DA10_2_5,DA10_2_6,DA10_2_7,
     &   DA10_2_8,DA10_2_9,
     &   DA10_3_1,DA10_3_2,DA10_3_3,DA10_3_4,DA10_3_5,DA10_3_6,DA10_3_7,
     &   DA10_3_8,DA10_3_9,
     &   DA10_4_1,DA10_4_2,DA10_4_3,DA10_4_4,DA10_4_5,DA10_4_6,DA10_4_7,
     &   DA10_4_8,DA10_4_9,
     &   DA10_5_1,DA10_5_2,DA10_5_3,DA10_5_4,DA10_5_5,DA10_5_6,DA10_5_7,
     &   DA10_5_8,DA10_5_9,
     &   DA10_6_1,DA10_6_2,DA10_6_3,DA10_6_4,DA10_6_5,DA10_6_6,DA10_6_7,
     &   DA10_6_8,DA10_6_9,
     &   DA10_7_1,DA10_7_2,DA10_7_3,DA10_7_4,DA10_7_5,DA10_7_6,DA10_7_7,
     &   DA10_7_8,DA10_7_9,
     &   DA10_8_1,DA10_8_2,DA10_8_3,DA10_8_4,DA10_8_5,DA10_8_6,DA10_8_7,
     &   DA10_8_8,DA10_8_9,
     &   DB10_1_1,DB10_1_2,DB10_1_3,DB10_1_4,DB10_1_5,DB10_1_6,DB10_1_7,
     &   DB10_1_8,DB10_1_9,
     &   DB10_2_1,DB10_2_2,DB10_2_3,DB10_2_4,DB10_2_5,DB10_2_6,DB10_2_7,
     &   DB10_2_8,DB10_2_9,
     &   DB10_3_1,DB10_3_2,DB10_3_3,DB10_3_4,DB10_3_5,DB10_3_6,DB10_3_7,
     &   DB10_3_8,DB10_3_9,
     &   DB10_4_1,DB10_4_2,DB10_4_3,DB10_4_4,DB10_4_5,DB10_4_6,DB10_4_7,
     &   DB10_4_8,DB10_4_9,
     &   DB10_5_1,DB10_5_2,DB10_5_3,DB10_5_4,DB10_5_5,DB10_5_6,DB10_5_7,
     &   DB10_5_8,DB10_5_9,
     &   DB10_6_1,DB10_6_2,DB10_6_3,DB10_6_4,DB10_6_5,DB10_6_6,DB10_6_7,
     &   DB10_6_8,DB10_6_9,
     &   DB10_7_1,DB10_7_2,DB10_7_3,DB10_7_4,DB10_7_5,DB10_7_6,DB10_7_7,
     &   DB10_7_8,DB10_7_9,
     &   DB10_8_1,DB10_8_2,DB10_8_3,DB10_8_4,DB10_8_5,DB10_8_6,DB10_8_7,
     &   DB10_8_8,DB10_8_9,
     &   A210_1_1,A210_1_2,A210_1_3,A210_1_4,A210_1_5,A210_1_6,A210_1_7,
     &   A210_1_8,A210_1_9,
     &   A210_2_1,A210_2_2,A210_2_3,A210_2_4,A210_2_5,A210_2_6,A210_2_7,
     &   A210_2_8,A210_2_9,
     &   A210_3_1,A210_3_2,A210_3_3,A210_3_4,A210_3_5,A210_3_6,A210_3_7,
     &   A210_3_8,A210_3_9,
     &   A210_4_1,A210_4_2,A210_4_3,A210_4_4,A210_4_5,A210_4_6,A210_4_7,
     &   A210_4_8,A210_4_9,
     &   A210_5_1,A210_5_2,A210_5_3,A210_5_4,A210_5_5,A210_5_6,A210_5_7,
     &   A210_5_8,A210_5_9,
     &   A210_6_1,A210_6_2,A210_6_3,A210_6_4,A210_6_5,A210_6_6,A210_6_7,
     &   A210_6_8,A210_6_9,
     &   A210_7_1,A210_7_2,A210_7_3,A210_7_4,A210_7_5,A210_7_6,A210_7_7,
     &   A210_7_8,A210_7_9,
     &   A210_8_1,A210_8_2,A210_8_3,A210_8_4,A210_8_5,A210_8_6,A210_8_7,
     &   A210_8_8,A210_8_9,
     &   B210_1_1,B210_2_1,B210_3_1,B210_4_1,B210_5_1,B210_6_1,B210_7_1,
     &   B210_8_1

      parameter(
     &  DA10_1_1    =   -1116185640857D0/134217728D4,
     &  DA10_1_2    =   268195813141D0/14680064D4,   
     &  DA10_1_3    =   -6468370013D0/33554432D2,
     &  DA10_1_4    =   +22440146897D0/12582912D3,
     &  DA10_1_5    =   -57418737523D0/402653184D2,
     &  DA10_1_6    =   +17489295313D0/2097152D4,
     &  DA10_1_7    =   -82495624981D0/25165824D4,
     &  DA10_1_8    =   +193947079D0/25165824D2,
     &  DA10_1_9    =   -15367501777D0/1879048192D3,
     &  DA10_2_1    =   -1408339461D0/1048576D3,   
     &  DA10_2_2    =   +133859921D0/114688D3,
     &  DA10_2_3    =   -21140373D0/32768D4,
     &  DA10_2_4    =   +15150619D0/16384D3,
     &  DA10_2_5    =   -67553873D0/524288D2,
     &  DA10_2_6    =   +44952727D0/49152D3,
     &  DA10_2_7    =   -26201627D0/65536D3,
     &  DA10_2_8    =   +8246437D0/8192D4,
     &  DA10_2_9    =   -246865651D0/22020096D3,
     &  DA10_3_1    =   -552495440623D0/805306368D3,
     &  DA10_3_2    =   -28197945999D0/29360128D3,
     &  DA10_3_3    =   +32441533071D0/16777216D3,
     &  DA10_3_4    =   -51906383051D0/6291456D4,
     &  DA10_3_5    =   +15467184783D0/134217728D2,
     &  DA10_3_6    =   -3942154371D0/4194304D3,
     &  DA10_3_7    =   +21125300879D0/50331648D3,
     &  DA10_3_8    =   -445562559D0/4194304D3,
     &  DA10_3_9    =   +111284620491D0/939524096D4,
     &  DA10_4_1    =   -19796057D0/2048D4,
     &  DA10_4_2    =   -535553D0/672D4,
     &  DA10_4_3    =   +63403D0/128D4,
     &  DA10_4_4    =   +281941D0/96D4,
     &  DA10_4_5    =   +50403D0/4096D1,  
     &  DA10_4_6    =   -219887D0/32D4,
     &  DA10_4_7    =   +754091D0/384D4,
     &  DA10_4_8    =   -39121D0/96D4,
     &  DA10_4_9    =   +581419D0/14336D4,
     &  DA10_5_1    =   -64468532577D0/536870912D2,  
     &  DA10_5_2    =   +16856063921D0/29360128D3,
     &  DA10_5_3    =   -47031869203D0/50331648D3,
     &  DA10_5_4    =   +4439913787D0/4194304D3,
     &  DA10_5_5    =   -6247096241D0/134217728D2,
     &  DA10_5_6    =   +1286778713D0/8388608D2,   
     &  DA10_5_7    =   -11512558907D0/16777216D3,
     &  DA10_5_8    =   +560510849D0/4194304D3,
     &  DA10_5_9    =   -75322449683D0/5637144576D3,
     &  DA10_6_1    =   -14564721451D0/1572864D4,
     &  DA10_6_2    =   -109580619D0/57344D4,
     &  DA10_6_3    =   +15285519D0/65536D3,
     &  DA10_6_4    =   -4234511D0/49152D3,
     &  DA10_6_5    =   -3573501D0/1048576D1,
     &  DA10_6_6    =   +84173553D0/8192D4,
     &  DA10_6_7    =   484815179D0/98304D4,
     &  DA10_6_8    =   -3694911D0/16384D3,
     &  DA10_6_9    =   +114744591D0/7340032D3,
     &  DA10_7_1    =   -10378315789D0/1073741824D1,
     &  DA10_7_2    =   -2207113739D0/12582912D3,
     &  DA10_7_3    =   +48853634347D0/8388608D4,
     &  DA10_7_4    =   -16701919087D0/12582912D3,
     &  DA10_7_5    =   +84761718349D0/402653184D2,
     &  DA10_7_6    =   -2054315207D0/8388608D2,
     &  DA10_7_7    =   +120669802351D0/50331648D3,
     &  DA10_7_8    =   -1748617243D0/2097152D4,
     &  DA10_7_9    =   -22368269479D0/268435456D3,
     &  DA10_8_1    =   -17349D0/16D3,   
     &  DA10_8_2    =   +1349D0/175D1,
     &  DA10_8_3    =   -9443D0/3D3,
     &  DA10_8_4    =   +9443D0/125D1,
     &  DA10_8_5    =   -9443D0/8D2,
     &  DA10_8_6    =   +9443D0/75D1,
     &  DA10_8_7    =   -9443D0/1D3,
     &  DA10_8_8    =   +1349D0/25D1,
     &  DA10_8_9    =   -466589D0/56D4)





      parameter(
     &  A210_1_1    =   -225991639143D0/134217728D4,
     &  A210_1_2    =   -121395173141D0/14680064D4,
     &  A210_1_3    =   +6468370013D0/33554432D2,
     &  A210_1_4    =   -22440146897D0/12582912D3,
     &  A210_1_5    =   +57418737523D0/402653184D2,
     &  A210_1_6    =   -17489295313D0/2097152D4,
     &  A210_1_7    =   +82495624981D0/25165824D4,
     &  A210_1_8    =   -193947079D0/25165824D2,
     &  A210_1_9    =   +15367501777D0/1879048192D3,
     &  A210_2_1    =   +359763461D0/1048576D3,
     &  A210_2_2    =   -133859921D0/114688D3,
     &  A210_2_3    =   +348820373D0/32768D4,
     &  A210_2_4    =   -15150619D0/16384D3,
     &  A210_2_5    =   +67553873D0/524288D2,
     &  A210_2_6    =   -44952727D0/49152D3,
     &  A210_2_7    =   +26201627D0/65536D3,
     &  A210_2_8    =   -8246437D0/8192D4,
     &  A210_2_9    =   +246865651D0/22020096D3,
     &  A210_3_1    =   -252810927377D0/805306368D3,
     &  A210_3_2    =   +28197945999D0/29360128D3,
     &  A210_3_3    =   -32441533071D0/16777216D3,
     &  A210_3_4    =   +114820943051D0/6291456D4,
     &  A210_3_5    =   -15467184783D0/134217728D2,
     &  A210_3_6    =   +3942154371D0/4194304D3,
     &  A210_3_7    =   -21125300879D0/50331648D3,
     &  A210_3_8    =   +445562559D0/4194304D3,
     &  A210_3_9    =   -111284620491D0/939524096D4,
     &  A210_4_1    =   -683943D0/2048D4,
     &  A210_4_2    =   +535553D0/672D4,
     &  A210_4_3    =   -63403D0/128D4,
     &  A210_4_4    =   -281941D0/96D4,
     &  A210_4_5    =   -9443D0/4096D1,
     &  A210_4_6    =   +219887D0/32D4,
     &  A210_4_7    =   -754091D0/384D4,
     &  A210_4_8    =   +39121D0/96D4,
     &  A210_4_9    =   -581419D0/14336D4,
     &  A210_5_1    =   +10781441377D0/536870912D2,
     &  A210_5_2    =   -16856063921D0/29360128D3,
     &  A210_5_3    =   +47031869203D0/50331648D3,
     &  A210_5_4    =   -4439913787D0/4194304D3,
     &  A210_5_5    =   +6247096241D0/134217728D2,
     &  A210_5_6    =   -447917913D0/8388608D2,
     &  A210_5_7    =   +11512558907D0/16777216D3,
     &  A210_5_8    =   -560510849D0/4194304D3,
     &  A210_5_9    =   +75322449683D0/5637144576D3,
     &  A210_6_1    =   -1163918549D0/1572864D4,
     &  A210_6_2    =   +109580619D0/57344D4,
     &  A210_6_3    =   -15285519D0/65536D3,
     &  A210_6_4    =   +4234511D0/49152D3,
     &  A210_6_5    =   +3573501D0/1048576D1,
     &  A210_6_6    =   -84173553D0/8192D4,
     &  A210_6_7    =   +498224821D0/98304D4,
     &  A210_6_8    =   +3694911D0/16384D3,
     &  A210_6_9    =   -114744591D0/7340032D3,
     &  A210_7_1    =   -359102451D0/1073741824D1,
     &  A210_7_2    =   +2207113739D0/12582912D3,
     &  A210_7_3    =   -48853634347D0/8388608D4,
     &  A210_7_4    =   +16701919087D0/12582912D3,
     &  A210_7_5    =   -84761718349D0/402653184D2,
     &  A210_7_6    =   +2054315207D0/8388608D2,
     &  A210_7_7    =   -120669802351D0/50331648D3,
     &  A210_7_8    =   +22720137243D0/2097152D4,
     &  A210_7_9    =   +22368269479D0/268435456D3,
     &  A210_8_1    =   +1349D0/16D3,
     &  A210_8_2    =   -1349D0/175D1,
     &  A210_8_3    =   +9443D0/3D3,
     &  A210_8_4    =   -9443D0/125D1,
     &  A210_8_5    =   +9443D0/8D2,
     &  A210_8_6    =   -9443D0/75D1,
     &  A210_8_7    =   +9443D0/1D3,
     &  A210_8_8    =   -1349D0/25D1,
     &  A210_8_9    =   +1026589D0/56D4)



      parameter(
     &  DB10_1_1    =   +761073292818799D0/2720626900992D3,
     &  DB10_1_2    =   +557937353D0/1297296D3,
     &  DB10_1_3    =   -214412651D0/2594592D2,
     &  DB10_1_4    =   +170436457D0/2594592D2,
     &  DB10_1_5    =   -7282853D0/2594592D1,
     &  DB10_1_6    =   -4668473D0/2594592D2,
     &  DB10_1_7    =   +22494019D0/2594592D2,
     &  DB10_1_8    =   -10434029D0/2594592D2,
     &  DB10_1_9    =   +3345851D0/5189184D2,
     &  DB10_2_1    =   +1183625282033D0/297568567296D1,
     &  DB10_2_2    =   +23930143D0/14189175D0,
     &  DB10_2_3    =   -213517489D0/162162D3,
     &  DB10_2_4    =   +3230893D0/2027025D0,
     &  DB10_2_5    =   -686209D0/405405D0,
     &  DB10_2_6    =   +2420083D0/2027025D0,
     &  DB10_2_7    =   -4396151D0/81081D2,
     &  DB10_2_8    =   +2038273D0/14189175D0,
     &  DB10_2_9    =   -1918153D0/1135134D2,
     &  DB10_3_1    =   +70530807725549D0/4232086290432D2,
     &  DB10_3_2    =   +330413059D0/2018016D2,
     &  DB10_3_3    =   -3305147D0/288288D2,
     &  DB10_3_4    =   +30733931D0/20592D3,
     &  DB10_3_5    =   -4993271D0/288288D1,
     &  DB10_3_6    =   +4741897D0/41184D2,
     &  DB10_3_7    =   -14566397D0/288288D2,
     &  DB10_3_8    =   +26359309D0/2018016D2,
     &  DB10_3_9    =   -6083071D0/4036032D2,
     &  DB10_4_1    =   +3304437553D0/1162377216D1,
     &  DB10_4_2    =   +21369776D0/14189175D0,
     &  DB10_4_3    =   +92756D0/289575D0,
     &  DB10_4_4    =   +3454736D0/2027025D0,
     &  DB10_4_5    =   -90373469D0/162162D3,
     &  DB10_4_6    =   +91376D0/2027025D0,
     &  DB10_4_7    =   +6956D0/289575D0,
     &  DB10_4_8    =   -251824D0/14189175D0,
     &  DB10_4_9    =   +46493D0/14189175D0,
     &  DB10_5_1    =   +6901687296484073D0/19044388306944D3,
     &  DB10_5_2    =   +115012739D0/72648576D0,
     &  DB10_5_3    =   +628493D0/10378368D0,
     &  DB10_5_4    =   +22584689D0/10378368D0,
     &  DB10_5_5    =   -4595D0/741312D0,
     &  DB10_5_6    =   +630493723D0/1297296D3,
     &  DB10_5_7    =   -3425557D0/10378368D0,
     &  DB10_5_8    =   +5553389D0/72648576D0,
     &  DB10_5_9    =   -1181891D0/145297152D0,
     &  DB10_6_1    =   +160618157669D0/635830272D3,
     &  DB10_6_2    =   +198307D0/121275D0,
     &  DB10_6_3    =   -9029D0/693D2,
     &  DB10_6_4    =   +44857D0/17325D0,
     &  DB10_6_5    =   -2131D0/3465D0,
     &  DB10_6_6    =   +37927D0/17325D0,
     &  DB10_6_7    =   -976837D0/1386D3,
     &  DB10_6_8    =   +11197D0/121275D0,
     &  DB10_6_9    =   -11197D0/9702D2,
     &  DB10_7_1    =   +153681594740771D0/5441253801984D2,
     &  DB10_7_2    =   +404839597D0/2594592D2,
     &  DB10_7_3    =   +5425339D0/370656D2,
     &  DB10_7_4    =   +73758847D0/370656D2,
     &  DB10_7_5    =   +878107D0/370656D1,
     &  DB10_7_6    =   +48743857D0/370656D2,
     &  DB10_7_7    =   +39269149D0/370656D2,
     &  DB10_7_8    =   -335961817D0/1297296D3,
     &  DB10_7_9    =   -7219753D0/5189184D2,
     &  DB10_8_1    =   +3956D0/14175D0,
     &  DB10_8_2    =   +23552D0/14175D0,
     &  DB10_8_3    =   -3712D0/14175D0,
     &  DB10_8_4    =   +41984D0/14175D0,
     &  DB10_8_5    =   -3632D0/2835D0,
     &  DB10_8_6    =   +41984D0/14175D0,
     &  DB10_8_7    =   -3712D0/14175D0,
     &  DB10_8_8    =   +23552D0/14175D0,
     &  DB10_8_9    =   -448403D0/1134D3,
     &  B210_1_1    =   +1037851801D0/33554432D3,
     &  B210_2_1    =   -3193083D0/262144D2,
     &  B210_3_1    =   +153048097D0/134217728D1,
     &  B210_4_1    =   +1349D0/1024D2,
     &  B210_5_1    =   -2509046919D0/33554432D3,
     &  B210_6_1    =   +3762361D0/131072D3,
     &  B210_7_1    =   +52838981D0/67108864D2,
     &  B210_8_1    =   0D0)


C ---------------------------------------------------------------------------------------

      goto(10,20) imas+1

10    continue
c     ODE case

c      Z=[Y0 Y]*(DA)'
      do i=1,m
            Z(i,1)=y0(i)* DA10_1_1+Y(i,1)*DA10_1_2+Y(i,2)*DA10_1_3+
     &           Y(i,3)*DA10_1_4+Y(i,4)*DA10_1_5+Y(i,5)*DA10_1_6+
     &           Y(i,6)*DA10_1_7+Y(i,7)*DA10_1_8+Y(i,8)*DA10_1_9
            Z(i,2)=y0(i)* DA10_2_1+Y(i,1)*DA10_2_2+Y(i,2)*DA10_2_3+
     &           Y(i,3)*DA10_2_4+Y(i,4)*DA10_2_5+Y(i,5)*DA10_2_6+
     &           Y(i,6)*DA10_2_7+Y(i,7)*DA10_2_8+Y(i,8)*DA10_2_9
            Z(i,3)=y0(i)* DA10_3_1+Y(i,1)*DA10_3_2+Y(i,2)*DA10_3_3+
     &           Y(i,3)*DA10_3_4+Y(i,4)*DA10_3_5+Y(i,5)*DA10_3_6+
     &           Y(i,6)*DA10_3_7+Y(i,7)*DA10_3_8+Y(i,8)*DA10_3_9
            Z(i,4)=y0(i)* DA10_4_1+Y(i,1)*DA10_4_2+Y(i,2)*DA10_4_3+
     &           Y(i,3)*DA10_4_4+Y(i,4)*DA10_4_5+Y(i,5)*DA10_4_6+
     &           Y(i,6)*DA10_4_7+Y(i,7)*DA10_4_8+Y(i,8)*DA10_4_9
            Z(i,5)=y0(i)* DA10_5_1+Y(i,1)*DA10_5_2+Y(i,2)*DA10_5_3+
     &           Y(i,3)*DA10_5_4+Y(i,4)*DA10_5_5+Y(i,5)*DA10_5_6+
     &           Y(i,6)*DA10_5_7+Y(i,7)*DA10_5_8+Y(i,8)*DA10_5_9
            Z(i,6)=y0(i)* DA10_6_1+Y(i,1)*DA10_6_2+Y(i,2)*DA10_6_3+
     &           Y(i,3)*DA10_6_4+Y(i,4)*DA10_6_5+Y(i,5)*DA10_6_6+
     &           Y(i,6)*DA10_6_7+Y(i,7)*DA10_6_8+Y(i,8)*DA10_6_9
            Z(i,7)=y0(i)* DA10_7_1+Y(i,1)*DA10_7_2+Y(i,2)*DA10_7_3+
     &           Y(i,3)*DA10_7_4+Y(i,4)*DA10_7_5+Y(i,5)*DA10_7_6+
     &           Y(i,6)*DA10_7_7+Y(i,7)*DA10_7_8+Y(i,8)*DA10_7_9
            Z(i,8)=y0(i)* DA10_8_1+Y(i,1)*DA10_8_2+Y(i,2)*DA10_8_3+
     &           Y(i,3)*DA10_8_4+Y(i,4)*DA10_8_5+Y(i,5)*DA10_8_6+
     &           Y(i,6)*DA10_8_7+Y(i,7)*DA10_8_8+Y(i,8)*DA10_8_9
      end do

c     Z=Z-h*[f0 F]*(dB)'
      do i=1,m
          Z(i,1) = Z(i,1)-
     &        h*(f0(i)* DB10_1_1+F(i,1)*DB10_1_2+F(i,2)*DB10_1_3+
     &           F(i,3)*DB10_1_4+F(i,4)*DB10_1_5+F(i,5)*DB10_1_6+
     &           F(i,6)*DB10_1_7+F(i,7)*DB10_1_8+F(i,8)*DB10_1_9)
          Z(i,2) = Z(i,2)-
     &        h*(f0(i)* DB10_2_1+F(i,1)*DB10_2_2+F(i,2)*DB10_2_3+
     &           F(i,3)*DB10_2_4+F(i,4)*DB10_2_5+F(i,5)*DB10_2_6+
     &           F(i,6)*DB10_2_7+F(i,7)*DB10_2_8+F(i,8)*DB10_2_9)
          Z(i,3) = Z(i,3)-
     &        h*(f0(i)* DB10_3_1+F(i,1)*DB10_3_2+F(i,2)*DB10_3_3+
     &           F(i,3)*DB10_3_4+F(i,4)*DB10_3_5+F(i,5)*DB10_3_6+
     &           F(i,6)*DB10_3_7+F(i,7)*DB10_3_8+F(i,8)*DB10_3_9)
          Z(i,4) = Z(i,4)-
     &        h*(f0(i)* DB10_4_1+F(i,1)*DB10_4_2+F(i,2)*DB10_4_3+
     &           F(i,3)*DB10_4_4+F(i,4)*DB10_4_5+F(i,5)*DB10_4_6+
     &           F(i,6)*DB10_4_7+F(i,7)*DB10_4_8+F(i,8)*DB10_4_9)
          Z(i,5) = Z(i,5)-
     &        h*(f0(i)* DB10_5_1+F(i,1)*DB10_5_2+F(i,2)*DB10_5_3+
     &           F(i,3)*DB10_5_4+F(i,4)*DB10_5_5+F(i,5)*DB10_5_6+
     &           F(i,6)*DB10_5_7+F(i,7)*DB10_5_8+F(i,8)*DB10_5_9)
          Z(i,6) = Z(i,6)-
     &        h*(f0(i)* DB10_6_1+F(i,1)*DB10_6_2+F(i,2)*DB10_6_3+
     &           F(i,3)*DB10_6_4+F(i,4)*DB10_6_5+F(i,5)*DB10_6_6+
     &           F(i,6)*DB10_6_7+F(i,7)*DB10_6_8+F(i,8)*DB10_6_9)
          Z(i,7) = Z(i,7)-
     &        h*(f0(i)* DB10_7_1+F(i,1)*DB10_7_2+F(i,2)*DB10_7_3+
     &           F(i,3)*DB10_7_4+F(i,4)*DB10_7_5+F(i,5)*DB10_7_6+
     &           F(i,6)*DB10_7_7+F(i,7)*DB10_7_8+F(i,8)*DB10_7_9)
          Z(i,8) = Z(i,8)-
     &        h*(f0(i)* DB10_8_1+F(i,1)*DB10_8_2+F(i,2)*DB10_8_3+
     &           F(i,3)*DB10_8_4+F(i,4)*DB10_8_5+F(i,5)*DB10_8_6+
     &           F(i,6)*DB10_8_7+F(i,7)*DB10_8_8+F(i,8)*DB10_8_9)
      end do

c     theta*Z=Z
      call sollu(m,theta,ldlu,Z(1,1),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,2),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,3),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,4),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,5),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,6),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,7),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,8),mljac,mujac,ipvt,ijob(1))

c      Z=Z+[y0 Y]*(A2)'
      do i=1,m
            Z(i,1)=Z(i,1)+
     &               y0(i)* A210_1_1+Y(i,1)*A210_1_2+Y(i,2)*A210_1_3+
     &              Y(i,3)* A210_1_4+Y(i,4)*A210_1_5+Y(i,5)*A210_1_6+
     &              Y(i,6)* A210_1_7+Y(i,7)*A210_1_8+Y(i,8)*A210_1_9
            Z(i,2)=Z(i,2)+
     &               y0(i)* A210_2_1+Y(i,1)*A210_2_2+Y(i,2)*A210_2_3+
     &              Y(i,3)* A210_2_4+Y(i,4)*A210_2_5+Y(i,5)*A210_2_6+
     &              Y(i,6)* A210_2_7+Y(i,7)*A210_2_8+Y(i,8)*A210_2_9
            Z(i,3)=Z(i,3)+
     &               y0(i)* A210_3_1+Y(i,1)*A210_3_2+Y(i,2)*A210_3_3+
     &              Y(i,3)* A210_3_4+Y(i,4)*A210_3_5+Y(i,5)*A210_3_6+
     &              Y(i,6)* A210_3_7+Y(i,7)*A210_3_8+Y(i,8)*A210_3_9
            Z(i,4)=Z(i,4)+
     &               y0(i)* A210_4_1+Y(i,1)*A210_4_2+Y(i,2)*A210_4_3+
     &              Y(i,3)* A210_4_4+Y(i,4)*A210_4_5+Y(i,5)*A210_4_6+
     &              Y(i,6)* A210_4_7+Y(i,7)*A210_4_8+Y(i,8)*A210_4_9
            Z(i,5)=Z(i,5)+
     &               y0(i)* A210_5_1+Y(i,1)*A210_5_2+Y(i,2)*A210_5_3+
     &              Y(i,3)* A210_5_4+Y(i,4)*A210_5_5+Y(i,5)*A210_5_6+
     &              Y(i,6)* A210_5_7+Y(i,7)*A210_5_8+Y(i,8)*A210_5_9
            Z(i,6)=Z(i,6)+
     &               y0(i)* A210_6_1+Y(i,1)*A210_6_2+Y(i,2)*A210_6_3+
     &              Y(i,3)* A210_6_4+Y(i,4)*A210_6_5+Y(i,5)*A210_6_6+
     &              Y(i,6)* A210_6_7+Y(i,7)*A210_6_8+Y(i,8)*A210_6_9
            Z(i,7)=Z(i,7)+
     &               y0(i)* A210_7_1+Y(i,1)*A210_7_2+Y(i,2)*A210_7_3+
     &              Y(i,3)* A210_7_4+Y(i,4)*A210_7_5+Y(i,5)*A210_7_6+
     &              Y(i,6)* A210_7_7+Y(i,7)*A210_7_8+Y(i,8)*A210_7_9
            Z(i,8)=Z(i,8)+
     &               y0(i)* A210_8_1+Y(i,1)*A210_8_2+Y(i,2)*A210_8_3+
     &              Y(i,3)* A210_8_4+Y(i,4)*A210_8_5+Y(i,5)*A210_8_6+
     &              Y(i,6)* A210_8_7+Y(i,7)*A210_8_8+Y(i,8)*A210_8_9
      end do

c     Z=Z-h*[f0 F]*(B2)'
      do i=1,m
          Z(i,1)=Z(i,1)-h*(f0(i)*B210_1_1 + F(i,1)*gamma)
          Z(i,2)=Z(i,2)-h*(f0(i)*B210_2_1 + F(i,2)*gamma)
          Z(i,3)=Z(i,3)-h*(f0(i)*B210_3_1 + F(i,3)*gamma)
          Z(i,4)=Z(i,4)-h*(f0(i)*B210_4_1 + F(i,4)*gamma)
          Z(i,5)=Z(i,5)-h*(f0(i)*B210_5_1 + F(i,5)*gamma)
          Z(i,6)=Z(i,6)-h*(f0(i)*B210_6_1 + F(i,6)*gamma)
          Z(i,7)=Z(i,7)-h*(f0(i)*B210_7_1 + F(i,7)*gamma)
          Z(i,8)=Z(i,8)-h*(f0(i)*B210_8_1 + F(i,8)*gamma)
      end do

c     theta*Z=Z
      call sollu(m,theta,ldlu,Z(1,1),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,2),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,3),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,4),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,5),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,6),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,7),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,8),mljac,mujac,ipvt,ijob(1))

      goto 30

20    continue
c     DAE case

c      Z=[Y0 Y]*(DA)'
      do i=1,m
            Z(i,1)=y0(i)* DA10_1_1+Y(i,1)*DA10_1_2+Y(i,2)*DA10_1_3+
     &           Y(i,3)*DA10_1_4+Y(i,4)*DA10_1_5+Y(i,5)*DA10_1_6+
     &           Y(i,6)*DA10_1_7+Y(i,7)*DA10_1_8+Y(i,8)*DA10_1_9
            Z(i,2)=y0(i)* DA10_2_1+Y(i,1)*DA10_2_2+Y(i,2)*DA10_2_3+
     &           Y(i,3)*DA10_2_4+Y(i,4)*DA10_2_5+Y(i,5)*DA10_2_6+
     &           Y(i,6)*DA10_2_7+Y(i,7)*DA10_2_8+Y(i,8)*DA10_2_9
            Z(i,3)=y0(i)* DA10_3_1+Y(i,1)*DA10_3_2+Y(i,2)*DA10_3_3+
     &           Y(i,3)*DA10_3_4+Y(i,4)*DA10_3_5+Y(i,5)*DA10_3_6+
     &           Y(i,6)*DA10_3_7+Y(i,7)*DA10_3_8+Y(i,8)*DA10_3_9
            Z(i,4)=y0(i)* DA10_4_1+Y(i,1)*DA10_4_2+Y(i,2)*DA10_4_3+
     &           Y(i,3)*DA10_4_4+Y(i,4)*DA10_4_5+Y(i,5)*DA10_4_6+
     &           Y(i,6)*DA10_4_7+Y(i,7)*DA10_4_8+Y(i,8)*DA10_4_9
            Z(i,5)=y0(i)* DA10_5_1+Y(i,1)*DA10_5_2+Y(i,2)*DA10_5_3+
     &           Y(i,3)*DA10_5_4+Y(i,4)*DA10_5_5+Y(i,5)*DA10_5_6+
     &           Y(i,6)*DA10_5_7+Y(i,7)*DA10_5_8+Y(i,8)*DA10_5_9
            Z(i,6)=y0(i)* DA10_6_1+Y(i,1)*DA10_6_2+Y(i,2)*DA10_6_3+
     &           Y(i,3)*DA10_6_4+Y(i,4)*DA10_6_5+Y(i,5)*DA10_6_6+
     &           Y(i,6)*DA10_6_7+Y(i,7)*DA10_6_8+Y(i,8)*DA10_6_9
            Z(i,7)=y0(i)* DA10_7_1+Y(i,1)*DA10_7_2+Y(i,2)*DA10_7_3+
     &           Y(i,3)*DA10_7_4+Y(i,4)*DA10_7_5+Y(i,5)*DA10_7_6+
     &           Y(i,6)*DA10_7_7+Y(i,7)*DA10_7_8+Y(i,8)*DA10_7_9
            Z(i,8)=y0(i)* DA10_8_1+Y(i,1)*DA10_8_2+Y(i,2)*DA10_8_3+
     &           Y(i,3)*DA10_8_4+Y(i,4)*DA10_8_5+Y(i,5)*DA10_8_6+
     &           Y(i,6)*DA10_8_7+Y(i,7)*DA10_8_8+Y(i,8)*DA10_8_9
      end do

c     MZ_i=M0*Z_i, i=1,2,...,8

      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,1),MZ(1,1),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,2),MZ(1,2),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,3),MZ(1,3),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,4),MZ(1,4),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,5),MZ(1,5),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,6),MZ(1,6),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,7),MZ(1,7),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,8),MZ(1,8),ijob(2))

c     MZ=MZ-h*[f0 F]*(dB)'
      do i=1,m
          MZ(i,1) = MZ(i,1)-
     &        h*(f0(i)* DB10_1_1+F(i,1)*DB10_1_2+F(i,2)*DB10_1_3+
     &           F(i,3)*DB10_1_4+F(i,4)*DB10_1_5+F(i,5)*DB10_1_6+
     &           F(i,6)*DB10_1_7+F(i,7)*DB10_1_8+F(i,8)*DB10_1_9)
          MZ(i,2) = MZ(i,2)-
     &        h*(f0(i)* DB10_2_1+F(i,1)*DB10_2_2+F(i,2)*DB10_2_3+
     &           F(i,3)*DB10_2_4+F(i,4)*DB10_2_5+F(i,5)*DB10_2_6+
     &           F(i,6)*DB10_2_7+F(i,7)*DB10_2_8+F(i,8)*DB10_2_9)
          MZ(i,3) = MZ(i,3)-
     &        h*(f0(i)* DB10_3_1+F(i,1)*DB10_3_2+F(i,2)*DB10_3_3+
     &           F(i,3)*DB10_3_4+F(i,4)*DB10_3_5+F(i,5)*DB10_3_6+
     &           F(i,6)*DB10_3_7+F(i,7)*DB10_3_8+F(i,8)*DB10_3_9)
          MZ(i,4) = MZ(i,4)-
     &        h*(f0(i)* DB10_4_1+F(i,1)*DB10_4_2+F(i,2)*DB10_4_3+
     &           F(i,3)*DB10_4_4+F(i,4)*DB10_4_5+F(i,5)*DB10_4_6+
     &           F(i,6)*DB10_4_7+F(i,7)*DB10_4_8+F(i,8)*DB10_4_9)
          MZ(i,5) = MZ(i,5)-
     &        h*(f0(i)* DB10_5_1+F(i,1)*DB10_5_2+F(i,2)*DB10_5_3+
     &           F(i,3)*DB10_5_4+F(i,4)*DB10_5_5+F(i,5)*DB10_5_6+
     &           F(i,6)*DB10_5_7+F(i,7)*DB10_5_8+F(i,8)*DB10_5_9)
          MZ(i,6) = MZ(i,6)-
     &        h*(f0(i)* DB10_6_1+F(i,1)*DB10_6_2+F(i,2)*DB10_6_3+
     &           F(i,3)*DB10_6_4+F(i,4)*DB10_6_5+F(i,5)*DB10_6_6+
     &           F(i,6)*DB10_6_7+F(i,7)*DB10_6_8+F(i,8)*DB10_6_9)
          MZ(i,7) = MZ(i,7)-
     &        h*(f0(i)* DB10_7_1+F(i,1)*DB10_7_2+F(i,2)*DB10_7_3+
     &           F(i,3)*DB10_7_4+F(i,4)*DB10_7_5+F(i,5)*DB10_7_6+
     &           F(i,6)*DB10_7_7+F(i,7)*DB10_7_8+F(i,8)*DB10_7_9)
          MZ(i,8) = MZ(i,8)-
     &        h*(f0(i)* DB10_8_1+F(i,1)*DB10_8_2+F(i,2)*DB10_8_3+
     &           F(i,3)*DB10_8_4+F(i,4)*DB10_8_5+F(i,5)*DB10_8_6+
     &           F(i,6)*DB10_8_7+F(i,7)*DB10_8_8+F(i,8)*DB10_8_9)
      end do

c     theta*MZ=MZ
      call sollu(m,theta,ldlu,MZ(1,1),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,2),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,3),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,4),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,5),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,6),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,7),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,8),mljac,mujac,ipvt,ijob(1))

c      MZ=MZ+[y0 Y]*(A2)'
      do i=1,m
            MZ(i,1)=MZ(i,1)+
     &               y0(i)* A210_1_1+Y(i,1)*A210_1_2+Y(i,2)*A210_1_3+
     &              Y(i,3)* A210_1_4+Y(i,4)*A210_1_5+Y(i,5)*A210_1_6+
     &              Y(i,6)* A210_1_7+Y(i,7)*A210_1_8+Y(i,8)*A210_1_9
            MZ(i,2)=MZ(i,2)+
     &               y0(i)* A210_2_1+Y(i,1)*A210_2_2+Y(i,2)*A210_2_3+
     &              Y(i,3)* A210_2_4+Y(i,4)*A210_2_5+Y(i,5)*A210_2_6+
     &              Y(i,6)* A210_2_7+Y(i,7)*A210_2_8+Y(i,8)*A210_2_9
            MZ(i,3)=MZ(i,3)+
     &               y0(i)* A210_3_1+Y(i,1)*A210_3_2+Y(i,2)*A210_3_3+
     &              Y(i,3)* A210_3_4+Y(i,4)*A210_3_5+Y(i,5)*A210_3_6+
     &              Y(i,6)* A210_3_7+Y(i,7)*A210_3_8+Y(i,8)*A210_3_9
            MZ(i,4)=MZ(i,4)+
     &               y0(i)* A210_4_1+Y(i,1)*A210_4_2+Y(i,2)*A210_4_3+
     &              Y(i,3)* A210_4_4+Y(i,4)*A210_4_5+Y(i,5)*A210_4_6+
     &              Y(i,6)* A210_4_7+Y(i,7)*A210_4_8+Y(i,8)*A210_4_9
            MZ(i,5)=MZ(i,5)+
     &               y0(i)* A210_5_1+Y(i,1)*A210_5_2+Y(i,2)*A210_5_3+
     &              Y(i,3)* A210_5_4+Y(i,4)*A210_5_5+Y(i,5)*A210_5_6+
     &              Y(i,6)* A210_5_7+Y(i,7)*A210_5_8+Y(i,8)*A210_5_9
            MZ(i,6)=MZ(i,6)+
     &               y0(i)* A210_6_1+Y(i,1)*A210_6_2+Y(i,2)*A210_6_3+
     &              Y(i,3)* A210_6_4+Y(i,4)*A210_6_5+Y(i,5)*A210_6_6+
     &              Y(i,6)* A210_6_7+Y(i,7)*A210_6_8+Y(i,8)*A210_6_9
            MZ(i,7)=MZ(i,7)+
     &               y0(i)* A210_7_1+Y(i,1)*A210_7_2+Y(i,2)*A210_7_3+
     &              Y(i,3)* A210_7_4+Y(i,4)*A210_7_5+Y(i,5)*A210_7_6+
     &              Y(i,6)* A210_7_7+Y(i,7)*A210_7_8+Y(i,8)*A210_7_9
            MZ(i,8)=MZ(i,8)+
     &               y0(i)* A210_8_1+Y(i,1)*A210_8_2+Y(i,2)*A210_8_3+
     &              Y(i,3)* A210_8_4+Y(i,4)*A210_8_5+Y(i,5)*A210_8_6+
     &              Y(i,6)* A210_8_7+Y(i,7)*A210_8_8+Y(i,8)*A210_8_9
      end do

c     Z_i=M0*MZ_i, i=1,2,...,8

      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,1),Z(1,1),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,2),Z(1,2),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,3),Z(1,3),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,4),Z(1,4),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,5),Z(1,5),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,6),Z(1,6),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,7),Z(1,7),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,8),Z(1,8),ijob(2))

c     Z=Z-h*[f0 F]*(B2)'
      do i=1,m
          Z(i,1)=Z(i,1)-h*(f0(i)*B210_1_1 + F(i,1)*gamma)
          Z(i,2)=Z(i,2)-h*(f0(i)*B210_2_1 + F(i,2)*gamma)
          Z(i,3)=Z(i,3)-h*(f0(i)*B210_3_1 + F(i,3)*gamma)
          Z(i,4)=Z(i,4)-h*(f0(i)*B210_4_1 + F(i,4)*gamma)
          Z(i,5)=Z(i,5)-h*(f0(i)*B210_5_1 + F(i,5)*gamma)
          Z(i,6)=Z(i,6)-h*(f0(i)*B210_6_1 + F(i,6)*gamma)
          Z(i,7)=Z(i,7)-h*(f0(i)*B210_7_1 + F(i,7)*gamma)
          Z(i,8)=Z(i,8)-h*(f0(i)*B210_8_1 + F(i,8)*gamma)
      end do

c     theta*Z=Z
      call sollu(m,theta,ldlu,Z(1,1),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,2),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,3),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,4),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,5),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,6),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,7),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,8),mljac,mujac,ipvt,ijob(1))

30    continue
c
c     Y = Y-Z
c
      do j=1,k
         do i=1,m
           Y(i,j) = Y(i,j) - Z(i,j)
         end do
      end do

      return
      end

C -------------------------------------------------------------------

      subroutine blendstep12(m,y0,f0,Y,F,h,theta,ipvt,Z,gamma,
     &                      ldlu,mljac,mujac,ijob,imas,
     &                                ldmas,mlmas,mumas,M0,MZ)
c
c     Blended iteration for the 12th order method
c

      implicit none

      integer k
      parameter (k=10)

c     Input parameters
      integer m,ipvt(m),ldlu,mljac,mujac,ijob(2),imas,ldmas,mlmas,mumas
      double precision  y0(m),f0(m),theta(ldlu,m),gamma,h,M0(ldmas,m)

c     I/O parameters
      double precision Y(m,k),F(m,k),MZ(m,k)

c     Output parameters
      double precision Z(m,k)

c     Local variables
      integer i,j

C ---------------------------------------------------------------------------------------
C 12th order BIM
C ---------------------------------------------------------------------------------------

      double precision
     &   DA12_1_1,DA12_1_2,DA12_1_3,DA12_1_4,DA12_1_5,DA12_1_6,DA12_1_7,
     &   DA12_1_8,DA12_1_9,DA12_1_10,DA12_1_11,
     &   DA12_2_1,DA12_2_2,DA12_2_3,DA12_2_4,DA12_2_5,DA12_2_6,DA12_2_7,
     &   DA12_2_8,DA12_2_9,DA12_2_10,DA12_2_11,
     &   DA12_3_1,DA12_3_2,DA12_3_3,DA12_3_4,DA12_3_5,DA12_3_6,DA12_3_7,
     &   DA12_3_8,DA12_3_9,DA12_3_10,DA12_3_11,
     &   DA12_4_1,DA12_4_2,DA12_4_3,DA12_4_4,DA12_4_5,DA12_4_6,DA12_4_7,
     &   DA12_4_8,DA12_4_9,DA12_4_10,DA12_4_11,
     &   DA12_5_1,DA12_5_2,DA12_5_3,DA12_5_4,DA12_5_5,DA12_5_6,DA12_5_7,
     &   DA12_5_8,DA12_5_9,DA12_5_10,DA12_5_11,
     &   DA12_6_1,DA12_6_2,DA12_6_3,DA12_6_4,DA12_6_5,DA12_6_6,DA12_6_7,
     &   DA12_6_8,DA12_6_9,DA12_6_10,DA12_6_11,
     &   DA12_7_1,DA12_7_2,DA12_7_3,DA12_7_4,DA12_7_5,DA12_7_6,DA12_7_7,
     &   DA12_7_8,DA12_7_9,DA12_7_10,DA12_7_11,
     &   DA12_8_1,DA12_8_2,DA12_8_3,DA12_8_4,DA12_8_5,DA12_8_6,DA12_8_7,
     &   DA12_8_8,DA12_8_9,DA12_8_10,DA12_8_11,
     &   DA12_9_1,DA12_9_2,DA12_9_3,DA12_9_4,DA12_9_5,DA12_9_6,DA12_9_7,
     &   DA12_9_8,DA12_9_9,DA12_9_10,DA12_9_11,
     &   DA12_10_1,DA12_10_2,DA12_10_3,DA12_10_4,DA12_10_5,DA12_10_6,
     &   DA12_10_7,DA12_10_8,DA12_10_9,DA12_10_10,DA12_10_11,
     &   DB12_1_1,DB12_1_2,DB12_1_3,DB12_1_4,DB12_1_5,DB12_1_6,DB12_1_7,
     &   DB12_1_8,DB12_1_9,DB12_1_10,DB12_1_11,
     &   DB12_2_1,DB12_2_2,DB12_2_3,DB12_2_4,DB12_2_5,DB12_2_6,DB12_2_7,
     &   DB12_2_8,DB12_2_9,DB12_2_10,DB12_2_11,
     &   DB12_3_1,DB12_3_2,DB12_3_3,DB12_3_4,DB12_3_5,DB12_3_6,DB12_3_7,
     &   DB12_3_8,DB12_3_9,DB12_3_10,DB12_3_11,
     &   DB12_4_1,DB12_4_2,DB12_4_3,DB12_4_4,DB12_4_5,DB12_4_6,DB12_4_7,
     &   DB12_4_8,DB12_4_9,DB12_4_10,DB12_4_11,
     &   DB12_5_1,DB12_5_2,DB12_5_3,DB12_5_4,DB12_5_5,DB12_5_6,DB12_5_7,
     &   DB12_5_8,DB12_5_9,DB12_5_10,DB12_5_11,
     &   DB12_6_1,DB12_6_2,DB12_6_3,DB12_6_4,DB12_6_5,DB12_6_6,DB12_6_7,
     &   DB12_6_8,DB12_6_9,DB12_6_10,DB12_6_11,
     &   DB12_7_1,DB12_7_2,DB12_7_3,DB12_7_4,DB12_7_5,DB12_7_6,DB12_7_7,
     &   DB12_7_8,DB12_7_9,DB12_7_10,DB12_7_11,
     &   DB12_8_1,DB12_8_2,DB12_8_3,DB12_8_4,DB12_8_5,DB12_8_6,DB12_8_7,
     &   DB12_8_8,DB12_8_9,DB12_8_10,DB12_8_11,
     &   DB12_9_1,DB12_9_2,DB12_9_3,DB12_9_4,DB12_9_5,DB12_9_6,DB12_9_7,
     &   DB12_9_8,DB12_9_9,DB12_9_10,DB12_9_11,
     &   DB12_10_1,DB12_10_2,DB12_10_3,DB12_10_4,DB12_10_5,DB12_10_6,
     &   DB12_10_7,DB12_10_8,DB12_10_9,DB12_10_10,DB12_10_11,
     &   A212_1_1,A212_1_2,A212_1_3,A212_1_4,A212_1_5,A212_1_6,A212_1_7,
     &   A212_1_8,A212_1_9,A212_1_10,A212_1_11,
     &   A212_2_1,A212_2_2,A212_2_3,A212_2_4,A212_2_5,A212_2_6,A212_2_7,
     &   A212_2_8,A212_2_9,A212_2_10,A212_2_11,
     &   A212_3_1,A212_3_2,A212_3_3,A212_3_4,A212_3_5,A212_3_6,A212_3_7,
     &   A212_3_8,A212_3_9,A212_3_10,A212_3_11,
     &   A212_4_1,A212_4_2,A212_4_3,A212_4_4,A212_4_5,A212_4_6,A212_4_7,
     &   A212_4_8,A212_4_9,A212_4_10,A212_4_11,
     &   A212_5_1,A212_5_2,A212_5_3,A212_5_4,A212_5_5,A212_5_6,A212_5_7,
     &   A212_5_8,A212_5_9,A212_5_10,A212_5_11,
     &   A212_6_1,A212_6_2,A212_6_3,A212_6_4,A212_6_5,A212_6_6,A212_6_7,
     &   A212_6_8,A212_6_9,A212_6_10,A212_6_11,
     &   A212_7_1,A212_7_2,A212_7_3,A212_7_4,A212_7_5,A212_7_6,A212_7_7,
     &   A212_7_8,A212_7_9,A212_7_10,A212_7_11,
     &   A212_8_1,A212_8_2,A212_8_3,A212_8_4,A212_8_5,A212_8_6,A212_8_7,
     &   A212_8_8,A212_8_9,A212_8_10,A212_8_11,
     &   A212_9_1,A212_9_2,A212_9_3,A212_9_4,A212_9_5,A212_9_6,A212_9_7,
     &   A212_9_8,A212_9_9,A212_9_10,A212_9_11,
     &   A212_10_1,A212_10_2,A212_10_3,A212_10_4,A212_10_5,A212_10_6,
     &   A212_10_7,A212_10_8,A212_10_9,A212_10_10,A212_10_11,
     &   B212_1_1,B212_2_1,B212_3_1,B212_4_1,B212_5_1,B212_6_1,B212_7_1,
     &   B212_8_1,B212_9_1,B212_10_1


      parameter(
     &  DA12_1_1    =   -2477009404842909D0/2D15,
     &  DA12_1_2    =   +144472720252463D0/45D12,
     &  DA12_1_3    =   -104423832752463D0/2D13,
     &  DA12_1_4    =   +9994398083607D0/125D10,
     &  DA12_1_5    =   -198624609755747D0/2D13,
     &  DA12_1_6    =   +578987204267241D0/625D11,
     &  DA12_1_7    =   -63206136585249D0/1D13,
     &  DA12_1_8    =   +3823813464403D0/125D10,
     &  DA12_1_9    =   -79610832752463D0/8D13,
     &  DA12_1_10   =   +976467842623D0/5D12,
     &  DA12_1_11   =   -78691832752463D0/45D14,
     &  DA12_2_1    =   -2813274212551D0/29296875D5,
     &  DA12_2_2    =   -1504054699D0/87890625D3,
     &  DA12_2_3    =   +334937703171D0/15625D7,
     &  DA12_2_4    =   -17253955733D0/732421875D1,
     &  DA12_2_5    =   +137190949259D0/5859375D4,
     &  DA12_2_6    =   -735194523679D0/3662109375D2,
     &  DA12_2_7    =   +153604208387D0/1171875D5,
     &  DA12_2_8    =   -4522677457D0/732421875D1,
     &  DA12_2_9    =   +23115833993D0/1171875D5,
     &  DA12_2_10   =   -1119363137D0/29296875D3,
     &  DA12_2_11   =   +119096296921D0/3515625D7,
     &  DA12_3_1    =   -4387223158569623D0/6D15,
     &  DA12_3_2    =   -14076290383771D0/15D12,
     &  DA12_3_3    =   +45847433651313D0/2D13,
     &  DA12_3_4    =   -8091971633771D0/375D10,
     &  DA12_3_5    =   +73204095186397D0/2D13,
     &  DA12_3_6    =   -244942223059191D0/625D11,
     &  DA12_3_7    =   +84461845186397D0/3D13,
     &  DA12_3_8    =   -1752429876253D0/125D10,
     &  DA12_3_9    =   +37162883651313D0/8D13,
     &  DA12_3_10   =   -461779273473D0/5D12,
     &  DA12_3_11   =   +12525477883771D0/15D14,
     &  DA12_4_1    =   -1861031701487D0/146484375D4,
     &  DA12_4_2    =   +163653743573D0/17578125D4,
     &  DA12_4_3    =   -167243587323D0/78125D6,
     &  DA12_4_4    =   +19779235397D0/48828125D2,
     &  DA12_4_5    =   -105419144767D0/29296875D3,
     &  DA12_4_6    =   +944544955011D0/244140625D3,
     &  DA12_4_7    =   -113325741529D0/390625D5,
     &  DA12_4_8    =   +21327766939D0/146484375D2,
     &  DA12_4_9    =   -9443080653D0/1953125D4,
     &  DA12_4_10   =   +5635786799D0/5859375D4,
     &  DA12_4_11   =   -152884212323D0/17578125D6,
     &  DA12_5_1    =   -209971213D0/2048D5,
     &  DA12_5_2    =   +375871D0/4608D3,
     &  DA12_5_3    =   -980573D0/6144D3,
     &  DA12_5_4    =   +76277D0/384D3,
     &  DA12_5_5    =   +456743D0/6144D3,
     &  DA12_5_6    =   +9236953D0/64D5,
     &  DA12_5_7    =   -2592499D0/3072D3,
     &  DA12_5_8    =   +116713D0/384D3,
     &  DA12_5_9    =   -2156893D0/24576D3,
     &  DA12_5_10   =   +8271D0/512D3,
     &  DA12_5_11   =   -640543D0/4608D5,
     &  DA12_6_1    =   -4727199180623D0/5859375D6,
     &  DA12_6_2    =   -38397576209D0/5859375D4,
     &  DA12_6_3    =   +227154597879D0/15625D7,
     &  DA12_6_4    =   -36961638709D0/146484375D2,
     &  DA12_6_5    =   +123083508919D0/390625D5,
     &  DA12_6_6    =   -625420975389D0/244140625D3,
     &  DA12_6_7    =   +360841627213D0/1171875D5,
     &  DA12_6_8    =   -7331573387D0/48828125D2,
     &  DA12_6_9    =   +34452186063D0/78125D6,
     &  DA12_6_10   =   -4904591801D0/5859375D4,
     &  DA12_6_11   =   +86487730543D0/1171875D7,
     &  DA12_7_1    =   -6527365601261023D0/6D15,
     &  DA12_7_2    =   +13219473015287D0/45D12,
     &  DA12_7_3    =   -12656585515287D0/2D13,
     &  DA12_7_4    =   +1312472696143D0/125D10,
     &  DA12_7_5    =   -72835248607009D0/6D13,
     &  DA12_7_6    =   +53134186107009D0/625D11,
     &  DA12_7_7    =   +663222376999D0/1D13,
     &  DA12_7_8    =   +3863677395041D0/375D10,
     &  DA12_7_9    =   -32920535515287D0/8D13,
     &  DA12_7_10   =   +906563815381D0/15D12,
     &  DA12_7_11   =   -21662785515287D0/45D14,
     &  DA12_8_1    =   -1017354514167D0/9765625D5,
     &  DA12_8_2    =   +3448272719D0/2197265625D1,
     &  DA12_8_3    =   -1517467504D0/3662109375D0,
     &  DA12_8_4    =   +3256140146D0/3662109375D0,
     &  DA12_8_5    =   -85969177441D0/5859375D4,
     &  DA12_8_6    =   +172929352099D0/91552734375D0,
     &  DA12_8_7    =   -14493362443D0/732421875D1,
     &  DA12_8_8    =   +7139238634D0/3662109375D0,
     &  DA12_8_9    =   5459539781D0/390625D5,
     &  DA12_8_10   =   -311981201D0/244140625D1,
     &  DA12_8_11   =   +3176506039D0/54931640625D1,
     &  DA12_9_1    =   -1883458029218509D0/2D15,
     &  DA12_9_2    =   -1274779497593D0/5D12,
     &  DA12_9_3    =   +16125452978337D0/2D13,
     &  DA12_9_4    =   -2480966997593D0/125D10,
     &  DA12_9_5    =   +72364256949453D0/2D13,
     &  DA12_9_6    =   -308280545848359D0/625D11,
     &  DA12_9_7    =   +51140018983151D0/1D13,
     &  DA12_9_8    =   -5198771570397D0/125D10,
     &  DA12_9_9    =   +239442452978337D0/8D13,
     &  DA12_9_10   =   -979844944177D0/5D12,
     &  DA12_9_11   =   -31292283002407D0/5D14,
     &  DA12_10_1   =   -106433D0/1D5,
     &  DA12_10_2   =   +6433D0/9D3,
     &  DA12_10_3   =   -57897D0/16D3,
     &  DA12_10_4   =   +2757D0/25D1,
     &  DA12_10_5   =   -45031D0/2D3,
     &  DA12_10_6   =   +405279D0/125D2,
     &  DA12_10_7   =   -135093D0/4D3,
     &  DA12_10_8   =   +6433D0/25D1,
     &  DA12_10_9   =   -57897D0/4D3,
     &  DA12_10_10  =   +6433D0/1D3,
     &  DA12_10_11  =   -3183139D0/36D5)



      parameter(
     &  A212_1_1    =   477009404842909D0/2D15,
     &  A212_1_2    =   -99472720252463D0/45D12,
     &  A212_1_3    =   104423832752463D0/2D13,
     &  A212_1_4    =   -9994398083607D0/125D10,
     &  A212_1_5    =   198624609755747D0/2D13,
     &  A212_1_6    =   -578987204267241D0/625D11,
     &  A212_1_7    =   63206136585249D0/1D13,
     &  A212_1_8    =   -3823813464403D0/125D10,
     &  A212_1_9    =   79610832752463D0/8D13,
     &  A212_1_10   =   -976467842623D0/5D12,
     &  A212_1_11   =   78691832752463D0/45D14,
     &  A212_2_1    =   -116413287449D0/29296875D5,
     &  A212_2_2    =   1504054699D0/87890625D3,
     &  A212_2_3    =   -178687703171D0/15625D7,
     &  A212_2_4    =   17253955733D0/732421875D1,
     &  A212_2_5    =   -137190949259D0/5859375D4,
     &  A212_2_6    =   735194523679D0/3662109375D2,
     &  A212_2_7    =   -153604208387D0/1171875D5,
     &  A212_2_8    =   4522677457D0/732421875D1,
     &  A212_2_9    =   -23115833993D0/1171875D5,
     &  A212_2_10   =   1119363137D0/29296875D3,
     &  A212_2_11   =   -119096296921D0/3515625D7,
     &  A212_3_1    =   -1612776841430377D0/6D15,
     &  A212_3_2    =   14076290383771D0/15D12,
     &  A212_3_3    =   -45847433651313D0/2D13,
     &  A212_3_4    =   11841971633771D0/375D10,
     &  A212_3_5    =   -73204095186397D0/2D13,
     &  A212_3_6    =   244942223059191D0/625D11,
     &  A212_3_7    =   -84461845186397D0/3D13,
     &  A212_3_8    =   1752429876253D0/125D10,
     &  A212_3_9    =   -37162883651313D0/8D13,
     &  A212_3_10   =   461779273473D0/5D12,
     &  A212_3_11   =   -12525477883771D0/15D14,
     &  A212_4_1    =   396187951487D0/146484375D4,
     &  A212_4_2    =   -163653743573D0/17578125D4,
     &  A212_4_3    =   167243587323D0/78125D6,
     &  A212_4_4    =   -19779235397D0/48828125D2,
     &  A212_4_5    =   134716019767D0/29296875D3,
     &  A212_4_6    =   -944544955011D0/244140625D3,
     &  A212_4_7    =   113325741529D0/390625D5,
     &  A212_4_8    =   -21327766939D0/146484375D2,
     &  A212_4_9    =   9443080653D0/1953125D4,
     &  A212_4_10   =   -5635786799D0/5859375D4,
     &  A212_4_11   =   152884212323D0/17578125D6,
     &  A212_5_1    =   5171213D0/2048D5,
     &  A212_5_2    =   -375871D0/4608D3,
     &  A212_5_3    =   980573D0/6144D3,
     &  A212_5_4    =   -76277D0/384D3,
     &  A212_5_5    =   -456743D0/6144D3,
     &  A212_5_6    =   -2836953D0/64D5,
     &  A212_5_7    =   2592499D0/3072D3,
     &  A212_5_8    =   -116713D0/384D3,
     &  A212_5_9    =   2156893D0/24576D3,
     &  A212_5_10   =   -8271D0/512D3,
     &  A212_5_11   =   640543D0/4608D5,
     &  A212_6_1    =   -1132175819377D0/5859375D6,
     &  A212_6_2    =   38397576209D0/5859375D4,
     &  A212_6_3    =   -227154597879D0/15625D7,
     &  A212_6_4    =   36961638709D0/146484375D2,
     &  A212_6_5    =   -123083508919D0/390625D5,
     &  A212_6_6    =   625420975389D0/244140625D3,
     &  A212_6_7    =   -243654127213D0/1171875D5,
     &  A212_6_8    =   7331573387D0/48828125D2,
     &  A212_6_9    =   -34452186063D0/78125D6,
     &  A212_6_10   =   4904591801D0/5859375D4,
     &  A212_6_11   =   -86487730543D0/1171875D7,
     &  A212_7_1    =   527365601261023D0/6D15,
     &  A212_7_2    =   -13219473015287D0/45D12,
     &  A212_7_3    =   12656585515287D0/2D13,
     &  A212_7_4    =   -1312472696143D0/125D10,
     &  A212_7_5    =   72835248607009D0/6D13,
     &  A212_7_6    =   -53134186107009D0/625D11,
     &  A212_7_7    =   -663222376999D0/1D13,
     &  A212_7_8    =   -113677395041D0/375D10,
     &  A212_7_9    =   32920535515287D0/8D13,
     &  A212_7_10   =   -906563815381D0/15D12,
     &  A212_7_11   =   21662785515287D0/45D14,
     &  A212_8_1    =   40792014167D0/9765625D5,
     &  A212_8_2    =   -3448272719D0/2197265625D1,
     &  A212_8_3    =   1517467504D0/3662109375D0,
     &  A212_8_4    =   -3256140146D0/3662109375D0,
     &  A212_8_5    =   85969177441D0/5859375D4,
     &  A212_8_6    =   -172929352099D0/91552734375D0,
     &  A212_8_7    =   14493362443D0/732421875D1,
     &  A212_8_8    =   -7139238634D0/3662109375D0,
     &  A212_8_9    =   33602960219D0/390625D5,
     &  A212_8_10   =   311981201D0/244140625D1,
     &  A212_8_11   =   -3176506039D0/54931640625D1,
     &  A212_9_1    =   -116541970781491D0/2D15,
     &  A212_9_2    =   1274779497593D0/5D12,
     &  A212_9_3    =   -16125452978337D0/2D13,
     &  A212_9_4    =   2480966997593D0/125D10,
     &  A212_9_5    =   -72364256949453D0/2D13,
     &  A212_9_6    =    308280545848359D0/625D11,
     &  A212_9_7    =   -51140018983151D0/1D13,
     &  A212_9_8    =   5198771570397D0/125D10,
     &  A212_9_9    =   -239442452978337D0/8D13,
     &  A212_9_10   =   5979844944177D0/5D12,
     &  A212_9_11   =   31292283002407D0/5D14,
     &  A212_10_1   =   6433D0/1D5,
     &  A212_10_2   =   -6433D0/9D3,
     &  A212_10_3   =   57897D0/16D3,
     &  A212_10_4   =   -2757D0/25D1,
     &  A212_10_5   =   45031D0/2D3,
     &  A212_10_6   =   -405279D0/125D2,
     &  A212_10_7   =   135093D0/4D3,
     &  A212_10_8   =   -6433D0/25D1,
     &  A212_10_9   =   57897D0/4D3,
     &  A212_10_10  =   -6433D0/1D3,
     &  A212_10_11  =   6783139D0/36D5)



      parameter(
     &  DB12_1_1    =   +11479451563170288931D0/2953665D13,
     &  DB12_1_2    =   +63885874501D0/9451728D4,
     &  DB12_1_3    =   -57527974223D0/352864512D2,
     &  DB12_1_4    =   +2126874941D0/88216128D1,
     &  DB12_1_5    =   -5375477389D0/19603584D2,
     &  DB12_1_6    =   +5059789129D0/22054032D2,
     &  DB12_1_7    =   -2704182763D0/19603584D2,
     &  DB12_1_8    =   +2542514071D0/44108064D2,
     &  DB12_1_9    =   -1112711029D0/705729024D1,
     &  DB12_1_10   =   +663357529D0/264648384D2,
     &  DB12_1_11   =   -36358589D0/2117187072D1,
     &  DB12_2_1    =   +43464304639852711D0/17306630859375D4,
     &  DB12_2_2    =   +443438449D0/24810786D1,
     &  DB12_2_3    =   -3078028627D0/172297125D1,
     &  DB12_2_4    =   +313340009D0/103378275D0,
     &  DB12_2_5    =   -3589321627D0/8270262D2,
     &  DB12_2_6    =   +61699387D0/1378377D1,
     &  DB12_2_7    =   -549273503D0/16540524D1,
     &  DB12_2_8    =   +178915433D0/103378275D0,
     &  DB12_2_9    =   -332030663D0/5513508D2,
     &  DB12_2_10   =   +156015413D0/12405393D2,
     &  DB12_2_11   =   -59216207D0/49621572D2,
     &  DB12_3_1    =   +38675312700454505333D0/20675655D13,
     &  DB12_3_2    =   +44161715831D0/264648384D2,
     &  DB12_3_3    =   -3283281181D0/117621504D2,
     &  DB12_3_4    =   +221229581497D0/11027016D4,
     &  DB12_3_5    =   -9249915641D0/352864512D1,
     &  DB12_3_6    =   +328488001D0/14702688D1,
     &  DB12_3_7    =   -24501545899D0/176432256D2,
     &  DB12_3_8    =   +2674826191D0/44108064D2,
     &  DB12_3_9    =   -2067552667D0/117621504D2,
     &  DB12_3_10   =   +160146433D0/529296768D1,
     &  DB12_3_11   =   -246451241D0/1058593536D2,
     &  DB12_4_1    =   +148916128991770691D0/40382138671875D4,
     &  DB12_4_2    =   +174240916D0/103378275D0,
     &  DB12_4_3    =   -2553949D0/6891885D0,
     &  DB12_4_4    =   +112512368D0/34459425D0,
     &  DB12_4_5    =   -4311283649D0/153153D4,
     &  DB12_4_6    =   +79999624D0/34459425D0,
     &  DB12_4_7    =   -5864678D0/3828825D0,
     &  DB12_4_8    =   +4859248D0/6891885D0,
     &  DB12_4_9    =   -7453148D0/34459425D0,
     &  DB12_4_10   =   +587116D0/14768325D0,
     &  DB12_4_11   =   -342427D0/103378275D0,
     &  DB12_5_1    =   +17710238627183D0/6351561216D4,
     &  DB12_5_2    =   +5562130625D0/3175780608D0,
     &  DB12_5_3    =   -935636375D0/1411458048D0,
     &  DB12_5_4    =   +2125610875D0/529296768D0,
     &  DB12_5_5    =   -6501442375D0/2117187072D0,
     &  DB12_5_6    =   +221113406161D0/5513508D4,
     &  DB12_5_7    =   -6539730625D0/2117187072D0,
     &  DB12_5_8    =   +829280125D0/529296768D0,
     &  DB12_5_9    =   -107876375D0/201636864D0,
     &  DB12_5_10   =   +349458875D0/3175780608D0,
     &  DB12_5_11   =   -130922555D0/12703122432D0,
     &  DB12_6_1    =   +84117653394000791D0/40382138671875D4,
     &  DB12_6_2    =   +706138327D0/4135131D2,
     &  DB12_6_3    =   -86081069D0/1837836D2,
     &  DB12_6_4    =   +24183767D0/6891885D0,
     &  DB12_6_5    =   -615563609D0/2756754D2,
     &  DB12_6_6    =   +94460189D0/2297295D1,
     &  DB12_6_7    =   -30823666891D0/1378377D4,
     &  DB12_6_8    =   +32702707D0/34459425D0,
     &  DB12_6_9    =   -1649311D0/525096D1,
     &  DB12_6_10   =   +25613911D0/4135131D2,
     &  DB12_6_11   =   -1836449D0/33081048D1,
     &  DB12_7_1    =   +9066850209038767819D0/2953665D13,
     &  DB12_7_2    =   +1273945301D0/75613824D1,
     &  DB12_7_3    =   -1831435457D0/50409216D2,
     &  DB12_7_4    =   +2028629887D0/6301152D2,
     &  DB12_7_5    =   -475125427D0/2800512D2,
     &  DB12_7_6    =   +213226433D0/6301152D1,
     &  DB12_7_7    =   -25983713D0/5601024D1,
     &  DB12_7_8    =   +8396816521D0/1575288D4,
     &  DB12_7_9    =   -1310408951D0/50409216D2,
     &  DB12_7_10   =   +175300267D0/37806912D2,
     &  DB12_7_11   =   -58184383D0/151227648D2,
     &  DB12_8_1    =   +49560837081155911D0/17306630859375D4,
     &  DB12_8_2    =   +76433504D0/44304975D0,
     &  DB12_8_3    =   -18894032D0/34459425D0,
     &  DB12_8_4    =   +385150592D0/103378275D0,
     &  DB12_8_5    =   -54124624D0/20675655D0,
     &  DB12_8_6    =   +6293440D0/1378377D0,
     &  DB12_8_7    =   -165253856D0/103378275D0,
     &  DB12_8_8    =   +250726016D0/103378275D0,
     &  DB12_8_9    =   -8968348841D0/1378377D4,
     &  DB12_8_10   =   +3948064D0/62026965D0,
     &  DB12_8_11   =   -1974032D0/310134825D0,
     &  DB12_9_1    =   +2538784308052589777D0/984555D13,
     &  DB12_9_2    =   +14976693257D0/88216128D2,
     &  DB12_9_3    =   -326374859D0/78414336D1,
     &  DB12_9_4    =   +4913986843D0/14702688D2,
     &  DB12_9_5    =   -10987159111D0/58810752D2,
     &  DB12_9_6    =   +862239689D0/2450448D2,
     &  DB12_9_7    =   -2973275233D0/58810752D2,
     &  DB12_9_8    =   +443339993D0/29405376D1,
     &  DB12_9_9    =   +4141950047D0/39207168D2,
     &  DB12_9_10   =   -52868905381D0/22054032D4,
     &  DB12_9_11   =   -411440411D0/352864512D2,
     &  DB12_10_1   =   +80335D0/299376D0,
     &  DB12_10_2   =   +132875D0/74844D0,
     &  DB12_10_3   =   -80875D0/99792D0,
     &  DB12_10_4   =   +28375D0/6237D0,
     &  DB12_10_5   =   -24125D0/5544D0,
     &  DB12_10_6   =   +89035D0/12474D0,
     &  DB12_10_7   =   -24125D0/5544D0,
     &  DB12_10_8   =   +28375D0/6237D0,
     &  DB12_10_9   =   -80875D0/99792D0,
     &  DB12_10_10  =   +132875D0/74844D0,
     &  DB12_10_11  =   -8769811D0/2338875D1,
     &  B212_1_1    =   -5169648083607D0/5D13,
     &  B212_2_1    =   +97697971D0/6103515625D0,
     &  B212_3_1    =   +4558075961257D0/5D13,
     &  B212_4_1    =   -8992156761D0/9765625D4,
     &  B212_5_1    =   -45031D0/512D4,
     &  B212_6_1    =   +1614856691D0/244140625D2,
     &  B212_7_1    =   -1513503946143D0/5D13,
     &  B212_8_1    =   -1333129889D0/9765625D4,
     &  B212_9_1    =   +872716997593D0/5D13,
     &  B212_10_1   =   +0D0)


C ---------------------------------------------------------------------------------------

      goto(10,20) imas+1

10    continue
c     ODE case

c      Z=[y0 Y]*(DA)'
      do i=1,m
         Z(i,1)=y0(i) *DA12_1_1  +Y(i,1) *DA12_1_2  +Y(i,2)*DA12_1_3+
     &          Y(i,3)*DA12_1_4  +Y(i,4) *DA12_1_5  +Y(i,5)*DA12_1_6+
     &          Y(i,6)*DA12_1_7  +Y(i,7) *DA12_1_8  +Y(i,8)*DA12_1_9+
     &          Y(i,9)*DA12_1_10 +Y(i,10)*DA12_1_11
         Z(i,2)=y0(i) *DA12_2_1  +Y(i,1) *DA12_2_2  +Y(i,2)*DA12_2_3+
     &          Y(i,3)*DA12_2_4  +Y(i,4) *DA12_2_5  +Y(i,5)*DA12_2_6+
     &          Y(i,6)*DA12_2_7  +Y(i,7) *DA12_2_8  +Y(i,8)*DA12_2_9+
     &          Y(i,9)*DA12_2_10 +Y(i,10)*DA12_2_11
         Z(i,3)=y0(i) *DA12_3_1  +Y(i,1) *DA12_3_2  +Y(i,2)*DA12_3_3+
     &          Y(i,3)*DA12_3_4  +Y(i,4) *DA12_3_5  +Y(i,5)*DA12_3_6+
     &          Y(i,6)*DA12_3_7  +Y(i,7) *DA12_3_8  +Y(i,8)*DA12_3_9+
     &            Y(i,9)*DA12_3_10 +Y(i,10)*DA12_3_11
         Z(i,4)=y0(i) *DA12_4_1  +Y(i,1) *DA12_4_2  +Y(i,2)*DA12_4_3+
     &          Y(i,3)*DA12_4_4  +Y(i,4) *DA12_4_5  +Y(i,5)*DA12_4_6+
     &          Y(i,6)*DA12_4_7  +Y(i,7) *DA12_4_8  +Y(i,8)*DA12_4_9+
     &          Y(i,9)*DA12_4_10 +Y(i,10)*DA12_4_11
         Z(i,5)=y0(i) *DA12_5_1  +Y(i,1) *DA12_5_2  +Y(i,2)*DA12_5_3+
     &          Y(i,3)*DA12_5_4  +Y(i,4) *DA12_5_5  +Y(i,5)*DA12_5_6+
     &          Y(i,6)*DA12_5_7  +Y(i,7) *DA12_5_8  +Y(i,8)*DA12_5_9+
     &          Y(i,9)*DA12_5_10 +Y(i,10)*DA12_5_11
         Z(i,6)=y0(i) *DA12_6_1  +Y(i,1) *DA12_6_2  +Y(i,2)*DA12_6_3+
     &          Y(i,3)*DA12_6_4  +Y(i,4) *DA12_6_5  +Y(i,5)*DA12_6_6+
     &          Y(i,6)*DA12_6_7  +Y(i,7) *DA12_6_8  +Y(i,8)*DA12_6_9+
     &          Y(i,9)*DA12_6_10 +Y(i,10)*DA12_6_11
         Z(i,7)=y0(i) *DA12_7_1  +Y(i,1) *DA12_7_2  +Y(i,2)*DA12_7_3+
     &          Y(i,3)*DA12_7_4  +Y(i,4) *DA12_7_5  +Y(i,5)*DA12_7_6+
     &          Y(i,6)*DA12_7_7  +Y(i,7) *DA12_7_8  +Y(i,8)*DA12_7_9+
     &          Y(i,9)*DA12_7_10 +Y(i,10)*DA12_7_11
         Z(i,8)=y0(i) *DA12_8_1  +Y(i,1) *DA12_8_2  +Y(i,2)*DA12_8_3+
     &          Y(i,3)*DA12_8_4  +Y(i,4) *DA12_8_5  +Y(i,5)*DA12_8_6+
     &          Y(i,6)*DA12_8_7  +Y(i,7) *DA12_8_8  +Y(i,8)*DA12_8_9+
     &          Y(i,9)*DA12_8_10 +Y(i,10)*DA12_8_11
         Z(i,9)=y0(i) *DA12_9_1  +Y(i,1) *DA12_9_2  +Y(i,2)*DA12_9_3+
     &          Y(i,3)*DA12_9_4  +Y(i,4) *DA12_9_5  +Y(i,5)*DA12_9_6+
     &          Y(i,6)*DA12_9_7  +Y(i,7) *DA12_9_8  +Y(i,8)*DA12_9_9+
     &          Y(i,9)*DA12_9_10 +Y(i,10)*DA12_9_11
         Z(i,10)=y0(i) *DA12_10_1 +Y(i,1) *DA12_10_2 +Y(i,2)*DA12_10_3+
     &         Y(i,3)*DA12_10_4 +Y(i,4) *DA12_10_5 +Y(i,5)*DA12_10_6+
     &         Y(i,6)*DA12_10_7 +Y(i,7) *DA12_10_8 +Y(i,8)*DA12_10_9+
     &         Y(i,9)*DA12_10_10+Y(i,10)*DA12_10_11
      end do

c     Z=Z-h*[f0 F]*(dB)'
      do i=1,m
         Z(i,1) = Z(i,1)-
     &      h*(f0(i) *DB12_1_1  +F(i,1) *DB12_1_2  +F(i,2)*DB12_1_3+
     &         F(i,3)*DB12_1_4  +F(i,4) *DB12_1_5  +F(i,5)*DB12_1_6+
     &         F(i,6)*DB12_1_7  +F(i,7) *DB12_1_8  +F(i,8)*DB12_1_9+
     &         F(i,9)*DB12_1_10 +F(i,10)*DB12_1_11)
         Z(i,2) = Z(i,2)-
     &      h*(f0(i) *DB12_2_1  +F(i,1) *DB12_2_2  +F(i,2)*DB12_2_3+
     &         F(i,3)*DB12_2_4  +F(i,4) *DB12_2_5  +F(i,5)*DB12_2_6+
     &         F(i,6)*DB12_2_7  +F(i,7) *DB12_2_8  +F(i,8)*DB12_2_9+
     &         F(i,9)*DB12_2_10 +F(i,10)*DB12_2_11)
         Z(i,3) = Z(i,3)-
     &      h*(f0(i) *DB12_3_1  +F(i,1) *DB12_3_2  +F(i,2)*DB12_3_3+
     &         F(i,3)*DB12_3_4  +F(i,4) *DB12_3_5  +F(i,5)*DB12_3_6+
     &         F(i,6)*DB12_3_7  +F(i,7) *DB12_3_8  +F(i,8)*DB12_3_9+
     &         F(i,9)*DB12_3_10 +F(i,10)*DB12_3_11)
         Z(i,4) = Z(i,4)-
     &      h*(f0(i) *DB12_4_1  +F(i,1) *DB12_4_2  +F(i,2)*DB12_4_3+
     &         F(i,3)*DB12_4_4  +F(i,4) *DB12_4_5  +F(i,5)*DB12_4_6+
     &         F(i,6)*DB12_4_7  +F(i,7) *DB12_4_8  +F(i,8)*DB12_4_9+
     &         F(i,9)*DB12_4_10 +F(i,10)*DB12_4_11)
         Z(i,5) = Z(i,5)-
     &      h*(f0(i) *DB12_5_1  +F(i,1) *DB12_5_2  +F(i,2)*DB12_5_3+
     &         F(i,3)*DB12_5_4  +F(i,4) *DB12_5_5  +F(i,5)*DB12_5_6+
     &         F(i,6)*DB12_5_7  +F(i,7) *DB12_5_8  +F(i,8)*DB12_5_9+
     &         F(i,9)*DB12_5_10 +F(i,10)*DB12_5_11)
         Z(i,6) = Z(i,6)-
     &      h*(f0(i) *DB12_6_1  +F(i,1) *DB12_6_2  +F(i,2)*DB12_6_3+
     &         F(i,3)*DB12_6_4  +F(i,4) *DB12_6_5  +F(i,5)*DB12_6_6+
     &         F(i,6)*DB12_6_7  +F(i,7) *DB12_6_8  +F(i,8)*DB12_6_9+
     &         F(i,9)*DB12_6_10 +F(i,10)*DB12_6_11)
         Z(i,7) = Z(i,7)-
     &      h*(f0(i) *DB12_7_1  +F(i,1) *DB12_7_2  +F(i,2)*DB12_7_3+
     &         F(i,3)*DB12_7_4  +F(i,4) *DB12_7_5  +F(i,5)*DB12_7_6+
     &         F(i,6)*DB12_7_7  +F(i,7) *DB12_7_8  +F(i,8)*DB12_7_9+
     &         F(i,9)*DB12_7_10 +F(i,10)*DB12_7_11)
         Z(i,8) = Z(i,8)-
     &      h*(f0(i) *DB12_8_1  +F(i,1) *DB12_8_2  +F(i,2)*DB12_8_3+
     &         F(i,3)*DB12_8_4  +F(i,4) *DB12_8_5  +F(i,5)*DB12_8_6+
     &         F(i,6)*DB12_8_7  +F(i,7) *DB12_8_8  +F(i,8)*DB12_8_9+
     &         F(i,9)*DB12_8_10 +F(i,10)*DB12_8_11)
         Z(i,9) = Z(i,9)-
     &      h*(f0(i) *DB12_9_1  +F(i,1) *DB12_9_2  +F(i,2)*DB12_9_3+
     &         F(i,3)*DB12_9_4  +F(i,4) *DB12_9_5  +F(i,5)*DB12_9_6+
     &         F(i,6)*DB12_9_7  +F(i,7) *DB12_9_8  +F(i,8)*DB12_9_9+
     &         F(i,9)*DB12_9_10 +F(i,10)*DB12_9_11)
         Z(i,10) = Z(i,10)-
     &      h*(f0(i) *DB12_10_1 +F(i,1) *DB12_10_2 +F(i,2)*DB12_10_3+
     &         F(i,3)*DB12_10_4 +F(i,4) *DB12_10_5 +F(i,5)*DB12_10_6+
     &         F(i,6)*DB12_10_7 +F(i,7) *DB12_10_8 +F(i,8)*DB12_10_9+
     &         F(i,9)*DB12_10_10+F(i,10)*DB12_10_11)
      end do

c     theta*Z=Z
      call sollu(m,theta,ldlu,Z(1,1),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,2),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,3),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,4),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,5),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,6),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,7),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,8),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,9),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,10),mljac,mujac,ipvt,ijob(1))

c      Z= Z+[y0 Y]*(A2)'
      do i=1,m
            Z(i,1)=Z(i,1)+
     &          y0(i)* A212_1_1  +Y(i,1) *A212_1_2  +Y(i,2)*A212_1_3+
     &         Y(i,3)* A212_1_4  +Y(i,4) *A212_1_5  +Y(i,5)*A212_1_6+
     &         Y(i,6)* A212_1_7  +Y(i,7) *A212_1_8  +Y(i,8)*A212_1_9+
     &         Y(i,9)* A212_1_10 +Y(i,10)*A212_1_11
            Z(i,2)=Z(i,2)+
     &          y0(i)* A212_2_1  +Y(i,1) *A212_2_2  +Y(i,2)*A212_2_3+
     &         Y(i,3)* A212_2_4  +Y(i,4) *A212_2_5  +Y(i,5)*A212_2_6+
     &         Y(i,6)* A212_2_7  +Y(i,7) *A212_2_8  +Y(i,8)*A212_2_9+
     &         Y(i,9)* A212_2_10 +Y(i,10)*A212_2_11
            Z(i,3)=Z(i,3)+
     &          y0(i)* A212_3_1  +Y(i,1) *A212_3_2  +Y(i,2)*A212_3_3+
     &         Y(i,3)* A212_3_4  +Y(i,4) *A212_3_5  +Y(i,5)*A212_3_6+
     &         Y(i,6)* A212_3_7  +Y(i,7) *A212_3_8  +Y(i,8)*A212_3_9+
     &         Y(i,9)* A212_3_10 +Y(i,10)*A212_3_11
            Z(i,4)=Z(i,4)+
     &          y0(i)* A212_4_1  +Y(i,1) *A212_4_2  +Y(i,2)*A212_4_3+
     &         Y(i,3)* A212_4_4  +Y(i,4) *A212_4_5  +Y(i,5)*A212_4_6+
     &         Y(i,6)* A212_4_7  +Y(i,7) *A212_4_8  +Y(i,8)*A212_4_9+
     &         Y(i,9)* A212_4_10 +Y(i,10)*A212_4_11
            Z(i,5)=Z(i,5)+
     &          y0(i)* A212_5_1  +Y(i,1) *A212_5_2  +Y(i,2)*A212_5_3+
     &         Y(i,3)* A212_5_4  +Y(i,4) *A212_5_5  +Y(i,5)*A212_5_6+
     &         Y(i,6)* A212_5_7  +Y(i,7) *A212_5_8  +Y(i,8)*A212_5_9+
     &         Y(i,9)* A212_5_10 +Y(i,10)*A212_5_11
            Z(i,6)=Z(i,6)+
     &          y0(i)* A212_6_1  +Y(i,1) *A212_6_2  +Y(i,2)*A212_6_3+
     &         Y(i,3)* A212_6_4  +Y(i,4) *A212_6_5  +Y(i,5)*A212_6_6+
     &         Y(i,6)* A212_6_7  +Y(i,7) *A212_6_8  +Y(i,8)*A212_6_9+
     &         Y(i,9)* A212_6_10 +Y(i,10)*A212_6_11
            Z(i,7)=Z(i,7)+
     &          y0(i)* A212_7_1  +Y(i,1) *A212_7_2  +Y(i,2)*A212_7_3+
     &         Y(i,3)* A212_7_4  +Y(i,4) *A212_7_5  +Y(i,5)*A212_7_6+
     &         Y(i,6)* A212_7_7  +Y(i,7) *A212_7_8  +Y(i,8)*A212_7_9+
     &         Y(i,9)* A212_7_10 +Y(i,10)*A212_7_11
            Z(i,8)=Z(i,8)+
     &          y0(i)* A212_8_1  +Y(i,1) *A212_8_2  +Y(i,2)*A212_8_3+
     &         Y(i,3)* A212_8_4  +Y(i,4) *A212_8_5  +Y(i,5)*A212_8_6+
     &         Y(i,6)* A212_8_7  +Y(i,7) *A212_8_8  +Y(i,8)*A212_8_9+
     &         Y(i,9)* A212_8_10 +Y(i,10)*A212_8_11
            Z(i,9)=Z(i,9)+
     &          y0(i)* A212_9_1  +Y(i,1) *A212_9_2  +Y(i,2)*A212_9_3+
     &         Y(i,3)* A212_9_4  +Y(i,4) *A212_9_5  +Y(i,5)*A212_9_6+
     &         Y(i,6)* A212_9_7  +Y(i,7) *A212_9_8  +Y(i,8)*A212_9_9+
     &         Y(i,9)* A212_9_10 +Y(i,10)*A212_9_11
            Z(i,10)=Z(i,10)+
     &         y0(i)* A212_10_1 +Y(i,1) *A212_10_2 +Y(i,2)*A212_10_3+
     &        Y(i,3)* A212_10_4 +Y(i,4) *A212_10_5 +Y(i,5)*A212_10_6+
     &        Y(i,6)* A212_10_7 +Y(i,7) *A212_10_8 +Y(i,8)*A212_10_9+
     &        Y(i,9)* A212_10_10+Y(i,10)*A212_10_11
      end do

c       Z=Z-h*[f0 F]*(B2)'
      do i=1,m
          Z(i,1)=Z(i,1) -h*(f0(i)*B212_1_1  + F(i,1)*gamma)
          Z(i,2)=Z(i,2) -h*(f0(i)*B212_2_1  + F(i,2)*gamma)
          Z(i,3)=Z(i,3) -h*(f0(i)*B212_3_1  + F(i,3)*gamma)
          Z(i,4)=Z(i,4) -h*(f0(i)*B212_4_1  + F(i,4)*gamma)
          Z(i,5)=Z(i,5) -h*(f0(i)*B212_5_1  + F(i,5)*gamma)
          Z(i,6)=Z(i,6) -h*(f0(i)*B212_6_1  + F(i,6)*gamma)
          Z(i,7)=Z(i,7) -h*(f0(i)*B212_7_1  + F(i,7)*gamma)
          Z(i,8)=Z(i,8) -h*(f0(i)*B212_8_1  + F(i,8)*gamma)
          Z(i,9)=Z(i,9) -h*(f0(i)*B212_9_1  + F(i,9)*gamma)
          Z(i,10)=Z(i,10)-h*(f0(i)*B212_10_1 + F(i,10)*gamma)
      end do

c     theta*Z=Z
      call sollu(m,theta,ldlu,Z(1,1),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,2),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,3),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,4),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,5),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,6),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,7),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,8),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,9),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,10),mljac,mujac,ipvt,ijob(1))

      goto 30

20    continue
c     DAE case

c      Z=[y0 Y]*(DA)'
      do i=1,m
         Z(i,1)=y0(i) *DA12_1_1  +Y(i,1) *DA12_1_2  +Y(i,2)*DA12_1_3+
     &          Y(i,3)*DA12_1_4  +Y(i,4) *DA12_1_5  +Y(i,5)*DA12_1_6+
     &          Y(i,6)*DA12_1_7  +Y(i,7) *DA12_1_8  +Y(i,8)*DA12_1_9+
     &          Y(i,9)*DA12_1_10 +Y(i,10)*DA12_1_11
         Z(i,2)=y0(i) *DA12_2_1  +Y(i,1) *DA12_2_2  +Y(i,2)*DA12_2_3+
     &          Y(i,3)*DA12_2_4  +Y(i,4) *DA12_2_5  +Y(i,5)*DA12_2_6+
     &          Y(i,6)*DA12_2_7  +Y(i,7) *DA12_2_8  +Y(i,8)*DA12_2_9+
     &          Y(i,9)*DA12_2_10 +Y(i,10)*DA12_2_11
         Z(i,3)=y0(i) *DA12_3_1  +Y(i,1) *DA12_3_2  +Y(i,2)*DA12_3_3+
     &          Y(i,3)*DA12_3_4  +Y(i,4) *DA12_3_5  +Y(i,5)*DA12_3_6+
     &          Y(i,6)*DA12_3_7  +Y(i,7) *DA12_3_8  +Y(i,8)*DA12_3_9+
     &            Y(i,9)*DA12_3_10 +Y(i,10)*DA12_3_11
         Z(i,4)=y0(i) *DA12_4_1  +Y(i,1) *DA12_4_2  +Y(i,2)*DA12_4_3+
     &          Y(i,3)*DA12_4_4  +Y(i,4) *DA12_4_5  +Y(i,5)*DA12_4_6+
     &          Y(i,6)*DA12_4_7  +Y(i,7) *DA12_4_8  +Y(i,8)*DA12_4_9+
     &          Y(i,9)*DA12_4_10 +Y(i,10)*DA12_4_11
         Z(i,5)=y0(i) *DA12_5_1  +Y(i,1) *DA12_5_2  +Y(i,2)*DA12_5_3+
     &          Y(i,3)*DA12_5_4  +Y(i,4) *DA12_5_5  +Y(i,5)*DA12_5_6+
     &          Y(i,6)*DA12_5_7  +Y(i,7) *DA12_5_8  +Y(i,8)*DA12_5_9+
     &          Y(i,9)*DA12_5_10 +Y(i,10)*DA12_5_11
         Z(i,6)=y0(i) *DA12_6_1  +Y(i,1) *DA12_6_2  +Y(i,2)*DA12_6_3+
     &          Y(i,3)*DA12_6_4  +Y(i,4) *DA12_6_5  +Y(i,5)*DA12_6_6+
     &          Y(i,6)*DA12_6_7  +Y(i,7) *DA12_6_8  +Y(i,8)*DA12_6_9+
     &          Y(i,9)*DA12_6_10 +Y(i,10)*DA12_6_11
         Z(i,7)=y0(i) *DA12_7_1  +Y(i,1) *DA12_7_2  +Y(i,2)*DA12_7_3+
     &          Y(i,3)*DA12_7_4  +Y(i,4) *DA12_7_5  +Y(i,5)*DA12_7_6+
     &          Y(i,6)*DA12_7_7  +Y(i,7) *DA12_7_8  +Y(i,8)*DA12_7_9+
     &          Y(i,9)*DA12_7_10 +Y(i,10)*DA12_7_11
         Z(i,8)=y0(i) *DA12_8_1  +Y(i,1) *DA12_8_2  +Y(i,2)*DA12_8_3+
     &          Y(i,3)*DA12_8_4  +Y(i,4) *DA12_8_5  +Y(i,5)*DA12_8_6+
     &          Y(i,6)*DA12_8_7  +Y(i,7) *DA12_8_8  +Y(i,8)*DA12_8_9+
     &          Y(i,9)*DA12_8_10 +Y(i,10)*DA12_8_11
         Z(i,9)=y0(i) *DA12_9_1  +Y(i,1) *DA12_9_2  +Y(i,2)*DA12_9_3+
     &          Y(i,3)*DA12_9_4  +Y(i,4) *DA12_9_5  +Y(i,5)*DA12_9_6+
     &          Y(i,6)*DA12_9_7  +Y(i,7) *DA12_9_8  +Y(i,8)*DA12_9_9+
     &          Y(i,9)*DA12_9_10 +Y(i,10)*DA12_9_11
         Z(i,10)=y0(i) *DA12_10_1 +Y(i,1) *DA12_10_2 +Y(i,2)*DA12_10_3+
     &         Y(i,3)*DA12_10_4 +Y(i,4) *DA12_10_5 +Y(i,5)*DA12_10_6+
     &         Y(i,6)*DA12_10_7 +Y(i,7) *DA12_10_8 +Y(i,8)*DA12_10_9+
     &         Y(i,9)*DA12_10_10+Y(i,10)*DA12_10_11
      end do

c     MZ_i=M0*Z_i, i=1,2,...,8

      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,1),MZ(1,1),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,2),MZ(1,2),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,3),MZ(1,3),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,4),MZ(1,4),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,5),MZ(1,5),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,6),MZ(1,6),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,7),MZ(1,7),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,8),MZ(1,8),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,9),MZ(1,9),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,10),MZ(1,10),ijob(2))

c     MZ=MZ-h*[f0 F]*(dB)'

      do i=1,m
         MZ(i,1) = MZ(i,1)-
     &      h*(f0(i) *DB12_1_1  +F(i,1) *DB12_1_2  +F(i,2)*DB12_1_3+
     &         F(i,3)*DB12_1_4  +F(i,4) *DB12_1_5  +F(i,5)*DB12_1_6+
     &         F(i,6)*DB12_1_7  +F(i,7) *DB12_1_8  +F(i,8)*DB12_1_9+
     &         F(i,9)*DB12_1_10 +F(i,10)*DB12_1_11)
         MZ(i,2) = MZ(i,2)-
     &      h*(f0(i) *DB12_2_1  +F(i,1) *DB12_2_2  +F(i,2)*DB12_2_3+
     &         F(i,3)*DB12_2_4  +F(i,4) *DB12_2_5  +F(i,5)*DB12_2_6+
     &         F(i,6)*DB12_2_7  +F(i,7) *DB12_2_8  +F(i,8)*DB12_2_9+
     &         F(i,9)*DB12_2_10 +F(i,10)*DB12_2_11)
         MZ(i,3) = MZ(i,3)-
     &      h*(f0(i) *DB12_3_1  +F(i,1) *DB12_3_2  +F(i,2)*DB12_3_3+
     &         F(i,3)*DB12_3_4  +F(i,4) *DB12_3_5  +F(i,5)*DB12_3_6+
     &         F(i,6)*DB12_3_7  +F(i,7) *DB12_3_8  +F(i,8)*DB12_3_9+
     &         F(i,9)*DB12_3_10 +F(i,10)*DB12_3_11)
         MZ(i,4) = MZ(i,4)-
     &      h*(f0(i) *DB12_4_1  +F(i,1) *DB12_4_2  +F(i,2)*DB12_4_3+
     &         F(i,3)*DB12_4_4  +F(i,4) *DB12_4_5  +F(i,5)*DB12_4_6+
     &         F(i,6)*DB12_4_7  +F(i,7) *DB12_4_8  +F(i,8)*DB12_4_9+
     &         F(i,9)*DB12_4_10 +F(i,10)*DB12_4_11)
         MZ(i,5) = MZ(i,5)-
     &      h*(f0(i) *DB12_5_1  +F(i,1) *DB12_5_2  +F(i,2)*DB12_5_3+
     &         F(i,3)*DB12_5_4  +F(i,4) *DB12_5_5  +F(i,5)*DB12_5_6+
     &         F(i,6)*DB12_5_7  +F(i,7) *DB12_5_8  +F(i,8)*DB12_5_9+
     &         F(i,9)*DB12_5_10 +F(i,10)*DB12_5_11)
         MZ(i,6) = MZ(i,6)-
     &      h*(f0(i) *DB12_6_1  +F(i,1) *DB12_6_2  +F(i,2)*DB12_6_3+
     &         F(i,3)*DB12_6_4  +F(i,4) *DB12_6_5  +F(i,5)*DB12_6_6+
     &         F(i,6)*DB12_6_7  +F(i,7) *DB12_6_8  +F(i,8)*DB12_6_9+
     &         F(i,9)*DB12_6_10 +F(i,10)*DB12_6_11)
         MZ(i,7) = MZ(i,7)-
     &      h*(f0(i) *DB12_7_1  +F(i,1) *DB12_7_2  +F(i,2)*DB12_7_3+
     &         F(i,3)*DB12_7_4  +F(i,4) *DB12_7_5  +F(i,5)*DB12_7_6+
     &         F(i,6)*DB12_7_7  +F(i,7) *DB12_7_8  +F(i,8)*DB12_7_9+
     &         F(i,9)*DB12_7_10 +F(i,10)*DB12_7_11)
         MZ(i,8) = MZ(i,8)-
     &      h*(f0(i) *DB12_8_1  +F(i,1) *DB12_8_2  +F(i,2)*DB12_8_3+
     &         F(i,3)*DB12_8_4  +F(i,4) *DB12_8_5  +F(i,5)*DB12_8_6+
     &         F(i,6)*DB12_8_7  +F(i,7) *DB12_8_8  +F(i,8)*DB12_8_9+
     &         F(i,9)*DB12_8_10 +F(i,10)*DB12_8_11)
         MZ(i,9) = MZ(i,9)-
     &      h*(f0(i) *DB12_9_1  +F(i,1) *DB12_9_2  +F(i,2)*DB12_9_3+
     &         F(i,3)*DB12_9_4  +F(i,4) *DB12_9_5  +F(i,5)*DB12_9_6+
     &         F(i,6)*DB12_9_7  +F(i,7) *DB12_9_8  +F(i,8)*DB12_9_9+
     &         F(i,9)*DB12_9_10 +F(i,10)*DB12_9_11)
         MZ(i,10) = MZ(i,10)-
     &      h*(f0(i) *DB12_10_1 +F(i,1) *DB12_10_2 +F(i,2)*DB12_10_3+
     &         F(i,3)*DB12_10_4 +F(i,4) *DB12_10_5 +F(i,5)*DB12_10_6+
     &         F(i,6)*DB12_10_7 +F(i,7) *DB12_10_8 +F(i,8)*DB12_10_9+
     &         F(i,9)*DB12_10_10+F(i,10)*DB12_10_11)
      end do

c     theta*MZ=MZ
      call sollu(m,theta,ldlu,MZ(1,1),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,2),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,3),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,4),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,5),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,6),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,7),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,8),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,9),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,MZ(1,10),mljac,mujac,ipvt,ijob(1))

c      MZ= MZ+[y0 Y]*(A2)'
      do i=1,m
            MZ(i,1)=MZ(i,1)+
     &          y0(i)* A212_1_1  +Y(i,1) *A212_1_2  +Y(i,2)*A212_1_3+
     &         Y(i,3)* A212_1_4  +Y(i,4) *A212_1_5  +Y(i,5)*A212_1_6+
     &         Y(i,6)* A212_1_7  +Y(i,7) *A212_1_8  +Y(i,8)*A212_1_9+
     &         Y(i,9)* A212_1_10 +Y(i,10)*A212_1_11
            MZ(i,2)=MZ(i,2)+
     &          y0(i)* A212_2_1  +Y(i,1) *A212_2_2  +Y(i,2)*A212_2_3+
     &         Y(i,3)* A212_2_4  +Y(i,4) *A212_2_5  +Y(i,5)*A212_2_6+
     &         Y(i,6)* A212_2_7  +Y(i,7) *A212_2_8  +Y(i,8)*A212_2_9+
     &         Y(i,9)* A212_2_10 +Y(i,10)*A212_2_11
            MZ(i,3)=MZ(i,3)+
     &          y0(i)* A212_3_1  +Y(i,1) *A212_3_2  +Y(i,2)*A212_3_3+
     &         Y(i,3)* A212_3_4  +Y(i,4) *A212_3_5  +Y(i,5)*A212_3_6+
     &         Y(i,6)* A212_3_7  +Y(i,7) *A212_3_8  +Y(i,8)*A212_3_9+
     &         Y(i,9)* A212_3_10 +Y(i,10)*A212_3_11
            MZ(i,4)=MZ(i,4)+
     &          y0(i)* A212_4_1  +Y(i,1) *A212_4_2  +Y(i,2)*A212_4_3+
     &         Y(i,3)* A212_4_4  +Y(i,4) *A212_4_5  +Y(i,5)*A212_4_6+
     &         Y(i,6)* A212_4_7  +Y(i,7) *A212_4_8  +Y(i,8)*A212_4_9+
     &         Y(i,9)* A212_4_10 +Y(i,10)*A212_4_11
            MZ(i,5)=MZ(i,5)+
     &          y0(i)* A212_5_1  +Y(i,1) *A212_5_2  +Y(i,2)*A212_5_3+
     &         Y(i,3)* A212_5_4  +Y(i,4) *A212_5_5  +Y(i,5)*A212_5_6+
     &         Y(i,6)* A212_5_7  +Y(i,7) *A212_5_8  +Y(i,8)*A212_5_9+
     &         Y(i,9)* A212_5_10 +Y(i,10)*A212_5_11
            MZ(i,6)=MZ(i,6)+
     &          y0(i)* A212_6_1  +Y(i,1) *A212_6_2  +Y(i,2)*A212_6_3+
     &         Y(i,3)* A212_6_4  +Y(i,4) *A212_6_5  +Y(i,5)*A212_6_6+
     &         Y(i,6)* A212_6_7  +Y(i,7) *A212_6_8  +Y(i,8)*A212_6_9+
     &         Y(i,9)* A212_6_10 +Y(i,10)*A212_6_11
            MZ(i,7)=MZ(i,7)+
     &          y0(i)* A212_7_1  +Y(i,1) *A212_7_2  +Y(i,2)*A212_7_3+
     &         Y(i,3)* A212_7_4  +Y(i,4) *A212_7_5  +Y(i,5)*A212_7_6+
     &         Y(i,6)* A212_7_7  +Y(i,7) *A212_7_8  +Y(i,8)*A212_7_9+
     &         Y(i,9)* A212_7_10 +Y(i,10)*A212_7_11
            MZ(i,8)=MZ(i,8)+
     &          y0(i)* A212_8_1  +Y(i,1) *A212_8_2  +Y(i,2)*A212_8_3+
     &         Y(i,3)* A212_8_4  +Y(i,4) *A212_8_5  +Y(i,5)*A212_8_6+
     &         Y(i,6)* A212_8_7  +Y(i,7) *A212_8_8  +Y(i,8)*A212_8_9+
     &         Y(i,9)* A212_8_10 +Y(i,10)*A212_8_11
            MZ(i,9)=MZ(i,9)+
     &          y0(i)* A212_9_1  +Y(i,1) *A212_9_2  +Y(i,2)*A212_9_3+
     &         Y(i,3)* A212_9_4  +Y(i,4) *A212_9_5  +Y(i,5)*A212_9_6+
     &         Y(i,6)* A212_9_7  +Y(i,7) *A212_9_8  +Y(i,8)*A212_9_9+
     &         Y(i,9)* A212_9_10 +Y(i,10)*A212_9_11
            MZ(i,10)=MZ(i,10)+
     &         y0(i)* A212_10_1 +Y(i,1) *A212_10_2 +Y(i,2)*A212_10_3+
     &        Y(i,3)* A212_10_4 +Y(i,4) *A212_10_5 +Y(i,5)*A212_10_6+
     &        Y(i,6)* A212_10_7 +Y(i,7) *A212_10_8 +Y(i,8)*A212_10_9+
     &        Y(i,9)* A212_10_10+Y(i,10)*A212_10_11
      end do

c     Z_i=M0*MZ_i, i=1,2,...,8

      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,1),Z(1,1),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,2),Z(1,2),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,3),Z(1,3),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,4),Z(1,4),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,5),Z(1,5),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,6),Z(1,6),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,7),Z(1,7),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,8),Z(1,8),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,9),Z(1,9),ijob(2))
      call MATVEC0(m,M0,ldmas,mlmas,mumas,MZ(1,10),Z(1,10),ijob(2))

c     Z=Z-h*[f0 F]*(B2)'
      do i=1,m
         Z(i,1)=Z(i,1) -h*(f0(i)*B212_1_1  + F(i,1)*gamma)
         Z(i,2)=Z(i,2) -h*(f0(i)*B212_2_1  + F(i,2)*gamma)
         Z(i,3)=Z(i,3) -h*(f0(i)*B212_3_1  + F(i,3)*gamma)
         Z(i,4)=Z(i,4) -h*(f0(i)*B212_4_1  + F(i,4)*gamma)
         Z(i,5)=Z(i,5) -h*(f0(i)*B212_5_1  + F(i,5)*gamma)
         Z(i,6)=Z(i,6) -h*(f0(i)*B212_6_1  + F(i,6)*gamma)
         Z(i,7)=Z(i,7) -h*(f0(i)*B212_7_1  + F(i,7)*gamma)
         Z(i,8)=Z(i,8) -h*(f0(i)*B212_8_1  + F(i,8)*gamma)
         Z(i,9)=Z(i,9) -h*(f0(i)*B212_9_1  + F(i,9)*gamma)
         Z(i,10)=Z(i,10)-h*(f0(i)*B212_10_1 + F(i,10)*gamma)
      end do

c     theta*Z=Z
      call sollu(m,theta,ldlu,Z(1,1),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,2),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,3),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,4),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,5),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,6),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,7),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,8),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,9),mljac,mujac,ipvt,ijob(1))
      call sollu(m,theta,ldlu,Z(1,10),mljac,mujac,ipvt,ijob(1))

30    continue

c
c     Y = Y-Z
c
      do j=1,k
         do i=1,m
           Y(i,j) = Y(i,j) - Z(i,j)
         end do
      end do

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     LOCALERR, ERRDOWN AND ERRUP
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine truncam(m,k,f0,F,h,Z,ord_ind)

c     Truncation Error estimate

      implicit none

c     Input parameters
      integer m,k,ord_ind
      double precision f0(m),F(m,k),h

c     Output parameter
      double precision Z(m)

c     Local variables
      integer i
      double precision PSI4_1,PSI4_2,
     &                 PSI6_1,PSI6_2,PSI6_3,
     &                 PSI8_1,PSI8_2,PSI8_3,PSI8_4,
     &                 PSI10_1,PSI10_2,PSI10_3,PSI10_4,PSI10_5,
     &                 PSI12_1,PSI12_2,PSI12_3,PSI12_4,PSI12_5,PSI12_6

      parameter( PSI4_1 = -1D0,
     &           PSI4_2 =  3D0,
     &           PSI6_1 =  1D0,
     &           PSI6_2 = -4D0,
     &           PSI6_3 =  6D0,
     &           PSI8_1 =  1D0,
     &           PSI8_2 = -6D0,
     &           PSI8_3 =  15D0,
     &           PSI8_4 = -20D0,
     &           PSI10_1 =  1D0,
     &           PSI10_2 = -8D0,
     &           PSI10_3 = 28D0,
     &           PSI10_4 =-56D0,
     &           PSI10_5 = 70D0,
     &           PSI12_1 =  1D0,
     &           PSI12_2 = -10D0,
     &           PSI12_3 =  45D0,
     &           PSI12_4 = -12D1,
     &           PSI12_5 =  21D1,
     &           PSI12_6 =-252D0)

      goto (10,20,30,40,50) ord_ind

10    continue
      do i=1,m
           Z(i)=h*(PSI4_1*f0(i) +PSI4_2*F(i,1)-
     &             PSI4_2*F(i,2)-PSI4_1*F(i,3))
      end do
      return

20    continue
      do i=1,m
           Z(i)=h*(PSI6_1*f0(i) +PSI6_2*F(i,1)+
     &             PSI6_3*F(i,2)+PSI6_2*F(i,3)+
     &             PSI6_1*F(i,4))
      end do
      return

30    continue
      do i=1,m
           Z(i)=h*(PSI8_1*f0(i) +PSI8_2*F(i,1)+
     &             PSI8_3*F(i,2)+PSI8_4*F(i,3)+
     &             PSI8_3*F(i,4)+PSI8_2*F(i,5)+
     &             PSI8_1*F(i,6))
      end do
      return

40    continue
      do i=1,m
           Z(i)=h*(PSI10_1*f0(i) +PSI10_2*F(i,1)+
     &             PSI10_3*F(i,2)+PSI10_4*F(i,3)+
     &             PSI10_5*F(i,4)+PSI10_4*F(i,5)+
     &             PSI10_3*F(i,6)+PSI10_2*F(i,7)+
     &             PSI10_1*F(i,8))
      end do
      return

50    continue
      do i=1,m
           Z(i)=h*(PSI12_1*f0(i) +PSI12_2*F(i,1)+
     &             PSI12_3*F(i,2)+PSI12_4*F(i,3)+
     &             PSI12_5*F(i,4)+PSI12_6*F(i,5)+
     &             PSI12_5*F(i,6)+PSI12_4*F(i,7)+
     &             PSI12_3*F(i,8)+PSI12_2*F(i,9)+
     &             PSI12_1*F(i,10))
      end do
      return

      end

C -------------------------------------------------------------------

      subroutine localerr4(m,f0,F,h,Z,scal,nerr,nerrup,nlinsys,
     &                     theta,vmax,ipvt,ldlu,mljac,mujac,ijob,
     &                     imas,ldmas,mlmas,mumas,M0,k,ord_ind,
     &                     index1,index2)
c
c     Local error estimate for the method of order 4
c
      implicit none

c     Input parameters
      integer m,ipvt(m),ldlu,mljac,mujac,ijob(2),imas,ldmas,mlmas,mumas,
     &        index1,index2,ord_ind,k
      double precision f0(m),F(m,3),theta(ldlu,m),h,scal(m),M0(ldmas,m),
     &                 vmax(3)

c     Output parameters
      double precision Z(m,3),nerr,nerrup

c     I/O parameters
      integer nlinsys

c     Local variables
      integer i

      call truncam(m,k,f0,F,h,Z,ord_ind)

      do i=1,m
         Z(i,2) = Z(i,1)
      end do

      call sollu(m,theta,ldlu,Z(1,2),mljac,mujac,ipvt,ijob(1))

      goto (10,20) imas+1

10    continue
c     ODE case

      do i=1,m
         Z(i,3)=Z(i,1)-Z(i,2)
      end do

      call sollu(m,theta,ldlu,Z(1,3),mljac,mujac,ipvt,ijob(1))

      do i=1,m
           Z(i,2)=vmax(1)*Z(i,2)
           Z(i,3)=vmax(2)*Z(i,3)
      end do

      call norm(m,2,scal,Z(1,2),nerr,nerrup)

      nlinsys = nlinsys + 2

      return

20    continue
c     DAE case

      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,2),Z(1,3),ijob(2))

      do i=1,m
          Z(i,3)=Z(i,1)-Z(i,3)
      end do

      call sollu(m,theta,ldlu,Z(1,3),mljac,mujac,ipvt,ijob(1))

      do i=1,index1
           Z(i,2)=vmax(1)*Z(i,2)
           Z(i,3)=vmax(2)*Z(i,3)
      end do
      do i =index1+1,index1+index2
           Z(i,2) = vmax(2)*Z(i,2)
           Z(i,3) = vmax(2)*Z(i,3)
      end do
      do i = index1+index2+1,m
           Z(i,2) = vmax(3)*Z(i,2)
           Z(i,3) = vmax(3)*Z(i,3)/2d0
      end do

      call norm(m,2,scal,Z(1,2),nerr,nerrup)

      nlinsys = nlinsys + 2

      return

      end

C -------------------------------------------------------------------

      subroutine localerr(m,f0,F,h,Z,scal,nerr,nerrup,nlinsys,
     &                    theta,vmax,ipvt,ldlu,mljac,mujac,ijob,
     &                    imas,ldmas,mlmas,mumas,M0,k,ord_ind,
     &                    index1,index2)
c
c     Local error estimate for the method of order 6
c
      implicit none

c     Input parameters
      integer m,ipvt(m),ldlu,mljac,mujac,ijob(2),imas,ldmas,mlmas,mumas,
     &        index1,index2,ord_ind,k
      double precision f0(m),F(m,k),theta(ldlu,m),h,scal(m),M0(ldmas,m),
     &                 vmax(3)

c     Output parameters
      double precision Z(m,4),nerr,nerrup

c     I/O parameters
      integer nlinsys

c     Local variables
      integer i

      call truncam(m,k,f0,F,h,Z,ord_ind)

      do i=1,m
         Z(i,2) = Z(i,1)
      end do

      call sollu(m,theta,ldlu,Z(1,2),mljac,mujac,ipvt,ijob(1))

      goto (10,20) imas+1

10    continue
c     ODE case

      do i=1,m
         Z(i,3)=2d0*Z(i,1)-Z(i,2)
      end do

      call sollu(m,theta,ldlu,Z(1,3),mljac,mujac,ipvt,ijob(1))

      do i=1,m
         Z(i,3)=Z(i,1)-Z(i,3)
      end do

      call sollu(m,theta,ldlu,Z(1,3),mljac,mujac,ipvt,ijob(1))

      do i=1,m
           Z(i,2)=vmax(1)*Z(i,2)
           Z(i,3)=vmax(2)*Z(i,3)
      end do

      call norm(m,2,scal,Z(1,2),nerr,nerrup)

      nlinsys = nlinsys + 3

      return

20    continue
c     DAE case

      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,2),Z(1,3),ijob(2))

      do i=1,m
         Z(i,3)=2d0*Z(i,1)-Z(i,3)
      end do

      call sollu(m,theta,ldlu,Z(1,3),mljac,mujac,ipvt,ijob(1))

      call MATVEC0(m,M0,ldmas,mlmas,mumas,Z(1,3),Z(1,4),ijob(2))

      do i=1,m
         Z(i,3)=Z(i,1)-Z(i,4)
      end do

      call sollu(m,theta,ldlu,Z(1,3),mljac,mujac,ipvt,ijob(1))

      do i=1,index1
           Z(i,2)=vmax(1)*Z(i,2)
           Z(i,3)=vmax(2)*Z(i,3)
      end do
      do i=index1+1,index1+index2
           Z(i,2)=vmax(2)*Z(i,2)
           Z(i,3)=vmax(2)*Z(i,3)
      end do
      do i=index1+index2+1,m
           Z(i,2)=vmax(3)*Z(i,2)
           Z(i,3)=vmax(3)*Z(i,3)/3d0
      end do
      call norm(m,2,scal,Z(1,2),nerr,nerrup)

      return
      end

C -------------------------------------------------------------------

      subroutine errdown(m,f0,F,h,Z,scal,nerrdown,nlinsys,vmax,
     &                   qinf,theta,ipvt,ldlu,mljac,mujac,ijob,
     &                   k,ord_ind,index1,index2)
c
c     Estimate of the error for the lower-order method,  k>=4
c
      implicit none

c     Input parameters
      integer m,ipvt(m),ldlu,mljac,mujac,ijob,index1,index2,k,ord_ind
      double precision f0(m),F(m,k-1),theta(ldlu,m),h,scal(m),vmax(3)
      logical qinf

c     Output parameters
      double precision Z(m),nerrdown

c     I/O parameters
      integer nlinsys

c     Local variables
      integer i
      double precision vmaxl

      call truncam(m,k-1,f0,F,h,Z,ord_ind-1)

      call sollu(m,theta,ldlu,Z,mljac,mujac,ipvt,ijob)

      if (qinf) then
         vmaxl=vmax(2)
      else
         vmaxl=vmax(1)
      end if

      do i=1,index1
           Z(i)=vmaxl*Z(i)
      end do
      do i =index1+1,index1+index2
           Z(i)=vmax(2)*Z(i)
      end do
      do i =index1+index2+1,m
           Z(i)=vmax(3)*Z(i)
      end do
      call norm(m,1,scal,Z,nerrdown,vmaxl)

      nlinsys = nlinsys + 1

      return

      end

C ---------------------- SUBROUTINE ERRUP ------------------------------

      subroutine errup(m,k,ord_ind,Z,h,h0,h00,vmax,nerrup,scal,
     &                 theta,ipvt,ldlu,mljac,mujac,ijob,
     &                 index1,index2)
c
c     Estimate of the higher-order method, in case of order reduction
c
      implicit none
c     Input parameter
      integer m,k,ipvt(m),ldlu,mljac,mujac,ijob,ord_ind,index1,index2
      double precision Z(m,k+2),theta(ldlu,m),scal(m),
     &                 h,h0,h00,vmax(2)

c     Output parameter
      double precision nerrup

c     Local variables
      integer i
      double precision dbk,dh,dh0,cp,cp0,cp1,cp2

      goto (10,20,20,20,20) ord_ind

10    continue
      cp  = 2d0/dble(k)*h/(h+h0)
      dbk = dble(k+1)
      dh  = (h/h0)**dbk
      do i=1,m
        Z(i,2) = cp *( Z(i,1) - Z(i,k+1)*dh )
      end do
      goto 30

20    continue
      dbk = dble(k+1)
      dh  = (h/h0)**dbk
      dh0 = (h/h00)**dbk
      cp0 = h0+h00
      cp2 = h+h0
      cp1 = cp0+cp2
      cp = cp0*cp1*cp2*dble(k*k)
      cp = 8d0*h*h/cp
      do i=1,m
         Z(i,2)= cp *(  cp0*Z(i,1)- cp1*Z(i,k+1)*dh +
     &                  cp2*Z(i,k+2)*dh0 )
      end do

30    continue
      call sollu(m,theta,ldlu,Z(1,2),mljac,mujac,ipvt,ijob)

      do i=1,index1+index2
         Z(i,2)=vmax(1)*Z(i,2)
      end do
      do i=index1+index2+1,m
         Z(i,2)=vmax(2)*Z(i,2)
      end do

      call norm(m,1,scal,Z(1,2),nerrup,dh)

      return

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     EXTRAPOLA, DIFFDIV, CONTSOL  AND NORM
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE EXTRAPOLA (M,K,KNEW,H0,H,YNEW,DD)
C
C     Initial profile by extrapolation.
C

      IMPLICIT NONE

C     Input parameters
      INTEGER M,K,KNEW
      DOUBLE PRECISION H0,H,DD(K+1,M)

C     Output parameters
      DOUBLE PRECISION YNEW(M,KNEW)

C     Local variables
      INTEGER I,J,L
      DOUBLE PRECISION DT,RATH

      RATH = (H/H0)

      DO I=1,M
C     Evaluation of the interpolating polynomial at the points of the new block
        DO L=1,KNEW
           DT  = dble(L)*RATH
           YNEW(I,L)=DD(K+1,I)
           DO J=K,1,-1
              DT = DT+1d0
              YNEW(I,L)=YNEW(I,L)*DT +DD(J,I)
           END DO
        END DO
      END DO

      RETURN
      END

      SUBROUTINE DIFFDIV(M,K,Y0,Y,DD)
C     Compute the divided differences of the interpolating polynomial

      IMPLICIT NONE

C     Input parameters
      INTEGER M,K
      DOUBLE PRECISION Y0(M),Y(M,K)

C     Output parameters
      DOUBLE PRECISION DD(K+1,M)

C     Local variables
      INTEGER I,J,L
      DOUBLE PRECISION DT

      DO I=1,M
          DD(1,I)=Y0(I)
          DO J=1,K
            DD(J+1,I)=Y(I,J)
          END DO
          DO J=2,K+1
            DT = dble(J-1)
            DO L=K+1,J,-1
               DD(L,I)=(DD(L,I)-DD(L-1,I))/DT
            END DO
          END DO
      END DO
      RETURN
      END


C Karline: this changed to a subroutine - directely all variables estimated
      SUBROUTINE contsolall(T,M,K,T0,TSTEP,DD, YOUT)
c
c     Function to be used if continuous output is desired.
c     It provides the value, at time T, of the polynomial 
c     interpolating the L-th component of the numerical solution 
c     obtained at the last successfully computed step

      IMPLICIT NONE

      INTEGER L,K,M
      DOUBLE PRECISION T,YC,YOUT(M),T0,TSTEP(K),DD(K+1,M)
      INTEGER I
      DOUBLE PRECISION DT      

      IF (TSTEP(1).LE.T0) THEN
        CALL Rprint('Warning in calling subroutine CONTSOL the input')
        CALL Rprint('par. T_0 must be strictly lower than TSTEP(1)')
        CALL Rprintd1('The solution is not returned at t = ', T)
         RETURN
      END IF
      
      DO L = 1, M
        DT = (T-TSTEP(K))/(TSTEP(1)-T0)
        YC = DD(K+1,L)

        DO I=K,1,-1
         DT = DT + 1D0
         YC = YC*DT + DD(I,L)
        END DO
      YOUT(L) = YC 
      ENDDO

      RETURN
      END


      SUBROUTINE NORM(M,K,SCAL,ERR,NERR,NERRUP)
C
C     Used norm
C
      IMPLICIT NONE

C     INPUT PARAMETER
      INTEGER M,K
      DOUBLE PRECISION SCAL(M),ERR(M,K)

C     OUTPUT PARAMETER
      DOUBLE PRECISION NERR,NERRUP

C     LOCAL VARIABLES
      DOUBLE PRECISION NERR0
      INTEGER I,J

      NERR = 0d0
      DO J=1,K-1
         NERR0=0d0
         DO I=1,M
           NERR0   = NERR0 + (ERR(I,J)*SCAL(I))**2
         END DO
         NERR = DMAX1(NERR,NERR0)
      END DO

      NERRUP=0d0
      DO I=1,M
         NERRUP  = NERRUP + (ERR(I,K)*SCAL(I))**2
      END DO
      NERR = DMAX1(NERR,NERRUP)

      NERR   = DSQRT(NERR/DBLE(M))
      NERRUP = DSQRT(NERRUP/DBLE(M))
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     LINEAR ALGEBRA
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C
C     SUBROUTINE DEC
C
      SUBROUTINE DEC (N, NDIM, A, IP, IER)
C VERSION REAL DOUBLE PRECISION
      INTEGER N,NDIM,IP,IER,NM1,K,KP1,M,I,J
      DOUBLE PRECISION A,T
      DIMENSION A(NDIM,N), IP(N)
C-----------------------------------------------------------------------
C  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION.
C  INPUT..
C     N = ORDER OF MATRIX.
C     NDIM = DECLARED DIMENSION OF ARRAY  A .
C     A = MATRIX TO BE TRIANGULARIZED.
C  OUTPUT..
C     A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U .
C     A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
C     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
C     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O .
C     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE
C           SINGULAR AT STAGE K.
C  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
C  DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N).
C  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO.
C
C  REFERENCE..
C     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
C     C.A.C.M. 15 (1972), P. 274.
C-----------------------------------------------------------------------
      IER = 0
      IP(N) = 1
      IF (N .EQ. 1) GO TO 70
      NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = K
        DO 10 I = KP1,N
          IF (DABS(A(I,K)) .GT. DABS(A(M,K))) M = I
 10     CONTINUE
        IP(K) = M
        T = A(M,K)
        IF (M .EQ. K) GO TO 20
        IP(N) = -IP(N)
        A(M,K) = A(K,K)
        A(K,K) = T
 20     CONTINUE
        IF (T .EQ. 0.D0) GO TO 80
        T = 1.D0/T
        DO 30 I = KP1,N
 30       A(I,K) = -A(I,K)*T
        DO 50 J = KP1,N
          T = A(M,J)
          A(M,J) = A(K,J)
          A(K,J) = T
          IF (T .EQ. 0.D0) GO TO 45
          DO 40 I = KP1,N
 40         A(I,J) = A(I,J) + A(I,K)*T
 45       CONTINUE
 50       CONTINUE
 60     CONTINUE
 70   K = N
      IF (A(N,N) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
C----------------------- END OF SUBROUTINE DEC -------------------------
      END

C
C     SUBROUTINE SOL
C
      SUBROUTINE SOL (N, NDIM, A, B, IP)
C VERSION REAL DOUBLE PRECISION
      INTEGER N,NDIM,IP,NM1,K,KP1,M,I,KB,KM1
      DOUBLE PRECISION A,B,T
      DIMENSION A(NDIM,N), B(N), IP(N)
C-----------------------------------------------------------------------
C  SOLUTION OF LINEAR SYSTEM, A*X = B .
C  INPUT..
C    N = ORDER OF MATRIX.
C    NDIM = DECLARED DIMENSION OF ARRAY  A .
C    A = TRIANGULARIZED MATRIX OBTAINED FROM DEC.
C    B = RIGHT HAND SIDE VECTOR.
C    IP = PIVOT VECTOR OBTAINED FROM DEC.
C  DO NOT USE IF DEC HAS SET IER .NE. 0.
C  OUTPUT..
C    B = SOLUTION VECTOR, X .
C-----------------------------------------------------------------------
      IF (N .EQ. 1) GO TO 50
      NM1 = N - 1
      DO 20 K = 1,NM1
        KP1 = K + 1
        M = IP(K)
        T = B(M)
        B(M) = B(K)
        B(K) = T
        DO 10 I = KP1,N
 10       B(I) = B(I) + A(I,K)*T
 20     CONTINUE
      DO 40 KB = 1,NM1
        KM1 = N - KB
        K = KM1 + 1
        B(K) = B(K)/A(K,K)
        T = -B(K)
        DO 30 I = 1,KM1
 30       B(I) = B(I) + A(I,K)*T
 40     CONTINUE
 50   B(1) = B(1)/A(1,1)
      RETURN
C----------------------- END OF SUBROUTINE SOL -------------------------
      END

      subroutine sollu(n,a,lda,b,ml,mu,ipvt,ijob)
      integer lda,n,ipvt(n),ijob
      double precision a(lda,n),b(n)

      goto (1,2) ijob

1     call sol(n,n,a,b,ipvt)
      return

2     call solb(n,lda,a,ml,mu,b,ipvt)
      return
      end

      subroutine declu(n,a,lda,ml,mu,ipvt,ijob,info)
      integer lda,n,ipvt(n),info,ijob
      double precision a(lda,n)

      goto(1,2) ijob

   1  call dec(n,n,a,ipvt,info)
      return

   2  call decb(n,lda,a,ml,mu,ipvt,info)
      return

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccca
c
C     SUBROUTINE DECB
C
      SUBROUTINE DECB (N, NDIM, A, ML, MU, IP, IER)
      REAL*8 A,T
      DIMENSION A(NDIM,N), IP(N)
C-----------------------------------------------------------------------
C  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A BANDED
C  MATRIX WITH LOWER BANDWIDTH ML AND UPPER BANDWIDTH MU
C  INPUT..
C     N       ORDER OF THE ORIGINAL MATRIX A.
C     NDIM    DECLARED DIMENSION OF ARRAY  A.
C     A       CONTAINS THE MATRIX IN BAND STORAGE.   THE COLUMNS
C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  A  AND
C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
C                ML+1 THROUGH 2*ML+MU+1 OF  A.
C     ML      LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C     MU      UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C  OUTPUT..
C     A       AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
C     IP      INDEX VECTOR OF PIVOT INDICES.
C     IP(N)   (-1)**(NUMBER OF INTERCHANGES) OR O .
C     IER     = 0 IF MATRIX A IS NONSINGULAR, OR  = K IF FOUND TO BE
C                SINGULAR AT STAGE K.
C  USE  SOLB  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
C  DETERM(A) = IP(N)*A(MD,1)*A(MD,2)*...*A(MD,N)  WITH MD=ML+MU+1.
C  IF IP(N)=O, A IS SINGULAR, SOLB WILL DIVIDE BY ZERO.
C
C  REFERENCE..
C     THIS IS A MODIFICATION OF
C     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
C     C.A.C.M. 15 (1972), P. 274.
C-----------------------------------------------------------------------
      IER = 0
      IP(N) = 1
      MD = ML + MU + 1
      MD1 = MD + 1
      JU = 0
      IF (ML .EQ. 0) GO TO 70
      IF (N .EQ. 1) GO TO 70
      IF (N .LT. MU+2) GO TO 7
      DO 5 J = MU+2,N
      DO 5 I = 1,ML
  5   A(I,J) = 0.D0
  7   NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = MD
        MDL = MIN(ML,N-K) + MD
        DO 10 I = MD1,MDL
          IF (DABS(A(I,K)) .GT. DABS(A(M,K))) M = I
 10     CONTINUE
        IP(K) = M + K - MD
        T = A(M,K)
        IF (M .EQ. MD) GO TO 20
        IP(N) = -IP(N)
        A(M,K) = A(MD,K)
        A(MD,K) = T
 20     CONTINUE
        IF (T .EQ. 0.D0) GO TO 80
        T = 1.D0/T
        DO 30 I = MD1,MDL
 30       A(I,K) = -A(I,K)*T
        JU = MIN0(MAX0(JU,MU+IP(K)),N)
        MM = MD
        IF (JU .LT. KP1) GO TO 55
        DO 50 J = KP1,JU
          M = M - 1
          MM = MM - 1
          T = A(M,J)
          IF (M .EQ. MM) GO TO 35
          A(M,J) = A(MM,J)
          A(MM,J) = T
 35       CONTINUE
          IF (T .EQ. 0.D0) GO TO 45
          JK = J - K
          DO 40 I = MD1,MDL
            IJK = I - JK
 40         A(IJK,J) = A(IJK,J) + A(I,K)*T
 45       CONTINUE
 50       CONTINUE
 55     CONTINUE
 60     CONTINUE
 70   K = N
      IF (A(MD,N) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
C----------------------- END OF SUBROUTINE DECB ------------------------
      END

C
C     SUBROUTINE SOLB
C
      SUBROUTINE SOLB (N, NDIM, A, ML, MU, B, IP)
      REAL*8 A,B,T
      DIMENSION A(NDIM,N), B(N), IP(N)
C-----------------------------------------------------------------------
C  SOLUTION OF LINEAR SYSTEM, A*X = B .
C  INPUT..
C    N      ORDER OF MATRIX A.
C    NDIM   DECLARED DIMENSION OF ARRAY  A .
C    A      TRIANGULARIZED MATRIX OBTAINED FROM DECB.
C    ML     LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C    MU     UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C    B      RIGHT HAND SIDE VECTOR.
C    IP     PIVOT VECTOR OBTAINED FROM DECB.
C  DO NOT USE IF DECB HAS SET IER .NE. 0.
C  OUTPUT..
C    B      SOLUTION VECTOR, X .
C-----------------------------------------------------------------------
      MD = ML + MU + 1
      MD1 = MD + 1
      MDM = MD - 1
      NM1 = N - 1
      IF (ML .EQ. 0) GO TO 25
      IF (N .EQ. 1) GO TO 50
      DO 20 K = 1,NM1
        M = IP(K)
        T = B(M)
        B(M) = B(K)
        B(K) = T
        MDL = MIN(ML,N-K) + MD
        DO 10 I = MD1,MDL
          IMD = I + K - MD
 10       B(IMD) = B(IMD) + A(I,K)*T
 20     CONTINUE
 25   CONTINUE
      DO 40 KB = 1,NM1
        K = N + 1 - KB
        B(K) = B(K)/A(MD,K)
        T = -B(K)
        KMD = MD - K
        LM = MAX0(1,KMD+1)
        DO 30 I = LM,MDM
          IMD = I - KMD
 30       B(IMD) = B(IMD) + A(I,K)*T
 40     CONTINUE
 50   B(1) = B(1)/A(MD,1)
      RETURN
C----------------------- END OF SUBROUTINE SOLB ------------------------
      END

C----------------------- SUBROUTINE FOR THE MATRIX-VECTOR PRODUCT  ------------------------

C
C     SUBROUTINE MATVEC0
C
      subroutine MATVEC0(m,M0,ldmas,mlmas,mumas,v,Mv,ijob)
C      Input
      integer m,ldmas,mlmas,mumas,ijob
      double precision M0(ldmas,m),v(m)
C      Output
      double precision Mv(m)
C     Local variables
      integer i,j

      goto(10,20) ijob

10    continue
c     Full matrix
      do i=1,m
            Mv(i)=0d0
            do j=1,m
                  Mv(i)=Mv(i) + M0(i,j)*v(j)
            end do
      end do

      return

20    continue
c     Banded matrix
      do i=1,m
         Mv(i)=0d0
         do j=MAX(1,i-mlmas),MIN(m,i+mumas)
            Mv(i)=Mv(i)+M0(i-j+mumas+1,j)*v(j)
         end do
      end do

      return
C----------------------- END OF SUBROUTINE MATVEC0 ------------------------
      end
