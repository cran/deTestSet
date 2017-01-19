c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     This file is part of the Test Set for IVP solvers
c     http://www.dm.uniba.it/~testset/
c
c        NAND gate
c        index 0 IDE of dimension 14
c
c     DISCLAIMER: see
c     http://www.dm.uniba.it/~testset/disclaimer.php
c
c     The most recent version of this source file can be found at
c     http://www.dm.uniba.it/~testset/src/problems/nand.f
c
c     This is revision
c     $Id: nand.F,v 1.2 2006/10/02 10:29:14 testset Exp $
c
c-----------------------------------------------------------------------

c----------------------------------------------------------------------
c     parameter common block initialisation
c----------------------------------------------------------------------

      SUBROUTINE nandpar(daeparms)
      EXTERNAL daeparms

      double precision parms(14)
      common / nandcom / parms

        CALL daeparms(14,parms)

      END SUBROUTINE nandpar

c-----------------------------------------------------------------------
c     residual function
c-----------------------------------------------------------------------

      subroutine nandres(t,y,yprime,cj,f,ierr,rpar,ipar)
      integer neqn,ierr,ipar(*)
      parameter (neqn=14)
      double precision t,y(neqn),yprime(neqn),f(neqn),rpar(*)

      integer i,j
      double precision am(14,14),fy(14),dum
      ierr = 0
      call   CAP(14,y,AM)
      call   nandfunc(14,t,y,fy,ierr, rpar, ipar)

      if(ierr.eq.-1)return

      do 20 i=1,14
         dum = -fy(i)
         do 10 j=1,14
            dum = dum+AM(i,j)*yprime(j)
   10    continue
         f(i) = dum
   20 continue

      return
      end

c-----------------------------------------------------------------------

      subroutine nandsoln(neqn, y)
      integer neqn
      double precision y(14)
c
c computed at Cray C90, using Cray double precision:
C Solving NAND gate using PSIDE
C
C User input:
C
C give relative error tolerance: 1d-16
C give absolute error tolerance: 1d-16
C
C
C Integration characteristics:
C
C    number of integration steps       22083
C    number of accepted steps          21506
C    number of f evaluations          308562
C    number of Jacobian evaluations      337
C    number of LU decompositions       10532
C
C CPU-time used:                         451.71 sec

      y(  1) =  0.4971088699385777d+1
      y(  2) =  0.4999752103929311d+1
      y(  3) = -0.2499998781491227d+1
      y(  4) = -0.2499999999999975d+1
      y(  5) =  0.4970837023296724d+1
      y(  6) = -0.2091214032073855d+0
      y(  7) =  0.4970593243278363d+1
      y(  8) = -0.2500077409198803d+1
      y(  9) = -0.2499998781491227d+1
      y( 10) = -0.2090289583878100d+0
      y( 11) = -0.2399999999966269d-3
      y( 12) = -0.2091214032073855d+0
      y( 13) = -0.2499999999999991d+1
      y( 14) = -0.2500077409198803d+1

      return
      end
c-----------------------------------------------------------------------
      SUBROUTINE nandfunc(N,T,Y,F,ierr, rpar, ipar)
C ---------------------------------------------------------------------
C
C Right-hand side f(X,t) for the network equation
C             C(Y) * Y' - f(Y,t) = 0
C describing the nand gate
C
C ---------------------------------------------------------------------
C
C Input parameters:
C          N......number of node potentials (14)
C          T......time point t
C          Y......node potentials at time point t
C Output parameter:
C          F......right-hand side f(Y,t)
C
C External reference:
C          nandIDS: Drain-source current
C          nandIBS: Nonlinear current characteristic for diode between
C               bulk and source
C          nandIBD: Nonlinear current characteristic for diode between
C               bulk and drain
C          nandPULSE: Input signal in pulse form
C
C ---------------------------------------------------------------------

      INTEGER N,ierr
      double precision T,Y(N),F(N),nandIDS,nandIBS,nandIBD,V1,V2,V1D,V2D
      EXTERNAL nandIDS, nandIBS, nandIBD, nandPULSE

      double precision  RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                  DELTA, CURIS, VTH, VDD, VBB, rpar(*)
      integer ipar(*)

      COMMON /nandcom/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *               DELTA, CURIS, VTH, VDD, VBB

      CALL nandPULSE(T,V1,V1D,0.D0,5.D0,5.D0,5.D0,5.D0,5.D0,20.D0)
      CALL nandPULSE(T,V2,V2D,0.D0,5.D0,15.D0,5.D0,15.D0,5.D0,40.D0)


      F(1)=-(Y(1)-Y(5))/RGS-nandIDS(1,Y(2)-Y(1),Y(5)-Y(1),Y(3)-Y(5),
     *       Y(5)-Y(2),Y(4)-VDD,ierr)
      F(2)=-(Y(2)-VDD)/RGD+nandIDS(1,Y(2)-Y(1),Y(5)-Y(1),Y(3)-Y(5),
     *       Y(5)-Y(2),Y(4)-VDD,ierr)
      F(3)=-(Y(3)-VBB)/RBS + nandIBS(Y(3)-Y(5))
      F(4)=-(Y(4)-VBB)/RBD + nandIBD(Y(4)-VDD)
      F(5)=-(Y(5)-Y(1))/RGS-nandIBS(Y(3)-Y(5))-(Y(5)-Y(7))/RGD-
     *       nandIBD(Y(9)-Y(5))
      F(6)=CGS*V1D-(Y(6)-Y(10))/RGS-
     *   nandIDS(2,Y(7)-Y(6),V1-Y(6),Y(8)-Y(10),V1-Y(7),Y(9)-Y(5),ierr)
      F(7)=CGD*V1D-(Y(7)-Y(5))/RGD+
     *   nandIDS(2,Y(7)-Y(6),V1-Y(6),Y(8)-Y(10),V1-Y(7),Y(9)-Y(5),ierr)
      F(8)=-(Y(8)-VBB)/RBS + nandIBS(Y(8)-Y(10))
      F(9)=-(Y(9)-VBB)/RBD + nandIBD(Y(9)-Y(5))
      F(10)=-(Y(10)-Y(6))/RGS-nandIBS(Y(8)-Y(10))-
     *         (Y(10)-Y(12))/RGD-nandIBD(Y(14)-Y(10))
      F(11)=CGS*V2D-Y(11)/RGS-nandIDS(2,Y(12)-Y(11),V2-Y(11),Y(13),
     *       V2-Y(12),Y(14)-Y(10),ierr)
      F(12)=CGD*V2D-(Y(12)-Y(10))/RGD+
     *   nandIDS(2,Y(12)-Y(11),V2-Y(11),Y(13),V2-Y(12),Y(14)-Y(10),ierr)
      F(13)=-(Y(13)-VBB)/RBS + nandIBS(Y(13))
      F(14)=-(Y(14)-VBB)/RBD + nandIBD(Y(14)-Y(10))

c      if(ierr.eq.-1)THEN
c         WRITE(MSG, *)"AN ERROR OCCURRED in NAND, at time", T
c         call rexit(MSG)
c      ENDIF
      RETURN
      END

      double precision FUNCTION nandIDS (NED,VDS,VGS,VBS,VGD,VBD,ierr)
C ---------------------------------------------------------------------------
C
C Function evaluating the drain-current due to the model of
C Shichman and Hodges
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   NED  Integer parameter for MOSFET-type
C   VDS  Voltage between drain and source
C   VGS  Voltage between gate and source
C   VBS  Voltage between bulk and source
C   VGD  Voltage between gate and drain
C   VBD  Voltage between bulk and drain
C
C External reference:
C   nandGDSP, nandGDSM Drain function for VDS > 0 gevalp. VDS < 0
C
C ---------------------------------------------------------------------------

      INTEGER NED,ierr
      double precision VDS, VGS, VBS, VGD, VBD,nandGDSP, nandGDSM
      EXTERNAL nandGDSP, nandGDSM

      IF ( VDS .GT. 0.D0 ) THEN
       nandIDS = nandGDSP (NED,VDS, VGS, VBS,ierr)
      ELSE IF ( VDS .EQ. 0.D0) THEN
       nandIDS = 0.D0
      ELSE IF ( VDS .LT. 0.D0) THE N
       nandIDS = nandGDSM (NED,VDS, VGD, VBD,ierr)
      END IF

      if(ierr.eq.-1)return

      RETURN
      END

      double precision FUNCTION nandGDSP (NED,VDS, VGS, VBS, ierr)
      integer NED,ierr
      double precision  VDS, VGS, VBS,VTE

      double precision  VT0, BETA, CGAMMA, PHI

      double precision  RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                  DELTA, CURIS, VTH, VDD, VBB

      COMMON /nandcom/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *               DELTA, CURIS, VTH, VDD, VBB

C


      IF(NED.EQ.1) THEN
C --- Depletion-type
        VT0=-2.43D0
        CGAMMA=.2D0
        PHI=1.28D0
        BETA=5.35D-4
      ELSE
C --- Enhancement-type
        VT0=.2D0
        CGAMMA=0.035D0
        PHI=1.01D0
        BETA=1.748D-3
      END IF
      if(phi-vbs.lt.0d0.or.phi.lt.0d0)then
         ierr=-1
!         return
         nandGDSP = 0.D0
         goto 100
      end if

      VTE = VT0 + CGAMMA * ( DSQRT(PHI-VBS) - DSQRT(PHI) )

      IF ( VGS-VTE .LE. 0.D0) THEN
       nandGDSP = 0.D0
      ELSE IF ( 0.D0 .LT. VGS-VTE .AND. VGS-VTE .LE. VDS ) THEN
       nandGDSP = - BETA * (VGS - VTE)**2.D0 * (1.D0 + DELTA*VDS)
      ELSE IF ( 0.D0 .LT. VDS .AND. VDS .LT. VGS-VTE ) THEN
       nandGDSP = - BETA * VDS * (2.D0*(VGS - VTE) - VDS) *
     *          (1.D0 + DELTA*VDS)
      END IF

100   RETURN
      END

      double precision FUNCTION nandGDSM (NED,VDS, VGD, VBD, ierr)
      integer NED,ierr
      double precision VDS, VGD, VBD,VTE


      double precision  VT0, BETA, CGAMMA, PHI

      double precision  RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                  DELTA, CURIS, VTH, VDD, VBB

      COMMON /nandcom/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *               DELTA, CURIS, VTH, VDD, VBB


      IF(NED.EQ.1) THEN
C --- Depletion-type
      VT0=-2.43D0
      CGAMMA=.2D0
      PHI=1.28D0
      BETA=5.35D-4
      ELSE
C --- Enhancement-type
      VT0=.2D0
      CGAMMA=0.035D0
      PHI=1.01D0
      BETA=1.748D-4
      END IF

      if(phi-vbd.lt.0d0.or.phi.lt.0d0)then
         ierr=-1
!         return
         nandGDSM = 0.D0
         goto 100
      end if

      VTE = VT0 + CGAMMA * ( DSQRT(PHI-VBD) - DSQRT(PHI) )

      IF ( VGD-VTE .LE. 0.D0) THEN
       nandGDSM = 0.D0
      ELSE IF ( 0.D0 .LT. VGD-VTE .AND. VGD-VTE .LE. -VDS ) THEN
       nandGDSM = BETA * (VGD - VTE)**2d0 * (1.D0 - DELTA*VDS)
      ELSE IF ( 0.D0 .LT. -VDS .AND. -VDS .LT. VGD-VTE ) THEN
       nandGDSM = - BETA * VDS * (2d0 *(VGD - VTE) + VDS) *
     *          (1.D0 - DELTA*VDS)
      END IF

100   RETURN
      END


      double precision FUNCTION nandIBS (VBS)
C ---------------------------------------------------------------------------
C
C Function evaluating the current of the pn-junction between bulk and
C source due to the model of Shichman and Hodges
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   VBS  Voltage between bulk and source
C
C ---------------------------------------------------------------------------

      double precision VBS

      double precision  RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                  DELTA, CURIS, VTH, VDD, VBB

      COMMON /nandcom/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *               DELTA, CURIS, VTH, VDD, VBB

C

C
C     nandIBS = GBS (VBS)
C

      IF ( VBS .LE. 0.D0 ) THEN
       nandIBS = - CURIS * ( DEXP( VBS/VTH ) - 1.D0 )
      ELSE
       nandIBS = 0.D0
      END IF

      RETURN
      END

      double precision FUNCTION nandIBD (VBD)
C ---------------------------------------------------------------------------
C
C Function evaluating the current of the pn-junction between bulk and
C drain  due to the model of Shichman and Hodges
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   VBS  Voltage between bulk and drain
C
C ---------------------------------------------------------------------------

      double precision VBD

      double precision  RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                  DELTA, CURIS, VTH, VDD, VBB

      COMMON /nandcom/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *               DELTA, CURIS, VTH, VDD, VBB


C

C
C     nandIBD = GBD (VBD)
C
      IF ( VBD .LE. 0.D0 ) THEN
       nandIBD = - CURIS * ( DEXP( VBD/VTH ) - 1.D0 )
      ELSE
       nandIBD = 0.D0
      END IF
      RETURN
      END

      SUBROUTINE nandPULSE(X,VIN,VIND,LOW,HIGH,DELAY,T1,T2,T3,PERIOD)
C ---------------------------------------------------------------------------
C
C Evaluating input signal at time point X
C
C Structure of input signal:
C
C                -----------------------                       HIGH
C               /                       \
C              /                         \
C             /                           \
C            /                             \
C           /                               \
C          /                                 \
C         /                                   \
C        /                                     \
C  ------                                       ---------      LOW
C
C |DELAY|   T1  |         T2           |   T3  |
C |          P     E     R     I     O     D            |
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   X                      Time-point at which input signal is evaluated
C   LOW                    Low-level of input signal
C   HIGH                   High-level of input signal
C   DELAY,T1,T2,T3, PERIOD Parameters to specify signal structure
C
C Output parameter:
C   VIN    Voltage of input signal at time point X
C   VIND   Derivative of VIN at time point X
C
C ---------------------------------------------------------------------------

      double precision X,VIN,VIND,LOW,HIGH,DELAY,T1,T2,T3,PERIOD,TIME

      TIME = DMOD(X,PERIOD)

      IF (TIME.GT.(DELAY+T1+T2+T3)) THEN
      VIN = LOW
      VIND= 0.D0
      ELSE IF (TIME.GT.(DELAY+T1+T2)) THEN
      VIN = ((HIGH-LOW)/T3)*(DELAY+T1+T2+T3-TIME) + LOW
      VIND= -((HIGH-LOW)/T3)
      ELSE IF (TIME.GT.(DELAY+T1)) THEN
      VIN = HIGH
      VIND= 0.D0
      ELSE IF (TIME.GT.DELAY) THEN
      VIN = ((HIGH-LOW)/T1)*(TIME-DELAY) + LOW
      VIND= ((HIGH-LOW)/T1)
      ELSE
      VIN = LOW
      VIND=0.D0
      END IF

      RETURN
      END

      SUBROUTINE CAP(N,Y,AM)
C ---------------------------------------------------------------------
C
C Voltage-dependent capacitance matrix C(Y) for the network equation
C             C(Y) * Y' - f(Y,t) = 0
C describing the nand gate
C
C ---------------------------------------------------------------------
C
C Input parameters:
C          N......number of node potentials (14)
C          Y......value of node potentials
C Output parameter:
C          AM.....voltage-dependent capacitance matrix
C
C External reference:
C          nandCBDBS: Voltage-dependent capacitance CBS(V) and CBD(V)
C
C ---------------------------------------------------------------------
      double precision nandCBDBS
      INTEGER N
      double precision Y(N), AM(N,N)

      double precision  RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                  DELTA, CURIS, VTH, VDD, VBB

      COMMON /nandcom/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *               DELTA, CURIS, VTH, VDD, VBB

      EXTERNAL nandCBDBS
      integer I,J

      DO 10 I=1,N
        DO 20 J=1,N
          AM(I,J)=0d0
 20     CONTINUE
 10   CONTINUE

      AM(1,1)=CGS
      AM(1,5)=-CGS
      AM(2,2)=CGD
      AM(2,5)=-CGD
      AM(3,3)=nandCBDBS(Y(3)-Y(5))
      AM(3,5)=-nandCBDBS(Y(3)-Y(5))
      AM(4,4)=nandCBDBS(Y(4)-VDD)
      AM(5,1)=-CGS
      AM(5,2)=-CGD
      AM(5,3)=-nandCBDBS(Y(3)-Y(5))
      AM(5,5)=CGS+CGD-AM(5,3)+
     *          nandCBDBS(Y(9)-Y(5))+C9
      AM(5,9)=-nandCBDBS(Y(9)-Y(5))
      AM(6,6)=CGS
      AM(7,7)=CGD
      AM(8,8)=nandCBDBS(Y(8)-Y(10))
      AM(8,10)=-nandCBDBS(Y(8)-Y(10))
      AM(9,5)=-nandCBDBS(Y(9)-Y(5))
      AM(9,9)=nandCBDBS(Y(9)-Y(5))
      AM(10,8)=-nandCBDBS(Y(8)-Y(10))
      AM(10,10)=-AM(8,10)+nandCBDBS(Y(14)-Y(10))+C9
      AM(10,14)=-nandCBDBS(Y(14)-Y(10))
      AM(11,11)=CGS
      AM(12,12)=CGD
      AM(13,13)=nandCBDBS(Y(13))
      AM(14,10)=-nandCBDBS(Y(14)-Y(10))
      AM(14,14)=nandCBDBS(Y(14)-Y(10))

      RETURN
      END

      double precision FUNCTION nandCBDBS (V)
C ---------------------------------------------------------------------------
C
C Function evaluating the voltage-dependent capacitance between bulk and
C drain gevalp. source  due to the model of Shichman and Hodges
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   V    Voltage between bulk and drain gevalp. source
C
C ---------------------------------------------------------------------------
      double precision V,PHIB

      double precision  RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                  DELTA, CURIS, VTH, VDD, VBB

      COMMON /nandcom/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *               DELTA, CURIS, VTH, VDD, VBB

      PHIB=0.87D0

      IF ( V .LE. 0.D0 ) THEN
       nandCBDBS = CBD/DSQRT(1.D0-V/PHIB)
      ELSE
       nandCBDBS = CBD*(1.D0+V/(2.D0*PHIB))
      END IF

      RETURN
      END

