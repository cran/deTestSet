      SUBROUTINE CASHKARP(N,FCN,X,Y,XEND,
     &                  RTOL,ATOL,ITOL,
     &                  SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID, VSTI)
C ----------------------------------------------------------
C     NUMERICAL SOLUTION OF A SYSTEM OF FIRST 0RDER
C     ORDINARY DIFFERENTIAL EQUATIONS  Y'=F(X,Y).
C     THIS IS AN EXPLICIT RUNGE-KUTTA METHOD OF ORDER  5(4)
C     DUE TO CASH AND KARP (WITH STEPSIZE CONTROL ,
C     DENSE OUTPUT and AUTOMATIC STIFFNESS DETECTION.)
C
C     AUTHOR : J. R. Cash
C              IMPERIAL COLLEGE LONDON
C              SOUTH KENSINGTON, LONDON S.W. 7, ENGLAND
C              E-MAIL:  j.cash@imperial.ac.uk
C
C              F. Mazzia
C              DIPARTMENTO DI MATEMATICA
C              UNIVERSITA' DI BARI, ITALY
C              E-MAIL: mazzia@dm.uniba.it
C
C     THIS CODE IS BASED ON THE CODE DOPRI5 WHICH IS DESCRIBED IN:
C         E. HAIRER, S.P. NORSETT AND G. WANNER, SOLVING ORDINARY
C         DIFFERENTIAL EQUATIONS I. NONSTIFF PROBLEMS. 2ND EDITION.
C         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
C         SPRINGER-VERLAG (1993)
C         The actual Code is given in J.R. CASH and Alan H. Karp,
C         A Variable Order Runge-Kutta Method for Initial Value
C         Problems With Rapidly Varying Right Hand Sides, ACM TOMS,
C         Vol 16, Sept 1990, pp. 201-222.
C
C         The stiffness detection algorithm based on conditioning is
C         given in:F. Mazzia, A.M. Nagy, Stiffness detection with
C         Runge-Kutta methods,
C         Rapporto del Dipartimento di Matematica
C         n.21 del 23 Novembre 2010.
C
C
C
C     INPUT PARAMETERS
C     ----------------
C     N           DIMENSION OF THE SYSTEM
C
C     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE
C                 VALUE OF F(X,Y):
C                    SUBROUTINE FCN(N,X,Y,F,IERR,RPAR,IPAR)
C                    DOUBLE PRECISION X,Y(N),F(N)
C                    F(1)=...   ETC.
C
C     karline: removed the IERR argument for consistency with other solvers...
C
C     X           INITIAL X-VALUE
C
C     Y(N)        INITIAL VALUES FOR Y
C
C     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE)
C
C     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY
C                 CAN BE BOTH SCALARS OR  BOTH VECTORS OF LENGTH N.
C
C     ITOL        SWITCH FOR RTOL AND ATOL:
C                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS.
C                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF
C                     Y(I) BELOW RTOL*ABS(Y(I))+ATOL
C                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS.
C                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW
C                     RTOL(I)*ABS(Y(I))+ATOL(I).
C
C     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE
C                 NUMERICAL SOLUTION DURING INTEGRATION.
C                 IF IOUT.GE.1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP.
C                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0.
C                 IT MUST HAVE THE FORM
C                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,N,CON,ICOMP,ND,
C                                       RPAR,IPAR,IRTRN)
C                    DIMENSION Y(N),CON(*),ICOMP(ND)
C                    ....
C                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH
C                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS
C                    THE FIRST GRID-POINT).
C                 "XOLD" IS THE PRECEEDING GRID-POINT.
C                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
C                    IS SET <0, CASHKARP WILL RETURN TO THE CALLING PROGRAM.
C                    IF THE NUMERICAL SOLUTION IS ALTERED IN SOLOUT,
C                    SET  IRTRN = 2
C
C          -----  CONTINUOUS OUTPUT: -----
C                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION
C                 FOR THE INTERVAL [XOLD,X] IS AVAILABLE THROUGH
C                 THE FUNCTION
C                        >>>   CONTCK(I,S,CON,ICOMP,ND)   <<<
C                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH
C                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE
C                 S SHOULD LIE IN THE INTERVAL [XOLD,X].
C
C     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT:
C                    IOUT=0: SUBROUTINE IS NEVER CALLED
C                    IOUT=1: SUBROUTINE IS USED FOR OUTPUT.
C                    IOUT=2: DENSE OUTPUT IS PERFORMED IN SOLOUT
C                            (IN THIS CASE WORK(5) MUST BE SPECIFIED)
C
C     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK".
C                 WORK(1),...,WORK(20) SERVE AS PARAMETERS FOR THE CODE.
C                 FOR STANDARD USE, SET THEM TO ZERO BEFORE CALLING.
C                 "LWORK" MUST BE AT LEAST  19*N+7*NRDENS+21
C                 WHERE  NRDENS = IWORK(5)
C
C     LWORK       DECLARED LENGTH OF ARRAY "WORK".
C
C     IWORK       INTEGER WORKING SPACE OF LENGTH "LIWORK".
C                 IWORK(1),...,IWORK(20) SERVE AS PARAMETERS FOR THE CODE.
C                 FOR STANDARD USE, SET THEM TO ZERO BEFORE CALLING.
C                 "LIWORK" MUST BE AT LEAST NRDENS+21 .
C
C     LIWORK      DECLARED LENGTH OF ARRAY "IWORK".
C
C     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH
C                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING
C                 PROGRAM AND THE FCN, JAC, MAS, SOLOUT SUBROUTINES.
C
C    VSTI IS AN ARRAY OF WORKSPACE OF LENGTH N. IT IS NOT NEEDED IF
C         AUTOMATIC STIFFNESS DETECTION IS NOT USED. now is in WORK
C
C     SOPHISTICATED SETTING OF PARAMETERS
C     -----------------------------------
C              SEVERAL PARAMETERS (WORK(1),..,IWORK(1),...) ALLOW us
C              TO ADAPT THE CODE TO THE PROBLEM AND TO THE NEEDS OF
C              THE USER. FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES.
C
C    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 2.3D-16.
C
C    WORK(2)   THE SAFETY FACTOR IN STEP SIZE PREDICTION,
C              DEFAULT 0.9D0.
C
C    WORK(3), WORK(4)   PARAMETERS FOR STEP SIZE SELECTION
C              THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION
C                 WORK(3) <= HNEW/HOLD <= WORK(4)
C              DEFAULT VALUES: WORK(3)=0.2D0, WORK(4)=10.D0
C
C    WORK(5)   IS THE "BETA" FOR STABILIZED STEP SIZE CONTROL
C              (SEE SECTION IV.2). LARGER VALUES OF BETA ( <= 0.1 )
C              MAKE THE STEP SIZE CONTROL MORE STABLE. CASHKARP NEEDS
C              A LARGER BETA THAN HIGHAM & HALL. NEGATIVE WORK(5)
C              CAUSES BETA=0.
C              THE DEFAULT VALUE IS  0.04D0.
C
C    WORK(6)   MAXIMAL STEP SIZE, DEFAULT XEND-X.
C
C    WORK(7)   INITIAL STEP SIZE, FOR WORK(7)=0.D0 AN INITIAL GUESS
C              IS COMPUTED WITH HELP OF THE FUNCTION HINIT
C
C    IWORK(1)  THIS IS THE MAXIMUM NUMBER OF ALLOWED STEPS.
C              THE DEFAULT VALUE (FOR IWORK(1)=0) IS 100000.
C
C    IWORK(2)  SWITCH FOR THE CHOICE OF THE COEFFICIENTS
C              IF IWORK(2).EQ.1  THE METHOD   USED IS CASHKARP
C              AT THE MOMENT THIS IS THE ONLY POSSIBLE CHOICE.
C              THE DEFAULT VALUE (FOR IWORK(2)=0) IS IWORK(2)=1.
C
C    IWORK(3)  SWITCH FOR PRINTING ERROR MESSAGES
C              IF IWORK(3).LT.0 NO MESSAGES ARE BEING PRINTED
C              DEFAULT VALUE (FOR IWORK(3)=0) IS IWORK(3)=6
C
C    IWORK(4)  TEST FOR STIFFNESS USING EIGENVALUE APPROXIMATION
C              IS ACTIVATED AFTER STEP 1
c             IWORK(4) = -1 NO TEST
C             IWORK(4) = 1 TEST NO STOP, WARNING
C             IWORK(4) = 2 TEST STOP IF STIFF
C             DEFAULT VALUE IS IWORK(4)=1
C
C    IWORK(5)  = NRDENS = NUMBER OF COMPONENTS, FOR WHICH DENSE OUTPUT
C              IS REQUIRED; DEFAULT VALUE IS IWORK(5)=0;
C              FOR   0 < NRDENS < N   THE COMPONENTS (FOR WHICH DENSE
C              OUTPUT IS REQUIRED) HAVE TO BE SPECIFIED IN
C              IWORK(21),...,IWORK(NRDENS+20);
C              FOR  NRDENS=N  THIS IS DONE BY THE CODE.
C
C    IWORK(6)  TEST FOR STIFFNESS BASED ON ERROR ESTIMATION
C             IWORK(6) = -1 NO TEST
C             IWORK(6) = 1 TEST NO STOP,WARNING
C             IWORK(6) = 2 TEST STOP IF STIFF
C             DEFAULT VALUE IS IWORK(6)=1
C    IWORK(7) TEST FOR STIFFNESS USING CONDITIONING
C              IWORK(7) = -1 NO TEST
C              IWORK(7) = 1 TEST NO STOP, WARNING
C              IWORK(7) = 2 TEST STOP IF STIFF
C              DEFAULT VALUE IS IWORK(7)=1
C----------------------------------------------------------------------
C
C     OUTPUT PARAMETERS
C     -----------------
C     X           X-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED
C                 (AFTER SUCCESSFUL RETURN X=XEND).
C
C     Y(N)        NUMERICAL SOLUTION AT X
C
C     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP
C
C     IDID        REPORTS ON SUCCESS UPON RETURN:
C                   IDID= 1  COMPUTATION SUCCESSFUL,
C                   IDID= 2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT)
C                   IDID=-1  INPUT IS NOT CONSISTENT,
C                   IDID=-2  LARGER NMAX IS NEEDED,
C                   IDID=-3  STEP SIZE BECOMES TOO SMALL.
C                   IDID=-4  PROBLEM IS PROBABLY STIFF (INTERRUPTED).
C
C   IWORK(17)  NFCN    NUMBER OF FUNCTION EVALUATIONS
C   IWORK(18)  NSTEP   NUMBER OF COMPUTED STEPS
C   IWORK(19)  NACCPT  NUMBER OF ACCEPTED STEPS
C   IWORK(20)  NREJCT  NUMBER OF REJECTED STEPS (DUE TO ERROR TEST),
C                      (STEP REJECTIONS IN THE FIRST STEP ARE NOT COUNTED)
C
C   WORK(7) = H (stepsize used in the last step)
C   IF NSTIFFCOND > 0
C   WORK(8) = KAPPA
C   WORK(9) = GAMMA
C   WORK(10) = SIGMA
C   WORK(11) = XST
C   WORK(12) = SIGMATOT
C-----------------------------------------------------------------------
C *** *** *** *** *** *** *** *** *** *** *** *** ***
C          DECLARATIONS
C *** *** *** *** *** *** *** *** *** *** *** *** ***
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(N),ATOL(*),RTOL(*),WORK(LWORK),IWORK(LIWORK)
      DIMENSION RPAR(*),IPAR(*),VSTI(N)
      DOUBLE PRECISION KAPPA, GAMMA,SIGMA,SIGMATOT
      LOGICAL ARRET
      EXTERNAL FCN
      EXTERNAL SOLOUT
C *** *** *** *** *** *** ***
C        SETTING THE PARAMETERS
C *** *** *** *** *** *** ***
      NFCN=0
      NSTEP=0
      NACCPT=0
      NREJCT=0
      ARRET=.FALSE.
C -------- IPRINT FOR MONITORING THE PRINTING
      IF(IWORK(3).EQ.0)THEN
         IPRINT=0
      ELSE
         IPRINT=IWORK(3)
      END IF
C -------- NMAX , THE MAXIMAL NUMBER OF STEPS -----
      IF(IWORK(1).EQ.0)THEN
         NMAX=100000
      ELSE
         NMAX=IWORK(1)
         IF(NMAX.LE.0)THEN
            IF (IPRINT.GT.0) THEN
        CALL Rprinti1('wrong input iwork(1) = ',IWORK(1))
            END IF
            ARRET=.TRUE.
         END IF
      END IF
C -------- METH   COEFFICIENTS OF THE METHOD
      IF(IWORK(2).EQ.0)THEN
         METH=1
      ELSE
         METH=IWORK(2)
         IF(METH.LE.0.OR.METH.GE.4)THEN
            IF (IPRINT.GT.0) THEN
        CALL Rprinti1('Curious input iwork(2) = ',IWORK(2))
            END IF
            ARRET=.TRUE.
         END IF
      END IF
C -------- NSTIFF   PARAMETER FOR STIFFNESS DETECTION

      NSTIFF=IWORK(4)
      IF (NSTIFF.EQ.0) NSTIFF=1
C -------- NRDENS   NUMBER OF DENSE OUTPUT COMPONENTS
      NRDENS=IWORK(5)
      IF(NRDENS.LT.0.OR.NRDENS.GT.N)THEN
         IF (IPRINT.GT.0) THEN
        CALL Rprinti1('Curious input iwork(5) = ',IWORK(5))
         END IF
         ARRET=.TRUE.
      ELSE
            IF(NRDENS.GT.0.AND.IOUT.LT.2)THEN
               IF (IPRINT.GT.0)  THEN
        CALL Rprint('Warning: put IOUT=2 for dense output ')
               END IF
            END IF
            IF (NRDENS.EQ.N) THEN
                DO 16 I=1,NRDENS
  16            IWORK(20+I)=I
            END IF
      END IF
      NSTIFFCOND=IWORK(7)
      IF (NSTIFFCOND.EQ.0) NSTIFFCOND=1
      NSTIFFERR = IWORK(6)
      IF (NSTIFFERR.EQ.0) NSTIFFERR=1
C -------- UROUND   SMALLEST NUMBER SATISFYING 1.D0+UROUND>1.D0

      IF(WORK(1).EQ.0.D0)THEN
         UROUND=2.3D-16
      ELSE
         UROUND=WORK(1)
         IF(UROUND.LE.1.D-35.OR.UROUND.GE.1.D0)THEN
            IF (IPRINT.GT.0) THEN
        CALL Rprintd1(
     &        ' Which machine do you have? Your uround was: ',WORK(1))
            END IF
            ARRET=.TRUE.
         END IF
      END IF
C -------  SAFETY FACTOR -------------
      IF(WORK(2).EQ.0.D0)THEN
         SAFE=0.9D0
      ELSE
         SAFE=WORK(2)
         IF(SAFE.GE.1.D0.OR.SAFE.LE.1.D-4)THEN
            IF (IPRINT.GT.0) THEN
        CALL Rprintd1(
     &          'Curious input for safety factor work(2) = ',WORK(2))
            END IF
            ARRET=.TRUE.
         END IF
      END IF
C -------  FAC1,FAC2     PARAMETERS FOR STEP SIZE SELECTION
      IF(WORK(3).EQ.0.D0)THEN
         FAC1=0.2D0
      ELSE
         FAC1=WORK(3)
      END IF
      IF(WORK(4).EQ.0.D0)THEN
         FAC2=10.D0
      ELSE
         FAC2=WORK(4)
      END IF
C --------- BETA FOR STEP CONTROL STABILIZATION -----------
      IF(WORK(5).EQ.0.D0)THEN
         BETA=0.04D0
      ELSE
         IF(WORK(5).LT.0.D0)THEN
            BETA=0.D0
         ELSE
            BETA=WORK(5)
            IF(BETA.GT.0.2D0)THEN
               IF (IPRINT.GT.0) THEN
        CALL Rprintd1('Curious input for beta: work(5) = ',WORK(5))
               END IF
            ARRET=.TRUE.
         END IF
         END IF
      END IF
C -------- MAXIMAL STEP SIZE
      IF(WORK(6).EQ.0.D0)THEN
         HMAX=XEND-X
      ELSE
         HMAX=WORK(6)
      END IF
C -------- INITIAL STEP SIZE
      H=WORK(7)
C ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK -----
      IEY1=21
      IEK1=IEY1+N
      IEK2=IEK1+N
      IEK3=IEK2+N
      IEK4=IEK3+N
      IEK5=IEK4+N
      IEK6=IEK5+N
      IEYS=IEK6+N
      IEVS=IEYS+N
      IEYSTI=IEVS+N
      IEYS=IEYSTI+N
      IEY1S=IEYS+N
      IEK1S=IEY1S+N
      IEK2S=IEK1S+N
      IEK3S=IEK2S+N
      IEK4S=IEK3S+N
      IEk5S=IEK4S+N
      IEK6S=IEK5S+N
      IEYSTIS=IEK6S+N
      IECO=IEYSTIS+N
C --- starting for IECO end after 7*NRDENS
C ------ TOTAL STORAGE REQUIREMENT -----------
      ISTORE=IEYSTIS+7*NRDENS-1
      IF(ISTORE.GT.LWORK)THEN
        IF (IPRINT.GT.0)  THEN
        CALL Rprinti1('Insufficient storage for work, min. = ',ISTORE)
        END IF
        ARRET=.TRUE.
      END IF
      ICOMP=21
      ISTORE=ICOMP+NRDENS-1
      IF(ISTORE.GT.LIWORK)THEN
        IF (IPRINT.GT.0) THEN
        CALL Rprinti1('Insufficient storage for iwork, min. = ',ISTORE)
        END IF
        ARRET=.TRUE.
      END IF
C ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1
      IF (ARRET) THEN
         IDID=-1
         RETURN
      END IF
C -------- CALL TO CORE INTEGRATOR ------------
      CALL CKCOR(N,FCN,X,Y,XEND,HMAX,H,RTOL,ATOL,ITOL,IPRINT,
     &   SOLOUT,IOUT,IDID,NMAX,UROUND,METH,NSTIFF,NSTIFFCOND,NSTIFFERR,
     &   SAFE,BETA,FAC1,FAC2,
     &   WORK(IEY1),WORK(IEK1),WORK(IEK2),WORK(IEK3),WORK(IEK4),
     &   WORK(IEK5),WORK(IEK6),WORK(IEYSTI),WORK(IECO),IWORK(ICOMP),
     &   NRDENS,RPAR,IPAR,NFCN,NSTEP,NACCPT,NREJCT,WORK(IEVS),
     &   WORK(IEYS),WORK(IEY1S),WORK(IEK1S),WORK(IEK2S),
     &   WORK(IEK3S),WORK(IEK4S),WORK(IEK5S),WORK(IEK6S),
     &   WORK(IEYSTIS),KAPPA,GAMMA,SIGMA,SIGMATOT,XST)
      WORK(7)=H
      IF (NSTIFFCOND .GT.  0) THEN
          WORK(8) = KAPPA
          WORK(9) = GAMMA
          WORK(10) = SIGMA
          WORK(11) = XST
          WORK(12) = SIGMATOT
      END IF
      IWORK(17)=NFCN
      IWORK(18)=NSTEP
      IWORK(19)=NACCPT
      IWORK(20)=NREJCT

C ----------- RETURN -----------
      RETURN
      END
C
C
C
C  ----- ... AND HERE IS THE CORE INTEGRATOR  ----------
C
      SUBROUTINE CKCOR(N,FCN,X,Y,XEND,HMAX,H,RTOL,ATOL,ITOL,IPRINT,
     &   SOLOUT,IOUT,IDID,NMAX,UROUND,METH,NSTIFF,NSTIFFCOND,NSTIFFERR,
     &   SAFE,BETA,FAC1,FAC2,
     &   Y1,K1,K2,K3,K4,K5,K6,YSTI,CONT,ICOMP,NRD,RPAR,IPAR,
     &   NFCN,NSTEP,NACCPT,NREJCT,VSTI,
     &   YS,Y1S,K1S,K2S,K3S,K4S,K5S,K6S,YSTIS,
     &   KAPPA,GAMMA,SIGMA,SIGMATOT,XST)
C ----------------------------------------------------------
C     CORE INTEGRATOR FOR CASHKARP
C     PARAMETERS SAME AS IN CASHKARP WITH WORKSPACE ADDED
C ----------------------------------------------------------
C         DECLARATIONS
C ----------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION K1(N),K2(N),K3(N),K4(N),K5(N),K6(N)
      DIMENSION Y(N),Y1(N),YSTI(N),ATOL(*),RTOL(*),RPAR(*),IPAR(*)
      DOUBLE PRECISION K1S(N),K2S(N),K3S(N),K4S(N),K5S(N),K6S(N),k4I
      DOUBLE PRECISION KAPPA, GAMMA,SIGMA,SIGMATOT
      DOUBLE PRECISION YS(N),Y1S(N),VSTI(N), YSTIS(N)
      DIMENSION CONT(*),ICOMP(NRD)
      DOUBLE PRECISION ETAVE(N)
      DOUBLE PRECISION  D1,D3,D4,D5,D6,D7
      LOGICAL REJECT,LAST
      EXTERNAL FCN
      LOGICAL HERMITE
      COMMON /CONTCKV/XOLD,HOUT,HERMITE
C *** *** *** *** *** *** ***
C  INITIALISATIONS
C *** *** *** *** *** *** ***
      IF (METH.EQ.1) CALL CKC(C2,C3,C4,C5,C6,E1,E3,E4,E5,E6,
     &               A21,A31,A32,A41,A42,A43,A51,A52,A53,A54,
     &               A61,A62,A63,A64,A65,A71,A72,A73,A74,A75,A76,
     &               D1,D3,D4,D5,D6,D7,ad1,ad2,ad3,ad4,ad5,ad6)

      FACOLD=1.0d-4
      EXPO1=0.2D0-beta*0.75d0
      FACC1=1.D0/FAC1
      FACC2=1.D0/FAC2
      POSNEG=SIGN(1.D0,XEND-X)
C --- INITIAL PREPARATIONS
      ATOLI=ATOL(1)
      RTOLI=RTOL(1)
      LAST=.FALSE.
      HLAMB=0.D0
      IASTI=0
      IASTIS=0
      NONSTI=0
      
cF we use the shampine algorithm to compute the contnuous
cF extension
      HERMITE = .true.
cF computation of eta vector in the direction of
cF dominant eigenvalue of the Jacobian for stiffness detection
      IERR = 0
cF used in the error check

      IF (NSTIFFCOND .GT. 0) THEN

        IF (N .eq. 1) THEN
cF scalar problems
          ETAVE(1) = 1.0d0
          IFLAGST = 0
        ELSE
cF vectorial problems
cF initialization of etave
         DO I=1,N
          ETAVE(I) = 0.0D0
         END DO
         XS = X
         XSTOP = XEND
cF initial stepsize for the computation of etave
         hst = MIN(1d-3,POSNEG*(XSTOP-XS))
         HST=SIGN(HST,POSNEG)
         IFLAGST = 1
         ITST = 0
         ITSTMAX = 10
         DO WHILE ((itst .LT. itstMAX) .AND. (IFLAGST .EQ.1))
         itst=itst+1

C removed ierr
         IERR = 0
         CALL FCN(N,X,Y,K1,RPAR,IPAR)
         IF (IERR .NE. 0 ) GOTO 333

         DO 220 I=1,N
  220      Y1(I)=Y(I)+HST*A21*K1(I)

C removed ierr
         IERR = 0
         CALL FCN(N,X+C2*HST,Y1,K2,RPAR,IPAR)
         IF (IERR .NE. 0 ) GOTO 333
         DO 230 I=1,N
  230      Y1(I)=Y(I)+HST*(A31*K1(I)+A32*K2(I))
C removed ierr
         IERR = 0
         CALL FCN(N,X+C3*HST,Y1,K3,RPAR,IPAR)
         IF (IERR .NE. 0 ) GOTO 333
         DO 240 I=1,N
  240     Y1(I)=Y(I)+HST*(A41*K1(I)+A42*K2(I)+A43*K3(I))
C removed ierr
         IERR = 0
         CALL FCN(N,X+C4*HST,Y1,K4,RPAR,IPAR)
         IF (IERR .NE. 0 ) GOTO 333
         DO 250 I=1,N
  250     YSTI(I)=Y(I)+HST*(A51*K1(I)+A52*K2(I)+A53*K3(I)+A54*K4(I))
C removed ierr
         IERR = 0
         CALL FCN(N,X+C5*HST,YSTI,K5,RPAR,IPAR)
         IF (IERR .NE. 0 ) GOTO 333
         DO 260 I=1,N
         Y1(I)=Y(I)+HST*(A61*K1(I)+A62*K2(I)+A63*K3(I)+A64*K4(I)+
     +   A65*K5(I))
  260    continue

C removed ierr
         IERR = 0
         CALL FCN(n,x+c6*hst,y1,k6,rpar,ipar)
         IF (IERR .NE. 0 ) GOTO 333
         do 2770 i=1,n
         y1(i)=y(i)+Hst*(a71*k1(i)+a72*k2(i)+a73*k3(i)+a74*k4(i)
     +   +a75*k5(i)+a76*k6(i))
 2770   continue

         DO I=1,N
           etave(I) =( (y1(I)) - ysti(I) )
         END DO
c  two norm of etave
         ENORMVE = 0.0d0
         DO I=1,N
          ENORMVE = MAX(ENORMVE,ABS(ETAVE(I)))
         END DO

         IF ( ENORMVE .GT. 0.0d0 .AND.
     +         (.NOT. (.NOT. ENORMVE .GT. 0.0d0))
     +     .AND.  (.NOT. ENORMVE .GT. 1d300)  ) THEN
           DO I=1,N
            ETAVE(I) = ETAVE(I)/ENORMVE
           END DO
           IFLAGST = 0
         ELSE IF  (.NOT. (.NOT. ENORMVE .GT. 0.0d0 )) THEN
cF we try with a bigger stepsize
           HST = 2*HST
         ENDIF
cF errors in function evaluation we reduce the stepsize
  333    IF (IERR .NE.0  .OR. .NOT. (ENORMVE .GT. 0.0d0)
     +                .OR. ENORMVE .GT. 1d300 )
     +             HST = HST/10.0d0
          NFCN=NFCN+6
       END DO
c ERROR ESTIMATION FOR THE COMPUTATION OF THE INITIAL STEPSIZE

C removed ierr
       IERR = 0
C Karline: CORRECTED function call fcn(x, ...
       call fcn(N,x+hst,y1,k2,rpar,ipar)
       NFCN=NFCN+1
       do 280 i=1,n
       k4(i)=(E1*k1(i)+E3*k3(i)+E4*k4(i)+e5*k5(i)+e6*k6(i))*hst
  280  continue

C --- ERROR ESTIMATION
        ERR=0.D0
        IF (ITOL.EQ.0) THEN
          DO  I=1,N
           SK=ATOLI+RTOLI*max(abs(y(i)),abs(y1(i)))
           ERR=ERR+(K4(I)/SK)**2.0d0
          END DO
        ELSE
          DO  I=1,N
           SK=ATOL(I)+RTOL(I)*MAX(ABS(Y(I)),ABS(Y1(I)))
           ERR=ERR+(K4(I)/SK)**2.0d0
          END DO
        END IF
        ERR=SQRT(ERR/N)
      ENDIF


       HMAX=ABS(HMAX)
       IORD=5

       IF (H.EQ.0.D0) H=HINITCKSTIFF(N,X,Y,XEND,POSNEG,K1,K2,Y1,IORD,
     &                HMAX,HST,ERR,EXPO1,ATOL,RTOL,ITOL,RPAR,IPAR)


      if (iprint .gt. 0 ) then
        CALL Rprintd1('Initial step h = ', H)
        end if
cF end computation of eta for stiffness detection


CF normalization of etave using the input tolerances
C
       SCAL_ETA = 1d-12
       SCAL_ATOL = 1d-4
cF two norm of y0
       ENORMY = 0.0d0
       DO I=1,N
         ENORMY = ENORMY + ABS(Y(I))**2.0d0
       END DO
       ENORMY = SQRT(ENORMY)
       ENORMVE = 0.0d0
       IF (ITOL.EQ.0) THEN
         IF (ENORMY .EQ. 0 ) THEN
           SK=MAX(ATOLI,SCAL_ETA)
         ELSE
            SK=MAX(RTOLI*ENORMY,SCAL_ETA)
         END IF
         DO  I=1,N
           YS(I) = ETAVE(I)*SK
           ENORMVE=ENORMVE+ABS(YS(I))**2.0d0
         END DO
       ELSE
         SK = SCAL_ETA
         IF (ENORMY .EQ. 0 ) THEN
            DO I=1,N
              SK=MAX(ATOL(I),SCAL_ETA)
            END DO
         ELSE
           DO I=1,N
             SK=MAX(RTOL(I)*ENORMY,SK)
           END DO
         END IF
         DO  I=1,N
          YS(I) = ETAVE(I)*SK
          ENORMVE=ENORMVE+ABS(YS(I))**2.0d0
         END DO
       END IF
       ENORMVE = SQRT(ENORMVE)

CF end normalization of eta
c initial value for perturbed solution
       DO I=1,N
          YS(I) = Y(I)+YS(I)
       END DO

C removed ierr
       IERR = 0
       CALL FCN(N,X,YS,K1S,RPAR,IPAR)
       NFCN=NFCN+1
       IF (IERR .NE. 0  .OR.
     *             (.NOT. ENORMVE .GT. 0.0d0)  ) THEN
        CALL Rprint
     *   ('Stiffness detection based on conditioning cannot be used')
          NSTIFFCOND = 0

       ELSE
cf  computation of the relative difference y-ys
        ERRZ = 0.0d0
        IF (ITOL.EQ.0) THEN
         DO I=1,N
          SK = ATOLI*SCAL_ATOL+RTOLI*ABS(Y(I))
          ERRZ = ERRZ + (ABS(Y(I)-YS(I))/SK)**2.0d0
         END DO
        ELSE
         DO I=1,N
          SK=ATOL(I)*SCAL_ATOL+RTOL(I)*(ABS(Y(I)))
          ERRZ = ERRZ + (abs(Y(I)-YS(I))/SK)**2.0d0
         END DO
        END IF
        ERRZ = SQRT(ERRZ/N)
c       write(*,*) 'ERRZ ', ERRZ
cf end of computation of the relative difference y-ys
cf initialization of conditioning parameters
        KAPPA   = ENORMVE
        ENORMZ1 = ENORMVE
        GAMMA   = 0d0
        SIGMA   = 0d0
        SIGMATOT = 0.0d0
        XST = X
        NRESTART = 0
CF Parameter for testing stiffness
        NST   = 0
        NEST  = 0
        NNST  = 0
        NNEST = 0
        MAXnst = 20
        MAXIaSti = 15
        sigmamMAX = 50d0
        sigmaMAX = 100d0
        nwarn1 = 0
        nwarn2 = 0
        nwarn3 = 0
        nwarn4 = 0
      END IF
cF end if id nstiffcond .gt. 0
      END IF
      IF (.NOT. NSTIFFCOND .GT. 0 ) THEN
cF stiffness is not used
C removed ierr
       IERR = 0
       CALL FCN(N,X,Y,K1,RPAR,IPAR)
       HMAX=ABS(HMAX)
       IORD=5
       IF (H.EQ.0.D0) H=HINITck(N,FCN,X,Y,XEND,POSNEG,K1,K2,K3,IORD,
     &                       HMAX,ATOL,RTOL,ITOL,RPAR,IPAR)
       NFCN=NFCN+2
      END IF

      IF (NSTIFFERR .GT.0) THEN
cF initialization of parameter for stiffness based on error
        NERRTA = 0
        NNERRTA = 0
      END IF

      REJECT=.FALSE.
      XOLD=X

      IF (IOUT.NE.0) THEN
          IRTRN=1
          HOUT=H
          CALL SOLOUT(NACCPT+1,XOLD,X,Y,N,CONT,ICOMP,NRD,
     &                RPAR,IPAR,IRTRN)
          IF (IRTRN.LT.0) GOTO 79
      ELSE
          IRTRN=0
      END IF

C --- BASIC INTEGRATION STEP
   1  CONTINUE

      IF (NSTEP.GT.NMAX) GOTO 78
      IF (0.1D0*ABS(H).LE.ABS(X)*UROUND)GOTO 77
      IF ((X+1.01D0*H-XEND)*POSNEG.GT.0.D0) THEN
         H=XEND-X
         LAST=.TRUE.
      END IF
      NSTEP=NSTEP+1
C --- THE FIRST 6 STAGES
      IF (IRTRN.GE.2) THEN
C removed ierr
         IERR = 0
         CALL FCN(N,X,Y,K1,RPAR,IPAR)
         NFCN=NFCN+1
         IF (IERR .NE. 0 ) GOTO 444
      END IF

      DO 22 I=1,N
  22  Y1(I)=Y(I)+H*A21*K1(I)
C removed ierr
      IERR = 0
      CALL FCN(N,X+C2*H,Y1,K2,RPAR,IPAR)
      IF (IERR .NE. 0 ) GOTO 444
      DO 23 I=1,N
  23  Y1(I)=Y(I)+H*(A31*K1(I)+A32*K2(I))
C removed ierr
      IERR = 0
      CALL FCN(N,X+C3*H,Y1,K3,RPAR,IPAR)
      IF (IERR .NE. 0 ) GOTO 444
      DO 24 I=1,N
  24  Y1(I)=Y(I)+H*(A41*K1(I)+A42*K2(I)+A43*K3(I))
C removed ierr
      IERR = 0
      CALL FCN(N,X+C4*H,Y1,K4,RPAR,IPAR)
      IF (IERR .NE. 0 ) GOTO 444
      DO 25 I=1,N
  25  YSTI(I)=Y(I)+H*(A51*K1(I)+A52*K2(I)+A53*K3(I)+A54*K4(I))
C removed ierr
      IERR = 0
      CALL FCN(N,X+C5*H,YSTI,K5,RPAR,IPAR)
      IF (IERR .NE. 0 ) GOTO 444
      DO 26 I=1,N
      Y1(I)=Y(I)+H*(A61*K1(I)+A62*K2(I)+A63*K3(I)+A64*K4(I)+
     +   A65*K5(I))
 26    continue
C removed ierr
      IERR = 0
      call fcn(n,x+c6*h,y1,k6,rpar,ipar)
      IF (IERR .NE. 0 ) GOTO 444
      do 277 i=1,n
      y1(i)=y(i)+h*(a71*k1(i)+a72*k2(i)+a73*k3(i)+a74*k4(i)
     +  +a75*k5(i)+a76*k6(i))
 277   continue

c      IF (NSTIFFERR .GT. 0) THEN
       do 704 i=1,n
        vsti(i)=h*(ad1*k1(i)+ad2*k2(i)+ad3*k3(i)+ad4*k4(i)
     +   +ad5*k5(i)+ad6*k6(i))
 704   continue
c      END IF

      XPH=X+h
C removed ierr
      IERR = 0
C KARLINE: was (x, => (n,
C      call fcn(x,xph,y1,k2,rpar,ipar)
      call fcn(n,xph,y1,k2,rpar,ipar)
      NFCN=NFCN+6
      IF (IERR .NE. 0 ) GOTO 444

C --- THE FIRST 6 STAGES FOR YS SOLUTION
C --- with initial condition y0+etave
      IF (NSTIFFCOND .GT.0) THEN

       DO  I=1,N
         Y1S(I)=YS(I)+H*A21*K1S(I)
       END DO

C removed ierr
       IERR = 0
       CALL FCN(N,X+C2*H,Y1S,K2S,RPAR,IPAR)
       IF (IERR .NE. 0 ) GOTO 444
       DO  I=1,N
           Y1S(I)=YS(I)+H*(A31*K1S(I)+A32*K2S(I))
       END DO
C removed ierr
       IERR = 0
       CALL FCN(N,X+C3*H,Y1S,K3S,RPAR,IPAR)
       IF (IERR .NE. 0 ) GOTO 444
       DO  I=1,N
         Y1S(I)=YS(I)+H*(A41*K1S(I)+A42*K2S(I)+A43*K3S(I))
       END DO
C removed ierr
       IERR = 0
       CALL FCN(N,X+C4*H,Y1S,K4S,RPAR,IPAR)
       IF (IERR .NE. 0 ) GOTO 444
       DO  I=1,N
          YSTIS(I)=YS(I)+H*(A51*K1S(I)+A52*K2S(I)+A53*K3S(I)+A54*K4S(I))
       END DO
C removed ierr
       IERR = 0
       CALL FCN(N,X+C5*H,YSTIS,K5S,RPAR,IPAR)
       IF (IERR .NE. 0 ) GOTO 444
       DO  I=1,N
       Y1S(I)=YS(I)+H*(A61*K1S(I)+A62*K2S(I)+A63*K3S(I)+A64*K4S(I)+
     +   A65*K5S(I))
       END DO
C removed ierr
       IERR = 0
       CALL fcn(n,x+c6*h,y1s,k6s,rpar,ipar)
       IF (IERR .NE. 0 ) GOTO 444
       do  i=1,n
       y1s(i)=ys(i)+h*(a71*k1s(i)+a72*k2s(i)+a73*k3s(i)+a74*k4s(i)
     +  +a75*k5s(i)+a76*k6s(i))
       END DO
C removed ierr
       IERR = 0
C Karline:corrected function call (fnc(x,...
       CALL fcn(N,xph,y1s,k2s,rpar,ipar)
       NFCN=NFCN+6
       IF (IERR .NE. 0 ) GOTO 444
      END IF
c end computation of basic step for the second solution
c with different initial condition

 444  IF (IERR .NE.0 ) THEN
cF we reduce the stepsize
        ERR = 10.0d0
      ELSE
cf set up continuous output
        IF (IOUT.GE.2)then
          IF (HERMITE) THEN
            DO 40 J=1,NRD
             I=ICOMP(J)
             CONT(4*NRD+J)=H*(D1*K1(I)+D3*K3(I)+D4*K4(I)+D5*K5(I)
     &                   +D6*K6(I)+D7*K2(I))
  40         CONTINUE
          ELSE
            DO  J=1,NRD
             I=ICOMP(J)
             CONT(J) = Y(I)
             CONT(NRD+J) = H*K1(I)
             CONT(2*NRD+J) = H*K3(I)
             CONT(3*NRD+J) = H*K4(I)
             CONT(4*NRD+J) = H*K5(I)
             CONT(5*NRD+J) = H*K6(I)
             CONT(6*NRD+J) = H*K2(I)
            END DO
          END IF
        END IF

C ------- STIFFNESS DETECTION  based on approximation of eigenvalues
         IF (NSTIFF .GT. 0) THEN
            STNUM=0.D0
            STDEN=0.D0
            DO 64 I=1,N
               STNUM=STNUM+(K2(I)-K5(I))**2
               STDEN=STDEN+(Y1(I)-YSTI(I))**2
 64         CONTINUE
            IF (STDEN.GT.0.D0) HLAMB=H*SQRT(STNUM/STDEN)

 188   format(1x,'The first est',i5,3g22.10)
            IF (HLAMB.GT.3.25D0 .OR. .NOT.(HLAMB .GT. 0.0d0) ) THEN
               NONSTI=0
               IASTI=IASTI+1
               IF (IASTI.EQ.15) THEN
C KARLINE: always print ...
C                  IF (IPRINT.GT.0 .OR. NSTIFF .EQ. 1) THEN
        CALL Rprintd1(
     &               'The problem seems to become stiff at x = ',X)
        CALL Rprint('using the standard check of eignevalues')
                    IF (NSTIFF .EQ. 2) GOTO 76
C                  END IF
               END IF
            ELSE
               NONSTI=NONSTI+1
               IF (NONSTI.EQ.6) IASTI=0
            END IF
         END IF

         IF (NSTIFFCOND .GT. 0) THEN
cF approximation of the eiganvalues using ys
            STNUM=0.D0
            STDEN=0.D0
            DO  I=1,N
               STNUM=STNUM+(K2(I)-K5(I)-K2S(I)+K5S(I))**2
               STDEN=STDEN+(Y1(I)-YSTI(I)-Y1S(I)+YSTIS(I))**2
            END DO
            IF (STDEN.GT.0.D0) HLAMB=H*SQRT(STNUM/STDEN)

            IF ((HLAMB .GT. 2.8D0  .AND. HLAMB .LT. 4.2d0)
     +               .OR. .NOT.(HLAMB .GT. 0.0d0) ) THEN
               NONSTIS=0
               IASTIS=IASTIS+1
               IF (IASTIS.EQ.15) THEN
C Karline: always print
C                  IF (IPRINT.GT.0 .OR. NSTIFFCOND.EQ. 1) THEN
        CALL Rprintd1(
     &               'The problem seems to become stiff at x = ',X)
        CALL Rprint('using the check of eigenvalues based on YS')
C                  END IF
                  IF (NSTIFFCOND .EQ. 2) GOTO 76
               END IF
            ELSE
               NONSTIS=NONSTIS+1
               IF (NONSTI.EQ.6) IASTIS=0
            END IF

         END IF

c  End of first stiffness detection looking for largest value of
c   eigenvalue of the jacobian.
C     We jump here if there is no testing
        IF (NSTIFFCOND .EQ. 0) THEN
         do 28 i=1,n
          k4(i)=(E1*k1(i)+E3*k3(i)+E4*k4(i)+e5*k5(i)+e6*k6(i))*h
   28    continue
        ELSE
cf estimation of the error using ys


         do  i=1,n
           k4i = k4(i)
           k4(i)=(E1*k1(i)+E3*k3(i)+E4*k4(i)+e5*k5(i)+e6*k6(i))*h
           YSTI(i)=(E1*k1s(i)+E3*k3s(i)+E4*k4s(i)+e5*k5s(i)+e6*k6s(i))*h
           YSTIS(i)=(E1*(k1(i)-k1s(i))+E3*(k3(i)-k3s(i))+E4*(k4i-k4s(i))
     &      +e5*(k5(i)-k5s(i))+e6*(k6(i)-k6s(i)) )*h
         end do
        END IF



C --- ERROR ESTIMATION
        ERRy=0.D0
        errta=0.0D+0
        IF (ITOL.EQ.0) THEN
          DO 41 I=1,N
           SK=ATOLI+RTOLI*max(abs(y(i)),abs(y1(i)))
           errta=errta+(VSTI(I)/sk)**2.0d0
           ERRy=ERRy+(K4(I)/SK)**2.0d0
 41      continue
        ELSE
          DO 42 I=1,N
           SK=ATOL(I)+RTOL(I)*MAX(ABS(Y(I)),ABS(Y1(I)))
           errta=errta+(VSTI(I)/sk)**2.0d0
  42       ERRy=ERRy+(K4(I)/SK)**2.0d0
        END IF
        ERRy=SQRT(ERRy/N)
        errta=sqrt(errta/n)
        ERR = ERRy

        IF (NSTIFFERR .GT. 0) THEN
          IF (ERRta .LT. 0.1d0*ERRY) THEN
            Nerrta = Nerrta+1
            NNerrta = 0
          ELSE
            NNerrta = NNerrta+1
            IF (NNerrta .EQ. 6) THEN
               Nerrta = 0
            END IF
          END IF
          IF (Nerrta .EQ. 15) THEN
        CALL Rprinti1(
     &               'The problem seems to become stiff at x = ',X)
        CALL Rprint('using stiffness detection based on error estimate')
C            END IF
            IF (NSTIFFERR .EQ. 2) THEN
               IDID = -5
               GOTO 760
            END IF
          END IF
        END IF


C --- ERROR ESTIMATION USING  YS
       IF (NSTIFFCOND .GT. 0 ) THEN
         ERRys=0.D0
         erryys=0.0D+0
         IF (ITOL.EQ.0) THEN
           DO  I=1,N
            SK=ATOLI+RTOLI*max(abs(ys(i)),abs(y1s(i)))
            errys=errys+(YSTI(I)/sk)**2.0d0
            SK=ATOLI*SCAL_ATOL +
     &           RTOLI*max(ABS(Y(I)-YS(I)),ABS(Y1(I)-Y1S(I)))
            ERRyys=ERRyys+(YSTIS(I)/SK)**2.0d0
           END DO
         ELSE
          DO  I=1,N
           SK=ATOL(I)+RTOL(I)*MAX(ABS(YS(I)),ABS(Y1S(I)))
           errys=errys+(YSTI(I)/sk)**2.0d0
           SK=SCAL_ATOL+
     &         RTOL(I)*max(ABS(Y(I)-YS(I)),ABS(Y1(I)-Y1S(I)))
           ERRyys=ERRyys+(YSTIS(I)/SK)**2.0d0
          END DO
         END IF
         ERRys=SQRT(ERRys/N)
         erryys=sqrt(erryys/n)
c the error is the maximum of the three error estimation
         IF (ERRys .GT. SCAL_ETA) THEN
            ERR = max(ERRy,ERRys,ERRyys)
         ELSE
            ERR = ERRy
         END IF

       END IF
cF end if 444
      END IF

C     We jump here if there is no testing

C --- COMPUTATION OF HNEW
      FAC11=ERR**EXPO1
C --- LUND-STABILIZATION
      FAC=FAC11/FACOLD**BETA
C --- WE REQUIRE  FAC1 <= HNEW/H <= FAC2
      FAC=MAX(FACC2,MIN(FACC1,FAC/SAFE))
      HNEW=H/FAC

      IF(ERR.LE.1.D0)THEN
C --- STEP IS ACCEPTED
         FACOLD=MAX(ERR,1.0D-4)
         NACCPT=NACCPT+1

c        This is the end of the second stiffness detector.


c        This is the end of the second stiffness detector
c        This is the end of the second stiffness detector
c set up continuous output
          IF (IOUT.GE.2 .AND. HERMITE) THEN
            DO 43 J=1,NRD
            I=ICOMP(J)
            YD0=Y(I)
            YDIFF=Y1(I)-YD0
            BSPL=H*K1(I)-YDIFF
            CONT(J)=Y(I)
            CONT(NRD+J)=YDIFF
            CONT(2*NRD+J)=BSPL
            CONT(3*NRD+J)=-H*K2(I)+YDIFF-BSPL
  43        CONTINUE
          END IF
c initialize the variables for the next step
          DO 44 I=1,N
            k1(i)=k2(i)
  44        Y(I)=Y1(I)

c initialize the variables related to ys for the next step
      IF (NSTIFFCOND .GT. 0) THEN
         DO  I=1,N
            k1s(i)=k2s(i)
            Ys(I)=Y1s(I)
         END DO

c computation of kappa gamma and sigma
         ENORMZ0 = ENORMZ1
         ENORMZ1 = 0.0d0
         DO I=1,N
           ENORMZ1 = ENORMZ1 + ABS(Y(I)-YS(I))**2.0d0
         END DO
         ENORMZ1 = SQRT(ENORMZ1)
         kappa = max( kappa, ENORMZ1)

         gamma = gamma + H*(ENORMZ0+ENORMZ1)/2.0d0

         sigma = (kappa/gamma)*(XPH-XST)


cf relation between erry and errys
         IF (ERRY .LT. 0.1d0*ERRYYS) THEN
            NST = NST+1
            NEST = 0
         ELSE
            NEST = NEST+1
            IF (NEST .EQ. 3) THEN
               NST = 0
            END IF
         END IF


cf  computation of errz the relative error between y and ys
        errzold = errz
        ERRZ = 0.0d0
        IF (ITOL.EQ.0) THEN
         DO  I=1,N
          SK=ATOLI*SCAL_ATOL+RTOLI*(ABS(Y(I)))
          ERRZ = ERRZ + (abs(Y(I)-YS(I))/SK)**2.0d0
         END DO
        ELSE
         DO I=1,N
           SK=ATOL(I)*SCAL_ATOL+RTOL(I)*(ABS(Y(I)))
          ERRZ = ERRZ + (abs(Y(I)-YS(I))/SK)**2.0d0
          END DO
        END IF
        ERRZ = SQRT(ERRZ/N)
cf end of computation of errz

        IF ((ABS(ERRZ-ERRZOLD) .LT. 1d-5).AND.(ERRZ .LT. 5d-3)) THEN
           NNST = NNST+1
           NNEST =0
        ELSE
            NNEST = NNEST+1
            IF (NNEST .EQ.6) THEN
               NNST = 0
            END IF
         END IF

c
cf  restarting
         IF (NNST .GT. 50 .AND. SIGMAtot .LT. 50) THEN

           nrestart = 1
           ENORMY = 0.0d0
           DO I=1,N
             ENORMY = ENORMY + Y(I)**2.0d0
           END DO
           ENORMY = SQRT(ENORMY)
           ENORMVE = 0.0d0
           IF (ITOL.EQ.0) THEN
             SK=MAX(ATOLI+RTOLI*ENORMY,SCAL_ETA)
             DO  I=1,N
               YS(I) = ETAVE(I)*SK
               ENORMVE=ENORMVE+(YS(I))**2.0d0
             END DO
           ELSE
             SK = SCAL_ETA
             DO I=1,N
               SK=MAX(ATOL(I)+RTOL(I)*ENORMY,SK)
             END DO
             DO  I=1,N
               YS(I) = ETAVE(I)*SK
               ENORMVE=ENORMVE+(YS(I))**2.0d0
             END DO
           END IF
           ENORMVE = SQRT(ENORMVE)
cF new starting vector for YS
           DO I=1,N
             YS(I) = Y(I) + YS(I)
           END DO
cF computation of the new value of ERRZ
           IERR = 0
           CALL FCN(N,XPH,YS,K1S,RPAR,IPAR)
           NFCN=NFCN+1
           NST = 0
           NNST = 0
            ERRZ = 0.0d0
            IF (ITOL.EQ.0) THEN
             DO  I=1,N
               SK=ATOLI*SCAL_ATOL+RTOLI*(ABS(Y(I)))
               ERRZ = ERRZ + (abs(Y(I)-YS(I))/SK)**2.0d0
             END DO
            ELSE
             DO I=1,N
               SK=ATOL(I)*SCAL_ATOL+RTOL(I)*(ABS(Y(I)))
               ERRZ = ERRZ + (abs(Y(I)-YS(I))/SK)**2.0d0
             END DO
            END IF

            ERRZ = SQRT(ERRZ/N)

cF initialization of the conditioning parameters
            ENORMZ1=ENORMVE
            KAPPA = ENORMVE
            GAMMA = 0d0
            XST = XPH
            SIGMATOT = MAX(SIGMA,SIGMATOT)

         END IF
         IF (NRESTART .EQ. 0) SIGMATOT=SIGMA

        IF (  ( nst .GT. MAXnst .OR. IaStiS .GT. MAXIaSti )
     &   .OR. ( (errz .LT. 1e-5) .AND. (sigmatot .GT. sigmamMAX) )
     &   .OR.  sigmatot .GT. sigmaMAX
     &   .OR. KAPPA/ENORMVE  .GT. 1d20 ) THEN

            IF  (sigmatot .GT. sigmaMAX) THEN
C Karline: always print if nwarn1=0
C              IF ( (IPRINT .GT. 0 .OR. NSTIFFCOND.EQ.1)
C     &          .AND. (NWARN1 .EQ. 0) ) THEN
               IF (NWARN1 .EQ. 0) THEN
        CALL Rprintd1(
     &               'The problem seems to become stiff at x = ',X)
        CALL Rprinti1('Sigma = ', SIGMATOT)
                 nwarn1 = 1
              END IF
               IDID=-6
               IF (NSTIFFCOND .EQ. 2)  GOTO 760
            endif
            IF ( (NST .GT. MAXnst .OR. IaStiS .GT. MAXIaSti)
     &          .AND. (NWARN2 .EQ. 0) ) THEN

c              IF (IPRINT .GT. 0 .OR. NSTIFFCOND.EQ.1 ) THEN
        CALL Rprintd1(
     &'The stepsize is restricted only by stability reason at x = ', X)
c              END IF
              nwarn2=1
              IDID=-7
              IF (NSTIFFCOND .EQ. 2)  GOTO 760
            endif

            if ((nst .GT. MAXnst .AND. kappa/ENORMVE .GT. 1.01)
     &          .AND. (NWARN3 .EQ. 0) ) THEN
c               IF (IPRINT .GT. 0  .OR. NSTIFFCOND.EQ.1 ) THEN
        CALL Rprintd2('Kappa > 1; x = , kappa = ', X, KAPPA/ENORMVE)
c               END IF
               NWARN3=1
               IDID=-8
               IF (NSTIFFCOND .EQ. 2)  GOTO 760
             endif
            if ( kappa/ENORMVE .GT. 1d20 .AND. (NWARN4.EQ.0)) THEN
c              IF (IPRINT .GT. 0 .OR. NSTIFFCOND.EQ.1 )  THEN
        CALL Rprintd1(
     &   'The numerical solution is unstable, kappa = ',KAPPA)
c               END IF
               NWARN4  = 1
             IDID=-9
             IF (NSTIFFCOND .EQ. 2)  GOTO 760

            endif


        endif

      END IF
c end computation of kappa gamma and sigma
           XOLD=X
           X=XPH
           IF (IOUT.NE.0) THEN
            HOUT=H
            CALL SOLOUT(NACCPT+1,XOLD,X,Y,N,CONT,ICOMP,NRD,
     &                  RPAR,IPAR,IRTRN)
            IF (IRTRN.LT.0) GOTO 79
         END IF
C ------- NORMAL EXIT
         IF (LAST) THEN
            H=HNEW
            IF (NSTIFFCOND .GT.0 ) THEN
             kappa = kappa/ENORMVE
             gamma =gamma/ENORMVE
             sigmatot=max(sigma,sigmatot)
             IF (IPRINT .GT. 0) THEN
        CALL Rprintd1('Time     = ',X)
        CALL Rprintd1('Kappa    = ',kappa/ENORMVE)
        CALL Rprintd1('Gamma    = ',(gamma/ENORMVE)/(XPH-XST))
        CALL Rprintd1('Sigma    = ',sigma)
        CALL Rprintd1('Sigmatot = ',sigmatot)
              END IF
            END IF
            IDID=1
            RETURN
         END IF
         IF(ABS(HNEW).GT.HMAX)HNEW=POSNEG*HMAX
         IF(REJECT)HNEW=POSNEG*MIN(ABS(HNEW),ABS(H))
         REJECT=.FALSE.
      ELSE
C --- STEP IS REJECTED
         HNEW=H/MIN(FACC1,FAC11/SAFE)
         REJECT=.TRUE.
         IF(NACCPT.GE.1)NREJCT=NREJCT+1
         LAST=.FALSE.
      END IF
      H=HNEW
      GOTO 1
C --- FAIL EXIT
  76  CONTINUE
      IF (NSTIFFCOND .GT.0 ) THEN
            kappa = kappa/ENORMVE
            gamma =gamma/ENORMVE
            sigmatot=max(sigma,sigmatot)
             IF (IPRINT .GT. 0) THEN
        CALL Rprintd1('Time     = ',X)
        CALL Rprintd1('Kappa    = ',kappa/ENORMVE)
        CALL Rprintd1('Gamma    = ',(gamma/ENORMVE)/(XPH-XST))
        CALL Rprintd1('Sigma    = ',sigma)
        CALL Rprintd1('Sigmatot = ',sigmatot)
         END IF
      END IF
      IDID=-4
      RETURN
  760 CONTINUE

      IF (NSTIFFCOND .GT.0 ) THEN
            kappa = kappa/ENORMVE
            gamma =gamma/ENORMVE
            sigmatot=max(sigma,sigmatot)
            IF (IPRINT .GT. 0) THEN
        CALL Rprintd1('Time     = ',X)
        CALL Rprintd1('Kappa    = ',kappa/ENORMVE)
        CALL Rprintd1('Gamma    = ',(gamma/ENORMVE)/(XPH-XST))
        CALL Rprintd1('Sigma    = ',sigma)
        CALL Rprintd1('Sigmatot = ',sigmatot)
         END IF
      END IF

      RETURN
  77  CONTINUE
      IF (IPRINT.GT.0) THEN
      CALL Rprintd1('Exit of CashKarp at x = ',X)
      END IF
      IF (IPRINT.GT.0) THEN
        CALL Rprintd1('Step size too small, h = ',H)
      END IF
      IDID=-3
      IF (NSTIFFCOND .GT.0 ) THEN
            kappa = kappa/ENORMVE
            gamma =gamma/ENORMVE
            sigmatot=max(sigma,sigmatot)
          IF (IPRINT .GT. 0) THEN
        CALL Rprintd1('Time     = ',X)
        CALL Rprintd1('Kappa    = ',kappa/ENORMVE)
        CALL Rprintd1('Gamma    = ',(gamma/ENORMVE)/(XPH-XST))
        CALL Rprintd1('Sigma    = ',sigma)
        CALL Rprintd1('Sigmatot = ',sigmatot)
         END IF
      END IF
      RETURN
  78  CONTINUE
      IF (IPRINT.GT.0) THEN
      CALL Rprintd1('Exit of CashKarp at x = ',X)
      END IF
      IF (IPRINT.GT.0) THEN
        CALL Rprinti1('More than nmax steps are needed: ',NMAX)
      END IF
      IDID=-2
      IF (NSTIFFCOND .GT.0 ) THEN
            kappa = kappa/ENORMVE
            gamma =gamma/ENORMVE
            sigmatot=max(sigma,sigmatot)
            IF (IPRINT .GT. 0) THEN
        CALL Rprintd1('Time     = ',X)
        CALL Rprintd1('Kappa    = ',kappa/ENORMVE)
        CALL Rprintd1('Gamma    = ',(gamma/ENORMVE)/(XPH-XST))
        CALL Rprintd1('Sigma    = ',sigma)
        CALL Rprintd1('Sigmatot = ',sigmatot)
            END IF
      END IF
      RETURN
  79  CONTINUE

      IF (NSTIFFCOND .GT.0 ) THEN
            kappa = kappa/ENORMVE
            gamma =gamma/ENORMVE
            sigmatot=max(sigma,sigmatot)
            IF (IPRINT .GT. 0) THEN
        CALL Rprintd1('Time     = ',X)
        CALL Rprintd1('Kappa    = ',kappa/ENORMVE)
        CALL Rprintd1('Gamma    = ',(gamma/ENORMVE)/(XPH-XST))
        CALL Rprintd1('Sigma    = ',sigma)
        CALL Rprintd1('Sigmatot = ',sigmatot)
            END IF
      END IF
      IF (IPRINT.GT.0) THEN
      CALL Rprintd1('Exit of CashKarp at x = ',X)
      END IF
      IDID=2
      RETURN
      END
C
      FUNCTION HINITck(N,FCN,X,Y,XEND,POSNEG,F0,F1,Y1,IORD,
     &                 HMAX,ATOL,RTOL,ITOL,RPAR,IPAR)
C ----------------------------------------------------------
C ----  COMPUTATION OF AN INITIAL STEP SIZE GUESS
C ----------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(N),Y1(N),F0(N),F1(N),ATOL(*),RTOL(*)
      DIMENSION RPAR(*),IPAR(*)
      EXTERNAL FCN
C ---- COMPUTE A FIRST GUESS FOR EXPLICIT EULER AS
C ----   H = 0.01 * NORM (Y0) / NORM (F0)
C ---- THE INCREMENT FOR EXPLICIT EULER IS SMALL
C ---- COMPARED TO THE SOLUTION
      DNF=0.0D0
      DNY=0.0D0
      ATOLI=ATOL(1)
      RTOLI=RTOL(1)
      IF (ITOL.EQ.0) THEN
        DO 10 I=1,N
        SK=ATOLI+RTOLI*ABS(Y(I))
        DNF=DNF+(F0(I)/SK)**2
  10    DNY=DNY+(Y(I)/SK)**2
      ELSE
        DO 11 I=1,N
        SK=ATOL(I)+RTOL(I)*ABS(Y(I))
        DNF=DNF+(F0(I)/SK)**2
  11    DNY=DNY+(Y(I)/SK)**2
      END IF
      IF (DNF.LE.1.D-10.OR.DNY.LE.1.D-10) THEN
         H=1.0D-6
      ELSE
         H=SQRT(DNY/DNF)*0.01D0
      END IF
      H=MIN(H,HMAX)
      IERR = 1.0d0
      DO WHILE (IERR .NE.0)
        H=SIGN(H,POSNEG)
C ---- PERFORM AN EXPLICIT EULER STEP
       DO 12 I=1,N
  12   Y1(I)=Y(I)+H*F0(I)
C removed ierr
       IERR = 0
       CALL FCN(N,X+H,Y1,F1,RPAR,IPAR)
       IF (IERR .NE.0) H=H/10.0d0
      END DO
C ---- ESTIMATE THE SECOND DERIVATIVE OF THE SOLUTION
      DER2=0.0D0
      IF (ITOL.EQ.0) THEN
        DO 15 I=1,N
        SK=ATOLI+RTOLI*ABS(Y(I))
  15    DER2=DER2+((F1(I)-F0(I))/SK)**2
      ELSE
        DO 16 I=1,N
        SK=ATOL(I)+RTOL(I)*ABS(Y(I))
  16    DER2=DER2+((F1(I)-F0(I))/SK)**2
      END IF
      DER2=SQRT(DER2)/H
C ---- STEP SIZE IS COMPUTED SUCH THAT
C ----  H**IORD * MAX ( NORM (F0), NORM (DER2)) = 0.01
      DER12=MAX(ABS(DER2),SQRT(DNF))
      IF (DER12.LE.1.D-15) THEN
         H1=MAX(1.0D-6,ABS(H)*1.0D-3)
      ELSE
         H1=(0.01D0/DER12)**(1.D0/IORD)
      END IF
      H=MIN(100*ABS(H),H1,HMAX)
      HINITck=SIGN(H,POSNEG)
      RETURN
      END


      FUNCTION HINITCKSTIFF(N,X,Y,XEND,POSNEG,F0,F1,Y1,IORD,
     &                 HMAX,HST,ERR,EXPO1,ATOL,RTOL,ITOL,RPAR,IPAR)
C ----------------------------------------------------------
C ----  COMPUTATION OF AN INITIAL STEP SIZE  USING THE
C ----  INFORMATION OBTAINED AFTER THE COMPUTATION OF ETA
C ----------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(N),Y1(N),F0(N),F1(N),ATOL(*),RTOL(*)
      DIMENSION RPAR(*),IPAR(*)

C ---- COMPUTE A FIRST GUESS
C ----   H = 0.01 * NORM (Y0) / NORM (F0)
      DNF=0.0D0
      DNY=0.0D0
      ATOLI=ATOL(1)
      RTOLI=RTOL(1)
      IF (ITOL.EQ.0) THEN
        DO 10 I=1,N
        SK=ATOLI+RTOLI*MAX(ABS(Y(I)),ABS(Y(I)))
        DNF=DNF+(F0(I)/SK)**2
  10    DNY=DNY+(Y(I)/SK)**2
      ELSE
        DO 11 I=1,N
        SK=ATOL(I)+RTOL(I)*MAX(ABS(Y(I)),ABS(Y(I)))
        DNF=DNF+(F0(I)/SK)**2
  11    DNY=DNY+(Y(I)/SK)**2
      END IF

      IF (DNF.LE.1.D-10.OR.DNY.LE.1.D-10) THEN
         H=1.0D-6
      ELSE
         H=SQRT(DNY/DNF)*0.01D0
      END IF
      H=MIN(H,HMAX)
      H=SIGN(H,POSNEG)


C ---- ESTIMATE THE SECOND DERIVATIVE OF THE SOLUTION
C ---- USING THE INFORMATION GIVEN BY THE COMPUTATION OF ETA
C ---- HST IS THE STEP USED FOR THE COMPUTATION OF ETA

      DER2=0.0D0
      IF (ITOL.EQ.0) THEN
        DO 15 I=1,N
        SK=ATOLI+RTOLI*MAX(ABS(Y(I)),ABS(Y(I)))
  15    DER2=DER2+((F1(I)-F0(I))/SK)**2
      ELSE
        DO 16 I=1,N
        SK=ATOL(I)+RTOL(I)*MAX(ABS(Y(I)),ABS(Y(I)))
  16    DER2=DER2+((F1(I)-F0(I))/SK)**2
      END IF


      DER2=SQRT(DER2)/HST
C ---- STEP SIZE IS COMPUTED SUCH THAT
C ----  H**IORD * MAX ( NORM (F0), NORM (DER2)) = 0.01
      DER12=MAX(ABS(DER2),SQRT(DNF))
      IF (DER12.LE.1.D-15) THEN
         H1=MAX(1.0D-6,ABS(H)*1.0D-3)
      ELSE
         H1=(0.01D0/DER12)**(1.D0/IORD)
      END IF

C --- THE STEPSIZE MUST BE GREATER THE 1.d-10
      IF (H1 .GT. 1.D-10) THEN
         H=MIN(100*ABS(H),H1,HMAX)
      ELSE
         H=MIN(  100*ABS(H),ABS(HST),HMAX)
      ENDIF


C --- WE USE THE ERROR ESTIMATION TO COMPUTE ANOTHER STEPSIZE
C      IF ERR .NE. NAN
      IF ( .NOT. ERR .GT. 0.0d0) THEN
        FAC=ERR**EXPO1
        FAC11 = FAC*10
        H2=ABS(HST)/FAC11

C --- THE FINAL STEP IS THE MINIMUM IF H2 > 1e-10
        IF (H2 .GT. 1e-10) H=MIN(H,H2,HMAX)
      END IF

      HINITCKSTIFF=SIGN(H,POSNEG)
      RETURN
      END
C
      FUNCTION CONTCK(II,X,CON,ICOMP,ND)
C ----------------------------------------------------------
C     THIS FUNCTION CAN BE USED FOR CONTINUOUS OUTPUT IN CONNECTION
C     WITH THE OUTPUT-SUBROUTINE FOR CASH-KARP. IT PROVIDES AN
C     APPROXIMATION TO THE II-TH COMPONENT OF THE SOLUTION AT X.
C ----------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION b1,b3,b4,b5,b6,b7
      DIMENSION CON(*),ICOMP(ND)
      LOGICAL HERMITE
      COMMON /CONTCKV/XOLD,H, HERMITE
      DOUBLE PRECISION  db31, db32,db33, db41,db42,db43,
     +                  db51,db52,db53, db61,db62,db63,
     +                  db71,db72,db73

      COMMON /CONCOEFF/ db31, db32,db33, db41,db42,db43,
     +                  db51,db52,db53, db61,db62,db63,
     +                  db71,db72,db73

C ----- COMPUTE PLACE OF II-TH COMPONENT
c      db31 = 500.0d0/161.0d0
c      db32 = -20000.0d0/4347.0d0
c      db33 = 2750.0d0/1449.0d0
c      db41 = 125.0d0/132.0d0
c      db42 = -625.0d0/594.0d0
c      db43 = 125.0d0/396.0d0
c      db51 = 15.0d0/28.0d0
c      db52 = -15.0d0/14.0d0
c      db53 = db51
c      db61 = - 6144.0d0/1771.0d0
c      db62 = 2048.0d0/253.0d0
c      db63 = - 7680.0d0/1771.0d0
c      db71 = 3.0d0/2.0d0
c      db72 = -4.0d0
c      db73 = 5.0d0/2.0d0

      I=0
      DO 5 J=1,ND
      IF (ICOMP(J).EQ.II) I=J
   5  CONTINUE
      IF (I.EQ.0) THEN
         CALL Rprinti1('No dense output available for comp. nr',II)
         RETURN
      END IF

      IF (hermite) THEN
       THETA=(X-XOLD)/H
       THETA1=1.D0-THETA
       CONTCK=CON(I)+THETA*(CON(ND+I)+THETA1*(CON(2*ND+I)+THETA*
     &            (CON(3*ND+I)+THETA1*CON(4*ND+I))))
      ELSE
       THETA=(X-XOLD)/H
       THETA2 = THETA*THETA
       b3 = THETA2*( db31 + THETA*(db32 + THETA*db33))
       b4 = THETA2*( db41 + THETA*(db42 + THETA*db43))
       b5 = THETA2*( db51 + THETA*(db52 + THETA*db53))
       b6 = THETA2*( db61 + THETA*(db62 + THETA*db63))
       b7 = THETA2*( db71 + THETA*(db72 + THETA*db73))
       b1 = THETA - (b3+b4+b5+b6+b7)


       CONTCK= CON(I) + ( b1*CON(ND+I) + b3*CON(2*ND+I) +
     &                     b4*CON(3*ND+I) + b5*CON(4*ND+I) +
     &                     b6*CON(5*ND+I) + b7*CON(6*ND+I) )

      END IF
      RETURN
      END

C
C
C ks: made this a subroutine  and loop over all states...
      SUBROUTINE CONTD5ck(NEQ,X,CON,ICOMP,ND,VAL)
C ----------------------------------------------------------
C     THIS FUNCTION CAN BE USED FOR CONTINUOUS OUTPUT IN CONNECTION
C     WITH THE OUTPUT-SUBROUTINE FOR DOPRI5. IT PROVIDES AN
C     APPROXIMATION TO THE SOLUTION AT X.
C ----------------------------------------------------------

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION b1,b3,b4,b5,b6,b7,Val(NEQ)
      DIMENSION CON(*),ICOMP(ND)
      LOGICAL HERMITE
      COMMON /CONTCKV/XOLD,H, HERMITE
      DOUBLE PRECISION  db31, db32,db33, db41,db42,db43,
     +                  db51,db52,db53, db61,db62,db63,
     +                  db71,db72,db73

      COMMON /CONCOEFF/ db31, db32,db33, db41,db42,db43,
     +                  db51,db52,db53, db61,db62,db63,
     +                  db71,db72,db73

      IF (hermite) THEN
        THETA=(X-XOLD)/H
        THETA1=1.D0-THETA
        DO I = 1, NEQ
          VAL(I)= CON(I)+THETA*(CON(ND+I)+THETA1*(CON(2*ND+I)+THETA*
     &            (CON(3*ND+I)+THETA1*CON(4*ND+I))))
        ENDDO
      ELSE
        THETA=(X-XOLD)/H
        THETA2 = THETA*THETA

        b3 = THETA2*( db31 + THETA*(db32 + THETA*db33))
        b4 = THETA2*( db41 + THETA*(db42 + THETA*db43))
        b5 = THETA2*( db51 + THETA*(db52 + THETA*db53))
        b6 = THETA2*( db61 + THETA*(db62 + THETA*db63))
        b7 = THETA2*( db71 + THETA*(db72 + THETA*db73))
        b1 = THETA - (b3+b4+b5+b6+b7)

        DO I = 1, NEQ
          VAL(I)= CON(I) + ( b1*CON(ND+I) + b3*CON(2*ND+I) +
     &                     b4*CON(3*ND+I) + b5*CON(4*ND+I) +
     &                     b6*CON(5*ND+I) + b7*CON(6*ND+I) )
        ENDDO
      END IF


      RETURN
      END


C
      SUBROUTINE CKC(C2,C3,C4,C5,C6,E1,E3,E4,E5,E6,
     &                    A21,A31,A32,A41,A42,A43,A51,A52,A53,A54,
     &                    A61,A62,A63,A64,A65,A71,A72,A73,A74,
     +                    A75,A76,
     &                    D1,D3,D4,D5,D6,D7
     +                    ,Ad1,Ad2,Ad3,Ad4,Ad5,Ad6)
C ----------------------------------------------------------
C     RUNGE-KUTTA COEFFICIENTS OF CASH KARP (1980)
C ----------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         DOUBLE PRECISION  db31, db32,db33, db41,db42,db43,
     +                  db51,db52,db53, db61,db62,db63,
     +                  db71,db72,db73

      COMMON /CONCOEFF/ db31, db32,db33, db41,db42,db43,
     +                  db51,db52,db53, db61,db62,db63,
     +                  db71,db72,db73
      DOUBLE PRECISION  D1,D3,D4,D5,D6,D7

      C2=0.2D0
      C3=0.3D0
      C4=0.6D0
      C5=1.0D+0
      C6=7.0D+0/8.0D+0
      A21=0.2D0
      A31=3.D0/40.D0
      A32=9.D0/40.D0
      A41=3.0D+0/10.0D+0
      A42=-9.0D+0/10.0D+0
      A43=6.0D+0/5.0D+0
      A51=-11.0D+0/54.0D+0
      A52=5.0D+0/2.0D+0
      A53=-70.0D+0/27.0D+0
      A54=35.0D+0/27.0D+0
      A61=1631.0D+0/55296.0D+0
      A62=175.0D+0/512.0D+0
      A63=575.0D+0/13824.0D+0
      A64=44275.0D+0/110592.0D+0
      A65=253.0D+0/4096.0D+0
c      A71=2825.0D+0/27648.0D+0
c      A72=0.0D+0
c      A73=18575.0D+0/48384.0D+0
c      A74=13525.0D+0/55296.0D+0
c      A75=277.0D+0/14336.0D+0
c      A76=1.0D+0/4.0D+0
cf the coefficients are of the order 5 method
      A71=37.0D+0/378.0D0
      A72 = 0.0D+0
      A73=250.0D+0/621.0d0
      A74=125.0D+0/594.0D0
      A75= 0.0d0
      A76=512.0D+0/1771.0D0
      E1=37.0D+0/378.0D+0 -2825.0D+0/27648.0D+0
      E2 = 0.0D+0
      E3=250.0D+0/621.0D+0 -18575.0D+0/48384.0D+0
      E4=125.0D+0/594.0D+0-13525.0D+0/55296.0D+0
      E5=-277.0D+0/14336.0D+0
      E6=512.0D+0/1771.0D+0-1.0D+0/4.0D+0
C ---- Coefficients for stiffness detection,page 21 of
C  Hairer-Wanner, Solving Ordinary Differential Equations 2,
C  second revised edition
      ad1=-0.08536D+0
      ad2=0.088d+0
      ad3=-0.0096D+0
      ad4=0.0052D+0
      ad5=0.00576D0
      ad6=-0.004d+0
C--- Coefficients for the continuous output
      db31 = 500.0d0/161.0d0
      db32 = -20000.0d0/4347.0d0
      db33 = 2750.0d0/1449.0d0
      db41 = 125.0d0/132.0d0
      db42 = -625.0d0/594.0d0
      db43 = 125.0d0/396.0d0
      db51 = 15.0d0/28.0d0
      db52 = -15.0d0/14.0d0
      db53 = db51
      db61 = - 6144.0d0/1771.0d0
      db62 = 2048.0d0/253.0d0
      db63 = - 7680.0d0/1771.0d0
      db71 = 3.0d0/2.0d0
      db72 = -4.0d0
      db73 = 5.0d0/2.0d0
C ---- DENSE OUTPUT OF SHAMPINE (1986)
C    D1 = 16*c^*_1 - 8*A71-2
C    DJ = 16*c^*_j - 8*A7j  j=3,4,5,6
C    D7 = 16*c^*_7 +2
C    c^*_j coefficients of the continuous polymonial
C    computed at 1/2
      D1 =-115.0d0/126.0d0
      D3 = 2750.0d0/1449.0d0
      D4 = 125.0d0/396.0d0
      D5 = 15.0d0/28.0d0
      D6 = -7680.0d0/1771.0d0
      D7 = 5.0d0/2.0d0


      RETURN
      END


