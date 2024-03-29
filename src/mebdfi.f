C*********************************************************************
C       MAIN mebdfi DRIVER
C*********************************************************************


      SUBROUTINE mebdfi(N,T0,HO,Y0,YPRIME,TOUT,TEND,MF,IDID,LWORK,
     +     WORK,LIWORK,IWORK,MBND,MAXDER,ITOL,RTOL,ATOL,RPAR,IPAR,
     +     pderv,resid,IERR)


c karline: LOUT removed, added npd to check bounds; BOUND CHECK FOR RWORK!
c ***********************************************************************
c ***********************************************************************
c    Written by T.J. Abdulla and J.R. Cash,
c    Department of Mathematics,
c    Imperial College,
c    London SW7 2AZ
c    England
c
c    t.abdulla@ic.ac.uk   or  j.cash@ic.ac.uk
c
c    The author would be pleased to receive any comments
c          good or bad!!
c ***********************************************************************
c ***********************************************************************
c
c
C     THIS IS THE SEPTEMBER 20th 1999 VERSION OF OVDRIV, A PACKAGE FOR
C     THE SOLUTION OF THE INITIAL VALUE PROBLEM FOR SYSTEMS OF
C     IMPLICIT DIFFERENTIAL ALGEBRAIC EQUATIONS
c     G(t,Y,Y')=0, Y=(Y(1),Y(2),Y(3),.....,Y(N)).
C     SUBROUTINE OVDRIV IS A DRIVER ROUTINE FOR THIS PACKAGE.
C
C                    REFERENCES
C
C     1.  J. R. CASH, THE INTEGRATION OF STIFF INITIAL VALUE PROBLEMS
C         IN O.D.E.S USING MODIFIED EXTENDED BACKWARD DIFFERENTIATION
C         FORMULAE, COMP. AND MATHS. WITH APPLICS., 9, 645-657, (1983).
C     2.  J.R. CASH AND S. CONSIDINE, AN MEBDF CODE FOR STIFF
C         INITIAL VALUE PROBLEMS, ACM TRANS MATH SOFTWARE, 142-158,
C         (1992).
C     3.  J.R. CASH, STABLE RECURSIONS WITH APPLICATIONS TO THE
C         NUMERICAL SOLUTION OF STIFF SYSTEMS, ACADEMIC PRESS,(1979).
C     4.  A.C. HINDMARSH, ODEPACK, A SYSTEMISED COLLECTION OF ODE
C         SOLVERS, in SCIENTIFIC COMPUTING, R.S. STEPLEMAN et. al.
C         (eds) North-Holland, AMSTERDAM, pp55-64 , (1983).
C     5.  E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
C         EQUATIONS II, STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS,
C         SPRINGER 1996, page 267.
C
C     ----------------------------------------------------------------
C     OVDRIV IS TO BE CALLED ONCE FOR EACH OUTPUT VALUE OF T, AND
C     IN TURN MAKES REPEATED CALLS TO THE CORE INTEGRATOR STIFF.
C
C     THE INPUT PARAMETERS ARE ..
C     N     =  THE NUMBER OF FIRST ORDER DIFFERENTIAL EQUATIONS.
C     T0    =  THE INITIAL VALUE OF T, THE INDEPENDENT VARIABLE
C              (USED ONLY ON THE FIRST CALL)
C     HO    =  THE NEXT STEP SIZE IN T (USED FOR INPUT ONLY ON THE
C              FIRST CALL)
C     Y0    =  A VECTOR OF LENGTH N CONTAINING THE INITIAL VALUES OF Y
C              (USED FOR INPUT ONLY ON FIRST CALL)
C     YPRIME   A VECTOR OF LENGTH N CONTAINING THE INITIAL VALUES OF
C              DY/DT
C     TOUT  =  THE VALUE OF T AT WHICH OUTPUT IS DESIRED NEXT.
C              INTEGRATION WILL NORMALLY GO SLIGHTLY BEYOND TOUT
C              AND THE PACKAGE WILL INTERPOLATE TO T = TOUT
C     TEND  =  END OF THE RANGE OF INTEGRATION.
C     MF    =  THE METHOD FLAG.  AT PRESENT MF=21,22,23 OR 24 IS
C              ALLOWED. THESE ARE EXTENDED BACKWARD DIFFERENTIATION
C              FORMULAE USING THE CHORD METHOD WITH ANALYTIC OR NUMERICAL
C              JACOBIAN FOR MF=21,22 RESPECTIVELY. MF=23/24 ARE  THE SAME
C              AS FOR 21/22 BUT THE JACOBIAN IS NOW BANDED.   THE USER
C              NEEDS TO SPECIFY SUBROUTINE PDERV IF MF=21 OR 23.
C     IDID   = THE INTEGER USED ON INPUT TO INDICATE THE TYPE OF CALL.
C              1   THIS IS THE FIRST CALL FOR THE PROBLEM.
C              0   THIS IS NOT THE FIRST CALL FOR THIS PROBLEM
C                  AND INTEGRATION IS TO CONTINUE.
C             -1   THIS IS NOT THE FIRST CALL FOR THE PROBLEM,
C                  AND THE USER HAS RESET N, RTOL, ATOL,H  AND/OR MF.
C              2   SAME AS 0 EXCEPT THAT TOUT HAS TO BE HIT
C                  EXACTLY (NO INTERPOLATION IS DONE).
C                  ASSUMES TOUT .GE. THE CURRENT T.
C              3   SAME AS 0 EXCEPT CONTROL RETURNS TO CALLING
C                  PROGRAM AFTER ONE STEP. TOUT IS IGNORED, UNTIL THE
C                  INTEGRATION REACHES TOUT OR BEYOND. IF IT PASSES TOUT
C                  THE PROGRAM INTERPOLATES THE SOLUTION VALUES AND
C                  RETURNS THE SOLUTION VALUE AT TOUT.
C              SINCE THE NORMAL OUTPUT VALUE OF IDID IS 0,
C              IT NEED NOT BE RESET FOR NORMAL CONTINUATION.
C              THE FIRST CALL TO THE DRIVER IS WITH IDID=1 AND FOR
C              A SUCCESSFUL STEP THE DRIVER RETURNS WITH IDID=1.THUS
C              THE CALL WITH IDID = 1 IS SIMPLY THE FIRST
C              INITIALISING STEP FOR THE CODE.  THE USER
C              THEN NEEDS TO CONTINUE WITH IDID=0,-1,2 OR 3 AS ABOVE.
C     LOUT   = THE LOGICAL OUTPUT CHANNEL FOR MESSAGE PASSING.
C              removed
C     MBND   = AN ARRAY OF DIMENSION 4 FOR USE WHEN THE NEWTON ITERATION
C              MATRIX IS BANDED.  IF THIS MATRIX HAS ML DIAGONALS
C              BELOW THE MAIN DIAGONAL AND MU DIAGONALS ABOVE THE
C              MAIN DIAGONAL THEN:
C              MBND(1) = ML
C              MBND(2) = MU
C              MBND(3) = MU + ML + 1
C              MBND(4) = 2*ML + MU + 1
C     MAXDER=  THE MAXIMUM ORDER IS MAXDER + 1.
C              THE VALUE OF MAXDER CANNOT EXCEED 7.  THIS IS THE
C              VALUE RECOMMENDED UNLESS IT IS BELIEVED THAT THERE
C              ARE SEVERE STABILITY PROBLEMS IN WHICH CASE MAXDER=3
C              OR 4 SHOULD BE TRIED INSTEAD.
C     ITOL  =  AN INDICATOR OF THE TYPE OF ERROR CONTROL. SEE
C              DESCRIPTION BELOW UNDER ATOL.
C     RTOL  =  A RELATIVE ERROR TOLERANCE PARAMETER. CAN BE EITHER A
C              SCALAR OR AN ARRAY OF LENGTH N.  SEE DESCRIPTION
C              BELOW UNDER ATOL.
C     ATOL  =  THE ABSOLUTE ERROR BOUND.
C              THE INPUT PARAMETERS ITOL, RTOL AND ATOL DETERMINE
C              THE ERROR CONTROL PERFORMED BY THE SOLVER.  THE
C              SOLVER WILL CONTROL THE VECTOR e = (e(i)) OF ESTIMATED
C              LOCAL ERRORS IN y ACCORDING TO AN INEQUALITY OF THE FORM
C                  RMS-NORM OF (e(i)/ewt(i)) .LE. 1
C              THE ROOT MEAN SQUARE NORM IS
C                   RMS-NORM(V) = SQRT((SUM v(i)**2)/N).  HERE
C                ewt = (ewt(i)) IS A VECTOR OF WEIGHTS WHICH MUST
C              ALWAYS BE POSITIVE, AND THE VALUES OF RTOL AND ATOL
C              SHOULD BE NON-NEGATIVE. IF ITOL = 1 THEN SINGLE STEP ERROR
C              ESTIMATES DIVIDED BY YMAX(I) WILL BE KEPT LESS THAN 1
C              IN ROOT-MEAN-SQUARE NORM.  THE VECTOR YMAX OF WEIGHTS IS
C              COMPUTED IN OVDRIV. INITIALLY YMAX(I) IS SET AS
C              THE MAXIMUM OF 1 AND ABS(Y(I)).  THEREAFTER YMAX(I) IS
C              THE LARGEST VALUE OF ABS(Y(I)) SEEN SO FAR, OR THE
C              INITIAL VALUE YMAX(I) IF THAT IS LARGER.
C              IF ITOL = 1 THE USER NEEDS TO SET ATOL = RTOL =
C              THE PRECISION REQUIRED.  THEN
C                         ewt(i) = RTOL(1)*YMAX(i)
C              IF ITOL IS GREATER THAN 1 THEN
C                 ewt(i) = rtol(i)*abs(y(i)) + atol(i)
C              THE FOLLOWING TABLE GIVES THE TYPES (SCALAR/ARRAY)
C              OF RTOL AND ATOL, AND THE CORRESPONDING FORM OF ewt(i)
C                  ITOL   RTOL      ATOL       ewt(i)
C                   2    SCALAR    SCALAR   rtol*abs(y(i))   + atol
C                   3    SCALAR    ARRAY    rtol*abs(y(i))   + atol(i)
C                   4    ARRAY     SCALAR   rtol(i)*abs(y(i))+ atol
C                   5    ARRAY     ARRAY    rtol(i)*abs(y(i))+ atol(i)
C              IF EITHER OF THESE PARAMETERS IS A SCALAR, IT NEED
C              NOT BE DIMENSIONED IN THE USER'S CALLING PROGRAM.
C     NIND1    = THE NUMBER OF VARIABLES OF INDEX 1,2,3 RESPECTIVELY.
C     NIND2,NIND3  THESE ARE SET IN IWORK(1),(2),(3).
C              THE EQUATIONS MUST BE DEFINED SO THAT THE INDEX 1
C              VARIABLES PRECEDE THE INDEX 2 VARIABLES WHICH IN
C              TURN PRECEDE THE INDEX 3 VARIABLES.
C     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS)
C              WHICH CAN BE USED FOR COMMUNICATION BETWEEN  THE USER'S
C              CALLING PROGRAM AND THE F AND PDERV SUBROUTINES.
C     IERR     IERR IS AN INTEGER FLAG WHICH IS ALWAYS EQUAL TO ZERO
C              ON INPUT.  SUBROUTINES F AND PDERV SHOULD ALTER
C              IERR ONLY IF ONE OF THEM ENCOUNTERS AN ILLEGAL OPERATION SUCH
C              AS THE SQUARE ROOT OF A NEGATIVE NUMBER OR EXPONENT
C              OVERFLOW. THE USER CAN THEN ALTER H AND CALL THE
C              SUBROUTINE AGAIN WITH IDID=-1 IF HE WISHES.
C
C     AFTER THE INITIAL CALL, IF A NORMAL RETURN OCCURED AND A NORMAL
C     CONTINUATION IS DESIRED, SIMPLY RESET TOUT AND CALL AGAIN.
C     ALL OTHER PARAMETERS WILL BE READY FOR THE NEXT CALL.
C     A CHANGE OF PARAMETERS WITH IDID = -1 CAN BE MADE AFTER
C     EITHER A SUCCESSFUL OR AN UNSUCCESSFUL RETURN.
C
C     THE OUTPUT PARAMETERS ARE..
C     T0    =  THE VALUE OF T WHICH RELATES TO THE CURRENT SOLUTION
C              POINT Y0()
C     HO    =  THE STEPSIZE H USED LAST, WHETHER SUCCESSFULLY OR NOT.
C     Y0    =  THE COMPUTED VALUES OF Y AT T = TOUT
C     YPRIME=  THE COMPUTED VALUES OF DY/DT AT T=TOUT.
C     TOUT  =  UNCHANGED FROM ITS INPUT VALUE.
C     IDID  =  INTEGER USED ON OUTPUT TO INDICATE RESULTS, WITH
C              THE FOLLOWING VALUES AND MEANINGS..
C
C      0   INTEGRATION WAS COMPLETED TO TOUT OR BEYOND.
C
C     -1   THE INTEGRATION WAS HALTED AFTER FAILING TO PASS THE
C          ERROR TEST EVEN AFTER REDUCING H BY A FACTOR OF
C          1.E10 FROM ITS INITIAL VALUE.
C
C     -2   AFTER SOME INITIAL SUCCESS, THE INTEGRATION WAS
C          HALTED EITHER BY REPEATED ERROR TEST FAILURES OR BY
C          A TEST ON RTOL/ATOL.  TOO MUCH ACCURACY HAS BEEN REQUESTED.
C
C     -3   THE INTEGRATION WAS HALTED AFTER FAILING TO ACHIEVE
C          CORRECTOR CONVERGENCE EVEN AFTER REDUCING H BY A
C          FACTOR OF 1.E10 FROM ITS INITIAL VALUE.
C
C     -4   IMMEDIATE HALT BECAUSE OF ILLEGAL VALUES OF INPUT
C          PARAMETERS.  SEE PRINTED MESSAGE.
C
C     -5   IDID WAS -1 ON INPUT, BUT THE DESIRED CHANGES OF
C          PARAMETERS WERE NOT IMPLEMENTED BECAUSE TOUT
C          WAS NOT BEYOND T.  INTERPOLATION AT T = TOUT WAS
C          PERFORMED AS ON A NORMAL RETURN.  TO TRY AGAIN,
C          SIMPLY CALL AGAIN WITH IDID = -1 AND A NEW TOUT.
C
C     -6   MAXIMUM ALLOWABLE NUMBER OF INTEGRATION STEPS EXCEEDED.
C          TO CONTINUE THE USER SHOULD RESET IWORK(14).
C
C
C     -7   STEPSIZE IS TOO SMALL (LESS THAN SQRT(UROUND)/100)
C
C
C     -11   INSUFFICIENT REAL WORKSPACE FOR THE INTEGRATION
C
C     -12   INSUFFICIENT INTEGER WORKSPACE FOR THE INTEGRATION
C
C
C     IN ADDITION TO OVDRIVE, THE FOLLOWING ROUTINES ARE PROVIDED
C     IN THE PACKAGE..
C
C     INTERP( - )   INTERPOLATES TO GET THE OUTPUT VALUES
C                   AT T=TOUT FROM THE DATA IN THE Y ARRAY.
C     STIFF( - )    IS THE CORE INTEGRATOR ROUTINE.  IT PERFORMS A
C                   SINGLE STEP AND ASSOCIATED ERROR CONTROL.
C     COSET( - )    SETS COEFFICIENTS FOR BACKWARD DIFFERENTIATION
C                   SCHEMES FOR USE IN THE CORE INTEGRATOR.
C     PSET( - )     COMPUTES AND PROCESSES THE NEWTON ITERATION
C                   MATRIX DG/DY + (1/(H*BETA))DG/DY'
C     DEC_mebdfi(-) PERFORMS AN LU DECOMPOSITION ON A MATRIX.
C     SOL_mebdfi(-) SOLVES LINEAR SYSTEMS A*X = B AFTER DEC_mebdfi
C                   HAS BEEN CALLED FOR THE MATRIX A
C     DGBFA ( - )   FACTORS A DOUBLE PRECISION BAND MATRIX BY
C                   ELIMINATION.
C     DGBSL ( - )   SOLVES A BANDED LINEAR SYSTEM A*x=b
C
C                   ALSO SUPPLIED ARE THE BLAS ROUTINES
C
C                   daxpy, dscal, idamax, ddot.
C
C
C     THE FOLLOWING ROUTINES ARE TO BE SUPPLIED BY THE USER AND
C                   SHOULD BE DECLARED AS EXTERNAL.
C
C     PDERV(T,Y,PD,N,YPRIME,MBND(4),CON,IPAR,RPAR,IERR)
C                         COMPUTES THE N*N NEWTON ITERATION MATRIX
C                         OF PARTIAL DERIVATIVES PD=DG/DY + c(DG/DY')
C                         WHERE c IS A SCALAR.  IF A NUMERICAL
C                         JACOBIAN IS REQUIRED (MF=22 OR 24) THEN THIS
C                         SUBROUTINE CAN BE A DUMMY ONE. THE ITERATION
C                         MATRIX IS STORED AS AN N BY N ARRAY IF THE
C                         MATRIX IS FULL.  IF THE ITERATION MATRIX IS
C                         BANDED THE ARRAY PD IS OF SIZE MBND(4)*N.
C                         IF THE ITERATION MATRIX IS FULL, PD(I,J) IS
C                         TO BE SET TO THE PARTIAL DERIVATIVE OF THE
C                         ith COMPONENT OF G(t,Y,Y') WITH
C                         RESPECT TO Y(J).  IF THE JACOBIAN IS BANDED
C                         WITH mu DIAGONALS ABOVE THE MAIN DIAGONAL
C                         THE PARTIAL DERIVATIVE OF THE ith COMPONENT
C                         OF G(t,Y,Y') WITH RESPECT TO Y(J) SHOULD BE
C                         PUT IN PD(i-j+mu+1,j). PDERV IS CALLED ONLY IF
C                         MITER = 1 OR 3.  OTHERWISE A DUMMY ROUTINE CAN
C                         BE SUBSTITUTED.
C
C KS: calling sequence made compatible with DASPK  - note: has NO IERR!
C     pderv(T, Y, YPRIME, PD, CJ, RPAR, IPAR)
C
C     THE DIMENSION OF PD  MUST BE AT LEAST N**2 FOR A FULL ITERATION
C     MATRIX AND MBND(4)*N FOR A BANDED ITERATION MATRIX
C     (MF=23 OR 24) .  THE DIMENSIONS
C     OF YMAX,ERROR,SAVE1,SAVE2,IPIV AND THE FIRST DIMENSION
C     OF Y SHOULD ALL BE AT LEAST N.
C
C     RESID(N,T,Y,DELTA,YPRIME,IPAR,RPAR,IERR)
C                          COMPUTES THE RESIDUAL VECTOR
C                          DELTA = G(T,Y,Y').
C KS: calling sequence made compatible with DASPK
C     resid(T,Y,Yprime, CJ, delta, ierr, rpar, ipar)
C
C
C     UROUND   THIS IS THE UNIT ROUNDOFF AND HAS TO BE SET AS
C                        UROUND = DLAMCH('Epsilon') - karline: has been removed
C     EPSJAC   = sqrt(UROUND).
C
C     HUSED  (=WORK(2))    LAST STEPSIZE SUCCESSFULLY USED BY THE INTEGRATOR
C     NQUSED (=IWORK(4))   LAST ORDER SUCCESSFULLY USED
C     NSTEP  (=IWORK(5))   NUMBER OF SUCCESSFUL STEPS TAKEN SO FAR
C     NFAIL  (=IWORK(6))   NUMBER OF FAILED STEPS
C     NRE    (=IWORK(7))   NUMBER OF RESIDUAL EVALUATIONS  SO FAR
C     NJE    (=IWORK(8))   NUMBER OF JACOBIAN EVALUATIONS  SO FAR
C     NDEC   (=IWORK(9))   NUMBER OF LU DECOMPOSITIONS  SO FAR
C     NBSOL  (=IWORK(10))   NUMBER OF 'BACKSOLVES'  SO FAR
C     NPSET  (=IWORK(11))   NUMBER OF TIMES A NEW COEFFICIENT MATRIX HAS BEEN
C               FORMED SO FAR
C     NCOSET (=IWORK(12))   NUMBER OF TIMES THE ORDER OF THE METHOD USED HAS
C               BEEN CHANGED SO FAR
C     MAXORD (=IWORK(13))   THE MAXIMUM ORDER USED SO FAR IN THE INTEGRATION
c     MAXSTP (=IWORK(14))   THE MAXIMUM ALLOWED NUMBER OF STEPS SET BY THE
C                           USER.
C
C    IF IT IS ONLY REQUIRED TO CONTROL THE ACCURACY IN THE
C    DIFFERENTIAL VARIABLES THEN THE USER SHOULD FIND THE
C    STRING 'AMMEND' AND MAKE CHANGES THERE
C

c

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     THIS SUBROUTINE IS FOR THE PURPOSE               *
C     OF SPLITTING UP THE WORK ARRAYS WORK AND IWORK   *
C     FOR USE INSIDE THE INTEGRATOR STIFF              *
C     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

C     .. SCALAR ARGUMENTS ..
      INTEGER IDID,LIWORK,LOUT,LWORK,MF,N,MAXDER,ITOL
C     ..
C     .. ARRAY ARGUMENTS ..
      DIMENSION  WORK(LWORK),Y0(N),RTOL(*),ATOL(*),RPAR(*),IPAR(*)
     +        ,YPRIME(N)
      INTEGER IWORK(LIWORK), MBND(4)
C     ..
C     .. LOCAL SCALARS ..
      INTEGER I1,I2,I3,I4,I5,I6,I7,I8,I9,I10, I11
C     COMMON BLOCKS
C     ..
C     .. EXTERNAL SUBROUTINES ..
      EXTERNAL OVDRIV,pderv,resid
C     ..
C FM added UROUND and EPSJAC in SAVE
C     .. SAVE STATEMENT ..
      SAVE  I1,I2,I3,I4,I5,I6,I7,I8,I9,I10, I11, UROUND, EPSJAC
C     ..

C KS: hard-coded LOUT = 0 (never used anymore...)
      LOUT = 0
      IF (IDID.EQ.1) THEN
         IF (N.LE.0) THEN
            CALL Rprinti1 ('Illegal value for number of equations ',N)
            call rexit('stopped')
            IDID = -4

         ELSE
            IF(MF.LT.23) MBND(4)=N
            I1 = N*12+ 3
            I2 = I1 + N*12
            I3 = I2 + N*2
            I4 = I3 + N
            I5 = I4 + N
            I6 = I5 + N
            I7 = I6 + N
            I8 = I7 + N
            I9 = I8 + N
            I10 = I9 + MBND(4)*N
            I11 = I10 + MBND(4)*N
c            UROUND = DLAMCH('Epsilon')     DID NOT WORK...
            UROUND = d1mach(3)
            WORK(1) = UROUND
            EPSJAC = SQRT(WORK(1))

c            IF (LWORK.LT.(I10+1)) THEN  KS: CHANGED THAT:
c Francesca added I11
           IF (LWORK.LT.(I11+1)) THEN
               IDID = -11
      CALL Rprint('Real workspace is insufficient ')
      CALL Rprinti1('Size of workspace must be at least ',I11 +1)
               call rexit('stopped')

            ENDIF

            IF (LIWORK.LT.N+14) THEN
               IDID = -12
      CALL Rprint('Integer workspace is insufficient ')
      CALL Rprinti1('Size of workspace must be at least ', N +14)
               call rexit('stopped')
            END IF


         END IF

         IF (IDID.LT.0) RETURN

      END IF

c
c    THE DIMENSION OF THE REAL WORKSPACE, WORK, HAS TO BE AT LEAST
c     (32 + 2*MBND(4))*N+2 WHILE THE DIMENSION OF THE INTEGER
c    WORKSPACE HAS TO BE AT LEAST N+14.
c

      IERR = 0
c
c     THE ERROR FLAG IS INITIALISED
      CALL OVDRIV(N,MBND(4),T0,HO,Y0,YPRIME,TOUT,TEND,MF,IDID,
     +    WORK(3),WORK(I1),WORK(I2),WORK(I3),WORK(I4),WORK(I5),WORK(I6),
     +     WORK(I7),WORK(I8),WORK(I9),WORK(I10),IWORK(15),
     +     MBND,IWORK(1),IWORK(2),IWORK(3),MAXDER,ITOL,RTOL,ATOL,RPAR,
     +     IPAR,pderv,resid,IWORK(4),IWORK(5),IWORK(6),IWORK(7),
     +     IWORK(8),IWORK(9),IWORK(10),IWORK(11),IWORK(12),IWORK(13),
     +     IWORK(14),WORK(1),WORK(2),EPSJAC,IERR)

C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     WORK() HOUSES THE FOLLOWING ARRAYS
C
C     Y(N,12)   , YHOLD(N,12) , YNHOLD(N,2) , YMAX(N)
C     ERRORS(N) , SAVE1(N) , SAVE2(N) , SCALE(N) , ARH(N) , PW(MBND(4)*N)
C     PWCOPY(MBND(4)*N)
C     IF THE BANDED OPTION IS NOT BEING USED THEN MBND(4)=N.
C     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      RETURN

      END
C--------------------------------------------------------------------------
C KS: added NPD to check pw and pwcopy
      SUBROUTINE OVDRIV(N,NPD,T0,HO,Y0,YPRIME,TOUT,TEND,MF,IDID,Y,
     +     YHOLD,YNHOLD,YMAX,ERRORS,SAVE1,SAVE2,SCALE,ARH,PW,PWCOPY,
     +     IPIV,MBND,NIND1,NIND2,NIND3,MAXDER,ITOL,RTOL,ATOL,RPAR,
     +     IPAR,pderv,resid,NQUSED,NSTEP,NFAIL,NRE,NJE,NDEC,NBSOL,
     +     NPSET,NCOSET,MAXORD,MAXSTP,UROUND,HUSED,EPSJAC,IERR)


      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     START OF PROGRAM PROPER
C
C     .. SCALAR ARGUMENTS ..
      INTEGER IDID,MF,N,NPD,MAXDER,ITOL
C     ..
C     .. ARRAY ARGUMENTS ..
      DIMENSION  ARH(N),ERRORS(N),PW(NPD*N),PWCOPY(NPD*N),SAVE1(N),
     +  SAVE2(N),Y(N,12),Y0(N),YHOLD(N,12),YMAX(N),YNHOLD(N,2),RTOL(*),
     +     ATOL(*),SCALE(N),RPAR(*),IPAR(*),YPRIME(N)
      INTEGER IPIV(N), MBND(4)
C     ..
C     .. LOCAL SCALARS ..
      INTEGER I,KGO,NHCUT
C     ..
C     .. EXTERNAL SUBROUTINES ..
      EXTERNAL INTERP,STIFF,pderv,resid
C     ..
C     .. INTRINSIC FUNCTIONS ..
      INTRINSIC DABS,DMAX1
C     ..
C     .. COMMON BLOCKS ..
      INTEGER JSTART,KFLAG,MAXORD,NBSOL,NCOSET,NDEC,
     +     NRE,NJE,NPSET,NQUSED,NSTEP
      SAVE T,H,HMIN,HMAX,KFLAG,JSTART
C     ..

      IF (IDID.EQ.0) THEN
C        I.E. NORMAL CONTINUATION OF INTEGRATION
         T0=T
CKS: hmax should become a parameter !
         HMAX = DABS(TEND-T0)*10.0D+0
         IF ((T-TOUT)*H.GE.0.0D+0) THEN
C           HAVE OVERSHOT THE OUTPUT POINT, SO INTERPOLATE
            CALL INTERP(N,JSTART,H,T,Y,TOUT,Y0)
            IDID = KFLAG
            T0 = TOUT
            HO = H
            RETURN

         END IF

      ELSE IF (IDID.EQ.2) THEN
C        I.E. CONTINUING INTEGRATION BUT WISH TO HIT TOUT
         T0 = T
         HMAX = DABS(TEND-T0)*10.0D+0
         IF (((T+H)-TOUT)*H.GT.0.0D+0) THEN
C           WE HAVE ALREADY OVERSHOT THE OUTPUT POINT OR WE WILL
C           DO SO ON THE NEXT STEP
            IF (((T-TOUT)*H.GE.0.0D+0) .OR. (DABS(T-TOUT).LE.
     +           100.0D+0*UROUND*HMAX)) THEN
C              HAVE OVERSHOT THE OUTPUT POINT, SO INTERPOLATE
               CALL INTERP(N,JSTART,H,T,Y,TOUT,Y0)
               T0 = TOUT
               HO = H
               IDID = KFLAG
               RETURN

            ELSE
C              WILL PASS TOUT ON NEXT STEP WITH CURRENT STEPSIZE
C              SO REDUCE STEPSIZE TO HIT TOUT 'EXACTLY'
               H = (TOUT-T)* (1.0D+0-4.0D+0*UROUND)
               JSTART = -1
            END IF

         END IF

      ELSE IF (IDID.EQ.-1) THEN
C        NOT FIRST CALL BUT PARAMETERS RESET
         H = HO
         IF(H.LT.EPSJAC/100.0D+0) THEN
      CALL Rprint('Stepsize is too small')
      call rexit('stopped')

            IDID = -7
            RETURN
         ENDIF
         T0 = T
         IF ((T-TOUT)*H.GE.0.0D+0) THEN
C           HAVE OVERSHOT TOUT
      CALL Rprintd1('IDID = -1 on input & (t-tout)*h .ge. 0. t= ', T)
      CALL Rprintd2('tout and h = ', TOUT, H)
      CALL Rprint(' interpolation was done as on normal return.')
      CALL Rprint(' desired parameter changes were not made.')
            call rexit('stopped')

            CALL INTERP(N,JSTART,H,T,Y,TOUT,Y0)
            HO = H
            T0 = TOUT
            IDID = -5
            RETURN

         ELSE
            JSTART = -1
         END IF

      ELSE IF (IDID.EQ.3) THEN
         T0 = T
         IF ((T-TOUT)*H.GE.0.0D+0) THEN
C           HAVE OVERSHOT,SO INTERPOLATE
            CALL INTERP(N,JSTART,H,T,Y,TOUT,Y0)
            IDID = KFLAG
            T0 = TOUT
            HO = H
            RETURN

         END IF

      ELSE
C        IDID SHOULD BE 1 AND THIS IS THE FIRST CALL FOR THIS PROBLEM
C        CHECK THE ARGUMENTS THAT WERE PASSED FOR CORRECTNESS
         IF (IDID.NE.1) THEN
C           VALUE OF IDID NOT ALLOWED
      CALL Rprinti1('Illegal input.. idid =', IDID)
            call rexit('stopped')

            IDID = -4
         END IF

         NN=N
         IF(ITOL.LE.3) NN = 1
         DO 1 I=1,NN
            IF (RTOL(I).LT.0.0D+0) THEN
C           ILLEGAL VALUE FOR RELATIVE ERROR TOLERENCE
      CALL Rprint('Illegal input.. rtol .le. 0.')
               call rexit('stopped')

               IDID = -4
            END IF
 1       CONTINUE

         NN=N
         IF(ITOL.EQ.1.OR.ITOL.EQ.2.OR.ITOL.EQ.4) NN=1
         DO 2 I=1,NN
            IF (ATOL(I).LT.0.0D+0) THEN
C           ILLEGAL ABSOLUTE ERROR TOLERANCE
      CALL Rprint('Illegal input.. atol .le. 0.')
               call rexit('stopped')

               IDID=-4
            ENDIF
 2       CONTINUE
         IF(ITOL.EQ.1.AND.RTOL(1).EQ.0) THEN
C           ILLEGAL ERROR TOLERANCE
      CALL Rprint('Illegal input.. rtol .le. 0.')
            call rexit('stopped')

            IDID = -4
         ENDIF



         IF(ITOL.NE.1) THEN
            VHOLD = 0.0D+0
            DO 3 I=1,N
               IF(ITOL.EQ.2) THEN
                  VHOLD = DMAX1(RTOL(1),ATOL(1))
               ELSE IF (ITOL.EQ.3) THEN
                  VHOLD = DMAX1(RTOL(1),ATOL(I))
               ELSE IF (ITOL.EQ.4) THEN
                  VHOLD = DMAX1(RTOL(I),ATOL(1))
               ELSE IF (ITOL.EQ.5) THEN
                  VHOLD = DMAX1(RTOL(I),ATOL(I))
               ENDIF
               IF(VHOLD.LE.0.0D+0) THEN
      CALL Rprint('Illegal input.. rtol .le. 0.')
                  call rexit('stopped')

                  IDID = -4
               ENDIF
 3          CONTINUE
         ENDIF

         IF (N.LE.0) THEN
C           ILLEGAL VALUE FOR THE NUMBER OF EQUATIONS
      CALL Rprint('Illegal input.. n .le. 0')
            call rexit('stopped')

            IDID = -4
         END IF

         IF ((T0-TOUT)*HO.GE.0.0D+0) THEN
C           PARAMETERS FOR INTEGRATION ARE ILLEGAL
      CALL Rprint('Illegal input.. (t0-tout)*h .ge. 0.')
            call rexit('stopped')

            IDID = -4
         END IF

         IF ((MF.NE.21).AND.(MF.NE.22).AND.(MF.NE.23).AND.(MF.NE.24))
     +   THEN
C           ILLEGAL VALUE FOR METHOD FLAG
      CALL Rprinti1('Illegal input.. method flag, mf, = ',MF)
      CALL Rprint  ('         allowed values are 21 or 22')

            IDID = -4
         END IF

         IF(ITOL.LT.1.OR.ITOL.GT.5) THEN
C           ILLEGAL VALUE FOR ERROR CONTROL PARAMETER
      CALL Rprint('Illegal value for itol')

            IDID=-4
         ENDIF

         IF(MAXDER.LT.1.OR.MAXDER.GT.7) THEN
C        ILLEGAL VALUE FOR MAXIMUM ORDER
      CALL Rprint('Illegal value for maxder')

            IDID = -4
         ENDIF

         IF(NIND1.EQ.0)  NIND1=N
         IF(NIND1 + NIND2 + NIND3 .NE. N) THEN
C    SUM OF VARIABLES OF DIFFERENT INDEX SHOULD BE N.
      CALL Rprint('Bad input for number of variables of index 1,2,3')

            IDID = -4
         ENDIF

         IF (IDID.NE.1) THEN
            RETURN

         ELSE
C           THE INITIAL PARAMETERS ARE O.K. SO INITIALISE EVERYTHING
C           ELSE NECESSARY FOR THE INTEGRATION.
C           IF VALUES OF YMAX OTHER THAN THOSE SET BELOW ARE DESIRED,
C           THEY SHOULD BE SET HERE. ALL YMAX(I) MUST BE POSITIVE. IF
C           VALUES FOR HMIN OR HMAX, THE BOUNDS ON DABS(H), OTHER THAN
C           THOSE BELOW ARE DESIRED, THEY SHOULD BE SET BELOW.

            IF(ITOL.EQ.1) THEN
               DO 10 I = 1,N
                  YMAX(I) = DABS(Y0(I))
                  YMAX(I)=DMAX1(YMAX(I),1.0D0)
 10            CONTINUE
            ENDIF
            DO 15 I=1,N
               Y(I,1)=Y0(I)
 15         CONTINUE
            T = T0
            H = HO
            HMIN = DABS(HO)
            HMAX = DABS(T0-TEND)*10.0D+0
            JSTART = 0
            NHCUT = 0
         END IF


      END IF
C     <<<<<<<<<<<<<<<<<
C     <  TAKE A STEP  >
C     <<<<<<<<<<<<<<<<<

 20   IF ((T+H).EQ.T) THEN
      CALL Rprint('Warning.. T + H = T on next step.')

      END IF


      CALL STIFF(H,HMAX,HMIN,JSTART,KFLAG,MF,MBND,
     +    NIND1,NIND2,NIND3,T,TOUT,TEND,Y,YPRIME,N,NPD,
     +    YMAX,ERRORS,SAVE1,SAVE2,SCALE,PW,PWCOPY,YHOLD,
     +    YNHOLD,ARH,IPIV,MAXDER,ITOL,RTOL,ATOL,RPAR,IPAR,
     +    pderv,resid,NQUSED,NSTEP,NFAIL,NRE,NJE,NDEC,NBSOL,
     +    NPSET,NCOSET,MAXORD,UROUND,EPSJAC,HUSED,IERR)

      KGO = 1 - KFLAG
      IF (KGO.EQ.1) THEN
C        NORMAL RETURN FROM STIFF
         GO TO 30

      ELSE IF (KGO.EQ.2) THEN
C        COULD NOT ACHIEVE REQUIRED PRECISION WITH HMIN
C        SO CHOP HMIN IF WE HAVEN'T DONE SO 10 TIMES
         GO TO 60

      ELSE IF (KGO.EQ.3) THEN
C        ERROR REQUIREMENT SMALLER THAN CAN BE HANDLED FOR THIS PROBLEM

      CALL Rprintd2('KFLAG = -2 at t  and h = ',T,H)
      CALL Rprint('The requested error is smaller than can be handled')

         GO TO 70

      ELSE IF (KGO.EQ.4) THEN
C        COULD NOT ACHIEVE CONVERGENCE WITH HMIN
      CALL Rprintd1('FLAG = -3 from integrator at t = ', T)
      CALL Rprint('corrector convergence could not be achieved')

         GO TO 60

      ELSE IF (KGO.EQ.6) THEN
C        PASSED TOUT
      CALL Rprintd2('Kflag = -5 at t and h = ',T, H)
      CALL Rprint('Overshoot Tout')

         GO TO 70

      ELSE IF (KGO.EQ.8) THEN
C        STEPSIZE TOO SMALL
      CALL Rprintd2('Kflag = -7 at t and h = ', T, H)
      CALL Rprint('Stepsize too small')

         GO TO 70

      END IF

 30   CONTINUE
C ---------------------------------------------------------------------
C     NORMAL RETURN FROM THE INTEGRATOR.
C
C     THE WEIGHTS YMAX(I) ARE UPDATED IF ITOL=1.
C     IF DIFFERENT VALUES ARE DESIRED, THEY SHOULD BE SET HERE.
C
C     ANY OTHER TESTS OR CALCULATIONS THAT ARE REQUIRED AFTER EVERY
C     STEP SHOULD BE INSERTED HERE.
C
C     IF IDID = 3, Y0 IS SET TO THE CURRENT Y VALUES ON RETURN.
C     IF IDID = 2, H IS CONTROLLED TO HIT TOUT (WITHIN ROUNDOFF
C     ERROR), AND THEN THE CURRENT Y VALUES ARE PUT IN Y0 ON RETURN.
C     FOR ANY OTHER VALUE OF IDID, CONTROL RETURNS TO THE INTEGRATOR
C     UNLESS TOUT HAS BEEN REACHED.  THEN INTERPOLATED VALUES OF Y ARE
C     COMPUTED AND STORED IN Y0 ON RETURN.
C     IF INTERPOLATION IS NOT DESIRED, THE CALL TO INTERP SHOULD BE
C     REMOVED AND CONTROL TRANSFERRED TO STATEMENT 500 INSTEAD OF 520.
C --------------------------------------------------------------------
      IF(NSTEP.GT.MAXSTP) THEN
         KGO=5
         KFLAG=-6
c   TOO MUCH WORK
      CALL Rprint('Number of steps exceeds maximum')

         IDID = -6
         GOTO 70
      END IF

      IF(ITOL.EQ.1) THEN
         D = 0.0D+0
         DO 40 I = 1,N
            AYI = DABS(Y(I,1))
            YMAX(I)=DMAX1(YMAX(I),AYI)
 40      CONTINUE
      ENDIF
      IF (IDID.EQ.3.OR.IDID.EQ.1) GO TO 70

      IF (DABS(T-TOUT).LE.DABS(15.0D+0*UROUND*TOUT)) THEN
C        EFFECTIVELY WE HAVE HIT TOUT
         IDID = KFLAG
         T0 = TOUT
         DO 50 I = 1,N
            Y0(I) = Y(I,1)
 50      CONTINUE
         HO = H
         RETURN

      END IF

      IF (IDID.EQ.2) THEN
C        CONTINUING INTEGRATION BUT MUST HIT TOUT EXACTLY
         IF (((T+H)-TOUT)*H.GT.0.0D+0) THEN
C           WE HAVE ALREADY OVERSHOT THE OUTPUT POINT OR WE WILL DO
C           SO ON THE NEXT STEP
            IF (((T-TOUT)*H.GE.0.0D+0) .OR. (DABS(T-TOUT).LE.
     +           100.0D+0*UROUND*HMAX)) THEN
C              HAVE OVERSHOT, SO INTERPOLATE
               CALL INTERP(N,JSTART,H,T,Y,TOUT,Y0)
               T0 = TOUT
               HO = H
               IDID = KFLAG
               RETURN

            ELSE
C              WILL PASS TOUT ON NEXT STEP WITH CURRENT STEPSIZE
C              SO REDUCE STEPSIZE TO HIT TOUT 'EXACTLY'
               H = (TOUT-T)* (1.0D+0-4.0D+0*UROUND)
               JSTART = -1
            END IF

         END IF

      ELSE IF ((T-TOUT)*H.GE.0.0D+0) THEN
C        HAVE OVERSHOT, SO INTERPOLATE
         CALL INTERP(N,JSTART,H,T,Y,TOUT,Y0)
         IDID = KFLAG
         HO = H
         T0 = TOUT
         RETURN

      END IF

      GO TO 20

C -------------------------------------------------------------------
C     ON AN ERROR RETURN FROM THE INTEGRATOR, AN IMMEDIATE RETURN OCCURS
C     IF KFLAG = -2, AND RECOVERY ATTEMPTS ARE MADE OTHERWISE.
C     H AND HMIN ARE REDUCED BY A FACTOR OF .1 UP TO 10 TIMES
C     BEFORE GIVING UP.
C --------------------------------------------------------------------
 60   CONTINUE
      IF (NHCUT.EQ.10) THEN
C        HAVE REDUCED H TEN TIMES
      CALL Rprint('Problem appears unsolvable with given input')
      CALL Rprint('         hmin reduced by a factor of 1.0e10')

         GO TO 70

      END IF

      NHCUT = NHCUT + 1
      HMIN = 0.1D+0*HMIN
      H = 0.1D+0*H
      JSTART = -1
      GO TO 20

 70   IF(DABS(T-TOUT).GT.1000.0D+0*UROUND) THEN
         DO 80 I = 1,N
            Y0(I) = Y(I,1)
 80      CONTINUE
         T0 = T

      ELSE
C        HAVE PASSED TOUT SO INTERPOLATE
         CALL INTERP(N,JSTART,H,T,Y,TOUT,Y0)
         T0 = TOUT
         IDID = KFLAG
      END IF
      HO = H
      IF(KFLAG.NE.0) IDID = KFLAG
      RETURN
C -------------------------- END OF SUBROUTINE OVDRIV -----------------


      END
C--------------------------------------------------------------------------
C
      SUBROUTINE INTERP(N,JSTART,H,T,Y,TOUT,Y0)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     .. SCALAR ARGUMENTS ..
      INTEGER JSTART,N
C     ..
C     .. ARRAY ARGUMENTS ..
      DIMENSION  Y(N,12),Y0(N)
C     ..
C     .. LOCAL SCALARS ..
      INTEGER I,J,L
C     ..
C     .. INTRINSIC FUNCTIONS ..
C     ..
      DO 10 I = 1,N
         Y0(I) = Y(I,1)
 10   CONTINUE
      L = JSTART + 2
      S = (TOUT-T)/H
      S1 = 1.0D+0
      DO 30 J = 2,L
         S1 = S1* (S+DBLE(FLOAT(J-2)))/DBLE(FLOAT(J-1))
         DO 20 I = 1,N
            Y0(I) = Y0(I) + S1*Y(I,J)
 20      CONTINUE
 30   CONTINUE
      RETURN
C -------------- END OF SUBROUTINE INTERP ---------------------------
      END
C
      SUBROUTINE COSET(NQ,EL,ELST,TQ,NCOSET,MAXORD)


      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C --------------------------------------------------------------------
C     COSET IS CALLED BY THE INTEGRATOR AND SETS THE COEFFICIENTS USED
C     BY THE CONVENTIONAL BACKWARD DIFFERENTIATION SCHEME AND THE
C     MODIFIED EXTENDED BACKWARD DIFFERENTIATION SCHEME.  THE VECTOR
C     EL OF LENGTH NQ+1 DETERMINES THE BASIC BDF METHOD WHILE THE VECTOR
C     ELST OF LENGTH NQ+2 DETERMINES THE MEBDF.  THE VECTOR TQ OF
C     LENGTH 4 IS INVOLVED IN ADJUSTING THE STEPSIZE IN RELATION TO THE
C     TRUNCATION ERROR.  ITS VALUES ARE GIVEN BY THE PERTST ARRAY.  THE
C     VECTORS EL AND TQ BOTH DEPEND ON METH AND NQ.  THE
C     COEFFICIENTS IN PERTST NEED TO BE GIVEN TO ONLY ABOUT ONE PERCENT
C     ACCURACY.  THE ORDER IN WHICH THE GROUPS APPEAR BELOW IS:
C     COEFFICIENTS FOR ORDER NQ-1, COEFFICIENTS FOR ORDER NQ,
C     COEFFICIENTS FOR ORDER NQ+1.
C -------------------------------------------------------------------
C     .. SCALAR ARGUMENTS ..
      INTEGER NQ
C     ..
C     .. ARRAY ARGUMENTS ..
      DIMENSION  EL(10),ELST(10),TQ(5)
C     ..
C     .. LOCAL SCALARS ..
      INTEGER K
C     ..
C     .. LOCAL ARRAYS ..
      DIMENSION  PERTST(8,3)
C     ..
C     .. INTRINSIC FUNCTIONS ..
C     ..
C     .. COMMON BLOCKS ..
      INTEGER NCOSET,MAXORD
C     ..
C     .. DATA STATEMENTS ..
      DATA  PERTST(1,1)/1./,PERTST(2,1)/2./,PERTST(3,1)/4.5/,
     +      PERTST(4,1)/7.333/,PERTST(5,1)/10.42/,PERTST(6,1)/13.7/,
     +      PERTST(7,1)/17.15/,PERTST(8,1)/20.74/
      DATA  PERTST(1,2)/2./,PERTST(2,2)/4.5/,PERTST(3,2)/7.333/,
     +      PERTST(4,2)/10.42/,PERTST(5,2)/13.7/,PERTST(6,2)/17.15/,
     +      PERTST(7,2)/20.74/,PERTST(8,2)/24.46/
      DATA  PERTST(1,3)/4.5/,PERTST(2,3)/7.333/,PERTST(3,3)/10.42/,
     +      PERTST(4,3)/13.7/,PERTST(5,3)/17.15/,PERTST(6,3)/20.74/,
     +      PERTST(7,3)/24.46/,PERTST(8,3)/1./
C     ..
C -------------------------------------------------------------------
C     THE FOLLOWING COEFFICIENTS SHOULD BE DEFINED TO MACHINE ACCURACY.
C     THEIR DERIVATION IS GIVEN IN REFERENCE 1.
C -------------------------------------------------------------------
      IF (NQ.GT.MAXORD) MAXORD = NQ
      NCOSET = NCOSET + 1
      IF (NQ .EQ. 1) THEN
        GOTO 10
      ELSE IF (NQ .EQ. 2) THEN
        GOTO 20
      ELSE IF (NQ .EQ. 3) THEN
        GOTO 30
      ELSE IF (NQ .EQ. 4) THEN
        GOTO 40
      ELSE IF (NQ .EQ. 5) THEN
        GOTO 50
      ELSE IF (NQ .EQ. 6) THEN
        GOTO 60
      ELSE IF (NQ .EQ. 7) THEN
        GOTO 70
      ENDIF        
C      GO TO (10,20,30,40,50,60,70) NQ

   10 EL(1) = 1.0D+0
      ELST(1) = 1.5D+0
      ELST(3) = -0.5D+0
      GO TO 80

   20 EL(1) = 6.6666666666667D-01
      EL(3) = 3.3333333333333D-01
      ELST(1) = 9.5652173913043D-01
      ELST(3) = 2.1739130434782D-01
      ELST(4) = -1.7391304347826D-01
      GO TO 80

   30 EL(1) = 5.4545454545455D-01
      EL(3) = 4.5454545454545D-01
      EL(4) = 1.8181818181818D-01
      ELST(1) = 7.6142131979695D-01
      ELST(3) = 3.2994923857868D-01
      ELST(4) = 8.6294416243654D-02
      ELST(5) = -9.1370558375634D-02
      GO TO 80

   40 EL(1) = 0.48D+0
      EL(3) = 0.52D+0
      EL(4) = 0.28D+0
      EL(5) = 0.12D+0
      ELST(1) = 6.5733706517393D-01
      ELST(3) = 4.0023990403838D-01
      ELST(4) = 1.5793682526989D-01
      ELST(5) = 4.4382247101159D-02
      ELST(6) = -5.7576969212315D-02
      GO TO 80

   50 EL(1) = 4.3795620437956D-01
      EL(3) = 5.62043795620436D-01
      EL(4) = 3.43065693430656D-01
      EL(5) = 1.97080291970802D-01
      EL(6) = 8.75912408759123D-02
      ELST(1) = 5.9119243917152D-01
      ELST(3) = 4.4902473356122D-01
      ELST(4) = 2.1375427307460D-01
      ELST(5) = 9.0421610027481503D-02
      ELST(6) = 2.6409276761177D-02
      ELST(7) = -4.0217172732757D-02
      GO TO 80

   60 EL(1) = 4.08163265306120D-01
      EL(3) = 5.91836734693874D-01
      EL(4) = 3.87755102040813D-01
      EL(5) = 2.51700680272107D-01
      EL(6) = 1.49659863945577D-01
      EL(7) = 6.80272108843534D-02
      ELST(1) = 5.4475876041119D-01
      ELST(3) = 4.8525549636077D-01
      ELST(4) = 2.5789750131312D-01
      ELST(5) = 1.3133738525800D-01
      ELST(6) = 5.7677396763462D-02
      ELST(7) = 1.7258197643881D-02
      ELST(8) = -3.0014256771967D-02
      GO TO 80

   70 EL(1) = 3.85674931129476D-01
      EL(3) = 6.14325068870521D-01
      EL(4) = 4.21487603305783D-01
      EL(5) = 2.9292929292929D-01
      EL(6) = 1.96510560146923D-01
      EL(7) = 1.19375573921028D-01
      EL(8) = 5.50964187327820D-02
      ELST(1) = 5.0999746293734D-01
      ELST(3) = 5.1345839935281D-01
      ELST(4) = 2.9364346131937D-01
      ELST(5) = 1.6664672120553D-01
      ELST(6) = 8.8013735242353D-02
      ELST(7) = 3.9571794884069D-02
      ELST(8) = 1.2039080338722D-02
      ELST(9) = -2.3455862290154D-02
   80 DO 90 K = 1,3
        TQ(K) = PERTST(NQ,K)
   90 CONTINUE
      TQ(4) = 0.5D+0*TQ(2)/DBLE(FLOAT(NQ))
      IF(NQ.NE.1) TQ(5)=PERTST(NQ-1,1)
      RETURN
C --------------------- END OF SUBROUTINE COSET ---------------------
      END

      SUBROUTINE PSET(Y,YPRIME,N,NPD,H,T,UROUND,EPSJAC,CON,MITER,MBND,
     +     IER,pderv,resid,NRENEW,YMAX,SAVE1,SAVE2,
     +     SAVE3,PW,PWCOPY,WRKSPC,IPIV,ITOL,RTOL,ATOL,NPSET,NJE,NRE,
     +     NDEC,IPAR,RPAR,IERR)


      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C -------------------------------------------------------------------
C     PSET IS CALLED BY STIFF TO COMPUTE AND PROCESS THE MATRIX
C     PD=DG/DY + (1/CON)DG/DY'. THIS MATRIX IS THEN SUBJECTED TO LU
C     DECOMPOSITION IN PREPARATION FOR LATER SOLUTION OF LINEAR SYSTEMS
C     OF ALGEBRAIC EQUATIONS WITH LU AS THE COEFFICIENT MATRIX.  THE
C     MATRIX PD IS FOUND BY THE USER-SUPPLIED ROUTINE pderv IF MITER=1
C     OR 3 OR BY FINITE DIFFERENCING IF MITER = 2 OR 4.
C     IN ADDITION TO VARIABLES DESCRIBED PREVIOUSLY, COMMUNICATION WITH
C     PSET USES THE FOLLOWING ..
C     EPSJAC = DSQRT(UROUND), USED IN NUMERICAL JACOBIAN INCREMENTS.
C *******************************************************************
C     THE ARGUMENT NRENEW IS USED TO SIGNAL WHETHER OR NOT
C     WE REQUIRE A NEW JACOBIAN TO BE CALCULATED.
C
C        IF NRENEW > 0 THEN WE REQUIRE A NEW J TO BE COMPUTED
C                  = 0 THEN USE A COPY OF THE LAST J COMPUTED
C *******************************************************************

C     .. SCALAR ARGUMENTS ..
      INTEGER IER,MITER,N,NRENEW,NPD
C     ..
C     .. ARRAY ARGUMENTS ..
      DIMENSION  PW(NPD*N),PWCOPY(NPD*N),SAVE2(N),WRKSPC(N),Y(N,12),
     + YMAX(N),SAVE1(N),SAVE3(N),RTOL(*),ATOL(*),IPAR(*),RPAR(*),
     + YPRIME(N)
      INTEGER IPIV(N), MBND(4)
C     ..
C     .. LOCAL SCALARS ..
      INTEGER I,J,J1,JJKK
C     ..
C     .. EXTERNAL SUBROUTINES ..
      EXTERNAL DEC_mebdfi,pderv,DGBFA,resid
C     ..
C     .. INTRINSIC FUNCTIONS ..
      INTRINSIC DABS,DMAX1,DSQRT
C     ..
C     .. COMMON BLOCKS ..
      INTEGER NDEC,NRE,NJE,NPSET

      NPSET = NPSET + 1
      ML = MBND(1)
      MU = MBND(2)
      IF (NRENEW.EQ.0) THEN
         IF (MITER.LT.3) THEN
            DO 10 I = 1,N*N
               PW(I) = PWCOPY(I)
 10         CONTINUE
         ELSE
            DO 15 I=0,N-1
               DO 18 J=MBND(1)+1,MBND(4)
                  PW(I*MBND(4)+J) = PWCOPY(I*MBND(4)+J)
 18            CONTINUE
 15         CONTINUE
         ENDIF
         GO TO 70
      ENDIF
C

C
C     PWCOPY IS INIZIALISED
C
      DO 19 I=1,NPD*N
        PWCOPY(I) = 0.0d0
 19   CONTINUE
C
      IF (MITER.EQ.2.OR.MITER.EQ.4) GO TO 30
C
      NJE = NJE + 1
      IF ( MITER.NE.3 ) THEN
C KS:         CALL PDERV(T,Y,PWCOPY,N,YPRIME,N,CON,IPAR,RPAR,IERR)
             CALL pderv(T, Y, YPRIME, PWCOPY, CON, RPAR, IPAR)
         DO 20 I=1,N*N
            PW(I)=PWCOPY(I)
 20      CONTINUE

      ELSE
C KS:         CALL PDERV(T,Y,PWCOPY,N,YPRIME,MBND(4),CON,IPAR,RPAR,IERR)
             CALL pderv (T, Y, YPRIME, PWCOPY, CON, RPAR, IPAR)
         DO 25 I=0,N-1
            ITMB=I*MBND(4)
            DO 28 J=MBND(1)+1,MBND(4)
               PW(ITMB+J) = PWCOPY(ITMB+J-MBND(1))
 28         CONTINUE
 25      CONTINUE
      ENDIF
      GOTO 70
C
 30   CONTINUE
      NJE = NJE + 1
      DO 35 I=1,N
         IF(ITOL.EQ.2) THEN
            YMAX(I)=Y(I,1)*RTOL(1)+ATOL(1)
         ELSE IF (ITOL.EQ.3) THEN
            YMAX(I)=Y(I,1)*RTOL(1)+ATOL(I)
         ELSE IF(ITOL.EQ.4) THEN
            YMAX(I)=Y(I,1)*RTOL(I)+ATOL(1)
         ELSE IF(ITOL.EQ.5) THEN
            YMAX(I)=Y(I,1)*RTOL(I)+ATOL(I)
         END IF
 35   CONTINUE
      D = 0.0D+0
      DO 40 I = 1,N
         D = D + (YPRIME(I)*(YMAX(I)))**2
 40   CONTINUE
      IF(ITOL.EQ.1) D=D*RTOL(1)**2
      D = DSQRT(D)/DBLE(N)
      R0 = DABS(H)*D*DBLE(N)*(1.0D+03)*UROUND
      IF(R0.EQ.0.0D+0) R0 = epsjac
      J1 = 0

      IF(MITER.EQ.4) GOTO 51
C KS: CALL RESID(N,T,Y,SAVE2,YPRIME,IPAR,RPAR,IERR)
      CALL resid(T,Y,Yprime, CON, SAVE2, ierr, rpar, ipar)
C FM: added check if ierr is OK
      if (ierr .ne. 0) return
      NRE=NRE+1
      DO 60 J = 1,N
         YJ = Y(J,1)
         YP=YPRIME(J)
         IF(ITOL.EQ.1) THEN
            R=DMAX1(EPSJAC*DABS(YJ),R0/(YMAX(J)*RTOL(1)))
         ELSE
            R=DMAX1(EPSJAC*DABS(YJ),R0/YMAX(J))
         ENDIF

         R = EPSJAC
         Y(J,1) = Y(J,1) + R
         YPRIME(J)=YPRIME(J)+R/CON

C KS:    CALL RESID(N,T,Y,WRKSPC,YPRIME,IPAR,RPAR,IERR)
         CALL resid (T, Y, YPRIME, CON, WRKSPC, IERR, RPAR, IPAR)
C FM: added check if ierr is OK
         if (ierr .ne. 0) return
         DO 50 I = 1,N
            JJKK = I + J1
            TEMPRY = (WRKSPC(I)-SAVE2(I))
            PWCOPY(JJKK) = TEMPRY/R
            PW(JJKK) =  PWCOPY(JJKK)
 50      CONTINUE
         Y(J,1) = YJ
         YPRIME(J)=YP
         J1 = J1 + N
 60   CONTINUE
      NRE = NRE + N


      GOTO 70

C
 51    CONTINUE

C KS:  CALL RESID(N,T,Y,SAVE2,YPRIME,IPAR,RPAR,IERR)
      CALL resid (T, Y, YPRIME, CON, SAVE2, IERR, RPAR, IPAR)
C FM: added check if ierr is OK
      if (ierr .ne. 0) return
      NRE = NRE+1
      MBA = min0(MBND(3),N)
      DO 61 J=1,MBA
         DO 161 I=J,N,MBND(3)
            SAVE1(I) = Y(I,1)
            SAVE3(I)=YPRIME(I)
            YI=Y(I,1)
c            YP=YPRIME(I)
            IF(ITOL.EQ.1) THEN
               R=DMAX1(EPSJAC*DABS(YI),R0/(YMAX(I)*RTOL(1)))
            ELSE
               R=DMAX1(EPSJAC*DABS(YI),R0/YMAX(I))
            ENDIF
            Y(I,1)=Y(I,1)+R
            YPRIME(I)=YPRIME(I)+R/CON
 161     CONTINUE
C KS:         CALL RESID(N,T,Y,WRKSPC,YPRIME,IPAR,RPAR,IERR)
         CALL resid (T, Y, YPRIME, CON, WRKSPC, IERR, RPAR, IPAR)
C FM: added check if ierr is OK
         if (ierr .ne. 0) return
         DO 261 JJ=J,N,MBND(3)
            Y(JJ,1)=SAVE1(JJ)
            YPRIME(JJ)=SAVE3(JJ)
            YJJ=Y(JJ,1)
c            YJP=YPRIME(JJ)
            IF (ITOL.EQ.1) THEN
               R=DMAX1(EPSJAC*DABS(YJJ),R0/(YMAX(JJ)*RTOL(1)))
            ELSE
               R=DMAX1(EPSJAC*DABS(YJJ),R0/YMAX(JJ))
            ENDIF
c           D=CON/R
            I1 = MAX0(JJ-MU,1)
            I2 = MIN0(JJ+ML,N)
            II = JJ*(MBND(4)-1)-ML
            DO 540 I = I1,I2
               TEMPRY = WRKSPC(I) - SAVE2(I)
               PWCOPY(II+I)=TEMPRY/R
               PW(II+I) = PWCOPY(II+I)
 540        CONTINUE
 261     CONTINUE
 61   CONTINUE
      NRE=NRE+MIN(MBND(3),N)
C


 70   IF (MITER.GT.2) THEN
         CALL DGBFA(PW,MBND(4),N,ML,MU,IPIV,IER)
         NDEC = NDEC + 1
      ELSE
         CALL DEC_mebdfi(N,N,PW,IPIV,IER)
         NDEC = NDEC + 1
      ENDIF
      RETURN
C ---------------------- END OF SUBROUTINE PSET ---------------------
      END
C
      SUBROUTINE DEC_mebdfi(N,NDIM,A,IP,IER)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C -------------------------------------------------------------------
C     MATRIX TRIANGULISATION BY GAUSSIAN ELIMINATION
C     INPUT..
C     N = ORDER OF MATRIX.
C     NDIM = DECLARED DIMENSION OF ARRAY A.
C     A = MATRIX TO BE TRIANGULARISED.
C     OUTPUT..
C     A(I,J),  I.LE.J = UPPER TRIANGULAR FACTOR, U.
C     A(I,J),  I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I-L.
C     IP(K), K.LT.N = INDEX OF KTH PIVOT ROW.
C     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR 0.
C     IER = 0 IF MATRIX IS NON-SINGULAR, OR K IF FOUND TO BE SINGULAR
C                  AT STAGE K.
C     USE SOL_mebdfi TO OBTAIN SOLUTION OF LINEAR SYSTEM.
C     DETERM(A) = IP(N)*A(1,1)*A(2,2)* . . . *A(N,N).
C     IF IP(N) = 0, A IS SINGULAR, SOL_mebdfi WILL DIVIDE BY ZERO.
C
C     REFERENCE.
C     C.B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, C.A.C.M
C     15 (1972), P.274.
C     ------------------------------------------------------------------
C     .. SCALAR ARGUMENTS ..
      INTEGER IER,N,NDIM
C     ..
C     .. ARRAY ARGUMENTS ..
      DIMENSION  A(NDIM,N)
      INTEGER IP(N)
C     ..
C     .. LOCAL SCALARS ..
      INTEGER I,J,K,KP1,M,NM1
C     ..
C     .. INTRINSIC FUNCTIONS ..
      INTRINSIC DABS
C     ..
C     .. COMMON BLOCKS ..
C     ..
      IER = 0
      IP(N) = 1
      IF (N.EQ.1) GO TO 70
      NM1 = N - 1
      DO 60 K = 1,NM1
         KP1 = K + 1
         M = K
         DO 10 I = KP1,N
            IF (DABS(A(I,K)).GT.DABS(A(M,K))) M = I
 10      CONTINUE
         IP(K) = M
         T = A(M,K)
         IF (M.EQ.K) GO TO 20
         IP(N) = -IP(N)
         A(M,K) = A(K,K)
         A(K,K) = T
 20      IF (T.EQ.0.0D+0) GO TO 80
         T = 1.0D+0/T
         DO 30 I = KP1,N
            A(I,K) = -A(I,K)*T
 30      CONTINUE
         DO 50 J = KP1,N
            T = A(M,J)
            A(M,J) = A(K,J)
            A(K,J) = T
            IF (T.EQ.0.0D+0) GO TO 50
            DO 40 I = KP1,N
               A(I,J) = A(I,J) + A(I,K)*T
 40         CONTINUE
 50      CONTINUE
 60   CONTINUE
 70   K = N
      IF (A(N,N).EQ.0.0D+0) GO TO 80
      RETURN

 80   IER = K
      IP(N) = 0
      RETURN
C--------------------- END OF SUBROUTINE DEC_mebdfi ----------------------
      END
C
      SUBROUTINE SOL_mebdfi(N,NDIM,A,B,IP)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     .. SCALAR ARGUMENTS ..
      INTEGER N,NDIM
C     ..
C     .. ARRAY ARGUMENTS ..
      DIMENSION  A(NDIM,N),B(N)
      INTEGER IP(N)
C     ..
C     .. LOCAL SCALARS ..
      INTEGER I,K,KB,KM1,KP1,M,NM1
C     ..
C     .. COMMON BLOCKS ..
C     ..

C     ------------------------------------------------------------------
C     SOLUTION OF LINEAR SYSTEM, A*X = B.
C     INPUT ..
C     N = ORDER OF MATRIX.
C     NDIM = DECLARED DIMENSION OF MATRIX A.
C     A = TRIANGULARISED MATRIX OBTAINED FROM DEC_mebdfi.
C     B = RIGHT HAND SIDE VECTOR.
C     IP = PIVOT VECTOR OBTAINED FROM DEC_mebdfi.
C     DO NOT USE IF DEC_mebdfi HAS SET IER .NE. 0
C     OUTPUT..
C     B = SOLUTION VECTOR, X.
C     ------------------------------------------------------------------
      IF (N.EQ.1) GO TO 50
      NM1 = N - 1
      DO 20 K = 1,NM1
         KP1 = K + 1
         M = IP(K)
         T = B(M)
         B(M) = B(K)
         B(K) = T
         DO 10 I = KP1,N
            B(I) = B(I) + A(I,K)*T
 10      CONTINUE
 20   CONTINUE
      DO 40 KB = 1,NM1
         KM1 = N - KB
         K = KM1 + 1
         B(K) = B(K)/A(K,K)
         T = -B(K)
         DO 30 I = 1,KM1
            B(I) = B(I) + A(I,K)*T
 30      CONTINUE
 40   CONTINUE
 50   B(1) = B(1)/A(1,1)
      RETURN
C------------------------- END OF SUBROUTINE SOL_mebdfi ------------------
      END
C
      subroutine dgbfa(abd,lda,n,ml,mu,ipvt,info)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      integer lda,n,ml,mu,ipvt(*),info
      double precision abd(lda,*)
c
c     dgbfa factors a double precision band matrix by elimination.
c
c     dgbfa is usually called by dgbco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c
c     on entry
c
c        abd     double precision(lda, n)
c                contains the matrix in band storage.  the columns
c                of the matrix are stored in the columns of  abd  and
c                the diagonals of the matrix are stored in rows
c                ml+1 through 2*ml+mu+1 of  abd .
c                see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array  abd .
c                lda must be .ge. 2*ml + mu + 1 .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c                0 .le. ml .lt. n .
c
c        mu      integer
c                number of diagonals above the main diagonal.
c                0 .le. mu .lt. n .
c                more efficient if  ml .le. mu .
c     on return
c
c        abd     an upper triangular matrix in band storage and
c                the multipliers which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgbsl will divide by zero if
c                     called.  use  rcond  in dgbco for a reliable
c                     indication of singularity.
c
c     band storage
c
c           if  a  is a band matrix, the following program segment
c           will set up the input.
c
c                   ml = (band width below the diagonal)
c                   mu = (band width above the diagonal)
c                   m = ml + mu + 1
c                   do 20 j = 1, n
c                      i1 = max0(1, j-mu)
c                      i2 = min0(n, j+ml)
c                      do 10 i = i1, i2
c                         k = i - j + m
c                         abd(k,j) = a(i,j)
c                10    continue
c                20 continue
c
c           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
c           in addition, the first  ml  rows in  abd  are used for
c           elements generated during the triangularization.
c           the total number of rows needed in  abd  is  2*ml+mu+1 .
c           the  ml+mu by ml+mu  upper left triangle and the
c           ml by ml  lower right triangle are not referenced.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c     fortran max0,min0
c
c     internal variables
c
      double precision t
      integer i,idamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1
c
c
      m = ml + mu + 1
      info = 0
c
c     zero initial fill-in columns
c
      j0 = mu + 2
      j1 = min0(n,m) - 1
      if (j1 .lt. j0) go to 30
      do 20 jz = j0, j1
         i0 = m + 1 - jz
         do 10 i = i0, ml
            abd(i,jz) = 0.0d0
 10      continue
 20   continue
 30   continue
      jz = j1
      ju = 0
c
c     gaussian elimination with partial pivoting
c
      nm1 = n - 1
      if (nm1 .lt. 1) go to 130
      do 120 k = 1, nm1
         kp1 = k + 1
c
c        zero next fill-in column
c
         jz = jz + 1
         if (jz .gt. n) go to 50
         if (ml .lt. 1) go to 50
         do 40 i = 1, ml
            abd(i,jz) = 0.0d0
 40      continue
 50      continue
c
c        find l = pivot index
c
         lm = min0(ml,n-k)
         l = idamax(lm+1,abd(m,k),1) + m - 1
         ipvt(k) = l + k - m
c
c        zero pivot implies this column already triangularized
c
         if (abd(l,k) .eq. 0.0d0) go to 100
c
c           interchange if necessary
c
         if (l .eq. m) go to 60
         t = abd(l,k)
         abd(l,k) = abd(m,k)
         abd(m,k) = t
 60      continue
c
c           compute multipliers
c
         t = -1.0d0/abd(m,k)
         call dscal(lm,t,abd(m+1,k),1)
c
c           row elimination with column indexing
c
         ju = min0(max0(ju,mu+ipvt(k)),n)
         mm = m
         if (ju .lt. kp1) go to 90
         do 80 j = kp1, ju
            l = l - 1
            mm = mm - 1
            t = abd(l,j)
            if (l .eq. mm) go to 70
            abd(l,j) = abd(mm,j)
            abd(mm,j) = t
 70         continue
            call daxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)
 80      continue
 90      continue
         go to 110
 100     continue
         info = k
 110     continue
 120  continue
 130  continue
      ipvt(n) = n
      if (abd(m,n) .eq. 0.0d0) info = n
      return
      end
C--------------------------------------------------------------------------
      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
         dy(iy) = dy(iy) + da*dx(ix)
         ix = ix + incx
         iy = iy + incy
 10   continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
 20   m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
         dy(i) = dy(i) + da*dx(i)
 30   continue
      if( n .lt. 4 ) return
 40   mp1 = m + 1
      do 50 i = mp1,n,4
         dy(i) = dy(i) + da*dx(i)
         dy(i + 1) = dy(i + 1) + da*dx(i + 1)
         dy(i + 2) = dy(i + 2) + da*dx(i + 2)
         dy(i + 3) = dy(i + 3) + da*dx(i + 3)
 50   continue
      return
      end
C---------------------------------------------------------------------------
      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified to correct problem with negative increment, 8/21/90.
c
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      double precision da,dx(*)
      integer i,incx,ix,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      do 10 i = 1,n
         dx(ix) = da*dx(ix)
         ix = ix + incx
 10   continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
 20   m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
         dx(i) = da*dx(i)
 30   continue
      if( n .lt. 5 ) return
 40   mp1 = m + 1
      do 50 i = mp1,n,5
         dx(i) = da*dx(i)
         dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
 50   continue
      return
      end
C--------------------------------------------------------------------------
      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified to correct problem with negative increment, 8/21/90.
c
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      double precision dx(*),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n .lt. 1 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      dmax = dabs(dx(ix))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
 10   continue
      return
c
c        code for increment equal to 1
c
 20   dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
 30   continue
      return
      end
C--------------------------------------------------------------------------
      subroutine dgbsl(abd,lda,n,ml,mu,ipvt,b,job)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      integer lda,n,ml,mu,ipvt(*),job
      double precision abd(lda,*),b(*)
c
c     dgbsl solves the double precision band system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgbco or dgbfa.
c
c     on entry
c
c        abd     double precision(lda, n)
c                the output from dgbco or dgbfa.
c
c        lda     integer
c                the leading dimension of the array  abd .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c
c        mu      integer
c                number of diagonals above the main diagonal.
c
c        ipvt    integer(n)
c                the pivot vector from dgbco or dgbfa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b , where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgbco has set rcond .gt. 0.0
c        or dgbfa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c     fortran min0
c
c     internal variables
c
      double precision ddot,t
      integer k,kb,l,la,lb,lm,m,nm1
c
      m = mu + ml + 1
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve l*y = b
c
      if (ml .eq. 0) go to 30
      if (nm1 .lt. 1) go to 30
      do 20 k = 1, nm1
         lm = min0(ml,n-k)
         l = ipvt(k)
         t = b(l)
         if (l .eq. k) go to 10
         b(l) = b(k)
         b(k) = t
 10      continue
         call daxpy(lm,t,abd(m+1,k),1,b(k+1),1)
 20   continue
 30   continue
c
c        now solve  u*x = y
c
      do 40 kb = 1, n
         k = n + 1 - kb
         b(k) = b(k)/abd(m,k)
         lm = min0(k,m) - 1
         la = m - lm
         lb = k - lm
         t = -b(k)
         call daxpy(lm,t,abd(la,k),1,b(lb),1)
 40   continue
      go to 100
 50   continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
      do 60 k = 1, n
         lm = min0(k,m) - 1
         la = m - lm
         lb = k - lm
         t = ddot(lm,abd(la,k),1,b(lb),1)
         b(k) = (b(k) - t)/abd(m,k)
 60   continue
c
c        now solve trans(l)*x = y
c
      if (ml .eq. 0) go to 90
      if (nm1 .lt. 1) go to 90
      do 80 kb = 1, nm1
         k = n - kb
         lm = min0(ml,n-k)
         b(k) = b(k) + ddot(lm,abd(m+1,k),1,b(k+1),1)
         l = ipvt(k)
         if (l .eq. k) go to 70
         t = b(l)
         b(l) = b(k)
         b(k) = t
 70      continue
 80   continue
 90   continue
 100  continue
      return
      end
C---------------------------------------------------------------------------
      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
         dtemp = dtemp + dx(ix)*dy(iy)
         ix = ix + incx
         iy = iy + incy
 10   continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
 20   m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
         dtemp = dtemp + dx(i)*dy(i)
 30   continue
      if( n .lt. 5 ) go to 60
 40   mp1 = m + 1
      do 50 i = mp1,n,5
         dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *        dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) +
     *        dx(i + 4)*dy(i + 4)
 50   continue
 60   ddot = dtemp
      return
      end
C---------------------------------------------------------------------------

      SUBROUTINE ERRORS(N,TQ,EDN,E,EUP,BND,EDDN)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     ***************************************************
C
C     THIS ROUTINE CALCULATES ERRORS USED IN TESTS
C     IN STIFF .
C
C     ***************************************************
C     .. SCALAR ARGUMENTS ..
      INTEGER N
C     ..
C     .. ARRAY ARGUMENTS ..
      DIMENSION  TQ(5)
C     ..
C     .. LOCAL SCALARS ..
C     ..
C     .. INTRINSIC FUNCTIONS ..
C     ..
      SQHOL = DBLE(FLOAT(N))
      EDN = TQ(1)*TQ(1)*SQHOL
C
C     ** ERROR ASSOCIATED WITH  METHOD OF ORDER ONE LOWER.
C
      E = TQ(2)*TQ(2)*SQHOL
C
C     ** ERROR ASSOCIATED WITH PRESENT ORDER
C
      EUP = TQ(3)*TQ(3)*SQHOL
C
C     ** ERROR ASSOCIATED WITH HIGHER ORDER METHOD
C
      BND = TQ(4)*TQ(4)*SQHOL*0.5D+0

      EDDN=TQ(5)*TQ(5)*SQHOL
C
C     ** ERROR ASSOCIATED WITH METHOD OF ORDER TWO LOWER.
      RETURN
      END
C--------------------------------------------------------------------------

      SUBROUTINE PRDICT(T,H,Y,L,N)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C **********************************************************************
C     PREDICTS A VALUE FOR Y AT (T+H) GIVEN THE HISTORY ARRAY AT T.
C **********************************************************************

C     .. SCALAR ARGUMENTS ..
      INTEGER L,N
C     ..
C     .. ARRAY ARGUMENTS ..
      DIMENSION  Y(N,12)
C     ..
C     .. LOCAL SCALARS ..
      INTEGER I,J2
C     ..
C     .. EXTERNAL SUBROUTINES ..
C     ..
C     .. COMMON BLOCKS ..
C     ..
      DO 10 I=1,N
         DO 20 J2 = 2,L
            Y(I,1) = Y(I,1) + Y(I,J2)
 20      CONTINUE
 10   CONTINUE
      T = T + H

      RETURN
      END
C------------------------------------------------------------------------

      SUBROUTINE ITRAT2(QQQ,Y,YPRIME,N,T,HBETA,ERRBND,ARH,CRATE,TCRATE
     +     ,M,WORKED,YMAX,ERROR,SAVE1,SAVE2,SCALE,PW,MF,MBND,NIND1,
     +     NIND2,NIND3,IPIV,LMB,ITOL,RTOL,ATOL,IPAR,RPAR,HUSED,NBSOL,
     +     NRE,NQUSED,resid,IERR)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     .. SCALAR ARGUMENTS ..
      INTEGER M,N,ITOL
      LOGICAL WORKED
C     ..
C     .. ARRAY ARGUMENTS ..
      DIMENSION ARH(N),ERROR(N),PW(*),SAVE1(N),SAVE2(N),Y(N,12),
     +     YMAX(N),RTOL(*),ATOL(*),SCALE(N),IPAR(*),RPAR(*),YPRIME(N)
      INTEGER IPIV(N),MBND(4)
C     ..
C     .. LOCAL SCALARS ..
      INTEGER I
C     ..
C     .. EXTERNAL SUBROUTINES ..
      EXTERNAL SOL_mebdfi,DGBSL,resid
C     ..
C     .. INTRINSIC FUNCTIONS ..
      INTRINSIC DMAX1,DMIN1
C     ..
C     .. COMMON BLOCKS ..
      INTEGER NBSOL,NRE

CKS: add this save statement?
      SAVE D1
C     ..
C     .. DATA STATEMENTS ..

      DATA  ZERO/0.0D+0/
C     ..



      DO 5 I=1,N
         AYI = DABS(Y(I,1))
         IF(ITOL.EQ.1) THEN
            SCALE(I) = YMAX(I)
         ELSE IF(ITOL.EQ.2) THEN
            SCALE(I) = RTOL(1)*AYI + ATOL(1)
         ELSE IF(ITOL.EQ.3) THEN
            SCALE(I) = RTOL(1)*AYI + ATOL(I)
         ELSE IF(ITOL.EQ.4) THEN
            SCALE(I) = RTOL(I)*AYI + ATOL(1)
         ELSE IF(ITOL.EQ.5) THEN
            SCALE(I) = RTOL(I)*AYI + ATOL(I)
         ENDIF
 5    CONTINUE

      IF(NIND2.NE.0) THEN
         DO 11 I = NIND1+1,NIND2+NIND1
            SCALE(I)=SCALE(I)/HUSED
 11      CONTINUE
      ENDIF
      IF(NIND3.NE.0) THEN
         DO 12 I = NIND1 +NIND2 + 1,NIND3+NIND2+NIND1
            SCALE(I)=SCALE(I)/(HUSED**2)
 12      CONTINUE
      ENDIF
C
      IF(LMB.EQ.1) GOTO 25
C
C KS:      call resid(n,t,y,save2,yprime,ipar,rpar,ierr)
      CALL resid (T, Y, YPRIME, hbeta, save2, IERR, RPAR, IPAR)
C FM: added check if ierr is OK
      if (ierr .ne. 0) return

      IF(MF.GE.23) THEN
         CALL DGBSL(PW,MBND(4),N,MBND(1),MBND(2),IPIV,SAVE2,0)
         NBSOL = NBSOL + 1
      ELSE
         CALL SOL_mebdfi(N,N,PW,SAVE2,IPIV)
         NBSOL = NBSOL + 1
      ENDIF
      D = ZERO
      DO 20 I = 1,N
         ERROR(I) = ERROR(I) - SAVE2(I)
         D = D + (SAVE2(I)/(SCALE(I)))**2
         SAVE1(I) = Y(I,1) + ERROR(I)
 20   CONTINUE
      IF(ITOL.EQ.1) D=D/(RTOL(1)**2)
      TCRATE = TCRATE + CRATE
      D1 = D
      M = 1
      do 1014 i=1,n
         yprime(i)=(save1(i)-arh(i))/qqq
 1014 continue
      NRE = NRE + 1
 25   CONTINUE
      WORKED = .TRUE.
 30   CONTINUE
C KS:      call resid(n,t,save1,save2,yprime,ipar,rpar,ierr)
      CALL resid (T, save1, YPRIME, hbeta, save2, IERR, RPAR, IPAR)
C FM: added check if ierr is OK
      if (ierr .ne. 0) return
      nre=nre+1
C
C     IF WE ARE HERE THEN PARTIALS ARE O.K.
C
      IF( MF.GE. 23) THEN
         CALL DGBSL(PW,MBND(4),N,MBND(1),MBND(2),IPIV,SAVE2,0)
         NBSOL=NBSOL + 1
      ELSE
         CALL SOL_mebdfi(N,N,PW,SAVE2,IPIV)
         NBSOL = NBSOL + 1
      ENDIF
C
C     WE NOW CALCULATE A WEIGHTED RMS TYPE NORM
C
      D = ZERO
      DO 50 I = 1,N
         ERROR(I) = ERROR(I) - SAVE2(I)
         D = D + (SAVE2(I)/(SCALE(I)))**2
         SAVE1(I) = Y(I,1) + ERROR(I)
 50   CONTINUE
      IF(ITOL.EQ.1) D=D/(RTOL(1)**2)
C -------------------------------------------------------------------
C     TEST FOR CONVERGENCE.  IF M.GT.0 , AN ESTIMATE OF THE CONVERGENCE
C     RATE CONSTANT IS STORED IN CRATE, AND THIS IS USED IN THE TEST.
C -------------------------------------------------------------------
      IF (M.NE.0) THEN
         IF (D1.NE.ZERO) CRATE = DMAX1(0.9D+0*CRATE,D/D1)
      END IF
      TCRATE = TCRATE + CRATE
      IF ((D*DMIN1(1.0D+0,2.0D+0*CRATE)).LT.ERRBND/DBLE(NQUSED))
     +     RETURN
      IF (M.NE.0) THEN
         IF (D.GT.D1) THEN
            WORKED = .FALSE.
            RETURN

         END IF

      END IF

      D1 = D
      IF (M.EQ.4) THEN
         WORKED = .FALSE.
         RETURN

      END IF

      M = M + 1
      do 40 i=1,n
         yprime(i)=(save1(i)-arh(i))/qqq
 40   continue
      GO TO 30

      END
C--------------------------------------------------------------------------

      SUBROUTINE STIFF(H,HMAX,HMIN,JSTART,KFLAG,MF,MBND,
     +     NIND1,NIND2,NIND3,T,TOUT,TEND,Y,YPRIME,N,NPD,
     +     YMAX,ERROR,SAVE1,SAVE2,SCALE,PW,PWCOPY,YHOLD,
     +     YNHOLD,ARH,IPIV,MAXDER,ITOL,RTOL,ATOL,RPAR,IPAR,
     +     pderv,resid,NQUSED,NSTEP,NFAIL,NRE,NJE,NDEC,NBSOL,NPSET,
     +     NCOSET,MAXORD,UROUND,EPSJAC,HUSED,IERR)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     ------------------------------------------------------------------
C     THE SUBROUTINE STIFF PERFORMS ONE STEP OF THE INTEGRATION OF AN
C     INITIAL VALUE PROBLEM FOR A SYSTEM OF
C     IMPLICIT DIFFERENTIAL ALGEBRAIC EQUATIONS.
C     COMMUNICATION WITH STIFF IS DONE WITH THE FOLLOWING VARIABLES..
C     Y      AN N BY LMAX+3 ARRAY CONTAINING THE DEPENDENT VARIABLES
C              AND THEIR BACKWARD DIFFERENCES.  MAXDER (=LMAX-1) IS THE
C              MAXIMUM ORDER AVAILABLE.  SEE SUBROUTINE COSET.
C              Y(I,J+1) CONTAINS THE JTH BACKWARD DIFFERENCE OF Y(I)
C     T      THE INDEPENDENT VARIABLE. T IS UPDATED ON EACH STEP TAKEN.
C     H      THE STEPSIZE TO BE ATTEMPTED ON THE NEXT STEP.
C              H IS ALTERED BY THE ERROR CONTROL ALGORITHM DURING
C              THE PROBLEM.  H CAN BE EITHER POSITIVE OR NEGATIVE BUT
C              ITS SIGN MUST REMAIN CONSTANT THROUGHOUT THE PROBLEM.
C     HMIN   THE MINIMUM AND MAXIMUM ABSOLUTE VALUE OF THE STEPSIZE
C     HMAX   TO BE USED FOR THE STEP.  THESE MAY BE CHANGED AT ANY
C              TIME BUT WILL NOT TAKE EFFECT UNTIL THE NEXT H CHANGE.
C     RTOL,ATOL  THE ERROR BOUNDS. SEE DESCRIPTION IN OVDRIV.
C     N      THE NUMBER OF FIRST ORDER DIFFERENTIAL EQUATIONS.
C     MF     THE METHOD FLAG.  MUST BE SET TO 21,22,23 OR 24 AT PRESENT
C     KFLAG  A COMPLETION FLAG WITH THE FOLLOWING MEANINGS..
C                  0  THE STEP WAS SUCCESSFUL
C                 -1  THE REQUESTED ERROR COULD NOT BE ACHIEVED
C                       WITH ABS(H) = HMIN.
C                 -2  THE REQUESTED ERROR IS SMALLER THAN CAN
C                       BE HANDLED FOR THIS PROBLEM.
C                 -3  CORRECTOR CONVERGENCE COULD NOT BE
C                       ACHIEVED FOR ABS(H)=HMIN.
C            ON A RETURN WITH KFLAG NEGATIVE, THE VALUES OF T AND
C            THE Y ARRAY ARE AS AT THE BEGINNING OF THE LAST
C            STEP ATTEMPTED, AND H IS THE LAST STEP SIZE ATTEMPTED.
C     JSTART  AN INTEGER USED ON INPUT AND OUTPUT.
C          ON INPUT IT HAS THE FOLLOWING VALUES AND MEANINGS..
C              0  PERFORM THE FIRST STEP.
C          .GT.0  TAKE A NEW STEP CONTINUING FROM THE LAST
C          .LT.0  TAKE THE NEXT STEP WITH A NEW VALUE OF H OR N.
C          ON EXIT, JSTART IS NQUSED, THE ORDER OF THE METHOD LAST USED.
C     YMAX     AN ARRAY OF N ELEMENTS WITH WHICH THE ESTIMATED LOCAL
C              ERRORS IN Y ARE COMPARED
C     ERROR    AN ARRAY OF N ELEMENTS.
C     SAVE1,2  TWO ARRAYS OF WORKING SPACE BOTH OF LENGTH N.
C     PW       A BLOCK OF LOCATIONS USED FOR PARTIAL DERIVATIVES
C     IPIV     AN INTEGER ARRAY OF LENGTH N USED FOR PIVOT INFORMATION.
C     JNEWIM   IS TO INDICATE IF PRESENT ITERATION MATRIX
C                WAS FORMED USING A NEW J OR OLD J.
C     JSNOLD   KEEPS TRACK OF NO. OF STEPS TAKEN WITH
C                PRESENT ITERATION MATRIX (BE IT FORMED BY
C                A NEW J OR NOT).
C     AVNEWJ   STORES VALUE FOR AVERAGE CRATE WHEN ITERATION
C                MATRIX WAS FORMED BY A NEW J.
C     AVOLDJ   STORES VALUE FOR AVERAGE CRATE WHEN ITERATION
C                MATRIX WAS FORMED BY AN OLD J.
C     NRENEW   FLAG THAT IS USED IN COMMUNICATION WITH SUBROUTINE PSET.
C                IF  NRENEW > 0  THEN FORM A NEW JACOBIAN BEFORE
C                                COMPUTING THE COEFFICIENT MATRIX FOR
C                                THE NEWTON-RAPHSON ITERATION
C                           = 0  FORM THE COEFFICIENT MATRIX USING A
C                                COPY OF AN OLD JACOBIAN
C     NEWPAR   FLAG USED IN THIS SUBROUTINE TO INDICATE IF A JACOBIAN
C              HAS BEEN EVALUATED FOR THE CURRENT STEP
C **********************************************************************

C     .. SCALAR ARGUMENTS ..
      INTEGER JSTART,KFLAG,MF,N,ITOL,NPD
C     ..
C     .. ARRAY ARGUMENTS ..
      DIMENSION HSTPSZ(2,14)
      DIMENSION  ARH(N),ERROR(N),PW(NPD*N),PWCOPY(NPD*N),SAVE1(N),
     + SAVE2(N),Y(N,12),YHOLD(N,12),YMAX(N),YNHOLD(N,2),RTOL(*),ATOL(*),
     +  SCALE(N),RPAR(*),IPAR(*),YPRIME(N)
      INTEGER IPIV(N), MBND(4)
C     ..
C     .. LOCAL SCALARS ..
      INTEGER I,IER,IITER,IITER2,IREDO,J,J1,J2,KFAIL,LL,LMAX,
     +        M3STEP,MAXDER,METH,MFOLD,MITER,NEWQ
      LOGICAL FINISH,OVRIDE,WORKED
C     ..
C     .. LOCAL ARRAYS ..
      DIMENSION  EL(10),ELST(10),TQ(5)
C     ..
C     .. EXTERNAL SUBROUTINES ..
      EXTERNAL COSET,CPYARY,ERRORS,HCHOSE,ITRAT2,
     +    PRDICT,PSET,RSCALE,SOL_mebdfi,DGBSL,pderv,resid
C     ..
C     .. INTRINSIC FUNCTIONS ..
      INTRINSIC DABS,DMAX1,DMIN1
C     ..
C     .. COMMON BLOCKS ..
      SAVE HSTPSZ
c      COMMON / STPSZE / HSTPSZ
      INTEGER IDOUB,ISAMP,IWEVAL,JCHANG,JSINUP,JSNOLD,L,M1,M2,MAXORD,
     +        MEQC1,MEQC2,MQ1TMP,MQ2TMP,NBSOL,NCOSET,NDEC,NEWPAR,
     +        NRE,NJE,NPSET,NQ,NQUSED,NRENEW,NSTEP
      LOGICAL CFAIL,JNEWIM,SAMPLE
C     ..
C     .. SAVE STATEMENT ..
      SAVE MFOLD,LMAX,MITER,IBND,KFAIL,JNEWIM,IEMB,
     + NEWPAR,NRENEW,JSINUP,JSNOLD,IDOUB,JCHANG,L,NQ,MEQC1,
     + MEQC2,MQ1TMP,MQ2TMP,ISAMP,IER,IREDO,IWEVAL
      SAVE EDN,EUP,BND,EDDN,EL,TQ,ELST,E,QI,QQQ,
     + RH,RMAX,TOLD,CRATE1,CRATE2,HOLD,
     + TCRAT1,TCRAT2,AVNEWJ,AVOLDJ,UPBND,
     + RC,VTOL,OLDLO,OVRIDE
      SAVE CFAIL,SAMPLE


C     ..
C     .. DATA STATEMENTS ..

      DATA  EL(2),ELST(2),OLDLO/3*1.0D+0/
      DATA  ZERO,ONE/0.0D+0,1.0D+0/
C     ..

 6000 TOLD = T
      KFLAG = 0
      CFAIL = .TRUE.
      DUP = ZERO    ! KARLINE ADDED THIS TO INITIALISE

      IF (JSTART.GT.0) GO TO 60
      IF (JSTART.NE.0) GO TO 30
C     ------------------------------------------------------------------
C     ON THE FIRST CALL, THE ORDER IS SET TO 1 AND THE INITIAL DERIVATIVE
C     IS CALCULATED OR GIVEN.  RMAX IS THE MAXIMUM RATIO BY WHICH H CAN BE
C     INCREASED IN A SINGLE STEP.  RMAX IS SET EQUAL TO 1.D4 INITIALLY
C     TO COMPENSATE FOR THE SMALL INITIAL H, BUT THEN IS NORMALLY = 10.
C     IF A FAILURE OCCURS (IN CORRECTOR CONVERGENCE OR ERROR TEST),
C     RMAX IS SET AT 2. FOR THE NEXT INCREASE.
C     ------------------------------------------------------------------
      DO 10 I = 1,N
         Y(I,2) = H*YPRIME(I)
 10   CONTINUE
      METH = 2
      MITER = MF - 10*METH
      IBND=5
      UPBND=0.2D+0
      IF(MF.GT.22) UPBND=0.1D+0
      NQ = 1
      NQUSED = NQ
      L = 2
      IDOUB = 3
      KFAIL = 0
      RMAX = 10000.0D+0
      IER=0
      RC = ZERO
      CRATE1 = 0.1D+0
      CRATE2 = 0.1D+0
      JSNOLD = 0
      JNEWIM = .TRUE.
      TCRAT1 = ZERO
      TCRAT2 = ZERO
      VTOL=DMAX1(RTOL(1),ATOL(1))/10.0D+0
      DO 15 I=1,12
         HSTPSZ(1,I)=1.0D+0
         HSTPSZ(2,I)=VTOL
 15   CONTINUE
      HOLD = H
      MFOLD = MF
      NSTEP = 0
      NRE = 1
      NJE = 0
      NDEC = 0
      NPSET = 0
      NCOSET = 0
      MAXORD = 1
      NFAIL = 0
      CFAIL = .TRUE.
      AVNEWJ = ZERO
      AVOLDJ = ZERO
c      AVNEW2 = ZERO
c      AVOLD2 = ZERO
      SAMPLE = .FALSE.
      ISAMP = 0
      IEMB=0
C     **************************************************
C     CFAIL=.TRUE. ENSURES THAT WE CALCULATE A NEW
C     J ON THE FIRST CALL.
C     **************************************************
      MEQC1 = 0
      MEQC2 = 0
      MQ1TMP = 0
      MQ2TMP = 0
      NBSOL = 0
      HUSED = H
C     -----------------------------------------------------------------
C     IF THE CALLER HAS CHANGED N , THE CONSTANTS E, EDN, EUP
C     AND BND MUST BE RESET.  E IS A COMPARISON FOR ERRORS AT THE
C     CURRENT ORDER NQ.  EUP IS TO TEST FOR INCREASING THE ORDER,
C     EDN FOR DECREASING THE ORDER.  BND IS USED TO TEST FOR CONVERGENCE
C     OF THE CORRECTOR ITERATES.   IF THE CALLER HAS CHANGED H, Y MUST
C     BE RE-SCALED.  IF H IS CHANGED, IDOUB IS SET TO L+1 TO PREVENT
C     FURTHER CHANGES IN H FOR THAT MANY STEPS.
C     -----------------------------------------------------------------
      CALL COSET(NQ,EL,ELST,TQ,NCOSET,MAXORD)
      LMAX = MAXDER + 1
      RC = RC*EL(1)/OLDLO
      OLDLO = EL(1)
      IWEVAL = MITER
      NRENEW = 1
      NEWPAR = 0
C     *****************************************************
C     NRENEW AND NEWPAR ARE TO INSTRUCT ROUTINE THAT
C     WE WISH A NEW J TO BE CALCULATED FOR THIS STEP.
C     *****************************************************
      CALL ERRORS(N,TQ,EDN,E,EUP,BND,EDDN)
      DO 20 I = 1,N
         ARH(I) = EL(2)*Y(I,1)
 20   CONTINUE
      CALL CPYARY(N*L,Y,YHOLD)
      QI = H*EL(1)
      QQ = ONE/QI
      CALL PRDICT(T,H,Y,L,N)
      IF(IERR.NE.0) THEN
         H=H/2
         IERR = 0
         GOTO 6000
      ENDIF
      DO 25 I=1,N
         YPRIME(I)=(Y(I,1)-ARH(I))/h
 25   CONTINUE
      GO TO 110
C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     DIFFERENT PARAMETERS ON THIS CALL        <
C     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 30   CALL CPYARY(N*L,YHOLD,Y)
      IF (MF.NE.MFOLD) THEN
         METH = MF/10
         MITER = MF - 10*METH
         MFOLD = MF
         IWEVAL = MITER
      END IF

      IF(NSTEP.GT.0) GOTO 35
      NJE = 0
      NRE = 1
      CFAIL = .TRUE.
      NEWPAR = 0
      MQ1TMP = 0
      MQ2TMP = 0
      MEQC1 = 0
      MEQC2 = 0
      TCRAT1 = 0.0D+0
      TCRAT2 = 0.0D+0
      CRATE1 = 1.0D-1
      CRATE2 = 1.0D-1
      NSTEP = 0
      NBSOL = 0
      NPSET = 0
      NCOSET = 0
      NDEC = 0
 35   CONTINUE
      IF (H.NE.HOLD) THEN
         RH = H/HOLD
         H = HOLD
         IREDO = 3
         GO TO 50

      ELSE
         GO TO 60

      END IF

C     *********************************************
C     RE-SCALE Y AFTER A CHANGE OF STEPSIZE   *
C     *********************************************
 40   RH = DMAX1(RH,HMIN/DABS(H))
 50   RH = DMIN1(RH,HMAX/DABS(H),RMAX)
      CALL RSCALE(N,L,RH,Y)
      RMAX = 10.0D+0
      JCHANG = 1
      H = H*RH
      RC = RC*RH
      IF (JSNOLD.GT.IBND) THEN
         CFAIL = .TRUE.
         NEWPAR = 0
         RC = ZERO
C **********************************************************************
C        CFAIL=TRUE AND NEWPAR=0 SHOULD FORCE A NEW J TO BE EVALUATED
C        AFTER 7 STEPS WITH AN OLD J, IF WE HAVE HAD A FAILURE OF ANY
C        KIND ON THE FIRST, SECOND OR THIRD STAGE OF THE CURRENT STEP
C **********************************************************************
      END IF

      IDOUB = L + 1
      CALL CPYARY(N*L,Y,YHOLD)

 60   IF (DABS(RC-ONE).GT.UPBND) IWEVAL = MITER
      HUSED = H
C     ------------------------------------------------------------------
C     THIS SECTION COMPUTES THE PREDICTED VALUES OF Y
C     AND THE RHS, ARH, FOR USE IN THE NEWTON ITERATION SCHEME.
C     RC IS THE RATIO OF THE NEW TO OLD VALUES OF THE COEFFICIENT
C     H*EL(1). WHEN RC DIFFERS FROM 1 BY MORE THAN 20 PERCENT, IWEVAL IS
C     SET TO MITER TO FORCE THE PARTIALS TO BE UPDATED.
C     ------------------------------------------------------------------
      QI = H*EL(1)
      QQ = ONE/QI
      DO 70 I = 1,N
         ARH(I) = EL(2)*Y(I,1)
 70   CONTINUE
      DO 90 J1 = 2,NQ
         JP1 =J1+1
         DO 80 I = 1,N
            ARH(I) = ARH(I) + EL(JP1)*Y(I,J1)
 80      CONTINUE
 90   CONTINUE
      IF (JCHANG.EQ.1) THEN
C        IF WE HAVE CHANGED STEPSIZE THEN PREDICT A VALUE FOR Y(T+H)
C        AND EVALUATE THE DERIVATIVE THERE (STORED IN SAVE2())
         CALL PRDICT(T,H,Y,L,N)
         IF(IERR.NE.0) GOTO 8000
         DO 95 I=1,N
            YPRIME(I)=(Y(I,1)-ARH(I))/QI
 95      CONTINUE

      ELSE
C        ELSE USE THE VALUES COMPUTED FOR THE SECOND BDF FROM THE LAST
C        STEP. Y( ,LMAX+2) HOLDS THE VALUE FOR THE DERIVATIVE AT (T+H)
C        AND Y( ,LMAX+3) HOLDS THE APPROXIMATION TO Y AT THIS POINT.
         LMP2=LMAX+2
         LMP3=LMAX+3
         DO 100 I = 1,N
            Y(I,1)=Y(I,LMP3)
            YPRIME(I) = (Y(I,1)-ARH(I))/QI
 100     CONTINUE
         T = T + H
      END IF

 110  IF (IWEVAL.LE.0) GO TO 120
C -------------------------------------------------------------------
C     IF INDICATED, THE MATRIX P = I/(H*EL(2)) - J IS RE-EVALUATED BEFORE
C     STARTING THE CORRECTOR ITERATION.  IWEVAL IS SET = 0 TO INDICATE
C     THAT THIS HAS BEEN DONE. P IS COMPUTED AND PROCESSED IN PSET.
C     THE PROCESSED MATRIX IS STORED IN PW
C -------------------------------------------------------------------
      IWEVAL = 0
      RC = ONE
      IITER = MEQC1 - MQ1TMP
      IITER2 = MEQC2 - MQ2TMP
      IF (JNEWIM) THEN
         IF (JSNOLD.GE.3) THEN
            AVNEWJ = TCRAT1/DBLE(FLOAT(IITER))
C FM            AVNEW2 = TCRAT2/DBLE(FLOAT(IITER2))

         ELSE
            AVNEWJ = ONE
c FM            AVNEW2 = ONE
         END IF

      ELSE
C
C          MATRIX P WAS FORMED WITH A COPY OF J
C
         IF (JSNOLD.GE.3) THEN
            AVOLDJ = TCRAT1/DBLE(FLOAT(IITER))
c            AVOLD2 = TCRAT2/DBLE(FLOAT(IITER2))
            IF (AVOLDJ.LT.AVNEWJ) THEN
               AVNEWJ = AVOLDJ

            ELSE IF (((DABS(AVOLDJ-AVNEWJ)).GT.0.3D+0) .OR.
     +          ((AVOLDJ.GT.0.85D+0).AND. (AVOLDJ.NE.ONE))) THEN
C
C              SINCE IN CERTAIN INSTANCES AVOLDJ WILL
C              BE 1.0 AND THERE WILL BE NO NEED TO
C              UPDATE J.
C
               CFAIL = .TRUE.
               CRATE1 = 0.1D+0
               CRATE2 = 0.1D+0
            END IF

         ELSE
            CFAIL = .TRUE.
            CRATE1 = 0.1D+0
            CRATE2 = 0.1D+0
C
C           *********************************************************
C           IF WE HAVE REACHED HERE THINGS MUST HAVE GONE WRONG
C           *********************************************************
C
         END IF

      END IF

      TCRAT1 = ZERO
      TCRAT2 = ZERO
      IF (CFAIL) THEN
         NRENEW = 1
         NEWPAR = 1
         JSINUP = -1
         JNEWIM = .TRUE.

      END IF

      CFAIL = .FALSE.
      JSNOLD = 0
      MQ1TMP = MEQC1
      MQ2TMP = MEQC2


      CALL PSET(Y,YPRIME,N,NPD,H,T,UROUND,EPSJAC,QI,MITER,MBND,
     +   IER,pderv,resid,NRENEW,YMAX,SAVE1,SAVE2,
     +   SCALE,PW,PWCOPY,ERROR,IPIV,ITOL,RTOL,ATOL,NPSET,NJE,NRE,NDEC
     +     ,IPAR,RPAR,IERR)
      IF(IERR.NE.0) GOTO 8000
      QQQ=QI
C
C     NOTE THAT ERROR() IS JUST BEING USED AS A WORKSPACE BY PSET
      IF (IER.NE.0) THEN
C     IF IER>0 THEN WE HAVE HAD A SINGULARITY IN THE ITERATION MATRIX
         IJUS=1
         RED=0.5D+0
         NFAIL = NFAIL + 1
         GO TO 450

      END IF


 120  DO 130 I = 1,N
         SAVE1(I) = Y(I,1)
         ERROR(I) = ZERO
 130  CONTINUE
      M1 = 0
C **********************************************************************
C     UP TO 4 CORRECTOR ITERATIONS ARE TAKEN.  A CONVERGENCE TEST IS MADE
C     ON THE R.M.S. NORM OF EACH CORRECTION ,USING BND, WHICH DEPENDS
C     ON ATOL AND RTOL.  THE SUM OF THE CORRECTIONS IS ACCUMULATED IN THE
C     VECTOR  ERROR(I).  THE Y ARRAY IS NOT ALTERED IN THE CORRECTOR
C     LOOP. THE UPDATED Y VECTOR IS STORED TEMPORARILY IN SAVE1.
C **********************************************************************
      IF (.NOT.SAMPLE) THEN

         CALL ITRAT2(QQQ,Y,YPRIME,N,T,QI,BND,ARH,CRATE1,TCRAT1,M1,
     +        WORKED,YMAX,ERROR,SAVE1,SAVE2,SCALE,PW,MF,MBND,
     +        NIND1,NIND2,NIND3,IPIV,1,ITOL,RTOL,ATOL,IPAR,RPAR,
     +        HUSED,NBSOL,NRE,NQUSED,resid,IERR)
         IF(IERR.NE.0) GOTO 8000

      ELSE
         CALL ITRAT2(QQQ,Y,YPRIME,N,T,QI,BND,ARH,CRATE1,TCRAT1,M1,
     +        WORKED,YMAX,ERROR,SAVE1,SAVE2,SCALE,PW,MF,MBND,
     +        NIND1,NIND2,NIND3,IPIV,0,ITOL,RTOL,ATOL,IPAR,RPAR,
     +        HUSED,NBSOL,NRE,NQUSED,resid,IERR)
         IF(IERR.NE.0) GOTO 8000
      END IF

      MEQC1 = MEQC1 + M1 + 1
C
C       NOW TEST TO SEE IF IT WAS SUCCESSFUL OR NOT
C
C
      IF (.NOT.WORKED) THEN
         NFAIL = NFAIL + 1
C **********************************************************************
C        THE CORRECTOR ITERATION FAILED TO CONVERGE IN 4 TRIES. IF
C        PARTIALS ARE NOT UP TO DATE, THEY ARE RE-EVALUATED FOR THE
C        NEXT TRY. OTHERWISE THE Y ARRAY IS REPLACED BY ITS VALUES
C        BEFORE PREDICTION AND H IS REDUCED IF POSSIBLE. IF NOT A
C        NON-CONVERGENCE EXIT IS TAKEN
C **********************************************************************
         IF (IWEVAL.EQ.-1) THEN
C           HAVE BEEN USING OLD PARTIALS, UPDATE THEM AND TRY AGAIN
            IWEVAL = MITER
            CFAIL = .TRUE.
            do 135 i=1,n
               yprime(i)=(y(I,1)-arh(i))*qq
 135        continue
            GO TO 110

         END IF

         IJUS=0
         RED=0.5D+0
C    ***    failed at step 1 because of Newton
         GO TO 450

      END IF

      IWEVAL = -1
      HUSED = H
      NQUSED = NQ
      DO 140 I = 1,N
         Y(I,1) = (SAVE1(I)-ARH(I))
 140  CONTINUE
      DO 145 I=1,N
         SAVE2(I) = Y(I,1)*QQ
         Y(I,1) = SAVE1(I)
 145  CONTINUE
C
C     UPDATE THE DIFFERENCES AT N+1
C
      DO 160 J = 2,L
         JM1 = J-1
         DO 150 I = 1,N
            Y(I,J) = Y(I,JM1) - YHOLD(I,JM1)
 150     CONTINUE
 160  CONTINUE
C
C     COMPUTE ERROR IN THE SOLUTION
C
      DO 161 I=1,N
         AYI = DABS(Y(I,1))
         IF(ITOL.EQ.1) THEN
            SCALE(I) = YMAX(I)
         ELSE IF(ITOL.EQ.2) THEN
            SCALE(I) = RTOL(1)*AYI + ATOL(1)
         ELSE IF(ITOL.EQ.3) THEN
            SCALE(I) = RTOL(1)*AYI + ATOL(I)
         ELSE IF(ITOL.EQ.4) THEN
            SCALE(I) = RTOL(I)*AYI + ATOL(1)
         ELSE IF(ITOL.EQ.5) THEN
            SCALE(I) = RTOL(I)*AYI + ATOL(I)
         ENDIF
 161  CONTINUE
      IF(NIND2.NE.0) THEN
         DO 162 I = NIND1+1,NIND2+NIND1
            SCALE(I)=SCALE(I)/HUSED
 162     CONTINUE
      ENDIF
      IF(NIND3.NE.0) THEN
         DO 163 I = NIND1 +NIND2 + 1,NIND1+NIND2+NIND3
            SCALE(I)=SCALE(I)/(HUSED**2)
 163     CONTINUE
      ENDIF
C ****
C ****  AMMEND
C ****  CHANGE 1,N BELOW TO 1,NVARS
C ****
      D = ZERO
      DO 170 I = 1,N
         D = D + ((Y(I,L)-YHOLD(I,L))/SCALE(I))**2
 170  CONTINUE
c
c    STORING Y FROM FIRST STEP FOR USE IN THIRD STEP.
C
      IF(ITOL .EQ. 1) D = D/(RTOL(1)**2)
      IF(D.GT.E) GOTO 330
      DO 180 I = 1,N
         YNHOLD(I,1) = Y(I,1)
         YNHOLD(I,2) = SAVE2(I)
 180  CONTINUE

      KFAIL = 0
      IREDO = 0
C ----------------------------------------------------------------------
      DO 190 I = 1,N
        ARH(I) = EL(2)*Y(I,1)
 190  CONTINUE
      DO 210 J1 = 2,NQ
         JP1 = J1+1
         DO 200 I = 1,N
            ARH(I) = ARH(I) + EL(JP1)*Y(I,J1)
 200     CONTINUE
 210  CONTINUE
      CALL PRDICT(T,H,Y,L,N)
      IF(IERR.NE.0) GOTO 8000

      DO 215 I=1,N
         YPRIME(I)=(Y(I,1)-ARH(I))/QQQ
 215  CONTINUE
      DO 220 I = 1,N
         SAVE1(I) = Y(I,1)
         ERROR(I) = ZERO
 220  CONTINUE
      M2 = 0
C
C     FOR NOW WILL ASSUME THAT WE DO NOT WISH TO SAMPLE
C     AT THE N+2 STEP POINT
C

      CALL ITRAT2(QQQ,Y,YPRIME,N,T,QI,BND,ARH,CRATE2,TCRAT2,M2,
     +     WORKED,YMAX,ERROR,SAVE1,SAVE2,SCALE,PW,MF,MBND,
     +     NIND1,NIND2,NIND3,IPIV,1,ITOL,RTOL,ATOL,IPAR,RPAR,
     +     HUSED,NBSOL,NRE,NQUSED,resid,IERR)
      IF(IERR.NE.0) GOTO 8000
      MEQC2 = MEQC2 + M2 + 1
C
C       NOW CHECK TO SEE IF IT WAS SUCCESSFUL OR NOT
C
      IF (.NOT.WORKED) THEN
         NFAIL = NFAIL + 1
         IJUS=0
         RED=0.5D+0
C ***have failed on step 2
         GOTO 450

      END IF
C
C        IF WE ARE DOWN TO HERE THEN THINGS MUST HAVE CONVERGED
C
      LMP2=LMAX+2
      LMP3=LMAX+3
      DO 230 I = 1,N
         Y(I,LMP3) = (SAVE1(I)-ARH(I))
 230  CONTINUE
      DO 233 I=1,N
         Y(I,LMP2) = Y(I,LMP3)*QQ
         Y(I,LMP3) = SAVE1(I)
 233  CONTINUE

C
C     WE ARE NOW COMPUTING THE THIRD STAGE
C
      LL = L + 1
      T = TOLD + H
      DELST = ELST(1)-EL(1)
      NQP2 = NQ+2
      LMP2=LMAX+2
      DO 280 I=1,N
         ARH(I) = H*(ELST(NQP2)*Y(I,LMP2)+DELST*YNHOLD(I,2))
         DO 270 J1 = 1,NQ
            ARH(I) = ARH(I) + ELST(J1+1)*YHOLD(I,J1)
 270     CONTINUE
 280  CONTINUE
      DO 290 I = 1,N
         SAVE2(I) = YNHOLD(I,2)
         Y(I,1) = YNHOLD(I,1)
 290  CONTINUE
      M3STEP = 0
 300  CONTINUE
      DO 310 I=1,N
         YPRIME(I)=(Y(I,1)-ARH(I))/QQQ
 310  CONTINUE
C KS:      CALL RESID(N,T,Y,SAVE1,YPRIME,IPAR,RPAR,IERR)
      CALL resid (T, Y, YPRIME, QI, SAVE1, IERR, RPAR, IPAR)
      if (ierr .ne. 0) go to 8000
      NRE=NRE+1
C
      IF (MF.GE. 23) THEN
         CALL DGBSL(PW,MBND(4),N,MBND(1),MBND(2),IPIV,SAVE1,0)
         NBSOL=NBSOL+1
      ELSE
         CALL SOL_mebdfi(N,N,PW,SAVE1,IPIV)
         NBSOL = NBSOL + 1
      ENDIF
      DO 321 I=1,N
         AYI = DABS(Y(I,1))
         IF(ITOL.EQ.1) THEN
            SCALE(I) = YMAX(I)
         ELSE IF(ITOL.EQ.2) THEN
            SCALE(I) = RTOL(1)*AYI + ATOL(1)
         ELSE IF(ITOL.EQ.3) THEN
            SCALE(I) = RTOL(1)*AYI + ATOL(I)
         ELSE IF(ITOL.EQ.4) THEN
            SCALE(I) = RTOL(I)*AYI + ATOL(1)
         ELSE IF(ITOL.EQ.5) THEN
            SCALE(I) = RTOL(I)*AYI + ATOL(I)
         ENDIF
 321  CONTINUE
      IF(NIND2.NE.0) THEN
         DO 322 I = NIND1+1,NIND2+NIND1
            SCALE(I)=SCALE(I)/HUSED
 322     CONTINUE
      ENDIF
      IF(NIND3.NE.0) THEN
         DO 133 I = NIND1+NIND2 + 1,NIND1+NIND2+NIND3
            SCALE(I)=SCALE(I)/(HUSED**2)
 133     CONTINUE
      ENDIF
      D = ZERO
      DO 320 I = 1,N
        D = D + (SAVE1(I)/SCALE(I))**2
        Y(I,1) = Y(I,1) - SAVE1(I)
 320  CONTINUE
      IF(ITOL .EQ. 1) D = D/(RTOL(1)**2)
      IF ((D*DMIN1(ONE,2.0D+0*CRATE1)).LE.BND) GO TO 360
      IF (M3STEP.EQ.4) THEN
         IJUS=1
         RED=0.5D+0
C    ****  step 3 fails
         NFAIL = NFAIL + 1
         GO TO 450
      END IF

      M3STEP = M3STEP + 1
c      IF(IERR.NE.0) GOTO 8000
      GO TO 300

 330  KFAIL = KFAIL - 1
C **********************************************************************
C     THE ERROR TEST FAILED. KFAIL KEEPS TRACK OF MULTIPLE FAILURES.
C     RESTORE T AND THE Y ARRAY TO THEIR PREVIOUS VALUES AND PREPARE TO
C     TRY THE STEP AGAIN. COMPUTE THE OPTIMAL STEP SIZE FOR THIS ORDER
C     AND ONE ORDER LOWER.
C **********************************************************************
C     ***  failed on step 1 because of accuracy
C     COMPUTE ERROR IN THE SOLUTION
C
      NFAIL = NFAIL + 1
      DO 561 I=1,N
         AYI = DABS(Y(I,1))
         IF(ITOL.EQ.1) THEN
            SCALE(I) = YMAX(I)
         ELSE IF(ITOL.EQ.2) THEN
            SCALE(I) = RTOL(1)*AYI + ATOL(1)
         ELSE IF(ITOL.EQ.3) THEN
            SCALE(I) = RTOL(1)*AYI + ATOL(I)
         ELSE IF(ITOL.EQ.4) THEN
            SCALE(I) = RTOL(I)*AYI + ATOL(1)
         ELSE IF(ITOL.EQ.5) THEN
            SCALE(I) = RTOL(I)*AYI + ATOL(I)
         ENDIF
 561  CONTINUE
      IF(NIND2.NE.0) THEN
         DO 562 I = NIND1+1,NIND2+NIND1
            SCALE(I)=SCALE(I)/HUSED
 562     CONTINUE
      ENDIF
      IF(NIND3.NE.0) THEN
         DO 563 I = NIND1+NIND2 + 1,NIND1+NIND2+NIND3
            SCALE(I)=SCALE(I)/(HUSED**2)
 563     CONTINUE
      ENDIF
      DDOWN = ZERO
      TWODWN = ZERO
      DO 1700 I=1,N
         DDOWN = DDOWN + ((Y(I,L))/SCALE(I))**2
         TWODWN = TWODWN + ((Y(I,l-1))/SCALE(I))**2
 1700 CONTINUE
      IF(ITOL .EQ. 1) D = D/(RTOL(1)**2)
      T = TOLD
      HOLD = H
C FM: added FFAIL=0.0d0 to avoid uninitialised value(s)
      FFAIL = 0.5d0
      FRFAIL = 0.5d0
      IF(NQ.GT.1) FFAIL = 0.5D+0/DBLE(FLOAT(NQ))
      IF(NQ.GT.2) FRFAIL = 0.5D+0/DBLE(FLOAT(NQ-1))
      EFAIL = 0.5D+0/DBLE(FLOAT(L))
      CALL CPYARY(N*L,YHOLD,Y)
      RMAX = 2.0D+0
      IF (DABS(H).LE.HMIN*1.00001D+0) THEN
C
C        REQUESTED ERROR NOT POSSIBLE WITH GIVEN HMIN
C
         KFLAG = -1
         HOLD = H
         RETURN

      END IF

      IF (KFAIL.LE.-3) GO TO 340
      IREDO = 2
C
C     PREDICTING A NEW H AFTER INSUFFICIENT ACCURACY
C
      PRFAIL = ((D/(0.2D+0*E))**EFAIL)*1.5D+0 + 1.6D-6
      PLFAIL = ((DDOWN/(0.2D+0*EDN))**FFAIL)*1.5D+0+1.7D-6
C FM:   added the following line to avoid  uninitialised value(s)
      IF(NQ.GT.2) THEN
         PLLFAL =((TWODWN/(0.2D+0*EDDN))**FRFAIL)*
     +     1.5D+0+1.7d-6
         IF(PLLFAL.GT.PLFAIL) PLFAIL=PLLFAL
      ENDIF
c      PLLFAL = PLFAIL
c      IF(NQ.GT.2) PLLFAL =((TWODWN/(0.2D+0*EDDN))**FRFAIL)*
c     +     1.5D+0+1.7d-6
C FM added NQ.GT.2 in the if
c      IF(PLLFAL.GT.PLFAIL) PLFAIL=PLLFAL
      IF(PLFAIL.LT.PRFAIL.AND.NQ.NE.1) THEN
         NEWQ=NQ-1
         NQ=NEWQ
         RH=ONE/(PLFAIL*DBLE(FLOAT(-KFAIL)))
         L=NQ+1
         CALL COSET(NQ,EL,ELST,TQ,NCOSET,MAXORD)
         RC=RC*EL(1)/OLDLO
         OLDLO=EL(1)
         CALL ERRORS(N,TQ,EDN,E,EUP,BND,EDDN)
      ELSE
         NEWQ = NQ
         RH = ONE/ (PRFAIL*DBLE(FLOAT(-KFAIL)))
      ENDIF
      GO TO 40
C **********************************************************************
C     CONTROL REACHES THIS STAGE IF 3 OR MORE FAILURES HAVE OCCURED.
C     IT IS ASSUMED THAT THE DERIVATIVES THAT HAVE ACCUMULATED IN THE Y
C     ARRAY HAVE ERRORS OF THE WRONG ORDER. HENCE THE FIRST DERIVATIVE
C     IS RE-COMPUTED, AND THE ORDER IS SET TO 1. THEN H IS REDUCED BY A
C     FACTOR OF 10, AND THE STEP IS RETRIED. AFTER A TOTAL OF 7
C     FAILURES AN EXIT IS TAKEN WITH KFLAG=-2.
C **********************************************************************
 340  IF (KFAIL.EQ.-7) THEN
C     ERROR SMALLER THAN CAN BE HANDLED FOR PROBLEM
         KFLAG = -2
         HOLD = H
         RETURN

      END IF
C     *********************************
C     START FROM ORDER 1 AGAIN    *
C     *********************************
      JCHANG = 1
      RH = DMAX1(HMIN/DABS(H),0.1D+0)
      CALL HCHOSE(RH,H,OVRIDE,HSTPSZ)
      H = H*RH
      DO 350 I = 1,N
         Y(I,1) = YHOLD(I,1)
C FM added this line to be equal to the testset code
         YHOLD(I,2) = YHOLD(I,2)*RH
         Y(I,2) = YHOLD(I,2)
 350  CONTINUE
      IWEVAL = MITER
      CFAIL = .TRUE.
C     SINCE WE HAVE HAD PROBLEMS PROCEED WITH THIS ORDER
C     FOR 10 STEPS (IF WE CAN)
      IDOUB = 10
      IF (NQ.EQ.1) GO TO 60
      NQ = 1
      L = 2
C     RESET ORDER, RECALCULATE ERROR BOUNDS
      CALL COSET(NQ,EL,ELST,TQ,NCOSET,MAXORD)
      LMAX = MAXDER + 1
      RC = RC*EL(1)/OLDLO
      OLDLO = EL(1)
      CALL ERRORS(N,TQ,EDN,E,EUP,BND,EDDN)
C     NOW JUMP TO NORMAL CONTINUATION POINT
      GO TO 60
C **********************************************************************
C     THE ITERATION FOR THE CORRECTED SOLUTION HAS CONVERGED.
C     UPDATE THE Y ARRAY.
C **********************************************************************
 360  CONTINUE
C   ****
C   **** AMMEND ****
C   **** CHANGE 1,N BELOW TO 1,NVARS
C   ****
      DEMB=0.0D+0
      DO 361 I=1,N
         DEMB=DEMB+((Y(I,1)-YNHOLD(I,1))/SCALE(I))**2
 361  CONTINUE
      IF(DEMB.GT.4.0D+0*DBLE(FLOAT(N))) THEN
         IEMB=1
         IJUS=1
         RED=0.5D+0
C ***  failed because of embedded error estimate
         NFAIL = NFAIL + 1
         GOTO 450
      ENDIF
      DO 380 J2 = 2,LL
         J2M1=J2-1
         DO 370 I = 1,N
            Y(I,J2) = Y(I,J2M1) - YHOLD(I,J2M1)
 370     CONTINUE
 380  CONTINUE
      do 385 i=1,n
         yprime(i)=(y(i,1)-arh(i))/qqq
 385  continue
C ---------------------------------------------------------------------
C     IF THE COUNTER IDOUB EQUALS 2 AND WE ARE NOT ALREADY USING THE
C     MAXIMUM ALLOWABLE ORDER , STORE Y(I,LMAX+4) WHICH IS USED IN
C     ASSESSING THE POSSIBILITY OF INCREASING THE ORDER. IF IDOUB = 0
C     CONTROL PASSES TO 480 WHERE AN ATTEMPT TO CHANGE THE STEPSIZE AND
C     ORDER IS MADE.
C ----------------------------------------------------------------------
C FM: moved this instruction outside the loop
      LMP4=LMAX+4
      IF (IDOUB.EQ.2.AND.L.NE.LMAX) THEN
         DO 390 I = 1,N
            Y(I,LMP4) = Y(I,LL)
 390     CONTINUE
      END IF
      IDOUB = IDOUB - 1
      TRANGE=(TEND-TOLD-H)*H
      IF(TRANGE.LT.0.0D+0) THEN
         IDOUB = IDOUB + 2
         GOTO 440
      ENDIF
      JCHANG = 0
      IF (IDOUB.EQ.0) THEN
         SAMPLE = .FALSE.
         ISAMP = ISAMP + 1
         IF (ISAMP.EQ.4) THEN
            SAMPLE = .TRUE.
            ISAMP = 0
         END IF

C **********************************************************************
C        NOW COMPUTE THE FACTORS PR1, PR2 AND PR3, BY WHICH
C        H COULD BE DIVIDED AT ORDER NQ-1, ORDER NQ AND ORDER NQ+1
C        RESPECTIVELY. THE SMALLEST OF THESE IS DETERMINED AND THE NEW
C        ORDER CHOSEN ACCORDINGLY. IF THE ORDER IS TO BE INCREASED WE
C        MUST COMPUTE ONE MORE BACKWARD DIFFERENCE.
C **********************************************************************
         PR3 = 1.D+20
         FAC = 1.5D+0
         IF(IEMB.EQ.1) FAC = 1.8D+0
         VHOLD = 0.0d0
         DO 400 I = 1,N
            AYI = DABS(Y(I,1))
            IF(ITOL.EQ.1) THEN
               VHOLD = YMAX(I)
            ELSE IF(ITOL.EQ.2) THEN
               VHOLD = RTOL(1)*AYI + ATOL(1)
            ELSE IF(ITOL.EQ.3) THEN
               VHOLD = RTOL(1)*AYI + ATOL(I)
            ELSE IF(ITOL.EQ.4) THEN
               VHOLD = RTOL(I)*AYI + ATOL(1)
            ELSE IF(ITOL.EQ.5) THEN
               VHOLD = RTOL(I)*AYI + ATOL(I)
            ENDIF
            SCALE(I)=VHOLD
 400     CONTINUE
         IF(NIND2.NE.0) THEN
            DO 4461 I=NIND1+1,NIND1+NIND2
               SCALE(I) = SCALE(I)/HUSED
 4461       CONTINUE
         ENDIF
         IF(NIND3.NE.0) THEN
            DO 4462 I=NIND1+NIND2+1,N
               SCALE(I) = SCALE(I)/(HUSED**2)
 4462       CONTINUE
         ENDIF
         IF(L.NE.LMAX) THEN
            LMP4 = LMAX + 4
            DUP = ZERO
            DO 401 I=1,N
               DUP = DUP + ((Y(I,LL)-Y(I,LMP4))/SCALE(I))**2
 401        CONTINUE
            IF(ITOL .EQ. 1) DUP = DUP/(RTOL(1)**2)
            ENQ3 = 0.5D+0/DBLE(FLOAT(L+1))
            PR3 = ((DUP/EUP)**ENQ3)*(FAC+0.2D+0) + 1.8D-6
         END IF

         ENQ2 = 0.5D+0/DBLE(FLOAT(L))
         D = ZERO
         DDOWN=ZERO

         DO 410 I = 1,N
            D = D + (Y(I,LL)/SCALE(I))**2
            DDOWN = DDOWN + (Y(I,LMP4)/SCALE(I))**2
 410     CONTINUE
         IF(ITOL .EQ.1) D = D/(RTOL(1)**2)
         PR2 = ((D/E)**ENQ2)*FAC + 1.6D-6
         PR1 = 1.D+20
         IF (NQ.GT.1) THEN
            DDDOWN=ZERO
            DDOWN = ZERO
            DO 57420 I = 1,N
               DDOWN = DDOWN + (Y(I,L)/SCALE(I))**2
               DDDOWN = DDDOWN + (Y(I,L-1)/SCALE(I))**2
57420       CONTINUE
            IF(ITOL .EQ. 1) DDOWN = DDOWN/(RTOL(1)**2)
            ENQ1 = 0.5D+0/DBLE(FLOAT(NQ))
            PR1 = ((DDOWN/EDN)**ENQ1)*(FAC+0.1D+0) + 1.7D-6
            IF(NQ.GT.2) THEN
               ENQ0 = 0.5D+0/DBLE(FLOAT(NQ-1))
               PR0 = ((DDDOWN/EDDN)**ENQ0)*(FAC+0.1D+0) + 1.7D-6
               IF(PR0.GT.PR1) PR1 = PR0
               IF(DDDOWN.LT.DDOWN) DDOWN = DDDOWN
            ENDIF
         END IF
         IF(L.EQ.LMAX) DUP = 0.0D+0
         IF(NQ.LE.1) GOTO 6578
         IF(DUP.GT.D.AND.D.GT.DDOWN) THEN
            PR2=1.0D+30
            PR3=1.0D+30
         ENDIF
 6578    CONTINUE
         IF (PR2.LE.PR3) THEN
            IF (PR2.GT.PR1) THEN
               NEWQ = NQ - 1
               RH = 1.0D+0/PR1

            ELSE
               NEWQ = NQ
               RH = 1.0D+0/PR2
            END IF

         ELSE IF (PR3.LT.PR1) THEN
            NEWQ = L
            RH = 1.0D+0/PR3

         ELSE
            NEWQ = NQ - 1
            RH = 1.0D+0/PR1
         END IF
         IEMB=0
         IF(RH.GT.1.0D+0.AND.RH.LT.1.1D+0) THEN
            IDOUB=10
            NQ=NQUSED
            L=NQ+1
            GOTO 440
         ENDIF
         RH = DMIN1(RH,RMAX)
         CALL HCHOSE(RH,H,OVRIDE,HSTPSZ)
         IF ((JSINUP.LE.20).AND.(KFLAG.EQ.0).AND.(RH.LT.1.1D+0)) THEN
C           WE HAVE RUN INTO PROBLEMS
            IDOUB = 10
            NQ = NQUSED
            L = NQ + 1
            GO TO 440

         END IF
C **********************************************************************
C        IF THERE IS A CHANGE IN ORDER, RESET NQ, L AND THE
C        COEFFICIENTS. IN ANY CASE H IS RESET  AND THE
C        Y ARRAY IS RE-SCALED
C **********************************************************************
         IF(NIND3.NE.0) THEN
            IF(NQ.LE.2.AND.PR3.LT.1.0D+0) NEWQ=NQ+1
         ENDIF
         IF (NEWQ.NE.NQ) THEN
            IF (NEWQ.GT.NQ) THEN
C              ADD AN EXTRA TERM TO THE HISTORY ARRAY
               DO 430 I = 1,N
                  Y(I,LL) = Y(I,L) - YHOLD(I,L)
 430           CONTINUE
            END IF

            NQ = NEWQ
            L = NQ + 1
C           RESET ORDER,RECALCULATE ERROR BOUNDS
            CALL COSET(NQ,EL,ELST,TQ,NCOSET,MAXORD)
            LMAX = MAXDER + 1
            RC = RC*EL(1)/OLDLO
            OLDLO = EL(1)
            CALL ERRORS(N,TQ,EDN,E,EUP,BND,EDDN)
         END IF

         RH = DMAX1(RH,HMIN/DABS(H))
         RH = DMIN1(RH,HMAX/DABS(H),RMAX)
         CALL RSCALE(N,L,RH,Y)
         RMAX = 10.0D+0
         JCHANG = 1
         H = H*RH
         RC = RC*RH
         IF(JSNOLD.GT.IBND) RC=ZERO

         IDOUB = L + 1

      END IF

 440  CONTINUE
C ----------------------------------------------------------------------
C     STORE THE Y ARRAY IN THE MATRIX YHOLD.  STORE IN THE Y ARRAY THE
C     INFORMATION NECESSARY TO PERFORM AN INTERPOLATION TO FIND THE
C     SOLUTION AT THE SPECIFIED OUTPUT POINT IF APPROPRIATE.
C ----------------------------------------------------------------------
      CALL CPYARY(N*L,Y,YHOLD)
      NSTEP = NSTEP + 1
      JSINUP = JSINUP + 1
      JSNOLD = JSNOLD + 1
      JSTART = NQUSED
      T = TOLD + HUSED
      HOLD = H
      KFAIL = 0
      NEWPAR = 0
      CFAIL = .FALSE.
      RETURN
 450  CONTINUE
      FINISH = .FALSE.
      T=TOLD
      RMAX=2.0D+0
      DO 465 J1=1,L
         DO 460 I=1,N
            Y(I,J1)=YHOLD(I,J1)
 460  CONTINUE
 465  CONTINUE
      IF(DABS(H).LE.HMIN*1.00001D+0) THEN
C
C   CORRECTOR CONVERGENCE COULD NOT BE ACHIEVED
C

         IF(NSTEP.EQ.0) THEN
            KFLAG=-1
         ELSE
            KFLAG=-3
         END IF
C
C    TO SUPPRESS ERROR MESSAGES AT START AS H MAY
C    HAVE BEEN TOO LARGE ON THE FIRST STEP.
C
         HOLD=H
         FINISH = .TRUE.
      END IF
      RH = RED
      IREDO=1
C
C     TRY AGAIN WITH UPDATED PARTIALS
C
 8000 if (ierr .ne. 0) then

      CALL Rprint('IERR is non-zero due to an illegal function call')
         h= h/2
         IF(H.LT.EPSJAC/100.0D+0) THEN

      CALL Rprint('Stepsize is too small')
            KFLAG = -7
            RETURN
         ENDIF
         T= TOLD
C FRANCESCA MAZZIA commented the following IF block
C   y0 H0 and T0 are not used in the function stiff
c    added in input
         IF ((T-TOUT)*H.GE.0.0D+0) THEN
C           HAVE OVERSHOT TOUT

c            CALL INTERP(N,JSTART,H,T,Y,TOUT,Y0)
c            HO = H
c            T0 = TOUT
            KFLAG = -5
            RETURN
         ENDIF
         IERR = 0
         jstart = -1
         goto 30
      endif
c
      IF(IJUS.EQ.0) CALL HCHOSE(RH,H,OVRIDE,HSTPSZ)
      IF(.NOT.FINISH) THEN
         GO TO 40
      ELSE
         RETURN
      END IF

      END
C ------------------- END OF SUBROUTINE STIFF --------------------------

      SUBROUTINE RSCALE(N,L,RH,Y)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     .. SCALAR ARGUMENTS ..
      INTEGER L,N
C     ..
C     .. ARRAY ARGUMENTS ..
      DIMENSION  Y(N,12)
C     ..
C     .. LOCAL SCALARS ..
      INTEGER I,J,J1
C     ..
C     .. LOCAL ARRAYS ..
      DIMENSION  DI(8,8)
C     ..
C     .. DATA STATEMENTS ..

C     SUBROUTINE IS FOR RESCALING THE HISTORY ARRAY AFTER A CHANGE IN
C     STEPSIZE
C
C     N      ORDER OF THE PROBLEM
C     L      NUMBER OF TERMS IN THE HISTORY ARRAY TO BE RESCALED
C     RH     RATIO OF THE STEPSIZE CHANGE (I.E. RH = HNEW/HOLD)
C     Y()    THE HISTORY ARRAY
C

      DATA  ZERO/0.0D+0/
C     ..

      DI(2,2) = RH
      IF (L.GT.2) THEN
         TA = RH*RH
         DI(2,3) = RH* (1.0D+0-RH)/2.0D+0
         DI(3,3) = TA
         IF (L.GT.3) THEN
            TB = TA*RH
            DI(2,4) = RH* ((RH-3.0D+0)*RH+2.0D+0)/6.0D+0
            DI(3,4) = TA* (1.0D+0-RH)
            DI(4,4) = TB
            IF (L.GT.4) THEN
               TC = TB*RH
               DI(2,5) = - (((RH-6.0D+0)*RH+11.0D+0)*RH-6.0D+0)*RH/
     +                      24.0D+0
               DI(3,5) = TA* ((7.0D+0*RH-18.0D+0)*RH+11.0D+0)/12.0D+0
               DI(4,5) = 1.5D+0*TB* (1.0D+0-RH)
               DI(5,5) = TC
               IF (L.GT.5) THEN
                  TD = TC*RH
                  DI(2,6) = ((((RH-10.0D+0)*RH+35.0D+0)*RH-50.0D+0)
     +                        *RH+24.0D+0)*RH/120.0D+0
                  DI(3,6) = - (((3.0D+0*RH-14.0D+0)*RH+21.0D+0)*RH
     +                       -10.0D+0)*TA/12.0D+0
                  DI(4,6) = ((5.0D+0*RH-12.0D+0)*RH+7.0D+0)*TB/4.0D+0
                  DI(5,6) = 2.0D+0*TC* (1.0D+0-RH)
                  DI(6,6) = TD
                  IF (L.GT.6) THEN
                     TE = TD*RH
                     DI(2,7) = -RH* (RH-1.0D+0)* (RH-2.0D+0)*
     +                          (RH-3.0D+0)*(RH-4.0D+0)*(RH-5.0D+0)/
     +                           720.0D+0
                     DI(3,7) = TA* ((((62.0D+0*RH-450.0D+0)*RH+
     +                         1190.0D+0)*RH-1350.0D+0)*RH+548.0D+0)
     +                         /720.0D+0
                     DI(4,7) = TB* (((-18.0D+0*RH+75.0D+0)*RH
     +                          -102.0D+0)*RH+45.0D+0)/24.0D+0
                     DI(5,7) = TC* ((13.0D+0*RH-30.0D+0)*RH+17.0D+0)
     +                         /6.0D+0
                     DI(6,7) = 2.5D+0*TD* (1.0D+0-RH)
                     DI(7,7) = TE
                     IF (L.GT.7) THEN
                        TF = TE*RH
                        DI(2,8) = RH*(RH-1.0D+0)*(RH-2.0D+0)*(RH
     +                            -3.0D+0)*(RH-4.0D+0)*(RH-5.0D+0)
     +                            *(RH-6.0D+0)/5040.0D+0
                        DI(3,8) = TA* ((((((-126.0D+0*RH)+1302.0D+0)*RH-
     +                            5250.0D+0)*RH+10290.0D+0)*RH-9744.0D+0
     +                            )*RH+3528.0D+0)/5040.0D+0
                        DI(4,8) = TB* ((((43.0D+0*RH-270.0D+0)*RH+
     +                            625.0D+0)*RH-630.0D+0)*RH+232.0D+0)
     +                           /120.0D+0
                        DI(5,8) = TC* (((-10.0D+0*RH+39.0D+0)*RH-
     +                            50.0D+0)*RH+21.0D+0)/6.0D+0
                        DI(6,8) = TD* ((20.0D+0*RH-45.0D+0)*RH+25.0D+0
     +                            )/6.0D+0
                        DI(7,8) = 3.0D+0*TE* (1.0D+0-RH)
                        DI(8,8) = TF
                     END IF

                  END IF

               END IF

            END IF

         END IF

      END IF

      DO 30 I = 1,N
        DO 20 J = 2,L
          ZZ = ZERO
          DO 10 J1 = J,L
            ZZ = ZZ + DI(J,J1)*Y(I,J1)
   10     CONTINUE
          Y(I,J) = ZZ
   20   CONTINUE
   30 CONTINUE
      RETURN

      END
C---------------------------------------------------------------------------

      SUBROUTINE CPYARY(NELEM,SOURCE,TARGET)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     COPIES THE ARRAY SOURCE() INTO THE ARRAY TARGET()
C
C     THIS SUBROUTINE COULD BE REPLACED BY THE BLAS ROUTINE SCOPY
C     (AFTER CHANGING THE ARGUMENT LIST APPROPRIATELY)
C
C     .. SCALAR ARGUMENTS ..
      INTEGER NELEM
C     ..
C     .. ARRAY ARGUMENTS ..
      DIMENSION  SOURCE(NELEM),TARGET(NELEM)
C     ..
C     .. LOCAL SCALARS ..
      INTEGER I
C     ..
      DO 10 I = 1,NELEM
         TARGET(I) = SOURCE(I)
 10   CONTINUE
      RETURN

      END
C----------------------------------------------------------------------------

      SUBROUTINE HCHOSE(RH,H,OVRIDE,HSTPSZ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION HSTPSZ(2,14)
c      COMMON / STPSZE / HSTPSZ
      LOGICAL OVRIDE
C
C     FIRST MOVE ALL ELEMENTS DOWN ONE PLACE
C
      IF (H.NE.HSTPSZ(2,1)) THEN
         DO 5 I=12,2,-1
            I2=I-1
            HSTPSZ(1,I)=HSTPSZ(1,I2)
            HSTPSZ(2,I)=HSTPSZ(2,I2)
 5       CONTINUE
C
C          NOW INSERT VALUE OF H USED BEFORE THIS CALL
C
         HSTPSZ(1,2)=H/HSTPSZ(2,1)
         HSTPSZ(2,1)=H
      END IF
C
C     NOW DECIDE ON THE NEW CHANGE
C
      IF (RH.GT.1.0D+0) THEN
         OVRIDE=.FALSE.
      ELSE IF (HSTPSZ(1,2).LE.1.0D+0) THEN
         OVRIDE=.FALSE.
      ELSE IF ((RH*H).LE.HSTPSZ(2,2)) THEN
         OVRIDE=.FALSE.
      ELSE
         RH=HSTPSZ(2,2)/H
         OVRIDE=.TRUE.
      END IF
      HSTPSZ(1,1)=RH

      RETURN
      END
C
C  ************************************************************
C      DOUBLE PRECISION FUNCTION DLAMCH( CMACH ) - removed

*
************************************************************************
*
      SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
*     ..
*
*  Purpose
*  =======
*
*  DLAMC1 determines the machine parameters given by BETA, T, RND, and
*  IEEE1.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  IEEE1   (output) LOGICAL
*          Specifies whether rounding appears to be done in the IEEE
*          'round to nearest' style.
*
*  Further Details
*  ===============
*
*  The routine is based on the routine  ENVRON  by Malcolm and
*  incorporates suggestions by Gentleman and Marovich. See
*
*     Malcolm M. A. (1972) Algorithms to reveal properties of
*        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
*
*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
*        that reveal properties of floating point arithmetic units.
*        Comms. of the ACM, 17, 276-277.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      DOUBLE PRECISION   A, B, C, F, ONE, QTR, SAVEC, T1, T2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Save statement ..
      SAVE        FIRST, LIEEE1, LBETA, LRND, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ONE = 1
*
*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
*        IEEE1, T and RND.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are  stored and not held in registers,  or
*        are not affected by optimizers.
*
*        Compute  a = 2.0**m  with the  smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.
*
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 10
         END IF
*+       END WHILE
*
*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.
*
         B = 1
         C = DLAMC3( A, B )
*
*+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = DLAMC3( A, B )
            GO TO 20
         END IF
*+       END WHILE
*
*        Now compute the base.  a and c  are neighbouring floating point
*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).
*
         QTR = ONE / 4
         SAVEC = C
         C = DLAMC3( C, -A )
         LBETA = INT(C + QTR)
*
*        Now determine whether rounding or chopping occurs,  by adding a
*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
*
         B = LBETA
         F = DLAMC3( B / 2, -B / 100 )
         C = DLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = DLAMC3( B / 2, B / 100 )
         C = DLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) )
     $      LRND = .FALSE.
*
*        Try and decide whether rounding is done in the  IEEE  'round to
*        nearest' style. B/2 is half a unit in the last place of the two
*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
*        A, but adding B/2 to SAVEC should change SAVEC.
*
         T1 = DLAMC3( B / 2, A )
         T2 = DLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
*
*        Now find  the  mantissa, t.  It should  be the  integer part of
*        log to the base beta of a,  however it is safer to determine  t
*        by powering.  So we find t as the smallest positive integer for
*        which
*
*           fl( beta**t + 1.0 ) = 1.0.
*
         LT = 0
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 30
         END IF
*+       END WHILE
*
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      RETURN
*
*     End of DLAMC1
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      DOUBLE PRECISION   EPS, RMAX, RMIN
*     ..
*
*  Purpose
*  =======
*
*  DLAMC2 determines the machine parameters specified in its argument
*  list.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  EPS     (output) DOUBLE PRECISION
*          The smallest positive number such that
*
*             fl( 1.0 - EPS ) .LT. 1.0,
*
*          where fl denotes the computed value.
*
*  EMIN    (output) INTEGER
*          The minimum exponent before (gradual) underflow occurs.
*
*  RMIN    (output) DOUBLE PRECISION
*          The smallest normalized number for the machine, given by
*          BASE**( EMIN - 1 ), where  BASE  is the floating point value
*          of BETA.
*
*  EMAX    (output) INTEGER
*          The maximum exponent before overflow occurs.
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest positive number for the machine, given by
*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
*          value of BETA.
*
*  Further Details
*  ===============
*
*  The computation of  EPS  is based on a routine PARANOIA by
*  W. Kahan of the University of California at Berkeley.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,
     $                   NGNMIN, NGPMIN
      DOUBLE PRECISION   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE,
     $                   SIXTH, SMALL, THIRD, TWO, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC1, DLAMC4, DLAMC5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Save statement ..
      SAVE       IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX,
     $                   LRMIN, LT
      SAVE     FIRST
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /

*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ZERO = 0
         ONE = 1
         TWO = 2
*
*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
*        BETA, T, RND, EPS, EMIN and RMIN.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are stored  and not held in registers,  or
*        are not affected by optimizers.
*
*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
*
         CALL DLAMC1( LBETA, LT, LRND, LIEEE1 )
*
*        Start to find EPS.
*
         B = LBETA
         A = B**( -LT )
         LEPS = A
*
*        Try some tricks to see whether or not this is the correct  EPS.
*
         B = TWO / 3
         HALF = ONE / 2
         SIXTH = DLAMC3( B, -HALF )
         THIRD = DLAMC3( SIXTH, SIXTH )
         B = DLAMC3( THIRD, -HALF )
         B = DLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS )
     $      B = LEPS
*
         LEPS = 1
*
*+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = DLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = DLAMC3( HALF, -C )
            B = DLAMC3( HALF, C )
            C = DLAMC3( HALF, -B )
            B = DLAMC3( HALF, C )
            GO TO 10
         END IF
*+       END WHILE
*
         IF( A.LT.LEPS )
     $      LEPS = A
*
*        Computation of EPS complete.
*
*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
*        Keep dividing  A by BETA until (gradual) underflow occurs. This
*        is detected when we cannot recover the previous A.
*
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = DLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE

         A = DLAMC3( ONE, SMALL )
         CALL DLAMC4( NGPMIN, ONE, LBETA )
         CALL DLAMC4( NGNMIN, -ONE, LBETA )
         CALL DLAMC4( GPMIN, A, LBETA )
         CALL DLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
*
         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
*            ( Non twos-complement machines, no gradual underflow;
*              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
*            ( Non twos-complement machines, with gradual underflow;
*              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
*            ( Twos-complement machines, no gradual underflow;
*              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND.
     $            ( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
*            ( Twos-complement machines with gradual underflow;
*              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
*         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
***
* Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
      CALL Rprinti1('Warning. The value emin may be incorrect: ',LEMIN)

         END IF
***
*
*        Assume IEEE arithmetic if we found denormalised  numbers above,
*        or if arithmetic seems to round in the  IEEE style,  determined
*        in routine DLAMC1. A true IEEE machine should have both  things
*        true; however, faulty machines may have one or the other.
*
         IEEE = IEEE .OR. LIEEE1
*
*        Compute  RMIN by successive division by  BETA. We could compute
*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
*        this computation.
*
         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = DLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
*
*        Finally, call DLAMC5 to compute EMAX and RMAX.
*
         CALL DLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
*
      RETURN
*
*
*     End of DLAMC2
*
      END
*
************************************************************************
*
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
*     ..
*
*  Purpose
*  =======
*
*  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
*  the addition of  A  and  B ,  for use in situations where optimizers
*  might hold one of these in a register.
*
*  Arguments
*  =========
*
*  A, B    (input) DOUBLE PRECISION
*          The values A and B.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      DLAMC3 = A + B
*
      RETURN
*
*     End of DLAMC3
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC4( EMIN, START, BASE )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      DOUBLE PRECISION   START
*     ..
*
*  Purpose
*  =======
*
*  DLAMC4 is a service routine for DLAMC2.
*
*  Arguments
*  =========
*
*  EMIN    (output) EMIN
*          The minimum exponent before (gradual) underflow, computed by
*          setting A = START and dividing by BASE until the previous A
*          can not be recovered.
*
*  START   (input) DOUBLE PRECISION
*          The starting point for determining EMIN.
*
*  BASE    (input) INTEGER
*          The base of the machine.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Executable Statements ..
*
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = DLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
*+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND.
     $    ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = DLAMC3( A / BASE, ZERO )
         C1 = DLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = DLAMC3( A*RBASE, ZERO )
         C2 = DLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
*+    END WHILE
*
      RETURN
*
*     End of DLAMC4
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      DOUBLE PRECISION   RMAX
*     ..
*
*  Purpose
*  =======
*
*  DLAMC5 attempts to compute RMAX, the largest machine floating-point
*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
*  approximately to a power of 2.  It will fail on machines where this
*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
*  EMAX = 28718).  It will also fail if the value supplied for EMIN is
*  too large (i.e. too close to zero), probably with overflow.
*
*  Arguments
*  =========
*
*  BETA    (input) INTEGER
*          The base of floating-point arithmetic.
*
*  P       (input) INTEGER
*          The number of base BETA digits in the mantissa of a
*          floating-point value.
*
*  EMIN    (input) INTEGER
*          The minimum exponent before (gradual) underflow.
*
*  IEEE    (input) LOGICAL
*          A logical flag specifying whether or not the arithmetic
*          system is thought to comply with the IEEE standard.
*
*  EMAX    (output) INTEGER
*          The largest exponent before overflow
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest machine floating-point number.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      DOUBLE PRECISION   OLDY, RECBAS, Y, Z
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     First compute LEXP and UEXP, two powers of 2 that bound
*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
*     approximately to the bound that is closest to abs(EMIN).
*     (EMAX is the exponent of the required number RMAX).
*
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
*
*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
*     than or equal to EMIN. EXBITS is the number of bits needed to
*     store the exponent.
*
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF
*
*     EXPSUM is the exponent range, approximately equal to
*     EMAX - EMIN + 1 .
*
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
*
*     NBITS is the total number of bits needed to store a
*     floating-point number.
*
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
*
*        Either there are an odd number of bits used to store a
*        floating-point number, which is unlikely, or some bits are
*        not used in the representation of numbers, which is possible,
*        (e.g. Cray machines) or the mantissa has an implicit bit,
*        (e.g. IEEE machines, DEC Vax machines), which is perhaps the
*        most likely. We have to assume the last alternative.
*        If this is true, then we need to reduce EMAX by one because
*        there must be some way of representing zero in an implicit-bit
*        system. On machines like Cray, we are reducing EMAX by one
*        unnecessarily.
*
         EMAX = EMAX - 1
      END IF
*
      IF( IEEE ) THEN
*
*        Assume we are on an IEEE machine which reserves one exponent
*        for infinity and NaN.
*
         EMAX = EMAX - 1
      END IF
*
*     Now create RMAX, the largest machine number, which should
*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
*
*     First compute 1.0 - BETA**(-P), being careful that the
*     result is less than 1.0 .
*
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      OLDY = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE )
     $      OLDY = Y
         Y = DLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE )
     $   Y = OLDY
*
*     Now multiply by BETA**EMAX to get RMAX.
*
      DO 30 I = 1, EMAX
         Y = DLAMC3( Y*BETA, ZERO )
   30 CONTINUE
*
      RMAX = Y
      RETURN
*
*     End of DLAMC5
*
      END
      

C----------------------------------------------------------------------------

C KS: instead of DLAMCH....
c This is the code that has gone into netlib as a replacement.
      DOUBLE PRECISION FUNCTION D1MACH(I)
      INTEGER I
C
C  DOUBLE-PRECISION MACHINE CONSTANTS
C
C  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C
C  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C
C  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C
C  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C
C  D1MACH( 5) = LOG10(B)
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.
C  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED.
C  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.)
C
C  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), ONE OF THE FIRST
C  TWO SETS OF CONSTANTS BELOW SHOULD BE APPROPRIATE.  IF YOU DO NOT
C  KNOW WHICH SET TO USE, TRY BOTH AND SEE WHICH GIVES PLAUSIBLE
C  VALUES.
C
C  WHERE POSSIBLE, DECIMAL, OCTAL OR HEXADECIMAL CONSTANTS ARE USED
C  TO SPECIFY THE CONSTANTS EXACTLY.  SOMETIMES THIS REQUIRES USING
C  EQUIVALENT INTEGER ARRAYS.  IF YOUR COMPILER USES HALF-WORD
C  INTEGERS BY DEFAULT (SOMETIMES CALLED INTEGER*2), YOU MAY NEED TO
C  CHANGE INTEGER TO INTEGER*4 OR OTHERWISE INSTRUCT YOUR COMPILER
C  TO USE FULL-WORD INTEGERS IN THE NEXT 5 DECLARATIONS.
C
C  COMMENTS JUST BEFORE THE END STATEMENT (LINES STARTING WITH *)
C  GIVE C SOURCE FOR D1MACH.
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
      INTEGER SC, CRAY1(38), J
      COMMON /D9MACH/ CRAY1
C/6S
C/7S
      COMMON/sizes/ SMALL, LARGE, RIGHT, DIVER, LOG10, SC
C/
      DOUBLE PRECISION DMACH(5)
C
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
      IF (SC .NE. 987) THEN
         DMACH(1) = 1.D13
         IF (      SMALL(1) .EQ. 1117925532
     *       .AND. SMALL(2) .EQ. -448790528) THEN
*           *** IEEE BIG ENDIAN ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2146435071
            LARGE(2) = -1
            RIGHT(1) = 1017118720
            RIGHT(2) = 0
            DIVER(1) = 1018167296
            DIVER(2) = 0
            LOG10(1) = 1070810131
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(2) .EQ. 1117925532
     *       .AND. SMALL(1) .EQ. -448790528) THEN
*           *** IEEE LITTLE ENDIAN ***
            SMALL(2) = 1048576
            SMALL(1) = 0
            LARGE(2) = 2146435071
            LARGE(1) = -1
            RIGHT(2) = 1017118720
            RIGHT(1) = 0
            DIVER(2) = 1018167296
            DIVER(1) = 0
            LOG10(2) = 1070810131
            LOG10(1) = 1352628735
         ELSE IF ( SMALL(1) .EQ. -2065213935
     *       .AND. SMALL(2) .EQ. 10752) THEN
*               *** VAX WITH D_FLOATING ***
            SMALL(1) = 128
            SMALL(2) = 0
            LARGE(1) = -32769
            LARGE(2) = -1
            RIGHT(1) = 9344
            RIGHT(2) = 0
            DIVER(1) = 9472
            DIVER(2) = 0
            LOG10(1) = 546979738
            LOG10(2) = -805796613
         ELSE IF ( SMALL(1) .EQ. 1267827943
     *       .AND. SMALL(2) .EQ. 704643072) THEN
*               *** IBM MAINFRAME ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2147483647
            LARGE(2) = -1
            RIGHT(1) = 856686592
            RIGHT(2) = 0
            DIVER(1) = 873463808
            DIVER(2) = 0
            LOG10(1) = 1091781651
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(1) .EQ. 1120022684
     *       .AND. SMALL(2) .EQ. -448790528) THEN
*           *** CONVEX C-1 ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2147483647
            LARGE(2) = -1
            RIGHT(1) = 1019215872
            RIGHT(2) = 0
            DIVER(1) = 1020264448
            DIVER(2) = 0
            LOG10(1) = 1072907283
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(1) .EQ. 815547074
     *       .AND. SMALL(2) .EQ. 58688) THEN
*           *** VAX G-FLOATING ***
            SMALL(1) = 16
            SMALL(2) = 0
            LARGE(1) = -32769
            LARGE(2) = -1
            RIGHT(1) = 15552
            RIGHT(2) = 0
            DIVER(1) = 15568
            DIVER(2) = 0
            LOG10(1) = 1142112243
            LOG10(2) = 2046775455
         ELSE
            DMACH(2) = 1.D27 + 1
            DMACH(3) = 1.D27
            LARGE(2) = LARGE(2) - RIGHT(2)
            IF (LARGE(2) .EQ. 64 .AND. SMALL(2) .EQ. 0) THEN
               CRAY1(1) = 67291416
               DO 10 J = 1, 20
                 CRAY1(J+1) = CRAY1(J) + CRAY1(J)
 10            CONTINUE
               CRAY1(22) = CRAY1(21) + 321322
               DO 20 J = 22, 37
                 CRAY1(J+1) = CRAY1(J) + CRAY1(J)
 20            CONTINUE
               IF (CRAY1(38) .EQ. SMALL(1)) THEN
*                  *** CRAY ***
*                 SMALL(1) = 2332160919536140288
                  SMALL(1) = 2332160
                  SMALL(1) = 1000000*SMALL(1) + 919536
                  SMALL(1) = 1000000*SMALL(1) + 140288
                  SMALL(2) = 0
*                 LARGE(1) = 6917247552664371199
                  LARGE(1) = 6917247
                  LARGE(1) = 1000000*LARGE(1) + 552664
                  LARGE(1) = 1000000*LARGE(1) + 371199
*                 LARGE(2) = 281474976710654
                  LARGE(2) = 28147497
                  LARGE(2) = 10000000*LARGE(2) + 6710654
*                 RIGHT(1) = 4585649583081652224
                  RIGHT(1) = 4585649
                  RIGHT(1) = 1000000*RIGHT(1) + 583081
                  RIGHT(1) = 1000000*RIGHT(1) + 652224
                  RIGHT(2) = 0
*                 DIVER(1) = 4585931058058362880
                  DIVER(1) = 4585931
                  DIVER(1) = 1000000*DIVER(1) + 058058
                  DIVER(1) = 1000000*DIVER(1) + 362880
                  DIVER(2) = 0
*                 LOG10(1) = 4611574008272714703
                  LOG10(1) = 4611574
                  LOG10(1) = 1000000*LOG10(1) +   8272
                  LOG10(1) = 1000000*LOG10(1) + 714703
*                 LOG10(2) = 272234615232940
                  LOG10(2) = 27223461
                  LOG10(2) = 10000000*LOG10(2) + 5232940
               ELSE
              CALL rexit('Stopped in D1mach')
C                  STOP 779
                  END IF
            ELSE
              CALL rexit('Stopped in D1mach')
C               STOP 779
               END IF
            END IF
         SC = 987
         END IF
      D1MACH = DMACH(I)
      RETURN
C/
C
      END
