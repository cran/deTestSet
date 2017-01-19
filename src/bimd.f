C karline: all write statements -> rprint and removed / from formats
C karline: removed argument ierr from call to func, jac and mas

C -----------------------------------------------------------------------------------
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


      SUBROUTINE BIMD(M,FCN,T0,TEND,Y0,H,
     &                RTOL,ATOL,ITOL,
     &                JAC,IJAC,MLJAC,MUJAC,
     &                MAS,IMAS,MLMAS,MUMAS,
     &                SOLOUT,IOUT,
     &                WORK,LWORK,IWORK,LIWORK,
     &                RPAR,IPAR,IDID)

C -----------------------------------------------------------------------------------
C -----------------------------------------------------------------------------------
C
C     PURPOSE:    BIMD SOLVES A (STIFF) DIFFERENTIAL ODE PROBLEM,
C     --------
C                               Y'    = F(T,Y),     T0<=T<=TEND,
C                               Y(T0) = Y0,
C
C                 OR A LINEARLY IMPLICIT DAE PROBLEM OF INDEX UP TO 3 WITH
C                 CONSTANT MASS MATRIX, NAMELY PROBLEM IN THE FORM
C
C                 M Y'  = F(T,Y),     T0<=T<=TEND,
C                 Y(T0) = Y0,
C
C                 WHERE M IS A POSSIBLY SINGULAR MATRIX.
C
C                 THE CODE IS BASED ON BLENDED IMPLICIT METHODS.
C                 BLENDED IMPLICIT METHODS ARE A CLASS OF BLOCK
C                 METHODS       PROVIDING A (RELATIVELY) EASY DEFINITION
C                 OF SUITABLE NONLINEAR SPLITTINGS FOR SOLVING THE
C                 CORRESPONDING DISCRETE PROBLEMS [1,5-7].
C                 THE CODE BIMD IMPLEMENTS A VARIABLE STEPSIZE-
C                 VARIABLE ORDER METHOD. ORDERS: 4-6-8-10-12.
C                 IMPLEMENTATION DETAILS ARE IN REFERENCES [1-5].
C
C
C
C     AUTHORS:    C.MAGHERINI,
C     --------    DIPARTIMENTO DI MATEMATICA APPLICATA "U.DINI"
C                 VIA BUONARROTI, 1/C
C                 56127 PISA
C                 ITALY
C                 E-MAIL: CECILIA.MAGHERINI@ING.UNIPI.IT
C
C                 L.BRUGNANO,
C                 DIPARTIMENTO DI MATEMATICA "U.DINI"
C                 VIALE MORGAGNI 67/A
C                 50134 FIRENZE
C                 ITALY
C                 E-MAIL: BRUGNANO@MATH.UNIFI.IT
C
C     CODE HOME PAGE:   http://www.math.unifi.it/~brugnano/BiM/index.html
C     ---------------
C
C     CODE:       THE CODE IS MADE UP OF TWO FILES:
C     -----        - BIMD.F     (I.E. THE PRESENT FILE) WHICH CONTAINS THE MAIN
C                    INTEGRATION PROCEDURE
C                  - BIMDA.F    CONTAINING ADDITIONAL AND LINEAR ALGEBRA
C                    PROCEDURES
C
C     CURRENT RELEASE:   1.1.1,  September, 2006.
C     ----------------
C
C     RELEASE HISTORY:   1.0,  October, 2005
C     ---------------    - first version released;
C
C                        1.1,  July, 2006
C                        Main features (wrt 1.0):
C                        - improved definition of the
C                          coefficients of the methods
C                        - the results described
C                          in reference [1], have been better
C                          exploited for the definition of
C                          the stopping criterion for
C                          the splitting Newton blended iteration
C                        - improved choice of the initial profile
C                          after a failure due to Newton convergence
C                        - possibility of solving the problem
C                          with vector-valued absolute input
C                          tolerances
C                        - new function CONTSOL, to be used when
C                          continuous output is desired
C                        - minor changes concerning the
C                          order variation strategy.
C
C                        1.1.1, September, 2006
C                        - some minor bugs fixed
C
C
C
C     REFERENCES:
C     -----------
C                 [1] L.BRUGNANO, C.MAGHERINI, F.MUGNAI.
C                     Blended Implicit Methods for the Numerical Solution
C                     of DAE problems.
C                     Jour. CAM 189 (2006) 34-50.
C
C                 [2] L.BRUGNANO, C.MAGHERINI
C                     The BiM code for the numerical solution of ODEs
C                     Jour. CAM 164-165 (2004) 145-158.
C
C                 [3] L.BRUGNANO, C.MAGHERINI
C                     Some Linear Algebra issues concerning the implementation
C                     of Blended Implicit Methods
C                     Numer. Linear Alg. Appl. 12 (2005) 305-314.
C
C                 [4] L.BRUGNANO, C.MAGHERINI
C                     Economical Error Estimates for Block Implicit Methods for
C                     ODEs via Deferred Correction.
C                     Appl. Numer. Math. 56 (2006) 608-617.
C
C                 [5] L.BRUGNANO, C.MAGHERINI
C                     Blended Implementation of Block Implicit Methods for ODEs
C                     Appl. Numer. Math. 42 (2002) 29-45.
C
C                 [6] L.BRUGNANO, D.TRIGIANTE
C                     Block Implicit Methods for ODEs
C                     in "Recent Trends in Numerical Analysis", D.Trigiante Ed.
C                     Nova Science Publ. Inc., New York, 2001, pp. 81-105.
C
C                 [7] L.BRUGNANO
C                     Blended Block BVMs (B$_3$VMs): a Family of economical
C                     implicit methods for ODEs
C                     Jour. CAM 116 (2000) 41-62.
C
C
C
C    REMARK:   The code BiMD has been written using a style very similar to the one
C    -------   used in the codes RADAU and GAM. Indeed, some subroutines and comments
C              have been imported from such codes. Moreover, the name and the meaning
C              of a number of input parameters and variables have been fully inherited
C              from them.
C
C -----------------------------------------------------------------------------------
C -----------------------------------------------------------------------------------
C
C     USAGE:
C     ------
C
C      CALL BIMD(M,FCN,T0,TEND,Y0,H,
C     &          RTOL,ATOL,ITOL,
C     &          JAC,IJAC,MLJAC,MUJAC,
C     &          MAS,IMAS,MLMAS,MUMAS,
C     &          SOLOUT,IOUT,
C     &          WORK,LWORK,IWORK,LIWORK,
C     &          RPAR,IPAR,IDID)
C
C     NOTE:   IN ORDER TO GAIN THE BEST PERFORMANCE, THE EXECUTABLE HAS TO
C     -----   BE CREATED WITH THE OPTION ALLOWING  TO CONTINUE THE EXECUTION
C             AFTER A FLOATING-POINT EXCEPTION  (E.G., BY USING THE OPTION
C             -FPE; SEE YOUR FORTRAN COMPILER REFERENCE MANUAL).
C                       THE ISNAN LOGICAL FUNCTION IS REQUIRED, TO RECOGNIZE NANs. IF
C                       NOT SUPPORTED BY YOUR COMPILER, A STANDARD ONE IS PROVIDED AT
C                       THE TOP OF THE SUBBIM.F FILE.  Karline: ISNANUM
C
C
C -----------------------------------------------------------------------------------
C           INPUT PARAMETERS
C -----------------------------------------------------------------------------------
C
C M             SIZE OF THE PROBLEM
C
C FCN           SUBROUTINE WITH THE FUNCTION F(T,Y) TO BE INTEGRATED.
C karline: changed this to fcn(m,t,y,dy,rpar,ipar)
C
C      subroutine fcn(m,t,y,dy,ierr,rpar,ipar)
C      double precision t,y,dy,rpar(*)
C      integer m,ierr,ipar(*)
C      dimension y(m),dy(m)
CC     m      size of the continuous problem
CC     t,y    is the point where f is evaluated
CC     dy     will contain the value of f(t,y)
CC     ierr   is a return code (0 means OK)
CC     rpar   possible external real parameters
CC     ipar   possible external integer parameters
C      ................
C      return
C      end
C
C T0-TEND       INTEGRATION INTERVAL
C
C Y0            INITIAL CONDITION
C
C H             INITIAL STEPSIZE
C
C RTOL-ATOL     RELATIVE AND ABSOLUTE TOLERANCES.
C               ATOL CAN BE EITHER SCALAR OR A
C               VECTOR OF LENGTH M.
C
C ITOL          SWITCH FOR ATOL:
C
C               ITOL = 0 --> ATOL IS SCALAR.
C                            THE CODE PROVIDES A NUMERICAL SOLUTION
C                            WITH THE LOCAL ERROR OF Y(I) ROUGHLY SMALLER
C                            THAN ATOL + RTOL*ABS(Y(I))
C
C               ITOL = 1 --> ATOL IS AN ARRAY OF LENGTH M
C                            THE LOCAL ERROR OF Y(I) IS KEPT
C                            BELOW ATOL(I) + RTOL*ABS(Y(I))
C
C JAC           SUBROUTINE EVALUATING THE JACOBIAN OF F (DUMMY, IF IJAC=0)
C karline: removed ierr from call
C francesca: added mu, ml
C      subroutine jac(m,t,y,ml,mu,jac,ldjac,rpar,ipar)
C      double precision t,y,jac,rpar(*)
C      integer m,ldjac,ierr,ipar(*),ml,mu
C      dimension y(m),jac(ldjac,m)
CC     m      size of the continuous problem
CC     t,y      is the point where the Jacobian is evaluated
CC     jac      will contain the value of the Jacobian at (t,y)
CC     ldjac  leading dimension of the array jac
CC     ierr     is a return code (0 means OK)
CC     rpar   possible external real parameters
CC     ipar     possible external integer parameters
C      ............
C      return
C      end
C
C IJAC          FLAG: 0=NUMERICAL JACOBIAN, ANALYTICAL OTHERWISE
C
C MLJAC-MUJAC   LOWER-UPPER BANDWIDTH OF THE JACOBIAN (MLJAC=M IF FULL JACOBIAN)
C
C MAS           SUBROUTINE EVALUATING THE MASS-MATRIX (DUMMY, IF IMAS=0)
C karline: remove ierr argument
C
C      subroutine mas(m,mas,ldmas,rpar,ipar)
C      double precision mas,rpar(*)
C      integer m,ldmas,ierr,ipar(*)
C      dimension mas(ldmas,m)
CC     m      size of the continuous problem
CC     mas    will contain the evaluated mass-matrix
CC     ldmas  leading dimension of the array mas
CC     ierr     is a return code (0 means OK)
CC     rpar   possible external real parameters
CC     ipar     possible external integer parameters
C      ..............
C      return
C      end
C
C
C IMAS            FLAG: 0=ODE, DAE OTHERWISE
C
C MLMAS-MUMAS     LOWER-UPPER BANDWIDTH OF THE MASS-MATRIX (MLMAS=M IF FULL MASS-MATRIX)
C                 MLMAS IS SUPPOSED TO BE .LE. MLJAC
C                 MUMAS IS SUPPOSED TO BE .LE. MUJAC.
C
C LWORK   LENGTH OF WORK   ( LWORK >= 14 +KMAX +9*M +5*KMAX*M +M*(LDJAC+LDLU+LDMAS),
C
C         WHERE:
C
C                            LDJAC = LDLU = M,            IN CASE OF A FULL JACOBIAN,
C            LDJAC = MLJAC+MUJAC+1,  LDLU = LDJAC+MLJAC,  IN CASE OF A BANDED JACOBIAN;
C
C                            LDMAS = M                    IN CASE OF A FULL MASS MATRIX,
C                            LDMAS = MLMAS+MUMAS+1        IN CASE OF A BANDED MASS MATRIX,
C                            LDMAS = 1                    IN THE ODE CASE (I.E. IMAS = 0)
C
C                             KMAX = ORDMAX-2,            IF ORDMAX>4,
C                                    3,                   IF ORDMAX=4. )
C
C WORK(1)   UROUND. MACHINE PRECISION. (DEFAULT = 1.D-16)
C
C WORK(2)   HMAX. MAXIMUM INTEGRATION STEP. (DEFAULT = (TEND-T0)/8)
C
C WORK(3)   FACTOL - SAFETY FACTOR FOR THE STOPPING CRITERION OF THE BLENDLED ITERATION.
C           METHOD OF ORDER 4. (DEFAULT = 1.D-1)
C
C WORK(4)   FACTOL - SAFETY FACTOR FOR THE STOPPING CRITERION OF THE BLENDLED ITERATION.
C           METHOD OF ORDER 6. (DEFAULT = 1.D-1)
C
C WORK(5)   FACTOL - SAFETY FACTOR FOR THE STOPPING CRITERION OF THE BLENDLED ITERATION.
C           METHOD OF ORDER 8. (DEFAULT = 1.D-1)
C
C WORK(6)   FACTOL - SAFETY FACTOR FOR THE STOPPING CRITERION OF THE BLENDLED ITERATION.
C           METHOD OF ORDER 10. (DEFAULT = 1.D-1)
C
C WORK(7)   FACTOL - SAFETY FACTOR FOR THE STOPPING CRITERION OF THE BLENDLED ITERATION.
C           METHOD OF ORDER 12. (DEFAULT = 1.D-1)
C
C WORK(8)   FACTOL - SAFETY FACTOR FOR THE STOPPING CRITERION OF THE BLENDLED ITERATION
C           IN CASE OF SMALL VALUES OF min(abs(y_0)), min(abs(f_0)) AND OF max(abs(f_0)).
C           (DEFAULT = 1D-2)
C
C WORK(9)   FACTOL - SAFETY FACTOR FOR THE STOPPING CRITERION OF THE BLENDLED ITERATION
C           IN CASE OF SLOWLY VARYING. (DEFAULT = 5D-2)
C
C WORK(10)-WORK(11)  FACL-FACR. THE NEW STEPSIZE MUST SATISFY FACL<=HNEW/HOLD<= FACR.
C           (DEFAULT: WORK(10)=1.2D-1, WORK(11)=1D1)
C
C WORK(12)  SFTY - SAFETY FACTOR FOR PREDICTING THE NEW STEPSIZE FOR THE CURRENT ORDER
C           METHOD. (DEFAULT = 1D0/2D1)
C
C WORK(13)  SFTYUP - SAFETY FACTOR FOR PREDICTING THE NEW STEPSIZE FOR THE HIGHER ORDER
C           METHOD. (DEFAULT = SFTY/2D0)
C
C WORK(14)  SFTYDN - SAFETY FACTOR FOR PREDICTING THE NEW STEPSIZE FOR THE LOWER ORDER
C           METHOD. (DEFAULT = SFTY)
C
C LIWORK    LENGTH OF IWORK  (LIWORK >= M+40)
C
C IWORK( 1) MAX NUMBER OF INTEGRATION STEPS (DEFAULT = 100000).
C
C IWORK( 2) ORDMIN, 4<=ORDMIN<=12. (DEFAULT = 4).
C
C IWORK( 3) ORDMAX, ORDMIN<=ORDMAX<=12. (DEFAULT = 12).
C
C IWORK( 4) MAX NUMBER OF BLENDED ITERATIONS PER INTEGRATION STEP, METHOD OF ORDER 4
C           (DEFAULT = 10).
C
C IWORK( 5) MAX NUMBER OF BLENDED ITERATIONS PER INTEGRATION STEP, METHOD OF ORDER 6
C           (DEFAULT = 12).
C
C IWORK( 6) MAX NUMBER OF BLENDED ITERATIONS PER INTEGRATION STEP, METHOD OF ORDER 8
C           (DEFAULT = 14).
C
C IWORK( 7) MAX NUMBER OF BLENDED ITERATIONS PER INTEGRATION STEP, METHOD OF ORDER 10
C           (DEFAULT = 16).
C
C IWORK( 8) MAX NUMBER OF BLENDED ITERATIONS PER INTEGRATION STEP, METHOD OF ORDER 12
C           (DEFAULT = 18).
C
C IWORK( 9) DIMENSION OF THE INDEX 1 VARIABLES (DEFAULT = M).
C           IT MUST BE GREATER THAN 0.
C
C IWORK(10) DIMENSION OF THE INDEX 2 VARIABLES (DEFAULT = 0).
C
C IWORK(11) DIMENSION OF THE INDEX 3 VARIABLES (DEFAULT = 0).
C
C REMARK: THE VARIABLES MUST BE SORTED BY INCREASING INDEX.
C -------
C
C RPAR,IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH
C            CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING
C            PROGRAM AND THE FCN, JAC, MAS AND SOLOUT SUBROUTINES.
C
C SOLOUT     NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE NUMERICAL
C            SOLUTION DURING INTEGRATION.
C            IF IOUT = 1, IT IS CALLED AFTER EACH SUCCESSFUL STEP.
C            SUPPLY A DUMMY SUBROUTINE IF IOUT = 0.
C            IT MUST HAVE THE FOLLOWING FORM:
C
C      subroutine solout(m,k,ord,t0,t,y,f,dd,rpar,ipar,irtrn)
C      integer m,k,ord,irtrn,ipar(*)
C      double precision t0,t,y,f,dd, rpar(*)
C      dimension t(k),y(m,k),f(m,k),dd(k+1,m)
C
CC     m                is the size of the problem
CC     k                is the block-size of the method
CC     ord              is the order of the method
CC     t0               is the starting time point of the step
CC     t                contains the (internal) mesh points of
CC                      the step
CC     y                is the current numerical solution
CC     f                contains the values of fcn(t,y)
CC     dd               contains the divided differences of y
CC                      over the internal mesh points of the step
CC                      (to be used, for example, if continuous
CC                      output is desired, see below)
CC     rpar             possible external real parameters
CC     ipar             possible external integer parameters
CC     irtrn            is a return code. If set <0, BiMD returns
CC                      to the calling program.
C
C      ................
C      return
C      end
C
C
C           CONTINUOUS OUTPUT:
C           ------------------
C
C           DURING CALLS TO SOLOUT, A CONTINUOUS SOLUTION
C           FOR THE INTERVAL [t0,t(k)] IS AVAILABLE THROUGH
C           THE FUNCTION
C
C               CONTSOL(I,T,M,K,T0,TSTEP,DD)
C
C           WHICH PROVIDES AN APPROXIMATION TO THE I-TH
C           COMPONENT OF THE SOLUTION AT THE TIME POINT T.
C
C
C IOUT      SWITCH FOR CALLING THE SUBROUTINE SOLOUT.
C
C           IOUT = 0, SOLOUT IS NEVER CALLED
C           IOUT = 1, SOLOUT IS CALLED AFTER EACH
C                     SUCCESSFULL STEP
C
C
C
C -----------------------------------------------------------------------------------
C           OUTPUT PARAMETERS
C -----------------------------------------------------------------------------------
C
C
C T0        VALUE OF T UP TO WHERE THE SOLUTION HAS BEEN COMPUTED
C           (IF THE INTEGRATION HAS BEEN SUCCESFULL,THEN T0=TEND)
C
C Y0        NUMERICAL SOLUTION AT T0
C
C IDID      RETURN CODE:
C              0  SUCCESFULL RUN
C             -1  WRONG INPUT PARAMETERS
C             -2  LARGER NMAX IS NEEDED
C             -3  STEPSIZE TOO SMALL
C             -4  REPEATEDLY SINGULAR MATRIX
C             -5  TOO MANY CONSECUTIVE NEWTON FAILURES
C             -6  ERROR CODE RETURNED BY THE JAC SUBROUTINE OR BY THE FCN SUBROUTINE
C                 AT THE STARTING POINT
C
C IWORK(12) NUMBER OF FUNCTION EVALUATIONS
C
C IWORK(13) NUMBER OF JACOBIAN EVALUATIONS
C
C IWORK(14) NUMBER OF LU DECOMPOSITION
C
C IWORK(15) NUMBER OF LINEAR SYSTEMS SOLVED
C
C IWORK(16)-IWORK(20) NUMBER OF BLENDED ITERATIONS PER METHOD
C
C IWORK(21)-IWORK(25) NUMBER OF STEP PER METHOD
C
C IWORK(26)-IWORK(30) NUMBER OF ACCEPTED STEP PER METHOD
C
C IWORK(31)-IWORK(35) NUMBER OF REFUSED STEP PER METHOD (ERROR TEST)
C
C IWORK(36)-IWORK(40) NUMBER OF REFUSED STEP PER METHOD (NEWTON'S CONVERGENCE)
C
C -----------------------------------------------------------------------------------
C -----------------------------------------------------------------------------------
Cc-----------------------------------------------------------------------
Cc     Sample driver for the code BIMD
Cc-----------------------------------------------------------------------
C      program chemakzod
C      implicit none
C      integer MMAX,lwork,liwork
C      parameter(MMAX=6,lwork=24+MMAX*(70+2*MMAX),liwork=MMAX+40)
C      double precision y(MMAX),work(lwork)
C      integer iwork(lwork),ijac,imas,iout
C      external feval, jeval,  solout, meval
C      character problm*8
C      double precision t0,tf,h0,h,rtol,atol,rpar(1)
C      integer neqn,mljac,mujac,mlmas,mumas,itol,i,ierr, ipar(1)
C      integer NSTEPS,NACCEPT,NFAILERR,NFAILNEWT,NITER,idid
C
C      double precision ks
C      parameter (ks   =115.83d0)
C
C      neqn  = 6
C      t0    = 0d0
C      tf    = 180d0
C      ijac  = 1
C      mljac = neqn
C      mujac = neqn
C      mlmas = 0
C      mumas = 0
C
C
C      y(1) = 0.444d0
C      y(2) = 0.00123d0
C      y(3) = 0d0
C      y(4) = 0.007d0
C      y(5) = 0d0
C      y(6) = ks*y(1)*y(4)
C
Cc-----------------------------------------------------------------------
Cc     read the tolerances and initial stepsize
Cc-----------------------------------------------------------------------
C      write(6,*) 'give the absolute tolerance'
C      read(5,*)   atol
C      write(6,*) 'give the relative tolerance'
C      read(5,*)   rtol
C      write(6,*) 'give the initial  stepsize '
C      read(5,*)   h0
C
C      h = h0
C      do i=1,8
C         iwork(i) = 0
C      end do
C      iwork(9)  = neqn
C      iwork(10) = 0
C      iwork(11) = 0
C      do i=1,14
C         work(i) = 0d0
C      end do
C
C      iout = 0
C      idid = 0
C      imas = 1
C      itol = 0
C
Cc-----------------------------------------------------------------------
Cc     call of the subroutine BIMD
Cc-----------------------------------------------------------------------
C      call BIMD(neqn,feval,t0,tf,y,h,rtol,atol,itol,
C     &          jeval,ijac,mljac,mujac,
C     &          meval,imas,mlmas,mumas,
C     &          solout,iout,
C     &          work,lwork,iwork,liwork,
C     &          rpar,ipar,idid)
C      if (idid.ne.0) then
C         write(6,*) 'ERROR: returned idid =', idid
C         goto 20
C      endif
Cc-----------------------------------------------------------------------
Cc     print final solution
Cc-----------------------------------------------------------------------
C      write(6,10)
C   10 format(//)
C
C      write(6,11) atol,rtol,h0
C   11 format(/,' we solved the problem with',//,
C     +       '       absolute tolerance = ',d10.4,',',/,
C     +       '       relative tolerance = ',d10.4,',',/,
C     +       '     and initial stepsize = ',d10.4,//)
Cc-----------------------------------------------------------------------
Cc     print error with respect to reference solution
Cc-----------------------------------------------------------------------
C      NSTEPS    = 0
C      NACCEPT   = 0
C      NFAILERR  = 0
C      NFAILNEWT = 0
C      DO I=1,5
C       NSTEPS    = NSTEPS    + iwork(I+19)
C       NACCEPT   = NACCEPT   + iwork(I+25)
C       NFAILERR  = NFAILERR  + iwork(I+30)
C       NFAILNEWT = NFAILNEWT + iwork(I+35)
C      END DO
C
C     write(6,41) NSTEPS,NACCEPT,NFAILNEWT,NFAILERR,
C    &            IWORK(12),IWORK(13),IWORK(14)
C
C  41 format(
C    +         ' # Steps              = ',i8,/
C    +         ' # Accept             = ',i8,/,
C    +         ' # Failnwt            = ',i8,/,
C    +         ' # Failerr            = ',i8,/,
C    +         ' # F-eval             = ',i8,/,
C    +         ' # Jac-eval           = ',i8,/,
C    +         ' # LU-decomp          = ',i8,/,/)
C
C     write(6,42) (y(i),i=1,6)
C  42 format(
C    &        'Numerical solution:',//
C    &        'y(1) = ',e22.15,/,
C    &        'y(2) = ',e22.15,/,
C    &        'y(3) = ',e22.15,/,
C    &        'y(4) = ',e22.15,/,
C    &        'y(5) = ',e22.15,/,
C    &        'y(6) = ',e22.15,/,/)
C
C20    continue
C     end
C
Cc-----------------------------------------------------------------------
Cc     AUXILIARY ROUTINES
C
C      subroutine feval(neqn,t,y,f,ierr,rpar,ipar)
C      integer neqn,ierr,ipar(*)
C      double precision t,y(neqn),f(neqn),rpar(*)
C
C      double precision k1,k2,k3,k4,kbig,kla,po2,hen,ks
C      parameter (
C     +   k1   = 18.7d0,
C     +   k2   = 0.58d0,
C     +   k3   = 0.09d0,
C     +   k4   = 0.42d0,
C     +   kbig = 34.4d0,
C     +   kla  = 3.3d0,
C     +   ks   = 115.83d0,
C     +   po2  = 0.9d0,
C     +   hen  = 737d0
C     +)
C      double precision r1,r2,r3,r4,r5,fin
C
C      if (y(2) .lt. 0d0) then
C         ierr = -1
C         return
C      endif
C
C      r1  = k1*(y(1)**4)*sqrt(y(2))
C      r2  = k2*y(3)*y(4)
C      r3  = k2/kbig*y(1)*y(5)
C      r4  = k3*y(1)*(y(4)**2)
C      r5  = k4*(y(6)**2)*sqrt(y(2))
C      fin = kla*(po2/hen-y(2))
C
C      f(1) =   -2d0*r1 +r2 -r3     -r4
C      f(2) = -0.5d0*r1             -r4     -0.5d0*r5 + fin
C      f(3) =        r1 -r2 +r3
C      f(4) =           -r2 +r3 -2d0*r4
C      f(5) =            r2 -r3         +r5
C      f(6) = ks*y(1)*y(4)-y(6)
C
C      return
C      end
Cc-----------------------------------------------------------------------
C      subroutine jeval(neqn,t,y,jac,ldjac,ierr,rpar,ipar)
C      integer ldjac,neqn,ierr,ipar(*)
C      double precision t,y(neqn),jac(ldjac,neqn),rpar(*)
C
C      integer i,j
C
C      double precision k1,k2,k3,k4,kbig,kla,ks
C      parameter (
C     +   k1   =18.7d0,
C     +   k2   =0.58d0,
C     +   k3   =0.09d0,
C     +   k4   =0.42d0,
C     +   kbig =34.4d0,
C     +   kla  =3.3d0,
C     +   ks   =115.83d0
C     +)
C      double precision r11,r12,r23,r24,r31,r35,r41,r44,r52,r56,fin2
C
C      if (y(2) .lt. 0d0) then
C         ierr = -1
C         return
C      endif
C
C      do 20 j=1,neqn
C         do 10 i=1,neqn
C            jac(i,j) = 0d0
C   10    continue
C   20 continue
C
C      r11  = 4d0*k1*(y(1)**3)*sqrt(y(2))
C      r12  = 0.5d0*k1*(y(1)**4)/sqrt(y(2))
C      r23  = k2*y(4)
C      r24  = k2*y(3)
C      r31  = (k2/kbig)*y(5)
C      r35  = (k2/kbig)*y(1)
C      r41  = k3*y(4)**2
C      r44  = 2d0*k3*y(1)*y(4)
C      r52  = 0.5d0*k4*(y(6)**2)/sqrt(y(2))
C      r56  = 2d0*k4*y(6)*sqrt(y(2))
C      fin2 = -kla
C
C      jac(1,1) = -2d0*r11-r31-r41
C      jac(1,2) = -2d0*r12
C      jac(1,3) = r23
C      jac(1,4) = r24-r44
C      jac(1,5) = -r35
C      jac(2,1) = -0.5d0*r11-r41
C      jac(2,2) = -0.5d0*r12-0.5d0*r52+fin2
C      jac(2,4) = -r44
C      jac(2,6) = -0.5d0*r56
C      jac(3,1) = r11+r31
C      jac(3,2) = r12
C      jac(3,3) = -r23
C      jac(3,4) = -r24
C      jac(3,5) = r35
C      jac(4,1) = r31-2d0*r41
C      jac(4,3) = -r23
C      jac(4,4) = -r24-2d0*r44
C      jac(4,5) = r35
C      jac(5,1) = -r31
C      jac(5,2) = r52
C      jac(5,3) = r23
C      jac(5,4) = r24
C      jac(5,5) = -r35
C      jac(5,6) = r56
C      jac(6,1) = ks*y(4)
C      jac(6,4) = ks*y(1)
C      jac(6,6) = -1d0
C
C      return
C      end
Cc-----------------------------------------------------------------------
C      subroutine meval(neqn,mas,ldmas,ierr,rpar,ipar)
C      integer ldmas,neqn,ierr,ipar(*)
C      double precision t,y(neqn),yprime(neqn),mas(ldmas,neqn),rpar(*)
C
C      integer i
C
C      do 10 i=1,neqn-1
C         mas(1,i)=1d0
C   10 continue
C
C      mas(1,neqn)=0d0
C
C      return
C      end
C
Cc-----------------------------------------------------------------------
C      subroutine solout(m,k,ord,t0,t,y,f,dd,rpar,ipar,irtrn)
C      implicit none
C      integer m,k,ord,irtrn,ipar(*)
C      double precision t0,t(k),y(m,k),f(m,k),dd(k+1,m),rpar(*)
C
Cc     dummy subroutine
C
C      return
C      end
C
C -----------------------------------------------------------------------------------
C -----------------------------------------------------------------------------------

      IMPLICIT NONE

      EXTERNAL FCN,JAC,MAS,SOLOUT
      INTEGER M,LWORK,LIWORK,IWORK(LIWORK),
     &        ITOL,IJAC,MLJAC,MUJAC,LDJAC,
     &        IMAS,MLMAS,MUMAS,LDMAS,IOUT,IDID,
     &        LDLU,IJOB(2),IPAR(*)
      LOGICAL JBAND,MBAND
      DOUBLE PRECISION T0,TEND,Y0(M),H,RTOL,ATOL(*),WORK(LWORK),RPAR(*)

      INTEGER NMETH,KMAX
      PARAMETER (NMETH=5,KMAX=10)

      INTEGER MAXSTEP,ORDMIN,ORDMAX,ITMAX(NMETH),STEP_ORD(NMETH)

      DOUBLE PRECISION UROUND, FACNEWTV(NMETH),FACNSMALL,
     &                 FACNRESTR,FACL,FACR,
     &                 SFTY, SFTYUP, SFTYDN, HMAX,
     &                 RHOMUV(NMETH),RHOMLV(NMETH),TOLESTRAPR

      INTEGER I,INDF0,INDT,INDIPVT,INDEJ0,INDY,INDF,INDTHETA,INDJ0,
     &        INDERR,INDSCAL,INDDELJ0,INDDELJ0OLD,INDTOLEXT,
     &        INDFJ0,INDORD,INDSCALEXT,IND_DD,INDM0,INDTEMP

      LOGICAL STOPINT

      STEP_ORD(1)= 3
      STEP_ORD(2)= 4
      STEP_ORD(3)= 6
      STEP_ORD(4)= 8
      STEP_ORD(5)= 10

      STOPINT = .FALSE.

C     INITIAL STEP-SIZE
      IF (H.EQ.0D0) THEN
          H=1.D-6
      ELSEIF(H.LT.0D0) THEN
        CALL Rprintd1('Wrong input H = ',H)
         STOPINT=.TRUE.
      END IF

C--------------------------------------------------
C     PARAMETERS INITIALIZATION
C--------------------------------------------------

      IF (IWORK(1).EQ.0) THEN
         MAXSTEP=100000
      ELSE
         MAXSTEP=IWORK(1)
         IF (MAXSTEP.LE.0) THEN
        CALL Rprinti1('Wrong input iwork(1) = ',IWORK(1))
             STOPINT=.TRUE.
         END IF
      ENDIF

      IF (IWORK(2).EQ.0) THEN
         ORDMIN = 4
      ELSE
         ORDMIN = IWORK(2)
         indord = ORDMIN/2-1
         IF ((indord.LE.0).OR.(indord.GT.NMETH)) THEN
        CALL Rprinti1('Wrong input iwork(2) = ',IWORK(2))
             STOPINT=.TRUE.
         END IF
         ORDMIN = 2*(indord+1)
      ENDIF

      IF (IWORK(3).EQ.0) THEN
         ORDMAX = 12
         indord = NMETH
      ELSE
         ORDMAX = IWORK(3)
         indord = ORDMAX/2 - 1
         IF ((indord.LE.0).OR.(indord.GT.NMETH)) THEN
        CALL Rprinti1('Wrong input iwork(3) = ',IWORK(3))
             STOPINT=.TRUE.
         END IF
         ORDMAX = 2*(indord+1)
      ENDIF

      IF (ORDMIN.GT.ORDMAX) THEN
        CALL Rprinti1('Invalid values for ORDMIN= ',IWORK(2))
        CALL Rprinti1('Invalid values for ORDMAX= ',IWORK(3))
        STOPINT=.TRUE.
      END IF

      IF (IWORK(4) .EQ. 0) THEN
         ITMAX(1) = 10
      ELSE
         ITMAX(1) = IWORK(4)
         IF (ITMAX(1).LE.0) THEN
          CALL Rprinti1('Wrong input iwork(4) = ',IWORK(4))
             STOPINT=.TRUE.
         END IF
      END IF

      IF (IWORK(5) .EQ. 0) THEN
         ITMAX(2) = 12
      ELSE
         ITMAX(2) = IWORK(5)
         IF (ITMAX(2).LE.0) THEN
        CALL Rprinti1('Wrong input iwork(5) = ',IWORK(5))
             STOPINT=.TRUE.
         END IF
      END IF

      IF (IWORK(6) .EQ. 0) THEN
         ITMAX(3) = 14
      ELSE
         ITMAX(3) = IWORK(6)
         IF (ITMAX(3).LE.0) THEN
        CALL Rprinti1('Wrong input iwork(6) = ',IWORK(6))
             STOPINT=.TRUE.
         END IF
      END IF

      IF (IWORK(7) .EQ. 0) THEN
         ITMAX(4) = 16
      ELSE
         ITMAX(4) = IWORK(7)
         IF (ITMAX(4).LE.0) THEN
        CALL Rprinti1('Wrong input iwork(7) = ',IWORK(7))
             STOPINT=.TRUE.
         END IF
      END IF

      IF (IWORK(8) .EQ. 0) THEN
         ITMAX(5) = 18
      ELSE
         ITMAX(5) = IWORK(8)
         IF (ITMAX(5).LE.0) THEN
        CALL Rprinti1('Wrong input iwork(8) = ',IWORK(8))
             STOPINT=.TRUE.
         END IF
      END IF

        IF ((IWORK(9)+IWORK(10)+IWORK(11)).EQ.0) THEN
                IWORK( 9) = M
                IWORK(10) = 0
                IWORK(11) = 0
        ELSE
          IF (IWORK(9).EQ.0) THEN
        CALL Rprinti1('Invalid value iwork(9), should be >0 ',iwork(9))
             STOPINT = .TRUE.
          END IF
          IF ((IWORK(9)+IWORK(10)+IWORK(11)).NE.M) THEN
           CALL Rprinti3('invalid values for iwork(9:11)', 
     &     IWORK(9), IWORK(10), IWORK(11))
              STOPINT=.TRUE.
          END IF
      END IF

      IF ((IMAS.NE.0).AND.((MLMAS.GT.MLJAC).OR.(MUMAS.GT.MUJAC)))
     &THEN
        CALL Rprint(
     &   'Bandwidth of MAS not smaller than bandwidth of JAC')
          STOPINT=.TRUE.
      END IF

      IF (WORK(1) .EQ. 0D0) THEN
         UROUND = 1.0D-16
      ELSE
         UROUND = WORK(1)
         IF ((UROUND.LE.0D0).OR.(UROUND.GE.1D0)) THEN
        CALL Rprintd1('Wrong input work(1) = ',WORK(1))
             STOPINT=.TRUE.
         END IF
      END IF

      IF (RTOL.LE.UROUND) THEN
        CALL Rprint('rtol is too small')
          STOPINT = .TRUE.
      END IF

      IF (ITOL.EQ.0) THEN
        IF (ATOL(1).LE.0D0) THEN
        CALL Rprint('atol is too small' )
          STOPINT = .TRUE.
        END IF
      ELSE
        DO I=1,M
          IF (ATOL(I).LE.0D0) THEN
        CALL Rprinti1('Atol is too small for index', I)
            STOPINT=.TRUE.
          END IF
        END DO
      END IF

      IF (WORK(2) .EQ. 0D0) THEN
         HMAX = (TEND-T0)/8d0
      ELSE
         HMAX = WORK(2)
         IF (HMAX.LT.0D0) HMAX=(TEND-T0)/8D0
         IF (HMAX.GT.(TEND-T0)) HMAX = TEND-T0
      END IF

      IF (WORK(3) .EQ. 0D0) THEN
         FACNEWTV(1) = 1D-1
      ELSE
         FACNEWTV(1) = WORK(3)
         IF ((FACNEWTV(1).LE.0D0).OR.(FACNEWTV(1).GE.1D0)) THEN
        CALL Rprintd1('Wrong input work(3)=',WORK(3))
             STOPINT=.TRUE.
         END IF
      END IF

      IF (WORK(4) .EQ. 0D0) THEN
         FACNEWTV(2) = 1D-1
      ELSE
         FACNEWTV(2) = WORK(4)
         IF ((FACNEWTV(2).LE.0D0).OR.(FACNEWTV(2).GE.1D0)) THEN
        CALL Rprintd1('Wrong input work(4) = ',WORK(4))
             STOPINT=.TRUE.
         END IF
      END IF

      IF (WORK(5) .EQ. 0D0) THEN
         FACNEWTV(3) = 1D-1
      ELSE
         FACNEWTV(3) = WORK(5)
         IF ((FACNEWTV(3).LE.0D0).OR.(FACNEWTV(3).GE.1D0)) THEN
        CALL Rprintd1('Wrong input work(5) = ',WORK(5))
             STOPINT=.TRUE.
         END IF
      END IF

      IF (WORK(6) .EQ. 0D0) THEN
         FACNEWTV(4) = 1D-1
      ELSE
         FACNEWTV(4) = WORK(6)
         IF ((FACNEWTV(4).LE.0D0).OR.(FACNEWTV(4).GE.1D0)) THEN
        CALL Rprintd1('Wrong input work(6) = ',WORK(6))
             STOPINT=.TRUE.
         END IF
      END IF

      IF (WORK(7) .EQ. 0D0) THEN
         FACNEWTV(5) = 1D-1
      ELSE
         FACNEWTV(5) = WORK(7)
         IF ((FACNEWTV(5).LE.0D0).OR.(FACNEWTV(5).GE.1D0)) THEN
        CALL Rprintd1('Wrong input work(7) = ',WORK(7))
             STOPINT=.TRUE.
         END IF
      END IF

      IF (WORK(8) .EQ. 0D0) THEN
         FACNSMALL = 1d-2
      ELSE
         FACNSMALL = WORK(8)
         IF ((FACNSMALL.LE.0D0).OR.(FACNSMALL.GE.1D0)) THEN
        CALL Rprintd1('Wrong input work(8) = ',WORK(8))
             STOPINT=.TRUE.
         END IF
      END IF

      IF (WORK(9) .EQ. 0D0) THEN
         FACNRESTR = 5d-2
      ELSE
         FACNRESTR = WORK(9)
         IF ((FACNRESTR.LE.0D0).OR.(FACNRESTR.GE.1D0)) THEN
        CALL Rprintd1('Wrong input work(9) = ',WORK(9))
             STOPINT=.TRUE.
         END IF
      END IF

      IF (WORK(10) .EQ. 0D0) THEN
         FACL = 1.2D-1
      ELSE
         FACL = WORK(10)
         IF (FACL.LT.0D0) THEN
        CALL Rprintd1('Wrong input work(10) = ',WORK(10))
             STOPINT=.TRUE.
         END IF
      END IF

      IF (WORK(11) .EQ. 0D0) THEN
                FACR = 10D0
      ELSE
         FACR = WORK(11)
         IF(FACR.LE.0D0) THEN
        CALL Rprintd1('Wrong input work(11) = ',WORK(11))
             STOPINT=.TRUE.
         END IF
         IF(FACL.GE.FACR) THEN
        CALL Rprintd1('Invalid values for work(10),FACL ', WORK(10))
        CALL Rprintd1('Invalid values for work(10),FACR ', WORK(11))
            STOPINT=.TRUE.
         END IF
      END IF

      IF (WORK(12) .EQ. 0D0) THEN
         SFTY = 1D0/20D0
      ELSE
         SFTY = WORK(12)
         IF(SFTY.LE.0D0) THEN
        CALL Rprintd1('Wrong input work(13) = ',WORK(12))
             STOPINT=.TRUE.
         END IF
      END IF

      IF (WORK(13) .EQ. 0D0) THEN
         SFTYUP = .5d0*SFTY
      ELSE
         SFTYUP = WORK(13)
         IF(SFTYUP.LE.0D0) THEN
        CALL Rprintd1('Wrong input work(13) = ',WORK(13))
             STOPINT=.TRUE.
         END IF
      END IF

      IF (WORK(14) .EQ. 0D0) THEN
         SFTYDN = SFTY
      ELSE
         SFTYDN = WORK(14)
         IF(SFTYDN.LE.0D0) THEN
        CALL Rprintd1('Wrong input work(14) = ',WORK(14))
             STOPINT=.TRUE.
         END IF
      END IF


      IF (STOPINT) THEN
C     INVALID INPUT PARAMETERS
         IDID = -1
         RETURN
      END IF

C---------------------------------------------------------
C     FIXED PARAMETERS
C---------------------------------------------------------

      RHOMUV(1) = 1d-2*DABS(DLOG10(DMIN1(RTOL,1D-1)))
      RHOMLV(1) = 5d-1
      DO I=2,NMETH
         RHOMUV(I) = RHOMUV(I-1)**(DBLE(STEP_ORD(I))
     &                             /DBLE(STEP_ORD(I-1)))
         RHOMLV(I) = RHOMLV(I-1)**(DBLE(STEP_ORD(I))
     &                             /DBLE(STEP_ORD(I-1)))
      END DO

C--------------------------------------------------------
C     BANDED MATRIX
C--------------------------------------------------------

        JBAND = (MLJAC .LT. M)
        IF (JBAND) THEN
                LDJAC   = MLJAC+MUJAC+1
                LDLU    = LDJAC+MLJAC
                IJOB(1) = 2
        ELSE
          LDJAC = M
                LDLU     = M
                IJOB(1)  = 1
        END IF

      IF (IMAS.EQ.0) THEN
           IWORK(9) = M
           MBAND     = .FALSE.
           LDMAS     = 1
      ELSE
             MBAND = (MLMAS .LT. M)
             IF (MBAND) THEN
               LDMAS   = MLMAS+MUMAS+1
               IJOB(2) = 2
           ELSE
               LDMAS   = M
               IJOB(2) = 1
           END IF
      END IF


C---------------------------------------------------------
C     COMPUTE THE VECTORS ENTRY-POINT IN IWORK AND WORK
C---------------------------------------------------------


      INDIPVT = 41
      IF ((INDIPVT + M-1) .GT. LIWORK) THEN
        CALL Rprinti1('Insuff. storage for iwork, min. = ',INDIPVT+M-1)
        IDID = -1
        RETURN
      END IF


      INDF0       = 15
      INDT        = INDF0       + M
      INDY        = INDT        + STEP_ORD(indord)
      INDF        = INDY        + M*STEP_ORD(indord)
      INDTHETA    = INDF        + M*STEP_ORD(indord)
      INDERR      = INDTHETA    + M*LDLU
      INDTEMP     = INDERR      + M*STEP_ORD(indord)
      INDSCAL     = INDTEMP     + M*STEP_ORD(indord)
      INDTOLEXT   = INDSCAL     + M
      INDSCALEXT  = INDTOLEXT   + M
      INDJ0       = INDSCALEXT  + M
      INDM0       = INDJ0       + M*LDJAC
      INDDELJ0    = INDM0       + M*LDMAS
      INDDELJ0OLD = INDDELJ0    + M
      IND_DD      = INDDELJ0OLD + M
      INDFJ0      = IND_DD      + M*(STEP_ORD(indord)+1)
      INDEJ0      = INDFJ0      + M

      IF ((INDEJ0 + M-1) .GT. LWORK) THEN
        CALL Rprinti1('Insuff. storage for work, min.=',INDEJ0+M-1)
        IDID = -1
        RETURN
      END IF

      TOLESTRAPR   = DMIN1(1D-2,1D2*RTOL)
      IF (ITOL.EQ.0) THEN
         DO I=1,M
            WORK(INDTOLEXT+I-1) = DMIN1(1D-2,1D2*ATOL(1))
         END DO
      ELSE
         DO I=1,M
            WORK(INDTOLEXT+I-1) = DMIN1(1D-2,1D2*ATOL(I))
         END DO
      END IF

      CALL  BIM0(M,FCN,JAC,NMETH,STEP_ORD(indord),Y0,WORK(INDF0),
     &           T0,TEND,H,RTOL,ATOL,ITOL,
     &           MAXSTEP,ORDMIN,ORDMAX,ITMAX,UROUND,HMAX,FACNEWTV,
     &           FACNSMALL,FACNRESTR,FACL,FACR,SFTY,SFTYUP,SFTYDN,
     &           RHOMUV,RHOMLV,
     &           IWORK(12),IWORK(13),IWORK(14),IWORK(15),IWORK(16),
     &           IWORK(21),IWORK(26),IWORK(31),IWORK(36),
     &           IWORK(INDIPVT),STEP_ORD,
     &           WORK(INDT),WORK(INDY),WORK(INDF),WORK(INDTHETA),
     &           WORK(INDERR),WORK(INDJ0),WORK(INDDELJ0),
     &           WORK(INDDELJ0OLD),WORK(INDFJ0),WORK(INDEJ0),
     &           WORK(INDSCAL),WORK(IND_DD),
     &           TOLESTRAPR,WORK(INDTOLEXT),
     &           WORK(INDSCALEXT),
     &           IJAC,MLJAC,MUJAC,LDJAC,LDLU,JBAND,IJOB,
     &           RPAR,IPAR,IOUT,SOLOUT,IDID,
     &           MAS,IMAS,MLMAS,MUMAS,LDMAS,MBAND,
     &           IWORK(9),IWORK(10),
     &           WORK(INDM0),WORK(INDTEMP))


      RETURN
      END

c------------------------------------------------------------------------------------------------------
C    CORE INTEGRATOR
C------------------------------------------------------------------------------------------------------

      SUBROUTINE  BIM0(M,FCN,JAC,NMETH,KMAX,Y0,F0,T0,TEND,
     &                 H,RTOL,ATOL,ITOL,
     &                 MAXSTEP,ORDMIN,ORDMAX,ITMAX,UROUND,HMAX,
     &                 FACNEWTV,FACNSMALL,FACNRESTR,FACL,FACR,SFTY,
     &                 SFTYUP,SFTYDN,RHOMUV,RHOMLV,
     &                 NFEVAL,NJEVAL,NLU,NLINSYS,NITER,NSTEP,NACCEPT,
     &                 NFAILERR,NFAILNEWT,
     &                 IPVT,STEP_ORD,T,Y,F,THETA,ERR,J0,
     &                 DELJ0,DELJ00,FJ0,EJ0,SCAL,
     &                 DD,TOLESTRAPR,TOLESTRAPA,SCALEXTRAP,
     &                 IJAC,MLJAC,MUJAC,LDJAC,LDLU,JBAND,IJOB,
     &                 RPAR,IPAR,IOUT,SOLOUT,IDID,
     &                 MAS,IMAS,MLMAS,MUMAS,LDMAS,MBAND,
     &                 INDEX1,INDEX2,M0,TEMP)

      IMPLICIT NONE
C
C INPUT PARAMETERS
C

      EXTERNAL FCN,JAC,MAS,SOLOUT
      INTEGER M,NMETH,KMAX,MAXSTEP,ORDMIN,ORDMAX,ITMAX(NMETH),
     &        IPVT(M),STEP_ORD(NMETH),INDEX1,INDEX2,
     &        MLJAC,MUJAC,LDJAC,LDLU,IJOB(2),IJAC,IPAR(*),
     &        ITOL,IMAS,MLMAS,MUMAS,LDMAS

      LOGICAL JBAND,MBAND

      DOUBLE PRECISION TEND,H,RTOL,ATOL(*),UROUND,HMAX,
     &                 FACNEWTV(NMETH),FACNSMALL,FACNRESTR,FACL,FACR,
     &                 SFTY,SFTYUP,SFTYDN,
     &                 T(KMAX),F0(M),THETA(LDLU,M),J0(LDJAC,M),DELJ0(M),
     &                 DELJ00(M),Y(M,KMAX),F(M,KMAX),ERR(M,KMAX),
     &                 RHOMUV(NMETH),RHOMLV(NMETH),SCAL(M),DD(KMAX+1,M),
     &                 TOLESTRAPR,TOLESTRAPA(M),FJ0(M),SCALEXTRAP(M),
     &                 EJ0(M),RPAR(*),M0(LDMAS,M),TEMP(M,KMAX)

C
C I/O PARAMETERS
C

      DOUBLE PRECISION T0,Y0(M)

C
C OUTPUT PARAMETERS
C

      INTEGER NFEVAL,NJEVAL,NLU,NLINSYS,NITER(NMETH),NSTEP(NMETH),
     &        NACCEPT(NMETH),NFAILERR(NMETH),NFAILNEWT(NMETH),
     &        IOUT,IDID


C
C LOCAL VARIABLES
C
      INTEGER I,J,IT,MAXIT,ORD,ORD_IND,ORDNEW,K,KOLD,KNEW,KUP,K0,
     &        NFAILCONV,NORD,INFO,IRTRN,IERR,IERR0,MUDIF,
     &        NFAILCONS,NSTEPS,NSING,NERROR,
     &        NFERRCONS,IT0,INDEXD,MINIT
      DOUBLE PRECISION NERR,NERRUP,NERRUP1,NERRDOWN,
     &                 NERR0,NERROLD,NERRSTOP,NERRLOC,NERRLOC0,
     &                 RHO,RHO0,RHOMU,RHOML,RHOT,RHOTUP,
     &                 RHOOLD,RHOI,RHOIUP,
     &                 HNEW,HNUP,HNDN,H0,H00,
     &                 RATH,RATRHO,FI,
     &                 ESP,ESPUP,ESPDN,VMAX(8),
     &                 FACNEWT,GAMMA,HGAMMA,MAXDELTA,NF0,NF,
     &                 MINY0,FMINY0,MAXF0,
     &                 NU, NUUP, NU1, NUUP1, FATERR,HNUP1,FI1,
     &                 SCALJ0_1,NJ0,NJ00,HJ0,RHOJ0,tolrhoJ0,
     &                 DJ0,FATDJ0,FATDJ0I,ABSY0,
     &                 DELTAH2,DELTAH1SF,DELTAH1,CFAT1,CFAT2,CFAT3,
     &                 RATHH,HFATT,ALFAFATT,RHOFATT,DISCR,
     &                 DELT,YSAFE,ITNEW,RHONEW,cscal,cscal0,
     &                 SISERR,SISERRUP,CSIS,CFACT,TNEXT



      LOGICAL JVAI,LAST,EXTRAP,EXTRAP0,EXTRAPS,
     &        CALJAC,CALFACT,SUCCESS,RESTRICT,
     &        TRUEJAC,STAGNA,
     &        LINEAR,CALJAC0,CALFACT0,
     &        QINF,QINFJ,QINFF,NQINF,SMALLM,
     &        NODJ0,NODJ00,ISNANUM,ERROR,ESTIM,ESTIM1

C
C  CONSTANT PARAMETERS
C

C  GAMMA
      double precision gamma4,gamma6,gamma8,gamma10,gamma12
      parameter(gamma4 =.7387D0,
     &          gamma6 =.8482D0,
     &          gamma8 =.7285D0,
     &          gamma10=.6745D0,
     &          gamma12=.6433D0)

C  RHOTILDE
      double precision rhot4,rhot6,rhot8,rhot10,rhot12
      parameter(rhot4 =.5021d0,
     &          rhot6 =.8975d0,
     &          rhot8 =.9178d0,
     &          rhot10=.9287d0,
     &          rhot12=.9361d0)

C RHOTILDEINF
      double precision rhoi4,rhoi6,rhoi8,rhoi10,rhoi12
      parameter(rhoi4 =  .9201d0,
     &          rhoi6 = 1.2475d0,
     &          rhoi8 = 1.7294d0,
     &          rhoi10= 2.0414d0,
     &          rhoi12= 2.2621d0)
C ERROR ESTIMATE

      DOUBLE PRECISION VMAX4_1,VMAX4_2,VMAX4_2_2,
     &                 VMAX6_1,VMAX6_2,VMAX6_2_2,
     &                 VMAX8_1,VMAX8_2,VMAX8_2_2,
     &                 VMAX10_1,VMAX10_2,VMAX10_2_2,
     &                 VMAX12_1,VMAX12_2,VMAX12_2_2

      parameter(VMAX4_1 = 1D0/15D0,VMAX4_2 = 1D0/4D0,
     &          VMAX4_2_2 = 2D0/3D0,
     &          VMAX6_1 =  4D0/45D0,VMAX6_2 =  1D0/5D0,
     &          VMAX6_2_2 =  5D0/6D0,
     &          VMAX8_1 = 81D0/28D2,VMAX8_2 = 1D0/7D0,
     &          VMAX8_2_2 = 7D0/1D1,
     &          VMAX10_1 = 73D0/4619D0,VMAX10_2 = 1D0/9D0,
     &          VMAX10_2_2 = 761D0/1260D0,
     &          VMAX12_1 = 62D0/9913D0,VMAX12_2 =  1D0/11D0,
     &          VMAX12_2_2 = 671D0/1260D0)

C JACOBIAN EVALUATION
      double precision fatdJ04,fatdJ06,fatdJ08,fatdJ010,fatdJ012,
     &                 fatdJ04i,fatdJ06i,fatdJ08i,fatdJ010i,fatdJ012i,
     &                 alfajac,scalJ0,
     &                 tolrhoJ0_1,toldJ0,
     &                 tolrhoJ4,tolrhoJ6,tolrhoJ8,tolrhoJ10,
     &                 tolrhoJ12
      integer  itmaxJ0,itmaxJ0_1
C     If hnew<h the Jacobian is not evaluated, if the estimated spectral
C     radius is less than alfajac.
      parameter( alfajac=.1d0)

C     Scaling factor used in the estimate of dJ0.
      parameter (scalJ0     =1d-3,
     &           tolrhoJ0_1 =1d-1,
     &           toldJ0     =1d-8,
     &           itmaxJ0    =3,
     &           itmaxJ0_1  =6)

      parameter(tolrhoJ4=5d-3,
     &          tolrhoJ6=4d-3,
     &          tolrhoJ8=3d-3,
     &          tolrhoJ10=2d-3,
     &          tolrhoJ12=1d-3)

C     If deltaJ0/J0 > fatdJ0 the Jacobian must be evaluated.
      parameter( fatdJ04 = 2d-2,
     &           fatdJ06 = 1d-2,
     &           fatdJ08 = 1d-3,
     &           fatdJ010= 2d-4,
     &           fatdJ012= 3d-5)

      parameter( fatdJ04i = 5d-2,
     &           fatdJ06i = 4d-2,
     &           fatdJ08i = 3d-2,
     &           fatdJ010i= 2d-2,
     &           fatdJ012i= 1d-2)

C FACTORIZATION
      double precision deltah2_4,deltah2_6,deltah2_8,deltah2_10,
     &     deltah2_12,
     &     deltah1_4sf,deltah1_6sf,deltah1_8sf,deltah1_10sf,
     &     deltah1_12sf
      double precision cfat4_1 ,cfat4_2, cfat6_1, cfat6_2,  cfat8_1,
     &                 cfat8_2, cfat10_1,cfat10_2,cfat12_1,cfat12_2

C     The factorization of theta is not computed if
C         max(deltah1,deltah1sf)<hnew/h<deltah2
      parameter (deltah2_4   =1.10d0,
     &           deltah2_6   =1.09d0,
     &           deltah2_8   =1.08d0,
     &           deltah2_10  =1.07d0,
     &           deltah2_12  =1.06d0,
     &           deltah1_4sf =.90d0,
     &           deltah1_6sf =.91d0,
     &           deltah1_8sf =.92d0,
     &           deltah1_10sf=.93d0,
     &           deltah1_12sf=.94d0)

C     Parameter cfat involved in computing deltah1 (see above).
      parameter (cfat4_1=-1.4487d0,  cfat4_2=2.3593d0,
     &           cfat6_1=-1.4983d0,  cfat6_2=3.1163d0,
     &           cfat8_1=-1.4662d0,  cfat8_2=3.5197d0,
     &           cfat10_1=-1.4290d0, cfat10_2=3.7538d0,
     &           cfat12_1=-1.3964d0, cfat12_2=3.9104d0)

C ORDER REDUCTION RECOVERY
      double precision faterr4,faterr6,faterr8,faterr10,
     &                 rath1,rath2,ratrho1,ratrho2

      parameter( rath1   = .95d0,   rath2 = 1.05d0,
     &           ratrho1 = .95d0, ratrho2 = 1.05d0)

      parameter (faterr4=7D0,faterr6=6D0,faterr8=5D0,
     &           faterr10 = 4D0)


C VARIOUS PARAMETERS
      double precision facu1,facu2,facnocon,faclro,
     &     rhobad,tolminy0,tolminf0,tolmaxf0
      integer flmx,flhlt

C     Interval for the new stepsize for which order can be augmented:
C              facu2*h<=hnew<=facu1*h
      parameter(facu1=1.25d0,facu2=.8d0)

C     Factors for determining the new stepsize in case of Newton
C     failure and large spectral radius, respectively.
C     (In both cases, the order is decreased).
      parameter(facnocon=.5d0, faclro=.5d0)

C     Max number of consecutive failures, after which a constant
C     initial guess is used.
      parameter(flmx=1)

C     Max number of consecutive failures, after which integration stops.
      parameter(flhlt=10)

C     Maximum value of the spectral radius, for having failure.
      parameter(rhobad =.99d0)

C     Parameters for the stop criterion, when the solution is 'small'.
      parameter(tolminy0=1d-2,
     &          tolminf0=1d-4,
     &          tolmaxf0=1d-3)

      double precision cscal4,cscal6,cscal8,cscal10,cscal12
      parameter( cscal4=16d0,
     &           cscal6=40d0,
     &           cscal8=16d1,
     &           cscal10=9d2,
     &           cscal12=7d3)


C-------------------------------------------------------------------------------------------------------------------
      IRTRN = 1

C STATISTICS
C
      NFEVAL = 0
      NJEVAL = 0
      NLU    = 0
      NLINSYS= 0
      NSTEPS = 0
      DO I=1,NMETH
         NITER(I)     = 0
         NSTEP(I)     = 0
         NACCEPT(I)   = 0
         NFAILERR(I)  = 0
         NFAILNEWT(I) = 0
      END DO

C-------------------------------------------------------------------------------------------------------------------
C     OTHER INITIALIZATIONS
CF
      qinfj = .FALSE.
      qinff = .FALSE.

      IF (INDEX1.EQ.M) THEN
         INDEXD = 1
      ELSEIF (INDEX1+INDEX2.EQ.M) THEN
         INDEXD = 2
      ELSE
         INDEXD = 3
      ENDIF

      IF (JBAND) THEN
         CSIS  = DBLE(2*M*(MLJAC+MUJAC))
         CFACT = DBLE(2*M*MLJAC*MUJAC)
         MUDIF = MUJAC -MUMAS
      ELSE
         CSIS  = DBLE(2*M*M)
         CFACT = DBLE(2*M*M*M)/3D0
      END IF
      SMALLM   = M.LE.5

C     VECTOR TO BE USED FOR THE ESTIMATE OF JACOBIAN VARIATION
      NERR     = 0.3141592654D0
      NERRUP   = 0D0
      DO I = 1,M
         NERR   = 4D0*NERR*(1D0-NERR)
         EJ0(I) = NERR
         IF (NERR.LT.5D-1) EJ0(I)=EJ0(I)+5D-1
         NERRUP = DMAX1(NERRUP,EJ0(I))
      END DO
      DO I = 1,M
         EJ0(I)=EJ0(I)/NERRUP
      END DO


      NFAILCONV  = 0
      NFERRCONS  = 0
      NFAILCONS  = 0

      ORD     = ORDMIN
      ORD_IND = INT(ORD/2) - 1
      K       = STEP_ORD(ORD_IND)
      KOLD    = 0
CF Initialize the parameters
      GAMMA  = 0.0d0
      RHOT   = 0.0d0
      RHOTUP = 0.0d0
      RHOI   = 0.0d0
      RHOIUP = 0.0d0
      RHOOLD = 0.0d0
      RHOMU  = 0.0d0
      RHOML  = 0.0d0
      RHOJ0  = 0.0d0
      RHOFATT = 0.0d0
      SISERR = 0
      SISERRUP  = 0
      FATERR = 0.0d0
      VMAX(1)  = 0.0d0
      VMAX(2)  = 0.0d0
      VMAX(3)  = 0.0d0
      VMAX(4)  = 0.0d0
      VMAX(5)  = 0.0d0
      VMAX(6)  = 0.0d0
      VMAX(7)  = 0.0d0
      VMAX(8)  = 0.0d0
      FATDJ0 = 0.0d0
      FATDJ0I = 0.0d0
      tolrhoJ0 = 0.0d0
      DELTAH2= 0.0d0
      DELTAH1SF = 0.0d0
      CFAT1  = 0.0d0
      CFAT2  = 0.0d0
      cscal = 0.0d0
      NORD     = 0
      NJ00 = 0
      NERROR = 0
      MINIT = 0
      MAXIT = 0
      MAXf0 = 0
      it0 = 0
      it = 0
      hj0 = 0
      fminy0 = 0.0d0
      extraps = .false.
      espup = 0.0d0
      espdn = 0.0d0
      esp = 0.0d0
      error = .false.
      mudif = 0
CF
      H       = DMIN1(H,HMAX)
C -------------------------------------------------------
C     INITIAL STEPSIZE TOO SMALL !!!!!
      IF (H .LT. 1D1*UROUND*T0) H = 2D1*UROUND*T0
C -------------------------------------------------------

      RHO      =  0D0
      MAXDELTA =  0D0

      LAST     = .FALSE.
      EXTRAP   = .FALSE.
      EXTRAP0  = .FALSE.
      RESTRICT = .FALSE.
      SUCCESS  = .FALSE.
      CALJAC   = .TRUE.
      CALJAC0  = .TRUE.
      CALFACT0 = .TRUE.
      LINEAR   = .FALSE.
      QINF     = .FALSE.

      NJ0     = 0D0
      NERRLOC = 0D0
      H0      = 0D0
      NODJ0   = .FALSE.
      DO J =1,KMAX
        DO I = 1,M
          ERR(I,J) = 0D0
        END DO
      END DO



      IERR = 0
      IF (IMAS.NE.0) THEN
C Karline: changed
          CALL MAS(M,M0,LDMAS,RPAR,IPAR)
          IF (IERR.NE.0) THEN
        CALL Rprinti1('Error Ierr returned by the subroutine MAS',IERR)
        CALL Rprintd1('Exit at t = ', T0)
             IDID = -6
             GOTO 800
          END IF
      END IF
C Karline: changed
      IERR = 0
      CALL FCN(M,T0,Y0,F0,RPAR,IPAR)
C      CALL FCN(M,T0,Y0,F0,IERR,RPAR,IPAR)

      IF (IERR.NE.0) THEN
        CALL Rprinti1('Error Ierr returned by the subroutine FCN',IERR)
        CALL Rprintd1('Exit at t = ', T0)
          IDID = -6
          GOTO 800
      END IF
      NFEVAL = NFEVAL + 1
      SCALJ0_1 = 0D0
      DO I=1,M
         ABSY0 = DABS(Y0(I))
         SCALEXTRAP(I)=1d0/(1d0+ABSY0)
         SCALJ0_1=DMAX1(SCALJ0_1,ABSY0)
      END DO

      HFATT=1D0

      MINY0 = TOLMINY0 + 1D0

c     MAIN LOOP
100   CONTINUE

      IF (K .EQ. KOLD) GOTO 140

C     THE ORDER OF THE METHOD HAS BEEN CHANGED
      ESP = 1D0/DBLE(K+1)
      IF (ORD .LT. ORDMAX)
     &    ESPUP=1D0/DBLE(STEP_ORD(ORD_IND+1)+1D0)
      IF (ORD .GT. ORDMIN)
     &    ESPDN=1D0/DBLE(STEP_ORD(ORD_IND-1)+1D0)
      MAXIT    = ITMAX(ORD_IND)
      RHOML    = RHOMLV(ORD_IND)
      RHOMU    = RHOMUV(ORD_IND)
      RHOOLD   = 0D0
      NORD     = 0
      NERROR   = 0
      ERROR    = .FALSE.
      MINIT    = INDEXD
C     DEFINE THE PARAMETERS DEPENDING ON THE ORDER OF THE METHOD
      if (ord_ind .eq. 1) then
        goto 105  
      else if (ord_ind .eq. 2) then
        goto 110  
      else if (ord_ind .eq. 3) then
        goto 115  
      else if (ord_ind .eq. 4) then
        goto 120  
      else if (ord_ind .eq. 5) then
        goto 125  
      endif  
C      GOTO(105,110,115,120,125) ORD_IND
105   GAMMA  = gamma4
      RHOT   = rhot4
      RHOTUP = rhot6
      RHOI   = rhoi4
      RHOIUP = rhoi6
      SISERR = 2
      SISERRUP  = 3
      FATERR = faterr4
      VMAX(1)  = VMAX4_1
      VMAX(2)  = gamma4*VMAX4_2
      VMAX(3)  = gamma4*gamma4*VMAX4_2_2
      VMAX(7)  = gamma4*VMAX6_2
      VMAX(8)  = gamma4*gamma4*VMAX6_2_2
      FATDJ0 = fatdj04
      FATDJ0I = fatdj04i
      tolrhoJ0 = tolrhoJ4
      DELTAH2= deltah2_4
      DELTAH1SF = deltah1_4sf
      CFAT1  = cfat4_1
      CFAT2  = cfat4_2
      cscal = cscal4
      GOTO 140
110   GAMMA  = gamma6
      RHOT   = rhot6
      RHOTUP = rhot8
      RHOI   = rhoi6
      RHOIUP = rhoi8
      SISERR = 3
      SISERRUP  = 3
      FATERR = faterr6
      VMAX(1)  = VMAX6_1
      VMAX(2)  = gamma6*VMAX6_2
      VMAX(3)  = gamma6*gamma6*VMAX6_2_2
      VMAX(4)  = VMAX4_1
      VMAX(5)  = gamma6*VMAX4_2
      VMAX(6)  = gamma6*gamma6*VMAX4_2_2
      VMAX(7)  = gamma6*VMAX8_2
      VMAX(8)  = gamma6*gamma6*VMAX8_2_2
      FATDJ0 = fatdj06
      FATDJ0I = fatdj06i
      tolrhoJ0 = tolrhoJ6
      DELTAH2= deltah2_6
      DELTAH1SF = deltah1_6sf
      CFAT1  = cfat6_1
      CFAT2  = cfat6_2
      cscal = cscal6
      GOTO 140
115   GAMMA  = gamma8
      RHOT   = rhot8
      RHOTUP = rhot10
      RHOI   = rhoi8
      RHOIUP = rhoi10
      SISERR = 3
      SISERRUP  = 3
      FATERR = faterr8
      VMAX(1)  = VMAX8_1
      VMAX(2)  = gamma8*VMAX8_2
      VMAX(3)  = gamma8*gamma8*VMAX8_2_2
      VMAX(4)  = VMAX6_1
      VMAX(5)  = gamma8*VMAX6_2
      VMAX(6)  = gamma8*gamma8*VMAX6_2_2
      VMAX(7)  = gamma8*VMAX10_2
      VMAX(8)  = gamma8*gamma8*VMAX10_2_2
      FATDJ0 = fatdj08
      FATDJ0i = fatdj08i
      tolrhoJ0 = tolrhoJ8
      DELTAH2= deltah2_8
      DELTAH1SF = deltah1_8sf
      CFAT1  = cfat8_1
      CFAT2  = cfat8_2
      cscal = cscal8
      GOTO 140
120   GAMMA  = gamma10
      RHOT   = rhot10
      RHOTUP = rhot12
      RHOI   = rhoi10
      RHOIUP = rhoi12
      SISERR = 3
      SISERRUP  = 3
      FATERR = faterr10
      VMAX(1)  = VMAX10_1
      VMAX(2)  = gamma10*VMAX10_2
      VMAX(3)  = gamma10*gamma10*VMAX10_2_2
      VMAX(4)  = VMAX8_1
      VMAX(5)  = gamma10*VMAX8_2
      VMAX(6)  = gamma10*gamma10*VMAX8_2_2
      VMAX(7)  = gamma10*VMAX12_2
      VMAX(8)  = gamma10*gamma10*VMAX12_2_2
      FATDJ0 = fatdj010
      FATDJ0i = fatdj010i
      tolrhoJ0 = tolrhoJ10
      DELTAH2= deltah2_10
      DELTAH1SF = deltah1_10sf
      CFAT1  = cfat10_1
      CFAT2  = cfat10_2
      cscal = cscal10
      GOTO 140
125   GAMMA  = gamma12
      RHOT   = rhot12
      RHOI   = rhoi12
      SISERR = 3
      VMAX(1)  = VMAX12_1
      VMAX(2)  = gamma12*VMAX12_2
      VMAX(3)  = gamma12*gamma12*VMAX12_2_2
      VMAX(4)  = VMAX10_1
      VMAX(5)  = gamma12*VMAX10_2
      VMAX(6)  = gamma12*gamma12*VMAX10_2_2
      FATDJ0 = fatdj012
      FATDJ0i = fatdj012i
      tolrhoJ0 = tolrhoJ12
      DELTAH2= deltah2_12
      DELTAH1SF = deltah1_12sf
      CFAT1  = cfat12_1
      CFAT2  = cfat12_2
      cscal = cscal12

140   CONTINUE

C-------------------------------------------------------------------------------------------------------------------
C     JACOBIAN EVALUATION
      LINEAR = .FALSE.
      ESTIM =(.NOT.SMALLM).AND.
     &       (SUCCESS.OR.(.NOT.SUCCESS.AND.CALJAC))
     &       .AND.(IMAS.EQ.0)
      IF (ESTIM) THEN
        SCALJ0_1 = 0D0
        DO I = 1,M
           SCALJ0_1 = DMAX1(SCALJ0_1,DABS(Y0(I)))
        END DO
        SCALJ0_1  = scalJ0*(1d0+SCALJ0_1)
        DO I=1,M
          DELJ0(I) = Y0(I) + EJ0(I)*SCALJ0_1
        END DO
        IERR=0
C Karline: changed
      IERR = 0
        CALL FCN(M,T0,DELJ0,FJ0,RPAR,IPAR)
        NFEVAL = NFEVAL + 1
        NODJ0 = (IERR.NE.0)
        NJ00 = NJ0
        IF (.NOT. NODJ0) THEN
            SCALJ0_1 = 1D0/SCALJ0_1
            NJ0 = 0D0
            DO I=1,M
              DELJ0(I) = (FJ0(I)-F0(I))*SCALJ0_1
              NJ0 = DMAX1(NJ0,DABS(DELJ0(I)))
            END DO
        END IF
      END IF

      IF (SUCCESS) THEN

        CALJAC=CALJAC0.OR.(IMAS.NE.0)
        IF (CALJAC) GOTO 150

        ESTIM1 = ESTIM.AND..NOT.NODJ00.AND..NOT.NODJ0
        IF (ESTIM1) THEN
           DJ0  = 0D0
           DO I=1,M
             DJ0   = DMAX1(DJ0,
     &                   DABS(DELJ0(I)-DELJ00(I)))
           END DO
           CALJAC=(DJ0.GT.1D2*UROUND*NJ00)
           LINEAR =.NOT.CALJAC
           IF (LINEAR) GOTO 150
        END IF

        CALJAC=(K.NE.KOLD).OR.(.NOT.EXTRAP)
        IF (CALJAC) GOTO 150

        CALJAC=(RHO.GE.tolrhoJ0).AND.(IT.GE.itmaxJ0)
        IF (.NOT.CALJAC.OR..NOT.ESTIM1) GOTO 150

        NQINF  = .NOT.(QINF.OR.QINFJ)
        CALJAC=.NOT.
     &      ( ((RHO.LT.5D-2) .OR. (IT .LT. 4)) .AND.
     &      ( (QINF.AND.(DJ0 .LE. fatdJ0i*NJ00)).OR.
     &        ((.NOT.QINF).AND.NQINF.AND.(RHOJ0.LT.1D0).AND.
     &         (HJ0*DJ0 .LE. fatdJ0*RHOJ0/RHOT)) ))
      END IF

150   CONTINUE

      TRUEJAC = CALJAC .OR.(.NOT.SUCCESS)
      IF (CALJAC) THEN
           IF (IJAC.EQ.0) THEN
C     NUMERICAL JACOBIAN
               DO I=1,M
                 YSAFE=Y0(I)
                 DELT = DSQRT(UROUND*DMAX1(1.D-5,DABS(YSAFE)))
                 Y0(I) = YSAFE+DELT
                 IERR = 0
C Karline: changed
      IERR = 0
                 CALL FCN(M,T0,Y0,FJ0,RPAR,IPAR)
                 IF (IERR.NE.0) THEN
        CALL Rprinti1('Error Ierr returned by the subroutine FNC',IERR)
        CALL Rprint('during the numerical evaluation of the Jacobian')
        CALL Rprintd1('Exit at t = ', T0)
                    GOTO 800
                 END IF
                 IF (JBAND) THEN
                    DO J=MAX(1,I-MUJAC),MIN(M,I+MLJAC)
                       J0(J-I+MUJAC+1,I)=(FJ0(J)-F0(J))/DELT
                    END DO
                 ELSE
                    DO J=1,M
                       J0(J,I)=(FJ0(J)-F0(J))/DELT
                    END DO
                 END IF
                 Y0(I)=YSAFE
               END DO
           ELSE
C     ANALYTICAL JACOBIAN
              IERR = 0
C karline: changed
              CALL JAC(M,T0,Y0,MLJAC,MUJAC,J0,LDJAC,RPAR,IPAR)
              IF (IERR.NE.0) THEN
        CALL Rprinti1('Error Ierr returned by the subroutine JAC',IERR)
        CALL Rprintd1('Exit at t = ', T0)
                 IDID = -6
                 GOTO 800
              END IF
           END IF

           NJEVAL = NJEVAL + 1
           IF (ESTIM) THEN
              DO I=1,M
                DELJ00(I) = DELJ0(I)
              END DO
           ENDIF
           NJ00   = NJ0
           NODJ00 = NODJ0

           IF ((SMALLM).AND.(IMAS.EQ.0)) THEN
            IF (JBAND) THEN
             NJ0 = 0D0
             DO I=1,M
                DELT=0D0
                DO J=MAX(1,I-MLJAC),MIN(M,MUJAC+I)
                   DELT=DELT+DABS(J0(MUJAC+1-J+I,J))
                END DO
                NJ0 = DMAX1(NJ0,DELT)
             END DO
            ELSE
              NJ0 = 0D0
              DO I=1,M
                 DELT = 0D0
                 DO J=1,M
                   DELT = DELT +DABS(J0(I,J))
                 END DO
                 NJ0 = DMAX1(NJ0,DELT)
              END DO
            END IF
            NJ0 = NJ0/DBLE(M)
           END IF
      END IF

C-------------------------------------------------------------------------------------------------------------------
c     COMPUTE AND FACTORIZE THE ITERATION MATRIX THETA

      RATHH = H/HFATT
      CALFACT=((KOLD.NE.K).OR.(.NOT.SUCCESS).OR.CALJAC .OR.
     &         (M.LT.2*K).OR.(IMAS.NE.0).OR.
     &         (QINF.AND.(DABS(RATHH-1D0).GT.fatdJ0i))
     &         .OR. ((.NOT.QINF).AND.
     &        ((RATHH.GT.DELTAH2).OR.(RATHH.LT.DELTAH1SF)))
     &        .OR. CALFACT0
     &        .OR.((RHO.GT.5D-2).AND.(IT.GE.4)) )
      IF (CALFACT.OR.QINF) GOTO 160

      CALFACT = QINFF.OR.(RHOFATT.GE.1D0)

      IF (CALFACT.OR.(RATHH.GE.1D0).OR.(IT.EQ.1)) GOTO 160

      RHONEW = RHO * H/H0
      ITNEW = DBLE(IT)*DLOG(RHO)/DLOG(RHONEW)
      ALFAFATT = (DBLE(M)/(DBLE(6*K)*ITNEW)) + 1D0
      CFAT3 = RHOT/(GAMMA*RHOFATT)*
     &            (DELTAH1SF*RHOFATT)**(1D0/ALFAFATT)
      CFAT3 = CFAT2-CFAT3*CFAT3
      DISCR = CFAT1*CFAT1 - CFAT3

      IF (DISCR .GE. 0D0) THEN
         DELTAH1 = -(CFAT1 + DSQRT(DISCR))
      ELSE
         DELTAH1 = 1D0
      ENDIF
      CALFACT=RATHH.LT.DELTAH1

160   CONTINUE

      IF (CALFACT) THEN
         NSING = 0
170      CONTINUE
         HGAMMA = H*GAMMA
         IF (IMAS.NE.0) THEN
           IF (JBAND) THEN
             DO J=1,M
               DO I=1,LDJAC
                   THETA(I+MLJAC,J)=-HGAMMA*J0(I,J)
               END DO
               DO I=1,LDMAS
                  THETA(I+MLJAC+MUDIF,J)=M0(I,J)+THETA(I+MLJAC+MUDIF,J)
               END DO
             END DO
           ELSEIF (MBAND) THEN
             DO J=1,M
               DO I=1,M
                   THETA(I,J) = -HGAMMA*J0(I,J)
               END DO
               DO I=MAX(1,J-MUMAS),MIN(M,J+MLMAS)
                   THETA(I,J)=THETA(I,J)+M0(I-J+MUMAS+1,J)
               END DO
             END DO
           ELSE
             DO J=1,M
                DO I=1,M
                   THETA(I,J) = M0(I,J)-HGAMMA*J0(I,J)
                END DO
             END DO
            END IF
         ELSE
           IF (JBAND) THEN
               DO J=1,M
                 DO I=1,LDJAC
                    THETA(I+MLJAC,J)=-HGAMMA*J0(I,J)
                 END DO
                 THETA(LDJAC,J)=1D0 +THETA(LDJAC,J)
               END DO
           ELSE
               DO J=1,M
                 DO I=1,M
                    THETA(I,J) = -HGAMMA*J0(I,J)
                 END DO
                 THETA(J,J) = 1D0 + THETA(J,J)
               END DO
           END IF
         END IF
         CALL DECLU(M,THETA,LDLU,MLJAC,MUJAC,IPVT,IJOB,INFO)
         NLU = NLU + 1
         IF (INFO.NE.0) THEN
               NSING = NSING + 1
               IF (NSING.GT.5) THEN
        CALL Rprinti1('Matrix is repeatedly singular, IER= ',INFO )
        CALL Rprintd1('Exit at t = ', T0)
                  IDID = -4
                  GOTO 800
               ELSE
                  H = H*.5D0
                  IF (.1d0*DABS(H) .LE. DABS(T0)*UROUND) THEN
        CALL Rprintd1('Stepsize too small, h = ',H)
        CALL Rprintd1('Exit at t = ', T0)
                    IDID = -3
                    GOTO 800
                  END IF
                  GOTO 170
               END IF
         END IF
      END IF

      CALJAC0 = .FALSE.
      CALFACT0 = .FALSE.

C-------------------------------------------------------------------------------------------------------------------
C     SCALING
      IF (ITOL.EQ.0) THEN
        DO I =1,INDEX1
          SCAL(I) = 1D0/(ATOL(1)+RTOL*DABS(Y0(I)))
        END DO
        cscal0 = DMIN1(1D0,RTOL/(UROUND*1D2*cscal))
        DO I =INDEX1+1,INDEX1+INDEX2
          SCAL(I) = DMIN1(1D0,cscal0*H)/(ATOL(1)+RTOL*DABS(Y0(I)))
        END DO
        cscal0 = DMIN1(1D0,RTOL/(UROUND*1D2*cscal*cscal))
        DO I =INDEX1+INDEX2+1,M
          SCAL(I) = DMIN1(1D0,cscal0*H*H)/
     &            (ATOL(1)+RTOL*DABS(Y0(I)))
        END DO
      ELSE
        DO I =1,INDEX1
          SCAL(I) = 1D0/(ATOL(I)+RTOL*DABS(Y0(I)))
        END DO
        cscal0 = DMIN1(1D0,RTOL/(UROUND*1D2*cscal))
        DO I =INDEX1+1,INDEX1+INDEX2
          SCAL(I) = DMIN1(1D0,cscal0*H)/
     &                   (ATOL(I)+RTOL*DABS(Y0(I)))
        END DO
        cscal0 = DMIN1(1D0,RTOL/(UROUND*1D2*cscal*cscal))
        DO I =INDEX1+INDEX2+1,M
          SCAL(I) = DMIN1(1D0,cscal0*H*H)/
     &            (ATOL(I)+RTOL*DABS(Y0(I)))
        END DO
      END IF

      DO I=1,K
         T(I)=T0+I*H
      END DO

C-------------------------------------------------------------------------------------------------------------------
C     STOPPING CRITERION WHEN THE SOLUTION HAS SMALL ENTRIES

      IF (SUCCESS) THEN
         IF (ITOL.EQ.0) THEN
           IF (IMAS.EQ.0) THEN
             MINY0  = DABS(Y0(1))
             FMINY0 = DABS(F0(1))*SCAL(1)*ATOL(1)
             MAXF0  = FMINY0
             DO I = 2,M
               IF (DABS(Y0(I)).LT.MINY0) THEN
                 MINY0 = DABS(Y0(I))
                 FMINY0 = DABS(F0(I))*SCAL(I)*ATOL(1)
               END IF
               MAXF0 = DMAX1(MAXF0,DABS(F0(I))*SCAL(I)*ATOL(1))
             END DO
           ELSE
             MINY0  = DABS(Y0(1))
             FMINY0 = DABS(Y0(1)-Y(1,K0-1))/H0*SCAL(1)*ATOL(1)
             MAXF0  = FMINY0
             DO I = 2,INDEX1
               FACNEWT = DABS(Y0(I)-Y(I,K0-1))/H0*SCAL(I)*ATOL(1)
               IF (DABS(Y0(I)).LT.MINY0) THEN
                  MINY0 = DABS(Y0(I))
                 FMINY0 = FACNEWT
               END IF
               MAXF0 = DMAX1(MAXF0,FACNEWT)
             END DO
           END IF
         ELSE
           IF (IMAS.EQ.0) THEN
             MINY0  = DABS(Y0(1))
             FMINY0 = DABS(F0(1))*SCAL(1)*ATOL(1)
             MAXF0  = FMINY0
             DO I = 2,M
               IF (DABS(Y0(I)).LT.MINY0) THEN
                 MINY0 = DABS(Y0(I))
                 FMINY0 = DABS(F0(I))*SCAL(I)*ATOL(I)
               END IF
               MAXF0 = DMAX1(MAXF0,DABS(F0(I))*SCAL(I)*ATOL(I))
             END DO
           ELSE
             MINY0  = DABS(Y0(1))
             FMINY0 = DABS(Y0(1)-Y(1,K0-1))/H0*SCAL(1)*ATOL(1)
             MAXF0  = FMINY0
             DO I = 2,INDEX1
               FACNEWT = DABS(Y0(I)-Y(I,K0-1))/H0*SCAL(I)*ATOL(I)
               IF (DABS(Y0(I)).LT.MINY0) THEN
                  MINY0 = DABS(Y0(I))
                 FMINY0 = FACNEWT
               END IF
               MAXF0 = DMAX1(MAXF0,FACNEWT)
             END DO
           END IF
         END IF
      END IF

      FACNEWT = FACNEWTV(ORD_IND)
      IF (((MINY0 .LT. TOLMINY0).AND.(FMINY0.LT.1d-1)))
     &THEN
          FACNEWT = DMIN1(1D-1,FACNEWT)
          IF ((FMINY0.LT.1D-5).AND.(MAXF0.LT.1D-2))
     &          FACNEWT=DMIN1(FACNEWT,FACNSMALL)
      END IF
      IF(RESTRICT) FACNEWT=DMIN1(FACNEWT,FACNRESTR)

C-------------------------------------------------------------------------------------------------------------------
C     SOLUTION INITIALIZATION

      IF ( (.NOT.EXTRAP).OR.(NFAILCONV.GT.flmx) ) THEN
C     THE SOLUTION IS INITIALIZED WITH THE CONSTANT PROFILE
         IF (NFAILCONV.GT.flhlt) THEN
        CALL Rprint('Too many consecutive Newton failures' )
        CALL Rprintd1('Exit at t = ', T0)
             IDID = -5
             GOTO 800
         END IF
         DO J = 1,K
            DO I=1,M
               Y(I,J)=Y0(I)
            END DO
         END DO
      ELSE
         CALL EXTRAPOLA(M,K0,K,H0,H,Y,DD)
      END IF


C-------------------------------------------------------------------------------------------------------------------
C     NEWTON ITERATION

      IT        = 0
      RHO       = 0D0
      NERR0     = 1D0
      NERRSTOP  = DMAX1(FACNEWT,UROUND/RTOL)

200   CONTINUE
      IERR0=0
      DO I=1,K
C karline: changed
        IERR = 0
        CALL FCN(M,T(I),Y(1,I),F(1,I),RPAR,IPAR)
        IERR0 = IERR0+IERR
      END DO
      NFEVAL = NFEVAL + K
      IF (IERR0.NE.0) THEN
         NERR = 2D0*NERRSTOP + 1D0
         EXTRAP0 = .FALSE.
         EXTRAPS = .FALSE.
         JVAI = .FALSE.
         GOTO 305
      END IF

      if (ord_ind .eq. 1) then
        goto 210  
      else if (ord_ind .eq. 2) then 
        goto 220  
      else if (ord_ind .eq. 3) then
        goto 230  
      else if (ord_ind .eq. 4) then
        goto 240  
      else if (ord_ind .eq. 5) then
        goto 250  
      endif  
C      GOTO (210,220,230,240,250) ORD_IND

210   CALL blendstep4(M,Y0,F0,Y,F,H,THETA,IPVT,ERR,
     &                GAMMA,LDLU,MLJAC,MUJAC,IJOB,IMAS,
     &                LDMAS,MLMAS,MUMAS,M0,TEMP)
      GOTO 300

220   CALL blendstep6(M,Y0,F0,Y,F,H,THETA,IPVT,ERR,
     &                GAMMA,LDLU,MLJAC,MUJAC,IJOB,IMAS,
     &                LDMAS,MLMAS,MUMAS,M0,TEMP)
      GOTO 300

230   CALL blendstep8(M,Y0,F0,Y,F,H,THETA,IPVT,ERR,
     &                GAMMA,LDLU,MLJAC,MUJAC,IJOB,IMAS,
     &                LDMAS,MLMAS,MUMAS,M0,TEMP)
      GOTO 300

240   CALL blendstep10(M,Y0,F0,Y,F,H,THETA,IPVT,ERR,
     &                 GAMMA,LDLU,MLJAC,MUJAC,IJOB,IMAS,
     &                 LDMAS,MLMAS,MUMAS,M0,TEMP)
      GOTO 300

250   CALL blendstep12(M,Y0,F0,Y,F,H,THETA,IPVT,ERR,
     &                 GAMMA,LDLU,MLJAC,MUJAC,IJOB,IMAS,
     &                 LDMAS,MLMAS,MUMAS,M0,TEMP)

300   CONTINUE

      IT = IT + 1
      CALL NORM(M,K,SCAL,ERR,NERR,NERRUP)
c CHECK FOR NANs
      IF (ISNANUM(NERR).OR.ISNANUM(NERRUP)) THEN
         NERR = 2D0*NERRSTOP + 1D0
         EXTRAP0 = .FALSE.
         EXTRAPS = .FALSE.
         JVAI = .FALSE.
         GOTO 305
      END IF

C SPECTRAL RADIUS ESTIMATE

      NERROLD = NERR0
      NERR0   = NERR
      RHO0    = RHO
      RHO     = NERR0/NERROLD
      IF (IT.GT.2) RHO = DSQRT(RHO0*RHO)

      JVAI = (NERR .GT. NERRSTOP) .AND. (IT .LE. MAXIT)
     &       .AND.(IT .LE. INDEXD+1 .OR. RHO .LE. RHOBAD)

      JVAI = JVAI .OR.
     &       (IMAS.NE.0 .AND. IT.LT.MINIT)
305   IF (JVAI)  GOTO 200

C END OF NEWTON'S ITERATION
C-------------------------------------------------------------------------------------------------------------------
      NLINSYS  = NLINSYS  + 2*IT*K
      NITER(ORD_IND) = NITER(ORD_IND) + IT
      NSTEP(ORD_IND) = NSTEP(ORD_IND) + 1

      IF (NERR.GT.NERRSTOP) GOTO 308

      IERR0 = 0
      DO I=1,K
C karline:changed
        IERR = 0
        CALL FCN(M,T(I),Y(1,I),F(1,I),RPAR,IPAR)
        IERR0 = IERR0 + IERR
      END DO

      NFEVAL = NFEVAL + K
      IF (IERR0.NE.0)  THEN
        NERR = 2D0*NERRSTOP + 1D0
        EXTRAP0 = .FALSE.
        EXTRAPS = .FALSE.
      END IF

308   CONTINUE

      IF (NERR .GT. NERRSTOP) THEN
C     NEWTON HAS FAILED
          NFAILNEWT(ORD_IND) = NFAILNEWT(ORD_IND) + 1
          NFAILCONV  = NFAILCONV  + 1
          NFAILCONS  = MAX(NFAILCONS + 1,2)
          NFERRCONS  = 0

          H       = facnocon * H
          KOLD    = K
          IF (ORD .GT. ORDMIN) THEN
              ORD     = ORD - 2
              ORD_IND = INT(ORD/2) - 1
              K       = STEP_ORD(ORD_IND)
          END IF
          SUCCESS  = .FALSE.
          LAST     = .FALSE.

          CALJAC   = .NOT.TRUEJAC
          TRUEJAC  = .TRUE.

          IF (NFAILCONV.EQ.1)
     &       EXTRAPS = EXTRAP0.AND.(IT.GT.MAXIT).AND.(RHO.LT.RHOBAD)

          EXTRAP   = EXTRAPS
          RESTRICT = .FALSE.

          MINIT    = 1

          GOTO 550
      END IF

      MINIT = INDEXD

C-------------------------------------------------------------------------------------------------------------------
C     LOCAL ERROR ESTIMATION

      NERRLOC0 = NERRLOC
      
      if (ord_ind .eq. 1) then
        goto 310  
      else if (ord_ind .le. 5) then 
        goto 320  
      endif  
C      GOTO (310,320,320,320,320) ORD_IND
310   CALL localerr4(M,F0,F,H,ERR,SCAL,NERR,NERRUP,NLINSYS,
     &               THETA,VMAX,IPVT,LDLU,MLJAC,MUJAC,IJOB,
     &               IMAS,LDMAS,MLMAS,MUMAS,M0,K,ORD_IND,
     &               INDEX1,INDEX2)
      GOTO 400

320   CALL localerr(M,F0,F,H,ERR,SCAL,NERR,NERRUP,NLINSYS,
     &              THETA,VMAX,IPVT,LDLU,MLJAC,MUJAC,IJOB,
     &              IMAS,LDMAS,MLMAS,MUMAS,M0,K,ORD_IND,
     &              INDEX1,INDEX2)

400   CONTINUE

      IF (ISNANUM(NERR).OR.ISNANUM(NERRUP)) THEN
          NERR = 2D0*NERRSTOP + 1D0
          EXTRAP0 = .FALSE.
          EXTRAPS = .FALSE.
          GOTO 308
      END IF

      NFAILCONV = 0
      NERRLOC = NERR

      IF ((IMAS.NE.0).AND.(NFERRCONS.GT.3).AND.
     &    (DABS(NERRLOC - NERRLOC0).LT.1D-4)
     &    .AND. (NERR.GE.1D0)) THEN
c         cwrite(MSG,*) 'WARNING: POSSIBLE INCONSISTENT INITIAL VALUE'
         NERR = .9D0
      END IF

      IF (NERR .GT. 0d0) THEN
          HNEW = H*(SFTY/NERR)**ESP
      ELSE
          HNEW = facr*H
      END IF

      IF (NERR .GE. 1D0) THEN
C     FAILURE DUE TO LOCAL ERROR TEST
          HNEW = H*(1D-1/NERR)**ESP
          NFAILERR(ORD_IND) = NFAILERR(ORD_IND) + 1
          NFAILCONS = MAX(NFAILCONS + 1,2)
          NFERRCONS = NFERRCONS + 1
          CALJAC   = .NOT. TRUEJAC
          TRUEJAC  = .TRUE.
          SUCCESS  = .FALSE.
          LAST     = .FALSE.
          KOLD     = K
          IF (ERROR .AND.(ORD.GT.ORDMIN)) THEN
              ORD = ORD - 2
              ORD_IND = INT(ORD/2) -1
              K = STEP_ORD(ORD_IND)
              H = DMAX1(HNEW/2D0,facl*H)
          ELSE
              H = DMAX1(HNEW,facl*H)
          END IF
          ERROR = .TRUE.
          GOTO 550
      END IF

C     STEP ACCEPTED
      NACCEPT(ORD_IND) = NACCEPT(ORD_IND) + 1
      NFERRCONS = 0

      SUCCESS   = .TRUE.
      NORD      =  NORD + 1
      NFAILCONS = MAX(NFAILCONS - 1,0)

      EXTRAP0 = .FALSE.
      DO I=1,INDEX1
         MAXDELTA = DABS(Y(I,K)-Y0(I))*SCALEXTRAP(I)
         EXTRAP0 = EXTRAP0 .OR.
     &      ((dabs(Y0(I)).LE.1D-1).AND.(MAXDELTA.GE.TOLESTRAPA(I))).OR.
     &      ((dabs(Y0(I)).GT.1D-1).AND.(MAXDELTA.GE.TOLESTRAPR))
      END DO
      RESTRICT = .NOT.EXTRAP0
      DO I=INDEX1+1,M
         MAXDELTA = DABS(Y(I,K)-Y0(I))*SCALEXTRAP(I)
         EXTRAP0 = EXTRAP0 .OR.
     &      ((dabs(Y0(I)).LE.1D-1).AND.(MAXDELTA.GE.TOLESTRAPA(I))).OR.
     &      ((dabs(Y0(I)).GT.1D-1).AND.(MAXDELTA.GE.TOLESTRAPR))
      END DO
      EXTRAP   =  EXTRAP0
      IF ((EXTRAP).OR.(IOUT.EQ.1)) CALL DIFFDIV(M,K,Y0,Y,DD)

      IF (IOUT.EQ.1) THEN
        CALL SOLOUT(M,K,ORD,T0,T,Y,F,DD,RPAR,IPAR,IRTRN)
        IF (IRTRN.LT.0) GOTO 800
      END IF

      IF (LAST) GOTO 600

C     ORDER VARIATION
      RATH   = HNEW/H
      IF (IT .LT. 2) RHO=1D-3
      RATRHO = RHOOLD/RHO
      KNEW   = K
      ORDNEW = ORD
      QINF   = (NERRUP.GE.NERR).AND.(NERR.NE.0D0)
      IF (TRUEJAC) THEN
         HJ0=H
         RHOJ0 = RHO
         IF (IT.LE.1) RHOJ0=1D0
         QINFJ = QINF
      END IF
      IF (CALFACT) THEN
         HFATT=H
         RHOFATT = RHO
         IF (IT.LE.1) RHOFATT=1D0
         QINFF   = QINF
      END IF

      IF ((ERROR).AND.(ORD.GT.ORDMIN)) THEN
         NERROR = MIN(NERROR+1,2)
         IF (NERROR.GT.1) THEN
             ORDNEW = ORD -2
             ORD_IND = ORD_IND - 1
             KNEW = STEP_ORD(ORD_IND)
             HNEW = DMIN1(H,HNEW)/2D0
             CALJAC = .TRUE.
             TRUEJAC = .TRUE.
             GOTO 500
         END IF
       ELSE
         NERROR = 0
       END IF
       ERROR = .FALSE.

      IF ( RHO.LE.RHOML) GOTO 445
c     SPECTRAL RADIUS TOO LARGE
      IF ((ORD.EQ.ORDMIN).OR.(IT.LE.INDEXD+2)) GOTO 500

      CALL errdown(M,F0,F,H,ERR(1,2),SCAL,NERRDOWN,NLINSYS,
     &             VMAX(4),QINF,THETA,IPVT,LDLU,MLJAC,MUJAC,IJOB,
     &             K,ORD_IND,INDEX1,INDEX2)


      IF (NERRDOWN .GT. 0D0) THEN
         HNDN = H*(SFTYDN/NERRDOWN)**ESPDN
      ELSE
         HNDN = HNEW
      END IF

      IF (QINF.AND.(INDEXD.LE.1).AND.(HNDN.LT.HNEW)) GOTO 500
      IF (.NOT. QINF) HNDN=DMIN1(HNDN,HNEW)
      ORDNEW  = ORD - 2
      ORD_IND = ORD_IND - 1
      KNEW    = STEP_ORD(ORD_IND)
      HNEW=HNDN
      GOTO 500

445   CONTINUE

      IF ((NFAILCONS.GT.0).OR.(ORD.GE.ORDMAX)
     &    .OR.(NORD.LT.3)
     &    .OR.(RATH.GE.FACU1).OR.(RATH.LE.facu2)
     &    .OR.(T(K)+DBLE(K)*HNEW.GE.TEND))
     &GOTO 500

      IF ((NERRUP.LT.1D1*UROUND).AND.(.NOT.QINF).AND.
     &    (IT.LT.3).AND.(RHO.LT.1D1*UROUND).AND.
     &    (IMAS.EQ.0))
     &THEN
C     THE PROBLEM IS A QUADRATURE
        CALL ERRUP(M,K,ORD_IND,ERR,H,H0,H00,VMAX(7),NERRUP,
     &             SCAL,THETA,IPVT,LDLU,MLJAC,MUJAC,IJOB,
     &             INDEX1,INDEX2)
        NLINSYS = NLINSYS + 1
        IF (NERRUP .GT. 0D0) THEN
           HNUP = H*(SFTYUP/NERRUP)**ESPUP
        ELSE
           HNUP = facr*H
        END IF
        IF (T(K)+DBLE(STEP_ORD(ORD_IND+1))*HNUP.GE.TEND)
     &  GOTO 500
        NU   = DBLE(IT)
        NUUP = DBLE(IT)
        KUP=STEP_ORD(ORD_IND+1)
        FI    = (HNEW/HNUP)*(DBLE(K)/DBLE(KUP))*
     &          (CSIS*(2D0*KUP*NUUP+SISERRUP)+CFACT)/
     &          (CSIS*(2D0*K*NU+SISERR)+CFACT)
        IF (FI.LT.1D0) THEN
         ORDNEW  = ORD + 2
         ORD_IND = ORD_IND + 1
         KNEW    = STEP_ORD(ORD_IND)
         HNEW    = HNUP
        END IF
        GOTO 500
       END IF

      IF (QINF) GOTO 450

      IF ((RHO.GE.RHOMU).AND.(IT.GT.INDEXD)) THEN
         IF (IMAS.NE.0) GOTO 450
         STAGNA=(DABS(RATH-1D0).LT.1D-1).AND.
     &          (DABS(RATRHO-1D0).LT.1D-1)
         IF (.NOT.STAGNA.OR.(T(K).EQ.0D0)) GOTO 450
         IERR=0
C karline: changed
         CALL FCN(M,T(K)*(1D0+1D-5),Y(1,K),FJ0,RPAR,IPAR)
         NFEVAL = NFEVAL + 1
         IF (IERR.NE.0) GOTO 450
         NF0= 0d0
         NF = 0d0
         DO I=1,M
            NF0 = DMAX1(NF0,DABS(FJ0(I)-F(I,K)))
            NF  = DMAX1(NF,DABS(F(I,K)))
         END DO
         NF0 = NF0/(1d-5*T(K))
         IF(NJ0*NF.GE.2D0*NF0) GOTO 450
      END IF

      IF (NERRUP .GT. 0D0) THEN
         HNUP = H*(SFTYUP/NERRUP)**ESPUP
      ELSE
         HNUP = facr*H
      END IF
      IF (T(K)+DBLE(STEP_ORD(ORD_IND+1))*HNUP.GE.TEND)
     &GOTO 500

      NU    = DMAX1(DBLE(INDEXD),(DBLE(IT) *DLOG(RHO))/DLOG(RHO*RATH))
      NUUP  = DMAX1(DBLE(INDEXD),
     &          (DBLE(IT) *DLOG(RHO))/DLOG(RHO*(RHOTUP/RHOT)*(HNUP/H)))

      KUP=STEP_ORD(ORD_IND+1)
      FI    = (HNEW/HNUP)*(DBLE(K)/DBLE(KUP))*
     &        (CSIS*(2D0*KUP*NUUP+SISERRUP)+CFACT)/
     &        (CSIS*(2D0*K*NU+SISERR)+CFACT)

      IF (FI.LT.1D0)
     &THEN
         ORDNEW  = ORD + 2
         ORD_IND = ORD_IND + 1
         KNEW    = STEP_ORD(ORD_IND)
         HNEW    = HNUP
         GOTO 500
      END IF

450   CONTINUE
C ---------------------------------------------------------------
C ORDER REDUCTION RECOVERY
C ---------------------------------------------------------------
      STAGNA = (RATH .GT. rath1).AND. (RATH .LT. rath2)
      IF ( (INDEXD.LE.1).OR.
     &     (IT.GT.INDEXD+1).OR.(IT0.GT.INDEXD+1) )
     &THEN
         STAGNA= STAGNA.AND.
     &       (RATRHO .GT. ratrho1).AND.(RATRHO.LT.ratrho2)
      END IF

      IF ( QINF .OR.
     &     (STAGNA.AND.(FATERR*NERRUP .GE. NERR))  )
     &THEN

             IF ((H/H0.GE.facu1).OR.(H/H0.LE.facu2) )
     &       GOTO 500
             IF ( (ORD_IND.GT.1).AND.
     &            ((H0/H00.GE.facu1).OR.(H0/H00.LE.facu2)) )
     &       GOTO 500

             IF ((.NOT.CALFACT.OR.(.NOT.TRUEJAC.AND..NOT.LINEAR)).AND.
     &          QINF) THEN
c IN THIS CASE, THE SPECTRAL RADIUS DOESN'T BEHAVES LIKE rho=rhoti/q.
                CALJAC0 = .NOT.LINEAR
                CALFACT0 = .TRUE.
                GOTO 500
             END IF

             IF (.NOT.QINF) THEN
                IF (IMAS.EQ.0) THEN
                    HNEW = HNEW *(vmax(1)/vmax(2))**ESP
                ELSE
                   DO I=1,INDEX1
                      ERR(I,2)=ERR(I,2)*vmax(2)/vmax(1)
                   END DO
                   CALL NORM(M,1,SCAL,ERR(1,2),NERR,NERROLD)
                   IF (NERR.GT.0) THEN
                      HNEW = H*(SFTY/NERR)**ESP
                   ELSE
                      HNEW = facr*H
                   END IF
                END IF
                RATH = HNEW/H
             END IF

             CALL ERRUP(M,K,ORD_IND,ERR,H,H0,H00,VMAX(7),NERRUP1,
     &                  SCAL,THETA,IPVT,LDLU,MLJAC,MUJAC,IJOB,
     &                  INDEX1,INDEX2)
             NLINSYS = NLINSYS + 1

             IF (NERRUP1 .GT. 0D0) THEN
                 HNUP1 = H*(SFTYUP/NERRUP1)**ESPUP
             ELSE
                 HNUP1 = facr*H
             END IF
             IF (T(K)+DBLE(STEP_ORD(ORD_IND+1))*HNUP1.GE.TEND)
     &       GOTO 500

             NU1    = DMAX1((DBLE(IT)*DLOG(RHO))/DLOG(RHO/RATH),
     &                       DBLE(INDEXD))
             NUUP1  = DMAX1(DBLE(INDEXD),(DBLE(IT)*DLOG(RHO))/
     &                DLOG(RHO*(H/HNUP1)*(RHOIUP/RHOI)) )


             KUP=STEP_ORD(ORD_IND+1)
             FI1    = (HNEW/HNUP1)*(DBLE(K)/DBLE(KUP))*
     &                (CSIS*(2D0*KUP*NUUP1+SISERRUP)+CFACT)/
     &                (CSIS*(2D0*K*NU1+SISERR)+CFACT)
             IF   (FI1 .LT. 1D0)
     &       THEN
                   ORDNEW  = ORD + 2
                   ORD_IND = ORD_IND + 1
                   KNEW    = STEP_ORD(ORD_IND)
                   HNEW    = HNUP1
             END IF
      END IF

500   CONTINUE
      DO I=1,M
         Y0(I)  =Y(I,K)
         ABSY0  =DABS(Y0(I))
         F0(I)  =F(I,K)
         SCALEXTRAP(I)=1d0/(1d0+ABSY0)
         IF (ORD.LT.ORDMAX) THEN
            ERR(I,K+2)=ERR(I,K+1)
            ERR(I,K+1)=ERR(I,1)
         END IF
      END DO

      T0  = T(K)

      K0   = K
      KOLD = K
      IT0  = IT

      K    = KNEW
      ORD  = ORDNEW

      H00 = H0
      H0  = H
      IF (NFAILCONS.GT.0) HNEW = DMIN1(H,HNEW)
      HNEW = DMAX1(HNEW,facl*H)
      H    = DMIN1(HNEW,facr*H,HMAX)

      RHOOLD    = RHO

550   CONTINUE

      NSTEPS=NSTEPS+1

      IF (.1d0*DABS(T0-TEND)/DBLE(K).GT.dabs(T0)*UROUND) THEN
         IF (.1d0*DABS(H) .LE. DABS(T0)*UROUND) THEN
        CALL Rprintd1('Stepsize too small, h = ',H)
        CALL Rprintd1('Exit at t = ', T0)
            IDID = -3
            GOTO 800
         END IF
         IF (NSTEPS .GE. MAXSTEP) THEN
        CALL Rprinti1('More than nmax steps are needed', MAXSTEP)
        CALL Rprintd1('Exit at t = ', T0)
            IDID = -2
            GOTO 800
         END IF

         TNEXT = T0 + DBLE(K)*H
         IF ( (IMAS.NE.0) .AND. (TNEXT.LT.TEND) .AND.
     &        ((TEND-TNEXT).LT.(DBLE(K)*H)) )
     &   H = (TEND-T0)/(2D0*DBLE(K))

         IF (TNEXT .GE. TEND) THEN
            H = (TEND-T0)/DBLE(K)
            LAST = .TRUE.
         END IF

         GOTO 100
C     END OF MAIN LOOP
      END IF

C     INTEGRATION END
600   CONTINUE
      IDID = 0
      DO I=1,M
          Y0(I)=Y(I,K)
          F0(I)=F(I,K)
      END DO
      T0 = T(K)


800   RETURN

      END
