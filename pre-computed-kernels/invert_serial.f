      PROGRAM invert
      
c     This program computes  the Laplace inverse for the kernels 
c     contained  in the  bimaterial  mode 1&2  code. It uses the 
c     input file 'invert.in'.
      
      IMPLICIT NONE

      REAL*8 alpha,relerr,cdcs,nu
      REAL*8 h11cut,dth11,k12cut,dtk12,h22cut,dth22
      INTEGER mm,kmax,nux 
      CHARACTER*160 word160
      CHARACTER*2 ii
      
      COMMON /para/cdcs,nu

      alpha = 0.d0              ! max value of the pole
      kmax = 200                ! max number of calculations
      
      OPEN(8,FILE='invert_serial.in',STATUS='unknown')
      
      READ(8,100) word160
      READ(8,*) h11cut          ! h11 kernel cutoff T
      READ(8,100) word160
      READ(8,*) dth11           ! h11 time spacing interval      
      READ(8,100) word160
      READ(8,*) k12cut          ! k12 kernel cutoff T
      READ(8,100) word160
      READ(8,*) dtk12           ! k12 time spacing interval
      READ(8,100) word160
      READ(8,*) h22cut          ! h22 kernel cutoff T
      READ(8,100) word160
      READ(8,*) dth22           ! h22 time spacing interval
      READ(8,100) word160       
      READ(8,*) nu              ! poission's ratio

      nux=INT(nu*1000.d0+0.5)/10
      CALL chartrans(ii,nux)
      OPEN(9000,FILE='nu_.'//ii//'_h11.dat',STATUS='unknown',
     $     FORM='unformatted')
      OPEN(9250,FILE='nu_.'//ii//'_k12.dat',STATUS='unknown',
     $     FORM='unformatted')
      OPEN(9500,FILE='nu_.'//ii//'_h22.dat',STATUS='unknown',
     $     FORM='unformatted')
c
c     Dilitational wave (cd) / Shear wave (cs)
c      
      cdcs=DSQRT(2.d0*(1.d0-nu)/(1.d0-2.d0*nu))

      DO mm=1,3
         IF(mm.EQ.1)THEN        ! H11 KERNEL
            relerr=.000001d0
            CALL kernel(dth11,relerr,h11cut,mm,kmax,alpha,9000)
            
         ELSE IF(mm.EQ.2)THEN   ! K12 KERNEL
            relerr=.00001d0
            CALL kernel(dtk12,relerr,k12cut,mm,kmax,alpha,9250)
            
         ELSE IF(mm.EQ.3)THEN  ! H22 KERNEL
            relerr=.000001d0  
            CALL kernel(dth22,relerr,k12cut,mm,kmax,alpha,9500)
         ENDIF
      ENDDO
      
 100  FORMAT(a160)
      
      END
c
c ------------------------------------------------------------- KERNEL
c
      SUBROUTINE kernel(dt,relerr,tcut,mm,kmax,alpha,iunit)

      IMPLICIT NONE
      
      INTEGER numt,kmax,mm,nn,iunit
      INTEGER k, ier
      REAL*8 x,tt,alpha,finv,estim,dt,relerr,tcut
      COMPLEX*16 h11,k12,h22
      PARAMETER (nn=100000)
      DIMENSION finv(nn),k(nn)
      DIMENSION estim(nn),ier(nn),tt(nn)
      
      EXTERNAL h11,k12,h22 

      numt=1
      x=dt
      DO WHILE(x.LE.(tcut+dt))
         numt=numt+1
         tt(numt)=x
         x=dt+x
      END DO

      IF(mm.EQ.1)THEN
         tt(1)=.00001d0
         CALL inverse(h11,numt,tt,alpha,relerr,kmax,finv,k,ier,estim)
         CALL outputres(finv,numt,tcut,dt,iunit,ier,relerr)
         PRINT*,'h11 kernel done.....phew'
      ELSE IF(mm.EQ.2)THEN
         tt(1)=.00001d0
         CALL inverse(k12,numt,tt,alpha,relerr,kmax,finv,k,ier,estim)
         CALL outputres(finv,numt,tcut,dt,iunit,ier,relerr)
         PRINT*,'k12 kernel done.....phew'
      ELSE
         tt(1)=.00001d0
         CALL inverse(h22,numt,tt,alpha,relerr,kmax,finv,k,ier,estim)
         CALL outputres(finv,numt,tcut,dt,iunit,ier,relerr)
         PRINT*,'h22 kernel done.....phew'
      ENDIF

      RETURN
      END
c     
c --------------------------------------------------------- FUNCTION h11
c
      COMPLEX*16 function h11(s)
      IMPLICIT NONE
      REAL*8 cdcs, nu
      COMPLEX*16 s
      COMMON /para/cdcs,nu
      
c     h11 kernel

      h11=(cdcs*cdcs+cdcs*SQRT(s*s+cdcs*cdcs)*SQRT(s*s+1)+cdcs*s*
     $     SQRT(s*s+cdcs*cdcs)-s*SQRT(1+s*s))/((SQRT(1+s*s)+s)*
     $     (s*s+1+cdcs*cdcs)) 
      
      RETURN
      END      
c
c --------------------------------------------------------- FUNCTION k12
c
      COMPLEX*16 function k12(s)
      IMPLICIT NONE
      REAL*8 cdcs,nu
      COMPLEX*16 s
      COMMON /para/cdcs,nu

c     k12 kernel
      
      k12 = cdcs*(s*s/(cdcs-SQRT(s*s+cdcs*cdcs)*SQRT(1+s*s))+1)
C
C      k12=n*(-n+((n*n+1)*s*s+n*n)/(s*s+SQRT(s*s+n*n)*SQRT(s*s+1)))/
C     $    (SQRT(s*s+n*n)*SQRT(s*s+1)-n)

      RETURN
      END      
      
c
c --------------------------------------------------------- FUNCTION h22
c      
      COMPLEX*16 function h22(s)
      IMPLICIT NONE
      REAL*8 cdcs,nu
      COMPLEX*16 s
      COMMON /para/cdcs,nu
      
c     h22 kernel
      
      h22=cdcs*s*(SQRT(1+s*s)*(-cdcs*cdcs)/(s+SQRT(s*s+cdcs*cdcs))+cdcs)
     $     /(SQRT(1+s*s)*SQRT(s*s+cdcs*cdcs)-cdcs)
      
      RETURN
      END      
c     
c ------------------------------------------------------------ OUTPUTRES
c     
      SUBROUTINE outputres(finv,ntim,tcut,dt,iunit,ier,relerr)
      
      IMPLICIT NONE
      
      REAL*8 finv(*),relerr,dt,tcut
      INTEGER i,j,ier(*),iunit,ntim
      real*8 cdcs, nu

      COMMON /para/cdcs,nu

      IF (ier(i).eq.22) THEN 
         PRINT*,' ERROR IN RESULTS, NOT ENOUGH ITERATIONS '
      ELSE IF (ier(i).eq.33) THEN 
         PRINT*,' ERROR IN RESULTS, PROGRAM WAS STOPED '
      ELSE IF (ier(i).eq.44) THEN
         PRINT*,' ERROR IN RESULTS, RELATIVE ERROR > ',relerr
      ELSE IF (ier(i).eq.11) THEN
         PRINT*,'ERROR IN RESULTS, NEED A LARGER MACHINE NUMBER(DMACH)'
      END IF

      WRITE(iunit) tcut,ntim,dt, cdcs, nu
      DO j=1,ntim
         WRITE(iunit) finv(j)
      END DO
      
      RETURN
      END
c
c -------------------------------------------------------------- INVERSE
c 
      SUBROUTINE inverse(f,n,t,alpha,relerr,kmax,finv,kl,ier,estim)

c     SPECIFICATIONS FOR ARGUMENTS
      
      INTEGER n, kmax
      INTEGER, DIMENSION(1:n) :: kl, ier
      DOUBLE PRECISION alpha,relerr,t(n),finv(n),estim(n)
      COMPLEX*16 f
      EXTERNAL  f
      
c     SPECIFICATIONS FOR LOCAL VARIABLES
      
      INTEGER  i, ibeg, iend, ih, ii, iter, itmax, j, jodd,
     &     k, kk, ll, mmmm, ncnt, newtry, nt,jmax
      DOUBLE PRECISION a, al, dela, er, erest, erj, expat, pik, qepal,
     &     QEPALJ, RATIO, REPAL, REPALI, REPALJ, SERR, T2, TCAP
      COMPLEX*16 DIFF, EPSUM(99), FK, FSUM, S(99), SI, SVI, Z, ZK
      
c     SPECIFICATIONS FOR SAVE VARIABLES
      
      DOUBLE PRECISION AX, DEFINT, HALF, ONE, PI, RPEXE, ZERO
      SAVE       AX, DEFINT, HALF, JMAX, ONE, PI, RPEXE, ZERO
c     was not used on HP ---->  IJ, (was in save)

c     SPECIFICATIONS FOR INTRINSICS
      
      INTRINSIC  DABS, DLOG, DMAX1, DMIN1, CDABS, CDEXP, DCMPLX, DEXP,
     &     MAX0, MOD, DREAL
      INTEGER    MAX0, MOD
      
c     SPECIFICATIONS FOR SUBROUTINES
c     SPECIFICATIONS FOR FUNCTIONS
      
      DOUBLE PRECISION DMACH
      
      DATA ZERO/0.0D0/, ONE/1.0D0/, HALF/.5D0/
      DATA RPEXE/-1.0D0/
      RPEXE = -1.D0
      DATA AX/0.1D0/
      AX = 0.D0
      DATA DEFINT/2.92D0/
      JMAX = 99
      PI = 4.0*datan(1.0d0)
      dmach = dtan(pi/2.d0)
      KK = 0
      ITMAX = 8
      IF (RPEXE .LT. 0.0D0) RPEXE=DLOG(DMACH) - 1.0D0
 
c     BEGIN MAIN LOOP
      DO 180  I=1, N
         NEWTRY = 0
         AL = ALPHA
         T2 = DEFINT*T(I)
         TCAP = T2*.5D0
         ITER = 1
 10      A = AL - DLOG(RELERR)/T2
         RATIO = DMIN1(ONE,T(I)/TCAP)
         Z = CDEXP(DCMPLX(ZERO,PI*T(I)/TCAP))
 20      ZK = Z
         K = 1
         FSUM = HALF*F(DCMPLX(A,ZERO))
         IF (A*T(I) .GT. RPEXE) THEN
            ier(i) = 11         ! this was IER=130
            GO TO 9000
         END IF
         EXPAT = DEXP(A*T(I))/TCAP
         JODD = 1
         J = 1
         S(1) = FSUM
         QEPAL = EXPAT*DREAL(FSUM)
         QEPALJ = QEPAL
         PIK = PI
 30      CONTINUE
         FK = F(DCMPLX(A,PIK/TCAP))
c     GENERATE PARTIAL SUM
         FSUM = FSUM + FK*ZK
         ZK = ZK*Z
         PIK = PIK + PI
         K = K + 1
         IF (J .EQ. 1) REPALJ = EXPAT*DREAL(FSUM)
c     line 163                         BEGIN EPSILON ALGORITHM
         JODD = -JODD
         J = J + 1
         S(J) = FSUM
         SVI = DCMPLX(ZERO,ZERO)
         II = J
         DO 40  ILOOP=2, J
            II = II - 1
            DIFF = S(II+1) - S(II)
            SI = SVI
            IF (CDABS(DIFF) .NE. ZERO) SI = SI + DCMPLX(ONE,ZERO)/DIFF
            SVI = S(II)
            S(II) = SI
 40      CONTINUE
         IF (JODD .LT. 0) GO TO 60
         REPALJ = EXPAT*DREAL(S(1))
         IF (J+1 .GE. JMAX) GO TO 50
c     CONVERGENCE CHECK WITHIN ONE
c     SAWTOOTH OF EPSILON ALGORITHM
         ERJ = DABS(REPALJ-QEPALJ)/DMAX1(AX,DABS(REPALJ))
         QEPALJ = REPALJ
         IF (ERJ .GT. RELERR*RATIO) GO TO 60
 50      J = 1
         S(1) = FSUM
         QEPALJ = EXPAT*DREAL(FSUM)
         REPAL = REPALJ
c     CONVERGENCE CHECK BETWEEN TWO
c     line 190                         SAWTEETH
         ER = DABS(REPAL-QEPAL)/DMAX1(AX,DABS(REPAL))
         IF (ER .LE. RATIO*RELERR) GO TO 70
         QEPAL = REPAL
 60      IF (K .LT. KMAX) GO TO 30
         ier(i) = 22            !                   this was 129 
         GO TO 170              !                    this was goto 9000
 70      CONTINUE
         KK = KK + K
         IF (MOD(ITER,2) .EQ. 0) GO TO 80
         ITER = ITER + 1
         A = A + DMAX1(1.15D0/TCAP,0.1D0*DABS(A))
         IF (A*T(I) .GT. RPEXE) THEN
            ier(i) = 33         !               this was 130           
            GO TO 9000
         END IF
         REPALI = REPAL
         GO TO 20
 80      CONTINUE
         ITER = ITER + 1
c     CONVERGENCE CHECK BETWEEN RESULTS
c     FROM THE USE OF TWO SUCCESSIVE A
         EREST = DABS(REPAL-REPALI)/DMAX1(AX,DABS(REPAL),DABS(REPALI))
         IF (EREST .LE. RELERR) GO TO 170
         DELA = 0.1D0*DABS(RELERR/EREST)
         IF (ITER .GT. 3) RATIO = RATIO*DELA
         IF (ITER .LE. ITMAX) GO TO 160
c     GET KMAX PARTIAL SUMS
         SERR = RELERR*AX
         AL = ALPHA
         T2 = DEFINT*T(I)
         TCAP = T2*.5D0
         A = AL - DLOG(SERR)/T2
         Z = CDEXP(DCMPLX(ZERO,PI*T(I)/TCAP))
 90      ZK = Z
         FSUM = HALF*F(DCMPLX(A,ZERO))
         EXPAT = DEXP(A*T(I))/TCAP
         PIK = PI
         JMAX = MAX0(99,KMAX)
         IF (MOD(JMAX,2) .EQ. 0) JMAX = JMAX - 1
         IBEG = JMAX - 98
         II = 1
         LL = JMAX/99
         DO 100  NT=1, JMAX
            FK = F(DCMPLX(A,PIK/TCAP))
            FSUM = FSUM + FK*ZK
            ZK = ZK*Z
            PIK = PIK + PI
            NCNT = NT/LL
            IF (NCNT .GT. 99) GO TO 110
            IF (NCNT .GT. 0) S(NCNT) = FSUM
 100     CONTINUE
c     INITIALIZE EPSUM TO ZERO AND SAVE
c     PSUMS
 110     DO 120  IH=1, 99
            EPSUM(IH) = ZERO
 120     CONTINUE
c     GET ACCELLARATED SUM
         MMMM = -1
         DO 140  IEND=2, 99
            DO 130  IH=99, IEND, MMMM
               DIFF = S(IH) - S(IH-1)
               S(IH) = S(IH-1)
               IF (CDABS(DIFF) .NE. ZERO) S(IH) = ONE/DIFF + EPSUM(IH)
               EPSUM(IH) = S(IH-1)
 130        CONTINUE
 140     CONTINUE
         EPSUM(1) = S(99)
         IF (NEWTRY .EQ. 0) REPAL = EXPAT*DREAL(EPSUM(1))
         IF (NEWTRY .EQ. 1) REPALJ = EXPAT*DREAL(EPSUM(1))
         IF (NEWTRY .EQ. 1) GO TO 150
         NEWTRY = NEWTRY + 1
         A = A + DMAX1(1.15D0/TCAP,0.1D0*DABS(A))
         GO TO 90
 150     EREST = DABS(REPAL-REPALJ)/DMAX1(AX,DABS(REPAL),DABS(REPALJ))
         IF (EREST .LE. RELERR) GO TO 170
         ier(i) = 44            !   this was 129 
         estim(i) = erest       
         GO TO 170
 160     AL = AL - DLOG(DELA)/T2
         GO TO 10
 170     FINV(I) = REPAL
         kl(i)=k
 180  CONTINUE
 9000 CONTINUE
      RETURN
      END
c
c ----------------------------------------------------------- CHARTRANS
c	
      SUBROUTINE chartrans(b,n)

      IMPLICIT NONE

      INTEGER n1,n2,n
      CHARACTER*2 b
      CHARACTER*1 a1,a2
c
c     This subroutine transforms an integer n (1<n<21) into its 
c     corresponding character*2 (e.g., 1 becomes '01', 12 becomes '12')
c
      n1=n/10
      n2=MOD(n,10)
      a1=CHAR(48+n1)
      a2=CHAR(48+n2)
      b=a1//a2

      RETURN
      END
