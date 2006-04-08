      SUBROUTINE GRAD(X,Y,M,N,H,ICEN,TCEN,XH,
     * R,TOL,IFLAG,WA,GUP)
C
C  X matrix (M BY N)
C  Y vector (M)
C  M = Number of Observations
C  N = Number of Parameters
C  H = basis, integer(N) vector 
C  ICEN (integer(M)) = 2 for deleted obs
C                    = 1 for crossed non-censored obs
C                    = 0 otherwise
C  TCEN(M) are the corresponding tau values where a censored obs is
C    first crossed (ICEN=1): weight = (tau - TCEN(I))/(1 - TCEN(I))
C    TCEN = 1 if censored obs is never crossed
C  XH = (N by N) X(H,) inverse matrix
C  R(M) = residuals
C  TOL = zero tolerance (1.D-10 by default)
C  IFLAG (M)  work vector
C  WA (M by N) work array
C  returns:
C    GUP(N)  = upper bounds on tau
C    IFLAG(1:N) = sign(grad denom)
C               = +1 if bound from lower grad bound
C               = -1 if bound from upper grad bound
C    
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER H(N),ICEN(M),IFLAG(M)
      DOUBLE PRECISION ONE
      DIMENSION X(M,N),Y(M),GUP(N),XH(N,N)
      DIMENSION R(M),TCEN(M),WA(M,N)
      DATA  ZERO/0.00D0/, ONE/1.0d0/
C
C GET X'(XH-INVERSE)
C
      DO 80 I = 1,M
      IF(ICEN(I) .EQ. 2) GO TO 80
      DO 70 J = 1,N
      A = ZERO
        DO 60 K = 1,N
  60    A = A + X(I,K)*XH(K,J)
  70  WA(I,J) = A
  80  CONTINUE
  
C
C  GET GRADIENT BOUNDS
C  FIRST SET IFLAG = 1 FOR BASIS INDICES (TEMPORARY)
C
      T = ZERO
      DO 90 I = 1,M
  90  IFLAG(I) = 0
      DO 95 J = 1,N
  95  IFLAG(H(J)) = 1
      DO 120 J=1,N
      A = ZERO
      B = ZERO
      C = ZERO
      D = ZERO
      DO 100 I = 1,M
        IF(ICEN(I) .EQ. 2) GO TO 100
        IF(ICEN(I) .EQ. 0) THEN
          IF(R(I) .GT. TOL)    A = A + WA(I,J)
          IF(R(I) .LT. - TOL)  B = B + WA(I,J)
          GO TO 100
        ENDIF
        IF(IFLAG(I).EQ.1) GO TO 100
C
C  IFLAG = 0: NOT BASIS  AND  ICEN = 1: CROSSED CENSORED
C
        IF(R(I) .LT. -TOL) THEN
           T = TCEN(I)/(ONE - TCEN(I))
           C = C - T*WA(I,J)
           GO TO 100
        ENDIF
        IF(R(I) .GT. TOL) THEN          
           D = D - WA(I,J)
        ENDIF
 100  CONTINUE
      D = D - C
      L = ICEN(H(J))
      IF(L .NE. 0) T = TCEN(H(J))/( ONE - TCEN(H(J)) )
      S = DFLOAT(L)*(T + ONE) - ONE
      E1 = A + B - D - S
      E2 = A + B - D + ONE
      IF(E1 .GT. ZERO) THEN
        GUP(J) = (B + C - S)/E1
        IFLAG(J+M) = 1
        GO TO 120
      ENDIF
      IF(E2 .LT. ZERO) THEN
        GUP(J) = (B + C)/E2
        IFLAG(J+M) = -1
        GO TO 120
      ENDIF
      GUP(J) = - ONE
 120  CONTINUE
      DO 130 J = 1,N
 130  IFLAG(J) = IFLAG(J+M)
      RETURN
      END
