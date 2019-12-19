module module_qlapack
  implicit none

  private
  public qzsteqr, lu_and_solve

  ! All subroutines and functions in this module is based on LAPACK (v3.4.2)
  ! provided at www.netlib.org.  ! The lapack copyright notice, terms
  ! and conditions, disclaimers are shown below:
  
  ! Copyright (c) 1992-2013 The University of Tennessee and The University
  !                         of Tennessee Research Foundation.  All rights
  !                         reserved.
  ! Copyright (c) 2000-2013 The University of California Berkeley. All
  !                         rights reserved.
  ! Copyright (c) 2006-2013 The University of Colorado Denver.  All rights
  !                         reserved.
  
  ! $COPYRIGHT$
  
  ! Additional copyrights may follow
  
  ! $HEADER$
  
  ! Redistribution and use in source and binary forms, with or without
  ! modification, are permitted provided that the following conditions are
  ! met:
  
  ! - Redistributions of source code must retain the above copyright
  !   notice, this list of conditions and the following disclaimer.
  
  ! - Redistributions in binary form must reproduce the above copyright
  !   notice, this list of conditions and the following disclaimer listed
  !   in this license in the documentation and/or other materials
  !   provided with the distribution.
  
  ! - Neither the name of the copyright holders nor the names of its
  !   contributors may be used to endorse or promote products derived from
  !   this software without specific prior written permission.
  
  ! The copyright holders provide no reassurances that the source code
  ! provided does not infringe any patent, copyright, or any other
  ! intellectual property rights of third parties.  The copyright holders
  ! disclaim any liability to any recipient for claims brought against
  ! recipient by any third party for infringement of that parties
  ! intellectual property rights.
  
  ! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  ! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  ! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  ! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  ! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  ! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  ! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  ! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  ! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  ! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  ! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  
contains

  SUBROUTINE QZSWAP(N,ZX,INCX,ZY,INCY)
    !     .. Scalar Arguments ..
    INTEGER INCX,INCY,N
    !     ..
    !     .. Array Arguments ..
    complex(16) ZX(*),ZY(*)
    !     ..
    !
    !  Purpose
    !  =======
    !
    !     interchanges two vectors.
    !     jack dongarra, 3/11/78.
    !     modified 12/3/93, array(1) declarations changed to array(*)
    !
    !
    !     .. Local Scalars ..
    complex(16) ZTEMP
    INTEGER I,IX,IY
    !     ..
    IF (N.LE.0) RETURN
    IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
    !
    !       code for unequal increments or equal increments not equal
    !         to 1
    !
    IX = 1
    IY = 1
    IF (INCX.LT.0) IX = (-N+1)*INCX + 1
    IF (INCY.LT.0) IY = (-N+1)*INCY + 1
    DO I = 1,N
       ZTEMP = ZX(IX)
       ZX(IX) = ZY(IY)
       ZY(IY) = ZTEMP
       IX = IX + INCX
       IY = IY + INCY
    end DO
    RETURN
    !
    !       code for both increments equal to 1
20  continue
    do I = 1,N
          ZTEMP = ZX(I)
          ZX(I) = ZY(I)
          ZY(I) = ZTEMP
       end do
    RETURN
  end SUBROUTINE QZSWAP

  SUBROUTINE QLASRT( ID, N, D, INFO )
    !
    !  -- LAPACK routine (version 3.0) --
    !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
    !     Courant Institute, Argonne National Lab, and Rice University
    !     September 30, 1994
    !
    !     .. Scalar Arguments ..
    CHARACTER          ID
    INTEGER            INFO, N
    !     ..
    !     .. Array Arguments ..
    real(16) D( * )
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  Sort the numbers in D in increasing order (if ID = 'I') or
    !  in decreasing order (if ID = 'D' ).
    !
    !  Use Quick Sort, reverting to Insertion sort on arrays of
    !  size <= 20. Dimension of STACK limits N to about 2**32.
    !
    !  Arguments
    !  =========
    !
    !  ID      (input) CHARACTER*1
    !          = 'I': sort D in increasing order;
    !          = 'D': sort D in decreasing order.
    !
    !  N       (input) INTEGER
    !          The length of the array D.
    !
    !  D       (input/output) DOUBLE PRECISION array, dimension (N)
    !          On entry, the array to be sorted.
    !          On exit, D has been sorted into increasing order
    !          (D(1) <= ... <= D(N) ) or into decreasing order
    !          (D(1) >= ... >= D(N) ), depending on ID.
    !
    !  INFO    (output) INTEGER
    !          = 0:  successful exit
    !          < 0:  if INFO = -i, the i-th argument had an illegal value
    !
    !  =====================================================================
    !
    !     .. Parameters ..
    INTEGER            SELECT
    PARAMETER          ( SELECT = 20 )
    !     ..
    !     .. Local Scalars ..
    INTEGER            DIR, ENDD, I, J, START, STKPNT
    real(16)   D1, D2, D3, DMNMX, TMP
    !     ..
    !     .. Local Arrays ..
    INTEGER            STACK( 2, 32 )
    !     ..
    !     .. External Functions ..
    LOGICAL            LSAME
    EXTERNAL           LSAME
    !     ..
    !     .. External Subroutines ..
    EXTERNAL           XERBLA
    !     ..
    !     .. Executable Statements ..
    !
    !     Test the input paramters.
    !
    INFO = 0
    DIR = -1
    IF( LSAME( ID, 'D' ) ) THEN
       DIR = 0
    ELSE IF( LSAME( ID, 'I' ) ) THEN
       DIR = 1
    END IF
    IF( DIR.EQ.-1 ) THEN
       INFO = -1
    ELSE IF( N.LT.0 ) THEN
       INFO = -2
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'DLASRT', -INFO )
       RETURN
    END IF
    !
    !     Quick return if possible
    !
    IF( N.LE.1 ) RETURN
    !
    STKPNT = 1
    STACK( 1, 1 ) = 1
    STACK( 2, 1 ) = N
10  CONTINUE
    START = STACK( 1, STKPNT )
    ENDD = STACK( 2, STKPNT )
    STKPNT = STKPNT - 1
    IF( ENDD-START.LE.SELECT .AND. ENDD-START.GT.0 ) THEN
       !
       !        Do Insertion sort on D( START:ENDD )
       !
       IF( DIR.EQ.0 ) THEN
          !
          !           Sort into decreasing order
          !
          DO I = START + 1, ENDD
             DO J = I, START + 1, -1
                IF( D( J ).GT.D( J-1 ) ) THEN
                   DMNMX = D( J )
                   D( J ) = D( J-1 )
                   D( J-1 ) = DMNMX
                ELSE
                   GO TO 30
                END IF
             end DO
30           continue
          end DO
          !
       ELSE
          !
          !           Sort into increasing order
          !
          DO I = START + 1, ENDD
             DO J = I, START + 1, -1
                IF( D( J ).LT.D( J-1 ) ) THEN
                   DMNMX = D( J )
                   D( J ) = D( J-1 )
                   D( J-1 ) = DMNMX
                ELSE
                   GO TO 50
                END IF
             end DO
50           continue
          end DO
          !
       END IF
       !
    ELSE IF( ENDD-START.GT.SELECT ) THEN
       !
       !        Partition D( START:ENDD ) and stack parts, largest one first
       !
       !        Choose partition entry as median of 3
       !
       D1 = D( START )
       D2 = D( ENDD )
       I = ( START+ENDD ) / 2
       D3 = D( I )
       IF( D1.LT.D2 ) THEN
          IF( D3.LT.D1 ) THEN
             DMNMX = D1
          ELSE IF( D3.LT.D2 ) THEN
             DMNMX = D3
          ELSE
             DMNMX = D2
          END IF
       ELSE
          IF( D3.LT.D2 ) THEN
             DMNMX = D2
          ELSE IF( D3.LT.D1 ) THEN
             DMNMX = D3
          ELSE
             DMNMX = D1
          END IF
       END IF
       !
       IF( DIR.EQ.0 ) THEN
          !
          !           Sort into decreasing order
          !
          I = START - 1
          J = ENDD + 1
60        CONTINUE
70        CONTINUE
          J = J - 1
          IF( D( J ).LT.DMNMX ) GO TO 70
80        CONTINUE
          I = I + 1
          IF( D( I ).GT.DMNMX ) GO TO 80
          IF( I.LT.J ) THEN
             TMP = D( I )
             D( I ) = D( J )
             D( J ) = TMP
             GO TO 60
          END IF
          IF( J-START.GT.ENDD-J-1 ) THEN
             STKPNT = STKPNT + 1
             STACK( 1, STKPNT ) = START
             STACK( 2, STKPNT ) = J
             STKPNT = STKPNT + 1
             STACK( 1, STKPNT ) = J + 1
             STACK( 2, STKPNT ) = ENDD
          ELSE
             STKPNT = STKPNT + 1
             STACK( 1, STKPNT ) = J + 1
             STACK( 2, STKPNT ) = ENDD
             STKPNT = STKPNT + 1
             STACK( 1, STKPNT ) = START
             STACK( 2, STKPNT ) = J
          END IF
       ELSE
    !
    !           Sort into increasing order
    !
          I = START - 1
          J = ENDD + 1
90        CONTINUE
100       CONTINUE
          J = J - 1
          IF( D( J ).GT.DMNMX ) GO TO 100
110       CONTINUE
          I = I + 1
          IF( D( I ).LT.DMNMX ) GO TO 110
          IF( I.LT.J ) THEN
             TMP = D( I )
             D( I ) = D( J )
             D( J ) = TMP
             GO TO 90
          END IF
          IF( J-START.GT.ENDD-J-1 ) THEN
             STKPNT = STKPNT + 1
             STACK( 1, STKPNT ) = START
             STACK( 2, STKPNT ) = J
             STKPNT = STKPNT + 1
             STACK( 1, STKPNT ) = J + 1
             STACK( 2, STKPNT ) = ENDD
          ELSE
             STKPNT = STKPNT + 1
             STACK( 1, STKPNT ) = J + 1
             STACK( 2, STKPNT ) = ENDD
             STKPNT = STKPNT + 1
             STACK( 1, STKPNT ) = START
             STACK( 2, STKPNT ) = J
          END IF
       END IF
    END IF
    IF( STKPNT.GT.0 )  GO TO 10
    RETURN
  end SUBROUTINE QLASRT

  SUBROUTINE QLARTG( F, G, CS, SN, R )
    !
    !  -- LAPACK auxiliary routine (version 3.1) --
    !     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
    !     November 2006
    !
    !     .. Scalar Arguments ..
    real(16) CS, F, G, R, SN
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  DLARTG generate a plane rotation so that
    !
    !     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
    !     [ -SN  CS  ]     [ G ]     [ 0 ]
    !
    !  This is a slower, more accurate version of the BLAS1 routine DROTG,
    !  with the following other differences:
    !     F and G are unchanged on return.
    !     If G=0, then CS=1 and SN=0.
    !     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
    !        floating point operations (saves work in DBDSQR when
    !        there are zeros on the diagonal).
    !
    !  If F exceeds G in magnitude, CS will be positive.
    !
    !  Arguments
    !  =========
    !
    !  F       (input) DOUBLE PRECISION
    !          The first component of vector to be rotated.
    !
    !  G       (input) DOUBLE PRECISION
    !          The second component of vector to be rotated.
    !
    !  CS      (output) DOUBLE PRECISION
    !          The cosine of the rotation.
    !
    !  SN      (output) DOUBLE PRECISION
    !          The sine of the rotation.
    !
    !  R       (output) DOUBLE PRECISION
    !          The nonzero component of the rotated vector.
    !
    !  This version has a few statements commented out for thread safety
    !  (machine parameters are computed on each entry). 10 feb 03, SJH.
    !
    !  =====================================================================
    !
    !     .. Parameters ..
    real(16),parameter :: zero=0.q0
    real(16),parameter :: one=1.q0
    real(16),parameter :: two=2.q0
    !     ..
    !     .. Local Scalars ..
    !     LOGICAL            FIRST
    INTEGER COUNT, I
    real(16) EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, INT, LOG, MAX, SQRT
    !     ..
    !     .. Save statement ..
    !     SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
    !     ..
    !     .. Data statements ..
    !     DATA               FIRST / .TRUE. /
    !     ..
    !     .. Executable Statements ..
    !
    !     IF( FIRST ) THEN
    SAFMIN = tiny(1.q0) !SAFMIN = DLAMCH( 'S' )
    EPS = epsilon(1.q0)/two !EPS = DLAMCH( 'E' )
    SAFMN2 = 2.q0**INT( LOG( SAFMIN / EPS ) / LOG( 2.q0 ) / TWO ) !SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) / LOG( DLAMCH( 'B' ) ) / TWO )
    SAFMX2 = ONE / SAFMN2 ! SAFMX2 = ONE / SAFMN2
    !        FIRST = .FALSE.
    !     END IF
    IF( G.EQ.ZERO ) THEN
       CS = ONE
       SN = ZERO
       R = F
    ELSE IF( F.EQ.ZERO ) THEN
       CS = ZERO
       SN = ONE
       R = G
    ELSE
       F1 = F
       G1 = G
       SCALE = MAX( ABS( F1 ), ABS( G1 ) )
       IF( SCALE.GE.SAFMX2 ) THEN
          COUNT = 0
10        CONTINUE
          COUNT = COUNT + 1
          F1 = F1*SAFMN2
          G1 = G1*SAFMN2
          SCALE = MAX( ABS( F1 ), ABS( G1 ) )
          IF( SCALE.GE.SAFMX2 ) GO TO 10
          R = SQRT( F1**2+G1**2 )
          CS = F1 / R
          SN = G1 / R
          DO I = 1, COUNT
             R = R*SAFMX2
          end DO
       ELSE IF( SCALE.LE.SAFMN2 ) THEN
          COUNT = 0
30        CONTINUE
          COUNT = COUNT + 1
          F1 = F1*SAFMX2
          G1 = G1*SAFMX2
          SCALE = MAX( ABS( F1 ), ABS( G1 ) )
          IF( SCALE.LE.SAFMN2 ) GO TO 30
          R = SQRT( F1**2+G1**2 )
          CS = F1 / R
          SN = G1 / R
          DO I = 1, COUNT
             R = R*SAFMN2
          end DO
       ELSE
          R = SQRT( F1**2+G1**2 )
          CS = F1 / R
          SN = G1 / R
       END IF
       IF( ABS( F ).GT.ABS( G ) .AND. CS.LT.ZERO ) THEN
          CS = -CS
          SN = -SN
          R = -R
       END IF
    END IF

    return
  end SUBROUTINE QLARTG

  real(16) function qlapy2( x, y )
    !
    !  -- lapack auxiliary routine (version 3.1) --
    !     univ. of tennessee, univ. of california berkeley and nag ltd..
    !     november 2006
    !
    !     .. scalar arguments ..
    real(16) x, y
    !     ..
    !
    !  purpose
    !  =======
    !
    !  dlapy2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
    !  overflow.
    !
    !  arguments
    !  =========
    !
    !  x       (input) double precision
    !  y       (input) double precision
    !          x and y specify the values x and y.
    !
    !  =====================================================================
    !
    !     .. parameters ..
    real(16),parameter :: zero=0.q0
    real(16),parameter :: one=1.q0
    !     ..
    !     .. local scalars ..
    real(16) w, xabs, yabs, z
    !     ..
    !     .. intrinsic functions ..
    intrinsic  abs, max, min, sqrt
    !     ..
    !     .. executable statements ..
    !
    xabs = abs( x )
    yabs = abs( y )
    w = max( xabs, yabs )
    z = min( xabs, yabs )
    if( z.eq.zero ) then
       qlapy2 = w
    else
       qlapy2 = w*sqrt( one+( z / w )**2 )
    end if
    return
  end function qlapy2

  subroutine qlae2( a, b, c, rt1, rt2 )
    !
    !  -- LAPACK auxiliary routine (version 3.1) --
    !     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
    !     November 2006
    !
    !     .. Scalar Arguments ..
    real(16) A, B, C, RT1, RT2
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix
    !     [  A   B  ]
    !     [  B   C  ].
    !  On return, RT1 is the eigenvalue of larger absolute value, and RT2
    !  is the eigenvalue of smaller absolute value.
    !
    !  Arguments
    !  =========
    !
    !  A       (input) DOUBLE PRECISION
    !          The (1,1) element of the 2-by-2 matrix.
    !
    !  B       (input) DOUBLE PRECISION
    !          The (1,2) and (2,1) elements of the 2-by-2 matrix.
    !
    !  C       (input) DOUBLE PRECISION
    !          The (2,2) element of the 2-by-2 matrix.
    !
    !  RT1     (output) DOUBLE PRECISION
    !          The eigenvalue of larger absolute value.
    !
    !  RT2     (output) DOUBLE PRECISION
    !          The eigenvalue of smaller absolute value.
    !
    !  Further Details
    !  ===============
    !
    !  RT1 is accurate to a few ulps barring over/underflow.
    !
    !  RT2 may be inaccurate if there is massive cancellation in the
    !  determinant A*C-B*B; higher precision or correctly rounded or
    !  correctly truncated arithmetic would be needed to compute RT2
    !  accurately in all cases.
    !
    !  Overflow is possible only if RT1 is within a factor of 5 of overflow.
    !  Underflow is harmless if the input data is 0 or exceeds
    !     underflow_threshold / macheps.
    !
    ! =====================================================================
    !
    !     .. Parameters ..
    real(16),parameter :: zero=0.q0
    real(16),parameter :: half=0.5q0
    real(16),parameter :: one=1.q0
    real(16),parameter :: two=2.q0
    !     ..
    !     .. Local Scalars ..
    real(16) AB, ACMN, ACMX, ADF, DF, RT, SM, TB
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, SQRT
    !     ..
    !     .. Executable Statements ..
    !
    !     Compute the eigenvalues
    !
    SM = A + C
    DF = A - C
    ADF = ABS( DF )
    TB = B + B
    AB = ABS( TB )
    IF( ABS( A ).GT.ABS( C ) ) THEN
       ACMX = A
       ACMN = C
    ELSE
       ACMX = C
       ACMN = A
    END IF
    IF( ADF.GT.AB ) THEN
       RT = ADF*SQRT( ONE+( AB / ADF )**2 )
    ELSE IF( ADF.LT.AB ) THEN
       RT = AB*SQRT( ONE+( ADF / AB )**2 )
    ELSE
       !
       !        Includes case AB=ADF=0
       !
       RT = AB*SQRT( TWO )
    END IF
    IF( SM.LT.ZERO ) THEN
       RT1 = HALF*( SM-RT )
       !
       !        Order of execution important.
       !        To get fully accurate smaller eigenvalue,
       !        next line needs to be executed in higher precision.
       !
       RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
    ELSE IF( SM.GT.ZERO ) THEN
       RT1 = HALF*( SM+RT )
       !
       !        Order of execution important.
       !        To get fully accurate smaller eigenvalue,
       !        next line needs to be executed in higher precision.
       !
       RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
    ELSE
       !
       !        Includes case RT1 = RT2 = 0
       !
       RT1 = HALF*RT
       RT2 = -HALF*RT
    END IF

    return
  end subroutine qlae2

  subroutine qzlasr( side, pivot, direct, m, n, c, s, a, lda )
    !
    !  -- LAPACK auxiliary routine (version 3.1) --
    !     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
    !     November 2006
    !
    !     .. Scalar Arguments ..
    character direct, pivot, side
    integer lda, m, n
    !     ..
    !     .. Array Arguments ..
    real(16) c( * ), s( * )
    complex(16) a( lda, * )
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  ZLASR applies a sequence of real plane rotations to a complex matrix
    !  A, from either the left or the right.
    !
    !  When SIDE = 'L', the transformation takes the form
    !
    !     A := P*A
    !
    !  and when SIDE = 'R', the transformation takes the form
    !
    !     A := A*P**T
    !
    !  where P is an orthogonal matrix consisting of a sequence of z plane
    !  rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R',
    !  and P**T is the transpose of P.
    !  
    !  When DIRECT = 'F' (Forward sequence), then
    !  
    !     P = P(z-1) * ... * P(2) * P(1)
    !  
    !  and when DIRECT = 'B' (Backward sequence), then
    !  
    !     P = P(1) * P(2) * ... * P(z-1)
    !  
    !  where P(k) is a plane rotation matrix defined by the 2-by-2 rotation
    !  
    !     R(k) = (  c(k)  s(k) )
    !          = ( -s(k)  c(k) ).
    !  
    !  When PIVOT = 'V' (Variable pivot), the rotation is performed
    !  for the plane (k,k+1), i.e., P(k) has the form
    !  
    !     P(k) = (  1                                            )
    !            (       ...                                     )
    !            (              1                                )
    !            (                   c(k)  s(k)                  )
    !            (                  -s(k)  c(k)                  )
    !            (                                1              )
    !            (                                     ...       )
    !            (                                            1  )
    !  
    !  where R(k) appears as a rank-2 modification to the identity matrix in
    !  rows and columns k and k+1.
    !  
    !  When PIVOT = 'T' (Top pivot), the rotation is performed for the
    !  plane (1,k+1), so P(k) has the form
    !  
    !     P(k) = (  c(k)                    s(k)                 )
    !            (         1                                     )
    !            (              ...                              )
    !            (                     1                         )
    !            ( -s(k)                    c(k)                 )
    !            (                                 1             )
    !            (                                      ...      )
    !            (                                             1 )
    !  
    !  where R(k) appears in rows and columns 1 and k+1.
    !  
    !  Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is
    !  performed for the plane (k,z), giving P(k) the form
    !  
    !     P(k) = ( 1                                             )
    !            (      ...                                      )
    !            (             1                                 )
    !            (                  c(k)                    s(k) )
    !            (                         1                     )
    !            (                              ...              )
    !            (                                     1         )
    !            (                 -s(k)                    c(k) )
    !  
    !  where R(k) appears in rows and columns k and z.  The rotations are
    !  performed without ever forming P(k) explicitly.
    !
    !  Arguments
    !  =========
    !
    !  SIDE    (input) CHARACTER*1
    !          Specifies whether the plane rotation matrix P is applied to
    !          A on the left or the right.
    !          = 'L':  Left, compute A := P*A
    !          = 'R':  Right, compute A:= A*P**T
    !
    !  PIVOT   (input) CHARACTER*1
    !          Specifies the plane for which P(k) is a plane rotation
    !          matrix.
    !          = 'V':  Variable pivot, the plane (k,k+1)
    !          = 'T':  Top pivot, the plane (1,k+1)
    !          = 'B':  Bottom pivot, the plane (k,z)
    !
    !  DIRECT  (input) CHARACTER*1
    !          Specifies whether P is a forward or backward sequence of
    !          plane rotations.
    !          = 'F':  Forward, P = P(z-1)*...*P(2)*P(1)
    !          = 'B':  Backward, P = P(1)*P(2)*...*P(z-1)
    !
    !  M       (input) INTEGER
    !          The number of rows of the matrix A.  If m <= 1, an immediate
    !          return is effected.
    !
    !  N       (input) INTEGER
    !          The number of columns of the matrix A.  If n <= 1, an
    !          immediate return is effected.
    !
    !  C       (input) DOUBLE PRECISION array, dimension
    !                  (M-1) if SIDE = 'L'
    !                  (N-1) if SIDE = 'R'
    !          The cosines c(k) of the plane rotations.
    !
    !  S       (input) DOUBLE PRECISION array, dimension
    !                  (M-1) if SIDE = 'L'
    !                  (N-1) if SIDE = 'R'
    !          The sines s(k) of the plane rotations.  The 2-by-2 plane
    !          rotation part of the matrix P(k), R(k), has the form
    !          R(k) = (  c(k)  s(k) )
    !                 ( -s(k)  c(k) ).
    !
    !  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
    !          The M-by-N matrix A.  On exit, A is overwritten by P*A if
    !          SIDE = 'R' or by A*P**T if SIDE = 'L'.
    !
    !  LDA     (input) INTEGER
    !          The leading dimension of the array A.  LDA >= max(1,M).
    !
    !  =====================================================================
    !
    !     .. Parameters ..
    real(16),parameter :: zero=0.q0
    real(16),parameter :: one=1.q0
    !     ..
    !     .. Local Scalars ..
    integer i, info, j
    real(16) ctemp, stemp
    complex(16) temp
    !     ..
    !     .. Intrinsic Functions ..
    intrinsic max
    !     ..
    !     .. External Functions ..
    logical lsame
    external lsame
    !     ..
    !     .. External Subroutines ..
    external xerbla
    !     ..
    !     .. Executable Statements ..
    !
    !     Test the input parameters
    !
    info = 0
    if( .not.( lsame( side, 'L' ) .or. lsame( side, 'R' ) ) ) then
       info = 1
    else if( .not.( lsame( pivot, 'V' ) .or. lsame( pivot, 'T' ) .or. lsame( pivot, 'B' ) ) ) then
       info = 2
    else if( .not.( lsame( direct, 'F' ) .or. lsame( direct, 'B' ) ) ) then
       info = 3
    else if( m.lt.0 ) then
       info = 4
    else if( n.lt.0 ) then
       info = 5
    else if( lda.lt.max( 1, m ) ) then
       info = 9
    end if
    if( info.ne.0 ) then
       call xerbla( 'QZLASR' , info )
       return
    end if
    !
    !     Quick return if possible
    !
    if( ( m.eq.0 ) .or. ( n.eq.0 ) ) return
    if( lsame( side, 'L' ) ) then
       !
       !        Form  P * A
       !
       IF( LSAME( PIVOT, 'V' ) ) THEN
          IF( LSAME( DIRECT, 'F' ) ) THEN
             DO J = 1, M - 1
                CTEMP = C( J )
                STEMP = S( J )
                IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                   DO I = 1, N
                      TEMP = A( J+1, I )
                      A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                      A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
                   end DO
                END IF
             end DO
          ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
             DO J = M - 1, 1, -1
                CTEMP = C( J )
                STEMP = S( J )
                IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                   DO I = 1, N
                      TEMP = A( J+1, I )
                      A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                      A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
                   end DO
                END IF
             end DO
          END IF
       ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
          IF( LSAME( DIRECT, 'F' ) ) THEN
             DO J = 2, M
                CTEMP = C( J-1 )
                STEMP = S( J-1 )
                IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                   DO I = 1, N
                      TEMP = A( J, I )
                      A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                      A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
                   end DO
                END IF
             end DO
          ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
             DO J = M, 2, -1
                CTEMP = C( J-1 )
                STEMP = S( J-1 )
                IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                   DO I = 1, N
                      TEMP = A( J, I )
                      A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                      A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
                   end DO
                END IF
             end DO
          END IF
       ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
          IF( LSAME( DIRECT, 'F' ) ) THEN
             DO J = 1, M - 1
                CTEMP = C( J )
                STEMP = S( J )
                IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                   DO I = 1, N
                      TEMP = A( J, I )
                      A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                      A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
                   end DO
                END IF
             end DO
          ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
             DO J = M - 1, 1, -1
                CTEMP = C( J )
                STEMP = S( J )
                IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                   DO I = 1, N
                      TEMP = A( J, I )
                      A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                      A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
                   end DO
                END IF
             end DO
          END IF
       END IF
    ELSE IF( LSAME( SIDE, 'R' ) ) THEN
       !
       !        Form A * P'
       !
       IF( LSAME( PIVOT, 'V' ) ) THEN
          IF( LSAME( DIRECT, 'F' ) ) THEN
             DO J = 1, N - 1
                CTEMP = C( J )
                STEMP = S( J )
                IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                   DO I = 1, M
                      TEMP = A( I, J+1 )
                      A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                      A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
                   end DO
                END IF
             end DO
          ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
             DO J = N - 1, 1, -1
                CTEMP = C( J )
                STEMP = S( J )
                IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                   DO I = 1, M
                      TEMP = A( I, J+1 )
                      A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                      A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
                   end DO
                END IF
             end DO
          END IF
       ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
          IF( LSAME( DIRECT, 'F' ) ) THEN
             DO J = 2, N
                CTEMP = C( J-1 )
                STEMP = S( J-1 )
                IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                   DO I = 1, M
                      TEMP = A( I, J )
                      A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                      A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
                   end DO
                END IF
             end DO
          ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
             DO J = N, 2, -1
                CTEMP = C( J-1 )
                STEMP = S( J-1 )
                IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                   DO I = 1, M
                      TEMP = A( I, J )
                      A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                      A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
                   end DO
                END IF
             end DO
          END IF
       ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
          IF( LSAME( DIRECT, 'F' ) ) THEN
             DO J = 1, N - 1
                CTEMP = C( J )
                STEMP = S( J )
                IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                   DO I = 1, M
                      TEMP = A( I, J )
                      A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                      A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
                   end DO
                END IF
             end DO
          ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
             DO J = N - 1, 1, -1
                CTEMP = C( J )
                STEMP = S( J )
                IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                   DO I = 1, M
                      TEMP = A( I, J )
                      A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                      A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
                   end DO
                END IF
             end DO
          END IF
       END IF
    end if
    return
  end subroutine qzlasr

  subroutine qlaev2( a, b, c, rt1, rt2, cs1, sn1 )
    !
    !  -- LAPACK auxiliary routine (version 3.1) --
    !     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
    !     November 2006
    !
    !     .. Scalar Arguments ..
    real(16) A, B, C, CS1, RT1, RT2, SN1
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
    !     [  A   B  ]
    !     [  B   C  ].
    !  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
    !  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
    !  eigenvector for RT1, giving the decomposition
    !
    !     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
    !     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].
    !
    !  Arguments
    !  =========
    !
    !  A       (input) DOUBLE PRECISION
    !          The (1,1) element of the 2-by-2 matrix.
    !
    !  B       (input) DOUBLE PRECISION
    !          The (1,2) element and the conjugate of the (2,1) element of
    !          the 2-by-2 matrix.
    !
    !  C       (input) DOUBLE PRECISION
    !          The (2,2) element of the 2-by-2 matrix.
    !
    !  RT1     (output) DOUBLE PRECISION
    !          The eigenvalue of larger absolute value.
    !
    !  RT2     (output) DOUBLE PRECISION
    !          The eigenvalue of smaller absolute value.
    !
    !  CS1     (output) DOUBLE PRECISION
    !  SN1     (output) DOUBLE PRECISION
    !          The vector (CS1, SN1) is a unit right eigenvector for RT1.
    !
    !  Further Details
    !  ===============
    !
    !  RT1 is accurate to a few ulps barring over/underflow.
    !
    !  RT2 may be inaccurate if there is massive cancellation in the
    !  determinant A*C-B*B; higher precision or correctly rounded or
    !  correctly truncated arithmetic would be needed to compute RT2
    !  accurately in all cases.
    !
    !  CS1 and SN1 are accurate to a few ulps barring over/underflow.
    !
    !  Overflow is possible only if RT1 is within a factor of 5 of overflow.
    !  Underflow is harmless if the input data is 0 or exceeds
    !     underflow_threshold / macheps.
    !
    ! =====================================================================
    !
    !     .. Parameters ..
    real(16),parameter :: zero=0.q0
    real(16),parameter :: half=0.5q0
    real(16),parameter :: one=1.q0
    real(16),parameter :: two=2.q0
    !     ..
    !     .. Local Scalars ..
    integer sgn1, sgn2
    real(16) ab, acmn, acmx, acs, adf, cs, ct, df, rt, sm, tb, tn
    !     ..
    !     .. Intrinsic Functions ..
    intrinsic abs, sqrt
    !     ..
    !     .. Executable Statements ..
    !
    !     Compute the eigenvalues
    !
    sm = a + c
    df = a - c
    adf = abs( df )
    tb = b + b
    ab = abs( tb )
    if( abs( a ).gt.abs( c ) ) then
       acmx = a
       acmn = c
    else
       acmx = c
       acmn = a
    end if
    if( adf.gt.ab ) then
       rt = adf*sqrt( one+( ab / adf )**2 )
    else if( adf.lt.ab ) then
       rt = ab*sqrt( one+( adf / ab )**2 )
    else
       !
       !        Includes case AB=ADF=0
       !
       rt = ab*sqrt( two )
    end if
    if( sm.lt.zero ) then
       rt1 = half*( sm-rt )
       sgn1 = -1
       !
       !        Order of execution important.
       !        To get fully accurate smaller eigenvalue,
       !        next line needs to be executed in higher precision.
       !
       rt2 = ( acmx / rt1 )*acmn - ( b / rt1 )*b
    else if( sm.gt.zero ) then
       rt1 = half*( sm+rt )
       sgn1 = 1
       !
       !        Order of execution important.
       !        To get fully accurate smaller eigenvalue,
       !        next line needs to be executed in higher precision.
       !
       rt2 = ( acmx / rt1 )*acmn - ( b / rt1 )*b
    else
       !
       !        Includes case RT1 = RT2 = 0
       !
       rt1 = half*rt
       rt2 = -half*rt
       sgn1 = 1
    end if
    !
    !     Compute the eigenvector
    !
    if( df.ge.zero ) then
       cs = df + rt
       sgn2 = 1
    else
       cs = df - rt
       sgn2 = -1
    end if
    acs = abs( cs )
    if( acs.gt.ab ) then
       ct = -tb / cs
       sn1 = one / sqrt( one+ct*ct )
       cs1 = ct*sn1
    else
       if( ab.eq.zero ) then
          cs1 = one
          sn1 = zero
       else
          tn = -cs / tb
          cs1 = one / sqrt( one+tn*tn )
          sn1 = tn*cs1
       end if
    end if
    if( sgn1.eq.sgn2 ) then
       tn = cs1
       cs1 = -sn1
       sn1 = tn
    end if
    return

  end subroutine qlaev2

  subroutine qlascl( type, kl, ku, cfrom, cto, m, n, a, lda, info )
    !
    !  -- LAPACK auxiliary routine (version 3.1) --
    !     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
    !     November 2006
    !
    !     .. Scalar Arguments ..
    character type
    integer info, kl, ku, lda, m, n
    real(16) cfrom, cto
    !     ..
    !     .. Array Arguments ..
    real(16) a( lda, * )
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  DLASCL multiplies the M by N real matrix A by the real scalar
    !  CTO/CFROM.  This is done without over/underflow as long as the final
    !  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
    !  A may be full, upper triangular, lower triangular, upper Hessenberg,
    !  or banded.
    !
    !  Arguments
    !  =========
    !
    !  TYPE    (input) CHARACTER*1
    !          TYPE indices the storage type of the input matrix.
    !          = 'G':  A is a full matrix.
    !          = 'L':  A is a lower triangular matrix.
    !          = 'U':  A is an upper triangular matrix.
    !          = 'H':  A is an upper Hessenberg matrix.
    !          = 'B':  A is a symmetric band matrix with lower bandwidth KL
    !                  and upper bandwidth KU and with the only the lower
    !                  half stored.
    !          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
    !                  and upper bandwidth KU and with the only the upper
    !                  half stored.
    !          = 'Z':  A is a band matrix with lower bandwidth KL and upper
    !                  bandwidth KU.
    !
    !  KL      (input) INTEGER
    !          The lower bandwidth of A.  Referenced only if TYPE = 'B',
    !          'Q' or 'Z'.
    !
    !  KU      (input) INTEGER
    !          The upper bandwidth of A.  Referenced only if TYPE = 'B',
    !          'Q' or 'Z'.
    !
    !  CFROM   (input) DOUBLE PRECISION
    !  CTO     (input) DOUBLE PRECISION
    !          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
    !          without over/underflow if the final result CTO*A(I,J)/CFROM
    !          can be represented without over/underflow.  CFROM must be
    !          nonzero.
    !
    !  M       (input) INTEGER
    !          The number of rows of the matrix A.  M >= 0.
    !
    !  N       (input) INTEGER
    !          The number of columns of the matrix A.  N >= 0.
    !
    !  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
    !          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
    !          storage type.
    !
    !  LDA     (input) INTEGER
    !          The leading dimension of the array A.  LDA >= max(1,M).
    !
    !  INFO    (output) INTEGER
    !          0  - successful exit
    !          <0 - if INFO = -i, the i-th argument had an illegal value.
    !
    !  =====================================================================
    !
    !     .. Parameters ..
    real(16),parameter :: zero=0.q0
    real(16),parameter :: one=1.q0
    !     ..
    !     .. Local Scalars ..
    logical done
    integer i, itype, j, k1, k2, k3, k4
    real(16) bignum, cfrom1, cfromc, cto1, ctoc, mul, smlnum
    !     ..
    !     .. External Functions ..
    logical lsame
    external lsame
    !     ..
    !     .. Intrinsic Functions ..
    intrinsic abs, max, min
    !     ..
    !     .. External Subroutines ..
    external xerbla
    !     ..
    !     .. Executable Statements ..
    !
    !     Test the input arguments
    !
    info = 0
    !
    if( lsame( type, 'G' ) ) then
       itype = 0
    else if( lsame( type, 'L' ) ) then
       itype = 1
    else if( lsame( type, 'U' ) ) then
       itype = 2
    else if( lsame( type, 'H' ) ) then
       itype = 3
    else if( lsame( type, 'B' ) ) then
       itype = 4
    else if( lsame( type, 'Q' ) ) then
       itype = 5
    else if( lsame( type, 'Z' ) ) then
       itype = 6
    else
       itype = -1
    end if
    !
    if( itype.eq.-1 ) then
       info = -1
    else if( cfrom.eq.zero ) then
       info = -4
    else if( m.lt.0 ) then
       info = -6
    else if( n.lt.0 .or. ( itype.eq.4 .and. n.ne.m ) .or. &
         ( itype.eq.5 .and. n.ne.m ) ) then
       info = -7
    else if( itype.le.3 .and. lda.lt.max( 1, m ) ) then
       info = -9
    else if( itype.ge.4 ) then
       if( kl.lt.0 .or. kl.gt.max( m-1, 0 ) ) then
          info = -2
       else if( ku.lt.0 .or. ku.gt.max( n-1, 0 ) .or. &
            ( ( itype.eq.4 .or. itype.eq.5 ) .and. kl.ne.ku ) ) then
          info = -3
       else if( ( itype.eq.4 .and. lda.lt.kl+1 ) .or. &
            ( itype.eq.5 .and. lda.lt.ku+1 ) .or. &
            ( itype.eq.6 .and. lda.lt.2*kl+ku+1 ) ) then
          info = -9
       end if
    end if
    !
    if( info.ne.0 ) then
       call xerbla( 'DLASCL', -info )
       return
    end if
    !
    !     Quick return if possible
    !
    if( n.eq.0 .or. m.eq.0 ) return
    !
    !     Get machine parameters
    !
    smlnum = tiny(1.q0)
    bignum = one / smlnum
    !
    cfromc = cfrom
    ctoc = cto
    !
10  continue
    cfrom1 = cfromc*smlnum
    cto1 = ctoc / bignum
    if( abs( cfrom1 ).gt.abs( ctoc ) .and. ctoc.ne.zero ) then
       mul = smlnum
       done = .false.
       cfromc = cfrom1
    else if( abs( cto1 ).gt.abs( cfromc ) ) then
       mul = bignum
       done = .false.
       ctoc = cto1
    else
       mul = ctoc / cfromc
       done = .true.
    end if
    !
    if( itype.eq.0 ) then
       !
       !        Full matrix
       !
       do j = 1, n
          do i = 1, m
             a( i, j ) = a( i, j )*mul
          end do
       end do
       ! 
    else if( itype.eq.1 ) then
       !
       !        Lower triangular matrix
       !
       do j = 1, n
          do i = j, m
             a( i, j ) = a( i, j )*mul
          end do
       end do
       !
    else if( itype.eq.2 ) then
       !
       !        Upper triangular matrix
       !
       do j = 1, n
          do i = 1, min( j, m )
             a( i, j ) = a( i, j )*mul
          end do
       end do
       !
    else if( itype.eq.3 ) then
       !
       !        Upper Hessenberg matrix
       !
       do j = 1, n
          do i = 1, min( j+1, m )
             a( i, j ) = a( i, j )*mul
          end do
       end do
       !
    else if( itype.eq.4 ) then
       !
       !        Lower half of a symmetric band matrix
       !
       k3 = kl + 1
       k4 = n + 1
       do j = 1, n
          do i = 1, min( k3, k4-j )
             a( i, j ) = a( i, j )*mul
          end do
       end do
       !
    else if( itype.eq.5 ) then
       !
       !        Upper half of a symmetric band matrix
       !
       k1 = ku + 2
       k3 = ku + 1
       do j = 1, n
          do i = max( k1-j, 1 ), k3
             a( i, j ) = a( i, j )*mul
          end do
       end do
       !
    else if( itype.eq.6 ) then
       !
       !        Band matrix
       !
       k1 = kl + ku + 2
       k2 = kl + 1
       k3 = 2*kl + ku + 1
       k4 = kl + ku + 1 + m
       do j = 1, n
          do i = max( k1-j, k2 ), min( k3, k4-j )
             a( i, j ) = a( i, j )*mul
          end do
       end do
       !
    end if
    !
    if( .not.done ) go to 10

    return
  end subroutine qlascl

  subroutine qlassq( n, x, incx, scale, sumsq )
    !
    !  -- LAPACK auxiliary routine (version 3.1) --
    !     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
    !     November 2006
    !
    !     .. Scalar Arguments ..
    integer incx, n
    real(16) scale, sumsq
    !     ..
    !     .. Array Arguments ..
    real(16) X( * )
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  DLASSQ  returns the values  scl  and  smsq  such that
    !
    !     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
    !
    !  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
    !  assumed to be non-negative and  scl  returns the value
    !
    !     scl = max( scale, abs( x( i ) ) ).
    !
    !  scale and sumsq must be supplied in SCALE and SUMSQ and
    !  scl and smsq are overwritten on SCALE and SUMSQ respectively.
    !
    !  The routine makes only one pass through the vector x.
    !
    !  Arguments
    !  =========
    !
    !  N       (input) INTEGER
    !          The number of elements to be used from the vector X.
    !
    !  X       (input) DOUBLE PRECISION array, dimension (N)
    !          The vector for which a scaled sum of squares is computed.
    !             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
    !
    !  INCX    (input) INTEGER
    !          The increment between successive values of the vector X.
    !          INCX > 0.
    !
    !  SCALE   (input/output) DOUBLE PRECISION
    !          On entry, the value  scale  in the equation above.
    !          On exit, SCALE is overwritten with  scl , the scaling factor
    !          for the sum of squares.
    !
    !  SUMSQ   (input/output) DOUBLE PRECISION
    !          On entry, the value  sumsq  in the equation above.
    !          On exit, SUMSQ is overwritten with  smsq , the basic sum of
    !          squares from which  scl  has been factored out.
    !
    ! =====================================================================
    !
    !     .. Parameters ..
    real(16),parameter :: zero=0.q0
    !     ..
    !     .. Local Scalars ..
    integer ix
    real(16) absxi
    !     ..
    !     .. Intrinsic Functions ..
    intrinsic abs
    !     ..
    !     .. Executable Statements ..
    !
    if( n.gt.0 ) then
       do ix = 1, 1 + ( n-1 )*incx, incx
          if( x( ix ).ne.zero ) then
             absxi = abs( x( ix ) )
             if( scale.lt.absxi ) then
                sumsq = 1 + sumsq*( scale / absxi )**2
                scale = absxi
             else
                sumsq = sumsq + ( absxi / scale )**2
             end if
          end if
       end do
    end if
    return
  end subroutine qlassq


  real(16) function qlanst( norm, n, d, e )
    !
    !  -- LAPACK auxiliary routine (version 3.1) --
    !     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
    !     November 2006
    !
    !     .. Scalar Arguments ..
    character norm
    integer n
    !     ..
    !     .. Array Arguments ..
    real(16) d( * ), e( * )
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  DLANST  returns the value of the one norm,  or the Frobenius norm, or
    !  the  infinity norm,  or the  element of  largest absolute value  of a
    !  real symmetric tridiagonal matrix A.
    !
    !  Description
    !  ===========
    !
    !  DLANST returns the value
    !
    !     DLANST = ( max(abs(A(i,j))), NORM = 'M' or 'm'
    !              (
    !              ( norm1(A),         NORM = '1', 'O' or 'o'
    !              (
    !              ( normI(A),         NORM = 'I' or 'i'
    !              (
    !              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
    !
    !  where  norm1  denotes the  one norm of a matrix (maximum column sum),
    !  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
    !  normF  denotes the  Frobenius norm of a matrix (square root of sum of
    !  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
    !
    !  Arguments
    !  =========
    !
    !  NORM    (input) CHARACTER*1
    !          Specifies the value to be returned in DLANST as described
    !          above.
    !
    !  N       (input) INTEGER
    !          The order of the matrix A.  N >= 0.  When N = 0, DLANST is
    !          set to zero.
    !
    !  D       (input) DOUBLE PRECISION array, dimension (N)
    !          The diagonal elements of A.
    !
    !  E       (input) DOUBLE PRECISION array, dimension (N-1)
    !          The (n-1) sub-diagonal or super-diagonal elements of A.
    !
    !  =====================================================================
    !
    !     .. Parameters ..
    real(16),parameter :: zero=0.q0
    real(16),parameter :: one=1.q0
    !     ..
    !     .. Local Scalars ..
    integer i
    real(16) anorm, scale, sum
    !     ..
    !     .. External Functions ..
    logical lsame
    external lsame
    !     ..
    !     .. Intrinsic Functions ..
    intrinsic abs, max, sqrt
    !     ..
    !     .. Executable Statements ..
    !
    if( n.le.0 ) then
       anorm = zero
    else if( lsame( norm, 'm' ) ) then
       !
       !        Find max(abs(A(i,j))).
       !
       anorm = abs( d( n ) )
       do i = 1, n - 1
          anorm = max( anorm, abs( d( i ) ) )
          anorm = max( anorm, abs( e( i ) ) )
       end do
    else if( lsame( norm, 'O' ) .or. norm.eq.'1' .or. lsame( norm, 'I' ) ) then
       !
       !        Find norm1(A).
       !
       if( n.eq.1 ) then
          anorm = abs( d( 1 ) )
       else
          anorm = max( abs( d( 1 ) )+abs( e( 1 ) ), abs( e( n-1 ) )+abs( d( n ) ) )
          do i = 2, n - 1
             anorm = max( anorm, abs( d( i ) )+abs( e( i ) )+ abs( e( i-1 ) ) )
          end do
       end if
    else if( ( lsame( norm, 'F' ) ) .or. ( lsame( norm, 'E' ) ) ) then
       !
       !        Find normF(A).
       !
       scale = zero
       sum = one
       if( n.gt.1 ) then
          call qlassq( n-1, e, 1, scale, sum )
          sum = 2*sum
       end if
       call qlassq( n, d, 1, scale, sum )
       anorm = scale*sqrt( sum )
    end if

    qlanst = anorm

    return
  end function qlanst

  subroutine qzlaset( uplo, m, n, alpha, beta, a, lda )
    !
    !  -- LAPACK auxiliary routine (version 3.1) --
    !     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
    !     November 2006
    !
    !     .. Scalar Arguments ..
    character uplo
    integer lda, m, n
    complex(16) alpha, beta
    !     ..
    !     .. Array Arguments ..
    COMPLEX(16) A( LDA, * )
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  ZLASET initializes a 2-D array A to BETA on the diagonal and
    !  ALPHA on the offdiagonals.
    !
    !  Arguments
    !  =========
    !
    !  UPLO    (input) CHARACTER*1
    !          Specifies the part of the matrix A to be set.
    !          = 'U':      Upper triangular part is set. The lower triangle
    !                      is unchanged.
    !          = 'L':      Lower triangular part is set. The upper triangle
    !                      is unchanged.
    !          Otherwise:  All of the matrix A is set.
    !
    !  M       (input) INTEGER
    !          On entry, M specifies the number of rows of A.
    !
    !  N       (input) INTEGER
    !          On entry, N specifies the number of columns of A.
    !
    !  ALPHA   (input) COMPLEX*16
    !          All the offdiagonal array elements are set to ALPHA.
    !
    !  BETA    (input) COMPLEX*16
    !          All the diagonal array elements are set to BETA.
    !
    !  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
    !          On entry, the m by n matrix A.
    !          On exit, A(i,j) = ALPHA, 1 <= i <= m, 1 <= j <= n, i.ne.j;
    !                   A(i,i) = BETA , 1 <= i <= min(m,n)
    !
    !  LDA     (input) INTEGER
    !          The leading dimension of the array A.  LDA >= max(1,M).
    !
    !  =====================================================================
    !
    !     .. Local Scalars ..
    integer i, j
    !     ..
    !     .. External Functions ..
    logical lsame
    external lsame
    !     ..
    !     .. Intrinsic Functions ..
    intrinsic min
    !     ..
    !     .. Executable Statements ..
    !
    if( lsame( uplo, 'U' ) ) then
       !
       !        Set the diagonal to BETA and the strictly upper triangular
       !        part of the array to ALPHA.
       !
       do j = 2, n
          do i = 1, min( j-1, m )
             a( i, j ) = alpha
          end do
       end do
       do i = 1, min( n, m )
          a( i, i ) = beta
       end do
       !
    else if( lsame( uplo, 'L' ) ) then
       !
       !        Set the diagonal to BETA and the strictly lower triangular
       !        part of the array to ALPHA.
       !
       do j = 1, min( m, n )
          do i = j + 1, m
             a( i, j ) = alpha
          end do
       end do
       do i = 1, min( n, m )
          a( i, i ) = beta
       end do
       !
    else
       !
       !        Set the array to BETA on the diagonal and ALPHA on the
       !        offdiagonal.
       !
       do j = 1, n
          do i = 1, m
             a( i, j ) = alpha
          end do
       end do
       do i = 1, min( m, n )
          a( i, i ) = beta
       end do
    end if

    return
  end subroutine qzlaset

  subroutine qzsteqr( compz, n, d, e, z, ldz, work, info )
    !
    !  -- LAPACK routine (version 3.1) --
    !     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
    !     November 2006
    !
    !     .. Scalar Arguments ..
    character compz
    integer info, ldz, n
    !     ..
    !     .. Array Arguments ..
    real(16) d( * ), e( * ), work( * )
    complex(16) z( ldz, * )
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  ZSTEQR computes all eigenvalues and, optionally, eigenvectors of a
    !  symmetric tridiagonal matrix using the implicit QL or QR method.
    !  The eigenvectors of a full or band complex Hermitian matrix can also
    !  be found if ZHETRD or ZHPTRD or ZHBTRD has been used to reduce this
    !  matrix to tridiagonal form.
    !
    !  Arguments
    !  =========
    !
    !  COMPZ   (input) CHARACTER*1
    !          = 'N':  Compute eigenvalues only.
    !          = 'V':  Compute eigenvalues and eigenvectors of the original
    !                  Hermitian matrix.  On entry, Z must contain the
    !                  unitary matrix used to reduce the original matrix
    !                  to tridiagonal form.
    !          = 'I':  Compute eigenvalues and eigenvectors of the
    !                  tridiagonal matrix.  Z is initialized to the identity
    !                  matrix.
    !
    !  N       (input) INTEGER
    !          The order of the matrix.  N >= 0.
    !
    !  D       (input/output) DOUBLE PRECISION array, dimension (N)
    !          On entry, the diagonal elements of the tridiagonal matrix.
    !          On exit, if INFO = 0, the eigenvalues in ascending order.
    !
    !  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
    !          On entry, the (n-1) subdiagonal elements of the tridiagonal
    !          matrix.
    !          On exit, E has been destroyed.
    !
    !  Z       (input/output) COMPLEX*16 array, dimension (LDZ, N)
    !          On entry, if  COMPZ = 'V', then Z contains the unitary
    !          matrix used in the reduction to tridiagonal form.
    !          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the
    !          orthonormal eigenvectors of the original Hermitian matrix,
    !          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
    !          of the symmetric tridiagonal matrix.
    !          If COMPZ = 'N', then Z is not referenced.
    !
    !  LDZ     (input) INTEGER
    !          The leading dimension of the array Z.  LDZ >= 1, and if
    !          eigenvectors are desired, then  LDZ >= max(1,N).
    !
    !  WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2))
    !          If COMPZ = 'N', then WORK is not referenced.
    !
    !  INFO    (output) INTEGER
    !          = 0:  successful exit
    !          < 0:  if INFO = -i, the i-th argument had an illegal value
    !          > 0:  the algorithm has failed to find all the eigenvalues in
    !                a total of 30*N iterations; if INFO = i, then i
    !                elements of E have not converged to zero; on exit, D
    !                and E contain the elements of a symmetric tridiagonal
    !                matrix which is unitarily similar to the original
    !                matrix.
    !
    !  =====================================================================
    !
    !     .. Parameters ..
    real(16),parameter :: zero=0.q0
    real(16),parameter :: one=1.q0
    real(16),parameter :: two=2.q0
    real(16),parameter :: three=3.q0
    complex(16),parameter :: czero=cmplx(0.q0,0.q0,kind(1.q0))
    complex(16),parameter :: cone=cmplx(1.q0,0.q0,kind(1.q0))
    integer,parameter :: maxit=30
    !     ..
    !     .. Local Scalars ..
    integer i, icompz, ii, iscale, j, jtot, k, l, l1, lend, &
         lendm1, lendp1, lendsv, lm1, lsv, m, mm, mm1, nm1, nmaxit
    real(16) anorm, b, c, eps, eps2, f, g, p, r, rt1, rt2, &
         s, safmax, safmin, ssfmax, ssfmin, tst
    !     ..
    !     .. External Functions ..
    logical lsame
    external lsame
    !     ..
    !     .. External Subroutines ..
    external xerbla
    !     ..
    !     .. Intrinsic Functions ..
    intrinsic abs, max, sign, sqrt
    !     ..
    !     .. Executable Statements ..
    !
    !     Test the input parameters.
    !
    info = 0
    !
    if( lsame( compz, 'N' ) ) then
       icompz = 0
    else if( lsame( compz, 'V' ) ) then
       icompz = 1
    else if( lsame( compz, 'I' ) ) then
       icompz = 2
    else
       icompz = -1
    end if
    if( icompz.lt.0 ) then
       info = -1
    else if( n.lt.0 ) then
       info = -2
    else if( ( ldz.lt.1 ) .or. ( icompz.gt.0 .and. ldz.lt.max( 1, n ) ) ) then
       info = -6
    end if
    if( info.ne.0 ) then
       call xerbla( 'ZSTEQR', -info )
       return
    end if
    !
    !     Quick return if possible
    !
    if( n.eq.0 ) return
    !
    if( n.eq.1 ) then
       if( icompz.eq.2 ) z( 1, 1 ) = cone
       return
    end if
    !
    !     Determine the unit roundoff and over/underflow thresholds.
    !
    eps = epsilon(1.q0)/2.q0
    eps2 = eps**2
    safmin = tiny(1.q0)
    safmax = one / safmin
    ssfmax = sqrt( safmax ) / three
    ssfmin = sqrt( safmin ) / eps2
    !
    !     Compute the eigenvalues and eigenvectors of the tridiagonal
    !     matrix.
    !
    if( icompz.eq.2 ) call qzlaset( 'Full', n, n, czero, cone, z, ldz )
    !
    nmaxit = n*maxit
    jtot = 0
    !
    !     Determine where the matrix splits and choose QL or QR iteration
    !     for each block, according to whether top or bottom diagonal
    !     element is smaller.
    !
    l1 = 1
    nm1 = n - 1

10  continue
    if( l1.gt.n ) go to 160
    if( l1.gt.1 ) e( l1-1 ) = zero
    if( l1.le.nm1 ) then
       do m = l1, nm1
          tst = abs( e( m ) )
          if( tst.eq.zero ) go to 30
          if( tst.le.( sqrt( abs( d( m ) ) )*sqrt( abs( d( m+1 ) ) ) )*eps ) then
             e( m ) = zero
             go to 30
          end if
       end do
    end if
    m = n
    !
30  continue
    l = l1
    lsv = l
    lend = m
    lendsv = lend
    l1 = m + 1
    if( lend.eq.l ) go to 10
    !
    !     Scale submatrix in rows and columns L to LEND
    !
    anorm = qlanst( 'I', lend-l+1, d( l ), e( l ) )
    iscale = 0
    if( anorm.eq.zero ) go to 10
    if( anorm.gt.ssfmax ) then
       iscale = 1
       call qlascl( 'G', 0, 0, anorm, ssfmax, lend-l+1, 1, d( l ), n, info )
       call qlascl( 'G', 0, 0, anorm, ssfmax, lend-l, 1, e( l ), n, info )
    else if( anorm.lt.ssfmin ) then
       iscale = 2
       call qlascl( 'G', 0, 0, anorm, ssfmin, lend-l+1, 1, d( l ), n, info )
       call qlascl( 'G', 0, 0, anorm, ssfmin, lend-l, 1, e( l ), n, info )
    end if
    !
    !     Choose between QL and QR iteration
    !
    if( abs( d( lend ) ).lt.abs( d( l ) ) ) then
       lend = lsv
       l = lendsv
    end if
    !
    if( lend.gt.l ) then
       !
       !        QL Iteration
       !
       !        Look for small subdiagonal element.
       !
40     CONTINUE
       if( l.ne.lend ) then
          lendm1 = lend - 1
          do m = l, lendm1
             tst = abs( e( m ) )**2
             if( tst.le.( eps2*abs( d( m ) ) )*abs( d( m+1 ) ) + safmin ) go to 60
          end do
       end if
       !
       m = lend
       !
60     continue
       if( m.lt.lend ) e( m ) = zero
       p = d( l )
       IF( M.EQ.L ) GO TO 80
       !
       !        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
       !        to compute its eigensystem.
       !
       IF( M.EQ.L+1 ) THEN
          IF( ICOMPZ.GT.0 ) THEN
             call qlaev2( d( l ), e( l ), d( l+1 ), rt1, rt2, c, s )
             work( l ) = c
             work( n-1+l ) = s
             call qzlasr( 'R', 'V', 'B', n, 2, work( l ), work( n-1+l ), z( 1, l ), ldz )
          else
             call qlae2( d( l ), e( l ), d( l+1 ), rt1, rt2 )
          end if
          d( l ) = rt1
          d( l+1 ) = rt2
          e( l ) = zero
          l = l + 2
          if( l.le.lend ) go to 40
          go to 140
       end if
       !
       if( jtot.eq.nmaxit ) go to 140
       jtot = jtot + 1
       !
       !        form shift.
       !
       g = ( d( l+1 )-p ) / ( two*e( l ) )
       r = qlapy2( g, one )
       g = d( m ) - p + ( e( l ) / ( g+sign( r, g ) ) )
       !
       s = one
       c = one
       p = zero
       !
       !        Inner loop
       !
       mm1 = m - 1
       do i = mm1, l, -1
          f = s*e( i )
          b = c*e( i )
          call qlartg( g, f, c, s, r )
          if( i.ne.m-1 ) e( i+1 ) = r
          g = d( i+1 ) - p
          r = ( d( i )-g )*s + two*c*b
          p = s*r
          d( i+1 ) = g + p
          g = c*r - b
          !
          !           if eigenvectors are desired, then save rotations.
          !
          if( icompz.gt.0 ) then
             work( i ) = c
             work( n-1+i ) = -s
          end if
          !
       end do
       !
       !        if eigenvectors are desired, then apply saved rotations.
       !
       if( icompz.gt.0 ) then
          mm = m - l + 1
          call qzlasr( 'R', 'V', 'B', n, mm, work( l ), work( n-1+l ), Z( 1, L ), LDZ )
       END IF
       !
       d( l ) = d( l ) - p
       e( l ) = g
       go to 40
       !
       !        Eigenvalue found.
       !
80     continue
       d( l ) = p
       !
       l = l + 1
       if( l.le.lend ) go to 40
       go to 140
       !
    ELSE
       !
       !        QR Iteration
       !
       !        Look for small superdiagonal element.
       !
90     CONTINUE
    IF( L.NE.LEND ) THEN
       LENDP1 = LEND + 1
       DO 100 M = L, LENDP1, -1
          TST = ABS( E( M-1 ) )**2
          IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M-1 ) )+ SAFMIN )GO TO 110
100       CONTINUE
       END IF
       !
       M = LEND
       !
110    CONTINUE
       IF( M.GT.LEND ) E( M-1 ) = ZERO
       P = D( L )
       IF( M.EQ.L ) GO TO 130
       !
       !        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
       !        to compute its eigensystem.
       !
       IF( M.EQ.L-1 ) THEN
          IF( ICOMPZ.GT.0 ) THEN
             CALL QLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
             WORK( M ) = C
             WORK( N-1+M ) = S
             CALL QZLASR( 'R', 'V', 'F', N, 2, WORK( M ), WORK( N-1+M ), Z( 1, L-1 ), LDZ )
          ELSE
             CALL QLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
          END IF
          D( L-1 ) = RT1
          D( L ) = RT2
          E( L-1 ) = ZERO
          L = L - 2
          IF( L.GE.LEND ) GO TO 90
          GO TO 140
       END IF
       !
       IF( JTOT.EQ.NMAXIT ) GO TO 140
       JTOT = JTOT + 1
       !
       !        Form shift.
       !
       G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
       R = QLAPY2( G, ONE )
       G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
       !
       S = ONE
       C = ONE
       P = ZERO
       !
       !        Inner loop
       !
       LM1 = L - 1
       DO I = M, LM1
          F = S*E( I )
          B = C*E( I )
          CALL QLARTG( G, F, C, S, R )
          IF( I.NE.M ) E( I-1 ) = R
          G = D( I ) - P
          R = ( D( I+1 )-G )*S + TWO*C*B
          P = S*R
          D( I ) = G + P
          G = C*R - B
          !
          !           If eigenvectors are desired, then save rotations.
          !
          IF( ICOMPZ.GT.0 ) THEN
             WORK( I ) = C
             WORK( N-1+I ) = S
          END IF

       end DO
       !
       !        If eigenvectors are desired, then apply saved rotations.
       !
       IF( ICOMPZ.GT.0 ) THEN
          MM = L - M + 1
          CALL QZLASR( 'R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ), Z( 1, M ), LDZ )
       END IF
       !
       D( L ) = D( L ) - P
       E( LM1 ) = G
       GO TO 90
       !
       !        Eigenvalue found.
       !
130    CONTINUE
       D( L ) = P
       !
       L = L - 1
       IF( L.GE.LEND ) GO TO 90
       GO TO 140
       !
    END IF
    !
    !     Undo scaling if necessary
    !
140 CONTINUE
    IF( ISCALE.EQ.1 ) THEN
       CALL QLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1, D( LSV ), N, INFO )
       CALL QLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ), N, INFO )
    ELSE IF( ISCALE.EQ.2 ) THEN
       CALL QLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1, D( LSV ), N, INFO )
       CALL QLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ), N, INFO )
    END IF
    !
    !     Check for no convergence to an eigenvalue after a total
    !     of N*MAXIT iterations.
    !
    IF( JTOT.EQ.NMAXIT ) THEN
       DO I = 1, N - 1
          IF( E( I ).NE.ZERO ) INFO = INFO + 1
       end DO
       RETURN
    END IF
    GO TO 10
    !
    !     Order eigenvalues and eigenvectors.
    !
160 continue
    IF( ICOMPZ.EQ.0 ) THEN
       !
       !        Use Quick Sort
       !
       CALL QLASRT( 'I', N, D, INFO )
       !
    ELSE
       !
       !        Use Selection Sort to minimize swaps of eigenvectors
       !
       DO II = 2, N
          I = II - 1
          K = I
          P = D( I )
          DO J = II, N
             IF( D( J ).LT.P ) THEN
                K = J
                P = D( J )
             END IF
          end DO
          IF( K.NE.I ) THEN
             D( K ) = D( I )
             D( I ) = P
             CALL QZSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
          END IF
       end DO
    END IF
    return
  end subroutine qzsteqr

  subroutine lu_and_solve(n,a,b,ip)
    integer,intent(in) :: n
    real(16),intent(inout) :: a(n,n), b(n)
    integer,intent(out) :: ip(n)

    integer :: i, j, k, im, itmp
    real(16) :: am, w, ar, y(n), s
    
    ip=(/(i,i=1,n)/)

    do k=1,n-1
       am=abs(a(k,k))
       im=k
       do i=k+1,n
          if(abs(a(i,k)).gt.am)then
             am=abs(a(i,k))
             im=i
          end if
       end do
       if(am==0.q0) stop "aho"
       if(im.ne.k)then
          itmp=ip(k)
          ip(k)=ip(im)
          ip(im)=itmp
          do j=1,n
             w=a(k,j)
             a(k,j)=a(im,j)
             a(im,j)=w
          end do
       end if
       ar=1.q0/a(k,k)
       do i=k+1,n
          a(i,k)=a(i,k)*ar
       end do
       do j=k+1,n
          do i=k+1,n
             a(i,j)=a(i,j)-a(i,k)*a(k,j)
          end do
       end do
    end do

    y(1)=b(ip(1))
    do i=2,n
       s=0.q0
       do k=1,i-1
          s=s+a(i,k)*y(k)
       end do
       y(i)=b(ip(i))-s
    end do

    b(n)=y(n)/a(n,n)
    do i=n-1,1,-1
       s=0.q0
       do k=i+1,n
          s=s+a(i,k)*b(k)
       end do
       b(i)=(y(i)-s)/a(i,i)
    end do
  end subroutine lu_and_solve
  
end module module_qlapack
