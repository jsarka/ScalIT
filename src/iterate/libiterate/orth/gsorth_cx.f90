!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         Complex version of Gram-Schimdit Algorithm          c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                            c
!c             Grant-Schimit Orthogonization                  c
!c                  Just for one step                         c
!c                                                            c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

LOGICAL FUNCTION GS_ORTH_CX(N, WJ, M, VJ)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: N, M
   DOUBLE COMPLEX, INTENT(IN)    :: WJ(N)
   DOUBLE COMPLEX, INTENT(INOUT) :: VJ(N,M+1)   

!ccccccccccccccccc
   DOUBLE COMPLEX   :: dotwv
   DOUBLE PRECISION :: normW
   INTEGER :: i, M0

   M0 = M + 1
   VJ(1:N, M0) = WJ(1:N)   

   DO I = 1, M
       dotwv = DOT_PRODUCT(VJ(1:N, I), WJ(1:N))
       VJ(1:N, M0) = VJ(1:N, M0) - dotwv * VJ(1:N, I)
   END DO
  
   normW = DSQRT(DBLE(DOT_PRODUCT(VJ(1:N, M0), VJ(1:N, M0))))

   IF (normW == 0.0D0) THEN
      GS_ORTH_CX = .FALSE.
   ELSE
      VJ(1:N, M0) = VJ(1:N, M0) / normW
      GS_ORTH_CX = .TRUE.
   END IF

END 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                            c
!c           Modified Grant-Schimit Orthogonization           c
!c                  Just for one step                         c
!c                                                            c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

LOGICAL FUNCTION MGS_ORTH_CX(N, WJ, M, VJ)
   IMPLICIT NONE   
   INTEGER, INTENT(IN) :: N, M
   DOUBLE COMPLEX, INTENT(IN)  :: WJ(N)
   DOUBLE COMPLEX, INTENT(INOUT) :: VJ(N,M+1)

!ccccccccccccccccc   
   DOUBLE COMPLEX   :: dotwv
   DOUBLE PRECISION :: normw
   INTEGER :: i, M0  

   M0 = M + 1
   VJ(1:N, M0) = WJ(1:N)   

   DO I = 1, M
       dotwv = DOT_PRODUCT( VJ(1:N,I),VJ(1:N, M0)) 
       VJ(1:N, M0) = VJ(1:N, M0) - dotwv * VJ(1:N, I)        
   END DO
  
   normW = DSQRT(DBLE(DOT_PRODUCT(VJ(1:N, M0), VJ(1:N, M0))))

   IF (normW /= 0.0D0) THEN      
      VJ(1:N, M0) = VJ(1:N, M0) / normW
      MGS_ORTH_CX = .TRUE.
   ELSE
      MGS_ORTH_CX = .FALSE.
   END IF

END FUNCTION MGS_ORTH_CX

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                            c
!c             Grant-Schimit Orthogonization                  c
!c                  for Lanczos Algorithm                     c
!c                                                            c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

LOGICAL FUNCTION LAN_GS_ORTH_CX(N, WJ, M, VJ, beta)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: N, M
   DOUBLE COMPLEX, INTENT(IN)  :: WJ(N)
   DOUBLE COMPLEX, INTENT(INOUT) :: VJ(N,M+1)
   double precision, intent(out)  :: beta

!ccccccccccccccccc
   DOUBLE COMPLEX   :: dotwv
   DOUBLE PRECISION :: normW
   INTEGER :: i, M0

   M0 = M + 1
   VJ(1:N, M0) = WJ(1:N)   

   DO I = 1, M
       dotwv = DOT_PRODUCT(VJ(1:N, I), WJ(1:N))
       VJ(1:N, M0) = VJ(1:N, M0) - dotwv * VJ(1:N, I)
   END DO
  
   normW = DSQRT(DBLE(DOT_PRODUCT(VJ(1:N, M0), VJ(1:N, M0))))
   beta  = normW

   IF (normW == 0.0D0) THEN
      LAN_GS_ORTH_CX = .FALSE.
   ELSE
      VJ(1:N, M0) = VJ(1:N, M0) / normW
      LAN_GS_ORTH_CX = .TRUE.
   END IF

END 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                            c
!c           Modified Grant-Schimit Orthogonization           c
!c                  Just for one step                         c
!c                                                            c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

LOGICAL FUNCTION LAN_MGS_ORTH_CX(N, WJ, M, VJ, beta)
   IMPLICIT NONE   
   INTEGER, INTENT(IN) :: N, M
   DOUBLE COMPLEX, INTENT(IN)  :: WJ(N)
   DOUBLE COMPLEX, INTENT(INOUT) :: VJ(N,M+1)
   double precision, intent(out) :: beta

!ccccccccccccccccc   
   DOUBLE COMPLEX   :: dotwv
   DOUBLE PRECISION :: normw
   INTEGER :: i, M0  

   M0 = M + 1
   VJ(1:N, M0) = WJ(1:N)   

   DO I = 1, M
       dotwv = DOT_PRODUCT( VJ(1:N,I),VJ(1:N, M0)) 
       VJ(1:N, M0) = VJ(1:N, M0) - dotwv * VJ(1:N, I)        
   END DO
  
   normW = DSQRT(DBLE(DOT_PRODUCT(VJ(1:N, M0), VJ(1:N, M0))))
   beta  = normw

   IF (normW /= 0.0D0) THEN      
      VJ(1:N, M0) = VJ(1:N, M0) / normW
      LAN_MGS_ORTH_CX = .TRUE.
   ELSE
      LAN_MGS_ORTH_CX = .FALSE.
   END IF

END 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                            c
!c             Grant-Schimit Orthogonization                  c
!c                  For the whole matrix                      c
!c                                                            c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

LOGICAL FUNCTION GS_FULL_ORTH_CX(nRow, nCol, Mat)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: nRow, nCol
   DOUBLE COMPLEX, INTENT(INOUT) :: Mat(nRow, nCol)   
                         
!ccccccccccccccccc
   DOUBLE COMPLEX     :: tmpV(nRow)
   DOUBLE PRECISION   :: normV
   LOGICAL            :: GS_ORTH_CX
   
   INTEGER :: i

   normV = DSQRT(DBLE(DOT_PRODUCT(Mat(1:nRow, 1), Mat(1:nRow, 1))))

   GS_FULL_ORTH_CX = .FALSE.
   IF (normV == 0.0D0)       RETURN

   Mat(1:nRow,1) = Mat(1:nRow, 1)/normV
   
   DO I = 2, nCol
       tmpV(1:nRow)  = Mat(1:nRow, I)
       IF ( .NOT. GS_ORTH_CX(nRow, tmpV, I-1, Mat) )    RETURN
   END DO

   GS_FULL_ORTH_CX = .TRUE.

END 
!****************************************************************


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                            c
!c       Modified Grant-Schimit Orthogonization               c
!c                  For the whole matrix                      c
!c                                                            c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

LOGICAL FUNCTION MGS_FULL_ORTH_CX(nRow, nCol, Mat)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: nRow, nCol
   DOUBLE COMPLEX, INTENT(INOUT) :: Mat(nRow, nCol)  
                         
!ccccccccccccccccc
   DOUBLE COMPLEX     :: tmpV(nRow)
   DOUBLE PRECISION   :: normV
   LOGICAL            :: MGS_ORTH_CX
   
   INTEGER :: i

   normV = DSQRT(DBLE(DOT_PRODUCT(Mat(1:nRow, 1), Mat(1:nRow, 1))))

   MGS_FULL_ORTH_CX = .FALSE.
   IF (normV == 0.0D0)      RETURN

   Mat(1:nRow,1) = Mat(1:nRow, 1)/normV

   DO I = 2, nCol       
       tmpV(1:nRow)  = Mat(1:nRow, I)
       IF ( .NOT. MGS_ORTH_CX(nRow, tmpV, I-1, Mat) )      RETURN
   END DO

   MGS_FULL_ORTH_CX = .TRUE.

END 
!****************************************************************
