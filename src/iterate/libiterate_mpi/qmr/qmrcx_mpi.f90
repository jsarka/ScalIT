!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             QMR Function in MPI Environment                      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 integer function QMR_CX_MPI(nTYPE, QMRMAX, QMRTOL, N, &
                 B, HX_MPI, PX_MPI, X, RES)
      implicit none
      include 'mpif.h'
      integer, intent(IN)          :: nTYPE, QMRMAX, N
      double precision, intent(IN) :: QMRTOL
      external :: HX_MPI, PX_MPI
      double complex,intent(IN)  :: B(N)
      double complex,intent(OUT) :: X(N) 
      double precision,intent(OUT)             :: RES    


!cccccccccccccc Temporary vectors    
      integer :: IERR  
 
      double complex :: R(N), UX(N), UR(N) 
      double complex :: Q(N), S(N),V(N), Z(N)   

      double precision  ::  NORM_CX_MPI
      double complex    ::  DOT_CX_MPI 
      double precision  ::  BNRM, ZNRM,VNRM, VNRMPO,ZVRAT         
      double complex    ::  QSDOT, VZDOT, BETA, ETA    
      double precision  ::  GAMAMO, GAMA, THETAMO, THETA  

      integer   :: i

      double precision :: RES1, RES2, RES3, XNRM  
      double complex,dimension(:),allocatable :: HX0  

!      INTEGER   :: MYID          ! my id in MPI
!
!*****************************************************
!      CALL MPI_COMM_RANK(MPI_COMM_WORLD, myID, IERR)

!cccccccccccccccccccccccccccccccccccccccccccc
      if (ABS(nTYPE) == 1) then    ! using |HX-B|/|HX| for the convergence
          allocate(HX0(N), STAT=QMR_CX_MPI)
          if (QMR_CX_MPI /= 0) then
              QMR_CX_MPI = -10
              return
          end if
      end if


!ccccccccccccccccccccccccccccccccccccccccc
!c     Initialization: R = V = B         c
!ccccccccccccccccccccccccccccccccccccccccc
      ZVRAT = 1.0d0
      if (nType<=0)  X(1:N) = 0.0D0 

      R (1:N)= B(1:N);    V(1:N) = R(1:N)
      BNRM = NORM_CX_MPI(MPI_COMM_WORLD, N, V, IERR)
      VNRM = BNRM    

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Solve Z = P^(-1) B                                   c
!c     This is the initial multiplication of vector B by    c
!c             inverse preconditioner.                      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call PX_MPI(N, B, Z, IERR)     
      ZNRM = NORM_CX_MPI(MPI_COMM_WORLD, N, Z, IERR)
      GAMAMO = 1.0D0;     ETA = -1.0D0

!ccccccccccccccccccccccccccccc      
!c     Main Loop             c
!ccccccccccccccccccccccccccccc
      QMR_CX_MPI = 0
      MAIN_LOOP: do i = 1,QMRMAX          

         if ( (VNRM == 0.0d0) .or. (ZNRM == 0.0d0)) then  
            QMR_CX_MPI = -1;      exit  
         end if

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  BP addition to reduce number of vectors and matrix-vector products.   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ZVRAT = ZVRAT * ZNRM/VNRM

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        Normalize V and Z, compute VZDOT              c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         V(1:N) = V(1:N) / VNRM  ;   Z(1:N) = Z(1:N) / ZNRM  

         VZDOT = DOT_CX_MPI(MPI_COMM_WORLD, N, V, Z, IERR)

         if ( abs(VZDOT)==0.0d0 )  then
            QMR_CX_MPI = -2;      exit
         end if

!cccccccccccccccccccccccccccc
!c        Update Q          c
!cccccccccccccccccccccccccccc
         if ( i == 1) then
            Q(1:N) = Z(1:N)
         else
           Q(1:N) = Z(1:N) - (VNRM*VZDOT/QSDOT) * Q(1:N) 
         end if
         
         
!cccccccccccccccccccccccccccccccccccccc         
!c        Update S = ZVRAT H Q        c
!cccccccccccccccccccccccccccccccccccccc

         call HX_MPI(N, Q, S, IERR)        
         S(1:N) = ZVRAT * S(1:N) 

!cccccccccccccccccccccccccccccccccccccccccccc
!c        Update QSDOT and BETA             c
!cccccccccccccccccccccccccccccccccccccccccccc

         QSDOT = DOT_CX_MPI(MPI_COMM_WORLD, N, Q, S, IERR)

         BETA = QSDOT/VZDOT

         if ( (abs(QSDOT)==0.0d0) .or. (abs(BETA)==0.0d0))  then    
            QMR_CX_MPI = -3;        exit   
         end if

!ccccccccccccccccccccccccccccccccccccc         
!c        Update V, VNRMPO           c
!ccccccccccccccccccccccccccccccccccccc
         V(1:N) = S(1:N) - BETA * V(1:N) 
         VNRMPO = NORM_CX_MPI(MPI_COMM_WORLD, N, V, IERR)  

!ccccccccccccccccccccccccccccccccccccccccccccccccc
!c        Update Z = P^(-1) V/ZVRAT              c
!ccccccccccccccccccccccccccccccccccccccccccccccccc

         call PX_MPI(N, V, Z, IERR)  
         Z(1:N) = Z(1:N)/ZVRAT  

!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        Update ZNRM, THETA, GAMA, ETA            c
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
      
         ZNRM  = NORM_CX_MPI(MPI_COMM_WORLD, N, Z, IERR)  

         THETA = VNRMPO/(GAMAMO*abs(BETA))
         GAMA = 1.0d0/DSQRT(1.0d0+THETA*THETA)
         ETA = -ETA*VNRM*GAMA*GAMA/(BETA*GAMAMO*GAMAMO)

         if (GAMA == 0.0d0) then
            QMR_CX_MPI = -4;      exit 
         end if

!cccccccccccccccccccccccccccccccccccccc         
!c        Update X and R              c
!cccccccccccccccccccccccccccccccccccccc
         if (i == 1) then
            UX(1:N) = (ETA*ZVRAT)*Q(1:N) 
            UR(1:N) = ETA * S(1:N)      
         else
            UX(1:N) = (ETA*ZVRAT)*Q(1:N) +((THETAMO*GAMA)**2)*UX(1:N)
            UR(1:N) =  ETA * S(1:N) + ((THETAMO*GAMA)**2)*UR(1:N)
         end if

         X(1:N) = X(1:N) + UX(1:N)  
         R(1:N) = R(1:N) - UR(1:N)  

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!c        Relative residual processing                  c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ! convergence criterion
         select case (ABS(nTYPE))   
         
         case (1)        ! APR, Not implemented
                call HX_MPI(N, X, HX0, IERR)    ! HX0 = H * X          
                XNRM = NORM_CX_MPI(MPI_COMM_WORLD, N, HX0,IERR)
                HX0(1:N) = HX0(1:N) - B(1:N)
                RES = NORM_CX_MPI(MPI_COMM_WORLD, N, HX0, IERR)/BNRM

         case (2)       ! Full
               RES1 = NORM_CX_MPI(MPI_COMM_WORLD, N, R, IERR)/BNRM  !|R|/BNRM     
               RES2 = NORM_CX_MPI(MPI_COMM_WORLD, N, UX,IERR)/   &
                      NORM_CX_MPI(MPI_COMM_WORLD, N, X, IERR)       !|UX|/|X|
               RES3 = NORM_CX_MPI(MPI_COMM_WORLD, N, UR,IERR)/   &
                      NORM_CX_MPI(MPI_COMM_WORLD, N, R, IERR)       !|UR|/|R|       
               RES = max(RES1, RES2, RES3)

         case DEFAULT       ! basic convergence criterion 
               RES = NORM_CX_MPI(MPI_COMM_WORLD, N, R, IERR)/BNRM  

         end select

!        CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
!         if (myid==0) then
!            print *, 'Iter:', i, 'Res:', RES
!            print *, 'X',X
!         end if

         if (DABS(RES) < QMRTOL) then  !  the  conv.
            QMR_CX_MPI = i;      exit  
         else
            QMR_CX_MPI = -i
         end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        Update variables on two levels                   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         VNRM = VNRMPO;     GAMAMO = GAMA
         THETAMO = THETA

      end do  MAIN_LOOP

      if (allocated(HX0))  deallocate(HX0)

  end 
!**************************************************************


