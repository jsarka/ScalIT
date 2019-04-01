!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             QMR Function in MPI Environment                      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 integer function QMR_MPI(nTYPE, QMRMAX, QMRTOL, N,  &
                 B, HX_MPI, PX_MPI, X, RES, myid)
      implicit none
      include 'mpif.h'
      integer, intent(IN) ::  nTYPE, QMRMAX , N
      double precision, intent(IN) ::QMRTOL
      external :: HX_MPI, PX_MPI
      double precision, intent(IN)  :: B(N)
      double precision, intent(OUT) :: X(N) 
      double precision, intent(OUT) :: RES    
      integer, intent(IN)   :: myid

!cccccccccccccc Temporary vectors  
      integer :: IERR
      double precision :: R(N), UX(N), UR(N) 
      double precision :: Q(N), S(N), V(N), Z(N)    

      double precision :: NORM_MPI, DOTPROD_MPI 
      double precision :: BNRM,ZNRM,VNRM,VNRMPO,ZVRAT               

      double precision :: QSDOT, VZDOT, BETA, ETA     
      double precision :: GAMAMO, GAMA, THETAMO, THETA  

      integer   :: i

      double precision :: RES1, RES2, RES3, XNRM  
      double precision,dimension(:),allocatable :: HX0                 

      integer :: rootID = 0
      double precision :: wtp1, wtp2, wtp3, wtp4
      double precision :: wtp5, wtp6, wtp7, wtp8
      double precision :: wtp9, wtp10, wtp11, wtp12
!
!      INTEGER   :: MYID          ! my id in MPI
!      
!      CALL MPI_COMM_RANK(MPI_COMM_WORLD, myID, IERR)
!
!cccccccccccccccccccccccccccccccccccccccccccc
      if (ABS(nTYPE) == 1) then    ! using |HX-B|/|HX| for conv.
          allocate(HX0(N), STAT=QMR_MPI)
          if (QMR_MPI /= 0) then
              QMR_MPI = -10;   return
          end if
      end if

!ccccccccccccccccccccccccccccccccccccccccc
!c     Initialization: R = V = B         c
!ccccccccccccccccccccccccccccccccccccccccc
      ZVRAT  = 1.0d0
      if (nType <= 0)   X(1:N) = 0.0D0   

      R(1:N) = B(1:N);  V(1:N) = R(1:N)

!cccccccccccccccccccccccccccccccccccccccccccccccccc
!c            BNRM = NORM(B,N) B = V              c
!c   It is curious that using NORM_MPI(MPI_COMM_  c
!c    WORLD) gives SIG 2/11 error                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccc
      BNRM = NORM_MPI(MPI_COMM_WORLD, N, V, IERR)
      VNRM = BNRM
    
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Solve Z = P^(-1) B                                   c
!c     This is the initial multiplication of vector B by    c
!c             inverse preconditioner.                      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call PX_MPI(N, B, Z, IERR)    

      ZNRM = NORM_MPI(MPI_COMM_WORLD, N, Z, IERR)

      GAMAMO = 1.0D0;     ETA = -1.0D0

!ccccccccccccccccccccccccccccc      
!c     Main Loop             c
!ccccccccccccccccccccccccccccc
      QMR_MPI = 0
      MAIN_LOOP: do i = 1,QMRMAX          

         if ( (VNRM == 0.0d0) .or. (ZNRM == 0.0d0)) then  
            QMR_MPI = -1;   exit   
         end if

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  BP addition to reduce number of vectors and matrix-vector products.   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ZVRAT = ZVRAT * ZNRM/VNRM

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        Normalize V and Z, compute VZDOT              c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         V(1:N) = V(1:N) / VNRM ;  Z(1:N) = Z(1:N) / ZNRM  
         VZDOT = DOTPROD_MPI(MPI_COMM_WORLD, N, V, Z, IERR)
         if ( VZDOT == 0.0d0 )  then
            QMR_MPI = -2;      exit
         end if

!          if (myid==0) then
!             print *, 'VZDOT=',VZDOT,' V(1:10)=',V(1:10), &
!                      ' Z(1:10)=',Z(1:10),' N=',N
!             print *, 'i=',i
!             print *
!          end if
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

         !if (i.eq.5) then
           wtp1 = MPI_WTime()
           call HX_MPI(N, Q, S, IERR)
           wtp2 = MPI_WTime()
           !if (myid==rootID) print *, ' QMR HX :',wtp2-wtp1
         !else
         !  call HX_MPI(N, Q, S, IERR)
         !end if

         !call HX_MPI(N, Q, S, IERR)         
         S(1:N) = ZVRAT * S(1:N) 

!cccccccccccccccccccccccccccccccccccccccccccc
!c        Update QSDOT and BETA             c
!cccccccccccccccccccccccccccccccccccccccccccc

         QSDOT = DOTPROD_MPI(MPI_COMM_WORLD, N, Q, S, IERR)
         BETA = QSDOT/VZDOT

         if ( (QSDOT == 0.0d0) .or. (BETA == 0.0d0))  then  
            !print *, 'VZDOT=',VZDOT,' QSDOT=',QSDOT,' BETA=',BETA  
            QMR_MPI = -3;      exit    
         end if

!ccccccccccccccccccccccccccccccccccccc         
!c        Update V, VNRMPO           c
!ccccccccccccccccccccccccccccccccccccc
          V(1:N) = S(1:N) - BETA * V(1:N)  
         VNRMPO = NORM_MPI(MPI_COMM_WORLD, N, V, IERR)  

!ccccccccccccccccccccccccccccccccccccccccccccccccc
!c        Update Z = P^(-1) V/ZVRAT              c
!ccccccccccccccccccccccccccccccccccccccccccccccccc
         !if (i.eq.5) then
           wtp3 = MPI_WTime()
           call PX_MPI(N, V, Z, IERR)
           wtp4 = MPI_WTime()
           !if (myid==rootID) print *, ' QMR PX :',wtp4-wtp3
         !else
         !  call PX_MPI(N, V, Z, IERR)
         !end if

         !call PX_MPI(N, V, Z, IERR)  
         Z(1:N) = Z(1:N)/ZVRAT    

!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        Update ZNRM, THETA, GAMA, ETA            c
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
      
         ZNRM  = NORM_MPI(MPI_COMM_WORLD, N, Z, IERR)  

         THETA = VNRMPO/(GAMAMO*DABS(BETA))
         GAMA = 1.0d0/DSQRT(1.0d0+THETA*THETA)
         ETA = -ETA*VNRM*GAMA*GAMA/(BETA*GAMAMO*GAMAMO)

         if (GAMA == 0.0d0) then
            QMR_MPI = -4;         exit   
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
            !if (i.eq.5) then
              wtp5 = MPI_WTime()

                call HX_MPI(N, X, HX0, IERR)    ! HX0 = H * X
                XNRM = NORM_MPI(MPI_COMM_WORLD, N, HX0, IERR)
                HX0(1:N) = HX0(1:N) - B(1:N)
                RES = NORM_MPI(MPI_COMM_WORLD, N, HX0, IERR)/XNRM

              wtp6 = MPI_WTime()
              !if (myid==rootID) print *, ' QMR CHK:',wtp6-wtp5
            !else
            !    call HX_MPI(N, X, HX0, IERR)    ! HX0 = H * X
            !    XNRM = NORM_MPI(MPI_COMM_WORLD, N, HX0, IERR)
            !    HX0(1:N) = HX0(1:N) - B(1:N)
            !    RES = NORM_MPI(MPI_COMM_WORLD, N, HX0, IERR)/XNRM
            !end if
                !call HX_MPI(N, X, HX0, IERR)    ! HX0 = H * X          
                !XNRM = NORM_MPI(MPI_COMM_WORLD, N, HX0, IERR)
                !HX0(1:N) = HX0(1:N) - B(1:N)
                !RES = NORM_MPI(MPI_COMM_WORLD, N, HX0, IERR)/XNRM


         case (2)        ! Full
               RES1 = NORM_MPI(MPI_COMM_WORLD, N, R, IERR)/BNRM    ! |R|/BNRM     
               RES2 = NORM_MPI(MPI_COMM_WORLD, N, UX, IERR)/        &
                      NORM_MPI(MPI_COMM_WORLD, N, X, IERR)         ! |UX|/|X|
               RES3 = NORM_MPI(MPI_COMM_WORLD, N, UR, IERR)/        &
                      NORM_MPI(MPI_COMM_WORLD, N, R,  IERR)         ! |UR|/|R|       
               RES = max(RES1, RES2, RES3)

         case DEFAULT       ! basic convergence criterion 
               RES = NORM_MPI(MPI_COMM_WORLD, N, R, IERR)/BNRM  
                      ! RES = NORM(R,N)/BNRM

         end select        

!         IF (myid==0) THEN
!            PRINT *, 'Iter:', i, 'Res:', RES
!         END IF

         if (DABS(RES) < QMRTOL) then 
            !if (myid==rootID) print *, ' SUM HX :',(wtp2-wtp1)*i
            !if (myid==rootID) print *, ' SUM PX :',(wtp4-wtp3)*i
            !if (myid==rootID) print *, ' SUM CHK:',(wtp6-wtp5)*i
            if (myid==rootID) print *, ' QMR HX :',wtp2-wtp1,(wtp2-wtp1)*i
            if (myid==rootID) print *, ' QMR PX :',wtp4-wtp3,(wtp4-wtp3)*i
            if (myid==rootID) print *, ' QMR CHK:',wtp6-wtp5,(wtp6-wtp5)*i
            QMR_MPI = i;      exit  
         else
            QMR_MPI = -i
         end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        Update variables on two levels                   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         VNRM = VNRMPO;    GAMAMO = GAMA
         THETAMO = THETA

      end do  MAIN_LOOP

      if (allocated(HX0))   deallocate(HX0)

 end 
!***************************************************************
