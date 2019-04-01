!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                         QMR Function                             c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function QMR(nTYPE, QMRMAX, QMRTOL, N, B, HX, PX, X, RES)
      implicit none
      integer, intent(IN)  :: nTYPE, QMRMAX, N
      double precision,intent(IN) :: QMRTOL  
      double precision,intent(IN)  :: B(N)         
      external :: HX, PX
      double precision, intent(OUT) :: X(N) 
      double precision, intent(OUT) :: RES    

!cccccccccccccc Temporary vectors  cccccccccccccccccccccccccccccc
      double precision :: R(N), UX(N), UR(N) 
      double precision :: Q(N), S(N), V(N), Z(N)   
      double precision :: BNRM, ZNRM, VNRM, VNRMPO, ZVRAT 
      double precision :: QSDOT, VZDOT, BETA, ETA    
      double precision :: GAMAMO, GAMA, THETAMO, THETA       

      integer ::  i, info 
      double precision :: RES1, RES2, RES3, XNRM  
      double precision,allocatable :: HX0 (:)                

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      QMR = -1
      if (ABS(nTYPE) == 1) then       ! use |HX-B|/|HX| for conv.
          allocate(HX0(N), STAT=info)
          if (info /= 0) return          
      end if

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Initialization: R=V=B:  nType < 0, use 0.0 as initial guess c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    

      ZVRAT = 1.0d0

      if (nType <= 0)  X(1:N) = 0.0D0

      R(1:N) = B(1:N) ;   V(1:N) = B(1:N) 
      BNRM   = DSQRT(dot_product(B, B))   
      VNRM   = BNRM

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    Z = P * B ~ H^-1 * B, the preconditioner         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call PX(N, B, Z)  
      ZNRM = DSQRT(dot_product(Z,Z))  
      GAMAMO = 1.0D0;       ETA = -1.0D0

!cccccccccccccccccccccccccccccccccccccccc      
!c     Main Loop for QMR iteration      c
!cccccccccccccccccccccccccccccccccccccccc

      QMR = 0
      MAIN_LOOP: do i = 1, QMRMAX          

         if ( (VNRM == 0.0d0) .OR. (ZNRM == 0.0d0))  exit

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  BP addition to reduce number of vectors and matrix-vector products.  c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         ZVRAT = ZVRAT * ZNRM/VNRM

!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        Normalize V and Z, compute VZDOT          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccc 
         V = V / VNRM ;    Z = Z / ZNRM  
         VZDOT = dot_product(V, Z)  
         if (VZDOT == 0.0D0)  exit

!cccccccccccccccccccccccccccccccccccccc
!c             Update Q               c
!cccccccccccccccccccccccccccccccccccccc

         if ( i == 1) then
            Q(1:N) = Z(1:N)  
         else
            Q(1:N) = Z(1:N) - (VNRM*VZDOT/QSDOT) * Q(1:N) 
         end if

!cccccccccccccccccccccccccccccccccccccc         
!c        Update S = ZVRAT H Q        c
!cccccccccccccccccccccccccccccccccccccc
         call HX(N, Q, S)  
         S(1:N) = ZVRAT * S(1:N)  

!cccccccccccccccccccccccccccccccccccccc
!c      Update QSDOT and BETA         c
!cccccccccccccccccccccccccccccccccccccc
         QSDOT = dot_product(Q, S)  
         BETA = QSDOT/VZDOT

         if ((QSDOT == 0.0D0) .OR. (BETA == 0.0d0)) exit

!ccccccccccccccccccccccccccccccccccccccc         
!c        Update V, VNRMPO             c
!ccccccccccccccccccccccccccccccccccccccc
         V(1:N) = S(1:N) - BETA * V(1:N)  
         VNRMPO = DSQRT(dot_product(V,V))  

!ccccccccccccccccccccccccccccccccccccccc
!c      Update Z = P^(-1) V/ZVRAT      c
!ccccccccccccccccccccccccccccccccccccccc
         call PX(N, V, Z)        
         Z(1:N) = Z(1:N)/ZVRAT   

!ccccccccccccccccccccccccccccccccccccccc
!c   Update ZNRM, THETA, GAMA, ETA     c
!ccccccccccccccccccccccccccccccccccccccc
         ZNRM  = DSQRT(dot_product(Z,Z))  
         THETA = VNRMPO/(GAMAMO*DABS(BETA))
         GAMA  = 1.0d0/DSQRT(1.0d0+THETA*THETA)
         ETA   = -ETA*VNRM*GAMA*GAMA/(BETA*GAMAMO*GAMAMO)

         if (GAMA == 0.0D0) exit   

!ccccccccccccccccccccccccccccccccccccccc         
!c        Update X and R               c
!ccccccccccccccccccccccccccccccccccccccc
         if (i == 1) then
            UX(1:N) = (ETA*ZVRAT)*Q(1:N) 
            UR(1:N) = ETA * S(1:N)       
         else
            UX(1:N) = (ETA*ZVRAT)*Q(1:N) +((THETAMO*GAMA)**2)*UX(1:N)
            UR(1:N) =  ETA * S(1:N) + ((THETAMO*GAMA)**2)*UR(1:N)
         end if

         X(1:N) = X(1:N) + UX(1:N)  
         R(1:N) = R(1:N) - UR(1:N)  

!cccccccccccccccccccccccccccccccccccccc
!c     Convergence Testing            c
!cccccccccccccccccccccccccccccccccccccc
         select case (ABS(nTYPE))             
         case (1)       ! B-H*X
               call HX(N, X, HX0)        ! HX0 = H * X          
               XNRM = DSQRT(dot_product(HX0, HX0))
               HX0(1:N) = HX0(1:N) - B(1:N)
               RES = DSQRT(dot_product(HX0, HX0))/XNRM

         case (2)       ! Full: |R|/|B|, |UX|/X, |UR|/|R|
               RES1 = DSQRT(dot_product(R, R))/BNRM                    
               RES2 = DSQRT(dot_product(UX,UX)/dot_product(X, X))  
               RES3 = DSQRT(dot_product(UR,UR)/dot_product(R, R))         
               RES = max(RES1, RES2, RES3)

         case DEFAULT   ! Simple: |R|/|B| 
               RES = DSQRT(dot_product(R, R))/BNRM   

         end select

         if (i/200*200 ==i)    write(*,10) i, RES

         if ( abs(RES) < QMRTOL) then
             QMR = i;  exit
         else
             QMR = -i
         end if         

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        Update variables on two levels                   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         VNRM = VNRMPO;          GAMAMO = GAMA
         THETAMO = THETA

      end do  MAIN_LOOP

      if (allocated(HX0))  deallocate(HX0)

 10   format('        QMR Iter:',I8,3x,'RES:',2X,E20.15)

end function QMR
!****************************************************************
