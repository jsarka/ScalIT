!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Subroutines to initialize HOSB and diagonize HOSB                c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   

 subroutine HosbDiag_Seq(BJMAX, BJTOL, BLK, MAK, NDK, HOSB, VOSB, EIGVAL, ID)
      implicit none
      include 'mpif.h'
      integer,intent(IN)              :: BJMAX, Blk, MAK, NDK, ID
      double precision, intent(IN)    :: BJTOL  
      double precision,intent(INOUT)  :: HOSB(NDK,MAK,MAK,BLK)
      double precision,intent(OUT)    :: VOSB(MAK*MAK,BLK), EIGVAL(NDK,MAK,BLK)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision :: eff(5), eff0
      double precision :: ye, initye,dye
      double precision :: sumOffDiag, sumDiag, PHI, MJACOBI 

      integer :: iter, i, j, k, rootID=0
      double precision :: ct1s, ct2s, ct1a, ct2a, ct3a, ct4a!, MPI_WTime

      if (id==rootID) print *

      EFF = 0.0D0
      MAIN_LOOP : do iter = 1, BLK                       
     
         initye = sumOffDiag(NDK, MAK, HOSB(1,1,1,iter))
         ye = initye;  dye = initye;     k = 1
         call VINIT(MAK, VOSB(1, iter))   ! V=I

         if (id==rootID) ct1s = MPI_WTIME()
         !if (id==rootID) print *

         do 10 while((DABS(dye) > initye * BJTOL).and.(k <= BJMAX))
         !if (id==rootID) ct1a = MPI_WTIME()

            do i = 1,(MAK-1)   
               do j = 1,(MAK-I) 
                  PHI = MJACOBI(j, j+i,MAK,NDK, HOSB(1,1,1,iter))

                  call VUPDATE(j, j+i, PHI, MAK, VOSB(1,iter))
               end do
            end do                         
            !if (id==rootID) ct2a = MPI_WTIME()
            !if (id==rootID) print *, 'MPI WTime Seq:', K, ct2a-ct1a

            dye = sumOffDiag(NDK, MAK, HOSB(1,1,1,iter)) - ye
            ye = ye + dye
            
            k = k + 1            
         
 10      end do  

         if (id==rootID) ct2s = MPI_WTIME()
         !if (id==rootID) print *
         if ((id==rootID).and.((iter.eq.1).or.(iter.eq.BLK))) print *, '   HOSB Seq Sum:',ITER, BLK, K-1, ct2s-ct1s
         !if (id==rootID) print *

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         eff0 = abs(1.0D0 - ye/initye)
         if (iter==1) then
             EFF(1) = eff0;     EFF(2) = eff0
         else
             if (eff0 < EFF(1) ) EFF(1) = eff0
             if (eff0 > EFF(2) ) EFF(2) = eff0
         end if

         EFF(3) = EFF(3) + eff0
         EFF(4) = EFF(4) + sumDiag(NDK, MAK, HOSB(1,1,1,iter))
         EFF(5) = EFF(5) + k
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         do i=1, MAK            
            EIGVAL(1:NDK,i,iter) = HOSB(1:NDK, i, i, iter)
         end do                    

         do i=1, MAK                        
            HOSB(1:NDK, i, i, iter) = 0.0D0
         end do
      end do MAIN_LOOP            
     
!      if(ID==rootID) then
!          print *
!          write(*, 20) EFF(5)/BLK,  EFF(3)/EFF(5)
!          write(*, 30) EFF(1), EFF(2), EFF(3)/BLK
!          print *, '-------------------------------------------------------------------'
!      end if

20        format('  Aver. # of Sweeps :', F10.5, ',   Aver. OffDiag/Diag :',F10.6)
30        format('  Jacobi Effi: Min :', F10.5, ',   Max :', F10.5,',  Mean :',F10.5)
 end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
