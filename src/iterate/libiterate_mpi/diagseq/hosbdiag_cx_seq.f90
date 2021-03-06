!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Subroutines to initialize HOSB and diagonize HOSB                c
!c                   HOSB is a Hermitian matrix                           c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   

 subroutine HosbDiagCX_Seq(BJMAX, BJTOL, BLK, MAK, NDK, HOSB, VOSB, EIGVAL)
      implicit none
      integer,intent(IN)           :: BJMAX, Blk, MAK, NDK           
      double precision, intent(IN) :: BJTOL  
      double complex,intent(INOUT) :: HOSB(NDK,MAK,MAK,BLK)
      double precision,intent(OUT) :: VOSB(MAK*MAK,BLK)
      double complex,intent(OUT)   :: EIGVAL(NDK,MAK,BLK)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision :: eff(5), eff0
      double precision :: ye, initye,dye
      double precision :: sumOffDiagCX, sumDiagCX, PHI, MJACOBI_CX 

      integer :: iter, i, j, k, rootID=0       

      eff = 0.0D0

      MAIN_LOOP : do iter = 1, BLK                       
     
         initye = sumOffDiagCX(NDK, MAK, HOSB(1,1,1,iter))
         ye = initye;  dye = initye;     k = 1

         call VINIT(MAK, VOSB(1, iter))   ! V=I

         do 10 while((DABS(dye) > initye * BJTOL).and.(k <= BJMAX))

            do i = 1,(MAK-1)   
               do j = 1,(MAK-I) 
                  PHI = MJACOBI_CX(j, j+i,MAK,NDK, HOSB(1,1,1,iter))

                  call VUPDATE(j, j+i, PHI, MAK, VOSB(1,iter))
               end do
            end do                         

            dye = sumOffDiagCX(NDK, MAK, HOSB(1,1,1,iter)) - ye
            ye = ye + dye
            
            k = k + 1            

 10      end do  

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         eff0 = abs(1.0D0 - ye/initye)
         if (iter==1) then
             eff(1) = eff0;     eff(2) = eff0
         else
             if (eff0 < eff(1) ) eff(1) = eff0
             if (eff0 > eff(2) ) eff(2) = eff0
         end if

         eff(3) = eff(3) + eff0
         eff(4) = eff(4) + sumDiagCX(NDK, MAK, HOSB(1,1,1,iter))
         eff(5) = eff(5) + k
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         do i=1, MAK            
            EIGVAL(1:NDK,i,iter) = HOSB(1:NDK, i, i, iter)
         end do                    

         do i=1, MAK                        
            HOSB(1:NDK, i, i, iter) = 0.0D0
         end do
      end do MAIN_LOOP            

      write(*, 20) eff(5)/BLK,  eff(3)/eff(5)
      write(*, 30) eff(1), eff(2), eff(3)/BLK

20    format(' # of Iter.:', F10.5, '.   OffDiag/Diag:',F10.6)
30    format(' Jacobi Effi.: Min:',F10.5,',   Max:',F10.5,'.  Mean:',F10.5)

 end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
