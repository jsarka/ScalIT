!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Get the information of the nodes to send/receive data           c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine getRecvNodesInd(id, nNode, sN1, sLen1, sPos1, ePos1, bNum1, locDim1, &
                                      sN2, sLen2, sPos2, ePos2, bNum2, locDim2, &
                                      nNum, nodeInd, lenInd, locInd, gridInd)
   implicit none
   include 'mpif.h'
   integer, intent(IN) :: id, nNode
   integer, intent(IN) :: sN1, sN2, bNum1(nNode), bNum2(nNode)
   integer(KIND=MPI_OFFSET_KIND),intent(IN) :: sLen1,sPos1(nNode),ePos1(nNode)
   integer(KIND=MPI_OFFSET_KIND),intent(IN) :: sLen2,sPos2(nNode),ePos2(nNode)
   integer, intent(IN)  :: locDim1(nNode), locDim2(nNode), nNum
   integer, intent(OUT) :: nodeInd(nNum), lenInd(nNum), locInd(nNum), gridInd(nNum)

   integer :: i, j, k, ind, myPos1, myPos2
   ! below change from 2**15 limits to 4096 cores -CP Mar 2016
   !integer :: SHIFTINT = 2**12;
   integer :: SHIFTINT = 2**9;

   integer(KIND=MPI_OFFSET_KIND) :: xPos1, xPos2,  yPos1, yPos2

   ind = 0; 
  
   if (bNum2(id+1)==0) then   !  grid <- grid/seq
       
       do i = 1, sN2
          xPos1 = sPos2((id+1))+(i-1)*sLen2
          xPos2 = ePos2((id+1))+(i-1)*sLen2

          do j = 1, nNode

             if (bNum1(j)==0) then   ! grid <- grid
                yPos1 = sPos1(j)
                yPos2 = ePos1(j)+(sN1-1)*sLen1
             else                    ! grid <- seq
                yPos1 = sPos1(j)
                yPos2 = ePos1(j)
             end if

             if (yPos2 < xPos1) cycle
             if (yPos1 > xPos2) exit
             
             if (bNum1(j)==0) then   ! grid <- grid
                do k=1, sN1
                   yPos1 = sPos1(j)+(k-1)*sLen1
                   yPos2 = ePos1(j)+(k-1)*sLen1
                   if ((yPos2<xPos1) .OR.(yPos1>xPos2)) then
                       cycle
                   else
                       ind = ind + 1
                       myPos1 = max(xPos1, yPos1); myPos2 = min(xPos2, yPos2) 
                       nodeInd(ind) = j-1; gridind(ind) = k*SHIFTINT + i
                       lenInd(ind)  = myPos2 - myPos1 + 1
                       locInd(ind)  = myPos1 - xPos1 + (i-1)*locDim2(id+1) + 1 
                   end if

                end do
             else
                yPos1 = sPos1(j)
                yPos2 = ePos1(j)
                if ((yPos2<xPos1) .OR.(yPos1>xPos2)) then
                     cycle
                else
                     ind = ind + 1
                     myPos1 = max(xPos1, yPos1); myPos2 = min(xPos2, yPos2)
                     nodeInd(ind) = j-1;  gridInd(ind) = i
                     lenInd(ind)  = myPos2 - myPos1 + 1
                     locInd(ind)  = myPos1 - xPos1 + (i-1)*locDim2(id+1) + 1
                end if
             end if
          end do
       end do
   else     ! seq <- grid/seq
       xPos1 = sPos2(id+1)
       xPos2 = ePos2(id+1)

       do j = 1, nNode

          if (bNum1(j)==0) then   ! seq<-grid
             yPos1 = sPos1(j)
             yPos2 = ePos1(j)+(sN1-1)*sLen1
          else                    ! seq->seq
             yPos1 = sPos1(j)
             yPos2 = ePos1(j)
          end if

          if (yPos2 < xPos1) cycle
          if (yPos1 > xPos2) exit
             
          if (bNum1(j)==0) then   ! seq<-grid
             do k=1, sN1
                yPos1 = sPos1(j)+(k-1)*sLen1
                yPos2 = ePos1(j)+(k-1)*sLen1
                if ((yPos2<xPos1) .OR.(yPos1>xPos2)) then
                   cycle
                else
                   ind = ind + 1
                   myPos1 = max(xPos1, yPos1); myPos2 = min(xPos2, yPos2)
                   nodeInd(ind) = j-1; gridInd(ind) = k*SHIFTINT
                   lenInd(ind)  = myPos2 - myPos1 + 1
                   locInd(ind)  = myPos1 - sPos2(id+1) + 1
                end if
             end do
          else
             yPos1 = sPos1(j)
             yPos2 = ePos1(j)
             if ((yPos2<xPos1) .OR.(yPos1>xPos2)) then
                 cycle
             else
                 ind = ind + 1
                 myPos1 = max(xPos1, yPos1); myPos2 = min(xPos2, yPos2)
                 nodeInd(ind) = j-1; gridInd(ind) = 0
                 lenInd(ind)  = myPos2 - myPos1 + 1
                 locInd(ind)  = myPos1 - sPos2(id+1) + 1
             end if
          end if
          
       end do 
 
   end if

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

