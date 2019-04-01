!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Get the information of the nodes to send/receive data      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getSendNodesInd(id, nNode, sN1, sLen1, sPos1, ePos1, bNum1, locDim1, &
                                      sN2, sLen2, sPos2, ePos2, bNum2, locDim2, &
                                      nNum, nodeInd, lenInd, locInd, tagInd, gridInd)
   implicit none
   include 'mpif.h'
   integer, intent(IN) :: id, nNode
   integer, intent(IN) :: sN1, sN2, bNum1(nNode), bNum2(nNode)
   integer(KIND=MPI_OFFSET_KIND),intent(IN) :: sLen1,sPos1(nNode),ePos1(nNode)
   integer(KIND=MPI_OFFSET_KIND),intent(IN) :: sLen2,sPos2(nNode),ePos2(nNode)
   integer, intent(IN)  :: locDim1(nNode), locDim2(nNode), nNum
   integer, intent(OUT) :: nodeInd(nNum),lenInd(nNum),locInd(nNum),tagInd(nNum),gridInd(nNum)
 
   integer :: i, j, k, ind, myPos1, myPos2

   integer(KIND=MPI_OFFSET_KIND) :: xPos1, xPos2,  yPos1, yPos2
   ! below change from 2**15 limits to 4096 cores -CP Mar 2016
   !integer,parameter :: SHIFTINT = 2**12;
   integer,parameter :: SHIFTINT = 2**9;

   ind = 0
  
   if (bNum1(id+1)==0) then   !  grid-> grid/seq
       
       do i = 1, sN1
          xPos1 = sPos1(id+1)+(i-1)*sLen1
          xPos2 = ePos1(id+1)+(i-1)*sLen1

          do j = 1, nNode

             if (bNum2(j)==0) then   ! grid->grid
                yPos1 = sPos2(j)
                yPos2 = ePos2(j)+(sN2-1)*sLen2
             else                    ! grid->seq
                yPos1 = sPos2(j)
                yPos2 = ePos2(j)
             end if

             if (yPos2 < xPos1) cycle
             if (yPos1 > xPos2) exit
             
             if (bNum2(j)==0) then   ! grid->grid
                do k=1, sN2
                   yPos1 = sPos2(j)+(k-1)*sLen2
                   yPos2 = ePos2(j)+(k-1)*sLen2

                   if ((xPos2<yPos1) .OR. (xPos1>yPos2)) then
                       cycle
                   else
                       ind = ind + 1
                       myPos1 = max(xPos1, yPos1); myPos2 = min(xPos2, yPos2)
                       nodeInd(ind) = j-1;  tagInd(ind)=i*SHIFTINT+k
                       lenInd(ind)  = myPos2 - myPos1 + 1; gridInd(ind)=i
                       locInd(ind)  = myPos1 - xPos1 + (i-1)*locDim1(id+1) + 1
                   end if
                end do
             else                 ! grid->seq
                yPos1 = sPos2(j)
                yPos2 = ePos2(j)

                if ((xPos2<yPos1) .OR. (xPos1>yPos2)) then
                   cycle
                else      
                   ind = ind + 1
                   myPos1 = max(xPos1, yPos1); myPos2 = min(xPos2, yPos2)
                   nodeInd(ind) = j-1;  tagInd(ind)=i*SHIFTINT
                   lenInd(ind)  = myPos2 - myPos1 + 1; gridInd(ind)=i
                   locInd(ind)  = myPos1 - xPos1 + (i-1)*locDim1(id+1) + 1                
                end if       
             end if
          end do
       end do
   else     ! seq->grid/seq
       xPos1 = sPos1(id+1)
       xPos2 = ePos1(id+1)

       do j = 1, nNode

          if (bNum2(j)==0) then   ! seq->grid
             yPos1 = sPos2(j)
             yPos2 = ePos2(j)+(sN2-1)*sLen2
          else                    ! seq->seq
             yPos1 = sPos2(j)
             yPos2 = ePos2(j)
          end if

          if (yPos2 < xPos1) cycle
          if (yPos1 > xPos2) exit
             
          if (bNum2(j)==0) then   ! seq->grid
             do k=1, sN2
                yPos1 = sPos2(j)+(k-1)*sLen2
                yPos2 = ePos2(j)+(k-1)*sLen2
                if ((xPos2<yPos1) .OR. (xPos1>yPos2)) then
                    cycle
                else
                    ind = ind + 1
                    myPos1 = max(xPos1, yPos1); myPos2 = min(xPos2, yPos2)
                    nodeInd(ind) = j-1; tagInd(ind) = k; gridInd(ind)=0
                    lenInd(ind)  = myPos2 - myPos1 + 1
                    locInd(ind)  = myPos1 - sPos1(id+1) + 1
                end if
             end do
          else                ! seq->seq
             yPos1 = sPos2(j)
             yPos2 = ePos2(j)
             if ((xPos2<yPos1) .OR. (xPos1>yPos2)) then
                 cycle
             else
                 ind = ind + 1
                 myPos1 = max(xPos1, yPos1); myPos2 = min(xPos2, yPos2)
                 nodeInd(ind) = j-1; tagInd(ind) = 0;gridInd(ind)=0
                 lenInd(ind)  = myPos2 - myPos1 + 1
                 locInd(ind)  = myPos1 - sPos1(id+1) + 1
             end if
          end if
       end do 
 
   end if

end 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


