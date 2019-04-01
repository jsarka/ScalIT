!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Get the information of the nodes to send/receive data      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function getSendNodesNum(id, nNode, sN1, sLen1, sPos1,     &
          ePos1, bNum1, sN2, sLen2, sPos2, ePos2, bNum2)
   implicit none
   include 'mpif.h'
   integer, intent(IN) :: id, nNode
   integer, intent(IN) :: sN1, sN2, bNum1(nNode),bNum2(nNode)
   integer(KIND=MPI_OFFSET_KIND),intent(IN) :: sLen1,sPos1(nNode),ePos1(nNode)
   integer(KIND=MPI_OFFSET_KIND),intent(IN) :: sLen2,sPos2(nNode),ePos2(nNode)

   integer :: i, j, k

   integer(KIND=MPI_OFFSET_KIND) :: xPos1, xPos2,  yPos1, yPos2

   getSendNodesNum = 0

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
                        getSendNodesNum = getSendNodesNum + 1
                     end if
                end do
             else                 ! grid->seq
                yPos1 = sPos2(j)
                yPos2 = ePos2(j)
                 if ((xPos2<yPos1) .OR. (xPos1>yPos2)) then
                     cycle
                 else
                     getSendNodesNum = getSendNodesNum + 1
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
             
          if (bNum2(j)==0) then   ! grid->grid
             do k=1, sN2
                yPos1 = sPos2(j)+(k-1)*sLen2
                yPos2 = ePos2(j)+(k-1)*sLen2
                 if ((xPos2<yPos1) .OR. (xPos1>yPos2)) then
                    cycle
                 else
                    getSendNodesNum = getSendNodesNum + 1
                 end if
             end do

          else
             yPos1 = sPos2(j)
             yPos2 = ePos2(j)
              if ((xPos2<yPos1) .OR. (xPos1>yPos2)) then
                   cycle
              else
                  getSendNodesNum = getSendNodesNum + 1
             end if
          end if
       end do 
 
   end if

end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

