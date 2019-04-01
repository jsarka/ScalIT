!
! Testing subroutine for vector_mpi functions
!
program vecmpi_cx
   implicit none
   include "mpif.h"
   integer, parameter :: N1 = 50
   integer, parameter :: N2 = 10
   double complex, dimension(N1, N2) :: locData
   double precision, dimension(N2)     :: normData
   double complex, dimension(N2, N2) :: dotData
   double complex, dimension(:, :), allocatable   :: glbData
   integer :: myID, node, ierr
   integer :: i, j
   double precision :: tmp1, tmp2, tmp3, tmp4
   double complex   :: cx1, cx2
   double precision ::  norm_svcx_mpi, norm2_svcx_mpi
   double complex   :: dotprod_svcx_mpi

   CALL MPI_INIT(ierr)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD, myID, ierr)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD, node, ierr)

   if (myid==0) THEN
        allocate(glbData(node*N1, N2))
   else
        allocate(glbData(1, N2))
   end if

   call randMat_cx(N1, N2, locData)

   do i = 1, N2
       CALL MPI_GATHER(locDATA(:, I), N1, MPI_DOUBLE_COMPLEX,   &
                 glbData(:, I), N1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, IERR )
   end do

   DO I = 1, N2
       if (myid == 0) then
            tmp1 = norm2_svcx_mpi(N1, locData(:, I), MPI_COMM_WORLD, IERR)
            tmp2 = norm_svcx_mpi(N1, locData(:, I), MPI_COMM_WORLD, IERR)
            tmp3 = DBLE(dot_product(glbData(:,I), glbData(:, I)))
            tmp4 = sqrt(tmp3)
            print *
            print *, 'Testing COMPLEX VERSION Norm Subroutines in MPI at ', I
            print *, tmp1, tmp2, tmp3, tmp4                                    
       else
            tmp1 = norm2_svcx_mpi(N1, locData(:, I), MPI_COMM_WORLD, IERR)
            tmp2 = norm_svcx_mpi(N1, locData(:, I), MPI_COMM_WORLD, IERR)
       end if
   END DO
   
   if (myid == 0) then
            print *
            print *, 'Multi Column NORM2 for complex version:'
            call norm2_mvcx_mpi(N1, N2, locData, MPI_COMM_WORLD, normData, IERR)             
            DO I = 1, N2 
                 tmp1 = DBLE(dot_product(glbDATA(:, I), glbData(:, I)))
                 print *, I, normDATA(I), tmp1, tmp1-normData(I)
            END DO  
            print *
            print *, 'Multi Column Norm for complex version'
            call norm_mvcx_mpi(N1, N2, locData, MPI_COMM_WORLD, normData, IERR)             
            DO I = 1, N2 
                 tmp1 = DSQRT(DBLE(dot_product(glbDATA(:, I), glbData(:, I))))
                 print *, I, normDATA(I), tmp1, tmp1-normData(I)
            END DO                                  
   else
            call norm2_mvcx_mpi(N1,N2, locData, MPI_COMM_WORLD, normData, IERR)
            call norm_mvcx_mpi(N1, N2,locData, MPI_COMM_WORLD, normData, IERR)
   end if

    DO I = 1, N2
       DO J = 1, N2
         if (myid == 0) then
            cx1 = dotprod_svcx_mpi(N1, locData(:, I),locData(:, J), MPI_COMM_WORLD, IERR) 
            cx2 = dot_product(glbData(:,I), glbData(:, J))           
            print *
            print *, 'Testing Dotprod Subroutines in MPI at ', I, J
            print *, cx1,  cx2, cx2-cx1               
         else
            cx1 = dotprod_svcx_mpi(N1, locData(:, I),locData(:, J), MPI_COMM_WORLD, IERR)           
         end if
       END DO
   END DO

   CALL dotprod_mvcx_mpi(N1, N2, locData, N2, locData, MPI_COMM_WORLD, dotDATA, IERR)

   if (myid == 0) then
       print *
       print *, 'Multi Column Dot_Product for complex version '
       DO I = 1, N2
            DO J = 1, N2
                 cx1 = dot_product(glbDATA(:, I), glbData(:, J))
                 print *, I, J, dotDATA(I, J), cx1, cx1-dotData(I, J)
            END DO
       END DO
      print *
   END IF

   deallocate(glbData)

   CALL MPI_FINALIZE(ierr)

end

