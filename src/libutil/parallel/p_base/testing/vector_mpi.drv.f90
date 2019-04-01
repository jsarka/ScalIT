!
! Testing subroutine for vector_mpi functions
!
program vecmpi
   include "mpif.h"
   integer, parameter :: N1 = 50
   integer, parameter :: N2 = 10
   double precision, dimension(N1, N2) :: locData
   double precision, dimension(N2)     :: normData
   double precision, dimension(N2, N2) :: dotData
   double precision, dimension(:, :), allocatable   :: glbData
   integer :: myID, node, ierr
   integer :: i, j
   double precision :: tmp1, tmp2, tmp3, tmp4
   double precision ::  norm_sv_mpi, norm2_sv_mpi, dotprod_sv_mpi

   CALL MPI_INIT(ierr)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD, myID, ierr)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD, node, ierr)

   if (myid==0) THEN
        allocate(glbData(node*N1, N2))
   else
        allocate(glbData(1, N2))
   end if

   call random_number(locData)

   do i = 1, N2
       CALL MPI_GATHER(locDATA(:, I), N1, MPI_DOUBLE_PRECISION,   &
                 glbData(:, I), N1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR )
   end do

   DO I = 1, N2
       if (myid == 0) then
            tmp1 = norm2_sv_mpi(N1, locData(:, I), MPI_COMM_WORLD, IERR)
            tmp2 = norm_sv_mpi(N1, locData(:, I), MPI_COMM_WORLD, IERR)
            tmp3 = dot_product(glbData(:,I), glbData(:, I))
            tmp4 = sqrt(tmp3)
            print *
            print *, 'Testing Norm Subroutines in MPI at ', I
            print *, tmp1, tmp2, tmp3, tmp4                                    
       else
            tmp1 = norm2_sv_mpi(N1, locData(:, I), MPI_COMM_WORLD, IERR)
            tmp2 = norm_sv_mpi(N1, locData(:, I), MPI_COMM_WORLD, IERR)
       end if
   END DO
   
   if (myid == 0) then
            print *
            print *, 'Multi Column NORM2:'
            call norm2_mv_mpi(N1, N2, locData, MPI_COMM_WORLD, normData, IERR)             
            DO I = 1, N2 
                 tmp1 = dot_product(glbDATA(:, I), glbData(:, I))
                 print *, I, normDATA(I), tmp1, tmp1-normData(I)
            END DO  
            print *
            print *, 'Multi Column Norm'
            call norm_mv_mpi(N1, N2, locData, MPI_COMM_WORLD, normData, IERR)             
            DO I = 1, N2 
                 tmp1 = DSQRT(dot_product(glbDATA(:, I), glbData(:, I)))
                 print *, I, normDATA(I), tmp1, tmp1-normData(I)
            END DO                                  
   else
            call norm2_mv_mpi(N1,N2, locData, MPI_COMM_WORLD, normData, IERR)
            call norm_mv_mpi(N1, N2,locData, MPI_COMM_WORLD, normData, IERR)
   end if

    DO I = 1, N2
       DO J = 1, N2
         if (myid == 0) then
            tmp1 = dotprod_sv_mpi(N1, locData(:, I),locData(:, J), MPI_COMM_WORLD, IERR)            
            tmp3 = dot_product(glbData(:,I), glbData(:, J))           
            print *
            print *, 'Testing Dotprod Subroutines in MPI at ', I, J
            print *, tmp1,  tmp3, tmp3-tmp1               
         else
            tmp1 = dotprod_sv_mpi(N1, locData(:, I),locData(:, J), MPI_COMM_WORLD, IERR)           
         end if
       END DO
   END DO

   CALL dotprod_mv_mpi(N1, N2, locData, N2, locData, MPI_COMM_WORLD, dotDATA, IERR)

   if (myid == 0) then
       print *
       print *, 'Multi Column Dot_Product'
       DO I = 1, N2
            DO J = 1, N2
                 tmp1 = dot_product(glbDATA(:, I), glbData(:, J))
                 print *, I, J, dotDATA(I, J), tmp1, tmp1-dotData(I, J)
            END DO
       END DO
      print *
   END IF

   deallocate(glbData)

   CALL MPI_FINALIZE(ierr)

end

