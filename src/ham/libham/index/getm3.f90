! Get Max, Min and Median
!
! Input:  j0[3], random integer
! Output: j1[3], [max, median, min]
!         ind[3], order of j1[3] in original j0
! Return:
!  .TRUE. (1),  even permutation
!  .FALSE.(-1), odd permutation

logical function getMMM(j0, j1, ind)
    implicit none
    integer,intent(IN)  :: j0(3)
    integer,intent(OUT) :: j1(3),ind(3)
     
    integer :: jmax, jmin, jmid

    getMMM = .TRUE. ;
    if (j0(1) > j0(3)) then
        jmax   = j0(1); jmin   = j0(3)
        ind(1) = 1    ; ind(3) = 3
    else
        jmax = j0(3);   jmin=j0(1);
        ind(1) = 3  ;   ind(3) = 1;
         getMMM = (.NOT.  getMMM);
    endif
      
    jmid = j0(2); ind(2) = 2;      
    if (jmid > jmax) then
        ind(2)  = ind(1); ind(1) = 2 ;
         getMMM = (.NOT. getMMM) ;
    else
        if (jmid < jmin) then
              ind(2)=ind(3);ind(3)=2;
              getMMM = (.NOT. getMMM);
        endif
     endif
     j1(1)=j0(ind(1)); 
     j1(2)=j0(ind(2));
     j1(3)=j0(ind(3));

end
