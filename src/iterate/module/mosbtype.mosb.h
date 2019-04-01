!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        Subroutines to calculate data structure for MOSB                c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine calViSeq(sF,sN, mySeq)
        integer, intent(IN) :: sF, sN(sF)
        TYPE(CSeqInfo),  intent(OUT):: mySeq

        myseq%mLevel = sF
        myseq%mSize(1:sF) = sN(1:sF)

        call calSeq(mySeq)
   
    end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine calH0Seq(sNDVR, sF, sN, mySeq)
        logical, intent(IN) :: sNDVR
        integer, intent(IN) :: sF, sN(sF)
        TYPE(CSeqInfo),  intent(OUT):: mySeq

        myseq%mLevel = sF
        myseq%mSize(1:sF) = sN(1:sF)**2

        if (sNDVR) myseq%mSize(sF)=0
    
	call calSeq(mySeq)

    end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine calHOSBGrid(sST, myGD, myData, myGrid)
        logical, intent(IN) :: sST 
        TYPE(GDataInfo), intent(IN) :: myGD
        TYPE(MDataInfo), intent(IN) :: myData
        TYPE(MGridInfo), intent(OUT):: myGrid

        integer ::  sF

        sF = myGD%sF;  

        myGrid%gSize(1:sF)=myGD%gN*myGD%sN(1:sF)
        
        myGrid%pSize(1:sF)=myData%pBlk(1:sF)*myData%pDim(1:sF)*myGD%sN(1:sF)**2

        call calMGrid(sF,myGrid)

        myGrid%gPos(1:sF) = myGrid%gStart(1:sF)+(myData%gBStart(1:sF)-1)*  &
 	         myGD%gDim(1:sF)*(myGD%sN(1:sF)**2)+myData%gDStart(1:sF)-1

        if (.NOT. sST) then
           myGrid%pStart(1:sF)=1
           myGrid%pEnd(1:sF)=myGrid%pSize(1:sF)
           myGrid%pLen = MaxVal(myGrid%pSize(1:sF))
        end if
    end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine calVOSBGrid(myGD, myData, myGrid)
        TYPE(GDataInfo), intent(IN) :: myGD
        TYPE(MDataInfo), intent(IN) :: myData
        TYPE(MGridInfo), intent(OUT):: myGrid

        integer :: i, sF

        sF = myGD%sF

        myGrid%gSize(1:sF) = myGD%gBlk(1:sF)*myGD%sN(1:sF)**2

        myGrid%pSize(1:sF) = myData%pBlk(1:sF)*myGD%sN(1:sF)**2

        call calMGrid(sF,myGrid)

        myGrid%gPos(1:sF) = myGrid%gStart(1:sF)+(myData%gBStart(1:sF)-1)*  &
 	           (myGD%sN(1:sF)**2)


    end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine calRESGrid(myGD, myData, myGrid)
        TYPE(GDataInfo), intent(IN) :: myGD
        TYPE(MDataInfo), intent(IN) :: myData
        TYPE(MGridInfo), intent(OUT):: myGrid

        integer :: sF

        sF = myGD%sF;  

        myGrid%gSize(1:sF) = myGD%gBlk(1:sF)*myGD%gDim(1:sF)*myGD%sN(1:sF)
        
        myGrid%pSize(1:sF) = myData%pBlk(1:sF)*myData%pDim(1:sF)*myGD%sN(1:sF)

        call calMGrid(sF,myGrid)

        myGrid%gLen = myGrid%gSize(1)

  	myGrid%gStart(1:sF)=1;	    myGrid%gEnd(1:sF)=myGrid%gSize(1:sF)
        myGrid%pStart(1:sF)=1;      myGrid%pEnd(1:sF)=myGrid%pSize(1:sF)
        myGrid%gPos(1:sF) = (myData%gBStart(1:sF)-1)*     &
 	             myGD%gDim(1:sF)*myGD%sN(1:sF)+myData%gDStart(1:sF)

        myGrid%pLen = MaxVal(myGrid%pSize(1:sF))

    end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine calDepGrid(sDep, myGD, myData,  myGrid)
        TYPE(GDataInfo), intent(IN) :: myGD
        logical, intent(IN)  :: sDep(myGD%sF)
        TYPE(MDataInfo), intent(IN) :: myData
        TYPE(MGridInfo), intent(OUT):: myGrid

        integer :: i, sF

        sF = myGD%sF
        do i = 1, sF
 	   if (sDep(i)) then
               myGrid%gSize(i) = myGD%gDim(i)
               myGrid%pSize(i) = myData%pDim(i)
           else
 	       myGrid%gSize(i) = 0
               myGrid%pSize(i) = 0
           end if
        end do        

        call calMGrid(sF,myGrid)

        do i = 1, sF          	 
	   if (sDep(i)) then
               myGrid%gPos(i) = myGrid%gStart(i)+myData%gDStart(i)-1
           else
	       myGrid%gPos(i)=myGrid%gStart(i)
           end if
        end do

    end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
