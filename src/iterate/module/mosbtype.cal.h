!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        Subroutines to calculate information of data structure          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     subroutine calGData(sF, sN, myData)
         integer, intent(IN) :: sF, sN(sF)
         TYPE(GDataInfo), intent(INOUT) :: myData
     
         integer :: i

         myData%sF = sF; 
         myData%sN(1:sF) = sN(1:sF)
       
         myData%gDim(1)=1
         do i = 1, sF-1
	    myData%gDim(i+1) = myData%gDim(i)*myData%sN(i)
         end do

         myData%gBlk(sF)=1
         do i = sF,2,-1
	    myData%gBlk(i-1) = myData%gBlk(i)*myData%sN(i)
         end do

         myData%gLen(1:sF) = myData%gDim(1:sF)*sN(1:sF)

         call getMaxLong(sF,myData%gLen, myData%gMaxLen)

         myData%gN = myData%gBlk(1)*myData%sN(1)

     end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     subroutine getMaxLong(sF, gLen, gMaxLen)
        integer, intent(IN) :: sF
	integer(kind=MPI_OFFSET_KIND), intent(IN) :: gLen(sF)
        integer(kind=MPI_OFFSET_KIND), intent(OUT) :: gMaxLen
   
        integer :: i     
     
	gMaxLen = gLen(1)

        do i = 2, sF
           if (gMaxLen<gLen(i)) gMaxLen = gLen(i)
        end do

     end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     subroutine calSeq(myseq)
         TYPE(CSeqInfo), intent(OUT) :: myseq

         call seqDataPos(myseq%mLevel, myseq%mSize, myseq%mLen,           &
                         myseq%mStart, myseq%mEnd)

         myseq%mMaxSize = MaxVal(myseq%mSize(1:myseq%mLevel))

     end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     subroutine calGDMNode (id, nNodes, myData, myNode)
         integer, intent(IN) :: id, nNodes
         Type(GDataInfo), intent(IN)  :: myData
         Type(MNodeInfo), intent(OUT) :: myNode

	call calMNode(id,nNodes,myData%sF, myData%gBlk, myNode)
     end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     subroutine calMNode (id, nNodes, sF, gBlk, myNode)
         integer, intent(IN)  :: id, nNodes, sF
 	 integer(kind=MPI_OFFSET_KIND) :: gblk(sF)
         Type(MNodeInfo), intent(OUT) :: myNode

         integer :: i, num, pnum, qnum

  	 myNode%id = id;  myNode%nNodes = nNodes
         myNode%nGroup(sF) = 1;              myNode%myID(sF) = myNode%id;
         myNode%nodNum(sF) = myNode%nNodes;  myNode%grpID(sF)= 0;       
         myNode%nodIDStart(sF) = 0;          myNode%spLevel = 0
         myNode%commID(1:sF) = MPI_COMM_WORLD
         
         do i = sF-1,1,-1

            if (gBlk(i) < myNode%nNodes) then
 	        myNode%nGroup(i) = gBlk(i)
                num  = myNode%nNodes / gBlk(i)
                pnum = myNode%nNodes - num * gBlk(i)
                qnum = (num+1)*pNum
    
                if (myNode%id < qNum) then
                   myNode%nodNum(i)  = num + 1
                   myNode%grpID(i) = myNode%id / (num+1)
                   myNode%myID(i)  = myNode%id - myNode%grpID(i)*(num+1)
                   myNode%nodIDStart(i)= myNode%grpID(i)*(num+1)
                else
                   myNode%nodNum(i)  = num
                   myNode%grpID(i) = (myNode%id-qnum)/num 
                   myNode%myID(i)  = (myNode%id-qnum)-myNode%grpID(i)*num
                   myNode%nodIDStart(i)= pNum*(num+1)+myNode%grpID(i)*num
                   myNode%grpID(i) = myNode%grpID(i)+pnum
                end if
            else
 	        if (myNode%spLevel==0) myNode%splevel=i+1

                myNode%nGroup(i) = myNode%nNodes
                myNode%myID(i)   = 0
                myNode%nodNum(i)   = 1
                myNode%grpID(i)  = myNode%id
                myNode%nodIDStart(i) = myNode%id
            end if
         end do

         i = myNode%spLevel-1
         pnum = gBlk(i) - (gBlk(i)/myNode%nNodes)*myNode%nNodes
         myNode%lbFlag = (pnum==0)

    end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine calMData(myGd,myNode,myData)
        TYPE(GDataInfo), intent(IN)  :: myGd
        TYPE(MNodeInfo), intent(IN)  :: myNode
        TYPE(MDataInfo), intent(OUT) :: myData

        integer :: i, sF, sP

        sF = myGd%sF
        sp = myNode%spLevel

        do i = 1, sP-1
           myData%pDim(i) = myGD%gDim(i)
           myData%gDStart(i) = 1
           call seqDataLenLong(myGD%gBlk(i), myNode%nNodes, myNode%ID,  &
                           myData%pBlk(i), myData%gBStart(i))

           myData%gBEnd(i) = myData%gBStart(i) + myData%pBlk(i) - 1
           myData%gDEnd(i) = myData%pDim(i)
        end do

        do i = sp, sF
           if (myNode%nodNum(i)==1) then
              myData%pBlk(i) = 1
	      myData%pDim(i) = myGD%gDim(i)
              myData%gDStart(i) = 1
           else
              myData%pBlk(i) = 1
              call seqDataLenLong(myGD%gDim(i),myNode%nodNum(i),myNode%myID(i), &
                              myData%pDim(i),myData%gDStart(i))
           end if
           myData%gBStart(i) = myNode%grpID(i) + 1
           myData%gDEnd(i) = myData%gDStart(i) + myData%pDim(i) - 1
           myData%gBEnd(i) = myData%gBStart(i)  
        end do

        myData%pLen(1:sF) = myData%pBlk(1:sF)*myData%pDim(1:sF)*myGD%sN(1:sF)

        myData%pMaxLen = MaxVal(myData%pLen(1:sF))

    end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine calMGrid(sF, myGrid)
        integer, intent(IN) :: sF
        TYPE(MGridInfo), intent(INOUT) :: myGrid

        integer :: i   

        call seqDataPosLong(sF,myGrid%gSize,myGrid%gLen,myGrid%gStart,myGrid%gEnd)
        call seqDataPos(sF,myGrid%pSize,myGrid%pLen,myGrid%pStart,myGrid%pEnd)       

        myGrid%pMaxSize = MaxVal(myGrid%pSize) 

        call getMaxLong(sF,myGrid%gSize,myGrid%gMaxSize)        

    end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

