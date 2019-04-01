!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        Subroutines to print information of data structure              c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     subroutine printSConv(myconv)
         TYPE(CConvSimple), intent(IN) :: myconv
         
         print *
         print *, ' ==============================================='
         print *, '    Information of Simple Convergence Testing   '
         print *, ' ==============================================='
         print *, '  Max. Number of iteration:', myconv%mMax
         print *, '  Convergence  Tolerance :', myconv%mTol
         print *, ' ========  End of Simple Convergence ==========='

     end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     subroutine printCConv(myconv)
         TYPE(CConv), intent(IN) :: myconv
         
         print *
         print *, ' ==============================================='
         print *, '       Information of Convergence Testing       '
         print *, ' ==============================================='
         write (*,100) myconv%mStart, myconv%mStep, myconv%mMax
         write (*,200) myconv%mNum, myconv%mGap
         write (*,300) myconv%mE0, myconv%mTol
         print *, ' ========  End of Convergence ============' 

         100 FORMAT('  Start # of Iter.:', I6, '. Inc. Step #: ',I6,        &
                     '. Max. # of Iter.:',I6)
         200 FORMAT('  # of Interest. Values:',I6,'. Gap for the Conv.:',I6)
         300 FORMAT('  Center of Interest. Values:',F15.9,'. Conv.  Tol :',F15.9)

     end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     subroutine printOsbw(myosbw)
         TYPE(COsbw), intent(IN) :: myosbw
         
         print *
         print *, ' ==============================================='
         print *, '     Information of OSBW/OSBD Preconditioner    '
         print *, ' ==============================================='
         write (*,100) myosbw%mE0, myosbw%mDe
         write (*,200) myosbw%mE0-myosbw%mDe, myosbw%mE0+myosbw%mDe
         write (*,300) myosbw%mBeta, myosbw%mCnt
         print *, ' =========    End of OSBD/OSBW     =============' 

         100 FORMAT('  Center of Window:',F15.9,'. Window Size:',F15.9) 
         200 FORMAT('  Window Range: [', F15.9, ',', F15.9, ']')
         300 FORMAT('  Param. for OSBD: Beta:', F15.9,'.  OSBW Window Size:', I10)

     end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     subroutine printSeqInfo(myseq)
         TYPE(CSeqInfo), intent(IN) :: myseq

         integer :: i
         
         print *
         print *, ' ==============================================='
         print *, '     Information of Multi-layer 1D Data         '
         print *, ' ==============================================='
         print *, ' Number of Layers:', myseq%mLevel
         print *, ' Total length of the 1D data:', myseq%mLen
         print *, ' Max. Size of data at 1 layer:', myseq%mMaxSize
         print *, ' -----------------------------------------------'
         print *, '  Layer        Start         End         Size   '
         print *, '------------------------------------------------'
         do i = 1, myseq%mLevel
     	    write (*,100) i, myseq%mStart(i), myseq%mEnd(i), myseq%mSize(i)
         end do
         print *, '=========  End of 1D data index info   ==============' 

         100 FORMAT(I8,3(2x,I10)) 

     end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine readGData(myData)
       TYPE(GDataInfo), intent(OUT) :: myData

       read(*,*) myData%sF
       read(*,*) myData%sN(1:myData%sF)

    end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine printGData(dInfo)
       TYPE(GDataInfo), intent(IN) :: dInfo

       integer :: i

       print *, ' ============================================================== '
       print *, '            Global  Layer  Configuration  Information           '
       print *, ' ============================================================== '
       print *, '    Number of Data Layers: ', dInfo%sF
       print *, '    Layer Config.:', dInfo%sN(1:dInfo%sF)
       print *, '    Total length:', dInfo%gN
       print *, '    Max. Length at one Layer:', dInfo%gMaxLen
       print *
       print *, ' ==========         Layer  Information:    ===================='
       print *, ' -----------------------------------------------------------------'
       print *, '    Layer | Grid Size |    Dim  Size   |  Block  Size  | Block Len'
       print *, ' -----------------------------------------------------------------'
       do i = 1, dInfo%sF
          write(*,100) i,dInfo%SN(i),dInfo%gDim(i),dInfo%gBlk(i), dInfo%gLen(i)
       end do
       print *, ' -----------------------------------------------------------------'
       print *, ' ============    Finish Layer Info   ============================='

      100 FORMAT(2X,I6, 1X, I6, 3(1X,I16))

    end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine printMNode(sF, nInfo)
       integer, intent(IN) :: sF
       TYPE(MNodeInfo), intent(IN) :: nInfo

       integer :: i

       print *
       print *, ' ============================================================= '
       print *, '               Layer  Information  of  MPI  Node               '
       print *, ' ============================================================= '
       print *
       print *, '     Node ID of current Node:', nInfo%id
       print *, '     Number of Computing Nodes:', nInfo%nNodes

       print *, '     Number of Layers:', sF
       print *, "     Split level: ", nInfo%spLevel
       print *, '     Load Balance Flag:',nInfo%lbFlag
       print *
       print *, " -------------------------------------------------------------------"
       print *, "     Layer    myID   GrpID     commID       #Grp      Num    nStart "
       print *, " -------------------------------------------------------------------"
       do i = 1, sF
          write(*,100) i, nInfo%myID(i),nInfo%grpID(i), nInfo%commID(i),   &
                        nInfo%nGroup(i),nInfo%nodNum(i),nInfo%nodIDStart(i)
       end do

       print *, ' ===========    Finish Layer Info of Nodes    ==================='

      100 FORMAT(2X,3(I7,1X),I14, 1X, 3(I7,1X))

    end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine printMData(sF,dInfo)
       integer, intent(IN) :: sF
       TYPE(MDataInfo), intent(IN) :: dInfo

       integer :: i

       print *, ' ============================================== '
       print *, '          Layer Information of MPI Data      '
       print *, ' ==============================================='
       print *, "    Number of Data Layers: ", sF
       print *, '    Max. Length of local data at one layer: ' , dInfo%pMaxLen 

       print *
       print *, ' ==========     Block Information:    ================='
       print *, " ------------------------------------------------------"
       print *, "    Layer        Start       End    pSize    pBlk*pDim "
       print *, " ------------------------------------------------------"
       do i = 1, sF
          write(*,100) i,dInfo%gBStart(i),dInfo%gBEnd(i),dInfo%pBlk(i),dInfo%pLen(i)
       end do

       print *
       print *, ' ===========   Dimenion Information:  ================='
       print *, " ------------------------------------------------------"
       print *, "   Layer       Start       End    pSize    pBlk*pDim "
       print *, " ------------------------------------------------------"
       do i = 1, sF
          write(*,100) i,dInfo%gDStart(i),dInfo%gDEnd(i),dInfo%pDim(i),dInfo%pLen(i)
       end do

       print *, ' ============    Finish Layer Info   =================='

      100 FORMAT(2X,I6, 1X, 2I12, I6, 2I12)

    end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine printMGrid(sF,dInfo)
       integer, intent(IN) :: sF
       TYPE(MGridInfo), intent(IN) :: dInfo

       integer :: i

       print *, ' =========================================================== '
       print *, '           Layer   Information   of   MPI   Grid   Data      '
       print *, ' =========================================================== '
       print *, '    Number of Data Layers:' , sF
       print *, '    Total Length of local data: ' , dInfo%pLen
       print *, '    Max. Size of local data at one layer: ' , dInfo%pMaxSize 

       print *, '    Total Length of global data: ' , dInfo%gLen 
       print *, '    Max. Size of global data at one layer:', dInfo%gMaxSize

       print *
       print *, ' ===========    Local Data Distribution:    ================ '
       print *, " ------------------------------------------------------------"
       print *, "   Layer     pStart      pEnd      pSize     gPosition   "
       print *, " ------------------------------------------------------------"
       do i = 1, sF
          write(*,100) i,dInfo%pStart(i),dInfo%pEnd(i),dInfo%pSize(i),dInfo%gPos(i)
       end do

       print *
       print *, ' ===========   Global Data Distribution :  ================='
       print *, " ------------------------------------------------------"
       print *, "   Layer          gStart         gEnd           gSize "
       print *, " ------------------------------------------------------"
       do i = 1, sF
          write(*,200) i,dInfo%gStart(i),dInfo%gEnd(i),dInfo%gSize(i)
       end do

       print *, ' ============    Finish Layer Info   =================='

      100 FORMAT(2X,I5,1x,3I10, I15)
      200 FORMAT(2X,I5,1x, 3I15)

    end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
