        Function rnfunc(numvar,input)

	Use rnvar

	Implicit None

	integer :: numvar,i,j,k
	double precision :: input(numinp), rnfunc,rninj(numneu),rnink(numneu)

	rnfunc=0.d0
	rninj=0.d0
	rnink=0.d0


	do k=1,numneu
		do j=1,numneu
			rninj(j)=dtanh(sum(input(1:numinp)*PL((j-1)*numinp+1:(j-1)*numinp+numinp))+PL(numw+j))
		enddo
		rnink(k)=dtanh(sum(rninj(:)*PL(numinp*numneu+(k-1)*numneu+1:numinp*numneu+(k-1)*numneu+numneu))+PL(numneu+numw+k))
	enddo

	rnfunc=rnfunc+sum(rnink(:)*PL(numneu**hidlay+numinp*numneu+1:numneu**hidlay+numinp*numneu+numneu))+PL(hidlay*numneu+numw+1)

!        if (rnfunc.gt.hardlim) then
!		rnfunc=hardlim
!        end if

	End function rnfunc

