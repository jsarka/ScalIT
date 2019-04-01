	Module rnvar

	Implicit none
	save
	public

	DOUBLE PRECISION, ALLOCATABLE :: PL(:)
	Integer :: numneu
	Integer, parameter :: hidlay=2, numinp=3
	Integer :: numw,numb
	Double Precision, dimension (:), allocatable :: Saida
	Double Precision, dimension (:,:), allocatable :: Entrada
        Double Precision, parameter :: hardlim = 9.5629979367616899d-2 

	End Module rnvar
