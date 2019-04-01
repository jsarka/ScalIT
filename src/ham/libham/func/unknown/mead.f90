subroutine minim(p, step, nop, func, maxfn, iprint, stopcr, nloop, iquad,  &
                 simp, var, functn, ifault)
!     A PROGRAM FOR FUNCTION MINIMIZATION USING THE SIMPLEX METHOD.
!     FOR DETAILS, SEE NELDER & MEAD, THE COMPUTER JOURNAL, JANUARY 1965
!     PROGRAMMED BY D.E.SHAW,
!     CSIRO, DIVISION OF MATHEMATICS & STATISTICS
!     P.O. BOX 218, LINDFIELD, N.S.W. 2070
!     WITH AMENDMENTS BY R.W.M.WEDDERBURN
!     ROTHAMSTED EXPERIMENTAL STATION
!     HARPENDEN, HERTFORDSHIRE, ENGLAND
!     Further amended by Alan Miller
!     CSIRO Division of Mathematics & Statistics
!     Private Bag 10, CLAYTON, VIC. 3169
!     Fortran 90 conversion by Alan Miller, June 1995
!     Latest revision - 14 September 1995
!     ARGUMENTS:-
!     P()     = INPUT, STARTING VALUES OF PARAMETERS
!               OUTPUT, FINAL VALUES OF PARAMETERS
!     STEP()  = INPUT, INITIAL STEP SIZES
!     NOP     = INPUT, NO. OF PARAMETERS, INCL. ANY TO BE HELD FIXED
!     FUNC    = OUTPUT, THE FUNCTION VALUE CORRESPONDING TO THE FINAL
!                 PARAMETER VALUES.
!     maxfn     = INPUT, THE MAXIMUM NO. OF FUNCTION EVALUATIONS ALLOWED.
!               Say, 20 times the number of parameters, NOP.
!     IPRINT  = INPUT, PRINT CONTROL PARAMETER
!                 < 0 NO PRINTING
!                 = 0 PRINTING OF PARAMETER VALUES AND THE FUNCTION
!                     VALUE AFTER INITIAL EVIDENCE OF CONVERGENCE.
!                 > 0 AS FOR IPRINT = 0 PLUS PROGRESS REPORTS AFTER
!                     EVERY IPRINT EVALUATIONS, PLUS PRINTING FOR THE
!                     INITIAL SIMPLEX.
!     STOPCR  = INPUT, STOPPING CRITERION.
!               The criterion is applied to the standard deviation of
!               the values of FUNC at the points of the simplex.
!     NLOOP   = INPUT, THE STOPPING RULE IS APPLIED AFTER EVERY NLOOP
!               FUNCTION EVALUATIONS.   Normally NLOOP should be slightly
!               greater than NOP, say NLOOP = 2*NOP.
!     IQUAD   = INPUT, = 1 IF FITTING OF A QUADRATIC SURFACE IS REQUIRED
!                      = 0 IF NOT
!               N.B. The fitting of a quadratic surface is strongly
!               recommended, provided that the fitted function is
!               continuous in the vicinity of the minimum.   It is often
!               a good indicator of whether a premature termination of
!               the search has occurred.
!     SIMP    = INPUT, CRITERION FOR EXPANDING THE SIMPLEX TO OVERCOME
!               ROUNDING ERRORS BEFORE FITTING THE QUADRATIC SURFACE.
!               The simplex is expanded so that the function values at
!               the points of the simplex exceed those at the supposed
!               minimum by at least an amount SIMP.
!     VAR()   = OUTPUT, CONTAINS THE DIAGONAL ELEMENTS OF THE INVERSE OF
!               THE INFORMATION MATRIX.
!     FUNCTN  = INPUT, NAME OF THE USER'S SUBROUTINE - ARGUMENTS (P,FUNC)
!               WHICH RETURNS THE FUNCTION VALUE FOR A GIVEN SET OF
!               PARAMETER VALUES IN ARRAY P.
!****     FUNCTN MUST BE DECLARED EXTERNAL IN THE CALLING PROGRAM.
!     IFAULT  = OUTPUT, = 0 FOR SUCCESSFUL TERMINATION
!                 = 1 IF MAXIMUM NO. OF FUNCTION EVALUATIONS EXCEEDED
!                 = 2 IF INFORMATION MATRIX IS NOT +VE SEMI-DEFINITE
!                 = 3 IF NOP < 1
!                 = 4 IF NLOOP < 1
!     N.B. P, STEP AND VAR (IF IQUAD = 1) MUST HAVE DIMENSION AT LEAST NOP
!          IN THE CALLING PROGRAM.
!*****************************************************************************
implicit double precision (a-h,o-z)
integer, intent(IN)             :: nop, maxfn, iprint, nloop, iquad
integer, intent(OUT)            :: ifault
double precision, intent(IN)    :: stopcr, simp
double precision, intent(INOUT) :: p(nop), step(nop)
double precision, intent(OUT)   :: var(nop), func
external functn
!     Local variables
double precision   :: g(nop+1,nop), h(nop+1), pbar(nop), pstar(nop),      &
                      pstst(nop), aval(nop), pmin(nop), temp(nop),        &
                      bmat(nop*(nop+1)/2), vc(nop*(nop+1)/2),             &
                      zero = 0.d0, half = 0.5d0, one = 1.d0, two = 2.d0
!     A = REFLECTION COEFFICIENT, B = CONTRACTION COEFFICIENT, AND
!     C = EXPANSION COEFFICIENT.
double precision, parameter :: a = 1.d0, b = 0.5d0, c = 2.d0
!     SET LOUT = LOGICAL UNIT NO. FOR OUTPUT
integer, parameter :: lout = 6
!     IF PROGRESS REPORTS HAVE BEEN REQUESTED, PRINT HEADING
if (iprint.gt.0) write (lout,5000) iprint
!     CHECK INPUT ARGUMENTS
ifault = 0
if (nop.le.0) ifault = 3
if (nloop.le.0) ifault = 4
if (ifault.ne.0) return
!     SET NAP = NO. OF PARAMETERS TO BE VARIED, I.E. WITH STEP.NE.0
nap = count(step.ne.zero)
neval = 0
loop = 0
iflag = 0
!     IF NAP = 0 EVALUATE FUNCTION AT THE STARTING POINT AND RETURN
if (nap.le.0) then
  call functn(p,func)
  return
end if
!     SET UP THE INITIAL SIMPLEX
20 g(1,:) = p
irow = 2
do i = 1, nop
  if (step(i).ne.zero) then
    g(irow,:) = p
    g(irow,i) = p(i) + step(i)
    irow = irow + 1
  end if
end do
np1 = nap + 1
do i = 1, np1
  p = g(i,:)
  call functn(p,h(i))
  neval = neval + 1
  if (iprint.gt.0) then
    write (lout,5100) neval, h(i), p
  end if
end do
!     START OF MAIN CYCLE.
!     FIND MAX. & MIN. VALUES FOR CURRENT SIMPLEX (HMAX & HMIN).
Main_loop: do
  loop = loop + 1
  imax = 1
  imin = 1
  hmax = h(1)
  hmin = h(1)
  do i = 2, np1
    if (h(i).gt.hmax) then
      imax = i
      hmax = h(i)
    else
      if (h(i).lt.hmin) then
        imin = i
        hmin = h(i)
      end if
    end if
  end do
!     FIND THE CENTROID OF THE VERTICES OTHER THAN P(IMAX)
  pbar = zero
  do i = 1, np1
    if (i.ne.imax) then
      pbar = pbar + g(i,:)
    end if
  end do
  pbar = pbar / FLOAT(nap)
!     REFLECT MAXIMUM THROUGH PBAR TO PSTAR,
!     HSTAR = FUNCTION VALUE AT PSTAR.
  pstar = a * (pbar - g(imax,:)) + pbar
  call functn(pstar,hstar)
  neval = neval + 1
  if (iprint.gt.0) then
    if (mod(neval,iprint).eq.0) write (lout,5100) neval, hstar, pstar
  end if
!     IF HSTAR < HMIN, REFLECT PBAR THROUGH PSTAR,
!     HSTST = FUNCTION VALUE AT PSTST.
  if (hstar.lt.hmin) then
    pstst = c * (pstar - pbar) + pbar
    call functn(pstst,hstst)
    neval = neval + 1
    if (iprint.gt.0) then
      if (mod(neval,iprint).eq.0) write (lout,5100) neval, hstst, pstst
    end if
!     IF HSTST < HMIN REPLACE CURRENT MAXIMUM POINT BY PSTST AND
!     HMAX BY HSTST, THEN TEST FOR CONVERGENCE.
    if (hstst.ge.hmin) then   ! REPLACE MAXIMUM POINT BY PSTAR & H(IMAX) BY HSTAR.
      g(imax,:) = pstar
      h(imax) = hstar
    else
      g(imax,:) = pstst
      h(imax) = hstst
    end if
    GO TO 250
  end if
!     HSTAR IS NOT < HMIN.
!     TEST WHETHER IT IS < FUNCTION VALUE AT SOME POINT OTHER THAN
!     P(IMAX).   IF IT IS REPLACE P(IMAX) BY PSTAR & HMAX BY HSTAR.
  do i = 1, np1
    if (i.ne.imax) then
      if (hstar.lt.h(i)) then  ! REPLACE MAXIMUM POINT BY PSTAR & H(IMAX) BY HSTAR.
        g(imax,:) = pstar
        h(imax) = hstar
        GO TO 250
      end if
    end if
  end do
!     HSTAR > ALL FUNCTION VALUES EXCEPT POSSIBLY HMAX.
!     IF HSTAR <= HMAX, REPLACE P(IMAX) BY PSTAR & HMAX BY HSTAR.
  if (hstar.le.hmax) then
    g(imax,:) = pstar
    hmax = hstar
    h(imax) = hstar
  end if
!     CONTRACTED STEP TO THE POINT PSTST,
!     HSTST = FUNCTION VALUE AT PSTST.
  pstst = b * g(imax,:) + (one-b) * pbar
  call functn(pstst,hstst)
  neval = neval + 1
  if (iprint.gt.0) then
    if (mod(neval,iprint).eq.0) write (lout,5100) neval, hstst, pstst
  end if
!     IF HSTST < HMAX REPLACE P(IMAX) BY PSTST & HMAX BY HSTST.
  if (hstst.le.hmax) then
    g(imax,:) = pstst
    h(imax) = hstst
    GO TO 250
  end if
!     HSTST > HMAX.
!     SHRINK THE SIMPLEX BY REPLACING EACH POINT, OTHER THAN THE CURRENT
!     MINIMUM, BY A POINT MID-WAY BETWEEN ITS CURRENT POSITION AND THE
!     MINIMUM.
  do i = 1, np1
    if (i.ne.imin) then
      do j = 1, nop
        if (step(j).ne.zero) g(i,j) = (g(i,j) + g(imin,j)) * half
        p(j) = g(i,j)
      end do
      call functn(p,h(i))
      neval = neval + 1
      if (iprint.gt.0) then
        if (mod(neval,iprint).eq.0) write (lout,5100) neval, h(i), p
      end if
    end if
  end do
!     IF LOOP = NLOOP TEST FOR CONVERGENCE, OTHERWISE REPEAT MAIN CYCLE.
  250 if (loop.lt.nloop) cycle Main_loop
!     CALCULATE MEAN & STANDARD DEVIATION OF FUNCTION VALUES FOR THE
!     CURRENT SIMPLEX.
  hmean = sum( h(1:np1) ) / FLOAT(np1)
  hstd = sum( (h(1:np1) - hmean) ** 2 )
  hstd = sqrt(hstd / FLOAT(np1))
!     IF THE RMS > STOPCR, SET IFLAG & LOOP TO ZERO AND GO TO THE
!     START OF THE MAIN CYCLE AGAIN.
  if (hstd.gt.stopcr .and. neval.le.maxfn) then
    iflag = 0
    loop = 0
    cycle Main_loop
  end if
!     FIND THE CENTROID OF THE CURRENT SIMPLEX AND THE FUNCTION VALUE THERE.
  do i = 1, nop
    if (step(i).ne.zero) then
      p(i) = sum( g(1:np1,i) ) / FLOAT(np1)
    end if
  end do
  call functn(p,func)
  neval = neval + 1
  if (iprint.gt.0) then
    if (mod(neval,iprint).eq.0) write (lout,5100) neval, func, p
  end if
!     TEST WHETHER THE NO. OF FUNCTION VALUES ALLOWED, maxfn, HAS BEEN
!     OVERRUN; IF SO, EXIT WITH IFAULT = 1.
  if (neval.gt.maxfn) then
    ifault = 1
    if (iprint.lt.0) return
    write (lout,5200) maxfn
    write (lout,5300) hstd
    write (lout,5400) p
    write (lout,5500) func
    return
  end if
!     CONVERGENCE CRITERION SATISFIED.
!     IF IFLAG = 0, SET IFLAG & SAVE HMEAN.
!     IF IFLAG = 1 & CHANGE IN HMEAN <= STOPCR THEN SEARCH IS COMPLETE.
  if (iprint.ge.0) then
    write (lout,5600)
    write (lout,5400) p
    write (lout,5500) func
  end if
  if (iflag.eq.0 .or. abs(savemn-hmean).ge.stopcr) then
    iflag = 1
    savemn = hmean
    loop = 0
  else
    exit Main_loop
  end if
end do Main_loop
if (iprint.ge.0) then
  write (lout,5700) neval
  write (lout,5800) p
  write (lout,5900) func
end if
if (iquad.le.0) return
!------------------------------------------------------------------
!     QUADRATIC SURFACE FITTING
if (iprint.ge.0) write (lout,6000)
!     EXPAND THE FINAL SIMPLEX, IF NECESSARY, TO OVERCOME ROUNDING
!     ERRORS.
hmin = func
nmore = 0
do i = 1, np1
  do
    test = abs(h(i)-func)
    if (test.lt.simp) then
      do j = 1, nop
        if (step(j).ne.zero) g(i,j) = (g(i,j)-p(j)) + g(i,j)
        pstst(j) = g(i,j)
      end do
      call functn(pstst,h(i))
      nmore = nmore + 1
      neval = neval + 1
      if (h(i).ge.hmin) cycle
      hmin = h(i)
      if (iprint.ge.0) write (lout,5100) neval, hmin, pstst
    else
      exit
    end if
  end do
end do
!     FUNCTION VALUES ARE CALCULATED AT AN ADDITIONAL NAP POINTS.
do i = 1, nap
  i1 = i + 1
  pstar = (g(1,:) + g(i1,:)) * half
  call functn(pstar,aval(i))
  nmore = nmore + 1
  neval = neval + 1
end do
!     THE MATRIX OF ESTIMATED SECOND DERIVATIVES IS CALCULATED AND ITS
!     LOWER TRIANGLE STORED IN BMAT.
a0 = h(1)
do i = 1, nap
  i1 = i - 1
  i2 = i + 1
  do j = 1, i1
    j1 = j + 1
    pstst = (g(i2,:) + g(j1,:)) * half
    call functn(pstst,hstst)
    nmore = nmore + 1
    neval = neval + 1
    l = i * (i-1) / 2 + j
    bmat(l) = two * (hstst + a0 - aval(i) - aval(j))
  end do
end do
l = 0
do i = 1, nap
  i1 = i + 1
  l = l + i
  bmat(l) = two * (h(i1) + a0 - two*aval(i))
end do
!     THE VECTOR OF ESTIMATED FIRST DERIVATIVES IS CALCULATED AND
!     STORED IN AVAL.
do i = 1, nap
  i1 = i + 1
  aval(i) = two * aval(i) - (h(i1) + 3.d0*a0) * half
end do
!     THE MATRIX Q OF NELDER & MEAD IS CALCULATED AND STORED IN G.
pmin = g(1,:)
do i = 1, nap
  i1 = i + 1
  g(i1,:) = g(i1,:) - g(1,:)
end do
do i = 1, nap
  i1 = i + 1
  g(i,:) = g(i1,:)
end do
!     INVERT BMAT
call syminv(bmat, nap, bmat, temp, nullty, ifault, rmax)
if (ifault.eq.0) then
  irank = nap - nullty
else                                 ! BMAT not +ve definite
                                     ! Resume search for the minimum
  if (iprint.ge.0) write (lout,6100)
  ifault = 2
  if (neval.gt.maxfn) return
  write (lout,6200)
  step = half * step
  GO TO 20
end if
!     BMAT*A/2 IS CALCULATED AND STORED IN H.
do i = 1, nap
  h(i) = zero
  do j = 1, nap
    if (j.le.i) then
      l = i * (i-1) / 2 + j
    else
      l = j * (j-1) / 2 + i
    end if
    h(i) = h(i) + bmat(l) * aval(j)
  end do
end do
!     FIND THE POSITION, PMIN, & VALUE, YMIN, OF THE MINIMUM OF THE
!     QUADRATIC.
ymin = dot_product( h(1:nap), aval(1:nap) )
ymin = a0 - ymin
do i = 1, nop
  pstst(i) = dot_product( h(1:nap), g(1:nap,i) )
end do
pmin = pmin - pstst
if (iprint.ge.0) then
  write (lout,6300) ymin, pmin
  write (lout,6400)
end if
!     Q*BMAT*Q'/2 IS CALCULATED & ITS LOWER TRIANGLE STORED IN VC
do i = 1, nop
  do j = 1, nap
    h(j) = zero
    do k = 1, nap
      if (k.le.j) then
        l = j * (j-1) / 2 + k
      else
        l = k * (k-1) / 2 + j
      end if
      h(j) = h(j) + bmat(l) * g(k,i) * half
    end do
  end do
  do j = i, nop
    l = j * (j-1) / 2 + i
    vc(l) = dot_product( h(1:nap), g(1:nap,j) )
  end do
end do
!     THE DIAGONAL ELEMENTS OF VC ARE COPIED INTO VAR.
j = 0
do i = 1, nop
  j = j + i
  var(i) = vc(j)
end do
if (iprint.lt.0) return
write (lout,6500) irank
call print_tri_matrix(vc, nop, lout)
write (lout,6600)
call syminv(vc, nap, bmat, temp, nullty, ifault, rmax)
!     BMAT NOW CONTAINS THE INFORMATION MATRIX
write (lout,6700)
call print_tri_matrix(bmat, nop, lout)
ii = 0
ij = 0
do i = 1, nop
  ii = ii + i
  if (vc(ii).gt.zero) then
    vc(ii) = one / sqrt(vc(ii))
  else
    vc(ii) = zero
  end if
  jj = 0
  do j = 1, i - 1
    jj = jj + j
    ij = ij + 1
    vc(ij) = vc(ij) * vc(ii) * vc(jj)
  end do
  ij = ij + 1
end do
write (lout,6800)
ii = 0
do i = 1, nop
  ii = ii + i
  if (vc(ii).ne.zero) vc(ii) = one
end do
call print_tri_matrix(vc, nop, lout)
!     Exit, on successful termination.
650 write (lout,6900) nmore
return
5000 format (' Progress Report every',i4,' function evaluations'/, &
             ' EVAL.   FUNC.VALUE.',10X,'PARAMETER VALUES')
5100 format (/1X,i4,2X,g12.5,2X,5G11.4,3(/21X,5G11.4))
5200 format (' No. of function evaluations > ',i5)
5300 format (' RMS of function values of last simplex =',g14.6)
5400 format (' Centroid of last simplex =',4(/1X,6G13.5))
5500 format (' Function value at centroid =',g14.6)
5600 format (/' EVIDENCE OF CONVERGENCE')
5700 format (/' Minimum found after',i5,' function evaluations')
5800 format (' Minimum at',4(/1X,6G13.6))
5900 format (' Function value at minimum =',g14.6)
6000 format (/' Fitting quadratic surface about supposed minimum'/)
6100 format (/' MATRIX OF ESTIMATED SECOND DERIVATIVES NOT +VE DEFN.'/ &
             ' MINIMUM PROBABLY NOT FOUND'/)
6200 format (/10X,'Search restarting'/)
6300 format (' Minimum of quadratic surface =',g14.6,' at',4(/1X,6G13.5))
6400 format (' IF THIS DIFFERS BY MUCH FROM THE MINIMUM ESTIMATED',1X,       &
             'FROM THE MINIMIZATION,'/' THE MINIMUM MAY BE FALSE &/OR THE '  &
             'INFORMATION MATRIX MAY BE',1X,'INACCURATE'/)
6500 format (' Rank of information matrix =',i3/ &
             ' Inverse of information matrix:-')
6600 format (/' If the function minimized was -LOG(LIKELIHOOD),'/         &
             ' this is the covariance matrix of the parameters.'/         &
             ' If the function was a sum of squares of residuals,'/       &
             ' this matrix must be multiplied by twice the estimated ',   &
             'residual variance'/' to obtain the covariance matrix.'/)
6700 format (' INFORMATION MATRIX:-'/)
6800 format (/' CORRELATION MATRIX:-')
6900 format (/' A further',i4,' function evaluations have been used'/)
end subroutine minim
subroutine syminv(a, n, c, w, nullty, ifault, rmax)
!     ALGORITHM AS7, APPLIED STATISTICS, VOL.17, 1968.
!     ARGUMENTS:-
!     A()    = INPUT, THE SYMMETRIC MATRIX TO BE INVERTED, STORED IN
!                LOWER TRIANGULAR FORM
!     N      = INPUT, ORDER OF THE MATRIX
!     C()    = OUTPUT, THE INVERSE OF A (A GENERALIZED INVERSE IF C IS
!                SINGULAR), ALSO STORED IN LOWER TRIANGULAR.
!                C AND A MAY OCCUPY THE SAME LOCATIONS.
!     W()    = WORKSPACE, DIMENSION AT LEAST N.
!     NULLTY = OUTPUT, THE RANK DEFICIENCY OF A.
!     IFAULT = OUTPUT, ERROR INDICATOR
!                 = 1 IF N < 1
!                 = 2 IF A IS NOT +VE SEMI-DEFINITE
!                 = 0 OTHERWISE
!     RMAX   = OUTPUT, APPROXIMATE BOUND ON THE ACCURACY OF THE DIAGONAL
!                ELEMENTS OF C.  E.G. IF RMAX = 1.E-04 THEN THE DIAGONAL
!                ELEMENTS OF C WILL BE ACCURATE TO ABOUT 4 DEC. DIGITS.
!     LATEST REVISION - 1 April 1985
!***************************************************************************
implicit double precision (a-h,o-z)
dimension a(*), c(*), w(n)
double precision :: zero = 0.d0, one = 1.d0
nrow = n
ifault = 1
if (nrow.gt.0) then
  ifault = 0
!     CHOLESKY FACTORIZATION OF A, RESULT IN C
  call chola(a, nrow, c, nullty, ifault, rmax, w)
  if (ifault.eq.0) then
!     INVERT C & FORM THE PRODUCT (CINV)'*CINV, WHERE CINV IS THE INVERSE
!     OF C, ROW BY ROW STARTING WITH THE LAST ROW.
!     IROW = THE ROW NUMBER, NDIAG = LOCATION OF LAST ELEMENT IN THE ROW.
    nn = nrow * (nrow+1) / 2
    irow = nrow
    ndiag = nn
    10 if (c(ndiag).ne.zero) then
      l = ndiag
      do i = irow, nrow
        w(i) = c(l)
        l = l + i
      end do
      icol = nrow
      jcol = nn
      mdiag = nn
      30 l = jcol
      x = zero
      if (icol.eq.irow) x = one / w(irow)
      k = nrow
      40 if (k.ne.irow) then
        x = x - w(k) * c(l)
        k = k - 1
        l = l - 1
        if (l.gt.mdiag) l = l - k + 1
        GO TO 40
      end if
      c(l) = x / w(irow)
      if (icol.eq.irow) GO TO 60
      mdiag = mdiag - icol
      icol = icol - 1
      jcol = jcol - 1
      GO TO 30
    end if ! (c(ndiag).NE.zero)
    l = ndiag
    do j = irow, nrow
      c(l) = zero
      l = l + j
    end do
    60 ndiag = ndiag - irow
    irow = irow - 1
    if (irow.ne.0) GO TO 10
  end if
end if
return
end subroutine syminv
subroutine chola(a, n, u, nullty, ifault, rmax, r)
!     ALGORITHM AS6, APPLIED STATISTICS, VOL.17, 1968, WITH
!     MODIFICATIONS BY A.J.MILLER
!     ARGUMENTS:-
!     A()    = INPUT, A +VE DEFINITE MATRIX STORED IN LOWER-TRIANGULAR
!                FORM.
!     N      = INPUT, THE ORDER OF A
!     U()    = OUTPUT, A LOWER TRIANGULAR MATRIX SUCH THAT U*U' = A.
!                A & U MAY OCCUPY THE SAME LOCATIONS.
!     NULLTY = OUTPUT, THE RANK DEFICIENCY OF A.
!     IFAULT = OUTPUT, ERROR INDICATOR
!                 = 1 IF N < 1
!                 = 2 IF A IS NOT +VE SEMI-DEFINITE
!                 = 0 OTHERWISE
!     RMAX   = OUTPUT, AN ESTIMATE OF THE RELATIVE ACCURACY OF THE
!                DIAGONAL ELEMENTS OF U.
!     R()    = OUTPUT, ARRAY CONTAINING BOUNDS ON THE RELATIVE ACCURACY
!                OF EACH DIAGONAL ELEMENT OF U.
!     LATEST REVISION - 1 April 1985
!***************************************************************************
implicit double precision (a-h,o-z)
dimension a(*), u(*), r(n)
!     ETA SHOULD BE SET EQUAL TO THE SMALLEST +VE VALUE SUCH THAT
!     1.D0 + ETA IS CALCULATED AS BEING GREATER THAN 1.D0 IN THE ACCURACY
!     BEING USED.
double precision :: eta = 1.d-16, zero = 0.d0
ifault = 1
if (n.gt.0) then
  ifault = 2
  nullty = 0
  rmax = eta
  r(1) = eta
  j = 1
  k = 0
!     FACTORIZE COLUMN BY COLUMN, ICOL = COLUMN NO.
  do 50 icol = 1, n
    l = 0
!     IROW = ROW NUMBER WITHIN COLUMN ICOL
    do 30 irow = 1, icol
      k = k + 1
      w = a(k)
      if (irow.eq.icol) rsq = (w*eta) ** 2
      m = j
      do 10 i = 1, irow
        l = l + 1
        if (i.eq.irow) exit
        w = w - u(l) * u(m)
        if (irow.eq.icol) rsq = rsq + (u(l)**2*r(i)) ** 2
        m = m + 1
      10 continue
      if (irow.eq.icol) GO TO 40
      if (u(l).ne.zero) then
        u(k) = w / u(l)
      else
        u(k) = zero
        if (abs(w).gt.abs(rmax*a(k))) GO TO 60
      end if
    30 continue
!     END OF ROW, ESTIMATE RELATIVE ACCURACY OF DIAGONAL ELEMENT.
    40 rsq = sqrt(rsq)
    if (abs(w).gt.5.*rsq) then
      if (w.lt.zero) return
      u(k) = sqrt(w)
      r(i) = rsq / w
      if (r(i).gt.rmax) rmax = r(i)
    else
      u(k) = zero
      nullty = nullty + 1
    end if
    j = j + icol
  50 continue
  ifault = zero
end if
60 return
end subroutine chola

subroutine print_tri_matrix(a, n, lout)
implicit none
integer, intent(IN)          :: n, lout
double precision, intent(IN) :: a(*)
!     Local variables
integer  :: i, ii, i1, i2, l
l = 1
do l = 1, n, 6
  ii = l * (l-1) / 2
  do i = l, n
    i1 = ii + l
    ii = ii + i
    i2 = min(ii,i1+5)
    write (lout,'(1X,6G13.5)') a(i1:i2)
  end do
end do
return
end subroutine print_tri_matrix

