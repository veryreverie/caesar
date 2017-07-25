module spline2_module
  use string_module
  use io_module
  use constants_module, only : dp
contains

subroutine spline2()
  implicit none
  
  ! independent data points
  real(dp),allocatable :: x1b(:),fb(:)
  real(dp),allocatable :: x2c(:),fc(:)
  ! coupled data points
  real(dp),allocatable :: x1a(:),x2a(:),f(:,:)
  ! second derivatives
  real(dp),allocatable :: f2(:,:),f2b(:),f2c(:)
  ! interpolated data
  real(dp) :: x1,x2,f_int
  real(dp) :: dx1,dx2
  integer :: mode1,mode2
  
  type(string), allocatable :: input_file(:)
  type(string), allocatable :: modes_file(:)
  type(string), allocatable :: mode1_file(:)
  type(string), allocatable :: mode2_file(:)
  type(string), allocatable :: energy_file(:)
  integer                   :: output_file
  
  ! Temporary variables.
  integer                   :: n,m,p
  integer                   :: i,j,ialloc
  type(string), allocatable :: line(:)
  
  ! read in number of integration points.
  input_file = read_lines('fit_input.dat')
  line = split(input_file(1))
  n = int(line(1))
  m = int(line(2))
  p = int(line(3))
  
  ! allocate arrays.
  allocate( x1a(n),  &
          & x2a(n),  &
          & f(n,n),  &
          & f2(n,n), &
          & x1b(p),  &
          & x2c(p),  &
          & fb(p),   &
          & fc(p),   &
          & f2b(p),  &
          & f2c(p),  &
          & stat=ialloc); call err(ialloc)

  
  ! read in input data.
  modes_file = read_lines('modes.dat')
  line = split(modes_file(1))
  mode1 = int(line(1))
  mode2 = int(line(2))
  
  mode1_file = read_lines('fit_energy_'//mode1//'.dat')
  do i=1,p
    line = split(mode1_file(i))
    x1b(i) = dble(line(1))
    fb(i)  = dble(line(2))
  enddo ! n
  
  mode2_file = read_lines('fit_energy_'//mode2//'.dat')
  do i=1,p
    line = split(mode2_file(i))
    x2c(i) = dble(line(1))
    fc(i)  = dble(line(2))
  enddo
  
  energy_file = read_lines('fit_energy.dat')
  do i=1,n
    do j=1,n
      if(i/=((n+1)/2).and.j/=((n+1)/2))then 
        line = split(energy_file((i-1)*n+j))
        x1a(i) = dble(line(1))
        x2a(j) = dble(line(2))
        f(i,j) = dble(line(3))
      endif
    enddo ! j
  enddo ! i
  
  x1a((n+1)/2) = 0.0_dp
  x2a((n+1)/2) = 0.0_dp
  f((n+1)/2,:) = 0.0_dp
  f(:,(n+1)/2) = 0.0_dp
  
  ! Run calculation.
  f2b = spline(x1b,fb,1.0e35_dp,1.0e35_dp)
  f2c = spline(x2c,fc,1.0e35_dp,1.0e35_dp)
  do i=1,n
    do j=1,n
      if(i/=((n+1)/2).and.j/=((n+1)/2))then
        f(i,j)=f(i,j)-splint(x1b,fb,f2b,x1a(i))-splint(x2c,fc,f2c,x2a(j))
      endif
    enddo ! j
  enddo ! i
  
  f2 = splie2(x1a,x2a,f)
  
  dx1=2*abs(x1a(1))/m
  dx2=2*abs(x2a(1))/m
  do i=1,m+1
    x1 = x1a(1) + (i-1)*dx1
    do j=1,m+1
      x2 = x2a(1) + (i-1)*dx2
      f_int = splin2(x1a,x2a,f,f2,x1,x2)
      call print_line(output_file, x1//' '//x2//' '//f_int)
    enddo ! j
  enddo ! i
end subroutine

! ----------------------------------------------------------------------
!     SPLINE use: given an 1D array of X data and an array of Y data,
!     both of length N, this routine computes the 2nd derivatives, Y2 at
!     each X data point.  The user needs to specify the values of YP1
!     and YP2, which flags how the Y2 are computed at the edges.  For
!     natural spline fitting (recommended), set YP1 and YPN to numbers
!     greater than 1.0E+30.

!     this routine called once, prior to using routine SPLINT, as a set
!     up for using routine SPLINT, which performs the actual
!     interpolation

!     IMPORTANT NOTE: the X data values in array X must be in ascending
!     order or the interpolation will fail
! ----------------------------------------------------------------------
function spline(x,y,yp1,ypn) result(y2)
  implicit none
  
  real(dp), intent(in)  :: x(:)
  real(dp), intent(in)  :: y(:)
  real(dp), intent(in)  :: yp1
  real(dp), intent(in)  :: ypn
  real(dp), allocatable :: y2(:)
  
  integer :: n
  
  real(dp) :: sig
  real(dp) :: p
  real(dp), allocatable :: u(:)
  
  real(dp) :: qn
  real(dp) :: un
  
  integer :: i,k,ialloc
  
  n = size(x)
  if (size(y)/=n) then
    call err()
  endif
  
  allocate(y2(n), stat=ialloc); call err(ialloc)
  
  allocate(u(n-1), stat=ialloc); call err(ialloc)
  
  ! If yp1>1.0e+30 use natural spline,
  !    otherwise estimate y2 at the first point.
  if (yp1 > 0.99e30_dp) then
    y2(1) = 0
    u(1) = 0
  else
    y2(1) = -0.5_dp
    u(1) = (3/(x(2)-x(1))) * ((y(2)-y(1))/(x(2)-x(1))-yp1)
  endif
  
  ! Store intermediate values of terms in the expansion series.
  do i=2,n-1
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*y2(i-1)+2
    y2(i)=(sig-1)/p
    u(i)=(6*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))&
        &/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  enddo
  
  ! If YPN>1.0E+30 use natural spline,
  !    otherwise estimate Y2 at the last point point.
  if (ypn > 0.99e30_dp) then
    qn = 0
    un = 0
  else
    qn = 0.5_dp
    un = (3/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  endif
  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
  
  ! Compute the Y2 from the 2nd order expansion series.
  do k=n-1,1,-1
    y2(k)=y2(k)*y2(k+1)+u(k)
  enddo
end function

! ----------------------------------------------------------------------
!     SPLINT use: given an 1D array of XA data, an array of YA data, and
!     an array of the 2nd derivatives Y2A, all of length N, this routine
!     performs cubic spline interpolation, returning the interpolated
!     value Y at the user input value X.  The Y2A are computed in
!     routine SPLINE, which is called once before calling SPLINT.

!     IMPORTANT NOTE: the X data values in array X must be in ascending
!     order or the interpolation will fail
! ----------------------------------------------------------------------
function splint(xa,ya,y2a,x) result(y)
  implicit none
  
  real(dp), intent(in) :: xa(:)
  real(dp), intent(in) :: ya(:)
  real(dp), intent(in) :: y2a(:)
  real(dp), intent(in) :: x
  real(dp)             :: y
  
  integer :: below
  integer :: above
  integer :: k
  
  real(dp) :: h
  real(dp) :: a
  real(dp) :: b
  
  integer :: n
  
  n = size(xa)
  if (size(ya)/=n .or. size(y2a)/=n) then
    call err()
  endif
  
  ! Use bisection to find the indices of array xa that bracket x.
  below=1
  above=n
  do while (above-below > 1)
    k = (above+below)/2
    if (xa(k) > x) then
      above = k
    else
      below = k
    endif
  enddo
  
  ! Determine the finite difference along the X dimension
  h = xa(above)-xa(below)
  if (h <= 0) then
    call err()
  endif

  ! Interpolate.
  a = (xa(above)-x)/h
  b = (x-xa(below))/h
  y = a*ya(below) &
  & + b*ya(above) &
  & + ( (a**3-a)*y2a(below) + (b**3-b)*y2a(above) ) * (h**2)/6.
end function

! ----------------------------------------------------------------------
!     SPLIE2 use: given an array of X1A data of length M, and an array
!     of X2A data of length N, this routine computes the 2nd
!     derivatives, Y2A, at each X1A,X2A data point.  Thus Y2A has
!     dimensions Y2A(M,N).  Natural spline fitting is assumed.

!     this routine called once, prior to using routine SPLIN2, as a set
!     up for using routine SPLIN2, which performs the actual
!     interpolation.

!     Uses routines: SPLINE

!     IMPORTANT NOTE: the X1A and X2A data values must both be in
!     ascending order or the interpolation will fail
! ----------------------------------------------------------------------
function splie2(x1a,x2a,ya) result(y2a)
  implicit none
  
  real(dp), intent(in)  :: x1a(:) ! m
  real(dp), intent(in)  :: x2a(:) ! n
  real(dp), intent(in)  :: ya(:,:) ! m,n
  real(dp), allocatable :: y2a(:,:) ! m,n
  
  integer :: m,n
  integer :: i,ialloc
  
  m = size(x1a)
  n = size(x2a)
  
  if (any(shape(ya)/=[m,n] .or. shape(y2a)/=[m,n])) then
    call err()
  endif
  
  allocate(y2a(m,n), stat=ialloc); call err(ialloc)
  
  do i=1,m
    y2a(i,:) = spline(x2a,ya(i,:),1.0e35_dp,1.0e35_dp)
  enddo
end function

! ----------------------------------------------------------------------
! SPLIN2 use: given an array of X1A data of length M, an array of
!    X2A data of length N, and an array of 2nd derivatives, Y2A at each
!    X1A,X2A data point, dimensioned Y2A(M,N), this routine performs 2D
!    interpolation, returning the interpolated value Y at user input
!    values X1 and X2.  Natural spline fitting is assumed.
!
! IMPORTANT NOTE: the X1A and X2A data values must both be in
!    ascending order or the interpolation will fail
! ----------------------------------------------------------------------
function splin2(x1a,x2a,ya,y2a,x1,x2) result(y)
  implicit none
  
  real(dp), intent(in) :: x1a(:)
  real(dp), intent(in) :: x2a(:)
  real(dp), intent(in) :: ya(:,:)
  real(dp), intent(in) :: y2a(:,:)
  real(dp), intent(in) :: x1
  real(dp), intent(in) :: x2
  real(dp)             :: y
  
  real(dp), allocatable :: yytmp(:)
  real(dp), allocatable :: y3tmp(:)
  
  integer :: m,n
  integer :: i,ialloc
  
  m = size(x1a)
  n = size(x2a)
  
  if (any(shape(ya)/=[m,n] .or. shape(y2a)/=[m,n])) then
    call err()
  endif
  
  allocate(yytmp(m), stat=ialloc); call err(ialloc)
  do i=1,m
    yytmp(i) = splint(x2a,ya(i,:),y2a(i,:),x2)
  enddo
  
  y3tmp = spline(x1a,yytmp,1.0e30_dp,1.0e30_dp)
  y = splint(x1a,yytmp,y3tmp,x1)
end function
end module
