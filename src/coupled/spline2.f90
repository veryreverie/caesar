module spline2_module
  use string_module
  use io_module
  use constants_module, only : dp
  
  interface spline
    module procedure spline_1d
    module procedure spline_2d
  end interface
  
  interface splint
    module procedure splint_1d
    module procedure splint_2d
  end interface
contains

subroutine spline2()
  use ifile_module
  use ofile_module
  implicit none
  
  ! Independent data points.
  real(dp),allocatable :: xb(:),fb(:)
  real(dp),allocatable :: xc(:),fc(:)
  ! Coupled data points.
  real(dp),allocatable :: x1a(:),x2a(:),f(:,:)
  ! Second derivatives.
  real(dp),allocatable :: d2f(:,:),d2fb(:),d2fc(:)
  ! Interpolated data.
  real(dp) :: u1,u2,f_int
  real(dp) :: dx1,dx2
  integer :: mode1,mode2
  
  type(IFile) :: input_file
  type(IFile) :: modes_file
  type(IFile) :: mode1_file
  type(IFile) :: mode2_file
  type(IFile) :: energy_file
  type(OFile) :: output_file
  
  ! Temporary variables.
  integer                   :: n,m,p
  integer                   :: i,j,ialloc
  type(string), allocatable :: line(:)
  
  ! Read in number of integration points.
  input_file = 'fit_input.dat'
  line = split(input_file%line(1))
  n = int(line(1))
  m = int(line(2))
  p = int(line(3))
  
  ! Allocate arrays.
  allocate( x1a(n),  &
          & x2a(n),  &
          & f(n,n),  &
          & d2f(n,n), &
          & xb(p),  &
          & xc(p),  &
          & fb(p),   &
          & fc(p),   &
          & d2fb(p),  &
          & d2fc(p),  &
          & stat=ialloc); call err(ialloc)

  
  ! Read in input data.
  modes_file = 'modes.dat'
  line = split(modes_file%line(1))
  mode1 = int(line(1))
  mode2 = int(line(2))
  
  mode1_file = 'fit_energy_'//mode1//'.dat'
  do i=1,p
    line = split(mode1_file%line(i))
    xb(i) = dble(line(1))
    fb(i)  = dble(line(2))
  enddo
  
  mode2_file = 'fit_energy_'//mode2//'.dat'
  do i=1,p
    line = split(mode2_file%line(i))
    xc(i) = dble(line(1))
    fc(i)  = dble(line(2))
  enddo
  
  energy_file = 'fit_energy.dat'
  do i=1,n
    do j=1,n
      if(i/=((n+1)/2).and.j/=((n+1)/2))then 
        line = split(energy_file%line((i-1)*n+j))
        x1a(i) = dble(line(1))
        x2a(j) = dble(line(2))
        f(i,j) = dble(line(3))
      endif
    enddo
  enddo
  
  x1a((n+1)/2) = 0
  x2a((n+1)/2) = 0
  f((n+1)/2,:) = 0
  f(:,(n+1)/2) = 0
  
  ! Run calculation.
  d2fb = spline(xb,fb)
  d2fc = spline(xc,fc)
  do i=1,n
    do j=1,n
      if(i/=((n+1)/2).and.j/=((n+1)/2))then
        f(i,j)=f(i,j)-splint(xb,fb,d2fb,x1a(i))-splint(xc,fc,d2fc,x2a(j))
      endif
    enddo
  enddo
  
  d2f = spline(x1a,x2a,f)
  
  dx1=2*abs(x1a(1))/m
  dx2=2*abs(x2a(1))/m
  do i=1,m+1
    u1 = x1a(1) + (i-1)*dx1
    do j=1,m+1
      u2 = x2a(1) + (i-1)*dx2
      f_int = splint(x1a,x2a,f,d2f,u1,u2)
      call output_file%print_line(u1//' '//u2//' '//f_int)
    enddo
  enddo
end subroutine

! ----------------------------------------------------------------------
! Given x and y, both of length n, spline computes d2y/dx2 at each x.
! Called once, to set up for splint, which performs the actual interpolation.
! x must be in ascending order.
! ----------------------------------------------------------------------
function spline_1d(x,y) result(output)
  implicit none
  
  real(dp), intent(in)  :: x(:)
  real(dp), intent(in)  :: y(:)
  real(dp), allocatable :: output(:)
  
  integer :: n
  
  real(dp)              :: sig
  real(dp)              :: p
  real(dp), allocatable :: u(:)
  
  integer :: i,ialloc
  
  n = size(x)
  if (size(y)/=n) then
    call err()
  endif
  
  allocate( output(n), &
          & u(n-1),    &
          & stat=ialloc); call err(ialloc)
  
  ! Store intermediate values of terms in the expansion series.
  u(1) = 0
  do i=2,n-1
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*output(i-1)+2
    u(i) = ( ( (y(i+1)-y(i))/(x(i+1)-x(i)) &
       &     - (y(i)-y(i-1))/(x(i)-x(i-1)) &
       &     )                             &
       &     * 6 / (x(i+1)-x(i-1))         &
       &     - sig*u(i-1)                  &
       &   ) / p
    output(i) = (sig-1) / p
  enddo
  
  ! Compute the derivitives from the 2nd order expansion series.
  output(n) = 0
  do i=n-1,2,-1
    output(i) = output(i)*output(i+1) + u(i)
  enddo
  output(0) = 0
end function

! ----------------------------------------------------------------------
! Given x1 of length m, x2a of length n, and y at each point in x1*x2,
!    calculates d2y/dx2 at each point.
! Called once to set up for splint, which performs the actual interpolation.
! x1 and x2 must be in ascending order.
! ----------------------------------------------------------------------
function spline_2d(x1,x2,y) result(output)
  implicit none
  
  real(dp), intent(in)  :: x1(:) ! m
  real(dp), intent(in)  :: x2(:) ! n
  real(dp), intent(in)  :: y(:,:) ! m,n
  real(dp), allocatable :: output(:,:) ! m,n
  
  integer :: m,n
  integer :: i,ialloc
  
  m = size(x1)
  n = size(x2)
  
  if (any(shape(y)/=[m,n])) then
    call err()
  endif
  
  allocate(output(m,n), stat=ialloc); call err(ialloc)
  do i=1,m
    output(i,:) = spline(x2,y(i,:))
  enddo
end function

! ----------------------------------------------------------------------
! Given x, y(x), and d2y(x), performs cubic interpolation.
! Returns y(u).
! x must be in ascending order.
! ----------------------------------------------------------------------
function splint_1d(x,y,d2y,u) result(output)
  implicit none
  
  real(dp), intent(in) :: x(:)
  real(dp), intent(in) :: y(:)
  real(dp), intent(in) :: d2y(:)
  real(dp), intent(in) :: u
  real(dp)             :: output
  
  integer :: below
  integer :: above
  
  real(dp) :: h
  real(dp) :: a
  real(dp) :: b
  
  integer :: n
  integer :: i
  
  n = size(x)
  if (size(y)/=n .or. size(d2y)/=n) then
    call err()
  endif
  
  ! Use bisection to find the indices of array x that bracket u.
  below=1
  above=n
  do while (above-below > 1)
    i = (above+below)/2
    if (x(i) > u) then
      above = i
    else
      below = i
    endif
  enddo
  
  ! Determine the finite difference along the X dimension
  h = x(above)-x(below)
  if (h <= 0) then
    call err()
  endif

  ! Interpolate.
  a = (x(above)-u)/h
  b = (u-x(below))/h
  output = a*y(below) &
       & + b*y(above) &
       & + ( (a**3-a)*d2y(below) + (b**3-b)*d2y(above) ) * (h**2)/6.
end function

! ----------------------------------------------------------------------
! Given x1, x2, y(x1,x2) and d2y(x1,x2),
!    performs 2-d cubic interpolation.
! Returns y(x1,x2)
! x1 and x2 must be in ascending order.
! ----------------------------------------------------------------------
function splint_2d(x1,x2,y,d2y,u1,u2) result(output)
  implicit none
  
  real(dp), intent(in) :: x1(:)
  real(dp), intent(in) :: x2(:)
  real(dp), intent(in) :: y(:,:)
  real(dp), intent(in) :: d2y(:,:)
  real(dp), intent(in) :: u1
  real(dp), intent(in) :: u2
  real(dp)             :: output
  
  real(dp), allocatable :: y_temp(:)
  
  integer :: m,n
  integer :: i,ialloc
  
  m = size(x1)
  n = size(x2)
  
  if (any(shape(y)/=[m,n] .or. shape(d2y)/=[m,n])) then
    call err()
  endif
  
  allocate(y_temp(m), stat=ialloc); call err(ialloc)
  do i=1,m
    y_temp(i) = splint(x2,y(i,:),d2y(i,:),u2)
  enddo
  
  output = splint(x1,y_temp,spline(x1,y_temp),u1)
end function
end module
