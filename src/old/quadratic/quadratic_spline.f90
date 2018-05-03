!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Program file name: spline.f90                                          !
!                                                                         !
!  © Tao Pang 2006                                                        !
!                                                                         !
!  Last modified: August 1, 2012                                          !
!                                                                         !
!  (1) This F90 program is created for the book, "An Introduction to      !
!      Computational Physics, 2nd Edition," written by Tao Pang and       !
!      published by Cambridge University Press on January 19, 2006.       !
!                                                                         !
!  (2) No warranties, express or implied, are made for this program.      !
!                                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module quadratic_spline_module
  use common_module
  implicit none
contains

! An example of creating cubic-spline approximation of
! a discrete function fi=f(xi).
pure function quadratic_spline(m,xifi) result(output)
  implicit none
  
  integer,               intent(in) :: m           ! no. points to spline to
  real(dp), allocatable, intent(in) :: xifi(:,:)   ! {(x,f)} at n points
  real(dp), allocatable             :: output(:,:) ! {(x,f)} at m points
  
  ! Input data
  integer               :: n     ! = size(xi) = size(xifi,2)
  real(dp), allocatable :: xi(:) ! = xifi(1,:)
  real(dp), allocatable :: fi(:) ! = xifi(2,:)
  
  integer :: p
  integer :: i, k
  real(dp) :: x, f, dx, h, alpha, beta, gamma, eta
  real(dp),allocatable  :: p2(:)

  n = size(xifi,2)
  
  allocate(xi(n))
  allocate(fi(n))
  allocate(p2(n))
  allocate(output(2,m))
  
  xi = xifi(1,:)
  fi = xifi(2,:)

  p2 = cubic_spline(n-1, xi, fi)

  ! find the approximation of the function
  
  output(1,1) = xi(1)
  output(2,1) = fi(1)
  
  h = (xi(n)-xi(1))/m 
  x = xi(1)
  
  p=m-1
  
  do i = 1, p
    x = x + h
    ! find the interval that x resides
    k = 1
    dx = x-xi(1)
    do while (dx .ge. 0)
      k = k + 1
      dx = x-xi(k)
    end do
    k = k - 1
  
    ! find the value of function f(x)
    dx = xi(k+1) - xi(k)
    alpha = p2(k+1)/(6*dx)
    beta = -p2(k)/(6*dx)
    gamma = fi(k+1)/dx - dx*p2(k+1)/6
    eta = dx*p2(k)/6 - fi(k)/dx
    f = alpha*(x-xi(k))*(x-xi(k))*(x-xi(k)) &
       +beta*(x-xi(k+1))*(x-xi(k+1))*(x-xi(k+1)) &
       +gamma*(x-xi(k))+eta*(x-xi(k+1))
    output(1,i+1) = x
    output(2,i+1) = f
  end do
end function

! Function to carry out the cubic-spline approximation
! with the second-order derivatives returned.
pure function cubic_spline (n, xi, fi) result(p2)
  implicit none
  
  integer,  intent(in) :: n
  real(dp), intent(in) :: xi(n+1)
  real(dp), intent(in) :: fi(n+1)
  real(dp)             :: p2(n+1)
  
  integer  :: i
  real(dp) :: g(n), h(n)
  real(dp) :: d(n-1), b(n-1), c(n-1)
  
  ! Assign the intervals and function differences
  do i = 1, n
    h(i) = xi(i+1) - xi(i)
    g(i) = fi(i+1) - fi(i)
  end do
  
  ! Evaluate the coefficient matrix elements
  do i = 1, n-1
    d(i) = 2*(h(i+1)+h(i))
    b(i) = 6*(g(i+1)/h(i+1)-g(i)/h(i))
    c(i) = h(i+1)
  end do
  
  ! Obtain the second-order derivatives
  g = tridiagonal_linear_eq (n-1, d, c, c, b)
  p2(1) = 0
  p2(n+1) = 0
  do i = 2, n 
    p2(i) = g(i-1)
  end do
end function

! Function to solve the tridiagonal linear equation set.
pure function tridiagonal_linear_eq (l, d, e, c, b) result(z)
  implicit none
  
  integer,  intent(in) :: l
  real(dp), intent(in) :: d(l), e(l), c(l), b(l)
  real(dp)             :: z(l)
  
  integer  :: i
  real(dp) :: y(l), w(l)
  real(dp) :: v(l-1), t(l-1)
  
  ! Evaluate the elements in the LU decomposition
  w(1) = d(1)
  v(1)  = c(1)
  t(1)  = e(1)/w(1)
  do i = 2, l - 1
    w(i) = d(i)-v(i-1)*t(i-1)
    v(i) = c(i)
    t(i) = e(i)/w(i)
  end do
  w(l) = d(l)-v(l-1)*t(l-1)
  
  ! Forward substitution to obtain y
  y(1) = b(1)/w(1)
  do i = 2, l
    y(i) = (b(i)-v(i-1)*y(i-1))/w(i)
  end do
  
  ! Backward substitution to obtain z
  z(l) = y(l)
  do i = l-1, 1, -1
    z(i) = y(i) - t(i)*z(i+1)
  end do
end function
end module
