! ======================================================================
! Subroutines for the calculation of minimum-image distances.
! ======================================================================
module min_images_module
  use constants_module, only : dp
  use string_module
  use io_module
  use linear_algebra_module
  implicit none
  
  type MinImages
    type(RealVector), allocatable :: images(:)
  end type
  
  interface MinImages
    module procedure new_MinImages
  end interface
  
  interface size
    module procedure size_MinImages
  end interface
contains

function new_MinImages(no_images) result(this)
  implicit none
  
  integer, intent(in) :: no_images
  type(MinImages)     :: this
  
  integer :: ialloc
  
  allocate(this%images(no_images), stat=ialloc); call err(ialloc)
end function

function size_MinImages(this) result(output)
  implicit none
  
  type(MinImages) :: this
  integer         :: output
  
  output = size(this%images)
end function

! Computes the minimum distance of the vector 'a' modulo the lattice vectors
!    of the supplied structure.
! Returns multiple vectors if they have similar lengths.
! 'a' and outputs are given in fractional lattice co-ordinates.
function min_images_brute_force(a,structure) result(output)
  use linear_algebra_module
  use structure_module
  implicit none
  
  type(RealVector),    intent(in) :: a
  type(StructureData), intent(in) :: structure
  type(MinImages)                 :: output
  
  type(RealVector) :: b(8)
  integer  :: nim
  
  real(dp)         :: mag_b_sq, dist2, tol_l2
  type(RealVector) :: delta
  integer          :: n(3)
  
  integer :: i,j,k
  
  ! Number of "shells" of lattice points to check.  Only used in setup, so
  ! may as well overkill.
  integer,parameter :: check_shell=3
  real(dp),parameter :: tol=1.0e-8_dp
  ! Maximum number of images
  integer, parameter :: maxim = 8
  
  tol_L2 = tol*l2_norm(structure%lattice*vec([1,0,0]))**2
  n = floor(dble(a))
  
  nim = 0
  mag_b_sq = -1.0_dp
  do i=n(1)-check_shell,n(1)+check_shell+1
    do j=n(2)-check_shell,n(2)+check_shell+1
      do k=n(3)-check_shell,n(3)+check_shell+1
        delta = transpose(structure%lattice) * (a-vec([i,j,k]))
        dist2 = delta * delta
        if(nim/=0 .and. abs(dist2-mag_b_sq)<=tol_L2)then
          nim = nim+1
          if (nim>maxim) then
            call print_line('Error: min_images_brute_force: maxim too small.')
            call err()
          endif
          b(nim) = a-vec([i,j,k])
        elseif(dist2<mag_b_sq.or.nim==0)then
          mag_b_sq = dist2
          nim = 1
          b(1) = a-vec([i,j,k])
        endif
      enddo
    enddo
  enddo
  
  if (nim==0) then
    call print_line('Error: bug in min_images_brute_force.')
    call err()
  endif
  
  output = MinImages(nim)
  output%images = b(1:nim)
end function
end module
