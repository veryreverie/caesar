! Subroutines for the calculation of minimum-image distances.
module min_images_module
  use constants, only : dp
  implicit none
  
  type MinImages
    real(dp), allocatable :: images(:,:)
  end type
  
  interface new
    module procedure new_MinImages
  end interface
  
  interface drop
    module procedure drop_MinImages
  end interface
  
  interface size
    module procedure size_MinImages
  end interface
contains

subroutine new_MinImages(this,no_images)
  implicit none
  
  type(MinImages), intent(out) :: this
  integer,         intent(in)  :: no_images
  
  allocate(this%images(3,no_images))
end subroutine

subroutine drop_MinImages(this)
  implicit none
  
  type(MinImages), intent(inout) :: this
  
  deallocate(this%images)
end subroutine

function size_MinImages(this) result(output)
  implicit none
  
  type(MinImages) :: this
  integer         :: output
  
  output = size(this%images,2)
end function

! Compute the minimum image vector(s) b of vector a with respect to the 
! lattice specified by the columns of lat_vec. rec_vec are the reciprocal 
! lattice vectors (w/o 2pi). -b is the vector from a to its closest lattice 
! point.  nim is the number of image vectors.
function min_images_brute_force(a,structure) result(output)
  use constants,      only : dp
  use linear_algebra, only : invert
  use string_module
  use structure_module
  implicit none
  
  real(dp),            intent(in) :: a(3)
  type(StructureData), intent(in) :: structure
  type(MinImages)                 :: output
  
  real(dp) :: b(3,8)
  integer  :: nim
  
  real(dp) :: mag_b_sq,dist2,tol_l2
  real(dp) :: delta(3)
  integer :: n(3),i,j,k
  
  ! Number of "shells" of lattice points to check.  Only used in setup, so
  ! may as well overkill.
  integer,parameter :: check_shell=3
  real(dp),parameter :: tol=1.d-8
  ! Maximum number of images
  integer, parameter :: maxim = 8
  
  tol_L2 = tol*dot_product(structure%lattice(:,1),structure%lattice(:,1))
  n = floor(matmul(structure%recip_lattice,a))
  
  nim = 0
  mag_b_sq = -1.0_dp
  do i=n(1)-check_shell,n(1)+check_shell+1
    do j=n(2)-check_shell,n(2)+check_shell+1
      do k=n(3)-check_shell,n(3)+check_shell+1
        delta = a-matmul(transpose(structure%lattice),(/i,j,k/))
        dist2 = dot_product(delta,delta)
        if(nim/=0 .and. abs(dist2-mag_b_sq)<=tol_L2)then
          nim = nim+1
          if (nim>maxim) then
            call print_line('Error: min_images_brute_force: maxim too small.')
            stop
          endif
          b(:,nim) = delta
        elseif(dist2<mag_b_sq.or.nim==0)then
          mag_b_sq = dist2
          nim = 1
          b(:,1) = delta
        endif
      enddo
    enddo
  enddo
  
  if (nim==0) then
    call print_line('Error: bug in min_images_brute_force.')
    stop
  endif
  
  call new(output,nim)
  output%images = b(:,1:nim)
end function
end module
