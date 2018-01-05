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
    type(IntVector), allocatable :: image_rvectors(:)
  end type
  
  interface size
    module procedure size_MinImages
  end interface
contains

function size_MinImages(this) result(output)
  implicit none
  
  type(MinImages) :: this
  integer         :: output
  
  output = size(this%image_rvectors)
end function

! ----------------------------------------------------------------------
! Computes the set of R-vectors such that the cartesian displacement
!    between a+R and b is minimised.
! Returns multiple R-vectors if they have similar lengths.
! ----------------------------------------------------------------------
function min_images_brute_force(supercell,a,b) result(output)
  use linear_algebra_module
  use structure_module
  implicit none
  
  type(StructureData), intent(in) :: supercell
  integer,             intent(in) :: a
  integer,             intent(in) :: b
  type(MinImages)                 :: output
  
  ! Parameters.
  ! Number of "shells" of supercell R-vectors to check.
  integer,parameter :: check_shell=3
  ! Maximum number of images.
  integer, parameter :: max_images = 8
  
  ! Positions in fractional supercell co-ordinates.
  type(RealVector) :: a_frac_pos
  type(RealVector) :: b_frac_pos
  
  ! The R-vector from a to b ignoring supercell R-vectors.
  type(IntVector)  :: rvector
  
  ! Tolerance parameters.
  real(dp), parameter :: tol=1.0e-8_dp
  real(dp)            :: tol_l2
  
  ! Distances in cartesian co-ordinates.
  type(RealVector) :: displacement
  real(dp)         :: l2_distance
  real(dp)         :: min_l2_distance
  
  ! Output data.
  integer                      :: no_images
  type(IntVector), allocatable :: image_rvectors(:)
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  rvector = supercell%rvectors(supercell%atoms(b)%rvec_id()) &
        & - supercell%rvectors(supercell%atoms(a)%rvec_id())
  
  tol_l2 = tol*supercell%volume**(2/3.0_dp)
  
  ! Get the fractional position of both atoms, translated into the
  !    supercell's first unit cell.
  a_frac_pos = supercell%atoms(a)%fractional_position()
  a_frac_pos = a_frac_pos-vec(floor(dble(a_frac_pos)))
  
  b_frac_pos = supercell%atoms(b)%fractional_position()
  b_frac_pos = b_frac_pos-vec(floor(dble(b_frac_pos)))
  
  ! Identify R-vectors.
  allocate(image_rvectors(max_images), stat=ialloc); call err(ialloc)
  min_l2_distance = huge(1.0_dp)
  no_images = 0
  do i=-check_shell,check_shell
    do j=-check_shell,check_shell
      do k=-check_shell,check_shell
        displacement = transpose(supercell%lattice) &
                   & * (b_frac_pos-a_frac_pos+vec([i,j,k]))
        l2_distance = displacement * displacement
        if (l2_distance < min_l2_distance-tol_l2) then
          ! The new R-vector makes a+R and b closer than any previous R-vector.
          min_l2_distance = l2_distance
          no_images = 1
          image_rvectors(no_images) = vec([i,j,k])
        elseif (l2_distance <= min_l2_distance+tol_l2) then
          ! The new R-vector makes a+R and b similarly close to
          !    previous R-vectors.
          no_images = no_images+1
          if (no_images>max_images) then
            call print_line(ERROR//': maxim too small.')
            call err()
          endif
          
          ! Convert the supercell R-vector from fractional supercell co-ords
          !    to fractional primitive co-ords,
          !    and add on the primitive R-vector.
          image_rvectors(no_images) = transpose(supercell%supercell) &
                                  & * vec([i,j,k])                   &
                                  & + rvector
        endif
      enddo
    enddo
  enddo
  
  if (no_images==0) then
    call print_line(ERROR//': bug in min_images_brute_force.')
    call err()
  endif
  
  output = MinImages(image_rvectors(:no_images))
end function
end module
