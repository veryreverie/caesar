! ======================================================================
! Subroutines for the calculation of minimum-image distances.
! ======================================================================
module caesar_min_images_module
  use caesar_utils_module
  
  use caesar_structure_module
  implicit none
  
  private
  
  public :: MinImages
  public :: size
  public :: calculate_min_images
  
  type MinImages
    type(IntVector), allocatable :: image_rvectors(:)
  end type
  
  interface size
    module function size_MinImages(this) result(output) 
      type(MinImages) :: this
      integer         :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns min_images(a,b).
    ! Computes the set of supercell R-vectors such that the cartesian displacement
    !    between a+R and b is minimised.
    ! Returns multiple R-vectors if they have similar lengths.
    ! Output given in fractional primitive co-ordinates, so R-vectors are integers.
    ! ----------------------------------------------------------------------
    module function calculate_min_images(supercell) result(output) 
      type(StructureData), intent(in) :: supercell
      type(MinImages), allocatable    :: output(:,:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Computes the set of R-vectors such that the cartesian displacement
    !    between a+R and b is minimised.
    ! Returns multiple R-vectors if they have similar lengths.
    ! ----------------------------------------------------------------------
    module function min_images_brute_force(supercell,a,b) result(output) 
      type(StructureData), intent(in) :: supercell
      integer,             intent(in) :: a
      integer,             intent(in) :: b
      type(MinImages)                 :: output
    end function
  end interface
end module
