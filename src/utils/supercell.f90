module supercell_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
  
  ! A class to hold supercell data.
  ! n.b. recip_supercell and gvectors are stored in a representation where
  ! The primitive IBZ is a cube with side length sc_size.
  type SupercellData
    ! The number of primitive cells in the supercell.
    integer              :: sc_size
    
    ! The supercell lattice and reciprocal lattice vectors.
    integer              :: supercell(3,3)
    integer              :: recip_supercell(3,3)
    
    ! The G-vectors of the supercell Brillouin Zone which are unique 
    !    in the primitive cell Brillouin Zone.
    ! N.B. this is the same as the R-vectors of the primitive cell which are
    !    unique in the supercell.
    integer, allocatable :: gvectors(:,:)
  end type
  
  interface new
    module procedure new_SupercellData
  end interface
  
contains

subroutine new_SupercellData(this,sc_size)
  implicit none
  
  type(SupercellData), intent(out) :: this
  integer,             intent(in)  :: sc_size
  
  this%sc_size = sc_size
  allocate(this%gvectors(3,sc_size))
end subroutine
end module
