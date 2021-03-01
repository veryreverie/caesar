! ======================================================================
! A minimal representation of the Structure class.
! ======================================================================
module caesar_basic_structure_module
  use caesar_utils_module
  
  use caesar_basic_atom_module
  implicit none
  
  private
  
  public :: BasicStructure
  public :: BasicSupercell
  
  type, extends(NoDefaultConstructor) :: BasicStructure
    type(RealMatrix)             :: lattice_matrix
    type(BasicAtom), allocatable :: atoms(:)
  contains
    procedure, public :: volume => volume_BasicStructure
  end type
  
  type, extends(NoDefaultConstructor) :: BasicSupercell
    type(IntMatrix)              :: supercell_matrix
    type(IntVector), allocatable :: rvectors(:)
    type(IntVector), allocatable :: gvectors(:)
    integer,         allocatable :: atom_rvector_ids(:)
    integer,         allocatable :: atom_prim_ids(:)
  end type
  
  interface BasicStructure
    module function new_BasicStructure(lattice_matrix,atoms) result(this) 
      type(RealMatrix), intent(in) :: lattice_matrix
      type(BasicAtom),  intent(in) :: atoms(:)
      type(BasicStructure)         :: this
    end function
  
    module function new_BasicStructure_constructed(lattice_matrix,species, &
       & masses,cartesian_positions) result(this) 
      type(RealMatrix), intent(in) :: lattice_matrix
      type(String),     intent(in) :: species(:)
      real(dp),         intent(in) :: masses(:)
      type(RealVector), intent(in) :: cartesian_positions(:)
      type(BasicStructure)         :: this
    end function
  end interface
  
  interface BasicSupercell
    module function new_BasicSupercell(supercell_matrix,rvectors,gvectors, &
       & atom_rvector_ids,atom_prim_ids) result(output) 
      type(IntMatrix), intent(in) :: supercell_matrix
      type(IntVector), intent(in) :: rvectors(:)
      type(IntVector), intent(in) :: gvectors(:)
      integer,         intent(in) :: atom_rvector_ids(:)
      integer,         intent(in) :: atom_prim_ids(:)
      type(BasicSupercell)        :: output
    end function
  end interface
  
  interface
    impure elemental module function volume_BasicStructure(this) &
       & result(output) 
      class(BasicStructure), intent(in) :: this
      real(dp)                          :: output
    end function
  end interface
end module
