! ======================================================================
! All data relating to a given atomic configuration.
! ======================================================================
module caesar_structure_data_module
  use caesar_utils_module
  
  use caesar_atom_module
  use caesar_spglib_module
  
  use caesar_basic_symmetry_module
  use caesar_symmetry_module
  implicit none
  
  private
  
  public :: StructureData
  public :: BasicStructure
  
  ! ----------------------------------------------------------------------
  ! The structure type.
  ! ----------------------------------------------------------------------
  type, extends(Stringsable) :: StructureData
    ! ------------------------------
    ! Lattice data
    ! ------------------------------
    type(RealMatrix) :: lattice
    type(RealMatrix) :: recip_lattice
    real(dp)         :: volume
    
    ! ------------------------------
    ! Atom data
    ! ------------------------------
    integer :: no_atoms
    integer :: no_atoms_prim
    integer :: no_modes
    integer :: no_modes_prim
    
    type(AtomData), allocatable :: atoms(:)
    
    ! ------------------------------
    ! Symmetry data (in fractional co-ordinates).
    ! ------------------------------
    real(dp)                            :: symmetry_precision
    type(String)                        :: space_group
    type(SymmetryOperator), allocatable :: symmetries(:)
    
    ! ------------------------------
    ! Superell data
    ! ------------------------------
    ! The number of primitive cells in the supercell.
    integer :: sc_size
    
    ! The lattice vectors of the supercell,
    !    in fractional primitive cell co-ordinates.
    type(IntMatrix) :: supercell
    
    ! invert(transpose(supercell))
    type(FractionMatrix) :: recip_supercell
    
    ! The R-vectors of the primitive cell which are not related by supercell
    !    lattice vectors.
    type(IntVector), allocatable :: rvectors(:)
    ! The G-vectors of the reciprocal supercell which are not related by
    !    primitive reciprocal cell vectors.
    type(IntVector), allocatable :: gvectors(:)
    
    ! --------------------------------------------------
    ! The IDs of paired R-vectors and G-vectors.
    ! Used to provide inverse and pair functions.
    ! --------------------------------------------------
    ! The ID of the R-vector, j, s.t. rvectors(i) + rvectors(j) = 0.
    integer, allocatable :: rvector_paired_ids_(:)
    ! The ID of the G-vector, j, s.t. gvectors(i) + gvectors(j) = 0.
    integer, allocatable :: gvector_paired_ids_(:)
  contains
    ! Snaps the structure to symmetry.
    procedure, public :: snap_to_symmetry => snap_to_symmetry_StructureData
    
    ! Calculate the symmetry operators of the structure.
    procedure, public :: calculate_symmetry
    
    ! Return inverse symmetries or paired R-vectors / G-vectors.
    procedure, public :: paired_rvector_id
    procedure, public :: paired_rvector
    procedure, public :: paired_gvectors
    
    ! Primitive lattice procedures.
    procedure, public :: prim_lattice
    procedure, public :: prim_recip_lattice
    procedure, public :: prim_volume
    
    ! I/O.
    procedure, public :: read  => read_StructureData
    procedure, public :: write => write_StructureData
  end type
  
  interface
    ! ----------------------------------------------------------------------
    ! Functions to provide paired R-vectors and G-vectors,
    !    and the groups corresponding to those objects.
    ! ----------------------------------------------------------------------
    ! Returns the pair of an R-vector.
    ! i.e. paired_rvectors(i) + rvectors(i) = 0, modulo supercell R-vectors.
    module function paired_rvector_id(this,rvector_id) result(output) 
      class(StructureData), intent(in) :: this
      integer,              intent(in) :: rvector_id
      integer                          :: output
    end function
  end interface
  
  interface
    module function paired_rvector(this,rvector_id) result(output) 
      class(StructureData), intent(in) :: this
      integer,              intent(in) :: rvector_id
      type(IntVector)                  :: output
    end function
  end interface
  
  interface
    ! Returns the pair of an G-vector.
    ! i.e. paired_gvectors(i) + gvectors(i) = 0, modulo primitive G-vectors.
    module function paired_gvectors(this,gvector_id) result(output) 
      class(StructureData), intent(in) :: this
      integer,              intent(in) :: gvector_id
      type(IntVector)                  :: output
    end function
  end interface
  
  interface StructureData
    ! ----------------------------------------------------------------------
    ! Allocates all arrays, and sets no_ variables.
    ! ----------------------------------------------------------------------
    module function new_StructureData(basic_structure,basic_supercell, &
       & recip_lattice_matrix,recip_supercell_matrix) result(this) 
      type(BasicStructure), intent(in)           :: basic_structure
      type(BasicSupercell), intent(in), optional :: basic_supercell
      type(RealMatrix),     intent(in), optional :: recip_lattice_matrix
      type(FractionMatrix), intent(in), optional :: recip_supercell_matrix
      type(StructureData)                        :: this
    end function
  end interface
  
  interface BasicStructure
    ! ----------------------------------------------------------------------
    ! Convert to a BasicStructure.
    ! ----------------------------------------------------------------------
    impure elemental module function new_BasicStructure_StructureData(this) &
       & result(output) 
      type(StructureData), intent(in) :: this
      type(BasicStructure)            :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Snaps the structure to symmetry.
    ! ----------------------------------------------------------------------
    module function snap_to_symmetry_StructureData(this,symmetry_precision) &
       & result(output) 
      class(StructureData), intent(in)           :: this
      real(dp),             intent(in), optional :: symmetry_precision
      type(StructureData)                        :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Calculate symmetries.
    ! ----------------------------------------------------------------------
    module subroutine calculate_symmetry(this,symmetry_precision,symmetries, &
       & loto_direction) 
      class(StructureData), intent(inout)        :: this
      real(dp),             intent(in)           :: symmetry_precision
      type(BasicSymmetry),  intent(in), optional :: symmetries(:)
      type(FractionVector), intent(in), optional :: loto_direction
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Calculate properties of the primitive lattice.
    ! ----------------------------------------------------------------------
    module function prim_lattice(this) result(output) 
      class(StructureData), intent(in) :: this
      type(RealMatrix)                 :: output
    end function
  end interface
  
  interface
    module function prim_recip_lattice(this) result(output) 
      class(StructureData), intent(in) :: this
      type(RealMatrix)                 :: output
    end function
  end interface
  
  interface
    module function prim_volume(this) result(output) 
      class(StructureData), intent(in) :: this
      real(dp)                         :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_StructureData(this,input) 
      class(StructureData), intent(out) :: this
      type(String),         intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_StructureData(this) result(output) 
      class(StructureData), intent(in) :: this
      type(String), allocatable        :: output(:)
    end function
  end interface
  
  interface StructureData
    module function new_StructureData_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(StructureData)      :: this
    end function
  
    impure elemental module function new_StructureData_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(StructureData)           :: this
    end function
  end interface
end module
