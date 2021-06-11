! ======================================================================
! A storage class, for all the variables which are passed to the
!    various potential methods.
! ======================================================================
module caesar_anharmonic_data_module
  use caesar_common_module
  
  use caesar_subspaces_module
  use caesar_stars_module
  
  use caesar_interpolated_supercell_module
  use caesar_max_displacement_module
  use caesar_degenerate_symmetry_module
  implicit none
  
  private
  
  public :: AnharmonicData
  
  type, extends(Stringsable) :: AnharmonicData
    type(StructureData)                    :: structure
    type(StructureData)                    :: anharmonic_supercell
    type(QpointData),          allocatable :: qpoints(:)
    type(ComplexMode),         allocatable :: complex_modes(:)
    type(RealMode),            allocatable :: real_modes(:)
    type(DegenerateSubspace),  allocatable :: degenerate_subspaces(:)
    type(DegenerateSymmetry),  allocatable :: degenerate_symmetries(:)
    type(SubspaceCoupling),    allocatable :: subspace_couplings(:)
    integer                                :: max_subspace_coupling
    integer                                :: max_qpoint_coupling
    integer                                :: potential_expansion_order
    logical                                :: vscf_basis_functions_only
    type(MaxDisplacement)                  :: max_displacement
    type(Group),               allocatable :: qpoint_symmetry_groups(:)
  contains
    procedure, public :: read  => read_AnharmonicData
    procedure, public :: write => write_AnharmonicData
  end type
  
  interface AnharmonicData
    module function new_AnharmonicData(structure,anharmonic_supercell,   &
       & qpoints,complex_modes,real_modes,degenerate_subspaces,          &
       & degenerate_symmetries,subspace_couplings,max_subspace_coupling, &
       & max_qpoint_coupling,potential_expansion_order,                  &
       & vscf_basis_functions_only,max_displacement) result(this)
      type(StructureData),      intent(in) :: structure
      type(StructureData),      intent(in) :: anharmonic_supercell
      type(QpointData),         intent(in) :: qpoints(:)
      type(ComplexMode),        intent(in) :: complex_modes(:)
      type(RealMode),           intent(in) :: real_modes(:)
      type(DegenerateSubspace), intent(in) :: degenerate_subspaces(:)
      type(DegenerateSymmetry), intent(in) :: degenerate_symmetries(:)
      type(SubspaceCoupling),   intent(in) :: subspace_couplings(:)
      integer,                  intent(in) :: max_subspace_coupling
      integer,                  intent(in) :: max_qpoint_coupling
      integer,                  intent(in) :: potential_expansion_order
      logical,                  intent(in) :: vscf_basis_functions_only
      type(MaxDisplacement),    intent(in) :: max_displacement
      type(AnharmonicData)                 :: this
    end function
  
    module function new_AnharmonicData_data(structure,                        &
       & interpolated_supercell,potential_expansion_order,                    &
       & max_subspace_coupling,max_qpoint_coupling,vscf_basis_functions_only, &
       & max_displacement) result(this) 
      type(StructureData),         intent(in) :: structure
      type(InterpolatedSupercell), intent(in) :: interpolated_supercell
      integer,                     intent(in) :: potential_expansion_order
      integer,                     intent(in) :: max_subspace_coupling
      integer,                     intent(in) :: max_qpoint_coupling
      logical,                     intent(in) :: vscf_basis_functions_only
      type(MaxDisplacement),       intent(in) :: max_displacement
      type(AnharmonicData)                    :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_AnharmonicData(this,input) 
      class(AnharmonicData), intent(out) :: this
      type(String),          intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_AnharmonicData(this) result(output) 
      class(AnharmonicData), intent(in) :: this
      type(String), allocatable         :: output(:)
    end function
  end interface
  
  interface AnharmonicData
    module function new_AnharmonicData_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(AnharmonicData)     :: this
    end function
  
    impure elemental module function new_AnharmonicData_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(AnharmonicData)          :: this
    end function
  end interface
end module
