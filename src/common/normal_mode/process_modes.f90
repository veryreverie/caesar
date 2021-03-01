! ======================================================================
! Take the modes at a given q-point, and process them w/r/t the q-point,
!    structure and symmetry data.
! ======================================================================
! Fills in:
!    - id
!    - paired_id
!    - qpoint_id
!    - paired_qpoint_id
!    - subspace_id
! In general, a the modes in a given subspace will be returned as a linear
!    combination of one another.
module caesar_process_modes_module
  use caesar_utils_module
  
  use caesar_structure_module
  
  use caesar_complex_mode_module
  use caesar_complex_mode_symmetry_module
  implicit none
  
  private
  
  public :: process_modes
  
  ! Types used for splitting degenerate modes. The type SplitModes is for
  !    modes split by commuting symmetries, and the type AntiSplitModes is for
  !    modes split by anti-commuting symmetries.
  type, extends(NoDefaultConstructor) :: SplitModes
    type(ComplexMode), allocatable :: modes(:)
    integer,           allocatable :: phases(:)
  end type
  
  type, extends(NoDefaultConstructor) :: AntiSplitModes
    type(ComplexMode), allocatable :: modes(:)
  end type
  
  interface SplitModes
    ! Constructors.
    module function new_SplitModes(modes,phases) result(this) 
      type(ComplexMode), intent(in) :: modes(:)
      integer,           intent(in) :: phases(:)
      type(SplitModes)              :: this
    end function
  end interface
  
  interface AntiSplitModes
    module function new_AntiSplitModes(modes) result(this) 
      type(ComplexMode), intent(in) :: modes(:)
      type(AntiSplitModes)          :: this
    end function
  end interface
  
  interface
    ! Process modes.
    module function process_modes(input,structure,qpoint,subspace_id) &
       & result(output) 
      type(ComplexMode),   intent(in) :: input(:)
      type(StructureData), intent(in) :: structure
      type(QpointData),    intent(in) :: qpoint
      integer,             intent(in) :: subspace_id
      type(ComplexMode), allocatable  :: output(:)
    end function
  end interface
  
  interface
    ! --------------------------------------------------
    ! Chooses the basis for the degenerate subspace using symmetry operators.
    ! --------------------------------------------------
    ! Symmetries must all take the q-point to itself.
    module function choose_basis_complex(input,structure,symmetries,qpoint) &
       & result(output) 
      type(ComplexMode),      intent(in)    :: input(:)
      type(StructureData),    intent(in)    :: structure
      type(SymmetryOperator), intent(in)    :: symmetries(:)
      type(QpointData),       intent(in)    :: qpoint
      type(ComplexMode), allocatable        :: output(:)
    end function
  end interface
  
  interface
    module function choose_basis_real(input,structure,symmetries,qpoint) &
       & result(output) 
      type(ComplexMode),      intent(in)    :: input(:)
      type(StructureData),    intent(in)    :: structure
      type(SymmetryOperator), intent(in)    :: symmetries(:)
      type(QpointData),       intent(in)    :: qpoint
      type(ComplexMode), allocatable        :: output(:)
    end function
  end interface
  
  interface SplitModes
    module function new_SplitModes_ComplexModes(input,symmetry,qpoint, &
       & positive_superposition) result(this) 
      type(ComplexMode),      intent(in)    :: input(:)
      type(SymmetryOperator), intent(in)    :: symmetry
      type(QpointData),       intent(in)    :: qpoint
      logical,                intent(in)    :: positive_superposition
      type(SplitModes)                      :: this
    end function
  end interface
  
  interface AntiSplitModes
    module function new_AntiSplitModes_ComplexModes(input,symmetries,qpoint) &
       & result(this) 
      type(ComplexMode),      intent(in) :: input(:)
      type(SymmetryOperator), intent(in) :: symmetries(:)
      type(QpointData),       intent(in) :: qpoint
      type(AntiSplitModes)               :: this
    end function
  end interface
end module
