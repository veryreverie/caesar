! ======================================================================
! Harmonic normal modes, in complex co-ordinates.
! ======================================================================
module caesar_complex_mode_module
  use caesar_utils_module
  
  use caesar_structure_module
  implicit none
  
  private
  
  public :: ComplexMode
  public :: conjg
  public :: transform
  public :: select_qpoint
  public :: select_qpoints
  public :: generate_translational_modes
  
  ! A normal mode in complex co-ordinates.
  type, extends(Stringsable) :: ComplexMode
    ! An id which is unique to each mode, and the unique id of the equivalent
    !    mode at the q-point q' s.t. q+q' is a G-vector.
    ! If 2q is a G-vector, then paired_id=id.
    integer :: id
    integer :: paired_id
    
    ! The frequency, and frequency-relevant information.
    real(dp) :: frequency
    real(dp) :: spring_constant    ! k from V(u) = k*u*u/2.
    logical  :: soft_mode          ! True if frequency < -epsilon.
    logical  :: translational_mode ! True if frequency=0 and at gamma.
    
    ! The unit vector along the normal mode, in mass-weighted co-ordinates.
    ! The unit vectors are stored as displacements from
    !    equilibrium for each atom in the primitive cell.
    ! Normal modes are orthonormal in mass-weighted co-ordinates.
    type(ComplexVector), allocatable :: unit_vector(:)
    
    ! The ID of the q-point at which this mode exists,
    !    and that at which the pair to this mode exists.
    integer :: qpoint_id
    integer :: paired_qpoint_id
    
    ! The ID of the subspace modes which are degenerate with this mode.
    integer :: subspace_id
  contains
    ! Calculates the stress prefactor associated with the mode.
    procedure, public :: stress_prefactor
    ! I/O.
    procedure, public :: read  => read_ComplexMode
    procedure, public :: write => write_ComplexMode
  end type
  
  interface ComplexMode
    ! ----------------------------------------------------------------------
    ! Basic constructor.
    ! ----------------------------------------------------------------------
    module function new_ComplexMode(id,paired_id,frequency,spring_constant, &
       & soft_mode,translational_mode,unit_vector,qpoint_id,                &
       & paired_qpoint_id,subspace_id) result(this) 
      integer,             intent(in) :: id
      integer,             intent(in) :: paired_id
      real(dp),            intent(in) :: frequency
      real(dp),            intent(in) :: spring_constant
      logical,             intent(in) :: soft_mode
      logical,             intent(in) :: translational_mode
      type(ComplexVector), intent(in) :: unit_vector(:)
      integer,             intent(in) :: qpoint_id
      integer,             intent(in) :: paired_qpoint_id
      integer,             intent(in) :: subspace_id
      type(ComplexMode)               :: this
    end function
  
    ! ----------------------------------------------------------------------
    ! Construct a complex mode, but without ids.
    ! ----------------------------------------------------------------------
    module function new_ComplexMode_unprocessed(frequency,unit_vector) &
       & result(this) 
      real(dp),            intent(in) :: frequency
      type(ComplexVector), intent(in) :: unit_vector(:)
      type(ComplexMode)               :: this
    end function
  
    ! ----------------------------------------------------------------------
    ! Construct a complex mode from the eigenstuff of the dynamical matrix.
    ! ----------------------------------------------------------------------
    impure elemental module function new_ComplexMode_HermitianEigenstuff(eigenstuff,structure) result(this) 
      type(HermitianEigenstuff), intent(in) :: eigenstuff
      type(StructureData),       intent(in) :: structure
      type(ComplexMode)                     :: this
    end function
  end interface
  
  interface conjg
    ! ----------------------------------------------------------------------
    ! Take the conjugate of the mode.
    ! This is equivalent to the transformation q -> -q.
    ! ----------------------------------------------------------------------
    impure elemental module function conjg_ComplexMode(input) result(output) 
      type(ComplexMode), intent(in) :: input
      type(ComplexMode)             :: output
    end function
  end interface
  
  interface transform
    ! ----------------------------------------------------------------------
    ! Transform a mode by a symmetry operator.
    ! ----------------------------------------------------------------------
    ! N.B. the ID and paired ID will only be correct if 'symmetry' is the
    !    symmetry used to construct the transformed modes from the input modes.
    impure elemental module function transform_ComplexMode(input,symmetry, &
       & qpoint_from,qpoint_to) result(output) 
      type(ComplexMode),      intent(in) :: input
      type(SymmetryOperator), intent(in) :: symmetry
      type(QpointData),       intent(in) :: qpoint_from
      type(QpointData),       intent(in) :: qpoint_to
      type(ComplexMode)                  :: output
    end function
  end interface
  
  interface select_qpoint
    ! ----------------------------------------------------------------------
    ! Select q-points corresponding to a given mode or modes.
    ! ----------------------------------------------------------------------
    module function select_qpoint_ComplexMode(mode,qpoints) result(output) 
      type(ComplexMode), intent(in) :: mode
      type(QpointData),  intent(in) :: qpoints(:)
      type(QpointData)              :: output
    end function
  end interface
  
  interface select_qpoints
    module function select_qpoints_ComplexModes(modes,qpoints) result(output) 
      type(ComplexMode), intent(in) :: modes(:)
      type(QpointData),  intent(in) :: qpoints(:)
      type(QpointData), allocatable :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Generate the three purely translational modes at the gamma point.
    ! ----------------------------------------------------------------------
    module function generate_translational_modes(structure,qpoints) &
       & result(output) 
      type(StructureData), intent(in) :: structure
      type(QpointData),    intent(in) :: qpoints(:)
      type(ComplexMode), allocatable  :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Generate the stress prefactor associated with the mode,
    !    or with two modes if a second mode is given.
    ! ----------------------------------------------------------------------
    impure elemental module function stress_prefactor(this,that) &
       & result(output) 
      class(ComplexMode), intent(in)           :: this
      class(ComplexMode), intent(in), optional :: that
      type(RealMatrix)                         :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_ComplexMode(this,input) 
      class(ComplexMode), intent(out) :: this
      type(String),       intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_ComplexMode(this) result(output) 
      class(ComplexMode), intent(in) :: this
      type(String), allocatable      :: output(:)
    end function
  end interface
  
  interface ComplexMode
    module function new_ComplexMode_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(ComplexMode)        :: this
    end function
  
    impure elemental module function new_ComplexMode_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(ComplexMode)             :: this
    end function
  end interface
end module
