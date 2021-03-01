! ======================================================================
! As caesar_complex_mode_module, but in real co-ordinates.
! ======================================================================
module caesar_real_mode_module
  use caesar_utils_module
  
  use caesar_structure_module
  
  use caesar_cartesian_displacement_module
  use caesar_cartesian_force_module
  use caesar_mass_weighted_displacement_module
  use caesar_mass_weighted_force_module
  implicit none
  
  private
  
  public :: RealMode
  public :: MassWeightedDisplacement
  public :: MassWeightedForce
  public :: CartesianDisplacement
  public :: CartesianForce
  public :: select_qpoint
  public :: select_qpoints
  
  ! A normal mode in real co-ordinates.
  type, extends(Stringsable) :: RealMode
    integer :: id
    integer :: paired_id
    
    ! The frequency, and frequency-relevant information.
    real(dp) :: frequency
    real(dp) :: spring_constant    ! k from v(u) = k*u*u/2.
    logical  :: soft_mode          ! True if frequency < -epsilon.
    logical  :: translational_mode ! True if frequency=0 and at gamma.
    
    ! Unit vectors along the normal mode, in mass-weighted co-ordinates.
    ! The unit vectors are stored as displacements from
    !    equilibrium for each atom in a single primitive cell.
    ! Normal modes are orthonormal in mass-weighted co-ordinates.
    ! To get displacements in a specific unit cell,
    !    the cos_vector is multiplied by cos(2*pi*q.R),
    !    and the sin_vector is multiplied by sin(2*pi*q.R).
    type(RealVector), allocatable :: cos_vector(:)
    type(RealVector), allocatable :: sin_vector(:)
    
    ! The IDs of the q-points at which this mode exists.
    ! N.B. unlike complex modes, real modes are not localised to a single
    !    q-point, but rather a superposition across +/- pair.
    integer :: qpoint_id_plus
    integer :: qpoint_id_minus
    
    ! The ID of the subspace of modes which are degenerate with this mode.
    integer :: subspace_id
  contains
    procedure, private :: construct_vector
    ! I/O.
    procedure, public :: read  => read_RealMode
    procedure, public :: write => write_RealMode
  end type
  
  interface RealMode
    ! ----------------------------------------------------------------------
    ! Basic constructor.
    ! ----------------------------------------------------------------------
    module function new_RealMode(id,paired_id,frequency,spring_constant,    &
       & soft_mode,translational_mode,cos_vector,sin_vector,qpoint_id_plus, &
       & qpoint_id_minus,subspace_id) result(this) 
      integer,          intent(in) :: id
      integer,          intent(in) :: paired_id
      real(dp),         intent(in) :: frequency
      real(dp),         intent(in) :: spring_constant
      logical,          intent(in) :: soft_mode
      logical,          intent(in) :: translational_mode
      type(RealVector), intent(in) :: cos_vector(:)
      type(RealVector), intent(in) :: sin_vector(:)
      integer,          intent(in) :: qpoint_id_plus
      integer,          intent(in) :: qpoint_id_minus
      integer,          intent(in) :: subspace_id
      type(RealMode)               :: this
    end function
  end interface
  
  interface MassWeightedDisplacement
    ! ----------------------------------------------------------------------
    ! Return the mode in the mass-weighted or cartesian co-ordinates of
    !    a given supercell.
    ! ----------------------------------------------------------------------
    impure elemental module function new_MassWeightedDisplacement_RealMode(this,structure,qpoint) result(output) 
      class(RealMode),     intent(in) :: this
      type(StructureData), intent(in) :: structure
      type(QpointData),    intent(in) :: qpoint
      type(MassWeightedDisplacement)  :: output
    end function
  end interface
  
  interface MassWeightedForce
    impure elemental module function new_MassWeightedForce_RealMode(this, &
       & structure,qpoint) result(output) 
      class(RealMode),     intent(in) :: this
      type(StructureData), intent(in) :: structure
      type(QpointData),    intent(in) :: qpoint
      type(MassWeightedForce)         :: output
    end function
  end interface
  
  interface CartesianDisplacement
    impure elemental module function new_CartesianDisplacement_RealMode(this,structure,qpoint) result(output) 
      class(RealMode),     intent(in) :: this
      type(StructureData), intent(in) :: structure
      type(QpointData),    intent(in) :: qpoint
      type(CartesianDisplacement)     :: output
    end function
  end interface
  
  interface CartesianForce
    impure elemental module function new_CartesianForce_RealMode(this, &
       & structure,qpoint) result(output) 
      class(RealMode),     intent(in) :: this
      type(StructureData), intent(in) :: structure
      type(QpointData),    intent(in) :: qpoint
      type(CartesianForce)            :: output
    end function
  end interface
  
  interface
    module function construct_vector(this,structure,qpoint) result(output) 
      class(RealMode),     intent(in) :: this
      type(StructureData), intent(in) :: structure
      type(QpointData),    intent(in) :: qpoint
      type(RealVector), allocatable   :: output(:)
    end function
  end interface
  
  interface select_qpoint
    ! ----------------------------------------------------------------------
    ! Select q-points corresponding to a given mode or modes.
    ! ----------------------------------------------------------------------
    module function select_qpoint_RealMode(mode,qpoints) result(output) 
      type(RealMode),   intent(in) :: mode
      type(QpointData), intent(in) :: qpoints(:)
      type(QpointData)             :: output
    end function
  end interface
  
  interface select_qpoints
    module function select_qpoints_RealModes(modes,qpoints) result(output) 
      type(RealMode),   intent(in)  :: modes(:)
      type(QpointData), intent(in)  :: qpoints(:)
      type(QpointData), allocatable :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_RealMode(this,input) 
      class(RealMode), intent(out) :: this
      type(String),    intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_RealMode(this) result(output) 
      class(RealMode), intent(in) :: this
      type(String), allocatable   :: output(:)
    end function
  end interface
  
  interface RealMode
    module function new_RealMode_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(RealMode)           :: this
    end function
  end interface
  
  interface RealMode
    impure elemental module function new_RealMode_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(RealMode)                :: this
    end function
  end interface
end module
