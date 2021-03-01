! ======================================================================
! As caesar_vscf_rvector_module, but produces a list of R-vector combinations
!    rather than the R-vectors for each subspace individually.
! ======================================================================
module caesar_vscf_rvectors_module
  use caesar_common_module
  
  use caesar_sampling_points_module
  use caesar_vscf_rvector_module
  implicit none
  
  private
  
  public :: VscfRvectors
  public :: size
  public :: construct_vscf_rvectors
  
  type, extends(Stringsable) :: VscfRvectors
    type(VscfRvector), allocatable :: vscf_rvectors(:)
  contains
    ! List R-vectors for given modes.
    procedure, public :: rvectors => rvectors_VscfRvectors
    
    ! Transform vectors in normal-mode co-ordinates.
    generic,   public  :: transform =>                       &
                        & transform_ComplexModeDisplacement, &
                        & transform_ComplexModeForce,        &
                        & transform_RealModeDisplacement,    &
                        & transform_RealModeForce
    procedure, private :: transform_ComplexModeDisplacement
    procedure, private :: transform_ComplexModeForce
    procedure, private :: transform_RealModeDisplacement
    procedure, private :: transform_RealModeForce
    
    ! The inverse operation to transform.
    generic,   public  :: inverse_transform =>                       &
                        & inverse_transform_ComplexModeDisplacement, &
                        & inverse_transform_ComplexModeForce,        &
                        & inverse_transform_RealModeDisplacement,    &
                        & inverse_transform_RealModeForce
    procedure, private :: inverse_transform_ComplexModeDisplacement
    procedure, private :: inverse_transform_ComplexModeForce
    procedure, private :: inverse_transform_RealModeDisplacement
    procedure, private :: inverse_transform_RealModeForce
    
    ! Helper functions for transform and inverse_transform.
    procedure, private :: transform_complex_magnitudes
    procedure, private :: transform_real_magnitudes
    
    ! I/O.
    procedure, public :: read  => read_VscfRvectors
    procedure, public :: write => write_VscfRvectors
  end type
  
  ! Helper type.
  type :: RvectorArray
    integer                      :: subspace_id
    type(IntVector), allocatable :: rvectors(:)
  end type
  
  interface VscfRvectors
    ! ----------------------------------------------------------------------
    ! Constructor and size() module function for VscfRvectors and RvectorArray.
    ! ----------------------------------------------------------------------
    module function new_VscfRvectors(vscf_rvectors) result(this) 
      type(VscfRvector), intent(in), optional :: vscf_rvectors(:)
      type(VscfRvectors)                      :: this
    end function
  end interface
  
  interface size
    module function size_VscfRvectors(this) result(output) 
      type(VscfRvectors), intent(in) :: this
      integer                        :: output
    end function
  end interface
  
  interface RvectorArray
    module function new_RvectorArray(subspace_id,rvectors) result(this) 
      integer,         intent(in) :: subspace_id
      type(IntVector), intent(in) :: rvectors(:)
      type(RvectorArray)          :: this
    end function
  end interface
  
  interface size
    module function size_RvectorArray(this) result(output) 
      type(RvectorArray), intent(in) :: this
      integer                        :: output
    end function
  end interface
  
  interface operator(//)
    ! ----------------------------------------------------------------------
    ! Concatenation of a VscfRvector to a VscfRvectors.
    ! ----------------------------------------------------------------------
    module function concatenate_VscfRvectors_VscfRvector(this,that) &
       & result(output) 
      type(VscfRvectors), intent(in) :: this
      type(VscfRvector),  intent(in) :: that
      type(VscfRvectors)             :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! List R-vectors for the given set of modes.
    ! ----------------------------------------------------------------------
    module function rvectors_VscfRvectors(this,real_modes) result(output) 
      class(VscfRvectors), intent(in) :: this
      type(RealMode),      intent(in) :: real_modes(:)
      type(IntVector), allocatable    :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Transform a displacement.
    ! ----------------------------------------------------------------------
    module function transform_ComplexModeDisplacement(this,displacement, &
       & modes,qpoints) result(output) 
      class(VscfRvectors),           intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      type(ComplexMode),             intent(in) :: modes(:)
      type(QpointData),              intent(in) :: qpoints(:)
      type(ComplexModeDisplacement)             :: output
    end function
  end interface
  
  interface
    module function transform_ComplexModeForce(this,force,modes,qpoints) &
       & result(output) 
      class(VscfRvectors),    intent(in) :: this
      type(ComplexModeForce), intent(in) :: force
      type(ComplexMode),      intent(in) :: modes(:)
      type(QpointData),       intent(in) :: qpoints(:)
      type(ComplexModeForce)             :: output
    end function
  end interface
  
  interface
    module function transform_RealModeDisplacement(this,displacement,modes, &
       & qpoints) result(output) 
      class(VscfRvectors),        intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      type(RealMode),             intent(in) :: modes(:)
      type(QpointData),           intent(in) :: qpoints(:)
      type(RealModeDisplacement)             :: output
    end function
  end interface
  
  interface
    module function transform_RealModeForce(this,force,modes,qpoints) &
       & result(output) 
      class(VscfRvectors), intent(in) :: this
      type(RealModeForce), intent(in) :: force
      type(RealMode),      intent(in) :: modes(:)
      type(QpointData),    intent(in) :: qpoints(:)
      type(RealModeForce)             :: output
    end function
  end interface
  
  interface
    module function inverse_transform_ComplexModeDisplacement(this, &
       & displacement,modes,qpoints) result(output) 
      class(VscfRvectors),           intent(in) :: this
      type(ComplexModeDisplacement), intent(in) :: displacement
      type(ComplexMode),             intent(in) :: modes(:)
      type(QpointData),              intent(in) :: qpoints(:)
      type(ComplexModeDisplacement)             :: output
    end function
  end interface
  
  interface
    module function inverse_transform_ComplexModeForce(this,force,modes, &
       & qpoints) result(output) 
      class(VscfRvectors),    intent(in) :: this
      type(ComplexModeForce), intent(in) :: force
      type(ComplexMode),      intent(in) :: modes(:)
      type(QpointData),       intent(in) :: qpoints(:)
      type(ComplexModeForce)             :: output
    end function
  end interface
  
  interface
    module function inverse_transform_RealModeDisplacement(this, &
       & displacement,modes,qpoints) result(output) 
      class(VscfRvectors),        intent(in) :: this
      type(RealModeDisplacement), intent(in) :: displacement
      type(RealMode),             intent(in) :: modes(:)
      type(QpointData),           intent(in) :: qpoints(:)
      type(RealModeDisplacement)             :: output
    end function
  end interface
  
  interface
    module function inverse_transform_RealModeForce(this,force,modes,qpoints) &
       & result(output) 
      class(VscfRvectors), intent(in) :: this
      type(RealModeForce), intent(in) :: force
      type(RealMode),      intent(in) :: modes(:)
      type(QpointData),    intent(in) :: qpoints(:)
      type(RealModeForce)             :: output
    end function
  end interface
  
  interface
    module function transform_complex_magnitudes(this,modes,magnitudes, &
       & qpoints,inverse) result(output) 
      class(VscfRvectors), intent(in)  :: this
      type(ComplexMode),   intent(in)  :: modes(:)
      complex(dp),         intent(in)  :: magnitudes(:)
      type(QpointData),    intent(in)  :: qpoints(:)
      logical,             intent(in)  :: inverse
      complex(dp), allocatable         :: output(:)
    end function
  end interface
  
  interface
    module function transform_real_magnitudes(this,modes,magnitudes,qpoints, &
       & inverse) result(output) 
      class(VscfRvectors),   intent(in) :: this
      type(RealMode),        intent(in) :: modes(:)
      real(dp),              intent(in) :: magnitudes(:)
      type(QpointData),      intent(in) :: qpoints(:)
      logical,               intent(in) :: inverse
      real(dp), allocatable             :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Construct VSCF R-vectors.
    ! ----------------------------------------------------------------------
    ! Constructs the list of all R-vectors for each subspace,
    !    and then calls a recursive helper module function to generate all permutations.
    module function construct_vscf_rvectors(sampling_point,supercell, &
       & real_modes,qpoints) result(output) 
      type(RealModeDisplacement), intent(in) :: sampling_point
      type(StructureData),        intent(in) :: supercell
      type(RealMode),             intent(in) :: real_modes(:)
      type(QpointData),           intent(in) :: qpoints(:)
      type(VscfRvectors), allocatable        :: output(:)
    end function
  end interface
  
  interface
    ! First helper module function for construct_vscf_rvectors.
    ! Constructs the R-vectors at which each subspace needs to be sampled.
    module function construct_rvector_arrays(sampling_point,supercell, &
       & real_modes,qpoints) result(output) 
      type(RealModeDisplacement), intent(in) :: sampling_point
      type(StructureData),        intent(in) :: supercell
      type(RealMode),             intent(in) :: real_modes(:)
      type(QpointData),           intent(in) :: qpoints(:)
      type(RvectorArray), allocatable        :: output(:)
    end function
  end interface
  
  interface
    ! Second helper module function for construct_vscf_rvectors.
    ! Recursively constructs the permutations of R-vectors across subspaces.
    ! N.B. rvectors_in should only be included when this module function calls itself. &
    recursive module function list_rvector_permutations(rvector_arrays, &
       & vscf_rvectors_in) result(output) 
      type(RvectorArray), intent(in)           :: rvector_arrays(:)
      type(VscfRvectors), intent(in), optional :: vscf_rvectors_in
      type(VscfRvectors), allocatable          :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_VscfRvectors(this,input) 
      class(VscfRvectors), intent(out) :: this
      type(String),        intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_VscfRvectors(this) result(output) 
      class(VscfRvectors), intent(in) :: this
      type(String), allocatable       :: output(:)
    end function
  end interface
  
  interface VscfRvectors
    module function new_VscfRvectors_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(VscfRvectors)       :: this
    end function
  
    impure elemental module function new_VscfRvectors_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(VscfRvectors)            :: this
    end function
  end interface
end module
