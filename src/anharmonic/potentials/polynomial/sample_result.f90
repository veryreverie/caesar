! ======================================================================
! The energy and force from a sampling_points calculation.
! ======================================================================
! The energy and force are averaged over several VSCF R-vectors.
! The energy is normalised to be per anharmonic supercell,
!    and the force is transformed into real mode co-ordinates.
! N.B. normal-mode forces are extensive.
module caesar_sample_result_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  
  use caesar_vscf_rvectors_module
  use caesar_sampling_points_module
  implicit none
  
  private
  
  public :: SampleResult
  public :: construct_sample_vector
  
  type, extends(NoDefaultConstructor) :: SampleResult
    real(dp)                               :: energy
    type(RealModeForce)                    :: force
    type(RealMatrix), allocatable, private :: stress_
  contains
    procedure, public :: has_stress => has_stress_SampleResult
    procedure, public :: stress => stress_SampleResult
  end type
  
  interface SampleResult
    ! Constructor.
    module function new_SampleResult(energy,force,stress) result(this) 
      real(dp),            intent(in)           :: energy
      type(RealModeForce), intent(in)           :: force
      type(RealMatrix),    intent(in), optional :: stress
      type(SampleResult)                        :: this
    end function
  end interface
  
  interface
    ! Getters for the stress.
    impure elemental module function has_stress_SampleResult(this) &
       & result(output) 
      class(SampleResult), intent(in) :: this
      logical                         :: output
    end function
  end interface
  
  interface
    impure elemental module function stress_SampleResult(this) result(output) 
      class(SampleResult), intent(in) :: this
      type(RealMatrix)                :: output
    end function
  end interface
  
  interface SampleResult
    ! Construct a SampleResult from an electronic structure.
    module function new_SampleResult_calculation(calculation,supercell, &
       & real_modes,qpoints,anharmonic_data) result(this) 
      type(ElectronicStructure), intent(in) :: calculation
      type(StructureData),       intent(in) :: supercell
      type(RealMode),            intent(in) :: real_modes(:)
      type(QpointData),          intent(in) :: qpoints(:)
      type(AnharmonicData),      intent(in) :: anharmonic_data
      type(SampleResult)                    :: this
    end function
  
    ! Construct a SampleResult from a set of VSCF R-vectors and the electronic
    !    structure at each.
    module function new_SampleResult_calculations(vscf_rvectors,    &
       & calculations,supercell,real_modes,qpoints,anharmonic_data) &
       & result(this) 
      type(VscfRvectors),        intent(in) :: vscf_rvectors(:)
      type(ElectronicStructure), intent(in) :: calculations(:)
      type(StructureData),       intent(in) :: supercell
      type(RealMode),            intent(in) :: real_modes(:)
      type(QpointData),          intent(in) :: qpoints(:)
      type(AnharmonicData),      intent(in) :: anharmonic_data
      type(SampleResult)                    :: this
    end function
  end interface
  
  interface
    ! Construct the vector neccesary for L2 fitting of basis function coefficients.
    ! The components of the vector are the measured energies and forces at a given
    !    sampling point, minus the energies and forces already accounted for
    !    by the potential (if given).
    module function construct_sample_vector(sampling_points,sample_results, &
       & potential,modes,energy_force_ratio,sample_weights) result(output) 
      type(RealModeDisplacement), intent(in)           :: sampling_points(:)
      type(SampleResult),         intent(in)           :: sample_results(:)
      class(PotentialData),       intent(in), optional :: potential
      type(RealMode),             intent(in)           :: modes(:)
      real(dp),                   intent(in)           :: energy_force_ratio
      real(dp),                   intent(in), optional :: sample_weights(:)
      real(dp), allocatable                            :: output(:)
    end function
  end interface
end module
