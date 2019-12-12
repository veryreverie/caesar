! ======================================================================
! The energy and force from a sampling_points calculation.
! ======================================================================
! The energy and force are averaged over several VSCF R-vectors.
! The energy is normalised to be per anharmonic supercell,
!    and the force is transformed into real mode co-ordinates.
! N.B. normal-mode forces are extensive.
module sample_result_module
  use common_module
  
  use anharmonic_common_module
  
  use vscf_rvectors_module
  use sampling_points_module
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
    module procedure new_SampleResult
    module procedure new_SampleResult_calculation
    module procedure new_SampleResult_calculations
  end interface
contains

! Constructor.
function new_SampleResult(energy,force,stress) result(this)
  implicit none
  
  real(dp),            intent(in)           :: energy
  type(RealModeForce), intent(in)           :: force
  type(RealMatrix),    intent(in), optional :: stress
  type(SampleResult)                        :: this
  
  this%energy = energy
  this%force  = force
  if (present(stress)) then
    this%stress_ = stress
  endif
end function

! Getters for the stress.
impure elemental function has_stress_SampleResult(this) result(output)
  implicit none
  
  class(SampleResult), intent(in) :: this
  logical                         :: output
  
  output = allocated(this%stress_)
end function

impure elemental function stress_SampleResult(this) result(output)
  implicit none
  
  class(SampleResult), intent(in) :: this
  type(RealMatrix)                :: output
  
  if (this%has_stress()) then
    output = this%stress_
  else
    call print_line(ERROR//': Sample result does not contain stress.')
    call err()
  endif
end function

! Construct a SampleResult from an electronic structure.
function new_SampleResult_calculation(calculation,supercell,real_modes, &
   & qpoints,anharmonic_data) result(this)
  implicit none
  
  type(ElectronicStructure), intent(in) :: calculation
  type(StructureData),       intent(in) :: supercell
  type(RealMode),            intent(in) :: real_modes(:)
  type(QpointData),          intent(in) :: qpoints(:)
  type(AnharmonicData),      intent(in) :: anharmonic_data
  type(SampleResult)                    :: this
  
  ! Output variables.
  real(dp)                      :: energy
  type(RealModeForce)           :: force
  type(RealMatrix), allocatable :: stress
  
  ! Normalise the energy to be per anharmonic supercell.
  energy = (calculation%energy() / supercell%sc_size) &
       & * anharmonic_data%anharmonic_supercell%sc_size
  
  ! Transform the forces into normal mode co-ordinates.
  force = RealModeForce( calculation%forces(),  &
      &                 supercell,              &
      &                 real_modes,             &
      &                 qpoints               ) &
      & * real(anharmonic_data%anharmonic_supercell%sc_size,dp)
  
  ! Construct output.
  if (calculation%has_stress()) then
    stress = calculation%stress()
  endif
  
  this = SampleResult(energy,force,stress)
end function

! Construct a SampleResult from a set of VSCF R-vectors and the electronic
!    structure at each.
function new_SampleResult_calculations(vscf_rvectors,calculations,supercell, &
   & real_modes,qpoints,anharmonic_data) result(this)
  implicit none
  
  type(VscfRvectors),        intent(in) :: vscf_rvectors(:)
  type(ElectronicStructure), intent(in) :: calculations(:)
  type(StructureData),       intent(in) :: supercell
  type(RealMode),            intent(in) :: real_modes(:)
  type(QpointData),          intent(in) :: qpoints(:)
  type(AnharmonicData),      intent(in) :: anharmonic_data
  type(SampleResult)                    :: this
  
  type(SampleResult), allocatable :: results(:)
  
  ! Output variables.
  real(dp)            :: energy
  type(RealModeForce) :: force
  type(RealMatrix)    :: stress
  
  ! Temporary variables.
  integer :: i,ialloc
  
  ! Normalise energies and change co-ordinates of forces.
  allocate(results(size(calculations)), stat=ialloc); call err(ialloc)
  do i=1,size(calculations)
    results(i) = SampleResult( calculations(i), &
                             & supercell,       &
                             & real_modes,      &
                             & qpoints,         &
                             & anharmonic_data  )
    
    ! Reverse the VSCF R-vector transformation.
    results(i)%force = vscf_rvectors(i)%inverse_transform( results(i)%force,  &
                                                         & real_modes,        &
                                                         & qpoints            )
  enddo
  
  ! Average over calculations.
  energy = sum(results%energy) / size(results)
  force = sum(results%force) / real(size(results),dp)
  if (all(results%has_stress())) then
    stress = sum(results%stress()) / size(results)
    this = SampleResult(energy, force, stress)
  else
    this = SampleResult(energy, force)
  endif
end function

! Construct the vector neccesary for L2 fitting of basis function coefficients.
! The components of the vector are the measured energies and forces at a given
!    sampling point, minus the energies and forces already accounted for
!    by the potential (if given).
function construct_sample_vector(sampling_points,sample_results,potential, &
   & modes,energy_force_ratio) result(output)
  implicit none
  
  type(RealModeDisplacement), intent(in)           :: sampling_points(:)
  type(SampleResult),         intent(in)           :: sample_results(:)
  class(PotentialData),       intent(in), optional :: potential
  type(RealMode),             intent(in)           :: modes(:)
  real(dp),                   intent(in)           :: energy_force_ratio
  real(dp), allocatable                            :: output(:)
  
  real(dp)            :: energy
  type(RealModeForce) :: forces
  
  integer :: dims
  
  integer :: i,ialloc
  
  dims = 1+size(modes)
  
  allocate( output(size(sampling_points)*dims), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(sampling_points)
    energy = sample_results(i)%energy
    forces = sample_results(i)%force
    
    ! Subtract the energy and forces from the existing potential.
    if (present(potential)) then
      energy = energy - potential%energy(sampling_points(i))
      forces = forces - potential%force(sampling_points(i))
    endif
    
    ! Construct the vector.
    output((i-1)*dims+1:i+dims) = make_sample_vector( energy,            &
                                                    & forces,            &
                                                    & modes,             &
                                                    & energy_force_ratio )
  enddo
end function
end module
