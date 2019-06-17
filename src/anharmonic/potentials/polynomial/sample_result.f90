! ======================================================================
! The energy and force from a sampling_points calculation.
! ======================================================================
! The energy and force are averaged over several VSCF R-vectors.
! The energy is normalised to be per primitive cell,
!    and the force is transformed into real mode co-ordinates.
module sample_result_module
  use common_module
  
  use vscf_rvectors_module
  implicit none
  
  private
  
  public :: SampleResult
  
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
   & qpoints) result(this)
  implicit none
  
  type(ElectronicStructure), intent(in) :: calculation
  type(StructureData),       intent(in) :: supercell
  type(RealMode),            intent(in) :: real_modes(:)
  type(QpointData),          intent(in) :: qpoints(:)
  type(SampleResult)                    :: this
  
  ! Output variables.
  real(dp)            :: energy
  type(RealModeForce) :: force
  type(RealMatrix)    :: stress
  
  ! Normalise the energy by the number of unit cells in the supercell.
  energy = calculation%energy() / supercell%sc_size
  
  ! Transform the forces into normal mode co-ordinates.
  force = RealModeForce( calculation%forces(), &
                       & supercell,            &
                       & real_modes,           &
                       & qpoints               )
  
  ! Construct output.
  if (calculation%has_stress()) then
    stress = calculation%stress()
    this = SampleResult(energy,force,stress)
  else
    this = SampleResult(energy,force)
  endif
end function

! Construct a SampleResult from a set of VSCF R-vectors and the electronic
!    structure at each.
function new_SampleResult_calculations(vscf_rvectors,calculations,supercell, &
   & real_modes,qpoints) result(this)
  implicit none
  
  type(VscfRvectors),        intent(in) :: vscf_rvectors(:)
  type(ElectronicStructure), intent(in) :: calculations(:)
  type(StructureData),       intent(in) :: supercell
  type(RealMode),            intent(in) :: real_modes(:)
  type(QpointData),          intent(in) :: qpoints(:)
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
                             & qpoints          )
    
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
end module
