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
    real(dp)            :: energy
    type(RealModeForce) :: force
  end type
  
  interface SampleResult
    module procedure new_SampleResult
    module procedure new_SampleResult_calculations
  end interface
contains

! Constructor.
function new_SampleResult(energy,force) result(this)
  implicit none
  
  real(dp),            intent(in) :: energy
  type(RealModeForce), intent(in) :: force
  type(SampleResult)              :: this
  
  this%energy = energy
  this%force  = force
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
  
  ! Output variables.
  real(dp)            :: energy
  type(RealModeForce) :: force
  
  ! Working variables for forces.
  type(RealModeForce), allocatable :: forces(:)
  
  ! Temporary variables.
  integer :: i,ialloc
  
  ! Average the energy over the calculations, and normalise to be per
  !    primitive cell.
  energy = sum(calculations%energy) / (size(calculations) * supercell%sc_size)
  
  ! Construct forces in the correct representation.
  
  allocate(forces(size(calculations)), stat=ialloc); call err(ialloc)
  do i=1,size(forces)
    ! Transform forces into real mode co-ordinates.
    forces(i) = RealModeForce( calculations(i)%forces, &
                             & supercell,              &
                             & real_modes,             &
                             & qpoints)
    ! Reverse the VSCF R-vector transformation.
    forces(i) = vscf_rvectors(i)%inverse_transform( forces(i),  &
                                                  & real_modes, &
                                                  & qpoints)
  enddo
  ! Average the force over the calculations.
  force = sum(forces) / real(size(calculations),dp)
  
  this = SampleResult(energy,force)
end function
end module
