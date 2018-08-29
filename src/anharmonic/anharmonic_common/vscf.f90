! ======================================================================
! Run VSCF on a potential to generate VSCF ground states.
! ======================================================================
module vscf_module
  use common_module
  
  use states_module
  
  use anharmonic_data_module
  use potential_module
  use potential_pointer_module
  implicit none
  
  private
  
  public :: run_vscf
contains

function run_vscf(potential,basis,energy_convergence,max_pulay_iterations, &
   & pre_pulay_iterations,pre_pulay_damping,anharmonic_data) result(output)
  implicit none
  
  class(PotentialData), intent(in)   :: potential
  type(SubspaceBasis),  intent(in)   :: basis(:)
  real(dp),             intent(in)   :: energy_convergence
  integer,              intent(in)   :: max_pulay_iterations
  integer,              intent(in)   :: pre_pulay_iterations
  real(dp),             intent(in)   :: pre_pulay_damping
  type(AnharmonicData), intent(in)   :: anharmonic_data
  type(VscfGroundState), allocatable :: output(:)
  
  type(VscfGroundState), allocatable :: old_guess(:)
  type(VscfGroundState), allocatable :: new_guess(:)
  
  type(RealVector), allocatable :: old_coefficients(:)
  type(RealVector), allocatable :: new_coefficients(:)
  
  old_guess = initial_ground_state(basis)
  new_guess = update(old_guess,basis,potential,anharmonic_data)
  
  output = new_guess
end function

! For each subspace, integrate the potential across the ground states of all
!    other subspaces to give a single-subspace potential, and then diagonalise
!    this potential to give a new ground state.
function update(states,basis,potential,anharmonic_data) result(output)
  implicit none
  
  type(VscfGroundState), intent(in)  :: states(:)
  type(SubspaceBasis),   intent(in)  :: basis(:)
  class(PotentialData),  intent(in)  :: potential
  type(AnharmonicData),  intent(in)  :: anharmonic_data
  type(VscfGroundState), allocatable :: output(:)
  
  type(PotentialPointer) :: new_potential
  
  integer :: i,j
  
  do i=1,size(states)
    new_potential = potential
    do j=1,size(states)
      if (j/=i) then
        call new_potential%braket( SumState(states(j),basis(j)), &
                                 & SumState(states(j),basis(j)), &
                                 & anharmonic_data               )
      endif
    enddo
    
    ! TODO
  enddo
  
  output = states
end function

! Concatenates the coefficients of all the states together,
!    to allow the Pulay scheme to be called.
function coefficients(states) result(output)
  implicit none
  
  type(VscfGroundState), intent(in) :: states(:)
  real(dp), allocatable             :: output(:)
  
  integer :: i
  
  if (any(.not.is_int(states%wavevector))) then
    call print_line(ERROR//': The ground state has gained a finite &
       &wavevector. Caesar is at present unable to handle this.')
    call err()
  endif
  
  output = [real(dp)::]
  
  do i=1,size(states)
    output = [output, states(i)%coefficients]
  enddo
end function
end module
