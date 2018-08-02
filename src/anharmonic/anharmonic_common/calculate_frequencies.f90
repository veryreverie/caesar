! ======================================================================
! Calculates the frequencies of harmonic ground-states along each mode
!    which minimise the total energy of that ground-state.
! ======================================================================
module calculate_frequencies_module
  use common_module
  
  use states_module
  
  use potential_module
  use potential_pointer_module
  use anharmonic_data_module
  implicit none
  
  private
  
  public :: calculate_frequencies
contains

! The output is give as an array indexed as
!    anharmonic_data%degenerate_subspaces.
! (Every mode in a given subspace has the same frequency.)
subroutine calculate_frequencies(potential,anharmonic_data)
  implicit none
  
  class(PotentialData), intent(in) :: potential
  type(AnharmonicData), intent(in) :: anharmonic_data
  
  type(DegenerateSubspace), allocatable :: subspaces(:)
  real(dp),                 allocatable :: frequencies(:)
  
  type(SubspaceState), allocatable :: subspace_states(:)
  type(SubspaceState), allocatable :: states(:)
  
  type(PotentialPointer) :: new_potential
  
  integer :: i,j,ialloc
  
  subspaces = anharmonic_data%degenerate_subspaces
  
  ! Generate the first guess at the frequencies:
  !    - harmonic frequency             if > frequency_of_max_displacement.
  !    - frequency_of_max_displacement  otherwise.
  allocate( frequencies(size(subspaces)), stat=ialloc); call err(ialloc)
  do i=1,size(subspaces)
    frequencies(i) = max( subspaces(i)%frequency,                       &
                        & anharmonic_data%frequency_of_max_displacement )
  enddo
  
  ! Generate ground states at first guess frequencies.
  allocate(states(size(subspaces)), stat=ialloc); call err(ialloc)
  do i=1,size(subspaces)
    subspace_states = generate_subspace_states( &
               & subspaces(i),                  &
               & frequencies(i),                &
               & anharmonic_data%complex_modes, &
               & 0                              )
    if (size(subspace_states)/=1) then
      call err()
    endif
    states(i) = subspace_states(1)
  enddo
  
  do
    do i=1,size(subspaces)
      ! Integrate the potential across all other subspaces.
      new_potential = potential
      do j=1,size(subspaces)
        if (j/=i) then
          call new_potential%braket(states(j),states(j),anharmonic_data)
        endif
      enddo
      
      ! TODO
      ! Find the next guess at the frequency, (using Newton or bisection).
    enddo
    
    ! TODO
    ! Use a Pulay scheme to update frequencies self-consistently.
    exit
  enddo
  
  ! TODO
  ! Generate basis states from frequencies.
  ! Return basis.
  
end subroutine
end module
