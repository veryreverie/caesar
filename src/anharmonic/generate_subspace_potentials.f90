! ======================================================================
! Takes a potential or stress, and a list of subspaces,
!    and generates the single-subspace potentials or stresses.
! ======================================================================
module generate_subspace_potentials_module
  use common_module
  
  use anharmonic_common_module
  implicit none
  
  private
  
  public :: generate_subspace_potentials
  public :: generate_subspace_stresses
contains

! ----------------------------------------------------------------------
! Takes a potential V and an array of subspace states {|i>}, and generates the
!    set of single-subspace potentials {V_i}, defined by
!    V_i = (prod_{j/=i}<j|)V(prod_{j/=i}|j>)
!        - (prod_i<i|)V(prod_i|i>) * (n-1)/n
! ----------------------------------------------------------------------
! The naive method of calculating {V_i} for n subspaces takes
!    n(n-1) operations.
! This can be accelerated using a bisection method, outlined below.
! 
! V0(1) = V, the input potential.
! 
! The first iteration splits the states into two intervals,
!    [1,s-1] and [s,n], where s=n/2, and two potentials are calculated:
!    - V1(1) = (<s|<s+1|...<n|)V0(1)(|s>|s+1>...|n>)
!    - V1(s) = (<1|<2|...<s-1|)V0(1)(|1>|2>...|s-1>)
! These intervals are recorded in terms of their min and max values:
!   mins = [1  , s]
!   maxs = [s-1, n]
!
! The next iteration splits each of the intervals into two intervals,
!    copies the potential to both intervals, and integrates the potential
!    corresponding to each interval over the states in the other interval.
! This method takes O(n.log(n)) operations.
function generate_subspace_potentials(potential,subspaces,subspace_bases, &
   & subspace_states,thermal_energy,anharmonic_data) result(output)
  implicit none
  
  class(PotentialData),     intent(in) :: potential
  type(DegenerateSubspace), intent(in) :: subspaces(:)
  class(SubspaceBasis),     intent(in) :: subspace_bases(:)
  class(BasisStates),       intent(in) :: subspace_states(:)
  real(dp),                 intent(in) :: thermal_energy
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(PotentialPointer), allocatable  :: output(:)
  
  ! The minimum and maximum indices in each interval.
  integer, allocatable :: mins_in(:)
  integer, allocatable :: maxs_in(:)
  integer, allocatable :: mins_out(:)
  integer, allocatable :: maxs_out(:)
  
  ! The half-way point of each interval.
  integer :: s
  
  ! The correction energy, -<V>(n-1)/n
  type(PotentialPointer) :: integrated_potential
  real(dp)               :: correction_energy
  
  integer :: i,j,ialloc
  
  if (size(subspace_states)==0) then
    call print_line(ERROR//': No states to integrate over.')
    call err()
  endif
  
  ! Initialise to a single interval spanning all states,
  !    and the un-integrated potential.
  mins_in = [1]
  maxs_in = [size(subspace_states)]
  allocate(output(size(subspace_states)), stat=ialloc); call err(ialloc)
  output(1) = PotentialPointer(potential)
  call output(1)%zero_energy()
  
  ! Loop over iteration until every interval contains exactly one subspace.
  do while (any(mins_in/=maxs_in))
    mins_out = [integer::]
    maxs_out = [integer::]
    
    ! Loop over intervals.
    do i=1,size(mins_in)
      if (mins_in(i)==maxs_in(i)) then
        ! The interval contains only one subspace; nothing needs doing.
        mins_out = [mins_out, mins_in(i)]
        maxs_out = [maxs_out, maxs_in(i)]
      else
        ! The interval contains more than one subspace;
        !    and generate two integrated potentials.
        
        ! Bisect the interval (or split it as evenly as possible).
        s = mins_in(i) + (maxs_in(i)-mins_in(i)+1)/2
        mins_out = [mins_out, mins_in(i), s     ]
        maxs_out = [maxs_out, s-1   , maxs_in(i)]
        
        ! Copy the potential from the first interval to the second.
        output(s) = output(mins_in(i))
        
        ! Integrate the first potential over all states in the second interval,
        !    and the second potential over all states in the first interval.
        do j=s,maxs_in(i)
          call output(mins_in(i))%braket(            &
             & states          = subspace_states(j), &
             & thermal_energy  = thermal_energy,     &
             & subspace        = subspaces(j),       &
             & subspace_basis  = subspace_bases(j),  &
             & anharmonic_data = anharmonic_data     )
        enddo
        do j=mins_in(i),s-1
          call output(s)%braket(                     &
             & states          = subspace_states(j), &
             & thermal_energy  = thermal_energy,     &
             & subspace        = subspaces(j),       &
             & subspace_basis  = subspace_bases(j),  &
             & anharmonic_data = anharmonic_data     )
        enddo
      endif
    enddo
    
    mins_in = mins_out
    maxs_in = maxs_out
  enddo
  
  ! Calculate the fully-integrated potential,
  !    <V> = (prod_i<i|)V(prod_i|i>).
  ! Use this to correct for the overcounting of terms.
  integrated_potential = output(1)
  
  call integrated_potential%braket( states          = subspace_states(1), &
                                  & thermal_energy  = thermal_energy,     &
                                  & subspace        = subspaces(1),       &
                                  & subspace_basis  = subspace_bases(1),  &
                                  & anharmonic_data = anharmonic_data     )
  
  correction_energy = integrated_potential%undisplaced_energy() &
                  & * (1.0_dp-size(subspaces))/size(subspaces)
  
  call output%add_constant(correction_energy)
  
  ! Process the subspace potentials if necessary.
  output = subspace_bases%process_subspace_potential( output,          &
                                                    & subspace_states, &
                                                    & subspaces,       &
                                                    & thermal_energy,  &
                                                    & anharmonic_data  )
  
  call output%finalise_subspace_potential(subspaces, anharmonic_data)
end function

! ----------------------------------------------------------------------
! Breaks a general stress into single-subspace stresses,
!    in a manner analogous to generate_subspace_potentials.
! ----------------------------------------------------------------------
function generate_subspace_stresses(stress,subspaces,subspace_bases, &
   & subspace_states,thermal_energy,anharmonic_data) result(output)
  implicit none
  
  class(StressData),        intent(in) :: stress
  type(DegenerateSubspace), intent(in) :: subspaces(:)
  class(SubspaceBasis),     intent(in) :: subspace_bases(:)
  class(BasisStates),       intent(in) :: subspace_states(:)
  real(dp),                 intent(in) :: thermal_energy
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(StressPointer), allocatable     :: output(:)
  
  ! The minimum and maximum indices in each interval.
  integer, allocatable :: mins_in(:)
  integer, allocatable :: maxs_in(:)
  integer, allocatable :: mins_out(:)
  integer, allocatable :: maxs_out(:)
  
  ! The half-way point of each interval.
  integer :: s
  
  ! The correction stress, -<S>(n-1)/n
  type(StressPointer) :: integrated_stress
  type(RealMatrix)    :: correction_stress
  
  integer :: i,j,ialloc
  
  if (size(subspace_states)==0) then
    call print_line(ERROR//': No states to integrate over.')
    call err()
  endif
  
  ! Initialise to a single interval spanning all states,
  !    and the un-integrated stress.
  mins_in = [1]
  maxs_in = [size(subspace_states)]
  allocate(output(size(subspace_states)), stat=ialloc); call err(ialloc)
  output(1) = StressPointer(stress)
  
  ! Loop over iteration until every interval contains exactly one subspace.
  do while (any(mins_in/=maxs_in))
    mins_out = [integer::]
    maxs_out = [integer::]
    
    ! Loop over intervals.
    do i=1,size(mins_in)
      if (mins_in(i)==maxs_in(i)) then
        ! The interval contains only one subspace; nothing needs doing.
        mins_out = [mins_out, mins_in(i)]
        maxs_out = [maxs_out, maxs_in(i)]
      else
        ! The interval contains more than one subspace;
        !    and generate two integrated stresses.
        
        ! Bisect the interval (or split it as evenly as possible).
        s = mins_in(i) + (maxs_in(i)-mins_in(i)+1)/2
        mins_out = [mins_out, mins_in(i), s     ]
        maxs_out = [maxs_out, s-1   , maxs_in(i)]
        
        ! Copy the stress from the first interval to the second.
        output(s) = output(mins_in(i))
        
        ! Integrate the first stress over all states in the second interval,
        !    and the second stress over all states in the first interval.
        do j=s,maxs_in(i)
          call output(mins_in(i))%braket(            &
             & states          = subspace_states(j), &
             & thermal_energy  = thermal_energy,     &
             & subspace        = subspaces(j),       &
             & subspace_basis  = subspace_bases(j),  &
             & anharmonic_data = anharmonic_data     )
        enddo
        do j=mins_in(i),s-1
          call output(s)%braket(                     &
             & states          = subspace_states(j), &
             & thermal_energy  = thermal_energy,     &
             & subspace        = subspaces(j),       &
             & subspace_basis  = subspace_bases(j),  &
             & anharmonic_data = anharmonic_data     )
        enddo
      endif
    enddo
    
    mins_in = mins_out
    maxs_in = maxs_out
  enddo
  
  ! Calculate the fully-integrated potential,
  !    <V> = (prod_i<i|)V(prod_i|i>).
  ! Use this to correct for the overcounting of terms.
  integrated_stress = output(1)
  
  call integrated_stress%braket( states          = subspace_states(1), &
                               & thermal_energy  = thermal_energy,     &
                               & subspace        = subspaces(1),       &
                               & subspace_basis  = subspace_bases(1),  &
                               & anharmonic_data = anharmonic_data     )
  
  correction_stress = integrated_stress%undisplaced_stress() &
                  & * (1.0_dp-size(subspaces))/size(subspaces)
  
  call output%add_constant(correction_stress)
  
  ! Process the subspace stresses if necessary.
  output = subspace_bases%process_subspace_stress( output,          &
                                                 & subspace_states, &
                                                 & subspaces,       &
                                                 & thermal_energy,  &
                                                 & anharmonic_data  )
end function
end module
