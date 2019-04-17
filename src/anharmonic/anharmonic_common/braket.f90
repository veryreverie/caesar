! ======================================================================
! Various integrals, of the form <i|i>, <i|j>, <i|X|i> and <i|X|j>,
!    where X is one of:
!       - the kinetic energy.
!       - the harmonic potential.
!       - a general potential.
!       - the kinetic stress (arising from the change in kinetic energy
!                             with strain).
!       - a general stress.
! ======================================================================
module braket_module
  use common_module
  
  use anharmonic_data_module
  use abstract_classes_module
  implicit none
  
  private
  
  public :: generate_subspace_potentials
  public :: braket
  public :: kinetic_energy
  public :: harmonic_potential_energy
  public :: potential_energy
  
  public :: generate_subspace_stresses
  public :: kinetic_stress
  public :: potential_stress
  
  interface braket
    module procedure braket_state
    module procedure braket_state_state
    module procedure braket_state_ComplexUnivariate
    module procedure braket_state_ComplexUnivariate_state
    module procedure braket_state_ComplexMonomial
    module procedure braket_state_ComplexMonomial_state
    module procedure braket_state_ComplexPolynomial
    module procedure braket_state_ComplexPolynomial_state
    module procedure braket_state_potential
    module procedure braket_state_potential_state
    module procedure braket_state_stress
    module procedure braket_state_stress_state
  end interface
  
  interface kinetic_energy
    module procedure kinetic_energy_state
    module procedure kinetic_energy_state_state
  end interface
  
  interface harmonic_potential_energy
    module procedure harmonic_potential_energy_state
    module procedure harmonic_potential_energy_state_state
  end interface
  
  interface potential_energy
    module procedure potential_energy_state
    module procedure potential_energy_state_state
  end interface
  
  interface kinetic_stress
    module procedure kinetic_stress_state
    module procedure kinetic_stress_state_state
  end interface
  
  interface potential_stress
    module procedure potential_stress_state
    module procedure potential_stress_state_state
  end interface
contains

! ----------------------------------------------------------------------
! Takes a potential V and an array of subspace states {|i>}, and generates the
!    set of single-subspace potentials {V_i}, defined by
!    V_i = (prod_{j/=i}<j|)V(prod_{j/=i}|j>).
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
   & subspace_states,anharmonic_data) result(output)
  implicit none
  
  class(PotentialData),     intent(in) :: potential
  type(DegenerateSubspace), intent(in) :: subspaces(:)
  class(SubspaceBasis),     intent(in) :: subspace_bases(:)
  class(SubspaceStates),    intent(in) :: subspace_states(:)
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(PotentialPointer), allocatable  :: output(:)
  
  ! The minimum and maximum indices in each interval.
  integer, allocatable :: mins_in(:)
  integer, allocatable :: maxs_in(:)
  integer, allocatable :: mins_out(:)
  integer, allocatable :: maxs_out(:)
  
  ! The half-way point of each interval.
  integer :: s
  
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
          output(mins_in(i)) = PotentialPointer(                   &
             & subspace_states(j)%integrate( output(mins_in(i)),   &
             &                               subspaces(j),         &
             &                               subspace_bases(j),    &
             &                               anharmonic_data     ) )
        enddo
        do j=mins_in(i),s-1
          output(s) = PotentialPointer(                           &
             & subspace_states(j)%integrate( output(s),           &
             &                               subspaces(j),        &
             &                               subspace_bases(j),   &
             &                               anharmonic_data    ) )
        enddo
      endif
    enddo
    
    mins_in = mins_out
    maxs_in = maxs_out
  enddo
  
  ! Set the constant term to zero.
  do i=1,size(output)
    call output(i)%zero_energy()
  enddo
end function

! ----------------------------------------------------------------------
! <i|i> or <i|j>.
! ----------------------------------------------------------------------
! Calculates <state|state>.
recursive function braket_state(state,subspace,subspace_basis, &
   & anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: state
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp)                             :: output
  
  output = state%braket( subspace        = subspace,       &
                       & subspace_basis  = subspace_basis, &
                       & anharmonic_data = anharmonic_data )
end function

! Calculates <bra|ket>.
recursive function braket_state_state(bra,ket,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: bra
  class(SubspaceState),     intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp)                             :: output
  
  output = bra%braket( ket             = ket,            &
                     & subspace        = subspace,       &
                     & subspace_basis  = subspace_basis, &
                     & anharmonic_data = anharmonic_data )
end function

! ----------------------------------------------------------------------
! <i|X|i> or <i|X|j>, where X is a potential energy type.
! N.B. since |i> may only span a subset of dimensions, these functions return
!    a partially integrated object. To get the total (scalar) integration,
!    the function potential_energy should be used.
! ----------------------------------------------------------------------
! Calculates <state|V|state> for a complex univariate.
recursive function braket_state_ComplexUnivariate(state,potential, &
   & subspace,subspace_basis,anharmonic_data,qpoint) result(output)
  implicit none
  
  class(SubspaceState),     intent(in)           :: state
  type(ComplexUnivariate),  intent(in)           :: potential
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(QpointData),         intent(in), optional :: qpoint
  type(ComplexMonomial)                          :: output
  
  output = state%braket( univariate      = potential,       &
                       & subspace        = subspace,        &
                       & subspace_basis  = subspace_basis,  &
                       & anharmonic_data = anharmonic_data, &
                       & qpoint          = qpoint           )
end function

! Calculates <bra|V|ket> for a complex univariate.
recursive function braket_state_ComplexUnivariate_state(bra,potential,ket, &
   & subspace,subspace_basis,anharmonic_data,qpoint) result(output)
  implicit none
  
  class(SubspaceState),     intent(in)           :: bra
  type(ComplexUnivariate),  intent(in)           :: potential
  class(SubspaceState),     intent(in)           :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(QpointData),         intent(in), optional :: qpoint
  type(ComplexMonomial)                          :: output
  
  output = bra%braket( univariate      = potential,       &
                     & ket             = ket,             &
                     & subspace        = subspace,        &
                     & subspace_basis  = subspace_basis,  &
                     & anharmonic_data = anharmonic_data, &
                     & qpoint          = qpoint           )
end function

! Calculates <state|V|state> for a complex monomial.
recursive function braket_state_ComplexMonomial(state,potential, &
   & subspace,subspace_basis,anharmonic_data,qpoint) result(output)
  implicit none
  
  class(SubspaceState),     intent(in)           :: state
  type(ComplexMonomial),    intent(in)           :: potential
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(QpointData),         intent(in), optional :: qpoint
  type(ComplexMonomial)                          :: output
  
  output = state%braket( monomial        = potential,       &
                       & subspace        = subspace,        &
                       & subspace_basis  = subspace_basis,  &
                       & anharmonic_data = anharmonic_data, &
                       & qpoint          = qpoint           )
end function

! Calculates <bra|V|ket> for a complex monomial.
recursive function braket_state_ComplexMonomial_state(bra,potential,ket, &
   & subspace,subspace_basis,anharmonic_data,qpoint) result(output)
  implicit none
  
  class(SubspaceState),     intent(in)           :: bra
  type(ComplexMonomial),    intent(in)           :: potential
  class(SubspaceState),     intent(in)           :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(QpointData),         intent(in), optional :: qpoint
  type(ComplexMonomial)                          :: output
  
  output = bra%braket( monomial        = potential,       &
                     & ket             = ket,             &
                     & subspace        = subspace,        &
                     & subspace_basis  = subspace_basis,  &
                     & anharmonic_data = anharmonic_data, &
                     & qpoint          = qpoint           )
end function

! Calculates <state|V|state> for a complex polynomial.
recursive function braket_state_ComplexPolynomial(state,potential, &
   & subspace,subspace_basis,anharmonic_data,qpoint) result(output)
  implicit none
  
  class(SubspaceState),     intent(in)           :: state
  type(ComplexPolynomial),  intent(in)           :: potential
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(QpointData),         intent(in), optional :: qpoint
  type(ComplexPolynomial)                        :: output
  
  output = state%braket( polynomial      = potential,       &
                       & subspace        = subspace,        &
                       & subspace_basis  = subspace_basis,  &
                       & anharmonic_data = anharmonic_data, &
                       & qpoint          = qpoint           )
end function

! Calculates <bra|V|ket> for a complex polynomial.
recursive function braket_state_ComplexPolynomial_state(bra,potential,ket, &
   & subspace,subspace_basis,anharmonic_data,qpoint) result(output)
  implicit none
  
  class(SubspaceState),     intent(in)           :: bra
  type(ComplexPolynomial),  intent(in)           :: potential
  class(SubspaceState),     intent(in)           :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(QpointData),         intent(in), optional :: qpoint
  type(ComplexPolynomial)                        :: output
  
  output = bra%braket( polynomial      = potential,       &
                     & ket             = ket,             &
                     & subspace        = subspace,        &
                     & subspace_basis  = subspace_basis,  &
                     & anharmonic_data = anharmonic_data, &
                     & qpoint          = qpoint           )
end function

! Calculates <state|V|state> for a general potential.
recursive function braket_state_potential(state,potential,subspace, &
   & subspace_basis,anharmonic_data,qpoint) result(output)
  implicit none
  
  class(SubspaceState),     intent(in)           :: state
  class(PotentialData),     intent(in)           :: potential
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(QpointData),         intent(in), optional :: qpoint
  type(PotentialPointer)                         :: output
  
  output = potential%braket( bra             = state,           &
                           & subspace        = subspace,        &
                           & subspace_basis  = subspace_basis,  &
                           & anharmonic_data = anharmonic_data, &
                           & qpoint          = qpoint           )
end function

! Calculates <bra|V|ket> for a general potential.
recursive function braket_state_potential_state(bra,potential,ket, &
   & subspace,subspace_basis,anharmonic_data,qpoint) result(output)
  implicit none
  
  class(SubspaceState),     intent(in)           :: bra
  class(PotentialData),     intent(in)           :: potential
  class(SubspaceState),     intent(in)           :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(QpointData),         intent(in), optional :: qpoint
  type(PotentialPointer)                         :: output
  
  output = potential%braket( bra             = bra,             &
                           & ket             = ket,             &
                           & subspace        = subspace,        &
                           & subspace_basis  = subspace_basis,  &
                           & anharmonic_data = anharmonic_data, &
                           & qpoint          = qpoint           )
end function

! ----------------------------------------------------------------------
! Calculates <i|X|i> or <i|X|j>, where X is a stress object.
! As with the potential variant, this may only integrate across a subset
!    of dimensions. potential_stress should be called instead if it is
!    desired to integrate across all dimensions.
! ----------------------------------------------------------------------
! Calculates <state|stress|state>.
recursive function braket_state_stress(state,stress,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: state
  class(StressData),        intent(in) :: stress
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(StressPointer)                  :: output
  
  output = stress%braket( bra             = state,          &
                        & subspace        = subspace,       &
                        & subspace_basis  = subspace_basis, &
                        & anharmonic_data = anharmonic_data )
end function

! Calculates <bra|stress|ket>.
recursive function braket_state_stress_state(bra,stress,ket, &
   & subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: bra
  class(StressData),        intent(in) :: stress
  class(SubspaceState),     intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(StressPointer)                  :: output
  
  output = stress%braket( bra             = bra,            &
                        & ket             = ket,            &
                        & subspace        = subspace,       &
                        & subspace_basis  = subspace_basis, &
                        & anharmonic_data = anharmonic_data )
end function

! ----------------------------------------------------------------------
! The bra-ket of the kinetic energy.
! ----------------------------------------------------------------------
! Calculates <state|T|state>.
recursive function kinetic_energy_state(state,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: state
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp)                             :: output
  
  output = state%kinetic_energy( subspace        = subspace,       &
                               & subspace_basis  = subspace_basis, &
                               & anharmonic_data = anharmonic_data )
end function

! Calculates <bra|T|ket>.
recursive function kinetic_energy_state_state(bra,ket,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: bra
  class(SubspaceState),     intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp)                             :: output
  
  output = bra%kinetic_energy( ket             = ket,            &
                             & subspace        = subspace,       &
                             & subspace_basis  = subspace_basis, &
                             & anharmonic_data = anharmonic_data )
end function

! ----------------------------------------------------------------------
! The bra-ket of the harmonic potential energy.
! ----------------------------------------------------------------------
! Calculates <state|Vh|state>, where Vh is the harmonic potential.
recursive function harmonic_potential_energy_state(state,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: state
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp)                             :: output
  
  output = state%harmonic_potential_energy( &
        & subspace        = subspace,       &
        & subspace_basis  = subspace_basis, &
        & anharmonic_data = anharmonic_data )
end function

! Calculates <bra|Vh|ket>, where Vh is the harmonic potential.
recursive function harmonic_potential_energy_state_state(bra,ket, &
   & subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: bra
  class(SubspaceState),     intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp)                             :: output
  
  output = bra%harmonic_potential_energy( ket             = ket,            &
                                        & subspace        = subspace,       &
                                        & subspace_basis  = subspace_basis, &
                                        & anharmonic_data = anharmonic_data )
end function

! ----------------------------------------------------------------------
! The bra-ket of an arbitrary potential.
! Integrates across all dimensions and returns a real scalar.
! ----------------------------------------------------------------------
! Calculates <state|V|state> as a constant.
recursive function potential_energy_state(state,potential,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: state
  class(PotentialData),     intent(in) :: potential
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp)                             :: output
  
  type(PotentialPointer), allocatable :: integrated_potential
  
  integrated_potential = braket( state,          &
                               & potential,      &
                               & subspace,       &
                               & subspace_basis, &
                               & anharmonic_data )
  output = integrated_potential%undisplaced_energy()
end function

! Calculates <bra|V|ket> as a constant.
recursive function potential_energy_state_state(bra,potential,ket,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: bra
  class(PotentialData),     intent(in) :: potential
  class(SubspaceState),     intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp)                             :: output
  
  type(PotentialPointer), allocatable :: integrated_potential
  
  integrated_potential = braket( bra,            &
                               & potential,      &
                               & ket,            &
                               & subspace,       &
                               & subspace_basis, &
                               & anharmonic_data )
  output = integrated_potential%undisplaced_energy()
end function

! ----------------------------------------------------------------------
! Breaks a general stress into single-subspace stresses,
!    in a manner analogous to generate_subspace_potentials.
! ----------------------------------------------------------------------
function generate_subspace_stresses(stress,subspaces,subspace_bases, &
   & subspace_states,anharmonic_data) result(output)
  implicit none
  
  class(StressData),        intent(in) :: stress
  type(DegenerateSubspace), intent(in) :: subspaces(:)
  class(SubspaceBasis),     intent(in) :: subspace_bases(:)
  class(SubspaceStates),    intent(in) :: subspace_states(:)
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(StressPointer), allocatable     :: output(:)
  
  ! The minimum and maximum indices in each interval.
  integer, allocatable :: mins_in(:)
  integer, allocatable :: maxs_in(:)
  integer, allocatable :: mins_out(:)
  integer, allocatable :: maxs_out(:)
  
  ! The half-way point of each interval.
  integer :: s
  
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
          output(mins_in(i)) = StressPointer(                      &
             & subspace_states(j)%integrate( output(mins_in(i)),   &
             &                               subspaces(j),         &
             &                               subspace_bases(j),    &
             &                               anharmonic_data     ) )
        enddo
        do j=mins_in(i),s-1
          output(s) = StressPointer(                              &
             & subspace_states(j)%integrate( output(s),           &
             &                               subspaces(j),        &
             &                               subspace_bases(j),   &
             &                               anharmonic_data    ) )
        enddo
      endif
    enddo
    
    mins_in = mins_out
    maxs_in = maxs_out
  enddo
  
  ! Set the constant term to zero.
  do i=1,size(output)
    call output(i)%zero_stress()
  enddo
end function

! ----------------------------------------------------------------------
! The bra-ket of the kinetic stress, which is the term in the stress
!    arising from the kinetic energy reducing as the volume increases.
! ----------------------------------------------------------------------
! Calculates <state|stress|state> for the kinetic stress.
recursive function kinetic_stress_state(state,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: state
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(RealMatrix)                     :: output
  
  output = state%kinetic_stress( subspace        = subspace,       &
                               & subspace_basis  = subspace_basis, &
                               & anharmonic_data = anharmonic_data )
end function

! Calculates <bra|stress|ket> for the kinetic stress.
recursive function kinetic_stress_state_state(bra,ket,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: bra
  class(SubspaceState),     intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(RealMatrix)                     :: output
  
  output = bra%kinetic_stress( ket             = ket,            &
                             & subspace        = subspace,       &
                             & subspace_basis  = subspace_basis, &
                             & anharmonic_data = anharmonic_data )
end function

! ----------------------------------------------------------------------
! The bra-ket of an arbitrary stress.
! Integrates across all dimensions and returns a constant tensor.
! ----------------------------------------------------------------------
! Calculates <state|stress|state> as a constant, for the potential stress.
recursive function potential_stress_state(state,stress,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: state
  class(StressData),        intent(in) :: stress
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(RealMatrix)                     :: output
  
  type(StressPointer), allocatable :: integrated_stress
  
  integrated_stress = braket( state,          &
                            & stress,         &
                            & subspace,       &
                            & subspace_basis, &
                            & anharmonic_data )
  output = integrated_stress%undisplaced_stress()
end function

! Calculates <bra|stress|ket> as a constant, for the potential stress.
recursive function potential_stress_state_state(bra,stress,ket,subspace, &
   & subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(SubspaceState),     intent(in) :: bra
  class(StressData),        intent(in) :: stress
  class(SubspaceState),     intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  class(SubspaceBasis),     intent(in) :: subspace_basis
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(RealMatrix)                     :: output
  
  type(StressPointer), allocatable :: integrated_stress
  
  integrated_stress = braket( bra,            &
                            & stress,         &
                            & ket,            &
                            & subspace,       &
                            & subspace_basis, &
                            & anharmonic_data )
  output = integrated_stress%undisplaced_stress()
end function
end module
