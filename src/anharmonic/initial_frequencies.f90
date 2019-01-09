! ======================================================================
! Finds the frequencies of the self-consistent harmonic potential which
!    minimises the free energy at zero temperature.
! ======================================================================
module initial_frequencies_module
  use common_module
  
  use states_module
  
  use anharmonic_common_module
  use potentials_module
  implicit none
  
  private
  
  public :: InitialFrequencies
  
  type, extends(Stringsable) :: InitialFrequencies
    integer,  allocatable, private :: subspace_ids_(:)
    real(dp), allocatable, private :: frequencies_(:)
  contains
    procedure, public :: frequency => frequency_InitialFrequencies
    
    ! I/O.
    procedure, public :: read  => read_InitialFrequencies
    procedure, public :: write => write_InitialFrequencies
  end type
  
  interface InitialFrequencies
    module procedure new_InitialFrequencies
    module procedure new_InitialFrequencies_PotentialData
    module procedure new_InitialFrequencies_Strings
    module procedure new_InitialFrequencies_StringArray
  end interface
  
  interface size
    module procedure size_InitialFrequencies
  end interface
contains

! Constructor and size function.
function new_InitialFrequencies(subspace_ids,frequencies) result(this)
  implicit none
  
  integer,  intent(in)     :: subspace_ids(:)
  real(dp), intent(in)     :: frequencies(:)
  type(InitialFrequencies) :: this
  
  if (size(subspace_ids)/=size(frequencies)) then
    call print_line(CODE_ERROR//': subspace_ids and frequencies do not match.')
    call err()
  endif
  
  this%subspace_ids_ = subspace_ids
  this%frequencies_ = frequencies
end function

function size_InitialFrequencies(this) result(output)
  implicit none
  
  type(InitialFrequencies), intent(in) :: this
  integer                              :: output
  
  output = size(this%frequencies_)
end function

! Returns the frequency corresponding to the given subspace id.
impure elemental function frequency_InitialFrequencies(this,subspace_id) &
   & result(output)
  implicit none
  
  class(InitialFrequencies), intent(in) :: this
  integer,                   intent(in) :: subspace_id
  real(dp)                              :: output
  
  output = this%frequencies_(first(this%subspace_ids_==subspace_id))
end function

! Calculate initial frequencies from a potential.
function new_InitialFrequencies_PotentialData(potential,anharmonic_data, &
   & frequency_convergence,max_pulay_iterations,pre_pulay_iterations,    &
   & pre_pulay_damping) result(output)
  implicit none
  
  class(PotentialData), intent(in) :: potential
  type(AnharmonicData), intent(in) :: anharmonic_data
  real(dp),             intent(in) :: frequency_convergence
  integer,              intent(in) :: max_pulay_iterations
  integer,              intent(in) :: pre_pulay_iterations
  real(dp),             intent(in) :: pre_pulay_damping
  type(InitialFrequencies)         :: output
  
  type(DegenerateSubspace), allocatable :: subspaces(:)
  real(dp),                 allocatable :: frequencies(:)
  
  type(MonomialState), allocatable :: states(:)
  type(MonomialState), allocatable :: subspace_states(:)
  
  type(RealVector), allocatable :: old_frequencies(:)
  type(RealVector), allocatable :: new_frequencies(:)
  
  type(RealVector) :: next_old_frequencies
  
  integer :: i,j,ialloc
  
  if (pre_pulay_damping<0_dp .or. pre_pulay_damping>1_dp) then
    call print_line(ERROR//': pre_pulay_damping must be between 0 and 1.')
    stop
  elseif (pre_pulay_iterations<0) then
    call print_line(ERROR//': pre_pulay_iterations must be >= 0.')
    stop
  elseif (max_pulay_iterations<0) then
    call print_line(ERROR//': max_pulay_iterations must be >= 0.')
    stop
  endif
  
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
    subspace_states = generate_monomial_states( &
               & subspaces(i),                  &
               & frequencies(i),                &
               & anharmonic_data%complex_modes, &
               & maximum_power = 0              )
    if (size(subspace_states)/=1) then
      call err()
    endif
    states(i) = subspace_states(1)
  enddo
  
  ! Find self-consistent frequencies which minimise the energy.
  old_frequencies = [vec(frequencies)]
  new_frequencies = [RealVector::]
  i = 1
  iter: do
    ! Update the states to have the new frequencies.
    states%frequency = dble(old_frequencies(i))
    
    ! Calculate mean-field potentials from the old frequencies, and use these
    !    mean-field potentials to calculate new frequencies.
    new_frequencies = [ new_frequencies,                              &
                      & optimise_frequencies( potential,              &
                      &                       states,                 &
                      &                       anharmonic_data,        &
                      &                       frequency_convergence ) ]
    
    ! Use a damped iterative scheme or a Pulay scheme to converge towards the
    !    self-consistent solution where new frequencies = old frequencies.
    if (i<=pre_pulay_iterations) then
      next_old_frequencies = (1-pre_pulay_damping) * old_frequencies(i) &
                         & + pre_pulay_damping * new_frequencies(i)
    else
      j = max(2, i-max_pulay_iterations+1)
      next_old_frequencies = pulay(old_frequencies(j:), new_frequencies(j:))
    endif
    
    old_frequencies = [old_frequencies, next_old_frequencies]
    
    ! Increment the loop counter.
    i = i+1
    
    ! Check whether the frequencies have converged.
    ! If the frequencies aren't converged, return to the top of the loop.
    if (i<=5) then
      cycle iter
    else
      do j=1,5
        if ( sum(abs(dble(old_frequencies(i)-old_frequencies(i-j)))) &
         & > frequency_convergence                                   ) then
          cycle iter
        endif
      enddo
    endif
    
    ! If the frequencies are converged, break out of the loop.
    frequencies = dble(old_frequencies(i))
    exit iter
  enddo iter
  
  output = InitialFrequencies(subspaces%id, frequencies)
end function

! For each subspace, integrate the potential across all other subspaces
!    to get a single-subspace mean-field potential.
! Then find the subspace frequency which minimises the single-subspace energy.
function optimise_frequencies(potential,states,anharmonic_data, &
   & frequency_convergence) result(output)
  implicit none
  
  class(PotentialData), intent(in) :: potential
  type(MonomialState),  intent(in) :: states(:)
  type(AnharmonicData), intent(in) :: anharmonic_data
  real(dp),             intent(in) :: frequency_convergence
  type(RealVector)                 :: output
  
  type(PotentialPointer), allocatable :: subspace_potentials(:)
  
  real(dp), allocatable  :: new_frequencies(:)
  
  integer :: i
  
  ! Calculate the single-subspace potentials, defined as
  !    V_i = (prod_{j/=i}<j|)V(prod_{j/=i}|j>).
  subspace_potentials = generate_subspace_potentials( potential,      &
                                                    & states,         &
                                                    & anharmonic_data )
  
  new_frequencies = [(0.0_dp,i=1,size(states))]
  
  ! Calculate update frequencies.
  do i=1,size(states)
    ! Find the frequency which minimises total energy.
    new_frequencies(i) = optimise_frequency(      &
       & subspace_potentials(i),                  &
       & states(i),                               &
       & anharmonic_data,                         &
       & anharmonic_data%degenerate_subspaces(i), &
       & frequency_convergence )
  enddo
  
  output = new_frequencies
end function

! Find the frequency which minimises energy in a single subspace.
function optimise_frequency(potential,state,anharmonic_data,subspace, &
   & frequency_convergence) result(output)
  implicit none
  
  class(PotentialData),     intent(in) :: potential
  type(MonomialState),      intent(in) :: state
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: frequency_convergence
  real(dp)                             :: output
  
  real(dp) :: frequencies(3)
  real(dp) :: energies(3)
  
  type(MonomialState)    :: new_state
  
  real(dp) :: first_derivative
  real(dp) :: second_derivative
  
  real(dp) :: old_frequency
  real(dp) :: new_frequency
  
  integer :: i
  
  old_frequency = state%frequency
  new_state = state
  
  do
    ! Calculate [w-dw, w, w+dw].
    frequencies = [ old_frequency - 0.01_dp*frequency_convergence, &
                  & old_frequency,                                 &
                  & old_frequency + 0.01_dp*frequency_convergence  ]
    
    ! Calculate [U(w-dw), U(w), U(w+dw)].
    do i=1,3
      new_state%frequency = frequencies(i)
      energies(i) = potential%potential_energy( new_state,               &
                &                               new_state,               &
                &                               anharmonic_data )        &
                & + kinetic_energy( new_state,                           &
                &                   new_state,                           &
                &                   subspace,                            &
                &                   anharmonic_data%anharmonic_supercell )
    enddo
    
    ! Calculate dU/dw = (U(w+dw)-U(w-dw))/(2dw).
    first_derivative = (energies(3)-energies(1)) &
                   & / (0.02_dp*frequency_convergence)
    
    ! Calculate d2U/dw2 = (U(w+dw)+U(w-dw)-2U(w))/(dw)^2.
    second_derivative = (energies(1)+energies(3)-2*energies(2)) &
                    & / (0.01_dp*frequency_convergence)**2
    
    ! Update the frequency, and check for convergence.
    ! At the extrema (w=w1), dU/dw=0. As such, w1 = w - (dU/dw)/(d2U/dw2).
    ! If |w1-w|>w/2, or if dU/dw<0 then cap |w1-w| at w/2.
    if (abs(old_frequency)*second_derivative<=abs(first_derivative)) then
      if (first_derivative>0) then
        new_frequency = 0.5_dp * old_frequency
      elseif (first_derivative<0) then
        new_frequency = 1.5_dp * old_frequency
      else
        new_frequency = old_frequency
        output = new_frequency
        return
      endif
    else
      new_frequency = old_frequency - first_derivative/second_derivative
      if (abs(new_frequency-old_frequency)<frequency_convergence) then
        output = new_frequency
        return
      else
      endif
    endif
    
    ! If the frequency hasn't converged, update the frequency,
    !    and return to the start of the loop.
    old_frequency = new_frequency
  enddo
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_InitialFrequencies(this,input)
  implicit none
  
  class(InitialFrequencies), intent(out) :: this
  type(String),              intent(in)  :: input(:)
  
  integer,  allocatable :: subspace_ids(:)
  real(dp), allocatable :: frequencies(:)
  
  type(String), allocatable :: line(:)
  
  integer :: i,ialloc
  
  select type(this); type is(InitialFrequencies)
    allocate( subspace_ids(size(input)), &
            & frequencies(size(input)),  &
            & stat=ialloc); call err(ialloc)
    do i=1,size(input)
      line = split_line(input(i))
      subspace_ids(i) = int(line(2))
      frequencies(i) = dble(line(6))
    enddo
    this = InitialFrequencies(subspace_ids, frequencies)
  class default
    call err()
  end select
end subroutine

function write_InitialFrequencies(this) result(output)
  implicit none
  
  class(InitialFrequencies), intent(in) :: this
  type(String), allocatable             :: output(:)
  
  type(String) :: max_id
  
  integer :: i,ialloc
  
  select type(this); type is(InitialFrequencies)
    max_id = str(maxval(this%subspace_ids_))
    
    allocate(output(size(this)), stat=ialloc); call err(ialloc)
    do i=1,size(this)
      output(i) = 'Subspace '//left_pad(this%subspace_ids_(i),max_id)//' : &
         &frequency = '//this%frequencies_(i)
    enddo
  class default
    call err()
  end select
end function

function new_InitialFrequencies_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(InitialFrequencies) :: this
  
  call this%read(input)
end function

impure elemental function new_InitialFrequencies_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(InitialFrequencies)      :: this
  
  this = InitialFrequencies(str(input))
end function
end module