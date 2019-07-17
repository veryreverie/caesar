! ======================================================================
! Finds the frequencies of the self-consistent harmonic potential which
!    minimises the free energy at zero temperature.
! ======================================================================
module initial_frequencies_module
  use common_module
  
  use states_module
  
  use anharmonic_common_module
  use potentials_module
  
  use generate_subspace_potentials_module
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
function new_InitialFrequencies_PotentialData(potential,anharmonic_data,   &
   & frequency_convergence,no_converged_calculations,max_pulay_iterations, &
   & pre_pulay_iterations,pre_pulay_damping) result(output)
  implicit none
  
  class(PotentialData), intent(in) :: potential
  type(AnharmonicData), intent(in) :: anharmonic_data
  real(dp),             intent(in) :: frequency_convergence
  integer,              intent(in) :: no_converged_calculations
  integer,              intent(in) :: max_pulay_iterations
  integer,              intent(in) :: pre_pulay_iterations
  real(dp),             intent(in) :: pre_pulay_damping
  type(InitialFrequencies)         :: output
  
  type(DegenerateSubspace), allocatable :: subspaces(:)
  real(dp),                 allocatable :: frequencies(:)
  type(PotentialPointer),   allocatable :: subspace_potentials(:)
  type(FullSubspaceBasis),  allocatable :: subspace_bases(:)
  type(WavevectorStates),   allocatable :: subspace_states(:)
  
  real(dp), allocatable :: max_differences(:)
  
  real(dp), allocatable :: input_frequency(:)
  real(dp), allocatable :: output_frequency(:)
  
  type(PulaySolver) :: solver
  
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
  
  ! For each subspace, generate a basis containing only the |0> state.
  ! Generate bases at first guess frequencies.
  allocate(subspace_bases(size(subspaces)), stat=ialloc); call err(ialloc)
  do i=1,size(subspaces)
    subspace_bases(i) = FullSubspaceBasis(                                 &
       & subspace                  = subspaces(i),                         &
       & frequency                 = frequencies(i),                       &
       & modes                     = anharmonic_data%complex_modes,        &
       & qpoints                   = anharmonic_data%qpoints,              &
       & supercell                 = anharmonic_data%anharmonic_supercell, &
       & maximum_power             = 0,                                    &
       & potential_expansion_order = 0,                                    &
       & symmetries                = anharmonic_data%structure%symmetries  )
  enddo
  
  ! Calculate the first guess ground states.
  subspace_states = [(                                           &
     & WavevectorStates(                                         &
     &    subspace_bases(i)%initial_states( subspaces(i),        &
     &                                      anharmonic_data ) ), &
     & i=1,                                                      &
     & size(subspaces)                                           )]
  
  ! Find self-consistent frequencies which minimise the energy.
  solver = PulaySolver( pre_pulay_iterations, &
                      & pre_pulay_damping,    &
                      & max_pulay_iterations, &
                      & frequencies,          &
                      & bound_at_zero=.true.  )
  i = 0
  allocate(max_differences(0), stat=ialloc); call err(ialloc)
  do
    input_frequency = solver%get_x()
    
    ! Calculate single-subspace mean-field potentials, defined as
    !    V_i = (prod_{j/=i}<j|)V(prod_{j/=i}|j>).
    call subspace_bases%set_frequency(input_frequency)
    
    subspace_potentials = PotentialPointer(               &
       & generate_subspace_potentials( potential,         &
       &                               subspaces,         &
       &                               subspace_bases,    &
       &                               subspace_states,   &
       &                               anharmonic_data  ) )
    
    ! Calculate updated frequencies, each of which minimises the free energy
    !    of the corresponding subspace potential.
    output_frequency = optimise_frequency( subspace_potentials,  &
                                         & subspaces,            &
                                         & subspace_bases,       &
                                         & subspace_states,      &
                                         & anharmonic_data,      &
                                         & frequency_convergence )
    
    ! Calculate the maximum error in self-consistency,
    !    for convergence purposes.
    max_differences = [ max_differences,                              &
                      & maxval(abs(output_frequency-input_frequency)) ]
    
    ! Increment the loop counter.
    i = i+1
    
    ! Check whether the frequencies have converged by the normal convergence
    !    condition.
    if (i>=no_converged_calculations) then
      if (all( max_differences(i-no_converged_calculations+1:) &
           & < frequency_convergence                           )) then
        frequencies = output_frequency
        exit
      endif
    endif
    
    call print_line('Frequency self-consistency step '//i// &
                   & '. Maximum error: '//                  &
                   & max_differences(i)                     &
                   & //' (Ha).'                             )
    
    call solver%set_f(output_frequency)
  enddo
  
  output = InitialFrequencies(subspaces%id, frequencies)
end function

! Returns whether or not two steps are converged relative to one another.
impure elemental function steps_converged(this,that,frequency_convergence) &
   & result(output)
  implicit none
  
  type(RealVector), intent(in) :: this
  type(RealVector), intent(in) :: that
  real(dp),         intent(in) :: frequency_convergence
  logical                      :: output
  
  output = all(abs(dble(this-that))<frequency_convergence)
end function

! Find the frequency which minimises energy in a single subspace.
impure elemental function optimise_frequency(potential,subspace,           &
   & subspace_basis,subspace_states,anharmonic_data,frequency_convergence) &
   & result(output)
  implicit none
  
  class(PotentialData),     intent(in) :: potential
  type(DegenerateSubspace), intent(in) :: subspace
  type(FullSubspaceBasis),  intent(in) :: subspace_basis
  type(WavevectorStates),   intent(in) :: subspace_states
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp),                 intent(in) :: frequency_convergence
  real(dp)                             :: output
  
  type(NewtonRaphson) :: solver
  
  real(dp) :: frequencies(3)
  real(dp) :: energies(3)
  
  type(FullSubspaceBasis) :: new_basis
  
  integer :: i
  
  new_basis = subspace_basis
  
  solver = NewtonRaphson(                                     &
     & starting_value        = subspace_basis%frequency,      &
     & finite_displacement   = 0.01_dp*frequency_convergence, &
     & convergence_threshold = frequency_convergence,         &
     & lower_bound           = 0.0_dp                         )
  do
    frequencies = solver%get_inputs()
    
    do i=1,3
      call new_basis%set_frequency(frequencies(i))
      energies(i) = potential_energy( subspace_states%states(1),    &
                &                     potential,                    &
                &                     subspace,                     &
                &                     new_basis,                    &
                &                     anharmonic_data            )  &
                & + new_basis%kinetic_energy(                       &
                &      bra             = subspace_states%states(1), &
                &      subspace        = subspace,                  &
                &      anharmonic_data = anharmonic_data            )
    enddo
    
    call solver%set_outputs(energies)
    
    if (solver%converged()) then
      output = solver%solution()
      exit
    endif
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
