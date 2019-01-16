! ======================================================================
! A wrapped polymorphic pointer to a potential.
! ======================================================================
! Wraps all of PotentialData's methods,
!    calling them on the pointed-to potential.
! See example module below for how to use this type.
module potential_pointer_module
  use common_module
  
  use states_module
  
  use anharmonic_common_module
  use polynomial_module
  implicit none
  
  private
  
  public :: PotentialPointer
  public :: generate_subspace_potentials
  
  type, extends(PotentialData) :: PotentialPointer
    type(String),                      private :: representation_
    class(PotentialData), allocatable, private :: potential_
  contains
    procedure, public :: generate_sampling_points => &
       & generate_sampling_points_PotentialPointer
    procedure, public :: generate_potential => &
       & generate_potential_PotentialPointer
    
    procedure, public :: zero_energy => zero_energy_PotentialPointer
    
    procedure, public :: energy_RealModeDisplacement => &
                       & energy_RealModeDisplacement_PotentialPointer
    procedure, public :: energy_ComplexModeDisplacement => &
                       & energy_ComplexModeDisplacement_PotentialPointer
    procedure, public :: force_RealModeDisplacement => &
                       & force_RealModeDisplacement_PotentialPointer
    procedure, public :: force_ComplexModeDisplacement => &
                       & force_ComplexModeDisplacement_PotentialPointer
    
    procedure, public :: braket => braket_PotentialPointer
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_PotentialPointer
    
    procedure, public :: iterate_damped => &
                       & iterate_damped_PotentialPointer
    procedure, public :: iterate_pulay => &
                       & iterate_pulay_PotentialPointer
    
    procedure, private :: check => check_PotentialPointer
    
    ! I/O.
    procedure, public :: read  => read_PotentialPointer
    procedure, public :: write => write_PotentialPointer
  end type
  
  interface PotentialPointer
    module procedure new_PotentialPointer
    module procedure new_PotentialPointer_Strings
    module procedure new_PotentialPointer_StringArray
  end interface
contains

! Construct a PotentialPointer from any type which extends PotentialData.
impure elemental function new_PotentialPointer(potential) result(this)
  implicit none
  
  class(PotentialData), intent(in) :: potential
  type(PotentialPointer)           :: this
  
  integer :: ialloc
  
  select type(potential); type is(PotentialPointer)
    this = potential
  type is(PolynomialPotential)
    this%representation_ = 'polynomial'
    allocate( this%potential_, source=potential, &
            & stat=ialloc); call err(ialloc)
  class default
    call err()
  end select
end function

! Checks that the pointer has been allocated before it is used.
impure elemental subroutine check_PotentialPointer(this)
  implicit none
  
  class(PotentialPointer), intent(in) :: this
  
  if (.not. allocated(this%potential_)) then
    call print_line(CODE_ERROR//': Trying to use a PotentialPointer before &
       &it has been allocated.')
    call err()
  endif
end subroutine

! Wrappers for all of PotentialData's methods.
subroutine generate_sampling_points_PotentialPointer(this,anharmonic_data, &
   & sampling_points_dir,calculation_writer,logfile)
  implicit none
  
  class(PotentialPointer), intent(inout) :: this
  type(AnharmonicData),    intent(in)    :: anharmonic_data
  type(String),            intent(in)    :: sampling_points_dir
  type(CalculationWriter), intent(inout) :: calculation_writer
  type(OFile),             intent(inout) :: logfile
  
  call this%check()
  
  call this%potential_%generate_sampling_points( anharmonic_data,     &
                                               & sampling_points_dir, &
                                               & calculation_writer,  &
                                               & logfile              )
end subroutine

subroutine generate_potential_PotentialPointer(this,anharmonic_data,     &
   & weighted_energy_force_ratio,sampling_points_dir,calculation_reader, &
   & logfile)
  implicit none
  
  class(PotentialPointer), intent(inout) :: this
  type(AnharmonicData),    intent(in)    :: anharmonic_data
  real(dp),                intent(in)    :: weighted_energy_force_ratio
  type(String),            intent(in)    :: sampling_points_dir
  type(CalculationReader), intent(inout) :: calculation_reader
  type(OFile),             intent(inout) :: logfile
  
  call this%check()
  
  call this%potential_%generate_potential( anharmonic_data,              &
                                         & weighted_energy_force_ratio,  &
                                         & sampling_points_dir,          &
                                         & calculation_reader,           &
                                         & logfile                       )
end subroutine

impure elemental subroutine zero_energy_PotentialPointer(this)
  implicit none
  
  class(PotentialPointer), intent(inout) :: this
  
  call this%check()
  
  call this%potential_%zero_energy()
end subroutine

impure elemental function energy_RealModeDisplacement_PotentialPointer(this, &
   & displacement) result(output)
  implicit none
  
  class(PotentialPointer),    intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  real(dp)                               :: output
  
  call this%check()
  
  output = this%potential_%energy(displacement)
end function

impure elemental function energy_ComplexModeDisplacement_PotentialPointer( &
   & this,displacement) result(output)
  implicit none
  
  class(PotentialPointer),       intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  complex(dp)                               :: output
  
  call this%check()
  
  output = this%potential_%energy(displacement)
end function

impure elemental function force_RealModeDisplacement_PotentialPointer(this, &
   & displacement) result(output)
  implicit none
  
  class(PotentialPointer),    intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  type(RealModeForce)                    :: output
  
  call this%check()
  
  output = this%potential_%force(displacement)
end function

impure elemental function force_ComplexModeDisplacement_PotentialPointer( &
   & this,displacement) result(output)
  implicit none
  
  class(PotentialPointer),       intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  type(ComplexModeForce)                    :: output
  
  call this%check()
  
  output = this%potential_%force(displacement)
end function

function braket_PotentialPointer(this,bra,ket,subspace,anharmonic_data) &
   & result(output)
  implicit none
  
  class(PotentialPointer),  intent(in) :: this
  class(SubspaceState),     intent(in) :: bra
  class(SubspaceState),     intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  class(PotentialData), allocatable    :: output
  
  call this%check()
  
  output = this%potential_%braket(bra,ket,subspace,anharmonic_data)
end function

function harmonic_expectation_PotentialPointer(this,frequency, &
   & thermal_energy,no_states,subspace,anharmonic_data) result(output)
  implicit none
  
  class(PotentialPointer),  intent(in) :: this
  real(dp),                 intent(in) :: frequency
  real(dp),                 intent(in) :: thermal_energy
  integer,                  intent(in) :: no_states
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp)                             :: output
  
  call this%check()
  
  output = this%potential_%harmonic_expectation( frequency,      &
                                               & thermal_energy, &
                                               & no_states,      &
                                               & subspace,       &
                                               & anharmonic_data )
end function

! ----------------------------------------------------------------------
! Generates the next iteration of the potential, either following a damped
!    iterative scheme or a pulay scheme.
! ----------------------------------------------------------------------
impure elemental function iterate_damped_PotentialPointer(this,new_potential, &
   & damping,anharmonic_data) result(output)
  implicit none
  
  class(PotentialPointer), intent(in) :: this
  class(PotentialData),    intent(in) :: new_potential
  real(dp),                intent(in) :: damping
  type(AnharmonicData),    intent(in) :: anharmonic_data
  class(PotentialData), allocatable   :: output
  
  class(PotentialData), allocatable :: potential
  
  call this%check()
  
  select type(new_potential); type is(PotentialPointer)
    call new_potential%check()
    potential = new_potential%potential_
  class default
    potential = new_potential
  end select
  
  output = this%potential_%iterate_damped( potential,      &
                                         & damping,        &
                                         & anharmonic_data )
end function

function iterate_pulay_PotentialPointer(this,input_potentials, &
   & output_potentials,anharmonic_data) result(output)
  implicit none
  
  class(PotentialPointer), intent(in) :: this
  class(PotentialData),    intent(in) :: input_potentials(:)
  class(PotentialData),    intent(in) :: output_potentials(:)
  type(AnharmonicData),    intent(in) :: anharmonic_data
  class(PotentialData), allocatable   :: output
  
  class(PotentialData), allocatable :: in_potentials(:)
  class(PotentialData), allocatable :: out_potentials(:)
  
  integer :: i
  
  call this%check()
  
  select type(input_potentials); type is(PotentialPointer)
    call input_potentials%check()
    in_potentials = [( input_potentials(i)%potential_, &
                     & i=1,                            &
                     & size(input_potentials)          )]
  class default
    in_potentials = input_potentials
  end select
  
  select type(output_potentials); type is(PotentialPointer)
    call output_potentials%check()
    out_potentials = [( output_potentials(i)%potential_, &
                      & i=1,                             &
                      & size(output_potentials)          )]
  class default
    out_potentials = output_potentials
  end select
  
  output = this%potential_%iterate_pulay( in_potentials,  &
                                        & out_potentials, &
                                        & anharmonic_data )
end function

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
! I/O.
! ----------------------------------------------------------------------
subroutine read_PotentialPointer(this,input)
  implicit none
  
  class(PotentialPointer), intent(out) :: this
  type(String),            intent(in)  :: input(:)
  
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  
  select type(this); type is(PotentialPointer)
    line = split_line(input(1))
    representation = line(3)
    if (representation=='polynomial') then
      this = PotentialPointer(PolynomialPotential(input(2:)))
    else
      call print_line( 'Unrecognised potential representation: '// &
                     & representation                              )
    endif
  class default
    call err()
  end select
end subroutine

function write_PotentialPointer(this) result(output)
  implicit none
  
  class(PotentialPointer), intent(in) :: this
  type(String), allocatable           :: output(:)
  
  select type(this); type is(PotentialPointer)
    output = [ 'Potential representation: '//this%representation_, &
             & str(this%potential_)                                ]
  end select
end function

function new_PotentialPointer_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(PotentialPointer)   :: this
  
  call this%read(input)
end function

impure elemental function new_PotentialPointer_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(PotentialPointer)        :: this
  
  this = PotentialPointer(str(input))
end function
end module

! ======================================================================
! An example module, showing how to use PotentialData and PotentialPointer.
! ======================================================================
module potential_example_module
  use common_module
  
  use anharmonic_common_module
  use states_module
  
  use anharmonic_data_module
  use potential_module
  use potential_pointer_module
  implicit none
  
  private
  
  public :: potential_example_subroutine
  
  type, extends(PotentialData) :: PotentialDataExample
    type(String) :: example_contents
  contains
    procedure, public :: generate_sampling_points => &
       & generate_sampling_points_PotentialDataExample
    procedure, public :: generate_potential => &
       & generate_potential_PotentialDataExample
    
    procedure, public :: zero_energy => zero_energy_PotentialDataExample
    
    procedure, public :: energy_RealModeDisplacement => &
                       & energy_RealModeDisplacement_PotentialDataExample
    procedure, public :: energy_ComplexModeDisplacement => &
                       & energy_ComplexModeDisplacement_PotentialDataExample
    procedure, public :: force_RealModeDisplacement => &
                       & force_RealModeDisplacement_PotentialDataExample
    procedure, public :: force_ComplexModeDisplacement => &
                       & force_ComplexModeDisplacement_PotentialDataExample
    
    procedure, public :: braket => braket_PotentialDataExample
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_PotentialDataExample
    
    procedure, public :: iterate_damped => &
                       & iterate_damped_PotentialDataExample
    procedure, public :: iterate_pulay => &
                       & iterate_pulay_PotentialDataExample
    
    procedure, public :: read  => read_PotentialDataExample
    procedure, public :: write => write_PotentialDataExample
  end type
  
  interface PotentialDataExample
    module procedure new_PotentialDataExample
    module procedure new_PotentialDataExample_Strings
    module procedure new_PotentialDataExample_StringArray
  end interface
contains

! --------------------------------------------------
! Constructor for example class.
! --------------------------------------------------
! This is where any PotentialDataExample-specific data is input.
function new_PotentialDataExample(example_contents) result(this)
  implicit none
  
  type(String), intent(in)   :: example_contents
  type(PotentialDataExample) :: this
  
  this%example_contents = example_contents
end function

! --------------------------------------------------
! Overloads of PotentialData's methods.
! --------------------------------------------------
subroutine generate_sampling_points_PotentialDataExample(this, &
   & anharmonic_data,sampling_points_dir,calculation_writer,logfile)
  implicit none
  
  class(PotentialDataExample), intent(inout) :: this
  type(AnharmonicData),        intent(in)    :: anharmonic_data
  type(String),                intent(in)    :: sampling_points_dir
  type(CalculationWriter),     intent(inout) :: calculation_writer
  type(OFile),                 intent(inout) :: logfile
  
  call print_line('PotentialDataExample: generating sampling points.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to generate sampling points goes here.
end subroutine

subroutine generate_potential_PotentialDataExample(this,anharmonic_data, &
   & weighted_energy_force_ratio,sampling_points_dir,calculation_reader, &
   & logfile)
  implicit none
  
  class(PotentialDataExample), intent(inout) :: this
  type(AnharmonicData),        intent(in)    :: anharmonic_data
  real(dp),                    intent(in)    :: weighted_energy_force_ratio
  type(String),                intent(in)    :: sampling_points_dir
  type(CalculationReader),     intent(inout) :: calculation_reader
  type(OFile),                 intent(inout) :: logfile
  
  call print_line('PotentialDataExample: generating potential.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to generate sampling points goes here.
end subroutine

impure elemental subroutine zero_energy_PotentialDataExample(this)
  implicit none
  
  class(PotentialDataExample), intent(inout) :: this
  
  call print_line('PotentialDataExample: zeroing energy.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to zero the energy (s.t. undisplaced_energy()=0) goes here.
end subroutine

impure elemental function energy_RealModeDisplacement_PotentialDataExample( &
   & this,displacement) result(output)
  implicit none
  
  class(potentialDataExample), intent(in) :: this
  type(RealModeDisplacement),  intent(in) :: displacement
  real(dp)                                :: output
  
  call print_line('PotentialDataExample: calculating energy.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to calculate energies at real displacements goes here.
end function

impure elemental function energy_ComplexModeDisplacement_PotentialDataExample(&
   & this,displacement) result(output)
  implicit none
  
  class(potentialDataExample),    intent(in) :: this
  type(ComplexModeDisplacement),  intent(in) :: displacement
  complex(dp)                                :: output
  
  call print_line('PotentialDataExample: calculating energy.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to calculate energies at complex displacements goes here.
end function

impure elemental function force_RealModeDisplacement_PotentialDataExample( &
   & this,displacement) result(output)
  implicit none
  
  class(potentialDataExample), intent(in) :: this
  type(RealModeDisplacement),  intent(in) :: displacement
  type(RealModeForce)                     :: output
  
  call print_line('PotentialDataExample: calculating force.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to calculate forces at real displacements goes here.
end function

impure elemental function force_ComplexModeDisplacement_PotentialDataExample( &
   & this,displacement) result(output)
  implicit none
  
  class(potentialDataExample),    intent(in) :: this
  type(ComplexModeDisplacement),  intent(in) :: displacement
  type(ComplexModeForce)                     :: output
  
  call print_line('PotentialDataExample: calculating force.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to calculate forces at complex displacements goes here.
end function

function braket_PotentialDataExample(this,bra,ket,subspace,anharmonic_data) &
   & result(output)
  implicit none
  
  class(PotentialDataExample), intent(in) :: this
  class(SubspaceState),        intent(in) :: bra
  class(SubspaceState),        intent(in) :: ket
  type(DegenerateSubspace),    intent(in) :: subspace
  type(AnharmonicData),        intent(in) :: anharmonic_data
  class(PotentialData), allocatable       :: output
  
  call print_line('PotentialDataExample: evaluating <bra|potential|ket>.')
  
  ! Code to integrate this potential between <bra| and |ket> goes here.
end function

function harmonic_expectation_PotentialDataExample(this,frequency, &
   & thermal_energy,no_states,subspace,anharmonic_data) result(output)
  implicit none
  
  class(PotentialDataExample), intent(in) :: this
  real(dp),                    intent(in) :: frequency
  real(dp),                    intent(in) :: thermal_energy
  integer,                     intent(in) :: no_states
  type(DegenerateSubspace),    intent(in) :: subspace
  type(AnharmonicData),        intent(in) :: anharmonic_data
  real(dp)                                :: output
  
  call print_line('PotentialDataExample: evaluating harmonic expectation of &
     &<V>.')
  
  ! Code to calculate thermal harmonic expectation goes here.
end function

impure elemental function iterate_damped_PotentialDataExample(this, &
   & new_potential,damping,anharmonic_data) result(output)
  implicit none
  
  class(PotentialDataExample), intent(in) :: this
  class(PotentialData),        intent(in) :: new_potential
  real(dp),                    intent(in) :: damping
  type(AnharmonicData),        intent(in) :: anharmonic_data
  class(PotentialData), allocatable       :: output
  
  ! Code to generate the potential which is equal to:
  !    (1-damping)*this + damping*new_potential
  ! goes here.
end function

function iterate_pulay_PotentialDataExample(this,input_potentials, &
   & output_potentials,anharmonic_data) result(output)
  implicit none
  
  class(PotentialDataExample), intent(in) :: this
  class(PotentialData),        intent(in) :: input_potentials(:)
  class(PotentialData),        intent(in) :: output_potentials(:)
  type(AnharmonicData),        intent(in) :: anharmonic_data
  class(PotentialData), allocatable       :: output
  
  ! Code to generate a new potential using a Pulay scheme goes here.
  ! See the PolynomialPotential Pulay scheme for an example.
end function

! --------------------------------------------------
! I/O.
! --------------------------------------------------
subroutine read_PotentialDataExample(this,input)
  implicit none
  
  class(PotentialDataExample), intent(out) :: this
  type(String),                intent(in)  :: input(:)
  
  ! Code to read potential from strings goes here.
end subroutine

function write_PotentialDataExample(this) result(output)
  implicit none
  
  class(PotentialDataExample), intent(in) :: this
  type(String), allocatable               :: output(:)
  
  ! Code to write potential to strings goes here.
end function

function new_PotentialDataExample_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)   :: input(:)
  type(PotentialDataExample) :: this
  
  call this%read(input)
end function

impure elemental function new_PotentialDataExample_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(PotentialDataExample)    :: this
  
  this = PotentialDataExample(str(input))
end function

! --------------------------------------------------
! The class in use.
! --------------------------------------------------
subroutine potential_example_subroutine()
  implicit none
  
  ! A polymorphic pointer, which can store an object of any type which
  !    extends PotentialData.
  type(PotentialPointer) :: potential
  
  ! An example variable with which to initialise a PotentialDataExample.
  type(String) :: example_contents
  
  ! Variables for generate_sampling points.
  type(AnharmonicData)    :: anharmonic_data
  type(String)            :: sampling_points_dir
  type(CalculationWriter) :: calculation_writer
  type(OFile)             :: logfile
  
  ! Variables for generate_potential.
  real(dp)                :: weighted_energy_force_ratio
  type(CalculationReader) :: calculation_reader
  
  ! Variables for energy and force.
  type(RealModeDisplacement) :: real_displacement
  real(dp)                   :: real_energy
  type(RealModeForce)        :: real_force
  type(ComplexModeDisplacement) :: complex_displacement
  complex(dp)                   :: complex_energy
  type(ComplexModeForce)        :: complex_force
  
  ! Variables for integrating potential.
  type(DegenerateSubspace) :: subspace
  type(MonomialState)      :: state_1
  type(MonomialState)      :: state_2
  
  ! Files.
  type(OFile) :: output_file
  type(IFile) :: input_file
  
  ! Set the pointer to point to a PotentialDataExample type.
  ! This is where any PotentialDataExample-specific data is input,
  !    in this case the variable example_contents.
  example_contents = 'example'
  potential = PotentialPointer(PotentialDataExample(example_contents))
  
  ! Now PotentialData's methods can be called.
  ! They will all be forwareded to the PotentialDataExample instance.
  
  ! Generates sampling points, in a manner specific to the representation.
  call potential%generate_sampling_points( anharmonic_data,     &
                                         & sampling_points_dir, &
                                         & calculation_writer,  &
                                         & logfile              )
  
  ! Code to run electronic structure goes here.
  
  ! Generates the potential, in a manner specific to the representation.
  call potential%generate_potential( anharmonic_data,             &
                                   & weighted_energy_force_ratio, &
                                   & sampling_points_dir,         &
                                   & calculation_reader,          &
                                   & logfile                      )
  
  ! Once the potential has been generated, it can be used to calculate
  !    energies and forces.
  real_energy = potential%energy(real_displacement)
  real_force  = potential%force(real_displacement)
  complex_energy = potential%energy(complex_displacement)
  complex_force  = potential%force(complex_displacement)
  
  ! The potential can also be integrated between two states.
  potential = PotentialPointer(                                   &
     & potential%braket(state_1,state_2,subspace,anharmonic_data) )
  
  ! The potential can be written to and read from file using the potential
  !    pointer's methods.
  output_file = OFile('example_potential.file')
  call output_file%print_lines(potential)
  
  input_file = IFile('example_potential.file')
  potential = PotentialPointer(input_file%lines())
end subroutine
end module
