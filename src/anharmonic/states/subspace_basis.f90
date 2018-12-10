! ======================================================================
! A basis of states which spans a subspace.
! ======================================================================
module subspace_basis_module
  use common_module
  
  use subspace_state_module
  use state_helper_module
  use monomial_state_module
  use harmonic_state_module
  use braket_module
  use state_conversion_module
  use wavevector_basis_module
  implicit none
  
  private
  
  public :: SubspaceBasis
  public :: size
  
  ! All states spanning the subspace.
  type, extends(Stringsable) :: SubspaceBasis
    ! The maximum power of the monomial states.
    ! This is also the maximum occupation of the harmonic basis states.
    integer  :: maximum_power
    ! The expansion order of the potential.
    ! This is also the limit on coupling between states.
    integer :: expansion_order
    ! The ID and frequency of the subspace.
    integer  :: subspace_id
    real(dp) :: frequency
    ! The states, wavevector by wavevector.
    ! N.B. this only includes one wavevector from each symmetry-related set.
    type(WavevectorBasis), allocatable :: wavevectors(:)
  contains
    ! Set the frequency of the basis.
    procedure, public :: set_frequency => set_frequency_SubspaceBasis
    
    ! Print the harmonic ground-state wavefunction of the basis.
    procedure, public :: ground_state_wavefunction
    
    ! I/O.
    procedure, public :: read  => read_SubspaceBasis
    procedure, public :: write => write_SubspaceBasis
  end type
  
  interface SubspaceBasis
    module procedure new_SubspaceBasis
    module procedure new_SubspaceBasis_subspace
    module procedure new_SubspaceBasis_Strings
    module procedure new_SubspaceBasis_StringArray
  end interface
  
  interface size
    module procedure size_SubspaceBasis
  end interface
contains

! Constructors and size functions.
function new_SubspaceBasis(maximum_power,expansion_order,subspace_id, &
   & frequency,wavevectors) result(this)
  implicit none
  
  integer,               intent(in) :: maximum_power
  integer,               intent(in) :: expansion_order
  integer,               intent(in) :: subspace_id
  real(dp),              intent(in) :: frequency
  type(WavevectorBasis), intent(in) :: wavevectors(:)
  type(SubspaceBasis)               :: this
  
  this%maximum_power   = maximum_power
  this%expansion_order = expansion_order
  this%subspace_id     = subspace_id
  this%frequency       = frequency
  this%wavevectors     = wavevectors
end function

function size_SubspaceBasis(this) result(output)
  implicit none
  
  type(SubspaceBasis), intent(in) :: this
  integer                         :: output
  
  output = size(this%wavevectors)
end function

! Set the frequency of the basis.
impure elemental subroutine set_frequency_SubspaceBasis(this,frequency)
  implicit none
  
  class(SubspaceBasis), intent(inout) :: this
  real(dp),             intent(in)    :: frequency
  
  this%frequency = frequency
  call this%wavevectors%set_frequency(frequency)
end subroutine

! Generates states up to a given power, spanning the whole subspace.
function new_SubspaceBasis_subspace(subspace,frequency,modes,qpoints, &
   & supercell,maximum_power,potential_expansion_order,symmetries)    &
   & result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: frequency
  type(ComplexMode),        intent(in) :: modes(:)
  type(QpointData),         intent(in) :: qpoints(:)
  type(StructureData),      intent(in) :: supercell
  integer,                  intent(in) :: maximum_power
  integer,                  intent(in) :: potential_expansion_order
  type(SymmetryOperator),   intent(in) :: symmetries(:)
  type(SubspaceBasis)                  :: output
  
  type(WavevectorBasis), allocatable :: wavevectors(:)
  
  wavevectors = WavevectorBasis( subspace,                  &
                               & frequency,                 &
                               & modes,                     &
                               & qpoints,                   &
                               & maximum_power,             &
                               & potential_expansion_order, &
                               & symmetries                 )
  
  output = SubspaceBasis( maximum_power,             &
                        & potential_expansion_order, &
                        & subspace%id,               &
                        & frequency,                 &
                        & wavevectors                )
end function

! Returns the harmonic ground-state wavefunction for the basis.
function ground_state_wavefunction(this,subspace,supercell) result(output)
  implicit none
  
  class(SubspaceBasis),     intent(in) :: this
  type(DegenerateSubspace), intent(in) :: subspace
  type(StructureData),      intent(in) :: supercell
  type(String)                         :: output
  
  real(dp) :: mass
  
  real(dp)                  :: coefficient
  type(String), allocatable :: terms(:)
  
  integer :: i
  
  ! Calculate the (geometric) average mass.
  mass = product(supercell%atoms%mass())
  mass = mass**(1.0_dp/size(supercell%atoms))
  
  ! Calculate the coefficient.
  coefficient = 1
  terms = [String::]
  do i=1,size(subspace)
    if (subspace%mode_ids(i)==subspace%paired_ids(i)) then
      ! |0_i> = sqrt(sqrt(m*w/pi)) exp(- 1/2 N w (u_i)^2 )
      coefficient = coefficient * (mass*this%frequency/PI)**0.25_dp
      terms = [terms, 'u_'//subspace%mode_ids(i)//'^2']
    elseif (subspace%mode_ids(i)<subspace%paired_ids(i)) then
      ! |0_i,0_j> = sqrt(2*m*w/pi) exp(- N w |u_i|^2 )
      coefficient = coefficient * (2*mass*this%frequency/PI)**0.5_dp
      terms = [ terms,                                                    &
              & '2*u_'//subspace%mode_ids(i)//'*'//subspace%paired_ids(i) ]
    endif
  enddo
  
  output = coefficient//'*e^('               // &
     & (-supercell%sc_size*this%frequency/2) // &
     & '*('//join(terms,' + ')//'))'
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_SubspaceBasis(this,input)
  implicit none
  
  class(SubspaceBasis), intent(out) :: this
  type(String),         intent(in)  :: input(:)
  
  integer                            :: maximum_power
  integer                            :: expansion_order
  integer                            :: subspace_id
  real(dp)                           :: frequency
  type(WavevectorBasis), allocatable :: wavevectors(:)
  
  integer, allocatable :: starting_lines(:)
  integer, allocatable :: ending_lines(:)
  
  type(String), allocatable :: line(:)
  type(String), allocatable :: wavevector_lines(:)
  
  integer :: i,ialloc
  
  select type(this); type is(SubspaceBasis)
    line = split_line(input(1))
    maximum_power = int(line(4))
    
    line = split_line(input(2))
    expansion_order = int(line(4))
    
    line = split_line(input(3))
    subspace_id = int(line(3))
    
    line = split_line(input(4))
    frequency = dble(line(3))
    
    starting_lines = [integer::]
    do i=5,size(input)
      line = split_line(input(i))
      if (size(line)>0) then
        if (line(1)=='Wavevector') then
          starting_lines = [starting_lines, i]
        endif
      endif
    enddo
    
    ending_lines = [starting_lines(2:)-1, size(input)]
    
    allocate(wavevectors(size(starting_lines)), stat=ialloc); call err(ialloc)
    do i=1,size(starting_lines)
      wavevector_lines = [input(:4), input(starting_lines(i):ending_lines(i))]
      wavevectors(i) = WavevectorBasis(wavevector_lines)
    enddo
    
    this = SubspaceBasis( maximum_power,   &
                        & expansion_order, &
                        & subspace_id,     &
                        & frequency,       &
                        & wavevectors      )
  class default
    call err()
  end select
end subroutine

function write_SubspaceBasis(this) result(output)
  implicit none
  
  class(SubspaceBasis), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  type(String), allocatable :: wavevector_strings(:)
  
  integer :: i
  
  select type(this); type is(SubspaceBasis)
    output = [ 'Maximum power   : '//this%maximum_power,   &
             & 'Expansion order : '//this%expansion_order, &
             & 'Subspace        : '//this%subspace_id,     &
             & 'Frequency       : '//this%frequency        ]
    do i=1,size(this%wavevectors)
      wavevector_strings = str(this%wavevectors(i))
      output = [ output, wavevector_strings(5:) ]
    enddo
  class default
    call err()
  end select
end function

function new_SubspaceBasis_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(SubspaceBasis)      :: this
  
  call this%read(input)
end function

impure elemental function new_SubspaceBasis_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(SubspaceBasis)           :: this
  
  this = SubspaceBasis(str(input))
end function
end module
