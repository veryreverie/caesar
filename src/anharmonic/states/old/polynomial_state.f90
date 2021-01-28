! ======================================================================
! A sum of monomial states and coefficients.
! ======================================================================
module caesar_polynomial_state_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  
  use caesar_monomial_state_module
  implicit none
  
  private
  
  public :: startup_polynomial_state
  
  public :: PolynomialState
  
  public :: size
  
  type, extends(BasisState) :: PolynomialState
    type(MonomialState), allocatable :: states(:)
    real(dp),            allocatable :: coefficients(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_PolynomialState
    
    procedure, public :: modes => modes_PolynomialState
    
    procedure, public :: inner_product => &
                       & inner_product_PolynomialState
    
    procedure, public :: braket_ComplexMonomial => &
                       & braket_ComplexMonomial_PolynomialState
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_PolynomialState
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_PolynomialState
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_PolynomialState
    
    procedure, public :: change_modes => change_modes_PolynomialState
    
    ! I/O.
    procedure, public :: read  => read_PolynomialState
    procedure, public :: write => write_PolynomialState
  end type
  
  interface PolynomialState
    module procedure new_PolynomialState
    module procedure new_PolynomialState_BasisState
    module procedure new_PolynomialState_Strings
    module procedure new_PolynomialState_StringArray
  end interface
  
  interface size
    module procedure size_PolynomialState
  end interface
contains

! Startup procedure.
subroutine startup_polynomial_state()
  implicit none
  
  type(PolynomialState) :: state
  
  call state%startup()
end subroutine

! Constructor and size function
function new_PolynomialState(subspace_id,states,coefficients) result(this)
  implicit none
  
  integer                          :: subspace_id
  type(MonomialState), allocatable :: states(:)
  real(dp),            allocatable :: coefficients(:)
  type(PolynomialState)            :: this
  
  if (size(states)==0) then
    call print_line(CODE_ERROR//': No states.')
    call err()
  elseif (size(states)/=size(coefficients)) then
    call print_line(CODE_ERROR//': States and coefficients do not match.')
    call err()
  endif
  
  this%subspace_id  = subspace_id
  this%states       = states
  this%coefficients = coefficients
end function

recursive function new_PolynomialState_BasisState(input) result(this)
  implicit none
  
  class(BasisState), intent(in) :: input
  type(PolynomialState)         :: this
  
  select type(input); type is(PolynomialState)
    this = input
  type is(BasisStatePointer)
    this = PolynomialState(input%state())
  class default
    call err()
  end select
end function

function size_PolynomialState(this) result(output)
  implicit none
  
  type(PolynomialState), intent(in) :: this
  integer                           :: output
  
  output = size(this%states)
end function

! ----------------------------------------------------------------------
! Type representation.
! ----------------------------------------------------------------------
impure elemental function representation_PolynomialState() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'polynomial'
end function

! ----------------------------------------------------------------------
! Returns the modes spanned by the state.
! ----------------------------------------------------------------------
function modes_PolynomialState(this) result(output)
  implicit none
  
  class(PolynomialState), intent(in) :: this
  integer, allocatable               :: output(:)
  
  output = this%states(1)%modes()
end function

! ----------------------------------------------------------------------
! BasisState methods.
! ----------------------------------------------------------------------
impure elemental function inner_product_PolynomialState(this, &
   & ket,subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(PolynomialState),   intent(in)           :: this
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(PolynomialState) :: polynomial_ket
  
  integer :: i,j
  
  if (present(ket)) then
    polynomial_ket = PolynomialState(ket)
    output = 0.0_dp
    do i=1,size(this)
      do j=1,size(polynomial_ket)
        output = output                                                    &
             & + this%states(i)%inner_product( polynomial_ket%states(j),   &
             &                                 subspace,                   &
             &                                 subspace_basis,             &
             &                                 anharmonic_data           ) &
             & * this%coefficients(i)                                      &
             & * polynomial_ket%coefficients(j)
      enddo
    enddo
  else
    output = 0.0_dp
    do i=1,size(this)
      do j=1,size(this)
        output = output                                          &
             & + this%states(i)%inner_product( this%states(j),   &
             &                                 subspace,         &
             &                                 subspace_basis,   &
             &                                 anharmonic_data ) &
             & * this%coefficients(i)                            &
             & * this%coefficients(j)
      enddo
    enddo
  endif
end function

impure elemental function braket_ComplexMonomial_PolynomialState(this, &
   & monomial,ket,subspace,subspace_basis,anharmonic_data,qpoint)      &
   & result(output)
  implicit none
  
  class(PolynomialState),   intent(in)           :: this
  type(ComplexMonomial),    intent(in)           :: monomial
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(QpointData),         intent(in), optional :: qpoint
  type(ComplexMonomial)                          :: output
  
  type(PolynomialState) :: polynomial_ket
  type(ComplexMonomial) :: integrated_monomial
  
  integer :: i,j
  
  if (present(ket)) then
    polynomial_ket = PolynomialState(ket)
    do i=1,size(this)
      do j=1,size(polynomial_ket)
        integrated_monomial = this%states(i)%braket(         &
                          &      monomial,                   &
                          &      polynomial_ket%states(j),   &
                          &      subspace,                   &
                          &      subspace_basis,             &
                          &      anharmonic_data,            &
                          &      qpoint                    ) &
                          & * this%coefficients(i)           &
                          & * polynomial_ket%coefficients(j)
        if (i==1 .and. j==1) then
          output = integrated_monomial
        else
          output%coefficient = output%coefficient &
                           & + integrated_monomial%coefficient
        endif
      enddo
    enddo
  else
    do i=1,size(this)
      do j=1,size(this)
        integrated_monomial = this%states(i)%braket( monomial,          &
                          &                          this%states(j),    &
                          &                          subspace,          &
                          &                          subspace_basis,    &
                          &                          anharmonic_data,   &
                          &                          qpoint           ) &
                          & * this%coefficients(i)                      &
                          & * this%coefficients(j)
        if (i==1 .and. j==1) then
          output = integrated_monomial
        else
          output%coefficient = output%coefficient &
                           & + integrated_monomial%coefficient
        endif
      enddo
    enddo
  endif
end function

impure elemental function kinetic_energy_PolynomialState(this,ket, &
   & subspace,subspace_basis,anharmonic_data,qpoint) result(output)
  implicit none
  
  class(PolynomialState),   intent(in)           :: this
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(QpointData), intent(in), optional :: qpoint
  real(dp)                                       :: output
  
  type(PolynomialState) :: polynomial_ket
  
  integer :: i,j
  
  if (present(ket)) then
    polynomial_ket = PolynomialState(ket)
    output = 0.0_dp
    do i=1,size(this)
      do j=1,size(polynomial_ket)
        output = output                                                     &
             & + this%states(i)%kinetic_energy( polynomial_ket%states(j),   &
             &                                  subspace,                   &
             &                                  subspace_basis,             &
             &                                  anharmonic_data,qpoint           ) &
             & * this%coefficients(i)                                       &
             & * polynomial_ket%coefficients(j)
      enddo
    enddo
  else
    output = 0.0_dp
    do i=1,size(this)
      do j=1,size(this)
        output = output                                           &
             & + this%states(i)%kinetic_energy( this%states(j),   &
             &                                  subspace,         &
             &                                  subspace_basis,   &
             &                                  anharmonic_data,qpoint ) &
             & * this%coefficients(i)                             &
             & * this%coefficients(j)
      enddo
    enddo
  endif
end function

impure elemental function harmonic_potential_energy_PolynomialState( &
   & this,ket,subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(PolynomialState),   intent(in)           :: this
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  type(PolynomialState) :: polynomial_ket
  
  integer :: i,j
  
  if (present(ket)) then
    polynomial_ket = PolynomialState(ket)
    output = 0.0_dp
    do i=1,size(this)
      do j=1,size(polynomial_ket)
        output = output                                    &
             & + this%states(i)%harmonic_potential_energy( &
             &                 polynomial_ket%states(j),   &
             &                 subspace,                   &
             &                 subspace_basis,             &
             &                 anharmonic_data           ) &
             & * this%coefficients(i)                      &
             & * polynomial_ket%coefficients(j)
      enddo
    enddo
  else
    output = 0.0_dp
    do i=1,size(this)
      do j=1,size(this)
        output = output                                    &
             & + this%states(i)%harmonic_potential_energy( &
             &                           this%states(j),   &
             &                           subspace,         &
             &                           subspace_basis,   &
             &                           anharmonic_data ) &
             & * this%coefficients(i)                      &
             & * this%coefficients(j)
      enddo
    enddo
  endif
end function

impure elemental function kinetic_stress_PolynomialState(this,ket, &
   & subspace,subspace_basis,stress_prefactors,anharmonic_data) result(output)
  implicit none
  
  class(PolynomialState),   intent(in)           :: this
  class(BasisState),        intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(StressPrefactors),   intent(in)           :: stress_prefactors
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(RealMatrix)                               :: output
  
  type(PolynomialState) :: polynomial_ket
  
  integer :: i,j
  
  if (present(ket)) then
    polynomial_ket = PolynomialState(ket)
    output = dblemat(zeroes(3,3))
    do i=1,size(this)
      do j=1,size(polynomial_ket)
        output = output                                                     &
             & + this%states(i)%kinetic_stress( polynomial_ket%states(j),   &
             &                                  subspace,                   &
             &                                  subspace_basis,             &
             &                                  stress_prefactors,          &
             &                                  anharmonic_data           ) &
             & * this%coefficients(i)                                       &
             & * polynomial_ket%coefficients(j)
      enddo
    enddo
  else
    output = dblemat(zeroes(3,3))
    do i=1,size(this)
      do j=1,size(this)
        output = output                                              &
             & + this%states(i)%kinetic_stress( this%states(j),      &
             &                                  subspace,            &
             &                                  subspace_basis,      &
             &                                  stress_prefactors,   &
             &                                  anharmonic_data    ) &
             & * this%coefficients(i)                                &
             & * this%coefficients(j)
      enddo
    enddo
  endif
end function

! ----------------------------------------------------------------------
! Change the modes of the state by the specified group.
! ----------------------------------------------------------------------
impure elemental function change_modes_PolynomialState(this,mode_group) &
   & result(output)
  implicit none
  
  class(PolynomialState), intent(in) :: this
  type(Group),            intent(in) :: mode_group
  type(PolynomialState)              :: output
  
  output = this
  output%states = output%states%change_modes(mode_group)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_PolynomialState(this,input)
  implicit none
  
  class(PolynomialState), intent(out) :: this
  type(String),           intent(in)  :: input(:)
  
  integer                          :: subspace_id
  real(dp),            allocatable :: coefficients(:)
  type(MonomialState), allocatable :: states(:)
  
  type(StringArray), allocatable :: sections(:)
  type(String),      allocatable :: line(:)
  
  integer :: i,ialloc
  
  select type(this); type is(PolynomialState)
    line = split_line(input(1))
    subspace_id = int(line(3))
    
    sections = split_into_sections(input(3:))
    allocate( coefficients(size(sections)), &
            & states(size(sections)),       &
            & stat=ialloc); call err(ialloc)
    do i=1,size(sections)
      line = split_line(sections(i)%strings(1))
      coefficients(i) = dble(line(3))
      states(i) = MonomialState(sections(i)%strings(2:))
    enddo
    
    this = PolynomialState( subspace_id  = subspace_id,  &
                          & coefficients = coefficients, &
                          & states       = states        )
  class default
    call err()
  end select
end subroutine

function write_PolynomialState(this) result(output)
  implicit none
  
  class(PolynomialState), intent(in) :: this
  type(String), allocatable          :: output(:)
  
  type(StringArray), allocatable :: states(:)
  
  integer :: i
  
  select type(this); type is(PolynomialState)
    states = [( StringArray( [ 'Coefficient : '//this%coefficients(i),     &
              &                str(this%states(i))                    ] ), &
              & i=1,                                                       &
              & size(this)                                                 )]
    output = [ 'Subspace : '//this%subspace_id,      &
             & str(''),                              &
             & str(join(states, separating_line='')) ]
  class default
    call err()
  end select
end function

function new_PolynomialState_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(PolynomialState)    :: this
  
  call this%read(input)
end function

impure elemental function new_PolynomialState_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(PolynomialState)         :: this
  
  this = PolynomialState(str(input))
end function
end module
