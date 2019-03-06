! ======================================================================
! A sum of monomial states and coefficients.
! ======================================================================
module polynomial_state_module
  use common_module
  
  use anharmonic_common_module
  
  use monomial_state_module
  implicit none
  
  private
  
  public :: startup_polynomial_state
  
  public :: PolynomialState
  
  public :: size
  
  type, extends(SubspaceState) :: PolynomialState
    type(MonomialState), allocatable :: states(:)
    real(dp),            allocatable :: coefficients(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_PolynomialState
    
    procedure, public :: braket_SubspaceState => &
                       & braket_SubspaceState_PolynomialState
    procedure, public :: braket_ComplexUnivariate => &
                       & braket_ComplexUnivariate_PolynomialState
    procedure, public :: braket_ComplexMonomial => &
                       & braket_ComplexMonomial_PolynomialState
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_PolynomialState
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_PolynomialState
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_PolynomialState
    
    ! I/O.
    procedure, public :: read  => read_PolynomialState
    procedure, public :: write => write_PolynomialState
  end type
  
  interface PolynomialState
    module procedure new_PolynomialState
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
! SubspaceState methods.
! ----------------------------------------------------------------------
impure elemental function braket_SubspaceState_PolynomialState(this, &
   & ket,subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(PolynomialState),   intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  integer :: i,j
  
  if (present(ket)) then
    select type(ket); type is(PolynomialState)
      output = 0.0_dp
      do i=1,size(this)
        do j=1,size(ket)
          output = output                    &
               & + braket( this%states(i),   &
               &           ket%states(j),    &
               &           subspace,         &
               &           subspace_basis,   &
               &           anharmonic_data ) &
               & * this%coefficients(i)      &
               & * ket%coefficients(j)
        enddo
      enddo
    class default
      call err()
    end select
  else
    output = 0.0_dp
    do i=1,size(this)
      do j=1,size(this)
        output = output                    &
             & + braket( this%states(i),   &
             &           this%states(j),   &
             &           subspace,         &
             &           subspace_basis,   &
             &           anharmonic_data ) &
             & * this%coefficients(i)      &
             & * this%coefficients(j)
      enddo
    enddo
  endif
end function

impure elemental function braket_ComplexUnivariate_PolynomialState(this, &
   & univariate,ket,subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(PolynomialState),   intent(in)           :: this
  type(ComplexUnivariate),  intent(in)           :: univariate
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(ComplexMonomial)                          :: output
  
  type(ComplexMonomial) :: integrated_univariate
  
  integer :: i,j
  
  if (present(ket)) then
    select type(ket); type is(PolynomialState)
      do i=1,size(this)
        do j=1,size(ket)
          integrated_univariate = braket( this%states(i),   &
                             &            univariate,       &
                             &            ket%states(j),    &
                             &            subspace,         &
                             &            subspace_basis,   &
                             &            anharmonic_data ) &
                             & * this%coefficients(i)       &
                             & * ket%coefficients(j)
          if (i==1 .and. j==1) then
            output = integrated_univariate
          else
            output%coefficient = output%coefficient &
                             & + integrated_univariate%coefficient
          endif
        enddo
      enddo
    class default
      call err()
    end select
  else
    do i=1,size(this)
      do j=1,size(this)
        integrated_univariate = braket( this%states(i),   &
                           &            univariate,       &
                           &            this%states(j),   &
                           &            subspace,         &
                           &            subspace_basis,   &
                           &            anharmonic_data ) &
                           & * this%coefficients(i)       &
                           & * this%coefficients(j)
        if (i==1 .and. j==1) then
          output = integrated_univariate
        else
          output%coefficient = output%coefficient &
                           & + integrated_univariate%coefficient
        endif
      enddo
    enddo
  endif
end function

impure elemental function braket_ComplexMonomial_PolynomialState(this, &
   & monomial,ket,subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(PolynomialState),   intent(in)           :: this
  type(ComplexMonomial),    intent(in)           :: monomial
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(ComplexMonomial)                          :: output
  
  type(ComplexMonomial) :: integrated_monomial
  
  integer :: i,j
  
  if (present(ket)) then
    select type(ket); type is(PolynomialState)
      do i=1,size(this)
        do j=1,size(ket)
          integrated_monomial = braket( this%states(i),   &
                            &           monomial,         &
                            &           ket%states(j),    &
                            &           subspace,         &
                            &           subspace_basis,   &
                            &           anharmonic_data ) &
                            & * this%coefficients(i)      &
                            & * ket%coefficients(j)
          if (i==1 .and. j==1) then
            output = integrated_monomial
          else
            output%coefficient = output%coefficient &
                             & + integrated_monomial%coefficient
          endif
        enddo
      enddo
    class default
      call err()
    end select
  else
    do i=1,size(this)
      do j=1,size(this)
        integrated_monomial = braket( this%states(i),   &
                          &           monomial,         &
                          &           this%states(j),   &
                          &           subspace,         &
                          &           subspace_basis,   &
                          &           anharmonic_data ) &
                          & * this%coefficients(i)      &
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
   & subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(PolynomialState),   intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  integer :: i,j
  
  if (present(ket)) then
    select type(ket); type is(PolynomialState)
      output = 0.0_dp
      do i=1,size(this)
        do j=1,size(ket)
          output = output                            &
               & + kinetic_energy( this%states(i),   &
               &                   ket%states(j),    &
               &                   subspace,         &
               &                   subspace_basis,   &
               &                   anharmonic_data ) &
               & * this%coefficients(i)              &
               & * ket%coefficients(j)
        enddo
      enddo
    class default
      call err()
    end select
  else
    output = 0.0_dp
    do i=1,size(this)
      do j=1,size(this)
        output = output                            &
             & + kinetic_energy( this%states(i),   &
             &                   this%states(j),   &
             &                   subspace,         &
             &                   subspace_basis,   &
             &                   anharmonic_data ) &
             & * this%coefficients(i)              &
             & * this%coefficients(j)
      enddo
    enddo
  endif
end function

impure elemental function harmonic_potential_energy_PolynomialState( &
   & this,ket,subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(PolynomialState),   intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  integer :: i,j
  
  
  if (present(ket)) then
    select type(ket); type is(PolynomialState)
      output = 0.0_dp
      do i=1,size(this)
        do j=1,size(ket)
          output = output                                       &
               & + harmonic_potential_energy( this%states(i),   &
               &                              ket%states(j),    &
               &                              subspace,         &
               &                              subspace_basis,   &
               &                              anharmonic_data ) &
               & * this%coefficients(i)                         &
               & * ket%coefficients(j)
        enddo
      enddo
    class default
      call err()
    end select
  else
    output = 0.0_dp
    do i=1,size(this)
      do j=1,size(this)
        output = output                                       &
             & + harmonic_potential_energy( this%states(i),   &
             &                              this%states(j),   &
             &                              subspace,         &
             &                              subspace_basis,   &
             &                              anharmonic_data ) &
             & * this%coefficients(i)                         &
             & * this%coefficients(j)
      enddo
    enddo
  endif
end function

impure elemental function kinetic_stress_PolynomialState(this,ket, &
   & subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(PolynomialState),   intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(RealMatrix)                               :: output
  
  integer :: i,j
  
  if (present(ket)) then
    select type(ket); type is(PolynomialState)
      output = dblemat(zeroes(3,3))
      do i=1,size(this)
        do j=1,size(ket)
          output = output                            &
               & + kinetic_stress( this%states(i),   &
               &                   ket%states(j),    &
               &                   subspace,         &
               &                   subspace_basis,   &
               &                   anharmonic_data ) &
               & * this%coefficients(i)              &
               & * ket%coefficients(j)
        enddo
      enddo
    class default
      call err()
    end select
  else
    output = dblemat(zeroes(3,3))
    do i=1,size(this)
      do j=1,size(this)
        output = output                            &
             & + kinetic_stress( this%states(i),   &
             &                   this%states(j),   &
             &                   subspace,         &
             &                   subspace_basis,   &
             &                   anharmonic_data ) &
             & * this%coefficients(i)              &
             & * this%coefficients(j)
      enddo
    enddo
  endif
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
