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
  public :: braket_PolynomialState
  public :: kinetic_energy_PolynomialState
  public :: harmonic_potential_energy_PolynomialState
  
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
                       & kinetic_energy_PolynomialState2
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_PolynomialState2
    
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
  
  interface braket_PolynomialState
    module procedure braket_PolynomialStates
    module procedure braket_PolynomialStates_ComplexUnivariate
    module procedure braket_PolynomialStates_ComplexMonomial
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
! Integrals of the form <bra|ket> and <bra|potential|ket>.
! ----------------------------------------------------------------------
impure elemental function braket_PolynomialStates(bra,ket) result(output)
  implicit none
  
  type(PolynomialState), intent(in) :: bra
  type(PolynomialState), intent(in) :: ket
  real(dp)                          :: output
  
  integer :: i,j
  
  output = 0.0_dp
  do i=1,size(bra)
    do j=1,size(ket)
      output = output                                            &
           & + braket_MonomialState(bra%states(i),ket%states(j)) &
           & * bra%coefficients(i)                               &
           & * ket%coefficients(j)
    enddo
  enddo
end function

impure elemental function braket_PolynomialStates_ComplexUnivariate(bra,ket, &
   & univariate,subspace,supercell) result(output)
  implicit none
  
  type(PolynomialState),    intent(in) :: bra
  type(PolynomialState),    intent(in) :: ket
  type(ComplexUnivariate),  intent(in) :: univariate
  type(DegenerateSubspace), intent(in) :: subspace
  type(StructureData),      intent(in) :: supercell
  type(ComplexMonomial)                :: output
  
  type(ComplexMonomial) :: integrated_univariate
  
  integer :: i,j
  
  do i=1,size(bra)
    do j=1,size(ket)
      integrated_univariate = braket_MonomialState( bra%states(i),   &
                         &                          ket%states(j),   &
                         &                          univariate,      &
                         &                          subspace,        &
                         &                          supercell      ) &
                         & * bra%coefficients(i)                     &
                         & * ket%coefficients(j)
      if (i==1 .and. j==1) then
        output = integrated_univariate
      else
        output%coefficient = output%coefficient &
                         & + integrated_univariate%coefficient
      endif
    enddo
  enddo
end function

impure elemental function braket_PolynomialStates_ComplexMonomial(bra,ket, &
   & monomial,subspace,supercell) result(output)
  implicit none
  
  type(PolynomialState),    intent(in) :: bra
  type(PolynomialState),    intent(in) :: ket
  type(ComplexMonomial),    intent(in) :: monomial
  type(DegenerateSubspace), intent(in) :: subspace
  type(StructureData),      intent(in) :: supercell
  type(ComplexMonomial)                :: output
  
  type(ComplexMonomial) :: integrated_monomial
  
  integer :: i,j
  
  do i=1,size(bra)
    do j=1,size(ket)
      integrated_monomial = braket_MonomialState( bra%states(i),   &
                        &                         ket%states(j),   &
                        &                         monomial,        &
                        &                         subspace,        &
                        &                         supercell      ) &
                        & * bra%coefficients(i)                    &
                        & * ket%coefficients(j)
      if (i==1 .and. j==1) then
        output = integrated_monomial
      else
        output%coefficient = output%coefficient &
                         & + integrated_monomial%coefficient
      endif
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Evaluates <bra|T|ket>, where T is the kinetic energy operator.
! Gives the result per primitive cell.
! ----------------------------------------------------------------------
function kinetic_energy_PolynomialState(bra,ket,subspace,supercell) &
   & result(output)
  implicit none
  
  type(PolynomialState),    intent(in) :: bra
  type(PolynomialState),    intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  type(StructureData),      intent(in) :: supercell
  real(dp)                             :: output
  
  integer :: i,j
  
  output = 0.0_dp
  do i=1,size(bra)
    do j=1,size(ket)
      output = output                                         &
           & + kinetic_energy_MonomialState( bra%states(i),   &
           &                                 ket%states(j),   &
           &                                 subspace,        &
           &                                 supercell      ) &
           & * bra%coefficients(i)                            &
           & * ket%coefficients(j)
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Evaluates <bra|V|ket>, where V is the harmonic potential energy operator.
! Gives the result per primitive cell.
! ----------------------------------------------------------------------
function harmonic_potential_energy_PolynomialState(bra,ket,subspace, &
   & supercell) result(output)
  implicit none
  
  type(PolynomialState),    intent(in) :: bra
  type(PolynomialState),    intent(in) :: ket
  type(DegenerateSubspace), intent(in) :: subspace
  type(StructureData),      intent(in) :: supercell
  real(dp)                             :: output
  
  integer :: i,j
  
  output = 0.0_dp
  do i=1,size(bra)
    do j=1,size(ket)
      output = output                                                    &
           & + harmonic_potential_energy_MonomialState( bra%states(i),   &
           &                                            ket%states(j),   &
           &                                            subspace,        &
           &                                            supercell      ) &
           & * bra%coefficients(i)                                       &
           & * ket%coefficients(j)
    enddo
  enddo
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
  
  if (present(ket)) then
    select type(ket); type is(PolynomialState)
      output = braket_PolynomialState(this,ket)
    class default
      call err()
    end select
  else
    output = braket_PolynomialState(this,this)
  endif
end function

impure elemental function braket_ComplexUnivariate_PolynomialState( &
   & this,univariate,ket,subspace,subspace_basis,anharmonic_data)        &
   & result(output)
  implicit none
  
  class(PolynomialState),   intent(in)           :: this
  type(ComplexUnivariate),  intent(in)           :: univariate
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(ComplexMonomial)                          :: output
  
  if (present(ket)) then
    select type(ket); type is(PolynomialState)
      output = braket_PolynomialState( this,                                &
                                     & ket,                                 &
                                     & univariate,                          &
                                     & subspace,                            &
                                     & anharmonic_data%anharmonic_supercell )
    class default
      call err()
    end select
  else
    output = braket_PolynomialState( this,                                &
                                   & this,                                &
                                   & univariate,                          &
                                   & subspace,                            &
                                   & anharmonic_data%anharmonic_supercell )
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
  
  if (present(ket)) then
    select type(ket); type is(PolynomialState)
      output = braket_PolynomialState( this,                                &
                                     & ket,                                 &
                                     & monomial,                            &
                                     & subspace,                            &
                                     & anharmonic_data%anharmonic_supercell )
    class default
      call err()
    end select
  else
    output = braket_PolynomialState( this,                                &
                                   & this,                                &
                                   & monomial,                            &
                                   & subspace,                            &
                                   & anharmonic_data%anharmonic_supercell )
  endif
end function

impure elemental function kinetic_energy_PolynomialState2(this,ket, &
   & subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(PolynomialState),   intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  if (present(ket)) then
    select type(ket); type is(PolynomialState)
      output = kinetic_energy_PolynomialState(  &
         & this,                                &
         & ket,                                 &
         & subspace,                            &
         & anharmonic_data%anharmonic_supercell )
    class default
      call err()
    end select
  else
    output = kinetic_energy_PolynomialState(  &
       & this,                                &
       & this,                                &
       & subspace,                            &
       & anharmonic_data%anharmonic_supercell )
  endif
end function

impure elemental function harmonic_potential_energy_PolynomialState2( &
   & this,ket,subspace,subspace_basis,anharmonic_data) result(output)
  implicit none
  
  class(PolynomialState),   intent(in)           :: this
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  real(dp)                                       :: output
  
  if (present(ket)) then
    select type(ket); type is(PolynomialState)
      output = harmonic_potential_energy_PolynomialState( &
                   & this,                                &
                   & ket,                                 &
                   & subspace,                            &
                   & anharmonic_data%anharmonic_supercell )
    class default
      call err()
    end select
  else
    output = harmonic_potential_energy_PolynomialState( &
                 & this,                                &
                 & this,                                &
                 & subspace,                            &
                 & anharmonic_data%anharmonic_supercell )
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
