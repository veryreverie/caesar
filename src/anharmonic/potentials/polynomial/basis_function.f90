! ======================================================================
! Basis functions for generating the Born-Oppenheimer surface mapping.
! ======================================================================
module basis_function_module
  use common_module
  
  use states_module
  use anharmonic_common_module
  implicit none
  
  private
  
  public :: BasisFunction
  public :: generate_basis_functions
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: operator(-)
  
  type, extends(Stringsable) :: BasisFunction
    ! The basis function in real co-ordinates.
    type(RealPolynomial), private :: real_representation_
    
    ! The basis function in complex co-ordinates.
    type(ComplexPolynomial), private :: complex_representation_
    
    ! A leading coefficient.
    real(dp), private :: coefficient_
  contains
    procedure, public :: real_representation
    procedure, public :: complex_representation
    
    procedure, public :: simplify => simplify_BasisFunction
    
    generic,   public  :: energy =>                                  &
                        & energy_RealModeDisplacement_BasisFunction, &
                        & energy_ComplexModeDisplacement_BasisFunction
    procedure, private :: energy_RealModeDisplacement_BasisFunction
    procedure, private :: energy_ComplexModeDisplacement_BasisFunction
    generic,   public  :: force =>                                  &
                        & force_RealModeDisplacement_BasisFunction, &
                        & force_ComplexModeDisplacement_BasisFunction
    procedure, private :: force_RealModeDisplacement_BasisFunction
    procedure, private :: force_ComplexModeDisplacement_BasisFunction
    
    generic,   public :: braket =>     &
                       & braket_state, &
                       & braket_states
    procedure, public :: braket_state  => braket_state_BasisFunction
    procedure, public :: braket_states => braket_states_BasisFunction
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_BasisFunction
    
    procedure, public :: undisplaced_energy => undisplaced_energy_BasisFunction
    
    procedure, private :: internal_coefficient
    procedure, public  :: coefficient => coefficient_BasisFunction
    procedure, public  :: set_coefficient => set_coefficient_BasisFunction
    
    procedure, public :: read  => read_BasisFunction
    procedure, public :: write => write_BasisFunction
  end type
  
  interface BasisFunction
    module procedure new_BasisFunction
    module procedure new_BasisFunction_Strings
    module procedure new_BasisFunction_StringArray
  end interface
  
  interface generate_basis_functions
    module procedure generate_basis_functions_SubspaceMonomial
  end interface
  
  interface operator(*)
    module procedure multiply_BasisFunction_real
    module procedure multiply_real_BasisFunction
  end interface
  
  interface operator(/)
    module procedure divide_BasisFunction_real
  end interface
  
  interface operator(+)
    module procedure add_BasisFunction_BasisFunction
  end interface
  
  interface operator(-)
    module procedure negative_BasisFunction
    
    module procedure subtract_BasisFunction_BasisFunction
  end interface
contains

! ----------------------------------------------------------------------
! Constructor.
! ----------------------------------------------------------------------
function new_BasisFunction(real_representation,complex_representation) &
   & result(this)
  implicit none
  
  type(RealPolynomial),    intent(in) :: real_representation
  type(ComplexPolynomial), intent(in) :: complex_representation
  type(BasisFunction)                 :: this
  
  this%real_representation_    = real_representation
  this%complex_representation_ = complex_representation
  this%coefficient_            = 1.0_dp
end function

! ----------------------------------------------------------------------
! Getters.
! ----------------------------------------------------------------------
impure elemental function real_representation(this) result(output)
  implicit none
  
  class(BasisFunction), intent(in) :: this
  type(RealPolynomial)             :: output
  
  output = this%coefficient_*this%real_representation_
end function

impure elemental function complex_representation(this) result(output)
  implicit none
  
  class(BasisFunction), intent(in) :: this
  type(ComplexPolynomial)          :: output
  
  output = this%coefficient_*this%complex_representation_
end function

! ----------------------------------------------------------------------
! Generate basis functions.
! ----------------------------------------------------------------------
function generate_basis_functions_SubspaceMonomial(subspace_monomial, &
   & structure,complex_modes,real_modes,qpoints,subspaces,            &
   & degenerate_symmetries,vscf_basis_functions_only,logfile) result(output)
  implicit none
  
  type(SubspaceMonomial),   intent(in)    :: subspace_monomial
  type(StructureData),      intent(in)    :: structure
  type(ComplexMode),        intent(in)    :: complex_modes(:)
  type(RealMode),           intent(in)    :: real_modes(:)
  type(QpointData),         intent(in)    :: qpoints(:)
  type(DegenerateSubspace), intent(in)    :: subspaces(:)
  type(DegenerateSymmetry), intent(in)    :: degenerate_symmetries(:)
  logical,                  intent(in)    :: vscf_basis_functions_only
  type(OFile),              intent(inout) :: logfile
  type(BasisFunction), allocatable        :: output(:)
  
  ! Monomials, in complex and real representations.
  type(ComplexMonomial), allocatable :: complex_monomials(:)
  type(RealMonomial),    allocatable :: real_monomials(:)
  
  ! The conversion from the complex monomial basis to the real monomial basis,
  !    and vice-versa.
  type(ComplexMatrix) :: complex_to_real_conversion
  type(ComplexMatrix) :: real_to_complex_conversion
  
  ! Symmetry data.
  type(ComplexMatrix) :: symmetry
  type(ComplexMatrix) :: projection
  
  ! Polynomial coefficients, in both bases.
  type(SymmetricEigenstuff), allocatable :: estuff(:)
  real(dp),                  allocatable :: real_coefficients(:)
  complex(dp),               allocatable :: complex_coefficients(:)
  
  ! Variables for constructing the output.
  type(RealPolynomial)    :: real_representation
  type(ComplexPolynomial) :: complex_representation
  
  ! Temporary variables.
  integer :: i,ialloc
  
  if (sum(subspace_monomial%powers)<2) then
    call print_line(CODE_ERROR//': Trying to generate basis functions with &
       &power less than 2.')
    call err()
  endif
  
  ! Generate the real monomials and complex monomials corresponding to the
  !    subspace monomial, with coefficients such that symmetries are unitary.
  real_monomials = generate_real_monomials( subspace_monomial, &
                                          & subspaces,         &
                                          & real_modes,        &
                                          & qpoints            )
  
  complex_monomials = generate_complex_monomials(            &
      & subspace_monomial,                                   &
      & subspaces,                                           &
      & complex_modes,                                       &
      & qpoints,                                             &
      & conserve_momentum=.true.,                            &
      & conserve_subspace_momentum=vscf_basis_functions_only )
  
  if (size(complex_monomials)==0) then
    output = [BasisFunction::]
    return
  endif
  
  ! Identify the mappings between complex monomials and real monomials.
  complex_to_real_conversion = basis_conversion_matrix( &
                          & real_monomials,             &
                          & complex_monomials,          &
                          & include_coefficients=.true. )
  
  real_to_complex_conversion = basis_conversion_matrix( &
                          & complex_monomials,          &
                          & real_monomials,             &
                          & include_coefficients=.true. )
  
  ! Construct projection matrix, which has allowed basis functions as
  !    eigenvectors with eigenvalue 1, and sends all other functions to 0.
  projection = cmplxmat(make_identity_matrix(size(complex_monomials)))
  do i=1,size(degenerate_symmetries)
    ! Constuct symmetry in complex monomial co-ordinates.
    symmetry = degenerate_symmetries(i)%calculate_symmetry( &
                              & complex_monomials,          &
                              & complex_modes,              &
                              & include_coefficients=.true. )
    call check_unitary(symmetry, 'symmetry in monomial basis', logfile)
    
    ! Construct the projection matrix for this symmetry,
    !    and multiply the total projection matrix by this.
    projection = projection                                                  &
             & * projection_matrix( symmetry,                                &
             &                      structure%symmetries(i)%symmetry_order() )
  enddo
  call check_hermitian( projection,               &
                      & 'projection matrix',      &
                      & logfile,                  &
                      & ignore_threshold=1e-10_dp )
  
  ! Transform the projection matrix into real co-ordinates,
  !    and check that it is real and symmetric.
  projection = complex_to_real_conversion &
           & * projection                 &
           & * real_to_complex_conversion
  call check_real(projection,'projection matrix',logfile)
  call check_symmetric( real(projection),         &
                      & 'projection_matrix',      &
                      & logfile,                  &
                      & ignore_threshold=1e-10_dp )
  
  ! Diagonalise the projection matrix,
  !    check its eigenvalues are either 0 or 1,
  !    and select only the eigenvectors with eigenvalue 1.
  estuff = diagonalise_symmetric(real(projection))
  if (any(abs(estuff%eval-1)>1e-2_dp .and. abs(estuff%eval)>1e-2_dp)) then
    call print_line(ERROR//': Projection matrix has eigenvalues which are &
       &neither 0 nor 1.')
    call print_line('Eigenvalues:')
    call print_line(estuff%eval)
    call err()
  endif
  estuff = estuff(filter(abs(estuff%eval-1)<1e-2_dp))
  
  ! Construct basis functions from the coefficients.
  allocate(output(size(estuff)), stat=ialloc); call err(ialloc)
  do i=1,size(estuff)
    real_coefficients = estuff(i)%evec
    complex_coefficients = cmplx( real_to_complex_conversion &
                              & * vec(real_coefficients)     )
    
    real_representation = RealPolynomial( real_coefficients &
                                      & * real_monomials    )
    complex_representation = ComplexPolynomial( complex_coefficients &
                                            & * complex_monomials    )
    
    output(i) = BasisFunction( real_representation,   &
                             & complex_representation )
    call output(i)%simplify()
  enddo
end function

! Given a unitary matrix U, s.t. U^n=I, returns the matrix
!    H = sum(j=0,n-1) U^j / n
! H is Hermitian, and is a projection matrix which projects out the
!    eigenvectors of U with eigenvalue 1.
! Uses Horner's rule to calculate the sum with minimal matrix multiplication.
function projection_matrix(input,order) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: input
  integer,             intent(in) :: order
  type(ComplexMatrix)             :: output
  
  type(ComplexMatrix) :: identity
  
  integer :: i
  
  if (order<1) then
    call print_line(CODE_ERROR//': symmetry order may not be < 1.')
    call err()
  elseif (order>6) then
    call print_line(CODE_ERROR//': symmetry order may not be > 6.')
    call err()
  endif
  
  identity = cmplxmat(make_identity_matrix(size(input,1)))
  
  output = identity
  do i=2,order
    output = input*output + identity
  enddo
  output = output/order
end function

! ----------------------------------------------------------------------
! Simplify the basis function.
! ----------------------------------------------------------------------
impure elemental subroutine simplify_BasisFunction(this)
  implicit none
  
  class(BasisFunction), intent(inout) :: this
  
  call this%real_representation_%simplify()
  call this%complex_representation_%simplify()
end subroutine

! ----------------------------------------------------------------------
! Evaluate the energy and forces due to the basis function.
! ----------------------------------------------------------------------
impure elemental function energy_RealModeDisplacement_BasisFunction(this, &
   & displacement) result(output)
  implicit none
  
  class(BasisFunction),        intent(in) :: this
  class(RealModeDisplacement), intent(in) :: displacement
  real(dp)                                :: output
  
  output = this%coefficient_ &
       & * this%real_representation_%energy(displacement)
end function

impure elemental function energy_ComplexModeDisplacement_BasisFunction(this, &
   & displacement) result(output)
  implicit none
  
  class(BasisFunction),           intent(in) :: this
  class(ComplexModeDisplacement), intent(in) :: displacement
  complex(dp)                                :: output
  
  output = this%coefficient_ &
       & * this%complex_representation_%energy(displacement)
end function

impure elemental function force_RealModeDisplacement_BasisFunction(this, &
   & displacement) result(output)
  implicit none
  
  class(BasisFunction),        intent(in) :: this
  class(RealModeDisplacement), intent(in) :: displacement
  type(RealModeForce)                     :: output
  
  output = this%coefficient_ &
       & * this%real_representation_%force(displacement)
end function

impure elemental function force_ComplexModeDisplacement_BasisFunction(this, &
   & displacement) result(output)
  implicit none
  
  class(BasisFunction),           intent(in) :: this
  class(ComplexModeDisplacement), intent(in) :: displacement
  type(ComplexModeForce)                     :: output
  
  output = this%coefficient_ &
       & * this%complex_representation_%force(displacement)
end function

! ----------------------------------------------------------------------
! Integrate the basis function between two states.
! ----------------------------------------------------------------------
impure elemental subroutine braket_state_BasisFunction(this,bra,ket,subspace, &
   & subspace_basis,anharmonic_data,qpoint)
  implicit none
  
  class(BasisFunction),     intent(inout)        :: this
  class(SubspaceState),     intent(in)           :: bra
  class(SubspaceState),     intent(in), optional :: ket
  type(DegenerateSubspace), intent(in)           :: subspace
  class(SubspaceBasis),     intent(in)           :: subspace_basis
  type(AnharmonicData),     intent(in)           :: anharmonic_data
  type(QpointData),         intent(in), optional :: qpoint
  
  type(ComplexMode),    allocatable :: subspace_modes(:)
  type(QpointData),     allocatable :: subspace_qpoints(:)
  integer,              allocatable :: integrated_modes(:)
  
  type(ComplexMatrix)               :: complex_to_real_conversion
  complex(dp),          allocatable :: complex_coefficients(:)
  real(dp),             allocatable :: real_coefficients(:)
  type(RealUnivariate), allocatable :: real_modes(:)
  logical,              allocatable :: mode_integrated(:)
  
  integer :: i,j
  
  subspace_modes = subspace%modes(anharmonic_data%complex_modes)
  subspace_qpoints = subspace%qpoints( anharmonic_data%complex_modes, &
                                     & anharmonic_data%qpoints        )
  if (present(qpoint)) then
    integrated_modes = subspace_modes(filter(         &
       & subspace_qpoints%id==qpoint%id .or.          &
       & subspace_qpoints%id==qpoint%paired_qpoint_id ))%id
  else
    integrated_modes = subspace_modes%id
  endif
  
  ! Generate conversion between complex and real representation.
  complex_to_real_conversion = coefficient_conversion_matrix( &
                        & this%real_representation_%terms,    &
                        & this%complex_representation_%terms, &
                        & include_coefficients = .false.      )
  
  ! Perform integration in complex co-ordinates.
  this%complex_representation_ = bra%braket( this%complex_representation_, &
                                           & ket,                          &
                                           & subspace,                     &
                                           & subspace_basis,               &
                                           & anharmonic_data,              &
                                           & qpoint                        )
  
  ! Use calculated complex coefficients and conversion to generate new
  !    coefficients for real representation.
  complex_coefficients = this%complex_representation_%terms%coefficient
  real_coefficients = real(cmplx( complex_to_real_conversion &
                              & * vec(complex_coefficients)  ))
  this%real_representation_%terms%coefficient = real_coefficients
  
  ! Remove modes in real representation which have been integrated over.
  do i=1,size(this%real_representation_)
    real_modes = this%real_representation_%terms(i)%modes()
    mode_integrated = [( any(real_modes(j)%id==integrated_modes), &
                       & j=1,                                     &
                       & size(real_modes)                         )]
    real_modes = real_modes(filter(.not.mode_integrated))
    this%real_representation_%terms(i) = RealMonomial(                &
       & modes       = real_modes,                                    &
       & coefficient = this%real_representation_%terms(i)%coefficient )
  enddo
end subroutine

impure elemental subroutine braket_states_BasisFunction(this,states,subspace, &
   & subspace_basis,anharmonic_data)
  implicit none
  
  class(BasisFunction),     intent(inout) :: this
  class(SubspaceStates),    intent(in)    :: states
  type(DegenerateSubspace), intent(in)    :: subspace
  class(SubspaceBasis),     intent(in)    :: subspace_basis
  type(AnharmonicData),     intent(in)    :: anharmonic_data
  
  type(ComplexMode), allocatable :: subspace_modes(:)
  integer,           allocatable :: integrated_modes(:)
  
  type(ComplexMatrix)               :: complex_to_real_conversion
  complex(dp),          allocatable :: complex_coefficients(:)
  real(dp),             allocatable :: real_coefficients(:)
  type(RealUnivariate), allocatable :: real_modes(:)
  logical,              allocatable :: mode_integrated(:)
  
  integer :: i,j
  
  subspace_modes = subspace%modes(anharmonic_data%complex_modes)
  integrated_modes = subspace_modes%id
  
  ! Generate conversion between complex and real representation.
  complex_to_real_conversion = coefficient_conversion_matrix( &
                        & this%real_representation_%terms,    &
                        & this%complex_representation_%terms, &
                        & include_coefficients = .false.      )
  
  ! Perform integration in complex co-ordinates.
  this%complex_representation_ = states%braket( this%complex_representation_, &
                                              & subspace,                     &
                                              & subspace_basis,               &
                                              & anharmonic_data               )
  
  ! Use calculated complex coefficients and conversion to generate new
  !    coefficients for real representation.
  complex_coefficients = this%complex_representation_%terms%coefficient
  real_coefficients = real(cmplx( complex_to_real_conversion &
                              & * vec(complex_coefficients)  ))
  this%real_representation_%terms%coefficient = real_coefficients
  
  ! Remove modes in real representation which have been integrated over.
  do i=1,size(this%real_representation_)
    real_modes = this%real_representation_%terms(i)%modes()
    mode_integrated = [( any(real_modes(j)%id==integrated_modes), &
                       & j=1,                                     &
                       & size(real_modes)                         )]
    real_modes = real_modes(filter(.not.mode_integrated))
    this%real_representation_%terms(i) = RealMonomial(                &
       & modes       = real_modes,                                    &
       & coefficient = this%real_representation_%terms(i)%coefficient )
  enddo
end subroutine

! ----------------------------------------------------------------------
! Returns the thermal expectation of the basis function.
! ----------------------------------------------------------------------
impure elemental function harmonic_expectation_BasisFunction(this,frequency, &
   & thermal_energy,no_states,subspace,anharmonic_data) result(output)
  implicit none
  
  class(BasisFunction),     intent(in) :: this
  real(dp),                 intent(in) :: frequency
  real(dp),                 intent(in) :: thermal_energy
  integer,                  intent(in) :: no_states
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  real(dp)                             :: output
  
  integer :: i
  
  output = 0
  do i=1,size(this%complex_representation_)
    output = output                                                       &
         & + harmonic_expectation( this%complex_representation_%terms(i), &
         &                         frequency,                             &
         &                         thermal_energy,                        &
         &                         no_states,                             &
         &                         subspace,                              &
         &                         anharmonic_data%anharmonic_supercell,  &
         &                         anharmonic_data                        )
  enddo
  output = output * this%coefficient_
end function

! ----------------------------------------------------------------------
! Returns the energy at zero displacement.
! ----------------------------------------------------------------------
impure elemental function undisplaced_energy_BasisFunction(this) result(output)
  implicit none
  
  class(BasisFunction), intent(in) :: this
  real(dp)                         :: output
  
  type(RealModeDisplacement) :: zero_displacement
  
  zero_displacement = RealModeDisplacement([RealSingleDisplacement::])
  
  output = this%energy(zero_displacement)
end function

! ----------------------------------------------------------------------
! Converts the basis function to and from a single coefficient.
! The coefficient is the sum of the absolute values of the
!    real representation coefficients, with the sign of the first non-zero
!    real representation coefficient.
! These methods allow for simpler linear algebra with basis functions.
! ----------------------------------------------------------------------
impure elemental function coefficient_BasisFunction(this) result(output)
  implicit none
  
  class(BasisFunction), intent(in) :: this
  real(dp)                         :: output
  
  output = this%coefficient_ * this%internal_coefficient()
end function

impure elemental subroutine set_coefficient_BasisFunction(this,coefficient)
  implicit none
  
  class(BasisFunction), intent(inout) :: this
  real(dp),             intent(in)    :: coefficient
  
  real(dp) :: internal_coefficient
  
  internal_coefficient = this%internal_coefficient()
  
  if (abs(internal_coefficient)>1.0e-300_dp) then
    this%coefficient_ = coefficient / this%internal_coefficient()
  endif
end subroutine

function internal_coefficient(this) result(output)
  implicit none
  
  class(BasisFunction), intent(in) :: this
  real(dp)                         :: output
  
  integer ::  i
  
  output = sum(abs(this%real_representation_%terms%coefficient))
  
  i = first(abs(this%real_representation_%terms%coefficient)>0, default=0)
  if (i/=0) then
    if (this%real_representation_%terms(i)%coefficient<0) then
      output = -output
    endif
  endif
end function

! ----------------------------------------------------------------------
! Arithmetic.
! ----------------------------------------------------------------------
impure elemental function multiply_BasisFunction_real(this,that) result(output)
  implicit none
  
  type(BasisFunction), intent(in) :: this
  real(dp),            intent(in) :: that
  type(BasisFunction)             :: output
  
  output = this
  output%coefficient_ = output%coefficient_ * that
end function

impure elemental function multiply_real_BasisFunction(this,that) result(output)
  implicit none
  
  real(dp),            intent(in) :: this
  type(BasisFunction), intent(in) :: that
  type(BasisFunction)             :: output
  
  output = that
  output%coefficient_ = this * output%coefficient_
end function

impure elemental function divide_BasisFunction_real(this,that) result(output)
  implicit none
  
  type(BasisFunction), intent(in) :: this
  real(dp),            intent(in) :: that
  type(BasisFunction)             :: output
  
  output = this
  output%coefficient_ = output%coefficient_ / that
end function

impure elemental function add_BasisFunction_BasisFunction(this,that) &
   & result(output)
  implicit none
  
  type(BasisFunction), intent(in) :: this
  type(BasisFunction), intent(in) :: that
  type(BasisFunction)             :: output
  
  output = BasisFunction(                                          &
     & this%real_representation()+that%real_representation(),      &
     & this%complex_representation()+that%complex_representation() )
end function

impure elemental function negative_BasisFunction(this) result(output)
  implicit none
  
  type(BasisFunction), intent(in) :: this
  type(BasisFunction)             :: output
  
  output = this
  output%coefficient_ = -output%coefficient_
end function

impure elemental function subtract_BasisFunction_BasisFunction(this,that) &
   & result(output)
  implicit none
  
  type(BasisFunction), intent(in) :: this
  type(BasisFunction), intent(in) :: that
  type(BasisFunction)             :: output
  
  output = BasisFunction(                                          &
     & this%real_representation()-that%real_representation(),      &
     & this%complex_representation()-that%complex_representation() )
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_BasisFunction(this,input)
  implicit none
  
  class(BasisFunction), intent(out) :: this
  type(String),         intent(in)  :: input(:)
  
  integer :: partition_line
  
  type(RealPolynomial)    :: real_representation
  type(ComplexPolynomial) :: complex_representation
  
  integer :: i
  
  select type(this); type is(BasisFunction)
    ! Locate the line between real terms and complex terms.
    do i=2,size(input)
      if (size(split_line(input(i)))>1) then
        partition_line = i
        exit
      endif
    enddo
    
    real_representation = RealPolynomial(                 &
       & join(input(2:partition_line-1), delimiter=' + ') )
    complex_representation = ComplexPolynomial(          &
       & join(input(partition_line+1:), delimiter=' + ') )
    
    this = BasisFunction( real_representation    = real_representation,   &
                        & complex_representation = complex_representation )
  class default
    call err()
  end select
end subroutine

function write_BasisFunction(this) result(output)
  implicit none
  
  class(BasisFunction), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  type(RealPolynomial)    :: real_representation
  type(ComplexPolynomial) :: complex_representation
  
  select type(this); type is(BasisFunction)
    real_representation = this%real_representation()
    complex_representation = this%complex_representation()
    output = [ str('Basis function in real co-ordinates:'),    &
             & str(real_representation%terms),                 &
             & str('Basis function in complex co-ordinates:'), &
             & str(complex_representation%terms)               ]
  class default
    call err()
  end select
end function

function new_BasisFunction_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(BasisFunction)      :: this
  
  call this%read(input)
end function

impure elemental function new_BasisFunction_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(BasisFunction)           :: this
  
  this = BasisFunction(str(input))
end function
end module
