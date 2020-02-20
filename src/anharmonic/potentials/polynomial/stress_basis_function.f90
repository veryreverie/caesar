! ======================================================================
! Basis functions for generating the stress tensor mapping.
! ======================================================================
module stress_basis_function_module
  use common_module
  
  use anharmonic_common_module
  
  use polynomial_interpolator_module
  use basis_function_module
  use polynomial_dynamical_matrices_module
  implicit none
  
  private
  
  public :: StressBasisFunction
  public :: generate_stress_basis_functions
  public :: operator(*)
  public :: operator(/)
  
  type, extends(Stringsable) :: StressBasisFunction
    ! This will always be 3x3, but is allocatable to avoid stack overflows.
    type(ComplexPolynomial), private, allocatable :: elements_(:,:)
    
    ! A leading coefficient.
    real(dp), private :: coefficient_
  contains
    procedure, public :: simplify => simplify_StressBasisFunction
    
    generic,   public  :: stress =>                                        &
                        & stress_RealModeDisplacement_StressBasisFunction, &
                        & stress_ComplexModeDisplacement_StressBasisFunction
    procedure, private :: stress_RealModeDisplacement_StressBasisFunction
    procedure, private :: stress_ComplexModeDisplacement_StressBasisFunction
    
    generic,   public :: braket =>              &
                       & braket_SubspaceBraKet, &
                       & braket_BasisState,     &
                       & braket_BasisStates
    procedure, public :: braket_SubspaceBraKet => &
                       & braket_SubspaceBraKet_StressBasisFunction
    procedure, public :: braket_BasisState => &
                       & braket_BasisState_StressBasisFunction
    procedure, public :: braket_BasisStates => &
                       & braket_BasisStates_StressBasisFunction
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_StressBasisFunction
    
    procedure, public :: undisplaced_stress => &
                       & undisplaced_stress_StressBasisFunction
    
    procedure, public :: interpolate_coefficients => &
                       & interpolate_coefficients_StressBasisFunction
    
    procedure, public :: calculate_dynamical_matrices => &
                       & calculate_dynamical_matrices_StressBasisFunction
    
    ! I/O.
    procedure, public :: read  => read_StressBasisFunction
    procedure, public :: write => write_StressBasisFunction
  end type
  
  interface StressBasisFunction
    module procedure new_StressBasisFunction
    module procedure new_StressBasisFunction_Strings
    module procedure new_StressBasisFunction_StringArray
  end interface
  
  interface generate_stress_basis_functions
    module procedure generate_stress_basis_functions_SubspaceMonomial
  end interface
  
  interface operator(*)
    module procedure multiply_StressBasisFunction_real
    module procedure multiply_real_StressBasisFunction
  end interface
  
  interface operator(/)
    module procedure divide_StressBasisFunction_real
  end interface
contains

! Constructor.
function new_StressBasisFunction(elements,coefficient) result(this)
  implicit none
  
  type(ComplexPolynomial), intent(in)           :: elements(3,3)
  real(dp),                intent(in), optional :: coefficient
  type(StressBasisFunction)                     :: this
  
  this%elements_ = elements
  if (present(coefficient)) then
    this%coefficient_ = coefficient
  else
    this%coefficient_ = 1
  endif
end function

! Generate basis functions.
function generate_stress_basis_functions_SubspaceMonomial(subspace_monomial, &
   & structure,complex_modes,real_modes,qpoints,subspaces,                   &
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
  type(StressBasisFunction), allocatable  :: output(:)
  
  ! Monomials, in complex and real representations.
  type(ComplexMonomial), allocatable :: complex_monomials(:)
  type(RealMonomial),    allocatable :: real_monomials(:)
  
  ! The conversion from the complex monomial basis to the real monomial basis,
  !    and vice-versa.
  type(ComplexMatrix) :: complex_to_real_conversion
  type(ComplexMatrix) :: real_to_complex_conversion
  
  ! Each symmetry, acting on scalar basis functions, the stress tensor,
  !    and the combined basis function.
  type(ComplexMatrix)   :: complex_scalar_symmetry
  real(dp), allocatable :: real_scalar_symmetry(:,:)
  real(dp)              :: tensor_symmetry(6,6)
  real(dp), allocatable :: symmetry(:,:)
  type(RealMatrix)      :: projection
  
  ! Polynomial coefficients, in both bases.
  type(SymmetricEigenstuff), allocatable :: estuff(:)
  real(dp),                  allocatable :: real_coefficients(:)
  complex(dp),               allocatable :: complex_coefficients(:)
  
  ! Variables for constructing the output.
  type(ComplexPolynomial) :: complex_representation
  type(ComplexPolynomial) :: elements(3,3)
  
  ! Mappings between 3x3 indices and 6 indices.
  integer :: x(6)
  integer :: y(6)
  
  ! Temporary variables.
  integer  :: i,j,k,l,m,n,o,ialloc
  real(dp) :: tensor(3,3)
  
  ! Initialise index mappings.
  x = [1,2,3,1,1,2]
  y = [1,2,3,2,3,3]
  
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
    allocate(output(0), stat=ialloc); call err(ialloc)
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
  projection = dblemat(make_identity_matrix(size(real_monomials)*6))
  allocate( symmetry(size(real_monomials)*6,size(real_monomials)*6), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(degenerate_symmetries)
    ! Construct the symmetry acting on the complex monomials.
    complex_scalar_symmetry = degenerate_symmetries(i)%calculate_symmetry( &
                                             & complex_monomials,          &
                                             & complex_modes,              &
                                             & include_coefficients=.true. )
    call check_unitary( complex_scalar_symmetry,              &
                      & 'symmetry in complex monomial basis', &
                      & logfile                               )
    
    ! Transform the symmetry into real co-ordinates,
    !    and check that it is real.
    complex_scalar_symmetry = complex_to_real_conversion &
                          & * complex_scalar_symmetry    &
                          & * real_to_complex_conversion
    call check_real( complex_scalar_symmetry,           &
                   & 'symmetry in real monomial basis', &
                   & logfile                            )
    real_scalar_symmetry = dble(real(complex_scalar_symmetry))
    
    ! Construct the symmetry acting on the tensor components.
    do j=1,6
      tensor = 0.0_dp
      if (x(j)==y(j)) then
        tensor(x(j),y(j)) = 1.0_dp
      else
        tensor(x(j),y(j)) = 1/sqrt(2.0_dp)
        tensor(y(j),x(j)) = 1/sqrt(2.0_dp)
      endif
      tensor = dble( structure%symmetries(i)%cartesian_tensor         &
                 & * mat(tensor)                                      &
                 & * invert(structure%symmetries(i)%cartesian_tensor) )
      do k=1,6
        if (x(k)==y(k)) then
          tensor_symmetry(k,j) = tensor(x(k),y(k))
        else
          tensor_symmetry(k,j) = (tensor(x(k),y(k)) + tensor(y(k),x(k))) &
                             & / sqrt(2.0_dp)
        endif
      enddo
    enddo
    call check_orthogonal( mat(tensor_symmetry),                     &
                         & 'symmetry in basis of tensor components', &
                         & logfile                                   )
    
    ! Construct the full symmetry, transforming both scalar and tensor
    !    compontents.
    l = 0
    do j=1,6
      do k=1,size(real_monomials)
        l = l+1
        
        o = 0
        do m=1,6
          do n=1,size(real_monomials)
            o = o+1
            
            symmetry(o,l) = tensor_symmetry(m,j) * real_scalar_symmetry(n,k)
          enddo
        enddo
      enddo
    enddo
    
    ! Construct the projection matrix for this symmetry,
    !    and multiply the total projection matrix by this.
    projection = projection                                                  &
             & * projection_matrix( mat(symmetry),                           &
             &                      structure%symmetries(i)%symmetry_order() )
  enddo
  call check_symmetric( projection,               &
                      & 'projection matrix',      &
                      & logfile,                  &
                      & ignore_threshold=1e-10_dp )
  
  ! Diagonalise the projection matrix,
  !    check its eigenvalues are either 0 or 1,
  !    and select only the eigenvectors with eigenvalue 1.
  estuff = diagonalise_symmetric(projection)
  if (any(abs(estuff%eval-1)>1e-2_dp .and. abs(estuff%eval)>1e-2_dp)) then
    call print_line(ERROR//': Projection matrix has eigenvalues which are &
       &neither 0 nor 1.')
    call print_line('Eigenvalues:')
    call print_line(estuff%eval)
    call err()
  endif
  estuff = estuff(filter(abs(estuff%eval-1)<1e-2_dp))
  
  ! Construct basis functions from coefficients.
  allocate(output(size(estuff)), stat=ialloc); call err(ialloc)
  do i=1,size(estuff)
    do j=1,6
      real_coefficients = estuff(i)%evec( size(real_monomials)*(j-1)+1 &
                                      & : size(real_monomials)*j       )
      
      complex_coefficients = cmplx( real_to_complex_conversion &
                                & * vec(real_coefficients)     )
      
      complex_representation = ComplexPolynomial( complex_coefficients &
                                              & * complex_monomials    )
      
      if (x(j)==y(j)) then
        elements(x(j),y(j)) = complex_representation
        call elements(x(j),y(j))%simplify()
      else
        elements(x(j),y(j)) = complex_representation/sqrt(2.0_dp)
        call elements(x(j),y(j))%simplify()
        elements(y(j),x(j)) = complex_representation/sqrt(2.0_dp)
        call elements(y(j),x(j))%simplify()
      endif
    enddo
    
    output(i) = StressBasisFunction(elements)
  enddo
end function

! Given an orthogonal matrix U, s.t. U^n=I, returns the matrix
!    H = sum(j=0,n-1) U^j / n
! H is symmetric, and is a projection matrix which projects out the
!    eigenvectors of U with eigenvalue 1.
! Uses Horner's rule to calculate the sum with minimal matrix multiplication.
function projection_matrix(input,order) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: input
  integer,          intent(in) :: order
  type(RealMatrix)             :: output
  
  type(RealMatrix) :: identity
  
  integer :: i
  
  if (order<1) then
    call print_line(CODE_ERROR//': symmetry order may not be < 1.')
    call err()
  elseif (order>6) then
    call print_line(CODE_ERROR//': symmetry order may not be > 6.')
    call err()
  endif
  
  identity = dblemat(make_identity_matrix(size(input,1)))
  
  output = identity
  do i=2,order
    output = input*output + identity
  enddo
  output = output/order
end function

! ----------------------------------------------------------------------
! Simplify the basis function.
! ----------------------------------------------------------------------
impure elemental subroutine simplify_StressBasisFunction(this)
  implicit none
  
  class(StressBasisFunction), intent(inout) :: this
  
  call this%elements_%simplify()
end subroutine

! ----------------------------------------------------------------------
! Evaluate the stress due to the basis function.
! ----------------------------------------------------------------------
impure elemental function stress_RealModeDisplacement_StressBasisFunction( &
   & this,displacement) result(output)
  implicit none
  
  class(StressBasisFunction),  intent(in) :: this
  class(RealModeDisplacement), intent(in) :: displacement
  type(RealMatrix)                        :: output
  
  real(dp) :: elements(3,3)
  
  integer :: i,j
  
  do i=1,3
    do j=1,3
      elements(j,i) = real(this%elements_(j,i)%energy(displacement))
    enddo
  enddo
  
  output = this%coefficient_ * mat(elements)
end function

impure elemental function stress_ComplexModeDisplacement_StressBasisFunction( &
   & this,displacement) result(output)
  implicit none
  
  class(StressBasisFunction),     intent(in) :: this
  class(ComplexModeDisplacement), intent(in) :: displacement
  type(ComplexMatrix)                        :: output
  
  complex(dp) :: elements(3,3)
  
  integer :: i,j
  
  do i=1,3
    do j=1,3
      elements(j,i) = this%elements_(j,i)%energy(displacement)
    enddo
  enddo
  
  output = this%coefficient_ * mat(elements)
end function

! ----------------------------------------------------------------------
! Integrate the basis function between two states.
! ----------------------------------------------------------------------
subroutine braket_SubspaceBraKet_StressBasisFunction(this,braket, &
   & anharmonic_data)
  implicit none
  
  class(StressBasisFunction), intent(inout) :: this
  class(SubspaceBraKet),      intent(in)    :: braket
  type(AnharmonicData),       intent(in)    :: anharmonic_data
  
  integer :: i,j
  
  do i=1,3
    do j=1,3
      call integrate( this%elements_(j,i)%terms, &
                    & braket,                    &
                    & anharmonic_data            )
    enddo
  enddo
end subroutine

subroutine braket_BasisState_StressBasisFunction(this,bra,ket,subspace, &
   & subspace_basis,anharmonic_data)
  implicit none
  
  class(StressBasisFunction), intent(inout)        :: this
  class(BasisState),          intent(in)           :: bra
  class(BasisState),          intent(in), optional :: ket
  type(DegenerateSubspace),   intent(in)           :: subspace
  class(SubspaceBasis),       intent(in)           :: subspace_basis
  type(AnharmonicData),       intent(in)           :: anharmonic_data
  
  integer :: i,j
  
  do i=1,3
    do j=1,3
      call integrate( this%elements_(j,i)%terms, &
                    & bra,                       &
                    & ket,                       &
                    & subspace,                  &
                    & subspace_basis,            &
                    & anharmonic_data            )
    enddo
  enddo
end subroutine

subroutine braket_BasisStates_StressBasisFunction(this,states,thermal_energy, &
   & subspace,subspace_basis,anharmonic_data)
  implicit none
  
  class(StressBasisFunction), intent(inout) :: this
  class(BasisStates),         intent(inout) :: states
  real(dp),                   intent(in)    :: thermal_energy
  type(DegenerateSubspace),   intent(in)    :: subspace
  class(SubspaceBasis),       intent(in)    :: subspace_basis
  type(AnharmonicData),       intent(in)    :: anharmonic_data
  
  integer :: i,j,k
  
  do i=1,3
    do j=1,3
      do k=1,size(this%elements_(j,i)%terms)
        call integrate( this%elements_(j,i)%terms(k), &
                      & states,                       &
                      & thermal_energy,               &
                      & subspace,                     &
                      & subspace_basis,               &
                      & anharmonic_data               )
      enddo
    enddo
  enddo
end subroutine

! ----------------------------------------------------------------------
! Returns the thermal expectation of the basis function.
! ----------------------------------------------------------------------
impure elemental function harmonic_expectation_StressBasisFunction(this, &
   & frequency,thermal_energy,supercell_size) result(output)
  implicit none
  
  class(StressBasisFunction), intent(in) :: this
  real(dp),                   intent(in) :: frequency
  real(dp),                   intent(in) :: thermal_energy
  integer,                    intent(in) :: supercell_size
  type(RealMatrix)                       :: output
  
  real(dp) :: elements(3,3)
  
  integer :: i,j
  
  do i=1,3
    do j=1,3
      elements(j,i) = this%elements_(j,i)%harmonic_expectation( &
                                              & frequency,      &
                                              & thermal_energy, &
                                              & supercell_size  )
    enddo
  enddo
  
  output = this%coefficient_ * mat(elements)
end function

! ----------------------------------------------------------------------
! Returns the stress at zero displacement.
! ----------------------------------------------------------------------
impure elemental function undisplaced_stress_StressBasisFunction(this) &
   & result(output)
  implicit none
  
  class(StressBasisFunction), intent(in) :: this
  type(RealMatrix)                       :: output
  
  type(RealSingleDisplacement) :: zero_displacement(0)
  
  output = this%stress(RealModeDisplacement(zero_displacement))
end function

! ----------------------------------------------------------------------
! Arithmetic.
! ----------------------------------------------------------------------
impure elemental function multiply_StressBasisFunction_real(this,that) &
   & result(output)
  implicit none
  
  type(StressBasisFunction), intent(in) :: this
  real(dp),                  intent(in) :: that
  type(StressBasisFunction)             :: output
  
  output = this
  output%coefficient_ = output%coefficient_ * that
end function

impure elemental function multiply_real_StressBasisFunction(this,that) &
   & result(output)
  implicit none
  
  real(dp),                  intent(in) :: this
  type(StressBasisFunction), intent(in) :: that
  type(StressBasisFunction)             :: output
  
  output = that
  output%coefficient_ = this * output%coefficient_
end function

impure elemental function divide_StressBasisFunction_real(this,that) &
   & result(output)
  implicit none
  
  type(StressBasisFunction), intent(in) :: this
  real(dp),                  intent(in) :: that
  type(StressBasisFunction)             :: output
  
  output = this
  output%coefficient_ = output%coefficient_ / that
end function

! Calculate the contribution to a given monomial from the interpolation of
!    this basis function.
! The result is given as a cartesian tensor.
impure elemental function interpolate_coefficients_StressBasisFunction(this, &
   & monomial,interpolator) result(output)
  implicit none
  
  class(StressBasisFunction),   intent(in) :: this
  type(ComplexMonomial),        intent(in) :: monomial
  type(PolynomialInterpolator), intent(in) :: interpolator
  type(ComplexMatrix)                      :: output
  
  complex(dp) :: elements(3,3)
  
  integer :: i,j
  
  do i=1,3
    do j=i,3
      elements(j,i) = interpolator%overlap(monomial, this%elements_(j,i)) &
                  & * this%coefficient_
    enddo
  enddo
  
  elements(1,2) = elements(2,1)
  elements(1,3) = elements(3,1)
  elements(2,3) = elements(3,2)
  
  output = mat(elements)
end function

! Calculate this basis function's contribution to the effective stress
!    dynamical matrix from which the stress can be interpolated in the
!    large-supercell limit.
function calculate_dynamical_matrices_StressBasisFunction(this,qpoints, &
   & thermal_energy,subspaces,subspace_bases,subspace_states,           &
   & subspaces_in_coupling,anharmonic_data) result(output) 
  implicit none
  
  class(StressBasisFunction), intent(in)    :: this
  type(QpointData),           intent(in)    :: qpoints(:)
  real(dp),                   intent(in)    :: thermal_energy
  type(DegenerateSubspace),   intent(in)    :: subspaces(:)
  class(SubspaceBasis),       intent(in)    :: subspace_bases(:)
  class(BasisStates),         intent(inout) :: subspace_states(:)
  integer,                    intent(in)    :: subspaces_in_coupling(:)
  type(AnharmonicData),       intent(in)    :: anharmonic_data
  type(StressDynamicalMatrix), allocatable  :: output(:)
  
  type(DynamicalMatrix), allocatable :: matrices(:)
  
  integer :: i,j,k
  
  output = [( StressDynamicalMatrix(anharmonic_data%structure%no_atoms), &
            & i=1,                                                       &
            & size(qpoints)                                              )]
  
  do i=1,3
    do j=1,3
      matrices = calculate_dynamical_matrices( this%elements_(j,i),   &
                                             & qpoints,               &
                                             & thermal_energy,        &
                                             & subspaces,             &
                                             & subspace_bases,        &
                                             & subspace_states,       &
                                             & subspaces_in_coupling, &
                                             & anharmonic_data        )
      do k=1,size(output)
        output(k)%elements(j,i) = matrices(k) * this%coefficient_
      enddo
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_StressBasisFunction(this,input)
  implicit none
  
  class(StressBasisFunction), intent(out) :: this
  type(String),               intent(in)  :: input(:)
  
  type(ComplexMonomial) :: no_monomials(0)
  
  type(ComplexPolynomial) :: elements(3,3)
  
  type(StringArray), allocatable :: sections(:)
  
  integer :: i,j,k
  
  select type(this); type is(StressBasisFunction)
    elements = ComplexPolynomial(no_monomials)
    
    sections = split_into_sections(input)
    do k=1,size(sections)
      i = int(token(sections(k)%strings(1),3))
      j = int(token(sections(k)%strings(1),5))
      elements(i,j) = ComplexPolynomial(join( sections(k)%strings(2:), &
                                            & delimiter=' + '          ))
    enddo
    
    this = StressBasisFunction(elements)
  class default
    call err()
  end select
end subroutine

function write_StressBasisFunction(this) result(output)
  implicit none
  
  class(StressBasisFunction), intent(in) :: this
  type(String), allocatable              :: output(:)
  
  character(1), parameter :: directions(3) = ['x', 'y', 'z']
  
  logical :: first_element
  
  integer :: i,j,ialloc
  
  select type(this); type is(StressBasisFunction)
    allocate(output(0), stat=ialloc); call err(ialloc)
    first_element = .true.
    do i=1,3
      do j=1,3
        if (size(this%elements_(i,j))>0) then
          if (.not. first_element) then
            output = [output, str('')]
          endif
          output = [                                                        &
             & output,                                                      &
             & str('Stress component '//directions(i)//directions(j)//':'), &
             & str(this%elements_(i,j)*this%coefficient_)                   ]
          first_element = .false.
        endif
      enddo
    enddo
  class default
    call err()
  end select
end function

function new_StressBasisFunction_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)  :: input(:)
  type(StressBasisFunction) :: this
  
  call this%read(input)
end function

impure elemental function new_StressBasisFunction_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(StressBasisFunction)     :: this
  
  this = StressBasisFunction(str(input))
end function
end module
