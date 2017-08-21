! ======================================================================
! The basis functions from which the Born-Oppenheimer potential is made.
! ======================================================================
module potential_module
  use constants_module, only : dp
  use string_module
  use io_module
  use linear_algebra_module
  
  private
  
  public :: calculate_potential
  
  ! A monomial, e.g.
  !    C * (u1)**a * (u2)**b * (u4)**d => coef=C, powers=[a,b,0,d]
  type :: Monomial
    real(dp)             :: coefficient
    integer, allocatable :: powers(:)
  contains
    procedure :: evaluate_energy => evaluate_energy_Monomial
    procedure :: evaluate_forces => evaluate_forces_Monomial
  end type
  
  ! A polynomial representation of a potential.
  type, public :: PolynomialPotential
    type(Monomial),   allocatable, private :: monomials(:)
    type(RealMatrix), allocatable, private :: harmonic_couplings(:,:)
  contains
    procedure, public :: evaluate_energy => evaluate_energy_PolynomialPotential
    procedure, public :: evaluate_forces => evaluate_forces_PolynomialPotential
    procedure, public :: simplify => simplify_PolynomialPotential
    procedure, public :: integrate_over_mode_average => &
                       & integrate_over_mode_average_PolynomialPotential
    procedure, public :: construct_hamiltonian => &
                       & construct_hamiltonian_PolynomialPotential
  end type
contains

! ----------------------------------------------------------------------
! Evaluates the energy of a Monomial at a given displacement.
! ----------------------------------------------------------------------
function evaluate_energy_Monomial(this,displacement) result(output)
  use sampling_points_module
  use normal_mode_module
  implicit none
  
  class(Monomial),  intent(in) :: this
  type(ModeVector), intent(in) :: displacement
  real(dp)                     :: output
  
  integer :: i
  
  output = this%coefficient
  do i=1,size(this%powers)
    output = output * displacement%vector(i)**this%powers(i)
  enddo
end function

! ----------------------------------------------------------------------
! Evaluates the force of a Monomial at a given displacement.
! Returns the result in normal mode co-ordinates.
! ----------------------------------------------------------------------
function evaluate_forces_Monomial(this,displacement) result(output)
  use sampling_points_module
  use normal_mode_module
  implicit none
  
  class(Monomial),  intent(in) :: this
  type(ModeVector), intent(in) :: displacement
  type(ModeVector)             :: output
  
  type(Monomial) :: derivative
  
  integer :: no_modes
  integer :: i,ialloc
  
  no_modes = size(this%powers)
  allocate(output%vector(no_modes), stat=ialloc); call err(ialloc)
  do i=1,no_modes
    derivative = this
    derivative%coefficient = derivative%coefficient &
                         & * derivative%powers(i)
    derivative%powers(i) = derivative%powers(i) - 1
    output%vector(i) = derivative%evaluate_energy(displacement)
  enddo
end function

! ----------------------------------------------------------------------
! Evaluates the energy of a potential at a given displacement.
! ----------------------------------------------------------------------
function evaluate_energy_PolynomialPotential(this,displacement) result(output)
  use sampling_points_module
  use normal_mode_module
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  type(ModeVector),           intent(in) :: displacement
  real(dp)                               :: output
  
  integer :: i
  
  output = 0
  do i=1,size(this%monomials)
    output = output + this%monomials(i)%evaluate_energy(displacement)
  enddo
end function

! ----------------------------------------------------------------------
! Evaluates the force of a potential at a given displacement.
! Returns the result in normal mode co-ordinates.
! ----------------------------------------------------------------------
function evaluate_forces_PolynomialPotential(this,displacement) result(output)
  use sampling_points_module
  use normal_mode_module
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  type(ModeVector),           intent(in) :: displacement
  type(ModeVector)                       :: output
  
  integer          :: no_modes
  type(ModeVector) :: temp
  
  integer :: i,ialloc
  
  no_modes = size(this%monomials(1)%powers)
  allocate(output%vector(no_modes), stat=ialloc); call err(ialloc)
  output%vector = 0
  do i=1,size(this%monomials)
    temp = this%monomials(i)%evaluate_forces(displacement)
    output%vector = output%vector + temp%vector
  enddo
end function

! ----------------------------------------------------------------------
! Calculates a representation of the Born-Oppenheimer surface using
!    samples from the surface.
! ----------------------------------------------------------------------
function calculate_potential(no_basis_functions,sampling,modes,qpoint, &
   & supercell,harmonic_states) result(output)
  use sampling_points_module
  use normal_mode_module
  use qpoints_module
  use structure_module
  use eigenstates_module
  implicit none
  
  integer,                intent(in) :: no_basis_functions
  type(CouplingSampling), intent(in) :: sampling(:)
  type(NormalMode),       intent(in) :: modes(:)
  type(QpointData),       intent(in) :: qpoint
  type(StructureData),    intent(in) :: supercell
  type(SingleModeState),  intent(in) :: harmonic_states(:,:)
  type(PolynomialPotential)          :: output
  
  integer :: no_modes
  integer :: no_states
  
  real(dp),                  allocatable :: gaussian_integrals(:)
  type(PolynomialPotential), allocatable :: basis_functions(:)
  real(dp),                  allocatable :: basis_coefficients(:)
  
  integer :: i,ialloc
  
  no_modes = size(harmonic_states,2)
  no_states = size(harmonic_states,1)
  
  ! Pre-process harmonic couplings (saves repeated calculation later).
  allocate( output%harmonic_couplings(no_basis_functions,no_modes), &
          & stat=ialloc); call err(ialloc)
  do i=1,no_modes
    ! Calculate gaussian integrals, I(n).
    ! I(n+1) = integral[ u^n e^(-freq*u*u) ].
    gaussian_integrals = calculate_gaussian_integrals( &
                     & harmonic_states(1,i)%frequency, &
                     & no_states,                      &
                     & no_basis_functions)
    
    ! Calculate coupling between harmonic states and u^n potentials, V_ij.
    ! V_ij(n+1) = <i|u^n|j>.
    output%harmonic_couplings(:,i) = generate_harmonic_couplings( &
                                          & harmonic_states(:,i), &
                                          & no_basis_functions,   &
                                          & gaussian_integrals)
  enddo
  
  ! Generate basis functions
  basis_functions = calculate_all_basis_functions(sampling,no_modes, &
     & no_basis_functions)
  
  ! Use linear least-squares regression to fit basis function coefficients
  !    to input energy and force data.
  basis_coefficients = fit_basis_functions(sampling,modes,qpoint,supercell, &
     & basis_functions)
  
  ! Combine basis functions into a single sum of monomials
  ! TODO
  
  ! Simplify output.
  call output%simplify()
end function

! ----------------------------------------------------------------------
! Generate all polynomial basis functions.
! At present, all basis functions are monomials with coefficient 1.
! ----------------------------------------------------------------------
! Basis functions are set up in a generalised octahedral manner, i.e.
!    sum(output(i)%monomials(j)%powers) in set [0,no_basis_functions]
! e.g. if no_basis_functions=3 then the functions u1^3, u1^2*u2 and u1*u2*u3
!    are all allowed, but e.g. u1^2*u2^2 is not.
! If there are n modes in the coupling, and no_basis_functions=k,
!    then the number of basis functions is given by the binomial coefficient
!  / n+k \    (n+k)!
! |       | = -----
!  \  k  /     n!k!

function calculate_all_basis_functions(sampling,no_modes,no_basis_functions) &
   & result(output)
  use utils_module, only : factorial
  use coupling_module
  use sampling_points_module
  implicit none
  
  type(CouplingSampling), intent(in)     :: sampling(:)
  integer,                intent(in)     :: no_modes
  integer,                intent(in)     :: no_basis_functions
  type(PolynomialPotential), allocatable :: output(:)
  
  integer              :: max_size
  integer, allocatable :: sizes(:)
  
  integer :: i,j,ialloc
  
  ! Find the size of the largest coupling.
  max_size = 0
  do i=1,size(sampling)
    max_size = max(max_size, size(sampling(i)%coupling))
  enddo
  
  ! Calculate how much space is needed and allocate output.
  allocate(sizes(size(sampling)), stat=ialloc); call err(ialloc)
  do i=1,size(sampling)
    sizes(i) = factorial(size(sampling(i)%coupling)+no_basis_functions) &
           & / ( factorial(size(sampling(i)%coupling))                  &
           &   * factorial(no_basis_functions) )
  enddo
  allocate(output(sum(sizes)), stat=ialloc); call err(ialloc)
  
  ! Calculate output, coupling by coupling.
  j = 0
  do i=1,size(sampling)
    output(j+1:j+sizes(i)) = calculate_basis_functions( sampling(i)%coupling%modes, &
                                                      & no_modes,           &
                                                      & no_basis_functions)
    j = j+sizes(i)
  enddo
end function

function calculate_basis_functions(coupling,no_modes,no_basis_functions) &
   & result(output)
  use grid_types_module
  implicit none
  
  integer, intent(in)                    :: coupling(:)
  integer, intent(in)                    :: no_modes
  integer, intent(in)                    :: no_basis_functions
  type(PolynomialPotential), allocatable :: output(:)
  
  integer, allocatable :: octahedral_points(:,:)
  
  integer :: i,j,ialloc
  
  ! Generate an octahedral grid of points.
  ! Each point will be mapped onto a single basis function.
  octahedral_points = generate_octahedral_grid( size(coupling), &
                                              & no_basis_functions)
  
  ! Generate output.
  allocate(output(size(octahedral_points,2)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    ! Each basis function consists of only one monomial.
    allocate(output(i)%monomials(1), stat=ialloc); call err(ialloc)
    output(i)%monomials(1)%coefficient = 1
    output(i)%monomials(1)%powers = int(zeroes(no_modes))
    do j=1,size(coupling)
      output(i)%monomials(1)%powers(coupling(j)) = octahedral_points(j,i)
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Simplifies a polynomial potential.
! Merges all equal terms, summing over their coefficients.
! ----------------------------------------------------------------------
subroutine simplify_PolynomialPotential(this)
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  
  type(Monomial), allocatable :: potential(:)
  integer,        allocatable :: unique_id(:)
  
  integer :: i,j,id,ialloc
  
  potential = this%monomials
  
  allocate(unique_id(size(potential)), stat=ialloc); call err(ialloc)
  id = 0
  do_i : do i=1,size(potential)
    do j=1,i-1
      if (all(potential(i)%powers==potential(j)%powers)) then
        unique_id(i) = unique_id(j)
        cycle do_i
      endif
    enddo
    
    id = id + 1
    unique_id(i) = id
  enddo do_i
  
  deallocate(this%monomials, stat=ialloc); call err(ialloc)
  allocate(this%monomials(id), stat=ialloc); call err(ialloc)
  do i=1,size(this%monomials)
    this%monomials(i)%coefficient = 0
  enddo
  
  do i=1,size(potential)
    id = unique_id(i)
    this%monomials(id)%coefficient = this%monomials(id)%coefficient &
                                 & + potential(i)%coefficient
    this%monomials(id)%powers = potential(i)%powers
  enddo
end subroutine

! ----------------------------------------------------------------------
! Fits a set of basis functions to the sampled B-O surface.
! ----------------------------------------------------------------------
function fit_basis_functions(sampling,modes,qpoint,supercell,basis_functions) &
   & result(output)
  use sampling_points_module
  use normal_mode_module
  use qpoints_module
  use structure_module
  use grid_types_module
  implicit none
  
  type(CouplingSampling),    intent(in) :: sampling(:)
  type(NormalMode),          intent(in) :: modes(:)
  type(QpointData),          intent(in) :: qpoint
  type(StructureData),       intent(in) :: supercell
  type(PolynomialPotential), intent(in) :: basis_functions(:)
  real(dp), allocatable                 :: output(:)
  
  ! Input array lengths.
  integer :: no_forces
  
  ! Independent sampling points.
  integer                          :: no_sampling_points
  type(SamplingPoint), allocatable :: sampling_points(:)
  real(dp),            allocatable :: energy(:)
  type(RealVector),    allocatable :: forces(:,:)
  
  ! The normal-mode displacement at each sampling point.
  type(ModeVector) :: displacement
  
  ! The basis-function force at each sampling point.
  type(ModeVector)              :: forces_mode
  type(RealVector), allocatable :: forces_cartesian(:)
  
  ! Linear least-squares matrix and vector, L=(a.x-b)**2.
  integer               :: m,n
  real(dp), allocatable :: a(:,:)
  real(dp), allocatable :: b(:)
  
  ! Temporary variables.
  integer :: i,j,k,l,ialloc
  
  no_forces = size(sampling(1)%forces)
  
  ! Concatenate all independent sampling points.
  no_sampling_points = 0
  do i=1,size(sampling)
    do j=1,size(sampling(i)%sampling_points)
      if (.not. sampling(i)%sampling_points(j)%duplicate) then
        no_sampling_points = no_sampling_points + 1
      endif
    enddo
  enddo
  k = 0
  allocate( sampling_points(no_sampling_points),  &
          & energy(no_sampling_points),           &
          & forces(no_forces,no_sampling_points), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(sampling)
    do j=1,size(sampling(i)%sampling_points)
      if (.not. sampling(i)%sampling_points(j)%duplicate) then
        k = k+1
        sampling_points(k) = sampling(i)%sampling_points(j)
        energy(k) = sampling(i)%energy(j)
        forces(:,k) = sampling(i)%forces(:,j)
      endif
    enddo
  enddo
  
  ! Construct a and b.
  m = no_sampling_points*(1+no_forces)
  n = size(basis_functions)
  allocate( a(m,n), &
          & b(m),   &
          & stat=ialloc); call err(ialloc)
  do i=1,no_sampling_points
    l = (i-1)*(1+no_forces) + 1
    displacement = sampling_points(i)%displacement
    
    do j=1,n
      a(l,j) = basis_functions(j)%evaluate_energy(displacement)
      forces_mode = basis_functions(j)%evaluate_forces(displacement)
      forces_cartesian = normal_mode_to_cartesian( forces_mode, &
                                                 & modes,       &
                                                 & qpoint,      &
                                                 & supercell)
      do k=1,size(forces_cartesian)
        a(l+3*k-2:l+3*k, j) = dble(forces_cartesian(k))
      enddo
    enddo
    
    b(l) = energy(i)
    do k=1,size(forces(:,i))
      b(l+3*k-2:l+3*k) = dble(forces(k,i))
    enddo
  enddo
  
  ! Perform linear least-squares optimisation.
  output = dble(linear_least_squares(a,b))
end function

! ----------------------------------------------------------------------
! Calculates integrals of the form I(n) = u^n*e^(-freq*u*u)
! ----------------------------------------------------------------------
! N.B. output(i) = I(i-1) because I(0) exists.
! Uses the relations:
!    I(n) = integral[ u^n*e^(-freq*u*u) ] from -infinity to infinity.
!    I(0) = sqrt(pi/freq)
!    I(1) = 0
!    I(n) = (n-1)/(2*freq) * I(n-2)
function calculate_gaussian_integrals(frequency,no_harmonic_states, &
   & no_basis_functions) result(output)
  use constants_module, only : pi
  implicit none
  
  real(dp), intent(in)  :: frequency
  integer,  intent(in)  :: no_harmonic_states
  integer,  intent(in)  :: no_basis_functions
  real(dp), allocatable :: output(:)
  
  integer :: i,ialloc
  
  allocate( output(2*no_harmonic_states+no_basis_functions), &
          & stat=ialloc); call err(ialloc)
  output(1) = sqrt(pi/frequency)
  output(2) = 0
  do i=3,size(output)
    output(i) = output(i-2)*i/(2*frequency)
  enddo
end function

! ----------------------------------------------------------------------
! Generates all elements <i|u^k|j> along a given mode.
! ----------------------------------------------------------------------
! output(k) is a matrix whose i,j element is <i-1|u^(k-1)|j-1>.
! The -1 offsets are due to |0> and u^0.
function generate_harmonic_couplings(harmonic_states,no_basis_functions, &
   & gaussian_integrals) result(output)
  use linear_algebra_module
  use eigenstates_module
  implicit none
  
  type(SingleModeState), intent(in) :: harmonic_states(:)
  integer,               intent(in) :: no_basis_functions
  real(dp),              intent(in) :: gaussian_integrals(:)
  type(RealMatrix), allocatable     :: output(:)
  
  real(dp), allocatable :: matrix(:,:)
  
  integer :: i,j,k,i2,j2,ialloc
  
  allocate( matrix(size(harmonic_states),size(harmonic_states)), &
          & output(no_basis_functions),                          &
          & stat=ialloc); call err(ialloc)
  
  do k=1,no_basis_functions
    ! Calculate output(k) = {<i|u^k|j>}
    matrix = 0
    do j=1,size(harmonic_states)
      do i=1,j
        do j2=1,size(harmonic_states(j)%coefficients)
          do i2=1,size(harmonic_states(i)%coefficients)
            ! N.B. the -2 is because i2 refers to u^(i2-1) etc.
            matrix(i,j) = matrix(i,j)                         &
                      & + harmonic_states(j)%coefficients(j2) &
                      & * harmonic_states(i)%coefficients(i2) &
                      & * gaussian_integrals(i2+j2+k-2)
          enddo
        enddo
        matrix(j,i) = matrix(i,j)
      enddo
    enddo
    output(k) = matrix
  enddo
end function

! ----------------------------------------------------------------------
! Takes a potential, a mode id, the harmonic couplings of that mode,
!    and the anharmonic states of that mode in terms of harmonic states.
! Integrates the potential between the average of the mode.
! ----------------------------------------------------------------------
! 
! f( V(u1,u2,u3,...), 1, {<u1|u^n|u1>}, {<U1|u1>} )   ->   V(u2,u3,...)
! f( V(u1,u2,u3,...), 2, {<u2|u^n|u2>}, {<U2|u2>} )   ->   V(u1,u3,...)
!
! V(u2,u3,...) = sum_i(<U1_i|V(u1,u2,u3,...)|U1_i>) / sum_i(1)
!              = sum_{i,j,k}(a_ij V_jk a_ik^T) / sum_i(1)
!
! where |u1> is a harmonic eigenstate, and |U1> is an anharmonic eigenstate.
!
! |U1_i> = a1_ij|u1_j>, where a1=eigenstuff%evecs=<U1_i|u1_j>
subroutine integrate_over_mode_average_PolynomialPotential(this,mode, &
   & eigenstuff)
  use linear_algebra_module
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  integer,                    intent(in)    :: mode
  type(RealEigenstuff),       intent(in)    :: eigenstuff
  
  integer :: i
  
  ! Integrate across each monomial.
  do i=1,size(this%monomials)
    this%monomials(i)%coefficient = this%monomials(i)%coefficient             &
          & + trace( mat(eigenstuff%evecs)                                    &
          &        * this%harmonic_couplings( this%monomials(i)%powers(mode), &
          &                                   mode)                           &
          &        * transpose(mat(eigenstuff%evecs)) )                       &
          & / size(eigenstuff)
    this%monomials(i)%powers(mode) = 0
  enddo
  
  ! Simplify across now identical monomials.
  call this%simplify()
end subroutine

function construct_hamiltonian_PolynomialPotential(this,mode) result(output)
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  integer,                    intent(in) :: mode
  type(RealMatrix)                       :: output
  
  integer :: no_states
  integer :: i
  
  no_states = size(this%harmonic_couplings(1,mode), 1)
  
  output = dble(int(zeroes(no_states,no_states)))
  do i=1,size(this%monomials)
    output = output &
         & + this%monomials(i)%coefficient &
         & * this%harmonic_couplings(this%monomials(i)%powers(mode),mode)
  enddo
end function
end module
