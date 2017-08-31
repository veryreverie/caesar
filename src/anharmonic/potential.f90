! ======================================================================
! The basis functions from which the Born-Oppenheimer potential is made.
! ======================================================================
module potential_module
  use constants_module, only : dp
  use string_module
  use io_module
  use linear_algebra_module
  implicit none
  
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
    integer,                       private :: cutoff
    type(Monomial),   allocatable, private :: monomials(:)
    real(dp),         allocatable, private :: gaussian_integrals(:,:)
  contains
    procedure, public  :: evaluate_energy => &
                        & evaluate_energy_PolynomialPotential
    procedure, public  :: evaluate_forces => &
                        & evaluate_forces_PolynomialPotential
    procedure, public  :: simplify => simplify_PolynomialPotential
    procedure, private :: calculate_gaussian_integrals => &
                        & calculate_gaussian_integrals_PolynomialPotential
    procedure, private :: calculate_power_coupling => &
                        & calculate_power_coupling_PolynomialPotential
    procedure, private :: generate_state_couplings => &
                        & generate_state_couplings_PolynomialPotential
    generic,   public  :: integrate =>                                   &
                        & integrate_PolynomialPotential_SingleModeState, &
                        & integrate_PolynomialPotential_ProductState
    procedure, private :: integrate_PolynomialPotential_SingleModeState
    procedure, private :: integrate_PolynomialPotential_ProductState
    procedure, public  :: integrate_to_constant => &
                        & integrate_to_constant_PolynomialPotential
    procedure, public  :: integrate_over_mode_average => &
                        & integrate_over_mode_average_PolynomialPotential
    procedure, public  :: construct_hamiltonian => &
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
  
  integer        :: no_modes
  type(Monomial) :: derivative
  
  integer :: i,ialloc
  
  no_modes = size(displacement%vector)
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
  
  no_modes = size(displacement%vector)
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
function calculate_potential(potential_basis_cutoff,sampling,modes,qpoint, &
   & supercell,harmonic_states_cutoff) result(output)
  use sampling_points_module
  use normal_mode_module
  use qpoints_module
  use structure_module
  implicit none
  
  integer,                intent(in) :: potential_basis_cutoff
  type(CouplingSampling), intent(in) :: sampling(:)
  type(NormalMode),       intent(in) :: modes(:)
  type(QpointData),       intent(in) :: qpoint
  type(StructureData),    intent(in) :: supercell
  integer,                intent(in) :: harmonic_states_cutoff
  type(PolynomialPotential)          :: output
  
  type(PolynomialPotential), allocatable :: basis_functions(:)
  real(dp),                  allocatable :: basis_coefficients(:)
  
  integer, allocatable :: no_monomials(:)
  
  integer :: i,j,k,ialloc
  
  output%cutoff = potential_basis_cutoff
  
  ! Calculate gaussian integrals, I(n).
  ! I(n+1) = integral[ u^n e^(-freq*u*u) ].
  call output%calculate_gaussian_integrals(modes, harmonic_states_cutoff)
  
  ! Generate basis functions
  basis_functions = calculate_basis_functions( sampling, &
                                             & size(modes), &
                                             & potential_basis_cutoff)
  
  ! Use linear least-squares regression to fit basis function coefficients
  !    to input energy and force data.
  basis_coefficients = fit_basis_functions(sampling,modes,qpoint,supercell, &
     & basis_functions)
  
  ! Calculate the number of monomials in each basis function.
  allocate(no_monomials(size(basis_functions)), stat=ialloc); call err(ialloc)
  do i=1,size(basis_functions)
    no_monomials(i) = size(basis_functions(i)%monomials)
  enddo
  
  ! Sum together basis functions * their coefficients.
  allocate(output%monomials(sum(no_monomials)), stat=ialloc); call err(ialloc)
  k = 0
  do i=1,size(basis_functions)
    do j=1,no_monomials(i)
      output%monomials(k+j) = basis_functions(i)%monomials(j)
      output%monomials(k+j)%coefficient = output%monomials(k+j)%coefficient &
                                      & * basis_coefficients(i)
    enddo
    k = k+no_monomials(i)
  enddo
  
  ! Simplify output, adding together all terms with equal powers.
  call output%simplify()
end function

! ----------------------------------------------------------------------
! Generate all polynomial basis functions.
! At present, all basis functions are monomials with coefficient 1.
! ----------------------------------------------------------------------
! Basis functions are set up in a generalised octahedral manner, i.e.
!    sum(output(i)%monomials(j)%powers) is in set [0,potential_basis_cutoff]
! e.g. if potential_basis_cutoff=3 then
!    u1^3, u1^2*u2 and u1*u2*u3 are allowed, but u1^2*u2^2 is not.
! These are generated by first generating an octahedral grid of points, and
!    then mapping the indices of each point to the powers of each monomial.
function calculate_basis_functions(sampling,no_modes, &
   & potential_basis_cutoff) result(output)
  use coupling_module
  use sampling_points_module
  use grid_types_module
  implicit none
  
  type(CouplingSampling), intent(in)     :: sampling(:)
  integer,                intent(in)     :: no_modes
  integer,                intent(in)     :: potential_basis_cutoff
  type(PolynomialPotential), allocatable :: output(:)
  
  integer              :: size_output
  integer, allocatable :: grid(:,:)
  
  integer :: i,j,k,l,ialloc
  
  ! Calculate how much space is needed and allocate output.
  size_output = 0
  do i=1,size(sampling)
    ! This gives the array length of 'grid' below.
    size_output = size_output &
              & + octahedral_grid_size( size(sampling(i)%coupling), &
              &                         potential_basis_cutoff-1,   &
              &                         include_negatives=.false.)
  enddo
  allocate(output(size_output), stat=ialloc); call err(ialloc)
  
  ! Calculate output, coupling by coupling.
  j = 0
  do i=1,size(sampling)
    ! Generate an octahedral grid of points.
    ! Each point will be mapped onto a single basis function.
    ! The -1 comes from not taking any u^0 terms, since they will be handled by
    !    subsidiary couplings.
    grid = generate_octahedral_grid( size(sampling(i)%coupling), &
                                   & potential_basis_cutoff-1,   &
                                   & include_negatives=.false.)
    do k=1,size(grid)
      ! Each basis function consists of only one monomial.
      allocate(output(j+k)%monomials(1), stat=ialloc); call err(ialloc)
      output(j+k)%monomials(1)%coefficient = 1
      output(j+k)%monomials(1)%powers = int(zeroes(no_modes))
      do l=1,size(sampling(i)%coupling)
        ! The +1 here is for the same reason as the -1 above.
        output(j+k)%monomials(1)%powers(sampling(i)%coupling%modes(l)) = &
           & grid(l,k) + 1
      enddo
    enddo
    j = j+size(grid)
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
! output(n+1) = I(n) = integral( u^n*e^(-freq*u*u) ) from u=-infty to infty.
! ----------------------------------------------------------------------
! I(0) = sqrt(pi/freq)
! I(1) = 0
! I(n) = (n-1)/(2*freq) * I(n-2)
subroutine calculate_gaussian_integrals_PolynomialPotential(this,modes, &
   & harmonic_states_cutoff)
  use constants_module, only : pi
  use normal_mode_module
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  type(NormalMode),           intent(in)    :: modes(:)
  integer,                    intent(in)    :: harmonic_states_cutoff
  
  integer :: no_integrals
  
  integer :: i,j,ialloc
  
  ! Integrals range from <u^0|u^0|u^0> to <u^h_s_c|u^cutoff|u^h_s_c>.
  !    i.e. integral[u^0e^~] to integral[u^(2*h_s_c+cutoff)e^~].
  no_integrals = 2*harmonic_states_cutoff+this%cutoff+1
  allocate( this%gaussian_integrals(no_integrals, size(modes)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(modes)
    this%gaussian_integrals(1,i) = sqrt(pi/modes(i)%frequency)
    this%gaussian_integrals(2,i) = 0
    do j=3,size(this%gaussian_integrals,1)
      this%gaussian_integrals(j,i) = this%gaussian_integrals(j-2,i) &
                                 & * j / (2*modes(i)%frequency)
    enddo
  enddo
end subroutine

! ----------------------------------------------------------------------
! Takes two single-mode states, |bra> and |ket>, and returns
!    output(n+1) = <bra|u^n|ket>.
! ----------------------------------------------------------------------
function calculate_power_coupling_PolynomialPotential(this,mode,bra,ket) &
   & result(output)
  use eigenstates_module
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  integer,                    intent(in) :: mode
  type(SingleModeState),      intent(in) :: bra
  type(SingleModeState),      intent(in) :: ket
  real(dp), allocatable                  :: output(:)
  
  integer :: i,j,k,ialloc
  
  allocate(output(this%cutoff+1), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    do j=1,size(bra%coefficients)
      do k=1,size(ket%coefficients)
        ! N.B. the -2 is because i refers to u^(i-1) etc.
        output(i) = output(i)           &
                & + bra%coefficients(j) &
                & * ket%coefficients(k) &
                & * this%gaussian_integrals(i+j+k-2,mode)
      enddo
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Generates all elements <i|u^k|j> along a given mode.
! ----------------------------------------------------------------------
! output(k) is a matrix whose i,j element is <i-1|u^(k-1)|j-1>.
! The -1 offsets are due to |0> and u^0.
function generate_state_couplings_PolynomialPotential(this,mode,states) &
   & result(output)
  use linear_algebra_module
  use eigenstates_module
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  integer,                    intent(in) :: mode
  type(SingleModeState),      intent(in) :: states(:)
  type(RealMatrix), allocatable          :: output(:)
  
  real(dp), allocatable :: temp(:,:,:)
  real(dp), allocatable :: coupling(:)
  
  integer :: i,j,ialloc
  
  allocate( temp(size(states),size(states),this%cutoff+1), &
          & output(this%cutoff+1),                         &
          & stat=ialloc); call err(ialloc)
  
  do i=1,size(states)
    do j=1,i
      coupling = this%calculate_power_coupling(mode,states(i),states(j))
      temp(i,j,:) = coupling
      temp(j,i,:) = coupling
    enddo
  enddo
  
  do i=1,size(temp,3)
    output(i) = mat(temp(:,:,i))
  enddo
end function

! ----------------------------------------------------------------------
! Takes two single-mode states, |bra> and |ket>, and forms <bra|V|ket>.
! ----------------------------------------------------------------------
subroutine integrate_PolynomialPotential_SingleModeState(this,mode,bra,ket)
  use eigenstates_module
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  integer,                    intent(in)    :: mode
  type(SingleModeState),      intent(in)    :: bra
  type(SingleModeState),      intent(in)    :: ket
  
  real(dp), allocatable :: couplings(:)
  
  integer :: i
  
  couplings = this%calculate_power_coupling(mode,bra,ket)
  
  do i=1,size(this%monomials)
    ! (a) * u1^n1*...um^nm*... -> (a*<bra|um^nm|ket>) * u1^n1*...*um^0*...
    this%monomials(i)%coefficient = this%monomials(i)%coefficient &
                                & * couplings(this%monomials(i)%powers(mode))
    this%monomials(i)%powers(mode) = 0
  enddo
  
  ! Simplify across now identical monomials.
  call this%simplify()
end subroutine

! ----------------------------------------------------------------------
! Takes two product states, |bra> and |ket>, and forms <bra|V|ket>
! ----------------------------------------------------------------------
subroutine integrate_PolynomialPotential_ProductState(this,bra,ket)
  use eigenstates_module
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  type(ProductState),         intent(in)    :: bra
  type(ProductState),         intent(in)    :: ket
  
  integer :: i
  
  if (size(ket%states)/=size(bra%states)) then
    call err()
  endif
  
  do i=1,size(ket%states)
    call this%integrate(i,bra%states(i),ket%states(i))
  enddo
  
  if (size(this%monomials)/=1) then
    call err()
  endif
end subroutine

! ----------------------------------------------------------------------
! As above, but returns the coefficient rather than modifying the potential.
! ----------------------------------------------------------------------
function integrate_to_constant_PolynomialPotential(this,bra,ket) result(output)
  use eigenstates_module
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  type(ProductState),         intent(in) :: bra
  type(ProductState),         intent(in) :: ket
  real(dp)                               :: output
  
  type(PolynomialPotential) :: temp
  
  temp = this
  
  call temp%integrate(bra,ket)
  
  output = temp%monomials(1)%coefficient
end function

! ----------------------------------------------------------------------
! Takes a set of single-mode states, {|i>}, along one mode, and forms
!    sum_i[<i|V|i>] / size({|i>}).
! ----------------------------------------------------------------------
subroutine integrate_over_mode_average_PolynomialPotential(this,mode,states)
  use linear_algebra_module
  use eigenstates_module
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  integer,                    intent(in)    :: mode
  type(SingleModeState),      intent(in)    :: states(:)
  
  ! couplings(i+1,j+1) = <j|u^i|j>
  real(dp), allocatable :: couplings(:,:)
  
  integer :: i,ialloc
  
  ! Generate couplings between states along the mode.
  allocate( couplings(this%cutoff+1,size(states)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(states)
    couplings(:,i) = this%calculate_power_coupling(mode,states(i),states(i))
  enddo
  
  ! Integrate across each monomial.
  do i=1,size(this%monomials)
    ! (a) * u1^n1*...um^nm*... -> (avg_j[a*<j|um^nm|j>]) * u1^n1*...*um^0*...
    this%monomials(i)%coefficient = this%monomials(i)%coefficient &
             & * sum(couplings(this%monomials(i)%powers(mode),:)) &
             & / size(couplings,2)
    this%monomials(i)%powers(mode) = 0
  enddo
  
  ! Simplify across now identical monomials.
  call this%simplify()
end subroutine

! ----------------------------------------------------------------------
! Construct the Hamiltonian between states along a single mode.
! ----------------------------------------------------------------------
! Should only be called once the potential has been integrated across all
!    other modes.
! output(i,j) = <i-1|H|j-1>
function construct_hamiltonian_PolynomialPotential(this,mode,states) &
   & result(output)
  use eigenstates_module
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  integer,                    intent(in) :: mode
  type(SingleModeState),      intent(in) :: states(:)
  type(RealMatrix)                       :: output
  
  type(RealMatrix), allocatable :: couplings(:)
  
  integer :: i
  
  ! Generate couplings between states along the mode.
  couplings = this%generate_state_couplings(mode,states)
  
  ! Sum over (coefficient * coupling) at each power.
  output = dble(int(zeroes(size(states),size(states))))
  do i=1,size(this%monomials)
    output = output &
         & + this%monomials(i)%coefficient &
         & * couplings(this%monomials(i)%powers(mode))
  enddo
end function
end module
