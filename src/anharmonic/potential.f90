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
  public :: print_line
  
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
  
  ! An interim type for fitting a polynomial potential.
  type :: CouplingPotential
    type(PolynomialPotential), allocatable :: basis_functions(:)
    real(dp),                  allocatable :: coefficients(:)
  end type
  
  interface print_line
    module procedure print_line_PolynomialPotential
  end interface
contains

! ----------------------------------------------------------------------
! Evaluates the energy of a Monomial at a given displacement.
! ----------------------------------------------------------------------
function evaluate_energy_Monomial(this,displacement) result(output)
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
  
  type(CouplingPotential), allocatable :: basis_functions(:)
  
  integer                                :: no_basis_functions
  integer,                   allocatable :: no_monomials(:)
  type(PolynomialPotential), allocatable :: fit(:)
  integer                                :: size_output
  
  integer :: i,j,k,l,ialloc
  
  output%cutoff = potential_basis_cutoff
  
  ! Calculate gaussian integrals, {I(n)}.
  ! gaussian_integrals(n+1) = I(n) = integral[ u^n e^(-freq*u*u) ].
  call output%calculate_gaussian_integrals(modes, harmonic_states_cutoff)
  
  ! Generate potential at each coupling.
  allocate( basis_functions(size(sampling)), &
          & fit(size(sampling)),             &
          & stat=ialloc); call err(ialloc)
  do i=1,size(sampling)
    ! Generate basis functions.
    basis_functions(i)%basis_functions = calculate_basis_functions( &
                                            & sampling(i)%coupling, &
                                            & size(modes),          &
                                            & potential_basis_cutoff)
    
    ! Use linear least-squares regression to fit basis function coefficients
    !    to input energy and force data.
    basis_functions(i)%coefficients = fit_basis_functions( &
                     & basis_functions(i)%basis_functions, &
                     & sampling(i),                        &
                     & fit(:i-1),                          &
                     & sampling(:i-1),                     &
                     & modes,                              &
                     & qpoint,                             &
                     & supercell)
    
    ! Make a single potential containing all monomials in all basis functions,
    !    with coefficients as dictated by the fitting procedure.
    no_basis_functions = size(basis_functions(i)%basis_functions)
    
    allocate(no_monomials(no_basis_functions), stat=ialloc); call err(ialloc)
    do j=1,no_basis_functions
      no_monomials(j) = size(basis_functions(i)%basis_functions(j)%monomials)
    enddo
    
    allocate( fit(i)%monomials(sum(no_monomials)), &
            & stat=ialloc); call err(ialloc)
    l = 0
    do j=1,no_basis_functions
      do k=1,no_monomials(j)
        fit(i)%monomials(l+k) = &
           & basis_functions(i)%basis_functions(j)%monomials(k)
        fit(i)%monomials(l+k)%coefficient = fit(i)%monomials(l+k)%coefficient &
                                        & * basis_functions(i)%coefficients(j)
      enddo
      l = l+no_monomials(j)
    enddo
    
    deallocate(no_monomials,stat=ialloc); call err(ialloc)
    
    ! Simplify the potential, adding together identical monomials.
    call fit(i)%simplify()
  enddo
  
  ! Combine the potentials from each coupling together.
  size_output = 0
  do i=1,size(fit)
    size_output = size_output + size(fit(i)%monomials)
  enddo
  allocate(output%monomials(size_output), stat=ialloc); call err(ialloc)
  k = 0
  do i=1,size(fit)
    do j=1,size(fit(i)%monomials)
      output%monomials(k+j) = fit(i)%monomials(j)
    enddo
    k = k + size(fit(i)%monomials)
  enddo
end function

! ----------------------------------------------------------------------
! Generate all polynomial basis functions.
! At present, all basis functions are monomials with coefficient 1.
! ----------------------------------------------------------------------
! Basis functions are set up in a generalised octahedral manner, i.e.
!    sum(powers) is in set [0,potential_basis_cutoff]
! e.g. if potential_basis_cutoff=3 then
!    u1^3, u1^2*u2 and u1*u2*u3 are allowed, but u1^2*u2^2 is not.
! These are generated by first generating an octahedral grid of points, and
!    then mapping the indices of each point to the powers of each monomial.
function calculate_basis_functions(coupling,no_modes,potential_basis_cutoff) &
   & result(output)
  use coupling_module
  use sampling_points_module
  use grid_types_module
  implicit none
  
  type(CoupledModes), intent(in)         :: coupling
  integer,            intent(in)         :: no_modes
  integer,            intent(in)         :: potential_basis_cutoff
  type(PolynomialPotential), allocatable :: output(:)
  
  integer, allocatable :: grid(:,:)
  integer              :: mode
  
  integer :: i,j,ialloc
  
  ! Generate an octahedral grid of points.
  ! Each point will be mapped onto a single basis function.
  ! The -1 comes from not taking any u^0 terms, since they will be handled by
  !    subsidiary couplings.
  grid = generate_octahedral_grid( size(coupling), &
                                 & potential_basis_cutoff-1,   &
                                 & include_negatives=.false.)
  allocate(output(size(grid,2)),stat=ialloc); call err(ialloc)
  do i=1,size(output)
    ! Each basis function consists of only one monomial.
    allocate( output(i)%monomials(1), &
            & stat=ialloc); call err(ialloc)
    output(i)%monomials(1)%coefficient = 1
    output(i)%monomials(1)%powers = int(zeroes(no_modes))
    do j=1,size(coupling)
      ! The +1 here is for the same reason as the -1 above.
      mode = coupling%modes(j)
      output(i)%monomials(1)%powers(mode) = grid(j,i) + 1
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
! Only forces in the hyperplane of the coupling are considered.
function fit_basis_functions(basis_functions,this_sampling,fit,prev_sampling, &
   & modes,qpoint,supercell) result(output)
  use sampling_points_module
  use normal_mode_module
  use qpoints_module
  use structure_module
  use grid_types_module
  implicit none
  
  type(PolynomialPotential), intent(in) :: basis_functions(:)
  type(CouplingSampling),    intent(in) :: this_sampling
  type(PolynomialPotential), intent(in) :: fit(:)
  type(CouplingSampling),    intent(in) :: prev_sampling(:)
  type(NormalMode),          intent(in) :: modes(:)
  type(QpointData),          intent(in) :: qpoint
  type(StructureData),       intent(in) :: supercell
  real(dp), allocatable                 :: output(:)
  
  ! Independent sampling points.
  integer                          :: no_sampling_points
  type(SamplingPoint), allocatable :: sampling_points(:)
  real(dp),            allocatable :: energy(:)
  type(ModeVector),    allocatable :: forces(:)
  
  ! The normal-mode displacement at each sampling point.
  type(ModeVector) :: displacement
  
  ! The basis-function force at each sampling point.
  real(dp)         :: basis_energy
  type(ModeVector) :: basis_forces
  
  ! Linear least-squares matrix and vector, L=(a.x-b)**2.
  integer               :: m,n
  real(dp), allocatable :: a(:,:)
  real(dp), allocatable :: b(:)
  
  ! Temporary variables.
  integer :: i,j,k,l,ialloc
  
  ! Concatenate all independent sampling points.
  no_sampling_points = 0
  do i=1,size(this_sampling%sampling_points)
    if (.not. this_sampling%sampling_points(i)%duplicate) then
      no_sampling_points = no_sampling_points + 1
    endif
  enddo
  j = 0
  allocate( sampling_points(no_sampling_points),  &
          & energy(no_sampling_points),           &
          & forces(no_sampling_points),           &
          & stat=ialloc); call err(ialloc)
  do i=1,size(this_sampling%sampling_points)
    if (.not. this_sampling%sampling_points(i)%duplicate) then
      j = j+1
      
      ! Read in energy and forces.
      sampling_points(j) = this_sampling%sampling_points(i)
      energy(j) = this_sampling%energy(i)
      forces(j) = cartesian_to_normal_mode( this_sampling%forces(:,i), &
                                          & modes,                     &
                                          & qpoint,                    &
                                          & supercell )
    endif
  enddo
  
  ! Subtract off the contributions from subsidiary couplings.
  do i=1,size(sampling_points)
    do j=1,size(fit)
      if ( prev_sampling(j)%coupling%is_subsidiary_of( &
         & this_sampling%coupling)) then
        energy(i) = energy(i) &
                & - fit(j)%evaluate_energy(sampling_points(i)%displacement)
        forces(i) = forces(i) &
                & - fit(j)%evaluate_forces(sampling_points(i)%displacement)
      endif
    enddo
  enddo
  
  ! Construct a and b.
  m = no_sampling_points*(1+size(this_sampling%coupling))
  n = size(basis_functions)
  allocate( a(m,n), &
          & b(m),   &
          & stat=ialloc); call err(ialloc)
  do i=1,no_sampling_points
    l = (i-1)*(1+size(this_sampling%coupling)) + 1
    displacement = sampling_points(i)%displacement
    
    do j=1,n
      basis_energy = basis_functions(j)%evaluate_energy(displacement)
      basis_forces = basis_functions(j)%evaluate_forces(displacement)
      
      a(l,j) = basis_energy
      do k=1,size(this_sampling%coupling)
        a(l+k,j) = basis_forces%vector(this_sampling%coupling%modes(k))
      enddo
    enddo
    
    b(l) = energy(i)
    do k=1,size(this_sampling%coupling)
      b(l+k) = forces(i)%vector(this_sampling%coupling%modes(k))
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
  output = 0
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
                                & * couplings(this%monomials(i)%powers(mode)+1)
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
             & * sum(couplings(this%monomials(i)%powers(mode)+1,:)) &
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
    output = output                        &
         & + this%monomials(i)%coefficient &
         & * couplings(this%monomials(i)%powers(mode)+1)
  enddo
end function

! ----------------------------------------------------------------------
! Print functions for potentials.
! ----------------------------------------------------------------------
subroutine print_line_PolynomialPotential(this)
  implicit none
  
  type(PolynomialPotential), intent(in) :: this
  
  type(String) :: line
  
  integer :: i,j
  
  do i=1,size(this%monomials)
    if (i==1) then
      line = '   '
    else
      line = ' + '
    endif
    line = line // this%monomials(i)%coefficient
    do j=1,size(this%monomials(i)%powers)
      if (this%monomials(i)%powers(j)/=0) then
        line = line//'*u'//j//'^'//this%monomials(i)%powers(j)
      endif
    enddo
    call print_line(line)
  enddo
end subroutine
end module
