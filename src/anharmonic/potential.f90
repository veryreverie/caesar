! ======================================================================
! The basis functions from which the Born-Oppenheimer potential is made.
! ======================================================================
module potential_module
  use constants_module, only : dp
  use string_module
  use stringable_module
  use io_module
  use linear_algebra_module
  implicit none
  
  private
  
  public :: calculate_potential
  public :: print_line
  
  ! A monomial, e.g.
  !    C * (u1)**a * (u2)**b * (u4)**d => coef=C, powers=[a,b,0,d]
  type, extends(Stringable) :: Monomial
    real(dp)             :: coefficient
    integer, allocatable :: powers(:)
  contains
    ! Evaluate energy or forces at a given value of u.
    procedure :: evaluate_energy => evaluate_energy_Monomial
    procedure :: evaluate_forces => evaluate_forces_Monomial
    
    ! I/O.
    procedure, pass(that) :: assign_String => assign_String_Monomial
  end type
  
  ! A polynomial representation of a potential.
  type, public :: PolynomialPotential
    integer,                       private :: cutoff
    type(Monomial),   allocatable, private :: monomials(:)
  contains
    procedure, public  :: evaluate_energy => &
                        & evaluate_energy_PolynomialPotential
    procedure, public  :: evaluate_forces => &
                        & evaluate_forces_PolynomialPotential
    procedure, public  :: simplify => simplify_PolynomialPotential
    generic,   public  :: integrate =>                                   &
                        & integrate_PolynomialPotential_HarmonicStates,  &
                        & integrate_PolynomialPotential_VscfStates,      &
                        & integrate_PolynomialPotential_ProductState
    procedure, private :: integrate_PolynomialPotential_HarmonicStates
    procedure, private :: integrate_PolynomialPotential_VscfStates
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
  
  ! I/O functions for potentials.
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
    if (derivative%powers(i)==0) then
      derivative%coefficient = 0
      derivative%coefficient = 0
    else
      derivative%coefficient = derivative%coefficient &
                           & * derivative%powers(i)
      derivative%powers(i) = derivative%powers(i) - 1
    endif
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
function calculate_potential(potential_basis_cutoff,sampling,energy_error, &
   & force_error,modes,qpoint,supercell) result(output)
  use sampling_points_module
  use normal_mode_module
  use qpoints_module
  use structure_module
  implicit none
  
  integer,                intent(in) :: potential_basis_cutoff
  type(CouplingSampling), intent(in) :: sampling(:)
  real(dp),               intent(in) :: energy_error
  real(dp),               intent(in) :: force_error
  type(NormalMode),       intent(in) :: modes(:)
  type(QpointData),       intent(in) :: qpoint
  type(StructureData),    intent(in) :: supercell
  type(PolynomialPotential)          :: output
  
  type(CouplingPotential), allocatable :: basis_functions(:)
  
  integer                                :: no_basis_functions
  integer,                   allocatable :: no_monomials(:)
  type(PolynomialPotential), allocatable :: fit(:)
  integer                                :: size_output
  
  integer :: i,j,k,l,ialloc
  
  output%cutoff = potential_basis_cutoff
  
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
                     & energy_error,                       &
                     & force_error,                        &
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
  ! Remove any constant terms (all powers=0), if present.
  size_output = 0
  do i=1,size(fit)
    size_output = size_output + size(fit(i)%monomials)
    do j=1,size(fit(i)%monomials)
      if (all(fit(i)%monomials(j)%powers==0)) then
        size_output = size_output-1
      endif
    enddo
  enddo
  allocate(output%monomials(size_output), stat=ialloc); call err(ialloc)
  k = 0
  do i=1,size(fit)
    do j=1,size(fit(i)%monomials)
      if (.not. all(fit(i)%monomials(j)%powers==0)) then
        k = k+1
        output%monomials(k) = fit(i)%monomials(j)
      endif
    enddo
  enddo
  if (k/=size_output) then
    call err()
  endif
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
function fit_basis_functions(basis_functions,this_sampling,energy_error, &
   & force_error,fit,prev_sampling,modes,qpoint,supercell) result(output)
  use sampling_points_module
  use normal_mode_module
  use qpoints_module
  use structure_module
  use grid_types_module
  implicit none
  
  type(PolynomialPotential), intent(in) :: basis_functions(:)
  type(CouplingSampling),    intent(in) :: this_sampling
  real(dp),                  intent(in) :: energy_error
  real(dp),                  intent(in) :: force_error
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
      
      a(l,j) = basis_energy / energy_error
      do k=1,size(this_sampling%coupling)
        a(l+k,j) = basis_forces%vector(this_sampling%coupling%modes(k)) &
               & / force_error
      enddo
    enddo
    
    b(l) = energy(i) / energy_error
    do k=1,size(this_sampling%coupling)
      b(l+k) = forces(i)%vector(this_sampling%coupling%modes(k)) &
           & / force_error
    enddo
  enddo
  
  ! Perform linear least-squares optimisation.
  ! Finds x s.t. (a.x-b)^2 is minimal.
  output = dble(linear_least_squares(a,b))
end function

! ----------------------------------------------------------------------
! Forms <bra|V|ket> between two harmonic states.
! ----------------------------------------------------------------------
!    a              * u1^p1 * ... * um^pm * ...
! -> a*<bra|um^nm|ket>) * u1^p1 * ... * um^0  * ...
subroutine integrate_PolynomialPotential_HarmonicStates(this,states,bra,ket)
  use harmonic_states_module
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  type(HarmonicStates),       intent(in)    :: states
  integer,                    intent(in)    :: bra
  integer,                    intent(in)    :: ket
  
  integer :: mode
  integer :: power
  
  integer :: i
  
  mode = states%mode
  
  do i=1,size(this%monomials)
    power = this%monomials(i)%powers(mode)
    this%monomials(i)%coefficient = this%monomials(bra)%coefficient &
                                & * states%integral(bra,ket,power)
    this%monomials(i)%powers(mode) = 0
  enddo
  
  ! Simplify across now identical monomials.
  call this%simplify()
end subroutine

! ----------------------------------------------------------------------
! Forms <bra|V|ket> between two VSCF states.
! ----------------------------------------------------------------------
!    a              * u1^p1 * ... * um^pm * ...
! -> a*<bra|um^nm|ket>) * u1^p1 * ... * um^0  * ...
subroutine integrate_PolynomialPotential_VscfStates(this,states,bra,ket)
  use vscf_states_module
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  type(VscfStates),           intent(in)    :: states
  integer,                    intent(in)    :: bra
  integer,                    intent(in)    :: ket
  
  integer :: mode
  integer :: power
  
  integer :: i
  
  mode = states%mode
  
  do i=1,size(this%monomials)
    power = this%monomials(i)%powers(mode)
    this%monomials(i)%coefficient = this%monomials(i)%coefficient &
                                & * states%integral(bra,ket,power)
    this%monomials(i)%powers(mode) = 0
  enddo
  
  ! Simplify across now identical monomials.
  call this%simplify()
end subroutine

! ----------------------------------------------------------------------
! Takes two product states, |bra> and |ket>, and forms <bra|V|ket>
! ----------------------------------------------------------------------
subroutine integrate_PolynomialPotential_ProductState(this,states,bra,ket)
  use product_states_module
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  type(ProductStates),        intent(in)    :: states
  integer,                    intent(in)    :: bra
  integer,                    intent(in)    :: ket
  
  integer :: mode
  integer :: i,j
  
  do mode=1,size(states%vscf_states)
    i = states%single_mode_state(bra,mode)
    j = states%single_mode_state(ket,mode)
    call this%integrate(states%vscf_states(mode),i,j)
  enddo
  
  if (size(this%monomials)/=1) then
    call err()
  endif
end subroutine

! ----------------------------------------------------------------------
! As above, but returns the coefficient rather than modifying the potential.
! ----------------------------------------------------------------------
function integrate_to_constant_PolynomialPotential(this,states,bra,ket) &
   & result(output)
  use product_states_module
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  type(ProductStates),        intent(in)    :: states
  integer,                    intent(in)    :: bra
  integer,                    intent(in)    :: ket
  real(dp)                                  :: output
  
  type(PolynomialPotential) :: temp
  
  temp = this
  
  call temp%integrate(states,bra,ket)
  
  output = temp%monomials(1)%coefficient
end function

! ----------------------------------------------------------------------
! Takes a set of single-mode states, {|i>}, along one mode, and forms
!    sum_i[<i|V|i>] / size({|i>}).
! ----------------------------------------------------------------------
subroutine integrate_over_mode_average_PolynomialPotential(this,vscf_states)
  use linear_algebra_module
  use vscf_states_module
  implicit none
  
  class(PolynomialPotential), intent(inout) :: this
  type(VscfStates),           intent(in)    :: vscf_states
  
  integer :: power
  
  integer :: i
  
  ! Integrate across each monomial.
  do i=1,size(this%monomials)
    !    a                      * u1^p1 * ... * ui^pi * ...
    ! -> a * avg_j[<j|ui^pi|j>] * u1^p1 * ... * ui^0  * ...
    power = this%monomials(i)%powers(vscf_states%mode)
    this%monomials(i)%coefficient = this%monomials(i)%coefficient &
                                & * vscf_states%mean_field(power)
    this%monomials(i)%powers(vscf_states%mode) = 0
  enddo
  
  ! Simplify across now identical monomials.
  call this%simplify()
end subroutine

! ----------------------------------------------------------------------
! Construct the Hamiltonian between basis states along a single mode.
! ----------------------------------------------------------------------
! Should only be called once the potential has been integrated across all
!    other modes.
! output(i+1,j+1) = <i|H|j> = <i|(T+V)|j>
function construct_hamiltonian_PolynomialPotential(this,states) &
   & result(output)
  use linear_algebra_module
  use harmonic_states_module
  implicit none
  
  class(PolynomialPotential), intent(in) :: this
  type(HarmonicStates),       intent(in) :: states
  type(RealMatrix)                       :: output
  
  real(dp) :: coefficient
  integer  :: power
  
  integer :: i
  
  ! Add <i|T|j>.
  output = states%kinetic_energies()
  
  ! Add <i|V|j> = sum(coefficient*<i|u^power|j>).
  do i=1,size(this%monomials)
    coefficient = this%monomials(i)%coefficient
    power = this%monomials(i)%powers(states%mode)
    output = output + coefficient*states%integrals(power)
  enddo
end function

! ----------------------------------------------------------------------
! I/O functions for potentials.
! ----------------------------------------------------------------------
subroutine assign_String_Monomial(this,that)
  implicit none
  
  type(String),    intent(inout) :: this
  class(Monomial), intent(in)    :: that
  
  integer :: i
  
  this = that%coefficient
  do i=1,size(that%powers)
    if (that%powers(i)/=0) then
      this = this//'*u'//i//'^'//that%powers(i)
    endif
  enddo
end subroutine

subroutine print_line_PolynomialPotential(this)
  implicit none
  
  type(PolynomialPotential), intent(in) :: this
  
  integer :: i
  
  do i=1,size(this%monomials)
    if (i==1) then
      call print_line('   '//this%monomials(i))
    else
      call print_line(' + '//this%monomials(i))
    endif
  enddo
end subroutine
end module
