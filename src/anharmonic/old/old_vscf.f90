module caesar_vscf_states_module
  use caesar_common_module
  
  use caesar_single_mode_states_module
  use caesar_harmonic_states_module
  implicit none
  
  private
  
  public :: construct_harmonic_vscf_states
  
  ! VSCF states, {|i>}, and their relation to the basis states.
  type, public :: VscfStates
    integer :: mode
    
    ! states_(i+1) = |i>.
    ! The +1 is so states_(1) = |0>.
    type(SingleModeState), private, allocatable :: states_(:)
    
    ! vscf_energies_(i+1) = <i|V|i>, where V is the VSCF potential.
    ! The +1 is so vscf_energies_(1) = <0|V|0>.
    real(dp), private, allocatable :: vscf_energies_(:)
    
    ! basis_(i+1,j+1) = <i|j>, where
    !    |i> is a VSCF state.
    !    |j> is a basis state.
    ! The +1s are so basis_(1,1) = <0|0>.
    real(dp), private, allocatable :: basis_(:,:)
    
    ! The basis states.
    type(HarmonicStates), private :: basis_states_
  contains
    procedure, public :: state
    procedure, public :: vscf_energy
    procedure, public :: kinetic_energy
    procedure, public :: cutoff
    procedure, public :: integral
    procedure, public :: mean_field
    procedure, public :: update_hamiltonian
  end type
contains

! ----------------------------------------------------------------------
! Returns |i>
! ----------------------------------------------------------------------
function state(this,i) result(output)
  implicit none
  
  class(VscfStates), intent(in) :: this
  integer,           intent(in) :: i
  type(SingleModeState)         :: output
  
  output = this%states_(i+1)
end function

! ----------------------------------------------------------------------
! Returns <i|H|i>. For the VSCF Hamiltonian.
! ----------------------------------------------------------------------
function vscf_energy(this,i) result(output)
  implicit none
  
  class(VscfStates), intent(in) :: this
  integer,           intent(in) :: i
  real(dp)                      :: output
  
  output = this%vscf_energies_(i+1)
end function

! ----------------------------------------------------------------------
! Returns <i|T|j>.
! ----------------------------------------------------------------------
function kinetic_energy(this,i,j) result(output)
  implicit none
  
  class(VscfStates), intent(in) :: this
  integer,           intent(in) :: i
  integer,           intent(in) :: j
  real(dp)                      :: output
  
  output = vec(this%basis_(i+1,:))               &
       & * this%basis_states_%kinetic_energies() &
       & * vec(this%basis_(j+1,:))
end function

! ----------------------------------------------------------------------
! Returns the number of states-1.
! ----------------------------------------------------------------------
function cutoff(this) result(output)
  implicit none
  
  class(VscfStates), intent(in) :: this
  integer                       :: output
  
  output = size(this%states_)-1
end function

! ----------------------------------------------------------------------
! Returns <i|u^power|j>.
! ----------------------------------------------------------------------
function integral(this,i,j,power) result(output)
  implicit none
  
  class(VscfStates), intent(in) :: this
  integer,           intent(in) :: i
  integer,           intent(in) :: j
  integer,           intent(in) :: power
  real(dp)                      :: output
  
  output = vec(this%basis_(i+1,:))             &
       & * this%basis_states_%integrals(power) &
       & * vec(this%basis_(:,j+1))
end function

! ----------------------------------------------------------------------
! Returns sum_i <i|u^power|i> / no_states.
! ----------------------------------------------------------------------
function mean_field(this,power) result(output)
  implicit none
  
  class(VscfStates), intent(in) :: this
  integer,           intent(in) :: power
  real(dp)                      :: output
  
  output = trace( mat(this%basis_)                    &
              & * this%basis_states_%integrals(power) &
              & * mat(transpose(this%basis_)))
  output = output/(this%cutoff()+1)
end function

! ----------------------------------------------------------------------
! Takes a Hamiltonian, and constructs the VSCF states which diagonalise it.
! ----------------------------------------------------------------------
subroutine update_hamiltonian(this,hamiltonian)
  implicit none
  
  class(VscfStates), intent(inout) :: this
  type(RealMatrix),  intent(in)    :: hamiltonian
  
  type(SymmetricEigenstuff) :: eigenstuff
  
  integer :: i,j
  
  eigenstuff = diagonalise_symmetric(hamiltonian)
  
  this%vscf_energies_ = eigenstuff%evals
  this%basis_ = eigenstuff%evecs
  
  do i=0,this%cutoff()
    this%states_(i+1) = this%state(i) * 0.0_dp
    do j=0,this%cutoff()
      this%states_(i+1) = this%state(i)        &
                      & + this%basis_(i+1,j+1) &
                      & * this%basis_states_%state(j)
    enddo
  enddo
end subroutine

! ----------------------------------------------------------------------
! Constructs a first guess from which the vscf states can be found.
! Each vscf state is initially just the corresponding harmonic state.
! ----------------------------------------------------------------------
function construct_harmonic_vscf_states(harmonic_states) result(output)
  implicit none
  
  type(HarmonicStates), intent(in) :: harmonic_states
  type(VscfStates)                 :: output
  
  integer :: i,ialloc
  
  output%mode = harmonic_states%mode
  allocate( output%states_(harmonic_states%cutoff()+1),        &
          & output%vscf_energies_(harmonic_states%cutoff()+1), &
          & stat=ialloc); call err(ialloc)
  do i=0,harmonic_states%cutoff()
    output%states_(i+1) = harmonic_states%state(i)
    output%vscf_energies_(i+1) = (i+0.5_dp)*output%states_(i+1)%frequency
  enddo
  output%basis_ = dble(make_identity_matrix(harmonic_states%cutoff()+1))
  output%basis_states_ = harmonic_states
end function
end module
