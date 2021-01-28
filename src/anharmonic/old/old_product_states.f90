! ======================================================================
! A product of VSCF states.
! ======================================================================
module caesar_product_states_module
  use caesar_common_module
  
  use caesar_vscf_states_module
  use caesar_single_mode_states_module
  use caesar_coupled_modes_module
  use caesar_grid_types_module
  implicit none
  
  type :: ProductStates
    type(VscfStates), allocatable :: vscf_states(:)
    
    ! A no_modes*no_states array of state indices.
    ! e.g. the row [0,2,0] referrs to the state |0>|2>|0>, with no_modes=3.
    integer, private, allocatable :: states_(:,:)
  contains
    procedure, public :: evaluate
    procedure, public :: single_mode_state
    procedure, public :: vscf_energy
    procedure, public :: kinetic_energy
    
    procedure, public :: state_as_bra_string
    procedure, public :: state_as_ket_string
  end type
  
  interface size
    module procedure size_ProductStates
  end interface
contains

! ----------------------------------------------------------------------
! Returns the number of product states contained.
! ----------------------------------------------------------------------
function size_ProductStates(this) result(output)
  implicit none
  
  class(ProductStates), intent(in) :: this
  integer                          :: output
  
  output = size(this%states_,2)
end function

! ----------------------------------------------------------------------
! Evaluates the state at a given displacement in normal mode co-ordinates.
! ----------------------------------------------------------------------
function evaluate(this,state,displacement) result(output)
  implicit none
  
  class(ProductStates), intent(in) :: this
  integer,              intent(in) :: state
  type(ModeVector),     intent(in) :: displacement
  real(dp)                         :: output
  
  type(SingleModeState) :: vscf_state
  
  integer :: i
  
  output = 1
  do i=1,size(this%states_,1)
    vscf_state = this%vscf_states(i)%state(this%states_(i,state))
    output = output * vscf_state%evaluate(displacement%vector(i))
  enddo
end function

! ----------------------------------------------------------------------
! If a product state is |i_1>|i_2>...|i_mode>...|i_no_modes>,
!    returns i_mode.
! ----------------------------------------------------------------------
function single_mode_state(this,state,mode) result(output)
  implicit none
  
  class(ProductStates), intent(in) :: this
  integer,              intent(in) :: state
  integer,              intent(in) :: mode
  integer                          :: output
  
  output = this%states_(mode,state)
end function

! ----------------------------------------------------------------------
! Returns the VSCF energy, <i1|V1|i1> + <i2|V2|i2> + ...
! ----------------------------------------------------------------------
function vscf_energy(this,state) result(output)
  implicit none
  
  class(ProductStates), intent(in) :: this
  integer,              intent(in) :: state
  real(dp)                         :: output
  
  integer :: mode
  integer :: i
  
  output = 0
  do mode=1,size(this%states_,1)
    i = this%single_mode_state(state,mode)
    output = output+this%vscf_states(mode)%vscf_energy(i)
  enddo
end function

! ----------------------------------------------------------------------
! Returns the kinetic energy, <bra|T|ket>.
! ----------------------------------------------------------------------
function kinetic_energy(this,bra,ket) result(output)
  implicit none
  
  class(ProductStates), intent(in) :: this
  integer,              intent(in) :: bra
  integer,              intent(in) :: ket
  real(dp)                         :: output
  
  integer :: mode
  integer :: bra_i,ket_i
  
  output = 0
  do mode=1,size(this%states_,1)
    bra_i = this%single_mode_state(bra,mode)
    ket_i = this%single_mode_state(ket,mode)
    output = output + this%vscf_states(mode)%kinetic_energy(bra_i,ket_i)
  enddo
end function

! ----------------------------------------------------------------------
! Prints the product state in the basis of single-mode states.
! ----------------------------------------------------------------------
function state_as_bra_string(this,state) result(output)
  implicit none
  
  class(ProductStates), intent(in) :: this
  integer,              intent(in) :: state
  type(String)                     :: output
  
  integer      :: i
  
  output = ''
  do i=1,size(this%states_,1)
    output = output//'<'//this%states_(i,state)//'|'
  enddo
end function

function state_as_ket_string(this,state) result(output)
  implicit none
  
  class(ProductStates), intent(in) :: this
  integer,              intent(in) :: state
  type(String)                     :: output
  
  integer      :: i
  
  output = ''
  do i=1,size(this%states_,1)
    output = output//'|'//this%states_(i,state)//'>'
  enddo
end function

! ----------------------------------------------------------------------
! Constructs all relevant product states from single-mode states.
! ----------------------------------------------------------------------
function construct_product_states(vscf_states,coupling) result(output)
  implicit none
  
  type(VscfStates),   intent(in) :: vscf_states(:)
  type(CoupledModes), intent(in) :: coupling(:)
  type(ProductStates)            :: output
  
  integer              :: no_modes
  integer              :: cutoff
  integer              :: output_size
  integer, allocatable :: grid(:,:)
  
  integer :: mode
  
  integer :: i,j,k,l,ialloc
  
  output%vscf_states = vscf_states
  
  no_modes = size(vscf_states)
  cutoff = vscf_states(1)%cutoff()
  do i=2,no_modes
    if (vscf_states(i)%cutoff()/=cutoff) then
      call err()
    endif
  enddo
  
  ! --------------------------------------------------
  ! Calculate how many states there will be.
  ! --------------------------------------------------
  output_size = 0
  do i=1,size(coupling)
    ! This gives the array length of 'grid' below.
    output_size = output_size &
              & + octahedral_grid_size( size(coupling(i)), &
              &                         cutoff-1,          &
              &                         include_negatives=.false.)
  enddo
  
  ! --------------------------------------------------
  ! Allocate states, and initialise each state to |0>|0>...|0>
  ! --------------------------------------------------
  allocate( output%states_(no_modes,output_size), &
          & stat=ialloc); call err(ialloc)
  output%states_ = 0
  
  ! --------------------------------------------------
  ! Calculate states, coupling by coupling.
  ! --------------------------------------------------
  j = 0
  do i=1,size(coupling)
    ! Generate an octahedral grid of points.
    ! Each poing will be mapped on to the states in each product.
    ! The -1 comes from not taking any |0> states, since they will be handled
    !    by subsidiary couplings.
    grid = generate_octahedral_grid( size(coupling(i)), &
                                   & cutoff-1,          &
                                   & include_negatives=.false.)
    do k=1,size(grid,2)
      ! Set all states in the coupling.
      do l=1,size(coupling(i))
        mode = coupling(i)%modes(l)
        ! The +1 is for the same reason as the -1 above.
        output%states_(mode,j+k) = grid(l,k)+1
      enddo
    enddo
    j = j+size(grid,2)
  enddo
  
  if (j/=output_size) then
    call err()
  endif
end function
end module
