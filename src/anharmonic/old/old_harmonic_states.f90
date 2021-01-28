! ======================================================================
! Basis states, {|i>}, and integrals, {<i|u^p|j>}.
! ======================================================================
module caesar_harmonic_states_module
  use caesar_common_module
  
  use caesar_single_mode_states_module
  implicit none
  
  private
  
  public :: calculate_harmonic_states
  public :: size

  type, public :: HarmonicStates
    integer :: mode
    ! N.B. All arrays are stored such that array(1) corresponds to |0>, so
    !    array(i+1) corresponds to |i>.
    
    ! states_(i+1) = |i>.
    type(SingleModeState), private, allocatable :: states_(:)
    
    ! integrals_(i+1,j+1,p+1) = <i|u^p|j>.
    real(dp), private, allocatable :: integrals_(:,:,:)
    
    ! kinetic_energy(i+1,j+1) = <i|T|j> = <i|-0.5(d2/du2)|j>.
    real(dp), private, allocatable :: kinetic_energies_(:,:)
  contains
    ! Getters for states, integrals and energy, taking +1s into account.
    procedure, public :: state
    procedure, public :: integral
    procedure, public :: integrals
    procedure, public :: kinetic_energy
    procedure, public :: kinetic_energies
    
    ! Returns the cutoff, = size(states_)-1
    procedure, public :: cutoff
  end type
contains

! ----------------------------------------------------------------------
! Returns |i> along the specified mode.
! ----------------------------------------------------------------------
function state(this,i) result(output)
  class(HarmonicStates), intent(in) :: this
  integer,               intent(in) :: i
  type(SingleModeState)             :: output
  
  output = this%states_(i+1)
end function

! ----------------------------------------------------------------------
! Returns <i|u^power|j>.
! ----------------------------------------------------------------------
function integral(this,i,j,power) result(output)
  implicit none
  
  class(HarmonicStates), intent(in) :: this
  integer,               intent(in) :: i
  integer,               intent(in) :: j
  integer,               intent(in) :: power
  real(dp)                          :: output
  
  output = this%integrals_(i+1,j+1,power+1)
end function

! ----------------------------------------------------------------------
! Returns the matrix of <i|u^power|j>.
! ----------------------------------------------------------------------
function integrals(this,power) result(output)
  implicit none
  
  class(HarmonicStates), intent(in) :: this
  integer,               intent(in) :: power
  type(RealMatrix)                  :: output
  
  output = this%integrals_(:,:,power+1)
end function

! ----------------------------------------------------------------------
! Returns <i|T|j>.
! ----------------------------------------------------------------------
function kinetic_energy(this,i,j) result(output)
  implicit none
  
  class(HarmonicStates), intent(in) :: this
  integer,               intent(in) :: i
  integer,               intent(in) :: j
  real(dp)                          :: output
  
  output = this%kinetic_energies_(i+1,j+1)
end function

! ----------------------------------------------------------------------
! Returns the matrix of <i|T|j>.
! ----------------------------------------------------------------------
function kinetic_energies(this) result(output)
  implicit none
  
  class(HarmonicStates), intent(in) :: this
  type(RealMatrix)                  :: output
  
  output = this%kinetic_energies_
end function

! ----------------------------------------------------------------------
! Returns the number of states-1.
! ----------------------------------------------------------------------
function cutoff(this) result(output)
  implicit none
  
  class(HarmonicStates), intent(in) :: this
  integer                           :: output
  
  output = size(this%states_)-1
end function

! ----------------------------------------------------------------------
! Calculates states {|i>} and integrals {<i|u^p|j>} along all modes.
! ----------------------------------------------------------------------
function calculate_harmonic_states(mode,frequency,state_cutoff, &
   & potential_cutoff) result(output)
  implicit none
  
  integer,          intent(in) :: mode
  real(dp),         intent(in) :: frequency
  integer,          intent(in) :: state_cutoff
  integer,          intent(in) :: potential_cutoff
  type(HarmonicStates)         :: output
  
  integer :: i,j,p,ialloc
  
  ! --------------------------------------------------
  ! Copy across mode id.
  ! --------------------------------------------------
  output%mode = mode
  
  ! --------------------------------------------------
  ! Calculate states.
  ! --------------------------------------------------
  output%states_ = generate_harmonic_basis(mode,frequency,state_cutoff)
  
  ! --------------------------------------------------
  ! Calculate integrals.
  ! --------------------------------------------------
  
  allocate( output%integrals_( state_cutoff+1,      &
          &                    state_cutoff+1,      &
          &                    potential_cutoff+1), &
          & stat=ialloc); call err(ialloc)
  output%integrals_ = 0
  
  ! <i|u^0|j> = 1 if j=k
  !             0 otherwise
  do i=0,state_cutoff
    output%integrals_(i+1,i+1,1) = 1
  enddo
  
  ! --------------------------------------------------
  ! Fill out all <0|u^p|0>
  ! --------------------------------------------------
  ! <0|u^1|0> = 0
  ! <0|u^p|0> = (p-1)/(2*frequency) * <0|u^(p-2)|0>
  do p=2,potential_cutoff
    if (mod(p,2)==0) then
      output%integrals_(1,1,p+1) = ((p-1)/(2*frequency)) &
                               & * output%integral(0,0,p-2)
    endif
  enddo
  
  do p=1,potential_cutoff
    ! --------------------------------------------------
    ! Fill out all <i|u^p|0>.
    ! --------------------------------------------------
    ! <i|u^p|0> = 0                                        if i>p
    !             p/sqrt(2*i*frequency) <i-1|u^(p-1)|0>    otherwise
    do i=1,state_cutoff
      if (i<=p) then
        output%integrals_(i+1,1,p+1) = (p/sqrt(2*i*frequency)) &
                                   & * output%integral(i-1,0,p-1)
      endif
    enddo
    
    ! --------------------------------------------------
    ! Fill out all <0|u^p|j>.
    ! --------------------------------------------------
    ! <0|u^p|j> = <j|u^p|0>
    output%integrals_(1,:,p+1) = output%integrals_(:,1,p+1)
    
    ! --------------------------------------------------
    ! Fill out all <i|u^1|j>.
    ! --------------------------------------------------
    ! <i|u^1|j> = sqrt(i/2*frequency) <i-1|u^0|j  >
    !           + sqrt(j/2*frequency) <i  |u^0|j-1>
    if (p==1) then
      do i=1,state_cutoff
        do j=1,i
          output%integrals_(i+1,j+1,p+1) =                            &
             &   sqrt(i/(2*frequency)) * output%integral(i-1,j  ,p-1) &
             & + sqrt(j/(2*frequency)) * output%integral(i  ,j-1,p-1)
          
          ! <i|u^1|j> = <j|u^1|i>
          output%integrals_(j+1,i+1,p+1) = output%integrals_(i+1,j+1,p+1)
        enddo
      enddo
    endif
    
    ! --------------------------------------------------
    ! Fill out all <i|u^1|j>.
    ! --------------------------------------------------
    ! <i|u^p|j> = sqrt(i/2*frequency) <i-1|u^(p-1)|j  >
    !           + sqrt(j/2*frequency) <i  |u^(p-1)|j-1>
    !           + (p-1)/(2*frequency) <i  |u^(p-2)|j  >
    if (p>1) then
      do i=1,state_cutoff
        do j=1,i
          output%integrals_(i+1,j+1,p+1) =                            &
             &   sqrt(i/(2*frequency)) * output%integral(i-1,j  ,p-1) &
             & + sqrt(j/(2*frequency)) * output%integral(i  ,j-1,p-1) &
             & + ((p-1)/(2*frequency)) * output%integral(i  ,j  ,p-2)
          
          ! <i|u^1|j> = <j|u^1|i>
          output%integrals_(j+1,i+1,p+1) = output%integrals_(i+1,j+1,p+1)
        enddo
      enddo
    endif
  enddo
  
  ! --------------------------------------------------
  ! Calculate kinetic energies.
  ! --------------------------------------------------
  ! T|j> = -0.5(d2/du2)|j>
  !      = ( frequency*(j+0.5) - 0.5*(frequency*u)^2 ) |j>
  output%kinetic_energies_ = dble(zeroes(state_cutoff+1,state_cutoff+1))
  
  ! Add <i|frequency*(j+0.5)|j>.
  do i=0,state_cutoff
    output%kinetic_energies_(i+1,i+1) = frequency*(i+0.5_dp)
  enddo
  
  ! Add 0.5*frequency^2*<i|u^2|j>.
  output%kinetic_energies_ = output%kinetic_energies_ &
                         & + 0.5_dp*frequency**2*output%integrals_(:,:,3)
end function
end module
