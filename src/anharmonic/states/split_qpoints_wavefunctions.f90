! ======================================================================
! Wavefunctions which treat q-points separately.
! ======================================================================
module caesar_split_qpoints_wavefunctions_module
  use caesar_common_module
  
  use caesar_subspace_wavefunctions_module
  implicit none
  
  private
  
  public :: startup_split_qpoints_wavefunctions
  
  public :: SplitQpointsWavefunctions
  
  type, extends(SubspaceWavefunctions) :: SplitQpointsWavefunctions
    integer                   :: subspace_id
    integer,      allocatable :: mode_ids(:)
    integer,      allocatable :: paired_mode_ids(:)
    type(String)              :: harmonic_ground_state
    real(dp),     allocatable :: energies(:)
    type(String), allocatable :: wavefunctions(:)
  contains
    ! Type representation.
    procedure, public, nopass :: representation => &
                               & representation_SplitQpointsWavefunctions
    
    ! I/O.
    procedure, public :: read  => read_SplitQpointsWavefunctions
    procedure, public :: write => write_SplitQpointsWavefunctions
  end type
  
  interface SplitQpointsWavefunctions
    module procedure new_SplitQpointsWavefunctions
    module procedure new_SplitQpointsWavefunctions_Strings
    module procedure new_SplitQpointsWavefunctions_StringArray
  end interface
contains

subroutine startup_split_qpoints_wavefunctions()
  implicit none
  
  type(SplitQpointsWavefunctions) :: wavefunctions
  
  call wavefunctions%startup()
end subroutine

function new_SplitQpointsWavefunctions(subspace_id,mode_ids,paired_mode_ids, &
   & harmonic_ground_state,energies,wavefunctions) result(this)
  implicit none
    
  integer,      intent(in)        :: subspace_id
  integer,      intent(in)        :: mode_ids(:)
  integer,      intent(in)        :: paired_mode_ids(:)
  type(String), intent(in)        :: harmonic_ground_state
  real(dp),     intent(in)        :: energies(:)
  type(String), intent(in)        :: wavefunctions(:)
  type(SplitQpointsWavefunctions) :: this
  
  if (size(energies)/=size(wavefunctions)) then
    call err()
  endif
  
  this%subspace_id = subspace_id
  this%mode_ids = mode_ids
  this%paired_mode_ids = paired_mode_ids
  this%harmonic_ground_state = harmonic_ground_state
  this%energies = energies
  this%wavefunctions = wavefunctions
end function

impure elemental function representation_SplitQpointsWavefunctions() &
   & result(output)
  implicit none
  
  type(String) :: output
  
  output = 'split_qpoints'
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_SplitQpointsWavefunctions(this,input)
  implicit none
  
  class(SplitQpointsWavefunctions), intent(out) :: this
  type(String),                     intent(in)  :: input(:)
  
  integer                   :: subspace_id
  integer,      allocatable :: mode_ids(:)
  integer,      allocatable :: paired_mode_ids(:)
  type(String)              :: harmonic_ground_state
  real(dp),     allocatable :: energies(:)
  type(String), allocatable :: wavefunctions(:)
  
  type(String), allocatable :: line(:)
  
  integer :: no_wavefunctions
  
  integer :: i,ialloc
  
  select type(this); type is(SplitQpointsWavefunctions)
    line = split_line(input(1))
    subspace_id = int(line(3))
    
    line = split_line(input(2))
    mode_ids = int(line(4:))
    
    line = split_line(input(3))
    paired_mode_ids = int(line(5:))
    
    line = split_line(input(4))
    harmonic_ground_state = join(line(6:))
    
    no_wavefunctions = (size(input)-5)/4
    allocate( energies(no_wavefunctions),      &
            & wavefunctions(no_wavefunctions), &
            & stat=ialloc); call err(ialloc)
    do i=1,no_wavefunctions
      line = split_line(input(5+4*(i-1)+2))
      energies(i) = dble(line(4))
      
      line = split_line(input(5+4*(i-1)+4))
      wavefunctions(i) = join(line(3:))
    enddo
    
    this = SplitQpointsWavefunctions( subspace_id,           &
                                    & mode_ids,              &
                                    & paired_mode_ids,       &
                                    & harmonic_ground_state, &
                                    & energies,              &
                                    & wavefunctions          )
  class default
    call err()
  end select
end subroutine

function write_SplitQpointsWavefunctions(this) result(output)
  implicit none
  
  class(SplitQpointsWavefunctions), intent(in) :: this
  type(String), allocatable                    :: output(:)
  
  integer :: i
  
  select type(this); type is(SplitQpointsWavefunctions)
    output = [ 'Subspace                   : '//this%subspace_id,           &
             & 'Mode IDs                   : '//this%mode_ids,              &
             & 'Paired Mode Ids            : '//this%paired_mode_ids,       &
             & 'Harmonic ground state, |0> : '//this%harmonic_ground_state, &
             & str('Wavefunctions')                                         ]
    do i=1,size(this%wavefunctions)
      output = [ output,                                      &
               & str(''),                                     &
               & 'State energy     : '//this%energies(i),     &
               & 'Wavefunction     : '//this%wavefunctions(i) ]
    enddo
  
  
  class default
    call err()
  end select
end function

function new_SplitQpointsWavefunctions_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)        :: input(:)
  type(SplitQpointsWavefunctions) :: this
  
  call this%read(input)
end function

impure elemental function new_SplitQpointsWavefunctions_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in)   :: input
  type(SplitQpointsWavefunctions) :: this
  
  this = SplitQpointsWavefunctions(str(input))
end function
end module
