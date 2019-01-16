! ======================================================================
! Wavefunctions spanning the full subspace.
! ======================================================================
module full_subspace_wavefunctions_module
  use common_module
  
  use subspace_wavefunctions_module
  implicit none
  
  private
  
  public :: FullSubspaceWavefunctions
  
  type, extends(SubspaceWavefunctions) :: FullSubspaceWavefunctions
    integer                   :: subspace_id
    integer,      allocatable :: mode_ids(:)
    integer,      allocatable :: paired_mode_ids(:)
    type(String)              :: harmonic_ground_state
    real(dp),     allocatable :: energies(:)
    integer,      allocatable :: degeneracies(:)
    type(String), allocatable :: wavefunctions(:)
  contains
    procedure, public :: read  => read_FullSubspaceWavefunctions
    procedure, public :: write => write_FullSubspaceWavefunctions
  end type
  
  interface FullSubspaceWavefunctions
    module procedure new_FullSubspaceWavefunctions
    module procedure new_FullSubspaceWavefunctions_Strings
    module procedure new_FullSubspaceWavefunctions_StringArray
  end interface
  
contains

function new_FullSubspaceWavefunctions(subspace_id,mode_ids,paired_mode_ids, &
   & harmonic_ground_state,energies,degeneracies,wavefunctions) result(this)
  implicit none
    
  integer,      intent(in)        :: subspace_id
  integer,      intent(in)        :: mode_ids(:)
  integer,      intent(in)        :: paired_mode_ids(:)
  type(String), intent(in)        :: harmonic_ground_state
  real(dp),     intent(in)        :: energies(:)
  integer,      intent(in)        :: degeneracies(:)
  type(String), intent(in)        :: wavefunctions(:)
  type(FullSubspaceWavefunctions) :: this
  
  if (size(energies)/=size(degeneracies)) then
    call err()
  elseif (size(energies)/=size(wavefunctions)) then
    call err()
  endif
  
  this%subspace_id = subspace_id
  this%mode_ids = mode_ids
  this%paired_mode_ids = paired_mode_ids
  this%harmonic_ground_state = harmonic_ground_state
  this%energies = energies
  this%degeneracies = degeneracies
  this%wavefunctions = wavefunctions
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_FullSubspaceWavefunctions(this,input)
  implicit none
  
  class(FullSubspaceWavefunctions), intent(out) :: this
  type(String),                     intent(in)  :: input(:)
  
  integer                   :: subspace_id
  integer,      allocatable :: mode_ids(:)
  integer,      allocatable :: paired_mode_ids(:)
  type(String)              :: harmonic_ground_state
  real(dp),     allocatable :: energies(:)
  integer,      allocatable :: degeneracies(:)
  type(String), allocatable :: wavefunctions(:)
  
  type(String), allocatable :: line(:)
  
  integer :: no_wavefunctions
  
  integer :: i,ialloc
  
  select type(this); type is(FullSubspaceWavefunctions)
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
            & degeneracies(no_wavefunctions),  &
            & wavefunctions(no_wavefunctions), &
            & stat=ialloc); call err(ialloc)
    do i=1,no_wavefunctions
      line = split_line(input(5+4*(i-1)+2))
      energies(i) = dble(line(4))
      
      line = split_line(input(5+4*(i-1)+3))
      degeneracies(i) = int(line(4))
      
      line = split_line(input(5+4*(i-1)+4))
      wavefunctions(i) = join(line(3:))
    enddo
    
    this = FullSubspaceWavefunctions( subspace_id,           &
                                    & mode_ids,              &
                                    & paired_mode_ids,       &
                                    & harmonic_ground_state, &
                                    & energies,              &
                                    & degeneracies,          &
                                    & wavefunctions          )
  class default
    call err()
  end select
end subroutine

function write_FullSubspaceWavefunctions(this) result(output)
  implicit none
  
  class(FullSubspaceWavefunctions), intent(in) :: this
  type(String), allocatable                    :: output(:)
  
  integer :: i
  
  select type(this); type is(FullSubspaceWavefunctions)
    output = [ 'Subspace                   : '//this%subspace_id,           &
             & 'Mode IDs                   : '//this%mode_ids,              &
             & 'Paired Mode Ids            : '//this%paired_mode_ids,       &
             & 'Harmonic ground state, |0> : '//this%harmonic_ground_state, &
             & str('Wavefunctions')                                         ]
    do i=1,size(this%wavefunctions)
      output = [ output,                                      &
               & str(''),                                     &
               & 'State energy     : '//this%energies(i),     &
               & 'State degeneracy : '//this%degeneracies(i), &
               & 'Wavefunction     : '//this%wavefunctions(i) ]
    enddo
  
  
  class default
    call err()
  end select
end function

function new_FullSubspaceWavefunctions_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)        :: input(:)
  type(FullSubspaceWavefunctions) :: this
  
  call this%read(input)
end function

impure elemental function new_FullSubspaceWavefunctions_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in)   :: input
  type(FullSubspaceWavefunctions) :: this
  
  this = FullSubspaceWavefunctions(str(input))
end function
end module