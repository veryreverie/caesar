! ======================================================================
! Parses VSCF states to allow wavefunctions to be printed.
! ======================================================================
module vscf_wavefunctions_module
  use common_module
  
  use subspace_basis_module
  use vscf_state_module
  implicit none
  
  private
  
  public :: VscfWavefunctions
  
  type, extends(Stringsable) :: VscfWavefunctions
    integer                   :: subspace_id
    integer,      allocatable :: mode_ids(:)
    integer,      allocatable :: paired_mode_ids(:)
    type(String)              :: harmonic_ground_state
    real(dp),     allocatable :: energies(:)
    integer,      allocatable :: degeneracies(:)
    type(String), allocatable :: wavefunctions(:)
  contains
    procedure, public :: read  => read_VscfWavefunctions
    procedure, public :: write => write_VscfWavefunctions
  end type
  
  interface VscfWavefunctions
    module procedure new_VscfWavefunctions
    module procedure new_VscfWavefunctions_VscfStates
    module procedure new_VscfWavefunctions_Strings
    module procedure new_VscfWavefunctions_StringArray
  end interface
  
contains

function new_VscfWavefunctions(subspace_id,mode_ids,paired_mode_ids,  &
   & harmonic_ground_state,energies,degeneracies,wavefunctions) result(this)
  implicit none
    
  integer,      intent(in) :: subspace_id
  integer,      intent(in) :: mode_ids(:)
  integer,      intent(in) :: paired_mode_ids(:)
  type(String), intent(in) :: harmonic_ground_state
  real(dp),     intent(in) :: energies(:)
  integer,      intent(in) :: degeneracies(:)
  type(String), intent(in) :: wavefunctions(:)
  type(VscfWavefunctions)  :: this
  
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

function new_VscfWavefunctions_VscfStates(states,subspace,basis,supercell) &
   & result(this)
  implicit none
  
  type(VscfState),          intent(in) :: states(:)
  type(DegenerateSubspace), intent(in) :: subspace
  type(SubspaceBasis),      intent(in) :: basis
  type(StructureData),      intent(in) :: supercell
  type(VscfWavefunctions)              :: this
  
  type(String)              :: ground_state
  type(String), allocatable :: wavefunctions(:)
  
  ground_state = basis%ground_state_wavefunction(subspace,supercell)
  wavefunctions = states%wavefunction(basis,supercell)
  
  this = VscfWavefunctions( subspace%id,         &
                          & subspace%mode_ids,   &
                          & subspace%paired_ids, &
                          & ground_state,        &
                          & states%energy,       &
                          & states%degeneracy,   &
                          & wavefunctions        )
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_VscfWavefunctions(this,input)
  implicit none
  
  class(VscfWavefunctions), intent(out) :: this
  type(String),             intent(in)  :: input(:)
  
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
  
  select type(this); type is(VscfWavefunctions)
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
    
    this = VscfWavefunctions( subspace_id,           &
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

function write_VscfWavefunctions(this) result(output)
  implicit none
  
  class(VscfWavefunctions), intent(in) :: this
  type(String), allocatable            :: output(:)
  
  integer :: i
  
  select type(this); type is(VscfWavefunctions)
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

function new_VscfWavefunctions_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(VscfWavefunctions)  :: this
  
  call this%read(input)
end function

impure elemental function new_VscfWavefunctions_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(VscfWavefunctions)       :: this
  
  this = VscfWavefunctions(str(input))
end function
end module
