submodule (caesar_split_qpoints_wavefunctions_module) caesar_split_qpoints_wavefunctions_submodule
  use caesar_states_module
contains

module procedure new_SplitQpointsWavefunctions
  if (size(energies)/=size(wavefunctions)) then
    call err()
  endif
  
  this%subspace_id = subspace_id
  this%mode_ids = mode_ids
  this%paired_mode_ids = paired_mode_ids
  this%harmonic_ground_state = harmonic_ground_state
  this%energies = energies
  this%wavefunctions = wavefunctions
end procedure

module procedure representation_SplitQpointsWavefunctions
  output = 'split_qpoints'
end procedure

module procedure read_SplitQpointsWavefunctions
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
end procedure

module procedure write_SplitQpointsWavefunctions
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
end procedure

module procedure new_SplitQpointsWavefunctions_Strings
  call this%read(input)
end procedure

module procedure new_SplitQpointsWavefunctions_StringArray
  this = SplitQpointsWavefunctions(str(input))
end procedure
end submodule
