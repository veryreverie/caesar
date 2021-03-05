submodule (caesar_subspace_wavefunctions_module) caesar_subspace_wavefunctions_submodule
  use caesar_anharmonic_common_module
contains

module procedure new_SubspaceWavefunctionsPointer
  integer :: ialloc
  
  select type(wavefunctions); type is (SubspaceWavefunctionsPointer)
    this = wavefunctions
  class default
    this%representation_ = wavefunctions%representation()
    allocate( this%wavefunctions_, source=wavefunctions, &
            & stat=ialloc); call err(ialloc)
  end select
end procedure

module procedure check_SubspaceWavefunctionsPointer
  if (.not. allocated(this%wavefunctions_)) then
    call print_line(CODE_ERROR//': Trying to use a &
       &SubspaceWavefunctionsPointer before it has been allocated.')
    call err()
  endif
end procedure

module procedure representation_SubspaceWavefunctionsPointer
  output = 'pointer'
end procedure

module procedure read_SubspaceWavefunctionsPointer
  type(SubspaceWavefunctionsPointer), allocatable :: types(:)
  
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  
  integer :: i
  
  select type(this); type is(SubspaceWavefunctionsPointer)
    line = split_line(input(1))
    representation = line(3)
    
    ! Identify which type corresponds to the representation.
    types = types_SubspaceWavefunctions()
    i = first([( types(i)%wavefunctions_%representation()==representation, &
               & i=1,                                                      &
               & size(types)                                               )])
    this = types(i)
    
    call this%wavefunctions_%read(input(2:))
  class default
    call err()
  end select
end procedure

module procedure write_SubspaceWavefunctionsPointer
  select type(this); type is(SubspaceWavefunctionsPointer)
    output = [ 'Wavefunction representation: '//this%representation_, &
             & str(this%wavefunctions_)                               ]
  end select
end procedure

module procedure new_SubspaceWavefunctionsPointer_Strings
  call this%read(input)
end procedure

module procedure new_SubspaceWavefunctionsPointer_StringArray
  this = SubspaceWavefunctionsPointer(str(input))
end procedure
end submodule
