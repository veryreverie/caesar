submodule (caesar_harmonic_state_real_module) caesar_harmonic_state_real_submodule
  use caesar_states_module
contains

module procedure prod_real
  output = HarmonicStateReal([lhs%modes_,rhs%modes_])
end procedure

module procedure startup_harmonic_state_real
  type(HarmonicStateReal) :: state
  
  call state%startup()
end procedure

module procedure new_HarmonicStateReal
  this%modes_ = modes
end procedure

module procedure new_HarmonicStateReal_SubspaceState
  select type(input); type is(HarmonicStateReal)
    this = input
  type is(SubspaceStatePointer)
    ! WORKAROUND: ifort doesn't recognise the interface to this module procedure
    !    from within this module procedure, so the full name is used instead.
    this = new_HarmonicStateReal_SubspaceState(input%state())
  class default
    call err()
  end select
end procedure

module procedure harmonic_state_real_pointer
  select type(input); type is(HarmonicStateReal)
    this => input
  type is(SubspaceStatePointer)
    this => harmonic_state_real_pointer(input%state_pointer())
  class default
    call err()
  end select
end procedure

module procedure representation_HarmonicStateReal
  output = 'harmonic real'
end procedure

module procedure mode_ids_HarmonicStateReal
  output = this%modes_%id()
end procedure

module procedure paired_mode_ids_HarmonicStateReal
  output = this%modes_%id()
end procedure

module procedure occupation_HarmonicStateReal
  output = sum(this%modes_%total_occupation())
end procedure

module procedure wavevector_HarmonicStateReal
  integer :: i
  
  ! Effectively sum(this%modes_(:)%wavevector(modes,qpoints)).
  output = this%modes_(1)%wavevector(modes,qpoints)
  do i=2,size(this%modes_)
    output = output + this%modes_(i)%wavevector(modes,qpoints)
  enddo
end procedure

module procedure wavefunction_HarmonicStateReal
  ! TODO
  call err()
end procedure

module procedure change_modes_HarmonicStateReal
  integer, allocatable :: ids(:)
  integer, allocatable :: occupations(:)
  integer, allocatable :: sort_key(:)
  
  ! Get the ids and occupations of the single-mode terms.
  ids = this%modes_%id()
  occupations = this%modes_%occupation()
  
  ! Change the ids according to the given group.
  ids = mode_group*ids
  
  ! Sort the modes by id.
  sort_key = sort(ids)
  ids = ids(sort_key)
  occupations = occupations(sort_key)
  
  ! Construct output using the new ids.
  output = HarmonicStateReal(modes = HarmonicState1D(ids,occupations))
end procedure

module procedure read_HarmonicStateReal
  type(HarmonicState1D), allocatable :: modes(:)
  
  type(String), allocatable :: line(:)
  
  integer :: i
  
  select type(this); type is(HarmonicStateReal)
    line = split_line(input(1),delimiter='>')
    line = [(line(i)//'>',i=1,size(line))]
    modes = HarmonicState1D(line)
    
    this = HarmonicStateReal(modes)
  class default
    call err()
  end select
end procedure

module procedure write_HarmonicStateReal
  select type(this); type is(HarmonicStateReal)
    output = [join(str(this%modes_), delimiter='')]
  class default
    call err()
  end select
end procedure

module procedure new_HarmonicStateReal_Strings
  call this%read(input)
end procedure

module procedure new_HarmonicStateReal_StringArray
  this = HarmonicStateReal(str(input))
end procedure
end submodule
