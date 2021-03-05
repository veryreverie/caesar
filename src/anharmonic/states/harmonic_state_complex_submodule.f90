submodule (caesar_harmonic_state_complex_module) caesar_harmonic_state_complex_submodule
  use caesar_states_module
contains

module procedure prod_complex
  output = HarmonicStateComplex([lhs%modes_,rhs%modes_])
end procedure

module procedure new_HarmonicStateComplex
  this%modes_ = modes
end procedure

module procedure new_HarmonicStateComplex_SubspaceState
  select type(input); type is(HarmonicStateComplex)
    this = input
  type is(SubspaceStatePointer)
    ! WORKAROUND: ifort doesn't recognise the interface to this module procedure
    !    from within this module procedure, so the full name is used instead.
    this = new_HarmonicStateComplex_SubspaceState(input%state())
  class default
    call err()
  end select
end procedure

module procedure harmonic_state_complex_pointer
  select type(input); type is(HarmonicStateComplex)
    this => input
  type is(SubspaceStatePointer)
    this => harmonic_state_complex_pointer(input%state_pointer())
  class default
    call err()
  end select
end procedure

module procedure representation_HarmonicStateComplex
  output = 'harmonic complex'
end procedure

module procedure mode_ids_HarmonicStateComplex
  output = this%modes_%id()
end procedure

module procedure paired_mode_ids_HarmonicStateComplex
  output = this%modes_%paired_id()
end procedure

module procedure occupation_HarmonicStateComplex
  output = sum(this%modes_%total_occupation())
end procedure

module procedure wavevector_HarmonicStateComplex
  integer :: i
  
  ! Effectively sum(this%modes_(:)%wavevector(modes,qpoints)).
  output = this%modes_(1)%wavevector(modes,qpoints)
  do i=2,size(this%modes_)
    output = output + this%modes_(i)%wavevector(modes,qpoints)
  enddo
end procedure

module procedure wavefunction_HarmonicStateComplex
  ! TODO
  call err()
end procedure

module procedure change_modes_HarmonicStateComplex
  integer, allocatable :: ids(:)
  integer, allocatable :: paired_ids(:)
  integer, allocatable :: occupations(:)
  integer, allocatable :: paired_occupations(:)
  integer, allocatable :: sort_key(:)
  
  ! Get the ids and occupations of the single-mode terms.
  ids = this%modes_%id()
  paired_ids = this%modes_%paired_id()
  occupations = this%modes_%occupation()
  paired_occupations = this%modes_%paired_occupation()
  
  ! Change the ids according to the given group.
  ids = mode_group*ids
  
  ! Sort the modes by id.
  sort_key = sort(ids)
  ids = ids(sort_key)
  paired_ids = paired_ids(sort_key)
  occupations = occupations(sort_key)
  paired_occupations = paired_occupations(sort_key)
  
  ! Construct output using the new ids.
  output = HarmonicStateComplex(                                         &
     & modes = HarmonicState2D( id                = ids,                 &
     &                          paired_id         = paired_ids,          &
     &                          occupation        = occupations,         &
     &                          paired_occupation = paired_occupations ) )
end procedure

module procedure read_HarmonicStateComplex
  type(HarmonicState2D), allocatable :: modes(:)
  
  type(String), allocatable :: line(:)
  
  integer :: i
  
  select type(this); type is(HarmonicStateComplex)
    line = split_line(input(1),delimiter='>')
    line = [(line(i)//'>',i=1,size(line))]
    modes = HarmonicState2D(line)
    
    this = HarmonicStateComplex(modes)
  class default
    call err()
  end select
end procedure

module procedure write_HarmonicStateComplex
  select type(this); type is(HarmonicStateComplex)
    output = [join(str(this%modes_), delimiter='')]
  class default
    call err()
  end select
end procedure

module procedure new_HarmonicStateComplex_Strings
  call this%read(input)
end procedure

module procedure new_HarmonicStateComplex_StringArray
  this = HarmonicStateComplex(str(input))
end procedure
end submodule
