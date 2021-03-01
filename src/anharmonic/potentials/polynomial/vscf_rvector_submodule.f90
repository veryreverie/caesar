submodule (caesar_vscf_rvector_module) caesar_vscf_rvector_submodule
  use caesar_polynomial_module
contains

module procedure new_VscfRvector
  this%subspace_id = subspace_id
  this%rvector     = rvector
end procedure

module procedure read_VscfRvector
  integer         :: subspace_id
  type(IntVector) :: rvector
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(VscfRvector)
    line = split_line(input)
    subspace_id = int(line(2))
    rvector = vec(int(line(4:6)))
    
    this = VscfRvector(subspace_id,rvector)
  end select
end procedure

module procedure write_VscfRvector
  select type(this); type is(VscfRvector)
    output = 'Subspace: '//this%subspace_id//' R-vector: '//this%rvector
  end select
end procedure

module procedure new_VscfRvector_String
  call this%read(input)
end procedure
end submodule
