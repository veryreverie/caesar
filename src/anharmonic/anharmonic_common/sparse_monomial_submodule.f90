submodule (caesar_sparse_monomial_module) caesar_sparse_monomial_submodule
  use caesar_anharmonic_common_module
contains

module procedure new_SparseMonomial
  this%modes = modes
end procedure

module procedure read_SparseMonomial
  type(ComplexUnivariate), allocatable :: modes(:)
  
  select type(this); type is(SparseMonomial)
    modes = ComplexUnivariate(tokens(input))
    this = SparseMonomial(modes)
  class default
    call err()
  end select
end procedure

module procedure write_SparseMonomial
  select type(this); type is(SparseMonomial)
    output = join(str(this%modes))
  class default
    call err()
  end select
end procedure

module procedure new_SparseMonomial_String
  call this%read(input)
end procedure
end submodule
