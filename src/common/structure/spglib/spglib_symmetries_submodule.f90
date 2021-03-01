submodule (caesar_spglib_symmetries_module) caesar_spglib_symmetries_submodule
  use caesar_spglib_module
contains

module procedure new_SpglibSymmetries
  if (size(tensors)/=n_operations) then
    call print_line(ERROR//': n_operations does not match the number of &
       &tensors.')
    call err()
  elseif (size(translations)/=n_operations) then
    call print_line(ERROR//': n_operations does not match the number of &
       &translations.')
    call err()
  endif
  
  this%spacegroup_number = spacegroup_number
  this%international_symbol = international_symbol
  this%transformation = transformation
  this%origin_shift = origin_shift
  this%n_operations = n_operations
  this%tensors = tensors
  this%translations = translations
  this%n_atoms = n_atoms
  this%pointgroup_symbol = pointgroup_symbol
end procedure

module procedure size_SpglibSymmetries
  output = size(this%tensors)
end procedure

module procedure write_SpglibSymmetries
  integer :: i
  
  select type(this); type is(SpglibSymmetries)
    output = [ 'Spacegroup Number    : '//this%spacegroup_number,    &
             & 'International Symbol : '//this%international_symbol, &
             & 'Pointgroup Symbol    : '//this%pointgroup_symbol,    &
             & str('Transformation       : '),                       &
             & str(this%transformation),                             &
             & 'Origin Shift         : '//this%origin_shift,         &
             & 'No. Atoms            : '//this%n_atoms,              &
             & 'No. Operations       : '//this%n_operations,         &
             & str('Operations           : ')                        ]
    do i=1,this%n_operations
      output = [ output,                   &
               & str(this%tensors(i)),     &
               & str(this%translations(i)) ]
    enddo
  class default
    call err()
  end select
end procedure
end submodule
