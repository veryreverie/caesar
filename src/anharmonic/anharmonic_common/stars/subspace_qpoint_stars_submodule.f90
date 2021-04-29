submodule (caesar_subspace_qpoint_stars_module) &
   & caesar_subspace_qpoint_stars_submodule
  use caesar_stars_module
contains

module procedure new_SubspaceQpointStars
  this%subspace_id = subspace_id
  this%powers = powers
end procedure

module procedure read_SubspaceQpointStars
  integer                        :: subspace_id
  type(QpointStars), allocatable :: powers(:)
  
  select type(this); type is(SubspaceQpointStars)
    subspace_id = int(token(input(1), 5))
    powers = QpointStars(split_into_sections( input(3:size(input)-1),        &
                                            & separating_line=repeat('-',50) ))
    this = SubspaceQpointStars(subspace_id, powers)
  class default
    call err()
  end select
end procedure

module procedure write_SubspaceQpointStars
  select type(this); type is(SubspaceQpointStars)
    output = [ 'q-point stars in subspace '//this%subspace_id//' :', &
             & str(repeat('-',50)),                                  &
             & str(this%powers, separating_line=repeat('-',50)),     &
             & str(repeat('-',50))                                   ]
  class default
    call err()
  end select
end procedure

module procedure new_SubspaceQpointStars_Strings
  call this%read(input)
end procedure

module procedure new_SubspaceQpointStars_StringArray
  call this%read(str(input))
end procedure
end submodule
