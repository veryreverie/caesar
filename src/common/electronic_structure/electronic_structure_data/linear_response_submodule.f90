submodule (caesar_linear_response_module) caesar_linear_response_submodule
  use caesar_electronic_structure_data_module
contains

module procedure new_LinearResponse
  this%permittivity = permittivity
  this%born_charges = born_charges
end procedure

module procedure read_LinearResponse
  type(RealMatrix)              :: permittivity
  type(RealMatrix), allocatable :: born_charges(:)
  
  select type(this); type is(LinearResponse)
    permittivity = RealMatrix(input(2:4))
    born_charges = RealMatrix(split_into_sections( &
                              & input(6:),         &
                              & separating_line='' ))
    
    this = LinearResponse(permittivity,born_charges)
  class default
    call err()
  end select
end procedure

module procedure write_LinearResponse
  select type(this); type is(LinearResponse)
    output = [ str('Permittivity'),                      &
             & str(this%permittivity),                   &
             & str('Born Effective Charges'),            &
             & str(this%born_charges,separating_line='') ]
  class default
    call err()
  end select
end procedure

module procedure new_LinearResponse_Strings
  call this%read(input)
end procedure

module procedure new_LinearResponse_StringArray
  this = LinearResponse(str(input))
end procedure
end submodule
