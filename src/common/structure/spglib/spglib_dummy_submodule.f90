submodule (caesar_spglib_dummy_module) caesar_spglib_dummy_submodule
  use caesar_spglib_module
contains

module procedure new_SpglibSymmetries_calculated
  call print_line(ERROR//': Cannot calculate symmetries because Caesar has &
     &not been linked against spglib. Please use the CMake flag &
     &-DLINK_TO_SPGLIB:LOGICAL=true to link against spglib.')
  call quit()
end procedure

module procedure snap_to_symmetry
  call print_line(ERROR//': Cannot snap to symmetry because Caesar has &
     &not been linked against spglib. Please use the CMake flag &
     &-DLINK_TO_SPGLIB:LOGICAL=true to link against spglib.')
  call quit()
end procedure
end submodule
