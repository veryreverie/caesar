submodule (caesar_harmonic_data_module) caesar_harmonic_data_submodule
  use caesar_harmonic_module
contains

module procedure new_HarmonicData
  this%structure = structure
  this%large_supercell = large_supercell
  this%qpoints = qpoints
end procedure

module procedure new_HarmonicData_arguments
  type(String)                  :: input_filename
  type(StructureData)           :: structure
  type(IntMatrix)               :: identity_matrix
  type(IntMatrix)               :: large_supercell_matrix
  type(StructureData)           :: large_supercell
  type(QpointData), allocatable :: qpoints(:)
  
  ! Read in structure
  input_filename = make_input_filename(file_type, seedname)
  structure = input_file_to_StructureData(file_type, input_filename)
  
  ! Snap structure to symmetry.
  if (snap_to_symmetry) then
    structure = structure%snap_to_symmetry(symmetry_precision)
  endif
  
  ! Generate symmetries of structure.
  call print_line('Calculating symmetry of input structure.')
  if (present(loto_direction)) then
    if (size(structure%symmetries)==0) then
      call structure%calculate_symmetry( symmetry_precision,             &
                                       & loto_direction = loto_direction )
    else
      if (any(loto_breaks_symmetry( structure%symmetries%tensor, &
                                  & loto_direction               ))) then
        call print_line(ERROR//': Symmetries have been specified which are &
           &broken by the LO/TO direction.')
        call quit()
      endif
    endif
  else
    if (size(structure%symmetries)==0) then
      call structure%calculate_symmetry(symmetry_precision)
    endif
  endif
  
  ! Check that the given structure is the primitive cell.
  identity_matrix = make_identity_matrix(3)
  if (count(structure%symmetries%tensor==identity_matrix)==0) then
    call print_line(ERROR//': The identity symmetry is not present in the &
       &symmetries of the input structure.')
    call err()
  elseif (count(structure%symmetries%tensor==identity_matrix)>1) then
    call print_line(ERROR//': The input structure has a purely translational &
       &symmetry. This is usually because it is not a primitive cell of the &
       &system.')
    call print_line('If the structure is known to be a primitive cell of &
       &the system, try lowering symmetry_precision.')
    call print_line('N.B. running phonon calculations on a system with a unit &
       &cell which is x*y*z primitive cells is equivalent to running the &
       &calculation on that primitive cell but with an x*y*z q-point grid. &
       &Please increase the q-point grid rather than running calculations on &
       &a supercell.')
    call err()
  endif
  
  ! Generate large supercell, for which all q-points are G-vectors.
  call print_line('Generating large supercell. &
     &(This may take a while if the q-point grid is large).')
  large_supercell_matrix = mat( [ grid(1), 0      , 0     ,    &
                              &   0      , grid(2), 0     ,    &
                              &   0      , 0      , grid(3) ], &
                              & 3,3)
  large_supercell = construct_supercell( structure,             &
                                       & large_supercell_matrix )
  
  ! Generate q-points in Monkhorst-Pack grid.
  qpoints = generate_qpoints(large_supercell)
  
  ! Construct output.
  this = HarmonicData(structure, large_supercell, qpoints)
end procedure
end submodule
