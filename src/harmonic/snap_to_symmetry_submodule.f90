submodule (caesar_snap_to_symmetry_module) caesar_snap_to_symmetry_submodule
  use caesar_harmonic_module
contains

module procedure snap_to_symmetry_mode
  output%mode_name = 'snap_to_symmetry'
  output%description = 'Uses spglib to snap a structure to a symmetry. The &
     &symmetry to which the structure is snapped is found by calling spglib &
     &with the given symmetry precision. The routine attempts to find the &
     &snapped structure which is as close to the input structure as possible.'
  output%keywords = [                                                         &
     & KeywordData( 'file_type',                                              &
     &              'file_type is the file type which will be used for &
     &single-point energy calculations. Settings are: "castep", &
     &"quantum_espresso", "caesar" and "xyz".',                               &
     &              default_value='castep'),                                  &
     & KeywordData( 'seedname',                                               &
     &              'seedname is the seedname from which file names are &
     &constructed.'),                                                         &
     & KeywordData( 'symmetry_precision',                                     &
     &              'symmetry_precision is the tolerance up to which a&
     &symmetry is accepted. In order for a symmetry to be accepted, it must &
     &transform the position of every atom to within symmetry_precision of an &
     &atom of the same element. There must be no other atom within &
     &symmetry_precision of this point. symmetry_precision should be given in &
     &Bohr. Ideally, symmetry_precision should be much smaller than the &
     &minimum inter-atomic distance, and much larger than the geometry &
     &optimisation tolerance.',                                               &
     &              default_value='0.1')                                      ]
  output%main_subroutine => snap_to_symmetry_subroutine
end procedure

module procedure snap_to_symmetry_subroutine
  ! User input variables.
  type(String) :: file_type
  type(String) :: seedname
  real(dp)     :: symmetry_precision
  
  ! The structure.
  type(StructureData) :: structure
  
  ! Electronic structure calculation writer.
  type(CalculationWriter) :: calculation_writer
  
  ! Variables for checking that the input is a primitive cell.
  type(IntMatrix) :: identity_matrix
  
  ! Files.
  type(String) :: input_filename
  
  ! Get settings from user, and check them.
  file_type = arguments%value('file_type')
  seedname = arguments%value('seedname')
  symmetry_precision = dble(arguments%value('symmetry_precision'))
  
  ! Initialise calculation writer.
  calculation_writer = CalculationWriter( file_type = file_type, &
                                        & seedname  = seedname   )
  
  ! Read in structure.
  input_filename = make_input_filename(file_type, seedname)
  structure = input_file_to_StructureData(file_type, input_filename)
  
  ! Snap structure to symmetry.
  structure = structure%snap_to_symmetry(symmetry_precision)
  
  ! Generate symmetries of structure.
  call structure%calculate_symmetry(symmetry_precision)
  
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
    call print_line('If the structure is known to be the primitive cell of &
       &the system, try lowering symmetry_precision.')
    call print_line('N.B. running phonon calculations on a system with a unit &
       &cell which is x*y*z primitive cells is equivalent to running the &
       &calculation on that primitive cell but with an x*y*z q-point grid. &
       &Please increase the q-point grid rather than running calculations on &
       &a supercell.')
    call err()
  endif
  
  ! Write out structure, q-point and supercell data.
  call calculation_writer%write_calculation( structure,              &
                                           & str('symmetry_snapped') )
end procedure
end submodule
