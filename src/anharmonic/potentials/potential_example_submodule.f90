submodule (caesar_potential_example_module) caesar_potential_example_submodule
  use caesar_potentials_module
contains

module procedure new_PotentialDataExample
  this%example_contents = example_contents
end procedure

module procedure representation_PotentialDataExample
  output = 'example'
end procedure

module procedure generate_sampling_points_PotentialDataExample
  call print_line('PotentialDataExample: generating sampling points.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to generate sampling points goes here.
end procedure

module procedure generate_potential_PotentialDataExample
  call print_line('PotentialDataExample: generating potential.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to generate sampling points goes here.
end procedure

module procedure generate_stress_PotentialDataExample
  call print_line('PotentialDataExample: generating stress.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to generate sampling points goes here.
end procedure

module procedure zero_energy_PotentialDataExample
  call print_line('PotentialDataExample: zeroing energy.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to zero the energy (s.t. undisplaced_energy()=0) goes here.
end procedure

module procedure add_constant_PotentialDataExample
  call print_line('PotentialDataExample: adding constant.')
  
  ! Code to add a constant to the potential goes here.
end procedure

module procedure energy_RealModeDisplacement_PotentialDataExample
  call print_line('PotentialDataExample: calculating energy.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to calculate energies at real displacements goes here.
end procedure

module procedure energy_ComplexModeDisplacement_PotentialDataExample
  call print_line('PotentialDataExample: calculating energy.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to calculate energies at complex displacements goes here.
end procedure

module procedure force_RealModeDisplacement_PotentialDataExample
  call print_line('PotentialDataExample: calculating force.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to calculate forces at real displacements goes here.
end procedure

module procedure force_ComplexModeDisplacement_PotentialDataExample
  call print_line('PotentialDataExample: calculating force.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to calculate forces at complex displacements goes here.
end procedure

module procedure braket_SubspaceBraKet_PotentialDataExample
  call print_line('PotentialDataExample: evaluating <bra|potential|ket>.')
  
  ! Code to integrate this potential between <bra| and |ket> goes here.
  ! This likely just involves calling braket on the constituent parts of this.
end procedure

module procedure braket_BasisState_PotentialDataExample
  call print_line('PotentialDataExample: evaluating <bra|potential|ket>.')
  
  ! Code to integrate this potential between <bra| and |ket> goes here.
  ! This likely just involves calling braket on the constituent parts of this.
end procedure

module procedure braket_BasisStates_PotentialDataExample
  call print_line('PotentialDataExample: evaluating <potential>.')
  
  ! Code to integrate this potential between the states goes here.
  ! This likely just involves calling braket on the constituent parts of this.
end procedure

module procedure harmonic_expectation_PotentialDataExample
  call print_line('PotentialDataExample: evaluating harmonic expectation of &
     &<V>.')
  
  ! Code to calculate thermal harmonic expectation goes here.
end procedure

module procedure coefficients_PotentialDataExample
  ! Code to convert the potential to an array of real coefficients goes here.
end procedure

module procedure set_coefficients_PotentialDataExample
  ! Code to convert the coefficients into the potential goes here.
end procedure

module procedure all_basis_functions_PotentialDataExample
  ! Code to return all basis functions goes here.
end procedure

module procedure variable_basis_functions_PotentialDataExample
  ! Code to return the basis functions corresponding to coefficients() goes
  !    here.
end procedure

module procedure read_PotentialDataExample
  ! Code to read potential from strings goes here.
end procedure

module procedure write_PotentialDataExample
  ! Code to write potential to strings goes here.
end procedure

module procedure new_PotentialDataExample_Strings
  call this%read(input)
end procedure

module procedure new_PotentialDataExample_StringArray
  this = PotentialDataExample(str(input))
end procedure

module procedure potential_example_subroutine
  ! A polymorphic pointer, which can store an object of any type which
  !    extends PotentialData.
  type(PotentialPointer) :: potential
  
  ! An example variable with which to initialise a PotentialDataExample.
  type(String) :: example_contents
  
  ! Variables for generate_sampling points.
  type(AnharmonicData)    :: anharmonic_data
  logical                 :: use_forces
  logical                 :: use_hessian
  logical                 :: calculate_stress
  type(String)            :: sampling_points_dir
  type(CalculationWriter) :: calculation_writer
  type(OFile)             :: logfile
  
  ! Variables for generate_potential.
  real(dp)                :: weighted_energy_force_ratio
  type(CalculationReader) :: calculation_reader
  
  ! Variables for energy and force.
  type(RealModeDisplacement)    :: real_displacement
  real(dp)                      :: real_energy
  type(RealModeForce)           :: real_force
  type(ComplexModeDisplacement) :: complex_displacement
  complex(dp)                   :: complex_energy
  type(ComplexModeForce)        :: complex_force
  
  ! Variables for integrating potential.
  type(HarmonicStateReal), target :: state_1
  type(HarmonicStateReal), target :: state_2
  type(HarmonicBraKetReal)        :: braket
  
  ! Files.
  type(OFile) :: output_file
  type(IFile) :: input_file
  
  ! Set the pointer to point to a PotentialDataExample type.
  ! This is where any PotentialDataExample-specific data is input,
  !    in this case the variable example_contents.
  example_contents = 'example'
  potential = PotentialPointer(PotentialDataExample(example_contents))
  
  ! Now PotentialData's methods can be called.
  ! They will all be forwareded to the PotentialDataExample instance.
  
  ! Generates sampling points, in a manner specific to the representation.
  call potential%generate_sampling_points( anharmonic_data,             &
                                         & use_forces,                  &
                                         & weighted_energy_force_ratio, &
                                         & use_hessian,                 &
                                         & calculate_stress,            &
                                         & sampling_points_dir,         &
                                         & calculation_writer,          &
                                         & logfile                      )
  
  ! Code to run electronic structure goes here.
  
  ! Generates the potential, in a manner specific to the representation.
  call potential%generate_potential( anharmonic_data,             &
                                   & weighted_energy_force_ratio, &
                                   & sampling_points_dir,         &
                                   & calculation_reader,          &
                                   & logfile                      )
  
  ! Once the potential has been generated, it can be used to calculate
  !    energies and forces.
  real_energy = potential%energy(real_displacement)
  real_force  = potential%force(real_displacement)
  complex_energy = potential%energy(complex_displacement)
  complex_force  = potential%force(complex_displacement)
  
  ! The potential can also be integrated between two states.
  ! If whole_subspace is .true., then the potential will be integrated across
  !    the whole subspace. whole_subspace default to true.
  braket%bra_ => state_1
  braket%ket_ => state_2
  call potential%braket( braket          = braket,         &
                       & whole_subspace  = .true.,         &
                       & anharmonic_data = anharmonic_data )
  
  ! The potential can be written to and read from file using the potential
  !    pointer's methods.
  output_file = OFile('example_potential.file')
  call output_file%print_lines(potential)
  
  input_file = IFile('example_potential.file')
  potential = PotentialPointer(input_file%lines())
end procedure
end submodule
