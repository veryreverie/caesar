! ======================================================================
! An example module, showing how to use PotentialData and PotentialPointer.
! ======================================================================
module potential_example_module
  use common_module
  
  use anharmonic_common_module
  use states_module
  
  use anharmonic_data_module
  implicit none
  
  private
  
  public :: startup_potential_example
  public :: potential_example_subroutine
  
  type, extends(PotentialData) :: PotentialDataExample
    type(String) :: example_contents
  contains
    procedure, public, nopass :: representation => &
                               & representation_PotentialDataExample
    
    procedure, public :: generate_sampling_points => &
       & generate_sampling_points_PotentialDataExample
    procedure, public :: generate_potential => &
       & generate_potential_PotentialDataExample
    procedure, public :: generate_stress => &
       & generate_stress_PotentialDataExample
    
    procedure, public :: zero_energy => zero_energy_PotentialDataExample
    
    procedure, public :: energy_RealModeDisplacement => &
                       & energy_RealModeDisplacement_PotentialDataExample
    procedure, public :: energy_ComplexModeDisplacement => &
                       & energy_ComplexModeDisplacement_PotentialDataExample
    procedure, public :: force_RealModeDisplacement => &
                       & force_RealModeDisplacement_PotentialDataExample
    procedure, public :: force_ComplexModeDisplacement => &
                       & force_ComplexModeDisplacement_PotentialDataExample
    
    procedure, public :: braket => braket_PotentialDataExample
    
    procedure, public :: harmonic_expectation => &
                       & harmonic_expectation_PotentialDataExample
    
    procedure, public :: iterate_damped => &
                       & iterate_damped_PotentialDataExample
    procedure, public :: iterate_pulay => &
                       & iterate_pulay_PotentialDataExample
    
    procedure, public :: read  => read_PotentialDataExample
    procedure, public :: write => write_PotentialDataExample
  end type
  
  interface PotentialDataExample
    module procedure new_PotentialDataExample
    module procedure new_PotentialDataExample_Strings
    module procedure new_PotentialDataExample_StringArray
  end interface
contains

! Startup procedure.
subroutine startup_potential_example()
  implicit none
  
  type(PotentialDataExample) :: potential
  
  call potential%startup()
end subroutine

! --------------------------------------------------
! Constructor for example class.
! --------------------------------------------------
! This is where any PotentialDataExample-specific data is input.
function new_PotentialDataExample(example_contents) result(this)
  implicit none
  
  type(String), intent(in)   :: example_contents
  type(PotentialDataExample) :: this
  
  this%example_contents = example_contents
end function

! Type representation.
impure elemental function representation_PotentialDataExample() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'example'
end function

! --------------------------------------------------
! Overloads of PotentialData's methods.
! --------------------------------------------------
subroutine generate_sampling_points_PotentialDataExample(this, &
   & anharmonic_data,sampling_points_dir,calculation_writer,logfile)
  implicit none
  
  class(PotentialDataExample), intent(inout) :: this
  type(AnharmonicData),        intent(in)    :: anharmonic_data
  type(String),                intent(in)    :: sampling_points_dir
  type(CalculationWriter),     intent(inout) :: calculation_writer
  type(OFile),                 intent(inout) :: logfile
  
  call print_line('PotentialDataExample: generating sampling points.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to generate sampling points goes here.
end subroutine

subroutine generate_potential_PotentialDataExample(this,anharmonic_data, &
   & weighted_energy_force_ratio,sampling_points_dir,calculation_reader, &
   & logfile)
  implicit none
  
  class(PotentialDataExample), intent(inout) :: this
  type(AnharmonicData),        intent(in)    :: anharmonic_data
  real(dp),                    intent(in)    :: weighted_energy_force_ratio
  type(String),                intent(in)    :: sampling_points_dir
  type(CalculationReader),     intent(inout) :: calculation_reader
  type(OFile),                 intent(inout) :: logfile
  
  call print_line('PotentialDataExample: generating potential.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to generate sampling points goes here.
end subroutine

function generate_stress_PotentialDataExample(this,anharmonic_data,       &
   & sampling_points_dir,stress_expansion_order,stress_subspace_coupling, &
   & vscf_basis_functions_only,calculation_reader,logfile) result(output)
  implicit none
  
  class(PotentialDataExample), intent(in)    :: this
  type(AnharmonicData),        intent(in)    :: anharmonic_data
  type(String),                intent(in)    :: sampling_points_dir
  integer,                     intent(in)    :: stress_expansion_order
  type(SubspaceCoupling),      intent(in)    :: stress_subspace_coupling(:)
  logical,                     intent(in)    :: vscf_basis_functions_only
  type(CalculationReader),     intent(inout) :: calculation_reader
  type(OFile),                 intent(inout) :: logfile
  type(StressPointer)                        :: output
  
  call print_line('PotentialDataExample: generating stress.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to generate sampling points goes here.
end function

impure elemental subroutine zero_energy_PotentialDataExample(this)
  implicit none
  
  class(PotentialDataExample), intent(inout) :: this
  
  call print_line('PotentialDataExample: zeroing energy.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to zero the energy (s.t. undisplaced_energy()=0) goes here.
end subroutine

impure elemental function energy_RealModeDisplacement_PotentialDataExample( &
   & this,displacement) result(output)
  implicit none
  
  class(potentialDataExample), intent(in) :: this
  type(RealModeDisplacement),  intent(in) :: displacement
  real(dp)                                :: output
  
  call print_line('PotentialDataExample: calculating energy.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to calculate energies at real displacements goes here.
end function

impure elemental function energy_ComplexModeDisplacement_PotentialDataExample(&
   & this,displacement) result(output)
  implicit none
  
  class(potentialDataExample),    intent(in) :: this
  type(ComplexModeDisplacement),  intent(in) :: displacement
  complex(dp)                                :: output
  
  call print_line('PotentialDataExample: calculating energy.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to calculate energies at complex displacements goes here.
end function

impure elemental function force_RealModeDisplacement_PotentialDataExample( &
   & this,displacement) result(output)
  implicit none
  
  class(potentialDataExample), intent(in) :: this
  type(RealModeDisplacement),  intent(in) :: displacement
  type(RealModeForce)                     :: output
  
  call print_line('PotentialDataExample: calculating force.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to calculate forces at real displacements goes here.
end function

impure elemental function force_ComplexModeDisplacement_PotentialDataExample( &
   & this,displacement) result(output)
  implicit none
  
  class(potentialDataExample),    intent(in) :: this
  type(ComplexModeDisplacement),  intent(in) :: displacement
  type(ComplexModeForce)                     :: output
  
  call print_line('PotentialDataExample: calculating force.')
  call print_line('Example contents = '//this%example_contents)
  
  ! Code to calculate forces at complex displacements goes here.
end function

function braket_PotentialDataExample(this,bra,ket,subspace,subspace_basis, &
   & anharmonic_data) result(output)
  implicit none
  
  class(PotentialDataExample), intent(in)           :: this
  class(SubspaceState),        intent(in)           :: bra
  class(SubspaceState),        intent(in), optional :: ket
  type(DegenerateSubspace),    intent(in)           :: subspace
  class(SubspaceBasis),        intent(in)           :: subspace_basis
  type(AnharmonicData),        intent(in)           :: anharmonic_data
  type(PotentialPointer)                            :: output
  
  call print_line('PotentialDataExample: evaluating <bra|potential|ket>.')
  
  ! Code to integrate this potential between <bra| and |ket> goes here.
end function

function harmonic_expectation_PotentialDataExample(this,frequency, &
   & thermal_energy,no_states,subspace,anharmonic_data) result(output)
  implicit none
  
  class(PotentialDataExample), intent(in) :: this
  real(dp),                    intent(in) :: frequency
  real(dp),                    intent(in) :: thermal_energy
  integer,                     intent(in) :: no_states
  type(DegenerateSubspace),    intent(in) :: subspace
  type(AnharmonicData),        intent(in) :: anharmonic_data
  real(dp)                                :: output
  
  call print_line('PotentialDataExample: evaluating harmonic expectation of &
     &<V>.')
  
  ! Code to calculate thermal harmonic expectation goes here.
end function

impure elemental function iterate_damped_PotentialDataExample(this, &
   & new_potential,damping,anharmonic_data) result(output)
  implicit none
  
  class(PotentialDataExample), intent(in) :: this
  class(PotentialData),        intent(in) :: new_potential
  real(dp),                    intent(in) :: damping
  type(AnharmonicData),        intent(in) :: anharmonic_data
  type(PotentialPointer)                  :: output
  
  ! Code to generate the potential which is equal to:
  !    (1-damping)*this + damping*new_potential
  ! goes here.
end function

function iterate_pulay_PotentialDataExample(this,input_potentials, &
   & output_potentials,anharmonic_data) result(output)
  implicit none
  
  class(PotentialDataExample), intent(in) :: this
  type(PotentialPointer),      intent(in) :: input_potentials(:)
  type(PotentialPointer),      intent(in) :: output_potentials(:)
  type(AnharmonicData),        intent(in) :: anharmonic_data
  type(PotentialPointer)                  :: output
  
  ! Code to generate a new potential using a Pulay scheme goes here.
  ! See the PolynomialPotential Pulay scheme for an example.
end function

! --------------------------------------------------
! I/O.
! --------------------------------------------------
subroutine read_PotentialDataExample(this,input)
  implicit none
  
  class(PotentialDataExample), intent(out) :: this
  type(String),                intent(in)  :: input(:)
  
  ! Code to read potential from strings goes here.
end subroutine

function write_PotentialDataExample(this) result(output)
  implicit none
  
  class(PotentialDataExample), intent(in) :: this
  type(String), allocatable               :: output(:)
  
  ! Code to write potential to strings goes here.
end function

function new_PotentialDataExample_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)   :: input(:)
  type(PotentialDataExample) :: this
  
  call this%read(input)
end function

impure elemental function new_PotentialDataExample_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(PotentialDataExample)    :: this
  
  this = PotentialDataExample(str(input))
end function

! --------------------------------------------------
! The class in use.
! --------------------------------------------------
subroutine potential_example_subroutine()
  implicit none
  
  ! A polymorphic pointer, which can store an object of any type which
  !    extends PotentialData.
  type(PotentialPointer) :: potential
  
  ! An example variable with which to initialise a PotentialDataExample.
  type(String) :: example_contents
  
  ! Variables for generate_sampling points.
  type(AnharmonicData)    :: anharmonic_data
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
  type(DegenerateSubspace)   :: subspace
  type(SubspaceBasisPointer) :: subspace_basis
  type(MonomialState)        :: state_1
  type(MonomialState)        :: state_2
  
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
  call potential%generate_sampling_points( anharmonic_data,     &
                                         & sampling_points_dir, &
                                         & calculation_writer,  &
                                         & logfile              )
  
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
  potential = braket( state_1,        &
                    & potential,      &
                    & state_2,        &
                    & subspace,       &
                    & subspace_basis, &
                    & anharmonic_data )
  
  ! The potential can be written to and read from file using the potential
  !    pointer's methods.
  output_file = OFile('example_potential.file')
  call output_file%print_lines(potential)
  
  input_file = IFile('example_potential.file')
  potential = PotentialPointer(input_file%lines())
end subroutine
end module