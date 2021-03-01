submodule (caesar_calculation_runner_module) caesar_calculation_runner_submodule
  use caesar_electronic_structure_module
contains

module procedure new_CalculationRunner
  integer :: ialloc
    
  this%file_type_           = file_type
  this%seedname_            = seedname
  this%run_script_          = run_script
  this%no_cores_            = no_cores
  this%no_nodes_            = no_nodes
  this%run_script_data_     = run_script_data
  this%calculation_type_    = calculation_type
  this%use_forces_          = use_forces
  this%use_hessians_        = use_hessians
  this%calculate_stress_    = calculate_stress
  this%exit_on_error_       = exit_on_error
  this%repeat_calculations_ = repeat_calculations
  allocate(this%directories_(0), stat=ialloc); call err(ialloc)
end procedure

module procedure directories_run
  output = this%directories_
end procedure

module procedure run_calculation
  type(String) :: working_directory
  
  integer                   :: result_code
  type(IFile)               :: structure_file
  type(StructureData)       :: structure
  type(ElectronicStructure) :: electronic_structure
  type(OFile)               :: electronic_structure_file
  
  working_directory = format_path('.')
  
  ! Check that this directory has not already been run.
  if (any(this%directories_==directory)) then
    call print_line(CODE_ERROR//': Trying to run an electronic structure &
       &calculation in a directory in which an electronic structure &
       &calculation has already been run by this execution of Caesar.')
    call print_line('Directory: '//directory)
    call err()
  endif
  
  ! Record the directory.
  this%directories_ = [this%directories_, directory]
  
  ! Check if the calculation has been run succesfully before.
  if (.not. this%repeat_calculations_) then
    if (file_exists(directory//'/electronic_structure.dat')) then
      call print_line('Skipping successful calculation in directory '// &
         & directory)
      return
    endif
  endif
  
  ! Run the calculation.
  call print_line('Running calculation in directory '//directory)
  result_code = system_call( 'cd '//working_directory//';' //' '// &
                           & this%run_script_              //' '// &
                           & this%file_type_               //' '// &
                           & directory                     //' '// &
                           & this%no_cores_                //' '// &
                           & this%no_nodes_                //' '// &
                           & this%seedname_                //' '// &
                           & this%run_script_data_                 )
  call print_line('Result code: '//result_code)
  
  ! Check the result code.
  if (result_code/=0) then
    if (this%exit_on_error_) then
      call quit()
    else
      return
    endif
  endif
  
  ! Convert the electronic structure result into an ElectronicStructure.
  structure_file = IFile(directory//'/structure.dat')
  structure = StructureData(structure_file%lines())
  electronic_structure = read_output_file( this%file_type_,        &
                                         & structure,              &
                                         & directory,              &
                                         & this%seedname_,         &
                                         & this%calculation_type_, &
                                         & this%use_forces_,       &
                                         & this%use_hessians_,     &
                                         & this%calculate_stress_  )
  
  ! Check that the objects calculated correspond to those requested.
  if (this%use_forces_) then
    if (.not. electronic_structure%has_forces()) then
      call print_line(ERROR//': Use of forces has been requested, but no &
         &forces are present in electronic structure.')
      call err()
    endif
  endif
  
  if (this%use_hessians_) then
    if (.not. electronic_structure%has_hessian()) then
      call print_line(ERROR//': Use of forces has been requested, but no &
         &forces are present in electronic structure.')
      call err()
    endif
  endif
  
  if (this%calculate_stress_) then
    if (.not. electronic_structure%has_stress()) then
      call print_line(ERROR//': Stress calculation has been requested, but no &
         &stress tensor or stress tensor is present in electronic structure.')
      call err()
    endif
  endif
  
  ! Print the electronic structure to file.
  electronic_structure_file = OFile(directory//'/electronic_structure.dat' )
  call electronic_structure_file%print_lines(electronic_structure)
end procedure

module procedure run_calculations
  integer :: i
  
  do i=1,size(directories)
    call print_line(i//' of '//size(directories)//':')
    call this%run_calculation(directories(i))
  enddo
end procedure
end submodule
