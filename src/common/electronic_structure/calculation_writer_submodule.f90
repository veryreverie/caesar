submodule (caesar_calculation_writer_module) caesar_calculation_writer_submodule
  use caesar_electronic_structure_module
contains

module procedure new_CalculationWriter
  integer :: ialloc
  
  this%file_type_      = file_type
  this%seedname_       = seedname
  this%input_filename_ = make_input_filename(file_type,seedname)
  allocate(this%directories_(0), stat=ialloc); call err(ialloc)
  
  ! Check that the input file exists.
  if (.not. file_exists(this%input_filename_)) then
    call print_line(ERROR//': The input file '//this%input_filename_// &
                   &' does not exist in the working directory.'        )
    call quit()
  endif
end procedure

module procedure directories_written
  output = this%directories_
end procedure

module procedure write_calculation
  type(String), allocatable :: temp(:)
  
  type(OFile) :: structure_file
  
  ! Check that the directory has not already been written to by this class.
  if (any(this%directories_==directory)) then
    call print_line(CODE_ERROR//': Trying to write an electronic structure &
       &calculation to a directory in which an electronic structure &
       &calculation has already been written by this execution of Caesar.')
    call print_line('Directory: '//directory)
    call err()
  endif
  
  ! Make the directory, and add a structure.dat file and
  !    an electronic structure input file.
  call mkdir(directory)
  structure_file = OFile(directory//'/structure.dat')
  call structure_file%print_lines(structure)
  call StructureData_to_input_file(              &
     & this%file_type_,                          &
     & structure,                                &
     & this%input_filename_,                     &
     & directory//'/'//this%input_filename_      )
  
  ! Record the directory.
  ! WORKAROUND: this is done via a temp array to avoid a compiler bug with
  !    ifort 19.0.4.
  temp = [this%directories_, directory]
  this%directories_ = temp
end procedure
end submodule
