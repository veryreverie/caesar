submodule (caesar_calculation_reader_module) caesar_calculation_reader_submodule
  use caesar_electronic_structure_module
contains

module procedure new_CalculationReader
  integer :: ialloc
  
  if (present(loto_direction)) then
    this%loto_direction_ = loto_direction
  endif
  allocate(this%directories_(0), stat=ialloc); call err(ialloc)
end procedure

module procedure directories_read
  output = this%directories_
end procedure

module procedure read_calculation
  type(IFile)  :: electronic_structure_file
  
  type(IFile)          :: structure_file
  type(StructureData)  :: structure
  type(LotoCorrection) :: loto_correction
  
  electronic_structure_file = IFile(directory//'/electronic_structure.dat')
  output = ElectronicStructure(electronic_structure_file%lines())
  
  if (allocated(this%loto_direction_)) then
    if (.not. output%has_linear_response()) then
      call print_line(ERROR//': LO/TO splitting requested, but linear &
         &response data is not present in electronic structure file in &
         &directory '//directory)
      call err()
    elseif (.not. present(displacement)) then
      call print_line(CODE_ERROR//': LO/TO splitting requested, but &
         &displacement has not been passed to calculation reader.')
      call err()
    endif
    
    structure_file = IFile(directory//'/structure.dat')
    structure = StructureData(structure_file%lines())
    loto_correction = LotoCorrection( output%linear_response(), &
                                    & this%loto_direction_,     &
                                    & displacement,             &
                                    & structure                 )
    output = calculate_loto_correction(output, loto_correction)
  endif
  
  ! Record the directory.
  this%directories_ = [this%directories_, directory]
end procedure

module procedure read_calculations
  type(ElectronicStructure), allocatable  :: output(:)
  
  integer :: i,ialloc
  
  ! Read the calculations.
  allocate(output(size(directories)), stat=ialloc); call err(ialloc)
  do i=1,size(directories)
    output(i) = this%read_calculation(directories(i))
  enddo
  
  ! Record the directories.
  this%directories_ = [this%directories_, directories]
end procedure
end submodule
