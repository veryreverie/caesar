! ======================================================================
! A class which reads the results of electronic structure calculations.
! ======================================================================
module calculation_reader_module
  use utils_module
  
  use structure_module
  use normal_mode_module
  
  use electronic_structure_data_module
  use electronic_structure_common_module
  
  use electronic_structure_file_module
  implicit none
  
  private
  
  public :: CalculationReader
  
  type, extends(NoDefaultConstructor) :: CalculationReader
    type(FractionVector), private, allocatable :: loto_direction_
    type(String),         private, allocatable :: directories_(:)
  contains
    procedure, public :: directories_read
    procedure, public :: read_calculation
    procedure, public :: read_calculations
  end type
  
  interface CalculationReader
    module procedure new_CalculationReader
  end interface
contains

! Constructor.
function new_CalculationReader(loto_direction) result(this)
  implicit none
  
  type(FractionVector), intent(in), optional :: loto_direction
  type(CalculationReader)                    :: this
  
  if (present(loto_direction)) then
    this%loto_direction_ = loto_direction
  endif
  this%directories_ = [String::]
end function

! Return a list of the directories from which calculations have been read.
function directories_read(this) result(output)
  implicit none
  
  class(CalculationReader), intent(in) :: this
  type(String), allocatable            :: output(:)
  
  output = this%directories_
end function

! Read the results of an electronic structure calculation from the given
!    directory, and record the directory.
function read_calculation(this,directory,displacement) result(output)
  implicit none
  
  class(CalculationReader),    intent(inout)        :: this
  type(String),                intent(in)           :: directory
  type(CartesianDisplacement), intent(in), optional :: displacement
  type(ElectronicStructure)                         :: output
  
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
end function

! Read the results of a set of electronic structure calculations from the given
!    directories, and record the directories.
subroutine read_calculations(this,directories)
  implicit none
  
  class(CalculationReader), intent(inout) :: this
  type(String),             intent(in)    :: directories(:)
  type(ElectronicStructure), allocatable  :: output(:)
  
  integer :: i,ialloc
  
  ! Read the calculations.
  allocate(output(size(directories)), stat=ialloc); call err(ialloc)
  do i=1,size(directories)
    output(i) = this%read_calculation(directories(i))
  enddo
  
  ! Record the directories.
  this%directories_ = [this%directories_, directories]
end subroutine
end module
