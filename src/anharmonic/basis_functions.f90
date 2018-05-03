! ======================================================================
! An array of type BasisFunction.
! ======================================================================
! Exists to allow heterogeneous storage of basis functions, as well as to
!    provide easy file I/O.
module basis_functions_module
  use common_module
  
  use degeneracy_module
  use subspace_monomial_module
  use coupled_modes_module
  use degenerate_symmetry_module
  use basis_function_module
  implicit none
  
  private
  
  public :: BasisFunctions
  public :: generate_basis_functions
  
  type :: BasisFunctions
    type(BasisFunction), allocatable :: functions(:)
  contains
    procedure, public :: write_file => write_file_BasisFunctions
    procedure, public :: read_file  => read_file_BasisFunctions
  end type
  
  interface size
    module procedure size_BasisFunctions
  end interface
  
  interface generate_basis_functions
    module procedure generate_basis_functions_SubspaceMonomials
  end interface
contains

function size_BasisFunctions(this) result(output)
  implicit none
  
  type(BasisFunctions), intent(in) :: this
  integer                          :: output
  
  output = size(this%functions)
end function

function generate_basis_functions_SubspaceMonomials(couplings,structure, &
   & complex_modes,real_modes,qpoints,subspaces,degenerate_symmetries,   &
   & vscf_basis_functions_only,logfile) result(output)
  implicit none
  
  type(SubspaceMonomial),   intent(in)    :: couplings(:)
  type(StructureData),      intent(in)    :: structure
  type(ComplexMode),        intent(in)    :: complex_modes(:)
  type(RealMode),           intent(in)    :: real_modes(:)
  type(QpointData),         intent(in)    :: qpoints(:)
  type(DegenerateModes),    intent(in)    :: subspaces(:)
  type(DegenerateSymmetry), intent(in)    :: degenerate_symmetries(:)
  logical,                  intent(in)    :: vscf_basis_functions_only
  type(OFile),              intent(inout) :: logfile
  type(BasisFunctions)                    :: output
  
  integer :: i
  
  output%functions = [BasisFunction::]
  do i=1,size(couplings)
    output%functions = [ output%functions, &
                     &   generate_basis_functions( couplings(i),              &
                     &                             structure,                 &
                     &                             complex_modes,             &
                     &                             real_modes,                &
                     &                             qpoints,                   &
                     &                             subspaces,                 &
                     &                             degenerate_symmetries,     &
                     &                             vscf_basis_functions_only, &
                     &                             logfile)                   &
                     & ]
  enddo
end function

! ----------------------------------------------------------------------
! File I/O.
! ----------------------------------------------------------------------
subroutine write_file_BasisFunctions(this,filename)
  implicit none
  
  class(BasisFunctions), intent(in) :: this
  type(String),          intent(in) :: filename
  
  type(OFile) :: file
  
  integer :: i,j
  
  file = OFile(filename)
  do i=1,size(this%functions)
    call file%print_line('Basis function '//i)
    call file%print_lines(this%functions(i))
    call file%print_line('')
  enddo
end subroutine

subroutine read_file_BasisFunctions(this,filename)
  implicit none
  
  class(BasisFunctions), intent(out) :: this
  type(String),          intent(in)  :: filename
  
  type(IFile)                    :: file
  type(StringArray), allocatable :: sections(:)
  
  integer :: i,ialloc
  
  file = IFile(filename)
  sections = split(file%lines())
  allocate(this%functions(size(sections)), stat=ialloc); call err(ialloc)
  do i=1,size(sections)
    this%functions(i) = BasisFunction(sections(i)%strings(2:))
  enddo
end subroutine
end module
