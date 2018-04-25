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
  
  file = filename
  
  do i=1,size(this%functions)
    call file%print_line('Basis function '//i)
    call file%print_line(this%functions(i))
    call file%print_line('')
  enddo
end subroutine

subroutine read_file_BasisFunctions(this,filename)
  implicit none
  
  class(BasisFunctions), intent(out) :: this
  type(String),          intent(in)  :: filename
  
  type(IFile) :: file
  
  integer, allocatable :: first_lines(:)
  integer, allocatable :: last_lines(:)
  
  integer :: no_basis_functions
  
  integer :: i,ialloc
  
  file = filename
  
  allocate( first_lines(size(file)), &
          & last_lines(size(file)),  &
          & stat=ialloc); call err(ialloc)
  
  ! Identify separate basis functions by the blank lines between them.
  no_basis_functions = 1
  first_lines(1) = 2
  do i=2,size(file)
    if (len(file%line(i))==0) then
      last_lines(no_basis_functions) = i-1
      no_basis_functions = no_basis_functions+1
      first_lines(no_basis_functions) = i+2
    endif
  enddo
  last_lines(no_basis_functions) = size(file)
  
  first_lines = first_lines(:no_basis_functions)
  last_lines = last_lines(:no_basis_functions)
  
  allocate(this%functions(no_basis_functions), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    this%functions(i) = BasisFunction(file%lines(first_lines(i),last_lines(i)))
  enddo
end subroutine

end module
