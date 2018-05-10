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
  
  type, extends(Stringsable) :: BasisFunctions
    type(BasisFunction), allocatable :: functions(:)
  contains
    procedure, public :: read  => read_BasisFunctions
    procedure, public :: write => write_BasisFunctions
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
! I/O.
! ----------------------------------------------------------------------
subroutine read_BasisFunctions(this,input)
  implicit none
  
  class(BasisFunctions), intent(out) :: this
  type(String),          intent(in)  :: input(:)
  
  type(StringArray), allocatable :: functions(:)
  
  integer :: i,ialloc
  
  select type(this); type is(BasisFunctions)
    functions = split_into_sections(input)
    allocate(this%functions(size(functions)), stat=ialloc); call err(ialloc)
    do i=1,size(functions)
      this%functions(i) = functions(i)
    enddo
  end select
end subroutine

function write_BasisFunctions(this) result(output)
  implicit none
  
  class(BasisFunctions), intent(in) :: this
  type(String), allocatable         :: output(:)
  
  select type(this); type is(BasisFunctions)
    output = str(this%functions, separating_line='')
  end select
end function
end module