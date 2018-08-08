! ======================================================================
! An array of type BasisFunction.
! ======================================================================
! Exists to allow heterogeneous storage of basis functions, as well as to
!    provide easy file I/O.
module basis_functions_module
  use common_module
  
  use anharmonic_common_module
  
  use basis_function_module
  implicit none
  
  private
  
  public :: BasisFunctions
  public :: generate_basis_functions
  
  type, extends(Stringsable) :: BasisFunctions
    type(BasisFunction), allocatable :: functions(:)
    type(RealMonomial),  allocatable :: unique_terms(:)
  contains
    procedure, public :: read  => read_BasisFunctions
    procedure, public :: write => write_BasisFunctions
  end type
  
  interface BasisFunctions
    module procedure new_BasisFunctions
    module procedure new_BasisFunctions_Strings
    module procedure new_BasisFunctions_StringArray
  end interface
  
  interface size
    module procedure size_BasisFunctions
  end interface
  
  interface generate_basis_functions
    module procedure generate_basis_functions_SubspaceMonomials
  end interface
contains

function new_BasisFunctions(functions,unique_terms) result(this)
  implicit none
  
  type(BasisFunction), intent(in) :: functions(:)
  type(RealMonomial),  intent(in) :: unique_terms(:)
  type(BasisFunctions)            :: this
  
  this%functions = functions
  this%unique_terms = unique_terms
end function

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
  type(DegenerateSubspace), intent(in)    :: subspaces(:)
  type(DegenerateSymmetry), intent(in)    :: degenerate_symmetries(:)
  logical,                  intent(in)    :: vscf_basis_functions_only
  type(OFile),              intent(inout) :: logfile
  type(BasisFunctions)                    :: output
  
  type(BasisFunction), allocatable :: basis_functions(:)
  
  integer                         :: unique_term_id
  type(RealMonomial)              :: unique_term
  type(RealMonomial)              :: matching_term
  
  integer :: i,j,k
  
  output%functions = [BasisFunction::]
  output%unique_terms = [RealMonomial::]
  do i=1,size(couplings)
    ! Generate all basis functions for the subspace monomial.
    basis_functions = generate_basis_functions( couplings(i),              &
                                              & structure,                 &
                                              & complex_modes,             &
                                              & real_modes,                &
                                              & qpoints,                   &
                                              & subspaces,                 &
                                              & degenerate_symmetries,     &
                                              & vscf_basis_functions_only, &
                                              & logfile                    )

    ! Take linear combinations of basis functions such that each basis function
    !    contains at least term which is in no other basis function.
    do j=1,size(basis_functions)
      ! Identify the largest term in basis function i.
      unique_term_id = maxloc(                                            &
         & abs(basis_functions(j)%real_representation%terms%coefficient), &
         & 1                                                              )
      unique_term = &
         & basis_functions(j)%real_representation%terms(unique_term_id)
      
      ! Subtract a multiple of basis function i from all other basis functions,
      !    such that the coefficient of unique_term_id(i) in all other basis
      !    functions is zero.
      do k=1,size(basis_functions)
        if (k/=j) then
          ! TODO
          matching_term = basis_functions(k)%real_representation%terms(first( &
             & basis_functions(k)%real_representation%terms,                  &
             & compare_real_monomials,                                        &
             & unique_term                                                   ))
          basis_functions(k) = basis_functions(k)        &
                           & - basis_functions(j)        &
                           & * matching_term%coefficient &
                           & / unique_term%coefficient
        endif
      enddo
      
      ! Append the unique term to the unique term array.
      output%unique_terms = [output%unique_terms, unique_term]
    enddo

    ! Append the basis functions to the basis function array.
    output%functions = [output%functions, basis_functions]
  enddo
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_BasisFunctions(this,input)
  implicit none
  
  class(BasisFunctions), intent(out) :: this
  type(String),          intent(in)  :: input(:)
  
  type(StringArray), allocatable :: sections(:)
  
  type(BasisFunction), allocatable :: functions(:)
  type(RealMonomial),  allocatable :: unique_terms(:)
  
  select type(this); type is(BasisFunctions)
    sections = split_into_sections(input)
    
    functions = BasisFunction(sections(:size(sections)-1))
    
    unique_terms = RealMonomial(str(sections(size(sections))))
    
    this = BasisFunctions(functions=functions, unique_terms=unique_terms)
  class default
    call err()
  end select
end subroutine

function write_BasisFunctions(this) result(output)
  implicit none
  
  class(BasisFunctions), intent(in) :: this
  type(String), allocatable         :: output(:)
  
  select type(this); type is(BasisFunctions)
    output = [ str(this%functions, separating_line=''), &
             & str(''),                                 &
             & str(this%unique_terms)                   ]
  class default
    call err()
  end select
end function

function new_BasisFunctions_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(BasisFunctions)     :: this
  
  call this%read(input)
end function

impure elemental function new_BasisFunctions_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(BasisFunctions)          :: this
  
  this = BasisFunctions(str(input))
end function
end module
