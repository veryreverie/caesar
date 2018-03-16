! ======================================================================
! Generates basis functions.
! ======================================================================
module basis_function_module
  use common_module
  
  use coupling_module
  use polynomial_module
  use degeneracy_module
  use degenerate_symmetry_module
  implicit none
  
  private
  
  public :: generate_basis_functions
contains

function generate_basis_functions(coupling,normal_modes,qpoints, &
   & subspaces,degenerate_symmetries,vscf_basis_functions_only) result(output)
  implicit none
  
  type(CoupledSubspaces),   intent(in) :: coupling
  type(ComplexMode),        intent(in) :: normal_modes(:)
  type(QpointData),         intent(in) :: qpoints(:)
  type(DegenerateModes),    intent(in) :: subspaces(:)
  type(DegenerateSymmetry), intent(in) :: degenerate_symmetries(:)
  logical,                  intent(in) :: vscf_basis_functions_only
  type(Polynomial), allocatable        :: output(:)
  
  type(DegenerateModes), allocatable :: coupled_subspaces(:)
  
  type(CoupledModes), allocatable :: coupled_modes_with_duplicates(:)
  type(CoupledModes), allocatable :: coupled_modes(:)
  integer,            allocatable :: no_duplicates(:)
  
  type(ComplexMatrix) :: symmetry
  
  type(DegenerateModes) :: subspace
  
  integer :: i,j,k,l,ialloc
  
  integer :: j2,k2
  
  integer :: m
  
  if (size(coupling)<2) then
    call print_line(CODE_ERROR//': Trying to generate basis functions with &
       &order less than 2.')
    call err()
  endif
  
  ! List the subspaces which are coupled together in this coupling.
  coupled_subspaces = coupling%coupled_subspaces(subspaces)
  
  ! Generate every allowed mode coupling within the coupled subspaces.
  coupled_modes_with_duplicates = generate_mode_coupling( &
                                     & coupled_subspaces, &
                                     & normal_modes,      &
                                     & qpoints,           &
                                     & vscf_basis_functions_only)
  
  ! De-duplicate the list of mode couplings, and record how many times each
  !    is duplicated.
  coupled_modes =                                                        &
     & coupled_modes_with_duplicates(set( coupled_modes_with_duplicates, &
     &                                    equality_coupled_modes))
  allocate( no_duplicates(size(coupled_modes)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(coupled_modes)
    no_duplicates(i) = count(coupled_modes_with_duplicates==coupled_modes(i))
  enddo
  
  if (sum(no_duplicates)/=size(coupled_modes_with_duplicates)) then
    call print_line(CODE_ERROR//': error de-duplicating mode couplings.')
    call err()
  endif
  
  do i=1,size(degenerate_symmetries)
    symmetry = degenerate_symmetries(i)%calculate_symmetry( coupled_modes, &
                                                          & no_duplicates)
    call print_line('')
    call print_line(sum_squares(real(symmetry*conjg(transpose(symmetry)))-make_identity_matrix(size(symmetry,1))))
  enddo
contains
  ! Lambda for comparing coupled_modes.
  function equality_coupled_modes(this,that) result(output)
    implicit none
    
    class(*), intent(in) :: this
    class(*), intent(in) :: that
    logical              :: output
    
    select type(this); type is(CoupledModes)
      select type(that); type is(CoupledModes)
        output = this==that
      end select
    end select
  end function
end function
end module
