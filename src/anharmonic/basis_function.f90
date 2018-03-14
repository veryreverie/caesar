! ======================================================================
! Generates basis functions.
! ======================================================================
module basis_function_module
  use common_module
  
  use coupling_module
  use polynomial_module
  use degeneracy_module
  implicit none
  
  private
  
  public :: generate_basis_functions
contains

function generate_basis_functions(coupling,normal_modes,qpoints, &
   & subspaces,symmetries) result(output)
  implicit none
  
  type(CoupledSubspaces), intent(in) :: coupling
  type(ComplexMode),      intent(in) :: normal_modes(:)
  type(QpointData),       intent(in) :: qpoints(:)
  type(DegenerateModes),  intent(in) :: subspaces(:)
  type(ComplexMatrix),    intent(in) :: symmetries(:,:)
  type(Polynomial), allocatable      :: output(:)
  
  type(DegenerateModes), allocatable :: coupled_subspaces(:)
  
  type(IntArray1D), allocatable :: mode_ids_with_duplicates(:)
  type(IntArray1D), allocatable :: mode_ids(:)
  integer,          allocatable :: no_duplicates(:)
  
  complex(dp), allocatable :: symm_mats(:,:,:)
  
  integer :: i,j,k,l,ialloc
  
  integer :: m
  
  if (size(coupling)<2) then
    call print_line(CODE_ERROR//': Trying to generate basis functions with &
       &order less than 2.')
    call err()
  endif
  
  ! List the subspaces which are coupled together in this coupling.
  coupled_subspaces = coupling%coupled_subspaces(subspaces)
  
  ! Generate every allowed mode coupling within the coupled subspaces.
  mode_ids_with_duplicates = generate_coupled_modes( coupled_subspaces, &
                                                   & normal_modes,      &
                                                   & qpoints)
  
  ! De-duplicate the list of mode couplings, and record how many times each
  !    is duplicated.
  mode_ids = mode_ids_with_duplicates(set( mode_ids_with_duplicates, &
                                         & equality_mode_ids))
  allocate( no_duplicates(size(mode_ids)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(mode_ids)
    no_duplicates(i) = 0
    do j=1,size(mode_ids_with_duplicates)
      if (all(mode_ids_with_duplicates(j)%i==mode_ids(i)%i)) then
        no_duplicates(i) = no_duplicates(i) + 1
      endif
    enddo
  enddo
  
  do i=1,size(mode_ids)
    call print_line(mode_ids(i))
  enddo
  
  if (sum(no_duplicates)/=size(mode_ids_with_duplicates)) then
    call print_line(CODE_ERROR//': error de-duplicating mode couplings.')
    call err()
  endif
  
  allocate( symm_mats(size(mode_ids),size(mode_ids),size(symmetries)), &
          & stat=ialloc); call err(ialloc)
  m = 0
  symm_mats = cmplx(1.0_dp,0.0_dp,dp)
  call print_line('No symms: '//size(symmetries,2))
  call print_line('Coupling size: '//size(coupling))
  do i=1,size(symmetries,2)
    do j=1,size(mode_ids)
      do k=1,size(mode_ids)
        do l=1,size(coupling)
          m = m+1
        enddo
      enddo
    enddo
  enddo
  
  call print_line(m)
contains
  ! Lambda for comparing mode_ids.
  function equality_mode_ids(this,that) result(output)
    implicit none
    
    class(*), intent(in) :: this
    class(*), intent(in) :: that
    logical              :: output
    
    select type(this); type is(IntArray1D)
      select type(that); type is(IntArray1D)
        output = all(this%i==that%i)
      end select
    end select
  end function
end function

! ----------------------------------------------------------------------
! Recursively generates all sets of coupled modes within a set of coupled
!    subspaces.
! Only returns couplings with sum(q)=0, modulo G-vectors.
! ----------------------------------------------------------------------
recursive function generate_coupled_modes(coupled_subspaces,normal_modes, &
   & qpoints,mode_ids_in,sum_q_in) result(output)
  implicit none
  
  type(DegenerateModes), intent(in)           :: coupled_subspaces(:)
  type(ComplexMode),     intent(in)           :: normal_modes(:)
  type(QpointData),      intent(in)           :: qpoints(:)
  type(IntArray1D),      intent(in), optional :: mode_ids_in
  type(FractionVector),  intent(in), optional :: sum_q_in
  type(IntArray1D), allocatable               :: output(:)
  
  type(IntArray1D)              :: mode_ids
  type(QpointData), allocatable :: subspace_qpoints(:)
  type(FractionVector)          :: sum_q
  
  integer :: i
  
  if (present(mode_ids_in) .neqv. present(sum_q_in)) then
    call print_line(CODE_ERROR//': generate_coupled_modes must be called with &
       &both mode_ids_in and sum_q_in or with neither.')
    call err()
  endif
  
  if (.not. present(mode_ids_in)) then
    subspace_qpoints = coupled_subspaces(1)%qpoints(qpoints)
    output = [IntArray1D::]
    do i=1,size(coupled_subspaces(1))
      mode_ids = [coupled_subspaces(1)%mode_ids(i)]
      sum_q = subspace_qpoints(i)%qpoint
      output = [ output,                                        &
             &   generate_coupled_modes( coupled_subspaces(2:), &
             &                           normal_modes,          &
             &                           qpoints,               &
             &                           mode_ids,              &
             &                           sum_q)                 &
             & ]
    enddo
  elseif (size(coupled_subspaces)==0) then
    if (is_int(sum_q_in)) then
      mode_ids = mode_ids_in%i(sort(mode_ids_in%i))
      output = [mode_ids]
    else
      output = [IntArray1D::]
    endif
  else
    subspace_qpoints = coupled_subspaces(1)%qpoints(qpoints)
    output = [IntArray1D::]
    do i=1,size(coupled_subspaces(1))
      mode_ids = [mode_ids_in%i, coupled_subspaces(1)%mode_ids(i)]
      sum_q = sum_q_in + subspace_qpoints(i)%qpoint
      output = [ output,                                        &
             &   generate_coupled_modes( coupled_subspaces(2:), &
             &                           normal_modes,          &
             &                           qpoints,               &
             &                           mode_ids,              &
             &                           sum_q)                 &
             & ]
    enddo
  endif
end function

end module
