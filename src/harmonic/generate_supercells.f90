! ======================================================================
! Generates the supercells needed to simulate harmonic phonons at all q-points.
! ======================================================================
module generate_supercells_module
  use common_module
  implicit none
  
  private
  
  public :: generate_supercells
contains

! ----------------------------------------------------------------------
! Generates supercells.
! ----------------------------------------------------------------------
function generate_supercells(structure,qpoints,symmetry_precision, &
   & loto_direction) result(output)
  implicit none
  
  ! Inputs.
  type(StructureData),  intent(in)           :: structure
  type(QpointData),     intent(in)           :: qpoints(:)
  real(dp),             intent(in)           :: symmetry_precision
  type(FractionVector), intent(in), optional :: loto_direction
  type(StructureData), allocatable           :: output(:)
  
  ! q-point variables.
  logical, allocatable :: accounted_for(:)
  integer, allocatable :: min_sc_sizes(:)
  type(QpointData)     :: transformed_qpoint
  
  ! Supercell variables.
  integer         :: no_supercells
  type(IntMatrix) :: supercell_matrix
  
  ! Temporary variables.
  integer :: i,j,jp,k,l,lp,ialloc
  
  allocate( accounted_for(size(qpoints)), &
          & min_sc_sizes(size(qpoints)),  &
          & output(size(qpoints)),        &
          & stat=ialloc); call err(ialloc)
  
  ! Find the minimum size of the supercell required to simulate each q-point.
  do i=1,size(qpoints)
    min_sc_sizes(i) = qpoints(i)%min_sc_size()
  enddo
  
  ! Find the minimal set of supercells required to simulate every q-point.
  accounted_for = .false.
  no_supercells = 0
  do while (any(.not. accounted_for))
    ! Find the q-point with the largest min supercell size
    !    which is not yet accounted for.
    i = maxloc(min_sc_sizes, dim=1, mask=.not. accounted_for)
    
    ! Create a supercell for simulating this q-point.
    no_supercells = no_supercells+1
    call print_line('Generating supercell '//no_supercells)
    supercell_matrix = construct_supercell_matrix(qpoints(i), structure)
    output(no_supercells) = construct_supercell( structure,       &
                                               & supercell_matrix )
    call output(no_supercells)%calculate_symmetry( &
                 & symmetry_precision,             &
                 & loto_direction = loto_direction )
    
    ! Find all q-points which can be simulated using this supercell.
    do j=1,size(qpoints)
      if (.not. accounted_for(j)) then
        if (is_int(output(no_supercells)%supercell * qpoints(j)%qpoint)) then
          ! Mark the q-point, j, and its paired q-point, jp, as accounted for.
          jp = first(qpoints%id==qpoints(j)%paired_qpoint_id)
          accounted_for(j)  = .true.
          accounted_for(jp) = .true.
          
          ! Find any symmetrically equivalent q-points,
          !    and mark them as accounted for.
          do k=1,size(structure%symmetries)
            transformed_qpoint = structure%symmetries(k) * qpoints(j)
            do l=1,size(qpoints)
              if (qpoints(l) == transformed_qpoint) then
                ! Mark the q-point, l, and its paired q-point, lp,
                !    as accounted for.
                lp = first(qpoints%id==qpoints(l)%paired_qpoint_id)
                accounted_for(l)  = .true.
                accounted_for(lp) = .true.
              endif
            enddo
          enddo
        endif
      endif
    enddo
    
    ! Check that this supercell does simulate this q-point.
    if (.not. accounted_for(i)) then
      call print_line(CODE_ERROR//': Generated supercell does not match &
         &q-point.')
      call err()
    endif
  enddo
  
  output = output(:no_supercells)
end function
end module
