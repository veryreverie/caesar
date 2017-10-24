! ======================================================================
! Constructs structure data given the primitive structure and a supercell.
! ======================================================================
module construct_supercell_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Calculates the set of vectors which are not related to one another by
!    lattice vectors.
! ----------------------------------------------------------------------
! This is used to either calculate the R-vectors of the primitive cell which
!    are unique in the supercell, or the G-vectors of the reciprocal supercell
!    which are unique in the reciprocal primitive cell.
function calculate_unique_vectors(lattice,centre_on_origin) result(output)
  use linear_algebra_module
  implicit none
  
  ! Lattice vectors are the rows of lattice.
  type(IntMatrix), intent(in)  :: lattice
  logical,         intent(in)  :: centre_on_origin
  type(IntVector), allocatable :: output(:)
  
  integer         :: lattice_size
  integer         :: no_vectors
  type(IntVector) :: frac_vec
  type(IntVector) :: prim_vec
  integer         :: i,j,k,ialloc
  
  if (size(lattice,1)/=3 .or. size(lattice,2)/=3) then
    call print_line('Error: lattice is not 3x3.')
    call err()
  endif
  
  lattice_size = abs(determinant(lattice))
  
  allocate(output(lattice_size), stat=ialloc); call err(ialloc)
  
  no_vectors = 0
  do i=0,lattice_size-1
    do j=0,lattice_size-1
      do k=0,lattice_size-1
        
        ! Construct vectors in scaled fractional primitive co-ordinates.
        ! (scaled by lattice_size, to preserve integer representation).
        if (centre_on_origin) then
          frac_vec = [i-lattice_size/2,j-lattice_size/2,k-lattice_size/2]
        else
          frac_vec = [i,j,k]
        endif
        
        ! Transform to scaled fractional supercell co-ordinates.
        prim_vec = transpose(lattice)*frac_vec
        
        ! Check if the scaled co-ordinate is scaling*(integer co-ordinate).
        if (all(modulo(int(prim_vec),lattice_size)==0)) then
          no_vectors = no_vectors+1
          output(no_vectors) = prim_vec / lattice_size
        endif
      enddo
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Main function.
! ----------------------------------------------------------------------
function construct_supercell(structure,supercell_matrix,calculate_symmetries, &
   & symmetry_precision) result(supercell)
  use linear_algebra_module
  use structure_module
  implicit none
  
  type(StructureData), intent(in)           :: structure
  type(IntMatrix),     intent(in)           :: supercell_matrix
  logical,             intent(in), optional :: calculate_symmetries
  real(dp),            intent(in), optional :: symmetry_precision
  type(StructureData)                       :: supercell
  
  ! Atomic positions.
  type(RealVector) :: atom_pos_prim ! Fractional primitive cell co-ordinates.
  type(RealVector) :: atom_pos_sc   ! Scaled fractional supercell co-ordinates.
  type(IntVector)  :: atom_floor_sc ! Scaled fractional supercell co-ordinates.
  type(RealVector) :: copy_pos_prim ! Fractional primitive cell co-ordinates.
  type(RealVector) :: copy_pos_cart ! Cartesian co-ordinates.
  
  ! R-vector positions.
  type(IntVector) :: rvector_sc   ! Scaled fractional supercell co-ordinates.
  type(IntVector) :: rvector_prim ! Fractional primitive cell co-ordinates.
  
  ! Temporary variables
  integer :: i,j
  integer :: atom_counter
  integer :: no_atoms_sc
  logical :: calculate_symmetries_flag
  integer :: sc_size
  
  if (present(calculate_symmetries)) then
    calculate_symmetries_flag = calculate_symmetries
  else
    calculate_symmetries_flag = .true.
  endif
  
  ! Generate supercell and lattice data.
  sc_size = abs(determinant(supercell_matrix))
  no_atoms_sc = structure%no_atoms*sc_size
  call new(supercell, no_atoms_sc, 0, sc_size)
  supercell%lattice = supercell_matrix * structure%lattice
  supercell%supercell = supercell_matrix
  supercell%rvectors = calculate_unique_vectors(supercell_matrix, .false.)
  supercell%gvectors = calculate_unique_vectors( transpose(supercell_matrix), &
                                               & .true.)
  call supercell%calculate_derived_quantities()
  
  ! Generate atomic positions.
  do i=1,structure%no_atoms
    ! N.B. fractional supercell coordinates are scaled so that the supercell
    !    is a cartesian cube with side length sc_size.
    
    ! Transform atom i into fractional primitive cell co-ordinates,
    !    and translate it into the first primitive cell.
    atom_pos_prim = modulo( dble(structure%recip_lattice*structure%atoms(i)), &
                          & 1.0_dp)
    
    ! Calculate the position of atom i in scaled fractional supercell co-ords.
    atom_pos_sc = supercell%recip_supercell * atom_pos_prim
    
    ! Calculate the R-vector corresponding to the first atom.
    atom_floor_sc = floor( dble(atom_pos_sc)/supercell%sc_size ) &
                & * supercell%sc_size
    
    ! Loop accross R-vectors in the supercell, creating a copy of the atom i in
    !    the lattice copy corresponding to each R-vector.
    do j=1,supercell%sc_size
      rvector_prim = supercell%rvectors(j)
      
      ! Convert the R-vector from fractional primitive cell co-ordinates into
      !    scaled fractional supercell co-ordinates.
      rvector_sc = supercell%recip_supercell * rvector_prim
      
      ! Translate the R-vector by supercell lattice vectors, s.t. the copy
      !    lies inside the primitive supercell.
      rvector_sc = modulo(int(rvector_sc+atom_floor_sc),supercell%sc_size) &
               & - int(atom_floor_sc)
      
      ! Convert the R-vector back into fractional primitive cell co-ordinates.
      rvector_prim = transpose(supercell%supercell) * rvector_sc &
                 & / supercell%sc_size
      
      ! Calculate the position of the copy in cartesian co-ordinates.
      copy_pos_prim = rvector_prim + atom_pos_prim
      
      copy_pos_cart = transpose(structure%lattice) * copy_pos_prim
      
      ! Add the copy to the supercell.
      atom_counter = supercell%rvec_and_prim_to_atom(i,j)
      supercell%atoms(atom_counter) = copy_pos_cart
      supercell%mass(atom_counter) = structure%mass(i)
      supercell%species(atom_counter) = structure%species(i)
    enddo
  enddo
  
  if (calculate_symmetries_flag) then
    if (.not. present(symmetry_precision)) then
      call print_line('Symmetry requested but no precision specified.')
      call err()
    endif
    call supercell%calculate_symmetry(symmetry_precision)
  endif
end function
end module
