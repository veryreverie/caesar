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
function calculate_unique_vectors(lattice,centre_on_origin) result(output)
  use linear_algebra_module
  implicit none
  
  ! Lattice vectors are the rows of lattice.
  type(IntMatrix), intent(in)  :: lattice
  logical,         intent(in)  :: centre_on_origin
  type(IntVector), allocatable :: output(:)
  
  integer         :: lattice_size
  type(IntMatrix) :: recip_lattice
  type(IntVector) :: vector
  
  integer         :: i,j,k,l,ialloc
  integer         :: no_vectors
  type(IntVector) :: delta
  integer         :: fractional(3)
  
  if (size(lattice,1)/=3 .or. size(lattice,2)/=3) then
    call print_line('Error: lattice is not 3x3.')
    call err()
  endif
  
  lattice_size = abs(determinant(lattice))
  recip_lattice = transpose(invert_int(lattice))
  
  allocate(output(lattice_size), stat=ialloc); call err(ialloc)
  
  ! Loop over vectors of up to lattice_size in each direction.
  no_vectors = 0
  do i=0,lattice_size-1
    do j=0,lattice_size-1
      do_k : do k=0,lattice_size-1
        vector = [k,j,i]
        
        ! Check if the vector has already been found.
        do l=1,no_vectors
          delta = recip_lattice * (vector-output(l))
          if (all(modulo(int(delta),lattice_size) == 0)) then
            cycle do_k
          endif
        enddo
        
        no_vectors = no_vectors+1
        if (no_vectors>lattice_size) then
          call print_line('Error: more unique vectors found than exist.')
          call err()
        endif
        output(no_vectors) = vector
        
      enddo do_k
    enddo
  enddo
  
  if (no_vectors/=lattice_size) then
    call print_line('Error: only '//no_vectors//' unique vectors were found, &
       &where there should be '//lattice_size//'.')
    call err()
  endif
  
  ! Translate vectors by reciprocal lattice vectors such that the elements
  !    of the fractional co-ordinate lie in [0,1) or [-0.5,0.5), depending
  !    on whether or not centre_on_origin is true.
  do i=1,lattice_size
    
    ! Transform the vector into scaled fractional co-ordinates.
    fractional = int(recip_lattice * output(i))
    
    ! Translate the fractional vector so its elements lie in [0,lattice_size).
    fractional = modulo(fractional, lattice_size)
    
    ! Translate the fractional vector so its elements lie in [-0.5,0.5).
    if (centre_on_origin) then
      do j=1,3
        if (fractional(j) > (lattice_size-1)/2) then
          fractional(j) = fractional(j) - lattice_size
        endif
      enddo
    endif
    
    ! Transform the fractional vector back into original co-ordinates.
    output(i) = transpose(lattice) * vec(fractional) / lattice_size
  enddo
end function

! ----------------------------------------------------------------------
! Main function.
! ----------------------------------------------------------------------
function construct_supercell(structure,supercell_matrix,calculate_symmetries) &
   & result(supercell)
  use linear_algebra_module
  use structure_module
  implicit none
  
  type(StructureData), intent(in)           :: structure
  type(IntMatrix),     intent(in)           :: supercell_matrix
  logical,             intent(in), optional :: calculate_symmetries
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
  call calculate_derived_supercell_quantities(supercell)
  call calculate_derived_atom_quantities(supercell)
  
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
    call calculate_symmetry(supercell)
  endif
end function
end module
