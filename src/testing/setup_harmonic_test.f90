! ======================================================================
! Runs unit tests on setup_harmonic.
! ======================================================================
module setup_harmonic_test_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function setup_harmonic_test_keywords() result(keywords)
  use help_module
  use setup_harmonic_module
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  type(KeywordData), allocatable :: keywords_setup_harmonic(:)
  type(KeywordData)              :: keywords_test(0)
  
  integer :: ialloc
  
  keywords_setup_harmonic = setup_harmonic_keywords()
  
  allocate( keywords(size(keywords_setup_harmonic)+size(keywords_test)), &
          & stat=ialloc); call err(ialloc)
  
  keywords(:size(keywords_setup_harmonic)) = keywords_setup_harmonic
  keywords(size(keywords_setup_harmonic)+1:) = keywords_test
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine setup_harmonic_test(arguments)
  use constants_module, only : identity
  use utils_module, only : make_dft_output_filename, mkdir, l2_norm
  use structure_module
  use unique_directions_module
  use dft_input_file_module
  use dft_output_file_module
  use group_module
  use atom_mapping_module
  use dictionary_module
  use setup_harmonic_module
  use structure_test_module
  use qpoints_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Previous settings.
  type(String)              :: dft_code
  type(String)              :: seedname
  
  type(String), allocatable :: no_supercells_file(:)
  integer                   :: no_supercells
  
  ! Directory and file names.
  type(String) :: wd
  type(String) :: sdir
  
  ! Structure information.
  type(StructureData) :: structure
  type(StructureData) :: structure_grid
  type(StructureData), allocatable :: supercells(:)
  type(StructureData) :: supercell
  
  ! q-point information.
  type(QpointData), allocatable :: qpoints_ibz(:)
  
  ! Grid information.
  integer :: grid(3)
  
  ! Unique direction information.
  type(UniqueDirections) :: unique_directions
  type(Group), allocatable :: symmetry_group(:)
  
  ! Temporary variables.
  integer :: i,j,k,l
  integer :: rotation(3,3)
  integer :: gvector_grid(3)
  integer :: gvector_sc(3)
  integer :: ialloc
  integer :: atom_1,atom_2
  real(dp) :: frac_pos(3)
  real(dp) :: frac_diff(3)
  real(dp), allocatable :: rotations_cart(:,:,:)
  real(dp) :: rotation_identity(3,3)
  
  ! --------------------------------------------------
  ! Read in arguments from user.
  ! --------------------------------------------------
  wd = item(arguments, 'working_directory')
  dft_code = item(arguments, 'dft_code')
  seedname = item(arguments, 'seedname')
  grid = int(split(item(arguments, 'q-point_grid')))
  
  ! --------------------------------------------------
  ! Copy out settings to setup_harmonic.used_settings
  ! --------------------------------------------------
  call write_dictionary_file(arguments, wd//'/setup_harmonic.used_settings')
  
  ! --------------------------------------------------
  ! Run new setup_harmonic.
  ! --------------------------------------------------
  call setup_harmonic(arguments)
  
  ! --------------------------------------------------
  ! Read in outputs.
  ! --------------------------------------------------
  no_supercells_file = read_lines(wd//'/no_sc.dat')
  no_supercells = int(no_supercells_file(1))
  
  structure = read_structure_file(wd//'/structure.dat')
  
  structure_grid = read_structure_file(wd//'/structure_grid.dat')
  qpoints_ibz = read_qpoints_file(wd//'/qpoints_ibz.dat')
  
  allocate(supercells(no_supercells), stat=ialloc); call err(ialloc)
  do i=1,no_supercells
    sdir = wd//'/Supercell_'//i
    supercells(i) = read_structure_file(sdir//'/structure.dat')
  enddo
  
  ! --------------------------------------------------
  ! Check results.
  ! --------------------------------------------------
  
  ! Check all q-points are accounted for, 
  !    and that the listed rotations are correct.
  do_i : do i=1,product(grid)
    do j=1,size(qpoints_ibz)
      gvector_grid = nint(matmul( structure_grid%supercell, &
                                & qpoints_ibz(j)%qpoint))
      do k=1,size(qpoints_ibz(j)%gvectors)
        if (qpoints_ibz(j)%gvectors(k)/=i) then
          cycle
        endif
        
        rotation = structure%rotations(:,:,qpoints_ibz(j)%rotations(k))
        if (any(                                                           &
           & matmul(rotation,structure_grid%gvectors(:,i)) /= gvector_grid &
           & )) then
          call print_line('Error: q-point '//i//' does not correctly map onto &
             &the IBZ by rotation '//qpoints_ibz(j)%rotations(k)//'.')
          call err()
        endif
        
        cycle do_i
      enddo
    enddo
    
    call print_line('Error: q-point '//i//' not accounted for in IBZ.')
    call err()
  enddo do_i
  
  call print_line('All q-points accounted for.')
  
  ! Check that all q-points in the IBZ have been associated with G-vectors in 
  !    supercells.
  do i=1,size(qpoints_ibz)
    gvector_grid = nint(matmul( structure_grid%supercell, &
                        & qpoints_ibz(i)%qpoint))
    
    if (any( abs( matmul( structure_grid%supercell, &
                &         qpoints_ibz(i)%qpoint)    &
                & - gvector_grid                    &
                & ) > 1.0e-10_dp)) then
      call print_line('Error: IBZ q-point '//i//' is not on the specified &
         &q-point grid.')
      call err()
    endif
    
    supercell = supercells(qpoints_ibz(i)%sc_id)
    
    gvector_sc = supercell%gvectors(:,qpoints_ibz(i)%gvector_id)
    
    if (any(                                                    &
       &   matmul(transpose(structure_grid%recip_supercell), gvector_grid) &
       & * supercell%sc_size                                 &
       & /=                                                     &
       &   matmul(transpose(supercell%recip_supercell), gvector_sc)     &
       & * structure_grid%sc_size                               &
       & )) then
      call print_line(gvector_grid)
      call print_line(matmul(transpose(structure_grid%recip_supercell), gvector_grid)*supercell%sc_size)
      call print_line(gvector_sc)
      call print_line(matmul(transpose(supercell%recip_supercell), gvector_grid)*structure_grid%sc_size)
      call print_line('Error: IBZ q-point '//i//' does not match its assigned &
         &supercell and G-vector.')
      call err()
    endif
  enddo
  
  call print_line('All q-points correctly assigned to supercells.')
  
  ! Check that all atoms are placed correctly in each supercell.
  do i=1,no_supercells
    call print_line('Checking atoms in supercell '//i//'.')
    supercell = supercells(i)
    
    ! Check that no atoms are on top of one another.
    do j=1,supercell%no_atoms
      do k=1,j-1
        if (l2_norm(supercell%atoms(:,j)-supercell%atoms(:,k))<0.5) then
          call print_line('Error: atoms '//j//' and '//k//' are within 0.5 &
             &bohr of one another.')
          call err()
        endif
      enddo
    enddo
    
    ! Check that all atoms are within the primitive supercell.
    do j=1,supercell%no_atoms
      frac_pos = matmul(supercell%recip_lattice, supercell%atoms(:,j)) &
             & / supercell%sc_size
      if (any(frac_pos<-1.0e-10_dp) .or. any(frac_pos>1+1.0e-10_dp)) then
        call print_line(frac_pos)
        call print_line('Error: atom '//j//' lies outside of the primitive &
           &supercell.')
        call err()
      endif
    enddo
    
    ! Check that all atoms related by G-vectors are actually related.
    do j=1,structure%no_atoms
      do k=1,supercell%sc_size
        do l=1,k-1
          atom_1 = supercell%rvec_and_prim_to_atom(j,k)
          atom_2 = supercell%rvec_and_prim_to_atom(j,l)
          frac_diff = matmul( structure%recip_lattice,       &
                            &   supercell%atoms(:,atom_1) &
                            & - supercell%atoms(:,atom_2))
          if (any( frac_diff-nint(frac_diff)>1.0e-10_dp )) then
            call print_line('Error: atoms '//atom_1//' and '//atom_2//&
               &' are not correctly related.')
            call err()
          endif
        enddo
      enddo
    enddo
    call print_line('Supercell '//i//' constructed correctly.')
    
    ! Check all rotations.
    rotations_cart = calculate_cartesian_rotations(supercell)
    
    do j=1,supercell%no_symmetries
      rotation_identity = matmul( rotations_cart(:,:,j), &
                                & transpose(rotations_cart(:,:,j)))
      if ( any(rotation_identity-identity>1.0e-5_dp) .and. &
         & any(rotation_identity+identity>1.0e-5_dp)) then
        call print_line('')
        call print_line('Error: rotation '//j//' not a rotation.')
        call print_line('Rotation (frac)')
        call print_line(supercell%rotations(1,:,j))
        call print_line(supercell%rotations(2,:,j))
        call print_line(supercell%rotations(3,:,j))
        call print_line('Translation (frac)')
        call print_line(supercell%translations(:,j))
        call print_line('Rotation (cart)')
        call print_line(rotations_cart(1,:,j))
        call print_line(rotations_cart(2,:,j))
        call print_line(rotations_cart(3,:,j))
        call print_line('Translation (cart)')
        call print_line(matmul( transpose(supercell%lattice), &
                              & supercell%translations(:,j)))
        call print_line('R.R^T')
        call print_line(rotation_identity(1,:))
        call print_line(rotation_identity(2,:))
        call print_line(rotation_identity(3,:))
        call print_line('Lattice^T')
        call print_line(supercell%lattice(:,1))
        call print_line(supercell%lattice(:,2))
        call print_line(supercell%lattice(:,3))
        call print_line('Recip lattice')
        call print_line(supercell%recip_lattice(1,:))
        call print_line(supercell%recip_lattice(2,:))
        call print_line(supercell%recip_lattice(3,:))
        call print_line('Atom positions (cart)')
        do k=1,supercell%no_atoms
          call print_line(supercell%atoms(:,k))
        enddo
        call print_line('Transformed atoms (cart)')
        do k=1,supercell%no_atoms
          call print_line( matmul(rotations_cart(:,:,j),supercell%atoms(:,k)) &
                       & + matmul( transpose(supercell%lattice),              &
                       &           supercell%translations(:,j)))
        enddo
        call print_line('Atom positions (frac)')
        do k=1,supercell%no_atoms
          call print_line(matmul(supercell%recip_lattice,supercell%atoms(:,k)))
        enddo
        call print_line('Transformed atoms (frac)')
        do k=1,supercell%no_atoms
          call print_line( matmul( supercell%recip_lattice, &
                       & matmul(rotations_cart(:,:,j),supercell%atoms(:,k)) &
                       & + matmul( transpose(supercell%lattice),            &
                       &           supercell%translations(:,j))))
        enddo
        call print_line('Transformed frac atoms')
        do k=1,supercell%no_atoms
          call print_line( matmul( supercell%rotations(:,:,j), &
                         &        matmul( supercell%recip_lattice, &
                         &                supercell%atoms(:,k))) &
                         & + supercell%translations(:,j))
        enddo
        call err()
      endif
    enddo
    
    call print_line('Supercell '//i//' rotations correct.')
    
  enddo
  
  call print_line('All supercells constructed correctly.')
  
  ! Check that unique directions are sufficient to map out energy surface.
  do i=1,no_supercells
    call print_line('Checking perturbations in supercell '//i//'.')
    sdir = wd//'/Supercell_'//i
    
    supercell = supercells(i)
    
    unique_directions = read_unique_directions_file( &
       & sdir//'/unique_directions.dat')
    
    symmetry_group = read_group_file(sdir//'/symmetry_group.dat')
    
    ! Check that all atoms are accounted for.
    do_j : do j=1,supercell%no_atoms
      do k=1,size(unique_directions)
        do l=1,size(symmetry_group)
          if (operate(symmetry_group(l),unique_directions%atoms(k))==j) then
            cycle do_j
          endif
        enddo
      enddo
      
      call print_line('Error: atom '//j//' not related to any of the unique &
         &directions.')
      call err()
    enddo do_j
    
    call print_line('Supercell '//i//' perturbations constructed correctly.')
  enddo
  
  call print_line('All supercell perturbations constructed correctly.')
  
end subroutine
end module
