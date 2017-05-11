! ======================================================================
! Runs unit tests on lte_harmonic.
! ======================================================================
module lte_harmonic_test_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function lte_harmonic_test_keywords() result(keywords)
  use help_module
  use lte_harmonic_module
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  type(KeywordData), allocatable :: keywords_lte_harmonic(:)
  type(KeywordData)              :: keywords_test(0)
  
  integer :: ialloc
  
  keywords_lte_harmonic = lte_harmonic_keywords()
  
  allocate( keywords(size(keywords_lte_harmonic)+size(keywords_test)), &
          & stat=ialloc); call err(ialloc)
  
  keywords(:size(keywords_lte_harmonic)) = keywords_lte_harmonic
  keywords(size(keywords_lte_harmonic)+1:) = keywords_test
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine lte_harmonic_test(arguments)
  use constants_module, only : pi
  use utils_module,     only : mkdir
  use structure_module
  use group_module
  use atom_mapping_module
  use lte_module
  use unique_directions_module
  use displacement_patterns_module
  use lte_harmonic_module
  use qpoints_module
  use dictionary_module
  use linear_algebra_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Previous user inputs.
  type(Dictionary) :: setup_harmonic_arguments
  
  ! Directories and files.
  type(String) :: new_wd
  type(String) :: old_wd
  type(String) :: new_sdir
  type(String) :: old_sdir
  integer      :: lte_file
  type(String) :: qdir
  
  ! Setup data
  type(String), allocatable :: no_supercells_file(:)
  integer                   :: no_supercells
  type(String)              :: dft_code
  type(String)              :: seedname
  
  ! Structure data.
  type(StructureData)              :: structure
  type(StructureData), allocatable :: supercells(:)
  type(StructureData)              :: supercell
  type(UniqueDirections)           :: unique_directions
  
  ! Symmetry data.
  type(Group),      allocatable :: symmetry_group(:)
  type(RealMatrix), allocatable :: rotations_cart(:)
  
  ! R-vector data.
  type(Group), allocatable  :: rvector_group(:)
  
  ! Force constant data.
  real(dp),     allocatable :: forces(:,:,:)
  real(dp),     allocatable :: old_force_constants(:,:,:)
  real(dp),     allocatable :: new_force_constants(:,:,:)
  type(String), allocatable :: old_force_constants_file(:)
  
  ! q-point data.
  type(StructureData)           :: structure_grid
  type(QpointData), allocatable :: qpoints_ibz(:)
  type(RealVector)              :: qpoint
  
  ! Dynamical matrix data.
  integer                  :: gvector
  complex(dp), allocatable :: dynamical_matrix(:,:)
  real(dp)                 :: exponent
  
  ! Displacement pattern variables.
  type(LteReturn), allocatable :: lte_results(:)
  integer :: mode,atom,atom_sc
  
  ! Fourier interpolation data.
  integer  :: grid(3)
  real(dp) :: temperature
  
  ! Atom, mode and R-vector IDs.
  integer :: atom_1,atom_2
  integer :: atom_1p,atom_2p
  integer :: mode_1,mode_2,mode_1p,mode_2p
  integer :: atom_1_sc,atom_2_sc
  integer :: atom_1p_sc,atom_2p_sc
  integer :: rvector,rvector_p,rvector_1p,rvector_2p
  
  ! Errors.
  real(dp) :: new_error, old_error, force
  type(RealMatrix) :: matrix_new, matrix_p_new
  type(RealMatrix) :: matrix_old, matrix_p_old
  type(RealVector) :: x, f, new_f, old_f
  real(dp) :: new_total, old_total
  
  ! Temporary variables.
  integer                   :: i,j,k,l,ialloc
  integer                   :: line_no
  type(String), allocatable :: line(:)
  
  ! --------------------------------------------------
  ! Read in arguments from user.
  ! --------------------------------------------------
  new_wd = item(arguments, 'working_directory')
  temperature = dble(item(arguments, 'temperature'))
  
  ! --------------------------------------------------
  ! Copy out settings to lte_harmonic.used_settings
  ! --------------------------------------------------
  call write_dictionary_file(arguments, new_wd//'/lte_harmonic.used_settings')
  
  ! --------------------------------------------------
  ! Run new lte_harmonic.
  ! --------------------------------------------------
  call lte_harmonic(arguments)
  
  ! --------------------------------------------------
  ! Read in previous arguments.
  ! --------------------------------------------------
  setup_harmonic_arguments = read_dictionary_file( &
     & new_wd//'/setup_harmonic.used_settings')
  dft_code = item(setup_harmonic_arguments, 'dft_code')
  seedname = item(setup_harmonic_arguments, 'seedname')
  grid = int(split(item(setup_harmonic_arguments, 'q-point_grid')))
  
  ! Read in setup data.
  no_supercells_file = read_lines(new_wd//'/no_sc.dat')
  no_supercells = int(no_supercells_file(1))
  
  ! Read in structure.
  structure = read_structure_file(new_wd//'/structure.dat')
  
  ! Read in supercells.
  allocate( supercells(no_supercells), &
          & lte_results(no_supercells), &
          & stat=ialloc); call err(ialloc)
  do i=1,no_supercells
    new_sdir = new_wd//'/Supercell_'//i
    supercells(i) = read_structure_file(new_sdir//'/structure.dat')
  enddo
  
  ! --------------------------------------------------
  ! Loop across supercells, testing each in turn.
  ! --------------------------------------------------
  old_wd = new_wd//'/old_caesar'
  call mkdir(old_wd)
  
  do i=1,no_supercells
    call print_line('')
    call print_line('Checking supercell '//i)
    new_sdir = new_wd//'/Supercell_'//i
    old_sdir = old_wd//'/Supercell_'//i
    call mkdir(old_sdir)
    
    supercell = supercells(i)
    
    rotations_cart = calculate_cartesian_rotations(supercell)
    
    ! Read in symmetry group and unique atoms.
    symmetry_group = read_group_file(new_sdir//'/symmetry_group.dat')
    unique_directions = read_unique_directions_file( &
       & new_sdir//'/unique_directions.dat')
    
    ! Read in forces.
    forces = read_forces(supercell, unique_directions, new_sdir, dft_code, &
       & seedname)
    
    ! Write lte.dat to run old lte.
    lte_file = open_write_file(old_sdir//'/lte.dat')
    call print_line(lte_file, 'Primitive lattice vectors (rows, in a.u.)')
    call print_line(lte_file, structure%lattice)
    call print_line(lte_file, 'Supercell lattice vectors (rows, in a.u.)')
    call print_line(lte_file, supercell%lattice)
    call print_line(lte_file, 'Number of atoms in supercell')
    call print_line(lte_file, supercell%no_atoms)
    call print_line(lte_file, 'Species ; mass (a.u.) ; &
       &position of atom in supercell (in terms of SC LVs)')
    do j=1,supercell%no_atoms
      call print_line(lte_file, supercell%species(j) //' '// &
                              & supercell%mass(j)    //' '// &
                              & supercell%recip_lattice * supercell%atoms(j))
    enddo
    call print_line(lte_file, 'Number of point-symmetry operations')
    call print_line(lte_file, supercell%no_symmetries)
    call print_line(lte_file, 'Rotation matrices (3 rows) &
       &and translations (1 row, in terms of SC LVs)')
    do j=1,supercell%no_symmetries
      call print_line(lte_file, rotations_cart(j))
      call print_line(lte_file, supercell%translations(j))
    enddo
    call print_line(lte_file, 'Number of force constants supplied')
    call print_line(lte_file, 3*supercell%no_atoms*size(unique_directions))
    call print_line(lte_file, 'Atom 1 ; Cartesian direction ; &
       &Atom 2 ; Cartesian direction ; force constants (a.u.)')
    do j=1,size(unique_directions)
      do k=1,supercell%no_atoms
        do l=1,3
          call print_line(lte_file, unique_directions%atoms(j)         //' '//&
                                  & unique_directions%directions_int(j)//' '//&
                                  & k                                  //' '//&
                                  & l                                  //' '//&
                                  & forces(l,k,j))
        enddo
      enddo
    enddo
    call print_line(lte_file, 'Program fn: (4) evaluate freqs on grid')
    call print_line(lte_file, 4)
    call print_line(lte_file, 'Temperature (K)')
    call print_line(lte_file, 0)
    call print_line(lte_file, 'Number of lines in k space to plot')
    call print_line(lte_file, 0)
    call print_line(lte_file, 'Points on journey through k space')
    call print_line(lte_file, '0 0 0')
    close(lte_file)
    
    ! --------------------------------------------------
    ! Run old lte.
    ! --------------------------------------------------
    call execute_old_code(old_sdir, str('lte_lower > lte.out'))
    
    ! --------------------------------------------------
    ! Check force constants.
    ! --------------------------------------------------
    ! Mass reduce forces.
    do j=1,size(forces,3)
      do k=1,size(forces,2)
        forces(:,k,j) = forces(:,k,j)            &
                    & / sqrt ( supercell%mass(k) &
                    &        * supercell%mass(unique_directions%atoms(j)))
      enddo
    enddo
    
    ! Calculate force constants.
    new_force_constants = construct_force_constants(forces,supercell,&
       & unique_directions,symmetry_group)
    
    ! Read in old force constants.
    allocate( old_force_constants( supercell%sc_size,   &
            &                      structure%no_modes,  &
            &                      structure%no_modes), &
            & stat=ialloc); call err(ialloc)
    old_force_constants_file = read_lines(old_sdir//'/force_constants.dat')
    if ( size(old_force_constants_file) &
       & /= (supercell%sc_size+2)*structure%no_modes**2) then
      call err()
    endif
    do j=1,structure%no_modes
      do k=1,structure%no_modes
        do l=1,supercell%sc_size
          line_no = (supercell%sc_size+2)*(structure%no_modes*(j-1)+k-1)+l+1
          line = split(old_force_constants_file(line_no))
          old_force_constants(l,k,j) = dble(line(5))
        enddo
      enddo
    enddo
    
    ! Calculate L2 error between forces and force constants.
    new_error = 0.0_dp
    old_error = 0.0_dp
    force = 0.0_dp
    do j=1,supercell%no_symmetries
      do k=1,size(unique_directions)
        
        if (unique_directions%directions_char(k)=='x') then
          x = [1.0_dp, 0.0_dp, 0.0_dp]
        elseif (unique_directions%directions_char(k)=='y') then
          x = [0.0_dp, 1.0_dp, 0.0_dp]
        else
          x = [0.0_dp, 0.0_dp, 1.0_dp]
        endif
        x = rotations_cart(j) * x
        
        do l=1,supercell%no_atoms
          atom_1 = operate(symmetry_group(j),unique_directions%atoms(k))
          if (supercell%atom_to_rvec(atom_1)/=1) then
            cycle
          endif
          atom_2 = operate(symmetry_group(j),l)
          mode_1 = (supercell%atom_to_prim(atom_1)-1)*3+1
          mode_2 = (supercell%atom_to_prim(atom_2)-1)*3+1
          rvector = supercell%atom_to_rvec(atom_2)
          f = rotations_cart(j) * vec(forces(:,l,k))
          new_f = mat( new_force_constants( rvector,          &
              &                             mode_2:mode_2+2,  &
              &                             mode_1:mode_1+2)) &
              & * x
          old_f = mat( old_force_constants( rvector,          &
              &                             mode_2:mode_2+2,  &
              &                             mode_1:mode_1+2)) &
              & * x
          new_error = new_error + (new_f-f)*(new_f-f)
          old_error = old_error + (old_f-f)*(old_f-f)
          force = force + f*f
        enddo
      enddo
    enddo
    new_error = sqrt(new_error/force)
    old_error = sqrt(old_error/force)
    call print_line( &
       & 'Fractional L2 difference between symmetrised and raw forces:')
    call print_line('New code: '//new_error)
    call print_line('Old code: '//old_error)
    
    rvector_group = calculate_rvector_group(supercell)
    
    ! Calculate L2 error between symmetries.
    new_error = 0.0_dp
    old_error = 0.0_dp
    new_total = 0.0_dp
    old_total = 0.0_dp
    do j=1,supercell%no_symmetries
      do atom_1=1,structure%no_atoms
        mode_1 = (atom_1-1)*3+1
        do atom_2=1,structure%no_atoms
          mode_2 = (atom_2-1)*3+1
          do rvector=1,supercell%sc_size
            atom_1_sc = supercell%rvec_and_prim_to_atom(atom_1,1)
            atom_2_sc = supercell%rvec_and_prim_to_atom(atom_2,rvector)
            atom_1p_sc = operate(symmetry_group(j),atom_1_sc)
            atom_2p_sc = operate(symmetry_group(j),atom_2_sc)
            atom_1p = supercell%atom_to_prim(atom_1p_sc)
            atom_2p = supercell%atom_to_prim(atom_2p_sc)
            mode_1p = (atom_1p-1)*3+1
            mode_2p = (atom_2p-1)*3+1
            rvector_1p = supercell%atom_to_rvec(atom_1p_sc)
            rvector_2p = supercell%atom_to_rvec(atom_2p_sc)
            rvector_p = operate( &
               & rvector_group(supercell%paired_rvec(rvector_1p)), rvector_2p)
            matrix_new = new_force_constants( rvector,         &
                                            & mode_2:mode_2+2, &
                                            & mode_1:mode_1+2)
            matrix_old = old_force_constants( rvector,         &
                                            & mode_2:mode_2+2, &
                                            & mode_1:mode_1+2)
            matrix_p_new = new_force_constants( rvector_p,         &
                                              & mode_2p:mode_2p+2, &
                                              & mode_1p:mode_1p+2)
            matrix_p_old = old_force_constants( rvector_p,         &
                                              & mode_2p:mode_2p+2, &
                                              & mode_1p:mode_1p+2)
            matrix_p_new =  transpose(rotations_cart(j)) &
                       & * matrix_p_new                  &
                       & * rotations_cart(j)
            matrix_p_old =  transpose(rotations_cart(j)) &
                       & * matrix_p_old                  &
                       & * rotations_cart(j)
            new_total = new_total + sum(dble(matrix_new)**2)
            old_total = old_total + sum(dble(matrix_old)**2)
            new_error = new_error + sum(dble(matrix_new-matrix_p_new)**2)
            old_error = old_error + sum(dble(matrix_old-matrix_p_old)**2)
          enddo
        enddo
      enddo
    enddo
    
    if (new_total < 1.0e-10_dp .and. new_error > 1.0e-10_dp) then
      call print_line('Error in symmetrisation of force constants.')
      call err()
    endif
    
    new_error = sqrt(new_error/new_total)
    old_error = sqrt(old_error/old_total)
    
    if (new_total > 1.0e-10_dp .and. new_error > 1.0e-10_dp) then
      call print_line('Error in symmetrisation of force constants.')
      call err()
    endif
    
    call print_line( &
       & 'Fractional L2 error in symmetrisation of force constants:')
    call print_line('New code: '//new_error)
    call print_line('Old code: '//old_error)
    
    deallocate(old_force_constants, stat=ialloc); call err(ialloc)
    
    ! Run lte on each supercell.
    lte_results(i) = evaluate_freqs_on_grid(supercell, new_force_constants)
    
    ! Check dynamical matrices.
    do gvector=1,supercell%sc_size
      qpoint = transpose(supercell%recip_supercell) &
           & * supercell%gvectors(gvector) / dble(supercell%sc_size)
      allocate( dynamical_matrix( supercell%no_modes_prim,  &
              &                   supercell%no_modes_prim), &
              & stat=ialloc); call err(ialloc)
      dynamical_matrix = cmplx(0.0_dp, 0.0_dp, dp)
      do rvector=1,supercell%sc_size
        exponent = -2*pi*qpoint*supercell%rvectors(rvector)
        dynamical_matrix = dynamical_matrix                 &
                       & + new_force_constants(rvector,:,:) &
                       & * cmplx(cos(exponent),sin(exponent),dp)
      enddo
      
      new_error = sum(abs ( dynamical_matrix &
                        & - lte_results(i)%dynamical_matrices(:,:,gvector))**2)
      if (new_error > 1.0e-10_dp) then
        new_error = sqrt(new_error / sum(abs(dynamical_matrix)**2))
        if (new_error > 1.0e-10_dp) then
          call print_line('G-vector '//gvector//' dynamical matrix incorrect.')
          call err()
        endif
      endif
      call print_line('G-vector '//gvector//' dynamical matrix correct.')
      
      ! Check polarisation vectors and frequencies.
      dynamical_matrix = cmplx(0.0_dp, 0.0_dp, dp)
      do mode=1,supercell%no_modes_prim
        do atom=1,supercell%no_atoms_prim
          do rvector=1,supercell%sc_size
            atom_sc = supercell%rvec_and_prim_to_atom(atom,rvector)
            ! TODO: fill this out.
          enddo
        enddo
      enddo
      
      deallocate(dynamical_matrix, stat=ialloc); call err(ialloc)
    enddo
  enddo
  
  call print_line('')
  call print_line('All force constants and dynamical matrices constructed.')
  
  ! --------------------------------------------------
  ! Loop across q-points, testing each in turn.
  ! --------------------------------------------------
  
  ! Read qpoint data.
  structure_grid = read_structure_file(new_wd//'/structure_grid.dat')
  qpoints_ibz = read_qpoints_file(new_wd//'/qpoints_ibz.dat')
  
  do i=1,size(qpoints_ibz)
    qdir = new_wd//'qpoint_'//i
    
    ! Read in frequencies.
    
    ! Read in displacements.
    
    ! Reconstruct real part of 
  enddo
end subroutine
end module
