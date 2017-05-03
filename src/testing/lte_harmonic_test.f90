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
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Previous user inputs.
  type(Dictionary) :: setup_harmonic_arguments
  
  ! Directories and files.
  type(String) :: new_wd
  type(String) :: old_wd
  type(String) :: new_sdir
  type(String) :: old_sdir
  type(String) :: old_dir
  type(String) :: lte_dir
  
  ! Setup data
  type(String), allocatable :: no_sc_file(:)
  integer                   :: no_sc
  type(String)              :: dft_code
  type(String)              :: seedname
  
  ! Structure data.
  type(StructureData)   :: structure
  type(StructureData)   :: supercell
  real(dp), allocatable :: rotations_cart(:,:,:)
  
  ! Symmetry data.
  type(Group),           allocatable :: symmetry_group(:)
  type(UniqueDirections)             :: unique_directions
  
  ! lte file data.
  integer                   :: lte_file
  real(dp),     allocatable :: forces(:,:,:)
  type(String), allocatable :: old_lte_file(:)
  integer                   :: new_lte_file
  integer                   :: atoms_start_line
  integer                   :: atoms_end_line
  integer                   :: force_start_line
  integer                   :: force_end_line
  integer                   :: line_no
  
  ! Atom mapping data.
  type(Group)               :: new_to_old
  type(Group)               :: old_to_new
  type(String), allocatable :: atom_file(:)
  integer,      allocatable :: atom(:)
  
  ! Force constant data.
  real(dp),     allocatable :: old_force_constants(:,:,:)
  real(dp),     allocatable :: new_force_constants(:,:,:)
  real(dp)                  :: average
  real(dp)                  :: average_err
  type(String), allocatable :: old_force_constants_file(:)
  integer                   :: new_force_constants_file
  
  ! Dynamical matrix variables.
  type(LteReturn)           :: lte_result
  complex(dp),  allocatable :: new_dyn_mats(:,:,:)
  type(String), allocatable :: old_dyn_mat_file(:)
  type(String), allocatable :: old_line(:)
  complex(dp)               :: old_element
  complex(dp)               :: new_element
  integer                   :: delta_gvec(3)
  real(dp)                  :: arg
  complex(dp)               :: phase
  
  ! Displacement pattern variables.
  type(DispPatterns) :: old_disp_patts
  
  ! Fourier interpolation data.
  integer               :: no_kspace_lines
  real(dp), allocatable :: disp_qpoints(:,:)
  integer               :: grid(3)
  real(dp)              :: temperature
  
  ! q-point data.
  type(StructureData)           :: structure_grid
  type(QpointData), allocatable :: qpoints_ibz(:)
  
  ! Temporary variables.
  integer                   :: i,j,k,l
  integer                   :: atom_1,atom_2
  integer                   :: mode_1,mode_2
  integer                   :: atom_1_sc,atom_2_sc
  integer                   :: atom_1p_sc,atom_2p_sc
  integer                   :: gvec_1p,gvec_2p
  integer                   :: start_line
  type(String), allocatable :: line(:)
  integer                   :: result_code
  
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
  no_sc_file = read_lines(new_wd//'/no_sc.dat')
  no_sc = int(no_sc_file(1))
  
  ! Read in structure.
  structure = read_structure_file(new_wd//'/structure.dat')
  
  ! Read qpoint data.
  structure_grid = read_structure_file(new_wd//'/structure_grid.dat')
  qpoints_ibz = read_qpoints_file(new_wd//'/qpoints_ibz.dat')
  
  ! --------------------------------------------------
  ! Loop across supercells, testing each in turn.
  ! --------------------------------------------------
  old_wd = new_wd//'/old_caesar'
  call mkdir(old_wd)
  
  do i=1,no_sc
    call print_line('')
    call print_line('Checking supercell '//i)
    new_sdir = new_wd//'/Supercell_'//i
    old_sdir = old_wd//'/Supercell_'//i
    call mkdir(old_sdir)
    
    ! Read in supercell structure.
    supercell = read_structure_file(new_sdir//'/structure.dat')
    
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
    do j=1,3
      call print_line(lte_file, structure%lattice(:,j))
    enddo
    call print_line(lte_file, 'Supercell lattice vectors (rows, in a.u.)')
    do j=1,3
      call print_line(lte_file, supercell%lattice(:,j))
    enddo
    call print_line(lte_file, 'Number of atoms in supercell')
    call print_line(lte_file, supercell%no_atoms)
    call print_line(lte_file, 'Species ; mass (a.u.) ; &
       &position of atom in supercell (in terms of SC LVs)')
    do j=1,supercell%no_atoms
      call print_line(lte_file, supercell%species(j)     //' '// &
                              & supercell%mass(j)        //' '// &
                              & matmul( supercell%recip_lattice, &
                              &         supercell%atoms(:,j)))
    enddo
    call print_line(lte_file, 'Number of point-symmetry operations')
    call print_line(lte_file, supercell%no_symmetries)
    call print_line(lte_file, 'Rotation matrices (3 rows) &
       &and translations (1 row, in terms of SC LVs)')
    do j=1,supercell%no_symmetries
      do k=1,3
        call print_line(lte_file, rotations_cart(:,k,j))
      enddo
      call print_line(lte_file, supercell%translations(:,j))
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
    
!    ! Read old lte.dat file.
!    old_lte_file = read_lines(old_dir//'/'//new_sdir//'/lte/lte.dat')
!    do j=1,size(old_lte_file)
!      line = split(old_lte_file(j))
!      if (size(line) >= 3) then
!        if (join(line(1:3))=='Species ; mass') then
!          atoms_start_line = j
!        elseif (join(line(1:3))=='Number of point') then
!          atoms_end_line = j
!        elseif (join(line(1:3))=='Atom 1 ;') then
!          force_start_line = j
!        elseif (join(line(1:3))=='Program function: (1)') then
!          force_end_line = j
!        endif
!      endif
!    enddo
!    
!    ! Write new lte.dat file, rearranged to new atom layout.
!    lte_dir = new_sdir//'/old_lte'
!    call mkdir(lte_dir)
!    new_lte_file = open_write_file(lte_dir//'/lte.dat')
!    
!    do j=1,atoms_start_line
!      call print_line(new_lte_file,old_lte_file(j))
!    enddo
!    
!    do j=atoms_start_line+1,atoms_end_line-1
!      line_no = operate(new_to_old, j-atoms_start_line) + atoms_start_line
!      call print_line(new_lte_file, old_lte_file(line_no))
!    enddo
!    
!    do j=atoms_end_line,force_start_line
!      call print_line(new_lte_file,old_lte_file(j))
!    enddo
!    
!    do j=force_start_line+1,force_end_line-1
!      line = split(old_lte_file(j))
!      line(1) = operate(old_to_new,int(line(1)))
!      line(3) = operate(old_to_new,int(line(3)))
!      call print_line(new_lte_file,join(line))
!    enddo
!    
!    do j=force_end_line,size(old_lte_file)
!      call print_line(new_lte_file,old_lte_file(j))
!    enddo
!    close(new_lte_file)
!    
!    ! Run old lte.
!    result_code = system_call('cd '//lte_dir//'; lte_lower > lte.out')
!    call err(result_code==0)
!    
!    ! Read in atom.dat.
!    atom_file = read_lines(lte_dir//'/atom.dat')
!    allocate(atom(size(atom_file)))
!    do j=1,size(atom_file)
!      atom(j) = int(atom_file(j))
!    enddo
!    
!    ! Read in force constants from old lte.
!    old_force_constants_file = read_lines(lte_dir//'/force_constants.dat')
!    
!    allocate(old_force_constants( structure%no_modes, &
!                                & structure%no_modes, &
!                                & supercell%sc_size))
!    do atom_1=1,structure%no_atoms
!      atom_1_sc = supercell%rvec_and_prim_to_atom(atom_1,1)
!      mode_1 = (atom_1-1)*3+1
!      do atom_2=1,structure%no_atoms
!        mode_2 = (atom_2-1)*3+1
!        do j=1,supercell%sc_size
!          atom_2_sc = supercell%rvec_and_prim_to_atom(atom_2,j)
!          start_line = ( (atom_1-1)*supercell%no_atoms &
!                   &   + (j-1)         &
!                   &   + (atom_2-1)*supercell%sc_size) &
!                   & * 5 + 1
!          
!          line = split(old_force_constants_file(start_line))
!          if (any(int(line)/=[atom_1,atom_2,j])) then
!            call print_line('Force constants file in unexpected order.')
!            call print_line('Line '//start_line//' is '//join(line))
!            call print_line('Expected: '//atom_1//' '//atom_2//' '//j)
!            call err()
!          endif
!          
!          do k=1,3
!            line = split(old_force_constants_file(start_line+k))
!            old_force_constants(mode_2:mode_2+2,mode_1+k-1,j) = dble(line)
!          enddo
!        enddo
!      enddo
!    enddo
!    
!    ! Generate new force constants
!  !  new_force_constants = calculate_force_constants(structure,supercell, &
!  !     & symmetry_group,unique_directions,new_sdir,dft_code,seedname)
!  
!    ! Write out force constants for debugging purposes.
!    call mkdir(new_sdir//'/new_lte')
!    new_force_constants_file = &
!       & open_write_file(new_sdir//'/new_lte/force_constants.dat')
!    do atom_1=1,structure%no_atoms
!      mode_1 = (atom_1-1)*3 + 1
!      do atom_2=1,structure%no_atoms
!        mode_2 = (atom_2-1)*3 + 1
!        do j=1,supercell%sc_size
!          call print_line(new_force_constants_file,atom_1//' '//atom_2//' '//j)
!          do k=1,3
!            call print_line( new_force_constants_file, &
!                           & new_force_constants(mode_2:mode_2+2,mode_1+k-1,j))
!          enddo
!          call print_line(new_force_constants_file,'')
!        enddo
!      enddo
!    enddo
!    
!    ! Check force constants.
!    average = 0.0_dp
!    average_err = 0.0_dp
!    do mode_1=1,structure%no_modes
!      do mode_2=1,structure%no_modes
!        do j=1,supercell%sc_size
!          k = supercell%paired_gvec(j)
!          
!          average = average + old_force_constants(mode_2,mode_1,j)**2
!          
!          average_err = average_err + ( old_force_constants(mode_2,mode_1,j) &
!                                    & - new_force_constants(mode_2,mode_1,j) &
!                                    & ) **2
!          
!          if (abs( old_force_constants(mode_2,mode_1,j) &
!               & - old_force_constants(mode_1,mode_2,k) ) > 1.0e-8_dp) then
!            call print_line('')
!            call print_line('Old force constants are not symmetric.')
!            call print_line('Modes: '//mode_1//' '//mode_2)
!            call print_line('Atoms: '//(mode_1-1)/3+1//' '// &
!               & (mode_2-1)/3+1//' '//j)
!            call print_line('Elements: '//                    &
!               & old_force_constants(mode_2,mode_1,j) //' '// &
!               & old_force_constants(mode_1,mode_2,k))
!            call err()
!          endif
!          
!          if (abs( new_force_constants(mode_2,mode_1,j) &
!               & - new_force_constants(mode_1,mode_2,k) ) > 1.0e-8_dp) then
!            call print_line('')
!            call print_line('New force constants are not symmetric.')
!            call print_line('Modes: '//mode_1//' '//mode_2)
!            call print_line('Atoms: '//(mode_1-1)/3+1//' '// &
!               & (mode_2-1)/3+1//' '//j)
!            call print_line('Elements: '//                    &
!               & new_force_constants(mode_2,mode_1,j) //' '// &
!               & new_force_constants(mode_1,mode_2,k))
!            call err()
!          endif
!          
!          if (abs( old_force_constants(mode_2,mode_1,j) &
!               & - new_force_constants(mode_2,mode_1,j) ) > 1.0e-4_dp) then
!            call print_line('')
!            call print_line('Old force constants do not match &
!               &new force constants.')
!            call print_line('Modes: '//mode_1//' '//mode_2)
!            call print_line('Atoms: '//(mode_1-1)/3+1//' '// &
!               & (mode_2-1)/3+1//' '//j)
!            call print_line('Elements: '//                    &
!               & old_force_constants(mode_2,mode_1,j) //' '// &
!               & new_force_constants(mode_2,mode_1,j))
!            call print_line('')
!          endif
!        enddo
!      enddo
!    enddo
!    
!    average = sqrt( average &
!                & / (structure%no_modes**2*supercell%sc_size))
!    average_err = sqrt( average_err &
!                    & / (structure%no_modes**2*supercell%sc_size))
!    
!    call print_line('Checking force constants.')
!    call print_line('L2 Average force constant: '//average)
!    call print_line('L2 Average error:          '//average_err)
!    
!    ! Run new lte, with old force constants.
!    lte_result  = evaluate_freqs_on_grid( structure,    &
!                                        & supercell, &
!                                        & old_force_constants)
!    
!    deallocate(old_force_constants)
!    deallocate(new_force_constants)
!    
!    allocate(new_dyn_mats( structure%no_modes, &
!                         & structure%no_modes, &
!                         & supercell%sc_size))
!    new_dyn_mats = lte_result%dynamical_matrices
!    
!    ! Check that dynamical matrices are the same.
!    do j=1,supercell%sc_size
!      call print_line('Checking dynamical matrix '//j)
!      old_dyn_mat_file = read_lines(lte_dir//'/dyn_mat.'//j//'.dat')
!      do atom_1=1,structure%no_atoms
!        atom_1_sc = supercell%rvec_and_prim_to_atom(atom_1,1)
!        atom_1p_sc = atom(atom_1_sc)
!        gvec_1p = supercell%atom_to_rvec(atom_1p_sc)
!        do atom_2=1,structure%no_atoms
!          atom_2_sc = supercell%rvec_and_prim_to_atom(atom_2,1)
!          atom_2p_sc = atom(atom_2_sc)
!          gvec_2p = supercell%atom_to_rvec(atom_2p_sc)
!          do k=1,3
!            mode_1 = (atom_1-1)*3 + k
!            do l=1,3
!              mode_2 = (atom_2-1)*3 + l
!              old_line = split(old_dyn_mat_file(   ((atom_1-1)*3+k-1)     &
!                                               & * (structure%no_atoms*3) &
!                                               & + ((atom_2-1)*3+l)))
!              
!              if (any(int(old_line(1:4))/=[atom_1,k,atom_2,l])) then
!                call print_line('Old dyn_mat file in unexpected order.')
!                call err()
!              endif
!              
!              new_element = new_dyn_mats(mode_2,mode_1,j)
!              old_element = cmplx(dble(old_line(5)),dble(old_line(6)),dp)
!              
!              ! Correct by relative phase if old lte misordered atoms.
!              delta_gvec = supercell%gvectors(:,gvec_1p) &
!                       & - supercell%gvectors(:,gvec_2p)
!              arg = dot_product(matmul( &
!                 & delta_gvec, &
!                 & supercell%recip_supercell), & ! TODO: this may be wrong.
!                 & supercell%gvectors(:,j)) &
!                 & *2*pi/supercell%sc_size
!              phase = cmplx(cos(arg),sin(arg),dp)
!              
!              if ( abs(new_element-old_element*phase) > 1.0e-11_dp .and. &
!                 & abs(new_element-old_element/phase) > 1.0e-11_dp ) then
!                call print_line('')
!                call print_line('Dynamical matrices disagree.')
!                call print_line('Supercell '//i)
!                call print_line('Dynamical matrix '//j)
!                call print_line('Atoms '//atom_1//' '//atom_2)
!                call print_line('Directions '//k//' '//l)
!                call print_line('Old values:')
!                call print_line(real(old_element)//' '//aimag(old_element))
!                call print_line( real(old_element*phase)//' '// &
!                               & aimag(old_element*phase))
!                call print_line( real(old_element/phase)//' '// &
!                               & aimag(old_element/phase))
!                call print_line('New values:')
!                call print_line(real(new_element)//' '//aimag(new_element))
!                call err()
!              endif
!            enddo
!          enddo
!        enddo
!      enddo
!    enddo
!    
!    ! Check that displacement patterns are the same.
!    call print_line('Checking displacement patterns.')
!    old_disp_patts = read_disp_patterns_file( &
!       & new_sdir//'/old_lte/disp_patterns.dat',  &
!       & structure%no_modes)
!    
!    do j=1,supercell%sc_size
!      do mode_1=1,structure%no_modes
!        if (abs( old_disp_patts%frequencies(mode_1,j) &
!               & - lte_result%frequencies(mode_1,j)) > 1.0e-7_dp) then
!          call print_line('Frequencies are different.')
!          call print_line('G-vector      : '//j)
!          call print_line('Mode          : '//mode_1)
!          call print_line('Old frequency : '// &
!             & old_disp_patts%frequencies(mode_1,j))
!          call print_line('New frequency : '// &
!             & lte_result%frequencies(mode_1,j))
!          call err()
!        endif
!        
!        do atom_1=1,supercell%no_atoms
!          if (abs( old_disp_patts%prefactors(atom_1,mode_1,j) &
!                 & - lte_result%prefactors(atom_1,mode_1,j)) > 1.0e-10_dp) then
!            call print_line('Prefactors are different.')
!            call print_line('G-vector      : '//j)
!            call print_line('Mode          : '//mode_1)
!            call print_line('Atom          : '//atom_1)
!            call print_line('Old frequency : '// &
!               & old_disp_patts%prefactors(atom_1,mode_1,j))
!            call print_line('New frequency : '// &
!               & lte_result%prefactors(atom_1,mode_1,j))
!            call err()
!          endif
!          
!        !  if (any(abs( old_disp_patts%disp_patterns(:,atom_1,mode_1,j) &
!        !     &       - lte_result%displacements(:,atom_1,mode_1,j) &
!        !     & ) > 1.0e-2_dp)) then
!        !    call print_line('Displacement patterns are different.')
!        !    call print_line('G-vector                 : '//j)
!        !    call print_line('Mode                     : '//mode_1)
!        !    call print_line('Atom                     : '//atom_1)
!        !    call print_line('Old displacement pattern : ')
!        !    call print_line(old_disp_patts%disp_patterns(:,atom_1,mode_1,j))
!        !    call print_line('New displacement pattern : ')
!        !    call print_line(lte_result%displacements(:,atom_1,mode_1,j))
!        !    call err()
!        !  endif
!        enddo
!      enddo
!    enddo
!    
!    deallocate(atom)
!    deallocate(new_dyn_mats)
  enddo
end subroutine
end module
