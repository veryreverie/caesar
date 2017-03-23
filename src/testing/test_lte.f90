module test_lte_module
contains
subroutine test_lte()
  use utils, only : mkdir
  use err_module
  use string_module
  use file_module
  use structure_module
  use group_module
  use atom_mapping_module
  use lte_module
  use fourier_interpolation_module
  use unique_directions_module
  
  use lte_harmonic_module, only : calculate_force_constants
  implicit none
  
  ! Directories and files.
  type(String)              :: sdir
  type(String)              :: copy_dir
  type(String)              :: lte_dir
  type(String), allocatable :: grid_file(:)
  type(String), allocatable :: ibz_file(:)
  
  ! Setup data
  type(String), allocatable :: no_sc_file(:)
  integer                   :: no_sc
  type(String), allocatable :: user_input_file(:)
  type(String)              :: dft_code
  type(String)              :: seedname
  
  ! Structure data.
  type(StructureData) :: structure
  type(StructureData) :: structure_sc
  type(StructureData) :: structure_sc_old
  
  ! Symmetry data.
  type(Group),           allocatable :: symmetry_group(:)
  type(UniqueDirections)             :: unique_directions
  
  ! lte file data.
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
  complex(dp),  allocatable :: dyn_mats_ibz(:,:,:)
  type(String), allocatable :: old_dyn_mat_file(:)
  type(String), allocatable :: old_line(:)
  complex(dp)               :: old_element
  complex(dp)               :: new_element
  integer                   :: delta_gvec(3)
  real(dp)                  :: arg
  complex(dp)               :: phase
  
  ! Fourier interpolation data.
  integer               :: no_kspace_lines
  real(dp), allocatable :: disp_kpoints(:,:)
  integer               :: grid(3)
  real(dp)              :: temperature
  integer               :: no_kpoints
  integer,  allocatable :: kpoints(:,:)
  integer,  allocatable :: multiplicity(:)
  integer,  allocatable :: sc_ids(:)
  integer,  allocatable :: gvector_ids(:)
  type(Group), allocatable :: prim_symmetry_group(:)
  
  ! Temporary variables.
  integer                   :: i,j,k,l
  integer                   :: atom_1,atom_2
  integer                   :: mode_1,mode_2
  integer                   :: atom_1_sc,atom_2_sc
  integer                   :: atom_1p_sc,atom_2p_sc
  integer                   :: gvec_1p,gvec_2p
  integer                   :: start_line
  type(String), allocatable :: line(:)
  
  ! ----------------------------------------------------------------------
  ! Read in settings from user.
  ! ----------------------------------------------------------------------
  ! Read in directory to compare against and copy from.
  call print_line('')
  call print_line('This test will check lte against a previous &
     &calculations.')
  call print_line('')
  call print_line('Where is the harmonic directory for comparison?')
  copy_dir = read_line_from_user()
  call print_line('')
  
  ! Get temperature from user.
  call print_line('What temperature (K)?')
  temperature = dble(read_line_from_user())
  call print_line('')
  
  ! ----------------------------------------------------------------------
  ! Read in initial data.
  ! ----------------------------------------------------------------------
  ! Read in setup data.
  no_sc_file = read_lines('no_sc.dat')
  no_sc = int(no_sc_file(1))
  
  user_input_file = read_lines('user_input.txt')
  dft_code = user_input_file(1)
  seedname = user_input_file(2)
  
  ! Read in structure.
  structure = read_structure_file('structure.dat')
  
  ! Read grid file
  grid_file = read_lines('grid.dat')
  grid = int(split(grid_file(1)))
  
  ! Read kpoints from ibz.dat
  ibz_file = read_lines('ibz.dat')
  no_kpoints = size(ibz_file)
  allocate(kpoints(3,no_kpoints))
  allocate(multiplicity(no_kpoints))
  allocate(sc_ids(no_kpoints))
  allocate(gvector_ids(no_kpoints))
  do i=1,no_kpoints
    line = split(ibz_file(i))
    kpoints(:,i) = int(line(1:3))
    multiplicity(i) = int(line(4))
    sc_ids(i) = int(line(5))
    gvector_ids(i) = int(line(6))
  enddo
  
  allocate(dyn_mats_ibz(structure%no_modes,structure%no_modes,no_kpoints))
  
  ! ----------------------------------------------------------------------
  ! Loop across supercells, testing each in turn.
  ! ----------------------------------------------------------------------
  do i=1,no_sc
    call print_line('')
    call print_line('Checking supercell '//i)
    sdir = 'Supercell_'//i
    
    ! Read in supercell structure.
    structure_sc = read_structure_file(sdir//'/structure.dat')
    structure_sc_old = read_structure_file(copy_dir//'/'//sdir//'/structure.dat')
    
    ! Read in symmetry group and unique atoms.
    symmetry_group = read_group_file(sdir//'/symmetry_group.dat')
    unique_directions = read_unique_directions_file( &
       & sdir//'/unique_directions.dat')
    
    ! Calculate mapping between structures.
    new_to_old = atom_mapping(structure_sc,structure_sc_old)
    old_to_new = atom_mapping(structure_sc_old,structure_sc)
    
    ! Read old lte.dat file
    old_lte_file = read_lines(copy_dir//'/'//sdir//'/lte/lte.dat')
    do j=1,size(old_lte_file)
      line = split(old_lte_file(j))
      if (size(line) >= 3) then
        if (join(line(1:3))=='Species ; mass') then
          atoms_start_line = j
        elseif (join(line(1:3))=='Number of point') then
          atoms_end_line = j
        elseif (join(line(1:3))=='Atom 1 ;') then
          force_start_line = j
        elseif (join(line(1:3))=='Program function: (1)') then
          force_end_line = j
        endif
      endif
    enddo
    
    ! Write new lte.dat file, rearranged to new atom layout.
    lte_dir = sdir//'/old_lte'
    call mkdir(lte_dir)
    new_lte_file = open_write_file(lte_dir//'/lte.dat')
    
    do j=1,atoms_start_line
      call print_line(new_lte_file,old_lte_file(j))
    enddo
    
    do j=atoms_start_line+1,atoms_end_line-1
      line_no = operate(new_to_old, j-atoms_start_line) + atoms_start_line
      call print_line(new_lte_file, old_lte_file(line_no))
    enddo
    
    do j=atoms_end_line,force_start_line
      call print_line(new_lte_file,old_lte_file(j))
    enddo
    
    do j=force_start_line+1,force_end_line-1
      line = split(old_lte_file(j))
      line(1) = operate(old_to_new,int(line(1)))
      line(3) = operate(old_to_new,int(line(3)))
      call print_line(new_lte_file,join(line))
    enddo
    
    do j=force_end_line,size(old_lte_file)
      call print_line(new_lte_file,old_lte_file(j))
    enddo
    close(new_lte_file)
    
    ! Run old lte.
    call system('cd '//lte_dir//'; lte_lower > lte.out')
    
    ! Read in atom.dat.
    atom_file = read_lines(lte_dir//'/atom.dat')
    allocate(atom(size(atom_file)))
    do j=1,size(atom_file)
      atom(j) = int(atom_file(j))
    enddo
    
    ! Read in force constants from old lte.
    old_force_constants_file = read_lines(lte_dir//'/force_constants.dat')
    
    allocate(old_force_constants( structure%no_modes, &
                                & structure%no_modes, &
                                & structure_sc%sc_size))
    do atom_1=1,structure%no_atoms
      atom_1_sc = structure_sc%gvec_and_prim_to_atom(atom_1,1)
      mode_1 = (atom_1-1)*3+1
      do atom_2=1,structure%no_atoms
        mode_2 = (atom_2-1)*3+1
        do j=1,structure_sc%sc_size
          atom_2_sc = structure_sc%gvec_and_prim_to_atom(atom_2,j)
          start_line = ( (atom_1-1)*structure_sc%no_atoms &
                   &   + (j-1)         &
                   &   + (atom_2-1)*structure_sc%sc_size) &
                   & * 5 + 1
          
          line = split(old_force_constants_file(start_line))
          if (any(int(line)/=(/atom_1,atom_2,j/))) then
            call print_line('Force constants file in unexpected order.')
            call print_line('Line '//start_line//' is '//join(line))
            call print_line('Expected: '//atom_1//' '//atom_2//' '//j)
            call err()
          endif
          
          do k=1,3
            line = split(old_force_constants_file(start_line+k))
            old_force_constants(mode_2:mode_2+2,mode_1+k-1,j) = dble(line)
          enddo
        enddo
      enddo
    enddo
    
    ! Generate new force constants
    new_force_constants = calculate_force_constants(structure,structure_sc, &
       & symmetry_group,unique_directions,sdir,dft_code,seedname)
  
    ! Write out force constants for debugging purposes.
    call mkdir(sdir//'/new_lte')
    new_force_constants_file = &
       & open_write_file(sdir//'/new_lte/force_constants.dat')
    do atom_1=1,structure%no_atoms
      mode_1 = (atom_1-1)*3 + 1
      do atom_2=1,structure%no_atoms
        mode_2 = (atom_2-1)*3 + 1
        do j=1,structure_sc%sc_size
          call print_line(new_force_constants_file,atom_1//' '//atom_2//' '//j)
          do k=1,3
            call print_line( new_force_constants_file, &
                           & new_force_constants(mode_2:mode_2+2,mode_1+k-1,j))
          enddo
          call print_line(new_force_constants_file,'')
        enddo
      enddo
    enddo
    
    ! Check force constants.
    average = 0.0_dp
    average_err = 0.0_dp
    do mode_1=1,structure%no_modes
      do mode_2=1,structure%no_modes
        do j=1,structure_sc%sc_size
          k = structure_sc%paired_gvec(j)
          
          average = average + old_force_constants(mode_2,mode_1,j)**2
          
          average_err = average_err + ( old_force_constants(mode_2,mode_1,j) &
                                    & - new_force_constants(mode_2,mode_1,j) &
                                    & ) **2
          
          if (abs( old_force_constants(mode_2,mode_1,j) &
               & - old_force_constants(mode_1,mode_2,k) ) > 1.0e-8_dp) then
            call print_line('')
            call print_line('Old force constants are not symmetric.')
            call print_line('Modes: '//mode_1//' '//mode_2)
            call print_line('Atoms: '//(mode_1-1)/3+1//' '// &
               & (mode_2-1)/3+1//' '//j)
            call print_line('Elements: '//                    &
               & old_force_constants(mode_2,mode_1,j) //' '// &
               & old_force_constants(mode_1,mode_2,k))
            call err()
          endif
          
          if (abs( new_force_constants(mode_2,mode_1,j) &
               & - new_force_constants(mode_1,mode_2,k) ) > 1.0e-8_dp) then
            call print_line('')
            call print_line('New force constants are not symmetric.')
            call print_line('Modes: '//mode_1//' '//mode_2)
            call print_line('Atoms: '//(mode_1-1)/3+1//' '// &
               & (mode_2-1)/3+1//' '//j)
            call print_line('Elements: '//                    &
               & new_force_constants(mode_2,mode_1,j) //' '// &
               & new_force_constants(mode_1,mode_2,k))
            call err()
          endif
          
          if (abs( old_force_constants(mode_2,mode_1,j) &
               & - new_force_constants(mode_2,mode_1,j) ) > 1.0e-4_dp) then
            call print_line('')
            call print_line('Old force constants do not match &
               &new force constants.')
            call print_line('Modes: '//mode_1//' '//mode_2)
            call print_line('Atoms: '//(mode_1-1)/3+1//' '// &
               & (mode_2-1)/3+1//' '//j)
            call print_line('Elements: '//                    &
               & old_force_constants(mode_2,mode_1,j) //' '// &
               & new_force_constants(mode_2,mode_1,j))
            call print_line('')
          endif
        enddo
      enddo
    enddo
    
    average = dsqrt( average &
                 & / (structure%no_modes**2*structure_sc%sc_size))
    average_err = dsqrt( average_err &
                     & / (structure%no_modes**2*structure_sc%sc_size))
    
    call print_line('Checking force constants.')
    call print_line('L2 Average force constant: '//average)
    call print_line('L2 Average error:          '//average_err)
    
    ! Run new lte, with old force constants.
    lte_result  = evaluate_freqs_on_grid( &
                       & structure,                           &
                       & structure_sc,                        &
                       & old_force_constants,                 &
                       & sdir//'/new_lte/freq_grids.dat',     &
                       & sdir//'/new_lte/disp_patterns.dat',  &
                       & sdir//'/new_lte/kdisp_patterns.dat', &
                       & sdir//'/new_lte/pol_vec.dat')
    
    deallocate(old_force_constants)
    deallocate(new_force_constants)
    
    allocate(new_dyn_mats( structure%no_modes, &
                         & structure%no_modes, &
                         & structure_sc%sc_size))
    new_dyn_mats = lte_result%dynamical_matrices
    
    ! Check that dynamical matrices are the same.
    do j=1,structure_sc%sc_size
      call print_line('Checking dynamical matrix '//j)
      old_dyn_mat_file = read_lines(lte_dir//'/dyn_mat.'//j//'.dat')
      do atom_1=1,structure%no_atoms
        atom_1_sc = structure_sc%gvec_and_prim_to_atom(atom_1,1)
        atom_1p_sc = atom(atom_1_sc)
        gvec_1p = structure_sc%atom_to_gvec(atom_1p_sc)
        do atom_2=1,structure%no_atoms
          atom_2_sc = structure_sc%gvec_and_prim_to_atom(atom_2,1)
          atom_2p_sc = atom(atom_2_sc)
          gvec_2p = structure_sc%atom_to_gvec(atom_2p_sc)
          do k=1,3
            mode_1 = (atom_1-1)*3 + k
            do l=1,3
              mode_2 = (atom_2-1)*3 + l
              old_line = split(old_dyn_mat_file(   ((atom_1-1)*3+k-1)     &
                                               & * (structure%no_atoms*3) &
                                               & + ((atom_2-1)*3+l)))
              
              if (any(int(old_line(1:4))/=(/atom_1,k,atom_2,l/))) then
                call print_line('Old dyn_mat file in unexpected order.')
                call err()
              endif
              
              new_element = new_dyn_mats(mode_2,mode_1,j)
              old_element = cmplx(dble(old_line(5)),dble(old_line(6)),dp)
              
              ! Correct by relative phase if old lte misordered atoms.
              delta_gvec = structure_sc%gvectors(:,gvec_1p) &
                       & - structure_sc%gvectors(:,gvec_2p)
              arg = dot_product( structure_sc%gvectors(:,j), &
                &                delta_gvec)                 &
                & * 2*pi/structure_sc%sc_size
              phase = cmplx(cos(arg),sin(arg),dp)
              
              if ( abs(new_element-old_element*phase) > 1.0e-8_dp .and. &
                 & abs(new_element-old_element/phase) > 1.0e-8_dp ) then
                call print_line('')
                call print_line('Dynamical matrices disagree.')
                call print_line('Supercell '//i)
                call print_line('Dynamical matrix '//j)
                call print_line('Atoms '//atom_1//' '//atom_2)
                call print_line('Directions '//k//' '//l)
                call print_line('Old values:')
                call print_line(real(old_element)//' '//imag(old_element))
                call print_line('New values:')
                call print_line(real(new_element)//' '//imag(new_element))
                call err()
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
    
    ! Check that displacement patterns are the same.
    call print_line('Checking displacement patterns.')
    
    deallocate(atom)
    
    ! Move dynamical matrices into dyn_mats_ibz.
    do j=1,no_kpoints
      if (sc_ids(j)/=i) then
        cycle
      endif
      
      dyn_mats_ibz(:,:,j) = new_dyn_mats(:,:,gvector_ids(j))
    enddo
    
    deallocate(new_dyn_mats)
  enddo
  
  ! Write path for fourier interpolation
  no_kspace_lines = 4
  allocate(disp_kpoints(3,no_kspace_lines+1))
  disp_kpoints(:,1) = (/ 0.0_dp, 0.0_dp, 0.0_dp /) ! GM
  disp_kpoints(:,2) = (/ 0.5_dp, 0.5_dp, 0.5_dp /) ! T
  disp_kpoints(:,3) = (/ 0.0_dp, 0.5_dp, 0.5_dp /) ! FB
  disp_kpoints(:,4) = (/ 0.0_dp, 0.0_dp, 0.0_dp /) ! GM
  disp_kpoints(:,5) = (/ 0.0_dp, 0.5_dp, 0.0_dp /) ! L
  
  ! Read in primitive symmetry group.
  prim_symmetry_group = read_group_file(str('Supercell_1/symmetry_group.dat'))
  call mkdir('new_lte')
  
  call print_line('')
  call print_line('Running fourier interpolation (this may take some time).')
  call fourier_interpolation(                      &
     & dyn_mats_ibz,                               &
     & structure,                                  &
     & grid,                                       &
     & temperature,                                &
     & kpoints,                                    &
     & disp_kpoints,                               &
     & symmetry_group,                             &
     & str('new_lte/phonon_dispersion_curve.dat'), &
     & str('new_lte/high_symmetry_points.dat'),    &
     & str('new_lte/free_energy.dat'),             &
     & str('new_lte/freq_dos.dat'))
end subroutine
end module
