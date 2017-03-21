module test_lte_module
contains
subroutine test_lte()
  use utils, only : mkdir
  use string_module
  use file_module
  use structure_module
  use group_module
  use atom_mapping_module
  use lte_module
  use fourier_interpolation_module
  implicit none
  
  ! Tolerance.
  real(dp), parameter :: tol = 1.0e-8_dp
  
  ! Directories and files.
  type(String)              :: sdir
  type(String), allocatable :: no_sc_file(:)
  type(String)              :: copy_dir
  type(String)              :: lte_dir
  type(String), allocatable :: force_constants_file(:)
  type(String), allocatable :: grid_file(:)
  type(String), allocatable :: ibz_file(:)
  
  ! Supercell data.
  integer :: no_sc
  
  ! Structure data.
  type(StructureData) :: structure
  type(StructureData) :: structure_sc
  type(StructureData) :: structure_sc_old
  
  ! lte file data.
  type(String), allocatable :: old_lte_file(:)
  integer                   :: new_lte_file
  integer                   :: atoms_start_line
  integer                   :: atoms_end_line
  integer                   :: force_start_line
  integer                   :: force_end_line
  integer                   :: line_no
  
  ! Atom mapping data.
  type(Group) :: new_to_old
  type(Group) :: old_to_new
  
  ! Force constant data.
  real(dp), allocatable :: force_constants(:,:,:,:)
  real(dp), allocatable :: force_constants_2(:,:,:)
  
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
  type(Group), allocatable :: symmetry_group(:)
  
  ! Dynamic matrix variables.
  type(String), allocatable :: new_dyn_mat_file(:)
  type(String), allocatable :: old_dyn_mat_file(:)
  type(String), allocatable :: new_line(:)
  type(String), allocatable :: old_line(:)
  complex(dp)               :: new_element
  complex(dp)               :: old_element
  
  ! Temporary variables.
  integer                   :: i,j,k,l
  integer                   :: atom_1,atom_2,atom_1p,atom_2p
  integer                   :: mode_1,mode_2
  integer                   :: start_line
  type(String), allocatable :: line(:)
  
  ! ----------------------------------------------------------------------
  ! Read in directory to compare against and copy from.
  ! ----------------------------------------------------------------------
  call print_line('')
  call print_line('This test will check lte against a previous &
     &calculations.')
  call print_line('')
  call print_line('Where is the harmonic directory for comparison?')
  copy_dir = read_line_from_user()
  
  ! ----------------------------------------------------------------------
  ! Get temperature from user
  ! ----------------------------------------------------------------------
  call print_line('What temperature (K)?')
  temperature = dble(read_line_from_user())
  
  ! ----------------------------------------------------------------------
  ! Read in initial data
  ! ----------------------------------------------------------------------
  ! Read in no_supercells.
  no_sc_file = read_lines('no_sc.dat')
  no_sc = int(no_sc_file(1))
  
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
  
  do i=1,no_sc
    call print_line('Checking supercell '//i)
    sdir = 'Supercell_'//i
    
    ! Read in supercell structure.
    structure_sc = read_structure_file(sdir//'/structure.dat')
    structure_sc_old = read_structure_file(copy_dir//'/'//sdir//'/structure.dat')
    
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
    ! Write everything above atoms.
    do j=1,atoms_start_line
      call print_line(new_lte_file,old_lte_file(j))
    enddo
    ! write atoms
    do j=atoms_start_line+1,atoms_end_line-1
      line_no = operate(new_to_old, j-atoms_start_line) + atoms_start_line
      call print_line(new_lte_file, old_lte_file(line_no))
    enddo
    ! Write everything between atoms and force constants.
    do j=atoms_end_line,force_start_line
      call print_line(new_lte_file,old_lte_file(j))
    enddo
    ! Write force constants
    do j=force_start_line+1,force_end_line-1
      line = split(old_lte_file(j))
      line(1) = operate(old_to_new,int(line(1)))
      line(3) = operate(old_to_new,int(line(3)))
      call print_line(new_lte_file,join(line))
    enddo
    ! Write everything below force constants.
    do j=force_end_line,size(old_lte_file)
      call print_line(new_lte_file,old_lte_file(j))
    enddo
    close(new_lte_file)
    
    ! Run old lte.
    call system('cd '//lte_dir//'; lte_lower > lte.out')
    
    ! Read in force constants from old lte.
    force_constants_file = read_lines(lte_dir//'/force_constants.dat')
    
    allocate(force_constants(3,3,structure_sc%no_atoms,structure_sc%no_atoms))
    do atom_1=1,structure_sc%no_atoms
      do atom_2=1,structure_sc%no_atoms
        start_line = 5*((atom_1-1)*structure_sc%no_atoms+(atom_2-1))+1
        line = split(force_constants_file(start_line))
        if (any(int(line)/=(/atom_1,atom_2/))) then
          call print_line('Force constants file in unexpected order.')
          stop
        endif
        do j=1,3
          line = split(force_constants_file(start_line+j))
          force_constants(j,:,atom_2,atom_1) = dble(line)
        enddo
      enddo
    enddo
    
    do j=1,3
      call print_line(force_constants(j,:,4,1))
    enddo
    
    do j=1,3
      call print_line(force_constants(j,:,1,4))
    enddo
    
    ! Convert force constants into Gvector-mode-mode representation.
    allocate(force_constants_2( structure%no_modes, &
                              & structure%no_modes, &
                              & structure_sc%sc_size))
    do j=1,structure_sc%sc_size
      do atom_1=1,structure%no_atoms
        atom_1p = structure_sc%gvec_and_prim_to_atom(atom_1,1)
        mode_1 = (atom_1-1)*3+1
        do atom_2=1,structure%no_atoms
          atom_2p = structure_sc%gvec_and_prim_to_atom(atom_2,j)
          mode_2 = (atom_2-1)*3+1
          force_constants_2(mode_2:mode_2+2,mode_1:mode_1+2,j) = &
             & transpose(force_constants(:,:,atom_2p,atom_1p))
        enddo
      enddo
    enddo
    
    ! Check force constants.
    do j=1,structure_sc%sc_size
      if (any(abs( transpose(force_constants_2(:,:,j)) &
         & - force_constants_2(:,:,structure_sc%paired_gvec(j))) >tol)) then
        call print_line('')
        call print_line('Old force constants are not symmetric.')
        call print_line('G-vector: '//j)
        call print_line('Paired G-vector: '//structure_sc%paired_gvec(j))
        
        call print_line('Atom mapping:')
        call print_line(new_to_old%operation)
        
        call print_line('')
        do atom_1=1,structure%no_atoms
          do atom_2=1,structure%no_atoms
            call print_line( atom_1                             //' '// &
                           & atom_2                             //' '// &
                           & force_constants_2(atom_2,atom_1,j) //' '// &
                           & force_constants_2( atom_1,                 &
                                              & atom_2,                 &
                                              & structure_sc%paired_gvec(j)))
          enddo
        enddo
        stop
      endif
    enddo
    
    ! Run new lte.
    call mkdir(sdir//'/lte')
    call lte_4( structure,                         &
              & structure_sc,                      &
              & force_constants_2,                 &
              & 0.0_dp,                            &
              & sdir//'/lte/freq_grids.dat',       &
              & sdir//'/lte/disp_patterns.dat',    &
              & sdir//'/lte/kdisp_patterns.dat',   &
              & sdir//'/lte/pol_vec.dat',          &
              & sdir//'/lte/error.txt',            &
              & sdir//'/lte/dyn_mat.')
    
    deallocate(force_constants)
    deallocate(force_constants_2)
    
    ! Check that dynamical matrices are the same.
    do j=1,structure_sc%sc_size
      call print_line('Checking dynamical matrix '//j)
      new_dyn_mat_file = read_lines(sdir//'/lte/dyn_mat.'//j//'.dat')
      old_dyn_mat_file = read_lines(lte_dir//'/dyn_mat.'//j//'.dat')
      do atom_1=1,structure%no_atoms
        do atom_2=1,structure_sc%sc_size
          do k=1,3
            do l=1,3
              new_line = split(new_dyn_mat_file(   ((atom_1-1)*3+k-1)     &
                                               & * (structure%no_atoms*3) &
                                               & + ((atom_2-1)*3+l)))
              old_line = split(old_dyn_mat_file(   ((atom_1-1)*3+k-1)     &
                                               & * (structure%no_atoms*3) &
                                               & + ((atom_2-1)*3+l)))
              
              if (any(int(new_line(1:4))/=(/atom_1,k,atom_2,l/))) then
                call print_line('New dyn_mat file in unexpected order.')
                stop
              endif
              
              if (any(int(old_line(1:4))/=(/atom_1,k,atom_2,l/))) then
                call print_line('Old dyn_mat file in unexpected order.')
                stop
              endif
              
              new_element = cmplx(dble(new_line(5)),dble(new_line(6)),dp)
              old_element = cmplx(dble(old_line(5)),dble(old_line(6)),dp)
              
              if ( abs(new_element-old_element) > tol) then
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
                stop
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  
  ! Write path for fourier interpolation
  no_kspace_lines = 4
  allocate(disp_kpoints(3,no_kspace_lines+1))
  disp_kpoints(:,1) = (/ 0.0_dp, 0.0_dp, 0.0_dp /) ! GM
  disp_kpoints(:,2) = (/ 0.5_dp, 0.5_dp, 0.5_dp /) ! T
  disp_kpoints(:,3) = (/ 0.0_dp, 0.5_dp, 0.5_dp /) ! FB
  disp_kpoints(:,4) = (/ 0.0_dp, 0.0_dp, 0.0_dp /) ! GM
  disp_kpoints(:,5) = (/ 0.0_dp, 0.5_dp, 0.0_dp /) ! L
  
  symmetry_group = read_group_file(str('Supercell_1/symmetry_group.dat'))
  call mkdir('lte')
!  call fourier_interpolation(                  &
!     & structure,                              &
!     & grid,                                   &
!     & temperature,                            &
!     & kpoints, sc_ids, gvector_ids,           &
!     & str('lte/dyn_mat.'),                    &! Supercell_*/lte/dyn_mat.*.dat
!     & disp_kpoints,                           &
!     & symmetry_group,                         &
!     & str('lte/phonon_dispersion_curve.dat'), &
!     & str('lte/high_symmetry_points.dat'),    &
!     & str('lte/free_energy.dat'),             &
!     & str('lte/freq_dos.dat'))
end subroutine
end module
