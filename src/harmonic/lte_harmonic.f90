module lte_harmonic_module

contains

! Program to construct and execute LTE
subroutine lte_harmonic()
  use constants, only : dp, eV_per_A_to_au
  use file_module
  use string_module
  use structure_module
  use dft_output_file_module
  use lte_module
  use fourier_interpolation_module
  implicit none
  
  ! User-input temperature
  real(dp) :: temperature
  
  ! Setup data
  integer             :: no_sc
  type(String)        :: dft_code
  type(String)        :: seedname
  type(StructureData) :: structure
  integer             :: grid(3)
  
  ! Supercell-specific setup data
  type(StructureData), allocatable :: structure_scs(:)
  
  ! Force constant data
  integer               :: no_force_constants
  integer,  allocatable :: atoms(:)
  integer,  allocatable :: displacements(:)
  real(dp), allocatable :: forces(:,:,:)
  type(DftOutputFile)   :: positive
  type(DftOutputFile)   :: negative
  
  ! kpoint data
  integer               :: no_kpoints
  integer,  allocatable :: kpoints(:,:)
  integer,  allocatable :: multiplicity(:)
  integer,  allocatable :: sc_ids(:)
  
  ! gvector data
  integer              :: no_gvectors
  integer              :: gvector_id
  integer              :: gvec(3)
  integer, allocatable :: gvector_ids(:)
  
  ! lte input data
  integer               :: no_kspace_lines
  real(dp), allocatable :: disp_kpoints(:,:)
  
  ! Temporary variables
  integer        :: i,j
  character(100) :: line
  type(String)   :: sdir
  type(String)   :: ddir
  
  ! File units
  integer :: user_input_file
  integer :: no_sc_file
  integer :: ibz_file
  integer :: gvectors_file
  integer :: list_file
  integer :: grid_file
  integer :: force_constants_file
  
  ! ----------------------------------------------------------------------
  ! Get temperature from user
  ! ----------------------------------------------------------------------
  write(*,"(a)") "What temperature (K)?"
  read(*,*) temperature
  
  ! ----------------------------------------------------------------------
  ! Read in initial data
  ! ----------------------------------------------------------------------
  user_input_file = open_read_file('user_input.txt')
  read(user_input_file,"(a)") line
  dft_code = trim(line)
  read(user_input_file,"(a)") line
  seedname = trim(line)
  close(user_input_file)
  
  no_sc_file = open_read_file('no_sc.dat')
  read(no_sc_file,*) no_sc
  close(no_sc_file)
  
  structure = read_structure_file('structure.dat')
  
  ! Read grid file
  grid_file = open_read_file('grid.dat')
  read(grid_file,*) grid
  close(grid_file)
  
  ! Read kpoints from ibz.dat
  no_kpoints = count_lines('ibz.dat')
  allocate(kpoints(3,no_kpoints))
  allocate(multiplicity(no_kpoints))
  allocate(sc_ids(no_kpoints))
  ibz_file = open_read_file('ibz.dat')
  do i=1,no_kpoints
    read(ibz_file,*) kpoints(:,i), multiplicity(i), sc_ids(i)
  enddo
  
  ! Read in supercell structures
  allocate(structure_scs(no_sc))
  do i=1,no_sc
    sdir = str('Supercell_')//i
    structure_scs(i) = read_structure_file(sdir//'/structure.dat')
  enddo
  
  ! ----------------------------------------------------------------------
  ! Make directories
  ! ----------------------------------------------------------------------
  call system('mkdir lte')
  do i=1,no_sc
    call system('mkdir '//sdir//'/lte')
  enddo
  
  ! ----------------------------------------------------------------------
  ! Loop over supercells
  ! ----------------------------------------------------------------------
  do i=1,no_sc
    sdir = str('Supercell_')//i
    
    ! Read in force constants
    no_force_constants = count_lines(sdir//'/force_constants.dat')
    allocate(atoms(no_force_constants))
    allocate(displacements(no_force_constants))
    force_constants_file = open_read_file(sdir//'/force_constants.dat')
    do j=1,no_force_constants
      read(force_constants_file,*) atoms(j), displacements(j)
    enddo
    close(force_constants_file)
    
    ! Read forces from DFT output files
    allocate(forces(3,structure_scs(i)%no_atoms,no_force_constants))
    do j=1,no_force_constants
      ddir = sdir//'/atom.'//atoms(j)//'/disp.'//displacements(j)
      
      positive = read_dft_output_file(dft_code,ddir//'/positive',seedname)
      negative = read_dft_output_file(dft_code,ddir//'/negative',seedname)
      
      forces(:,:,j) = (positive%forces-negative%forces)/(0.02d0*eV_per_A_to_au)
    enddo
    
    call lte_4( 1e-5_dp,                           &
              & 1e-5_dp,                           &
              & 1e-2_dp,                           &
              & structure,                         &
              & structure_scs(i),                  &
              & atoms,                             &
              & displacements,                     &
              & forces,                            &
              & 0.0_dp,                            &
              & sdir//'/lte/kpairs.dat',           &
              & sdir//'/lte/freq_grids.dat',       &
              & sdir//'/lte/disp_patterns.dat',    &
              & sdir//'/lte/kdisp_patterns.dat',   &
              & sdir//'/lte/pol_vec.dat',          &
              & sdir//'/lte/gvectors.dat',         &
              & sdir//'/lte/error.txt',            &
              & sdir//'/lte/dyn_mat.',             &
              & sdir//'/lte/atoms_in_primitive_cell.dat')
    
    if (file_exists(sdir//'/lte/error.txt')) then
      write(*,"(a)") "There is an error in lte: check error.txt file."
      stop
    endif
    
    deallocate(atoms)
    deallocate(displacements)
    deallocate(forces)
  enddo
  
  ! Locate the corresponding gvector for each kpoint
  allocate(gvector_ids(no_kpoints))
  do i=1,no_kpoints
    sdir = str('Supercell_')//sc_ids(i)
    gvectors_file = open_read_file(sdir//'/lte/gvectors.dat')
    read(gvectors_file,*) no_gvectors
    do j=1,no_gvectors
      read(gvectors_file,*) gvector_id, gvec
      if (all(nint(matmul(gvec,structure%lattice))==kpoints(:,i))) then
        gvector_ids(i) = gvector_id
      endif
    enddo
    close(gvectors_file)
  enddo
  
  ! Write list.dat, a mapping between kpoints, supercells and gvectors
  list_file = open_write_file('list.dat')
    write(list_file,*) kpoints(:,i),    &
                     & multiplicity(i), &
                     & sc_ids(i),       &
                     & gvector_ids(i)
  close(list_file)
  
  ! Write path for fourier interpolation
  no_kspace_lines = 4
  allocate(disp_kpoints(3,0:no_kspace_lines))
  disp_kpoints(:,0) = (/ 0.0_dp, 0.0_dp, 0.0_dp /) ! GM
  disp_kpoints(:,1) = (/ 0.5_dp, 0.5_dp, 0.5_dp /) ! T
  disp_kpoints(:,2) = (/ 0.0_dp, 0.5_dp, 0.5_dp /) ! FB
  disp_kpoints(:,3) = (/ 0.0_dp, 0.0_dp, 0.0_dp /) ! GM
  disp_kpoints(:,4) = (/ 0.0_dp, 0.5_dp, 0.0_dp /) ! L
    
  call fourier_interpolation(                  &
     & structure,                              &
     & grid,                                   &
     & temperature,                            &
     & kpoints, sc_ids, gvector_ids,           &
     & str('lte/atoms_in_primitive_cell.'),    &
     & str('lte/dyn_mat.'),                    &! Supercell_*/lte/dyn_mat.*.dat
     & disp_kpoints,                           &
     & str('lte/phonon_dispersion_curve.dat'), &
     & str('lte/high_symmetry_points.dat'),    &
     & str('lte/free_energy.dat'),             &
     & str('lte/freq_dos.dat'))
  
end subroutine
end module
