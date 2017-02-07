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
  
  ! Parameters
  real(dp),     parameter :: tol = 1.0e-10_dp
  
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
  real(dp), allocatable :: kpoints(:,:)
  integer,  allocatable :: sc_ids(:)
  
  ! gvector data
  integer              :: no_gvectors
  integer              :: gvector_id
  real(dp)             :: gvec_frac(3)
  integer, allocatable :: gvector_ids(:)
  
  ! lte input data
  integer               :: no_kspace_lines
  real(dp), allocatable :: disp_kpoints(:,:)
  
  ! Temporary variables
  integer        :: i,j
  character(100) :: line
  type(String)   :: sdir
  type(String)   :: ddir
  character(100) :: dump
  
  ! File units
  integer :: no_sc_file
  integer :: dft_code_file
  integer :: seedname_file
  integer :: ibz_file
  integer :: gvectors_frac_file
  integer :: list_file
  integer :: grid_file
  
  ! ----------------------------------------------------------------------
  ! Get temperature from user
  ! ----------------------------------------------------------------------
  write(*,"(a)") "What temperature (K)?"
  read(*,*) temperature
  
  ! ----------------------------------------------------------------------
  ! Read in initial data
  ! ----------------------------------------------------------------------
  no_sc_file = open_read_file('no_sc.dat')
  read(no_sc_file,*) no_sc
  close(no_sc_file)
  
  dft_code_file = open_read_file('code.txt')
  read(dft_code_file,"(a)") line
  close(dft_code_file)
  dft_code = line
  
  seedname_file = open_read_file('seedname.txt')
  read(seedname_file,"(a)") line
  close(seedname_file)
  seedname = line
  
  structure = read_structure_file('structure.dat')
  
  ! Read grid file
  grid_file = open_read_file('grid.dat')
  read(grid_file,*) grid
  close(grid_file)
  
  ! Read kpoints from ibz.dat
  no_kpoints = count_lines('ibz.dat')
  allocate(kpoints(3,no_kpoints))
  allocate(sc_ids(no_kpoints))
  ibz_file = open_read_file('ibz.dat')
  do i=1,no_kpoints
    read(ibz_file,*) kpoints(:,i), dump, sc_ids(i)
  enddo
  
  ! ----------------------------------------------------------------------
  ! Read in supercell structures
  ! ----------------------------------------------------------------------
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
    
    ! Read forces from DFT output files
    no_force_constants = count_lines(sdir//'force_constants.dat')
    allocate(atoms(no_force_constants))
    allocate(displacements(no_force_constants))
    
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
              & sdir//'/lte/gvectors_frac.dat',    &
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
  
  ! Write list.dat, a mapping between kpoints, gvectors and supercells
  list_file = open_write_file('list.dat')
  allocate(gvector_ids(no_kpoints))
  do i=1,no_kpoints
    sdir = str('Supercell_')//sc_ids(i)
    gvectors_frac_file = open_read_file(sdir//'/lte/gvectors_frac.dat')
    read(gvectors_frac_file,*) no_gvectors
    do j=1,no_gvectors
      read(gvectors_frac_file,*) gvector_id, gvec_frac
      if (all(abs(gvec_frac-kpoints(:,i))<tol)) then
        gvector_ids(i) = gvector_id
        write(list_file,*) i, gvector_id, sc_ids(i)
        call system('cp ' &
           & //sdir//'/lte/dyn_mat.'//gvector_id//'.dat ' &
           & //'lte/dyn_mat.'//i//'.dat')
      endif
    enddo
    close(gvectors_frac_file)
  enddo
  close(list_file)
  
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
