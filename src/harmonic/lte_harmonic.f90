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
  character(3), parameter :: s = "(a)"        ! String format parameter
  
  ! User-input temperature
  real(dp) :: temperature
  
  ! Setup data
  integer             :: no_sc
  type(String)        :: dft_code
  type(String)        :: seedname
  type(StructureData) :: structure
  
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
  integer               :: no_gvectors
  integer,  allocatable :: gvector_ids(:)
  real(dp), allocatable :: gvec_frac(:,:)
  
  ! lte input data
  integer               :: no_kspace_lines
  real(dp), allocatable :: disp_kpoints(:,:)
  
  ! Temporary variables
  integer        :: i,j,k
  character(100) :: line
  type(String)   :: sdir
  type(String)   :: ddir
  character(100) :: dump
  
  ! File units
  integer :: no_sc_file
  integer :: dft_code_file
  integer :: seedname_file
  integer :: f
  integer :: ibz_file
  integer :: gvectors_frac_file
  integer :: list_file
  
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
  
  do i=1,no_sc
    sdir = str('Supercell_')//i
    
    ! --------------------------------------------------
    ! Read forces
    ! --------------------------------------------------
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
    
    no_kspace_lines = 4
    allocate(disp_kpoints(3,0:no_kspace_lines))
    disp_kpoints(:,0) = (/ 0.0_dp, 0.0_dp, 0.0_dp /) ! GM
    disp_kpoints(:,1) = (/ 0.5_dp, 0.5_dp, 0.5_dp /) ! T
    disp_kpoints(:,2) = (/ 0.0_dp, 0.5_dp, 0.5_dp /) ! FB
    disp_kpoints(:,3) = (/ 0.0_dp, 0.0_dp, 0.0_dp /) ! GM
    disp_kpoints(:,4) = (/ 0.0_dp, 0.5_dp, 0.0_dp /) ! L
    
    call lte( 4,                                 &
            & 1e-5_dp,                           &
            & 1e-5_dp,                           &
            & 1e-2_dp,                           &
            & structure,                         &
            & structure_scs(i),                  &
            & atoms,                             &
            & displacements,                     &
            & forces,                            &
            & 0.0_dp,                            &
            & no_kspace_lines,                   &
            & disp_kpoints,                      &
            & sdir//'/lte/freq_dos.dat',         &
            & sdir//'/lte/tdependence1.dat',     &
            & sdir//'/lte/tdependence2.dat',     &
            & sdir//'/lte/dispersion_curve.dat', &
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
    
    ! read gvectors_frac.dat
    gvectors_frac_file = open_read_file(sdir//'/lte/gvectors_frac.dat')
    read(gvectors_frac_file,*) no_gvectors
    allocate(gvector_ids(no_gvectors))
    allocate(gvec_frac(3,no_gvectors))
    do j=1,no_gvectors
      read(gvectors_frac_file,*) gvector_ids(j), gvec_frac(:,j)
    enddo
    close(gvectors_frac_file)
    
    list_file = open_append_file('list.dat')
    do j=1,no_gvectors
      do k=1,no_kpoints
        if (sc_ids(k)/=i) cycle ! skip kpoints not in this unit cell
        if (all(abs(gvec_frac(:,j)-kpoints(:,k))<tol)) then
          write(list_file,*) k, gvector_ids(j), i
          call system('cp ' &
             & //sdir//'/lte/dyn_mat.'//gvector_ids(i)//'.dat ' &
             & //'lte/dyn_mat.'//k//'.dat')
        endif
      enddo
    enddo
    close(list_file)
    
    deallocate(gvector_ids)
    deallocate(gvec_frac)
    
  enddo
  
  f = open_write_file('lte/path.dat')
  write(f,s) "0.000000 0.000000 0.000000  # GM"
  write(f,s) "0.500000 0.500000 0.500000  # T"
  write(f,s) "0.000000 0.500000 0.500000  # FB"
  write(f,s) "0.000000 0.000000 0.000000  # GM"
  write(f,s) "0.000000 0.500000 0.000000  # L"
  close(f)
  
  call fourier_interpolation(             &
     & str('structure.dat'),                   &
     & str('lte/phonon_dispersion_curve.dat'), &
     & str('lte/high_symmetry_points.dat'),    &
     & temperature,                       &
     & str('lte/free_energy.dat'),             &
     & str('lte/freq_dos.dat'),                &
     & str('grid.dat'),                        &
     & str('ibz.dat'),                         &
     & str('lte/atoms_in_primitive_cell.'),    &
     & str('lte/dyn_mat.'),                    &
     & str('lte/path.dat'))
  
end subroutine
end module
