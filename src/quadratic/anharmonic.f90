! ------------------------------------------------------------
! program to calculate anharmonic 1-dimensional correction
! ------------------------------------------------------------

module anharmonic_module
contains

subroutine anharmonic()
  use constants, only : dp
  use utils,     only : i2s
  use file_module
  
  use mapping_module
  use structure_module
  use string_module
  use dft_output_file_module
  use displacement_patterns_module
  
  use generate_amplitudes_module,  only : generate_amplitudes
  use calculate_anharmonic_module, only : calculate_anharmonic
  use quadratic_spline_module,     only : quadratic_spline
  use vscf_1d_module,              only : VscfReturn, vscf_1d, drop
  
  implicit none
  
  ! ----------------------------------------
  ! Parameters
  ! ----------------------------------------
  integer, parameter :: integration_points = 5000
  integer, parameter :: Nbasis = 20
  
  ! ----------------------------------------
  ! Working variables
  ! ----------------------------------------
  integer               :: no_supercells   ! no. of supercells
  integer, allocatable  :: no_atoms_sc(:)  ! no. atoms in supercell
  integer, allocatable  :: no_cells(:)     ! no_atoms_sc/no_atoms
  type(String)          :: castep          ! seedname.castep
  type(MappingData)     :: mapping         ! mapping.dat
  
  ! kpoint data
  integer, allocatable  :: kpoints(:)      ! kpoint ids
  integer, allocatable  :: gvectors(:)     ! gvector ids
  integer, allocatable  :: sc_ids(:)       ! supercell ids
  integer               :: no_kpoints
  integer, allocatable  :: sizes(:)
  
  ! ibz data
  integer, allocatable  :: multiplicity(:) ! fourth column of ibz.dat
  
  real(dp), allocatable :: energies(:,:,:)
  real(dp), allocatable :: static_energies(:)
  real(dp), allocatable :: frequencies(:,:) ! harmonic frequencies
  real(dp), allocatable :: amplitudes(:,:)
  real(dp), allocatable :: spline(:,:)
  type(VscfReturn)      :: vscf
  logical, allocatable  :: sc_acoustic(:)  ! if Supercell_i/acoustic.dat exists
  
  type(String)          :: sdir          ! Supercell_*          directory name
  type(String)          :: kpoint_dir    ! Supercell_*/kpoint.* directory name
  
  real(dp), allocatable :: eigenvals(:,:,:)
  real(dp), allocatable :: harmonic(:,:,:)
  
  
  type(StructureData)   :: structure     ! the contents of structure.dat
  type(StructureData)   :: structure_sc
  type(DispPatterns)    :: disp_patterns
  
  type(String)          :: harmonic_path ! the path to the harmonic directory
  type(String)          :: filename
  type(String)          :: ibz_filename
  
  type(DftOutputFile)   :: dft_output_file
  
  integer               :: i,j,k
  character(100)        :: temp_char
  character(100)        :: dump
  
  ! ----------------------------------------
  ! File units
  ! ----------------------------------------
  integer :: harmonic_path_file
  integer :: no_sc_file
  integer :: seedname_file
  integer :: list_file
  integer :: ibz_file
  integer :: result_file
  
  ! ----------------------------------------
  ! Read in data
  ! ----------------------------------------
  
  ! read in harmonic_path
  harmonic_path_file = open_read_file('harmonic_path.dat')
  read(harmonic_path_file,"(a)") temp_char
  close(harmonic_path_file)
  harmonic_path = trim(temp_char)
  
  ! read the number of Supercell_* directories into no_supercells
  no_sc_file = open_read_file(harmonic_path//'/no_sc.dat')
  read(no_sc_file,*) no_supercells
  close(no_sc_file)
  
  ! allocate arrays of size no_supercells
  allocate(sc_acoustic(no_supercells))
  allocate(no_atoms_sc(no_supercells))
  allocate(no_cells(no_supercells))
  allocate(static_energies(no_supercells))
  
  ! read structure data
  structure = read_structure_file(harmonic_path//'/structure.dat')
  
  ! read the castep seedname into castep variable
  seedname_file = open_read_file('seedname.txt')
  read(seedname_file,"(a)") temp_char
  close(seedname_file)
  castep = trim(temp_char)//'.castep'
  
  ! read sampling data from mapping.dat
  mapping = read_mapping_file('mapping.dat')
  
  ! check for Supercell_*/acoustic.dat
  do i=1,no_supercells
    sc_acoustic(i) = file_exists(str('Supercell_')//i//'/acoustic.dat')
  enddo
  
  ! Read kpoints
  no_kpoints = count_lines(harmonic_path//'/list.dat')
  allocate(kpoints(no_kpoints))
  allocate(gvectors(no_kpoints))
  allocate(sc_ids(no_kpoints))
  list_file = open_read_file(harmonic_path//'/list.dat')
  do i=1,no_kpoints
    read(list_file,*) kpoints(i),gvectors(i),sc_ids(i)
  enddo
  close(list_file)
  
  allocate(sizes(no_kpoints))
  
  ! read no_atoms_sc and no_cells
  do i=1,no_supercells
    if (.not. sc_acoustic(i)) then
      sdir = str('Supercell_')//i
      filename = harmonic_path//'/'//sdir//'/structure.dat' 
      structure_sc = read_structure_file(filename)
      no_atoms_sc(i) = structure_sc%no_atoms
      no_cells(i) = no_atoms_sc(i)/structure%no_atoms
    endif
  enddo
  
  ! allocate arguments for calculate_anharmonic
  allocate(energies(no_kpoints,structure%no_modes,mapping%count))
  allocate(frequencies(no_kpoints,structure%no_modes))
  allocate(harmonic(no_kpoints,structure%no_modes,Nbasis))
  allocate(eigenvals(no_kpoints,structure%no_modes,Nbasis))
  
  ! read multiplicity from ibz.dat
  ibz_filename = harmonic_path//'/ibz.dat'
  allocate(multiplicity(count_lines(ibz_filename)))
  ibz_file = open_read_file(ibz_filename)
  do i=1,size(multiplicity)
    read(ibz_file,*) dump, dump, dump, multiplicity(i)
  enddo
  close(ibz_file)
  
  ! read data from supercells
  do i=1,no_supercells
    if (.not. sc_acoustic(i)) then
      filename = str('Supercell_')//i//'/static/'//castep
      dft_output_file = read_castep_output_file(filename)
      static_energies(i) = dft_output_file%energy
    endif
  enddo
  
  ! Loop over kpoints
  do i=1,no_kpoints
    if (.not. sc_acoustic(sc_ids(i))) then
      kpoint_dir = str('kpoint.')//kpoints(i)
      
      ! set sizes
      sizes(i) = no_cells(sc_ids(i))
      
      ! Read frequencies
      filename = harmonic_path//'/Supercell_'//sc_ids(i)//'/lte/disp_patterns.dat'
      disp_patterns = read_disp_patterns_file(filename,structure_sc%no_modes)
      frequencies(i,:) = disp_patterns%frequencies(:,gvectors(i))
      
      do j=1,structure%no_modes
        ! read energies
        do k=mapping%first,mapping%last
          filename = kpoint_dir// &
                   & '/configurations/mode.'//j//'.'//k//'/'//castep
          if (file_exists(filename)) then
            dft_output_file = read_castep_output_file(filename)
            energies(i,j,k) = dft_output_file%energy
          else
            energies(i,j,k) = static_energies(sc_ids(i))
          endif
        enddo
      enddo
    endif
  enddo
  
  ! ----------------------------------------
  ! Process data.
  ! ----------------------------------------
  
  ! Calculate anharmonic 1-dimensional correction
  do i=1,no_kpoints
    if (kpoints(i)/=1 .or. .not. any(sc_acoustic)) then
      do j=1,structure%no_modes
        
        ! generate amplitudes
        ! generate potential at {q} defined by map
        amplitudes = generate_amplitudes(mapping, energies(i,j,:),&
          &frequencies(i,j),sizes(i))
        
        ! fit splines
        ! interpolate potential onto integration_points points
        ! min(q) and max(q) are unchanged
        spline = quadratic_spline(integration_points, amplitudes)
        
        ! calculate 1-d anharmonic energy
        vscf = vscf_1d(frequencies(i,j), spline, Nbasis)
        
        harmonic(i,j,:) = vscf%harmonic
        eigenvals(i,j,:) = vscf%eigenvals
        
        deallocate(amplitudes)
        deallocate(spline)
        call drop(vscf)
        
      enddo
    endif
  enddo
  
  ! ----------------------------------------
  ! Write output.
  ! ----------------------------------------
  
  ! copy acoustic.dat files
  call system('mkdir anharmonic')
  do i=1,no_supercells
    if (sc_acoustic(i)) then
      call system('cp Supercell_'//trim(i2s(i))//'/acoustic.dat anharmonic')
    endif
  enddo
  
  ! calculate free energy, F(T), for harmonic and anharmonic cases
  ! write output to anharmonic_correction.dat
  result_file = open_write_file('anharmonic/anharmonic_correction.dat')
  call calculate_anharmonic(multiplicity,structure%no_modes,Nbasis,harmonic,&
    &eigenvals,result_file)
  close(result_file)
end subroutine

end module
