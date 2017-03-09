! ------------------------------------------------------------
! program to calculate anharmonic 1-dimensional correction
! ------------------------------------------------------------

module anharmonic_module
contains

subroutine anharmonic()
  use constants, only : dp, eV
  use file_module
  
  use mapping_module
  use structure_module
  use string_module
  use dft_output_file_module
  use displacement_patterns_module
  use supercell_module
  
  use calculate_anharmonic_module
  use quadratic_spline_module
  use vscf_1d_module
  
  implicit none
  
  ! ----------------------------------------
  ! Parameters
  ! ----------------------------------------
  integer, parameter :: integration_points = 5000
  integer, parameter :: Nbasis = 20
  
  ! ----------------------------------------
  ! Input variables
  ! ----------------------------------------
  type(String) :: dft_code      ! The dft code name (castep,vasp,qe)
  type(String) :: seedname      ! The dft input file seedname
  type(String) :: harmonic_path ! The path to the harmonic directory
  
  ! ----------------------------------------
  ! Working variables
  ! ----------------------------------------
  integer               :: no_supercells   ! no. of supercells
  type(MappingData)     :: mapping         ! mapping.dat
  
  type(StructureData)              :: structure
  type(SupercellData), allocatable :: supercells(:)
  type(StructureData), allocatable :: structure_scs(:)
  type(DispPatterns)    :: disp_patterns
  
  ! kpoint data
  integer               :: no_kpoints
  integer, allocatable  :: multiplicity(:)
  integer, allocatable  :: sc_ids(:)       ! supercell ids
  integer, allocatable  :: gvectors(:)     ! gvector ids
  
  integer, allocatable  :: sizes(:)
  
  real(dp), allocatable :: energies(:,:,:)
  real(dp), allocatable :: static_energies(:)
  real(dp), allocatable :: frequencies(:,:) ! harmonic frequencies
  real(dp)              :: amplitude
  real(dp), allocatable :: amplitudes(:,:)
  real(dp), allocatable :: spline(:,:)
  type(VscfData)        :: vscf
  logical, allocatable  :: sc_acoustic(:)  ! if Supercell_i/acoustic.dat exists
  
  ! Directory names
  type(String) :: sdir       ! Supercell_*
  type(String) :: kpoint_dir ! kpoint.*
  type(String) :: ddir       ! kpoint.*/mode.*.*
  
  real(dp), allocatable :: eigenvals(:,:,:)
  real(dp), allocatable :: harmonic(:,:,:)
  
  type(String)          :: filename
  
  type(DftOutputFile)   :: dft_output_file
  
  ! ----------------------------------------
  ! Temporary variables
  ! ----------------------------------------
  integer               :: i,j,k
  character(100)        :: dump
  
  ! ----------------------------------------
  ! File contents
  ! ----------------------------------------
  type(String), allocatable :: user_inputs(:)
  
  ! ----------------------------------------
  ! File units
  ! ----------------------------------------
  integer :: no_sc_file
  integer :: list_file
  integer :: result_file
  
  ! ----------------------------------------
  ! Read in data
  ! ----------------------------------------
  ! Read in previous user inputs
  user_inputs = read_lines('user_input.txt')
  dft_code = user_inputs(1)
  seedname = user_inputs(2)
  harmonic_path = user_inputs(3)
  
  ! read the number of Supercell_* directories into no_supercells
  no_sc_file = open_read_file(harmonic_path//'/no_sc.dat')
  read(no_sc_file,*) no_supercells
  close(no_sc_file)
  
  ! allocate arrays of size no_supercells
  allocate(sc_acoustic(no_supercells))
  allocate(static_energies(no_supercells))
  
  ! read structure data
  structure = read_structure_file( harmonic_path//'/structure.dat', &
                                 & identity_supercell())
  
  ! read sampling data from mapping.dat
  mapping = read_mapping_file('mapping.dat')
  
  ! check for Supercell_*/acoustic.dat
  do i=1,no_supercells
    sc_acoustic(i) = file_exists('Supercell_'//i//'/acoustic.dat')
  enddo
  
  ! Read kpoints
  no_kpoints = count_lines(harmonic_path//'/list.dat')
  allocate(multiplicity(no_kpoints))
  allocate(gvectors(no_kpoints))
  allocate(sc_ids(no_kpoints))
  list_file = open_read_file(harmonic_path//'/list.dat')
  do i=1,no_kpoints
    read(list_file,*) dump,dump,dump,multiplicity(i),sc_ids(i),gvectors(i)
  enddo
  close(list_file)
  
  ! Read supercell structures
  supercells = read_supercells_file(str('supercells.dat'))
  allocate(structure_scs(no_supercells))
  do i=1,no_supercells
    sdir = 'Supercell_'//i
    filename = harmonic_path//'/'//sdir//'/structure.dat' 
    structure_scs(i) = read_structure_file(filename,supercells(i))
  enddo
  
  ! read data from supercells
  do i=1,no_supercells
    sdir = 'Supercell_'//i
    if (.not. sc_acoustic(i)) then
      dft_output_file = read_dft_output_file(dft_code,sdir//'/static',seedname)
      static_energies(i) = dft_output_file%energy
    endif
  enddo
  
  ! Loop over kpoints
  allocate(sizes(no_kpoints))
  allocate(frequencies(structure%no_modes,no_kpoints))
  allocate(energies(mapping%count,structure%no_modes,no_kpoints))
  do i=1,no_kpoints
    if (.not. sc_acoustic(sc_ids(i))) then
      kpoint_dir = 'kpoint.'//i
      
      ! set sizes
      sizes(i) = structure_scs(sc_ids(i))%no_atoms/structure%no_atoms
      
      ! Read frequencies
      filename = harmonic_path//'/Supercell_'//sc_ids(i)// &
         & '/lte/disp_patterns.dat'
      disp_patterns = read_disp_patterns_file( filename, &
         & structure_scs(sc_ids(i))%no_modes)
      frequencies(:,i) = disp_patterns%frequencies(:,gvectors(i))
      
      do j=1,structure%no_modes
        ! read energies
        do k=mapping%first,mapping%last
          ddir = kpoint_dir//'/mode.'//j//'.'//k
          
          if (dft_code=='castep') then
            filename = ddir//'/'//seedname//'.castep'
          elseif (dft_code=='qe') then
            filename = ddir//'/'//seedname//'.out'
          endif
          
          if (file_exists(filename)) then
            dft_output_file = read_dft_output_file(dft_code,ddir,seedname)
            energies(k,j,i) = dft_output_file%energy
          else
            energies(k,j,i) = static_energies(sc_ids(i))
          endif
        enddo
      enddo
    endif
  enddo
  
  ! ----------------------------------------
  ! Process data.
  ! ----------------------------------------
  
  ! Calculate anharmonic 1-dimensional correction
  allocate(harmonic(Nbasis,structure%no_modes,no_kpoints))
  allocate(eigenvals(Nbasis,structure%no_modes,no_kpoints))
  do i=1,no_kpoints
    if (i/=1 .or. .not. any(sc_acoustic)) then
      do j=1,structure%no_modes
        
        ! generate amplitudes, {(x,V(x))}
        ! generate potential at {q} defined by map
        allocate(amplitudes(2,mapping%count))
        amplitude = -mapping%max*eV/(2*dabs(frequencies(j,i)))
        do k=1,mapping%count
          amplitudes(1,k) = amplitude+(k-1)*dabs(amplitude/mapping%first)
          amplitudes(2,k) = (energies(k,j,i)-energies(mapping%mid,j,i)) &
                        & / sizes(i)
        enddo
        
        ! fit splines
        ! interpolate potential onto integration_points points
        ! min(q) and max(q) are unchanged
        spline = quadratic_spline(integration_points, amplitudes)
        
        ! calculate 1-d anharmonic energy
        vscf = vscf_1d(frequencies(j,i), spline, Nbasis)
        
        harmonic(:,j,i) = vscf%harmonic
        eigenvals(:,j,i) = vscf%eigenvals
        
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
      call system('cp Supercell_'//i//'/acoustic.dat anharmonic')
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
