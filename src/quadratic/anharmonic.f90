! ------------------------------------------------------------
! program to calculate anharmonic 1-dimensional correction
! ------------------------------------------------------------

module anharmonic_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

subroutine anharmonic(wd)
  use utils_module, only : mkdir, make_dft_output_filename
  use mapping_module
  use structure_module
  use dft_output_file_module
  use supercell_module
  use kpoints_module
  use calculate_anharmonic_module
  use quadratic_spline_module
  use vscf_1d_module
  implicit none
  
  ! Working directory.
  type(String), intent(in) :: wd
  
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
  type(StructureData), allocatable :: structure_scs(:)
  
  ! K-point data.
  type(StructureData)           :: structure_grid
  type(KpointData), allocatable :: kpoints(:)
  
  real(dp), allocatable :: energies(:,:,:)
  real(dp), allocatable :: static_energies(:)
  type(String), allocatable :: frequencies_file(:)
  real(dp), allocatable :: frequencies(:,:) ! harmonic frequencies
  real(dp)              :: amplitude
  real(dp), allocatable :: amplitudes(:,:)
  real(dp), allocatable :: spline(:,:)
  type(VscfData), allocatable :: vscf
  logical, allocatable  :: sc_acoustic(:)  ! if Supercell_i/acoustic.dat exists
  
  ! Directory names
  type(String) :: sdir       ! Supercell_*
  type(String) :: ddir       ! kpoint_*/mode_*/amplitude_*
  
  real(dp), allocatable :: eigenvals(:,:,:)
  real(dp), allocatable :: harmonic(:,:,:)
  
  type(String)          :: filename
  type(String)          :: dft_output_filename
  
  type(DftOutputFile)   :: dft_output_file
  
  ! ----------------------------------------
  ! Temporary variables
  ! ----------------------------------------
  integer               :: i,j,k
  
  ! ----------------------------------------
  ! Files.
  ! ----------------------------------------
  type(String), allocatable :: user_input_file(:)
  type(String), allocatable :: no_sc_file(:)
  integer                   :: result_file
  
  ! ----------------------------------------
  ! Read in data
  ! ----------------------------------------
  ! Read in previous user inputs
  user_input_file = read_lines(wd//'/user_input.txt')
  dft_code = user_input_file(1)
  seedname = user_input_file(2)
  harmonic_path = user_input_file(3)
  
  ! read the number of Supercell_* directories into no_supercells
  no_sc_file = read_lines(harmonic_path//'/no_sc.dat')
  no_supercells = int(no_sc_file(1))
  
  ! allocate arrays of size no_supercells
  allocate(sc_acoustic(no_supercells))
  allocate(static_energies(no_supercells))
  
  ! read structure data
  structure = read_structure_file(harmonic_path//'/structure.dat')
  
  ! read sampling data from mapping.dat
  mapping = read_mapping_file(wd//'/mapping.dat')
  
  ! check for Supercell_*/acoustic.dat
  do i=1,no_supercells
    sc_acoustic(i) = file_exists(wd//'/Supercell_'//i//'/acoustic.dat')
  enddo
  
  ! Read kpoints
  structure_grid = read_structure_file(harmonic_path//'/structure_grid.dat')
  kpoints = read_kpoints_file(harmonic_path//'/kpoints_ibz.dat')
  
  ! Read supercell structures
  allocate(structure_scs(no_supercells))
  do i=1,no_supercells
    sdir = wd//'/Supercell_'//i
    filename = harmonic_path//'/'//sdir//'/structure.dat' 
    structure_scs(i) = read_structure_file(filename)
  enddo
  
  ! Read data from supercells.
  do i=1,no_supercells
    if (.not. sc_acoustic(i)) then
      sdir = wd//'/Supercell_'//i
      dft_output_filename = make_dft_output_filename(dft_code,seedname)
      dft_output_filename = sdir//'/'//dft_output_filename
      dft_output_file = read_dft_output_file(dft_code,dft_output_filename)
      static_energies(i) = dft_output_file%energy
    endif
  enddo
  
  ! Loop over K-points.
  allocate(frequencies(structure%no_modes,size(kpoints)))
  allocate(energies(mapping%count,structure%no_modes,size(kpoints)))
  do i=1,size(kpoints)
    if (.not. sc_acoustic(kpoints(i)%sc_id)) then
      ! Read frequencies
      frequencies_file = read_lines( &
         & harmonic_path//'/kpoint_'//i//'/frequencies.dat')
      frequencies(:,i) = dble(frequencies_file)
      
      do j=1,structure%no_modes
        ! read energies
        do k=mapping%first,mapping%last
          ddir = wd//'/kpoint_'//i//'/mode_'//j//'/amplitude_'//k
          
          dft_output_filename = make_dft_output_filename(dft_code,seedname)
          dft_output_filename = ddir//'/'//dft_output_filename
          
          if (file_exists(dft_output_filename)) then
            dft_output_file = read_dft_output_file(dft_code,dft_output_filename)
            energies(k,j,i) = dft_output_file%energy
          else
            energies(k,j,i) = static_energies(kpoints(i)%sc_id)
          endif
        enddo
      enddo
    endif
  enddo
  
  ! ----------------------------------------
  ! Process data.
  ! ----------------------------------------
  
  ! Calculate anharmonic 1-dimensional correction
  allocate(harmonic(Nbasis,structure%no_modes,size(kpoints)))
  allocate(eigenvals(Nbasis,structure%no_modes,size(kpoints)))
  do i=1,size(kpoints)
    if (i/=1 .or. .not. any(sc_acoustic)) then
      do j=1,structure%no_modes
        
        ! generate amplitudes, {(x,V(x))}
        ! generate potential at {q} defined by map
        allocate(amplitudes(2,mapping%count))
        amplitude = -mapping%max/(2*abs(frequencies(j,i)))
        do k=1,mapping%count
          amplitudes(1,k) = amplitude+(k-1)*abs(amplitude/mapping%first)
          amplitudes(2,k) = (energies(k,j,i)-energies(mapping%mid,j,i)) &
                        & / structure_scs(kpoints(i)%sc_id)%sc_size
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
        deallocate(vscf)
        
      enddo
    endif
  enddo
  
  ! ----------------------------------------
  ! Write output.
  ! ----------------------------------------
  
  ! copy acoustic.dat files
  call mkdir(wd//'/anharmonic')
  do i=1,no_supercells
    if (sc_acoustic(i)) then
      call system_call('cp '//wd//'/Supercell_'//i//'/acoustic.dat '// &
         & wd//'/anharmonic')
    endif
  enddo
  
  ! calculate free energy, F(T), for harmonic and anharmonic cases
  ! write output to anharmonic_correction.dat
  result_file = open_write_file(wd//'/anharmonic/anharmonic_correction.dat')
  call calculate_anharmonic(structure,structure_grid,kpoints,Nbasis, &
     & harmonic,eigenvals,result_file)
  close(result_file)
end subroutine

end module
