! ======================================================================
! Program to calculate anharmonic 1-dimensional correction.
! ======================================================================
module anharmonic_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function anharmonic_keywords() result(keywords)
  use help_module
  implicit none
  
  type(KeywordData) :: keywords(2)
  
  keywords = [ &
  & make_keyword('integration_points', '5000', 'integration_points is the &
     &number of points onto which the 1-d potential is interpolated before &
     &VSCF is performed.'),                                                &
  & make_keyword('basis_size', '20', 'basis_size is the number of states &
     &used to construct a basis when performing VSCF.')                    ]
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine anharmonic(arguments)
  use utils_module, only : mkdir, make_dft_output_filename
  use mapping_module
  use structure_module
  use dft_output_file_module
  use supercell_module
  use qpoints_module
  use calculate_anharmonic_module
  use quadratic_spline_module
  use vscf_1d_module
  use dictionary_module
  implicit none
  
  type(Dictionary) :: arguments
  
  ! Working directory.
  type(String) :: wd
  
  ! ----------------------------------------
  ! Input variables
  ! ----------------------------------------
  type(Dictionary) :: setup_quadratic_arguments
  type(String)     :: dft_code
  type(String)     :: seedname
  type(String)     :: harmonic_path
  integer          :: integration_points
  integer          :: basis_size
  
  ! ----------------------------------------
  ! Working variables
  ! ----------------------------------------
  integer               :: no_supercells   ! no. of supercells
  type(MappingData)     :: mapping         ! mapping.dat
  
  type(StructureData)              :: structure
  type(StructureData), allocatable :: structure_scs(:)
  
  ! q-point data.
  type(StructureData)           :: structure_grid
  type(QpointData), allocatable :: qpoints(:)
  
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
  type(String) :: ddir       ! qpoint_*/mode_*/amplitude_*
  
  real(dp), allocatable :: eigenvals(:,:,:)
  real(dp), allocatable :: harmonic(:,:,:)
  
  type(String)          :: filename
  type(String)          :: dft_output_filename
  
  type(DftOutputFile)   :: dft_output_file
  
  ! ----------------------------------------
  ! Temporary variables
  ! ----------------------------------------
  integer :: i,j,k
  integer :: result_code
  
  ! ----------------------------------------
  ! Files.
  ! ----------------------------------------
  type(String), allocatable :: no_sc_file(:)
  integer                   :: result_file
  
  
  ! ----------------------------------------
  ! Read in data
  ! ----------------------------------------
  ! Read arguments
  wd = item(arguments, 'working_directory')
  integration_points = int(item(arguments, 'integration_points'))
  basis_size = int(item(arguments, 'basis_size'))
  
  ! Read in previous user inputs
  setup_quadratic_arguments = read_dictionary_file( &
     & wd//'/setup_quadratic.used_settings')
  dft_code = item(setup_quadratic_arguments, 'dft_code')
  seedname = item(setup_quadratic_arguments, 'seedname')
  harmonic_path = item(setup_quadratic_arguments, 'harmonic_path')
  
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
  
  ! Read qpoints
  structure_grid = read_structure_file(harmonic_path//'/structure_grid.dat')
  qpoints = read_qpoints_file(harmonic_path//'/qpoints_ibz.dat')
  
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
  
  ! Loop over q-points.
  allocate(frequencies(structure%no_modes,size(qpoints)))
  allocate(energies(mapping%count,structure%no_modes,size(qpoints)))
  do i=1,size(qpoints)
    if (.not. sc_acoustic(qpoints(i)%sc_id)) then
      ! Read frequencies
      frequencies_file = read_lines( &
         & harmonic_path//'/qpoint_'//i//'/frequencies.dat')
      frequencies(:,i) = dble(frequencies_file)
      
      do j=1,structure%no_modes
        ! read energies
        do k=mapping%first,mapping%last
          ddir = wd//'/qpoint_'//i//'/mode_'//j//'/amplitude_'//k
          
          dft_output_filename = make_dft_output_filename(dft_code,seedname)
          dft_output_filename = ddir//'/'//dft_output_filename
          
          if (file_exists(dft_output_filename)) then
            dft_output_file = read_dft_output_file(dft_code,dft_output_filename)
            energies(k,j,i) = dft_output_file%energy
          else
            energies(k,j,i) = static_energies(qpoints(i)%sc_id)
          endif
        enddo
      enddo
    endif
  enddo
  
  ! ----------------------------------------
  ! Process data.
  ! ----------------------------------------
  
  ! Calculate anharmonic 1-dimensional correction
  allocate(harmonic(basis_size,structure%no_modes,size(qpoints)))
  allocate(eigenvals(basis_size,structure%no_modes,size(qpoints)))
  do i=1,size(qpoints)
    if (i/=1 .or. .not. any(sc_acoustic)) then
      do j=1,structure%no_modes
        
        ! generate amplitudes, {(x,V(x))}
        ! generate potential at {q} defined by map
        allocate(amplitudes(2,mapping%count))
        amplitude = -mapping%max/(2*abs(frequencies(j,i)))
        do k=1,mapping%count
          amplitudes(1,k) = amplitude+(k-1)*abs(amplitude/mapping%first)
          amplitudes(2,k) = (energies(k,j,i)-energies(mapping%mid,j,i)) &
                        & / structure_scs(qpoints(i)%sc_id)%sc_size
        enddo
        
        ! fit splines
        ! interpolate potential onto integration_points points
        ! min(q) and max(q) are unchanged
        spline = quadratic_spline(integration_points, amplitudes)
        
        ! calculate 1-d anharmonic energy
        vscf = vscf_1d(frequencies(j,i), spline, basis_size)
        
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
      result_code = system_call('cp '//             &
         & wd//'/Supercell_'//i//'/acoustic.dat '// &
         & wd//'/anharmonic')
      call err(result_code==0)
    endif
  enddo
  
  ! calculate free energy, F(T), for harmonic and anharmonic cases
  ! write output to anharmonic_correction.dat
  result_file = open_write_file(wd//'/anharmonic/anharmonic_correction.dat')
  call calculate_anharmonic(structure,structure_grid,qpoints,basis_size, &
     & harmonic,eigenvals,result_file)
  close(result_file)
end subroutine

end module
