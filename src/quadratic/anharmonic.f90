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
  use keyword_module
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [                                                                &
  & make_keyword( 'integration_points',                                       &
  &               'integration_points is the number of points onto which the &
  &1-d potential is interpolated before VSCF is performed.',                  &
  &               default_value='5000'),                                      &
  & make_keyword( 'basis_size',                                               &
  &               'basis_size is the number of states used to construct a &
  &basis when performing VSCF.',                                              &
  &               default_value='20') ]
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine anharmonic(arguments)
  use utils_module, only : mkdir
  use ifile_module
  use ofile_module
  use structure_module
  use dft_output_file_module
  use qpoints_module
  use calculate_anharmonic_correction_module
  use quadratic_spline_module
  use vscf_1d_module
  use dictionary_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Working directory.
  type(String) :: wd
  
  ! User inputs.
  integer          :: integration_points
  integer          :: basis_size
  
  ! Previous user inputs.
  type(Dictionary) :: setup_quadratic_arguments
  type(String)     :: dft_code
  type(String)     :: seedname
  type(String)     :: harmonic_path
  integer          :: no_samples
  real(dp)         :: displacement
  
  ! Working variables.
  integer               :: no_supercells   ! no. of supercells
  
  type(StructureData)              :: structure
  type(StructureData), allocatable :: supercells(:)
  
  ! q-point data.
  type(StructureData)           :: structure_grid
  type(QpointData), allocatable :: qpoints_ibz(:)
  
  real(dp), allocatable :: energies(:,:,:)
  real(dp), allocatable :: static_energies(:)
  type(IFile)           :: frequencies_file
  real(dp), allocatable :: frequencies(:,:) ! harmonic frequencies
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
  
  ! Files.
  type(IFile) :: no_supercells_file
  type(OFile) :: result_file
  
  ! Temporary variables.
  integer :: i,j,k
  integer :: result_code
  
  ! --------------------------------------------------
  ! Read in data.
  ! --------------------------------------------------
  ! Read arguments
  wd = arguments%value('working_directory')
  integration_points = int(arguments%value('integration_points'))
  basis_size = int(arguments%value('basis_size'))
  
  ! Read in previous user inputs
  call setup_quadratic_arguments%read_file( &
     & wd//'/setup_quadratic.used_settings')
  dft_code = setup_quadratic_arguments%value('dft_code')
  seedname = setup_quadratic_arguments%value('seedname')
  harmonic_path = setup_quadratic_arguments%value('harmonic_path')
  no_samples = int(setup_quadratic_arguments%value('no_samples'))
  displacement = dble(setup_quadratic_arguments%value('displacement'))
  
  ! read the number of Supercell_* directories into no_supercells
  no_supercells_file = harmonic_path//'/no_supercells.dat'
  no_supercells = int(no_supercells_file%line(1))
  
  ! allocate arrays of size no_supercells
  allocate(sc_acoustic(no_supercells))
  allocate(static_energies(no_supercells))
  
  ! read structure data
  structure = read_structure_file(harmonic_path//'/structure.dat')
  
  ! check for Supercell_*/acoustic.dat
  do i=1,no_supercells
    sc_acoustic(i) = file_exists(wd//'/Supercell_'//i//'/acoustic.dat')
  enddo
  
  ! Read qpoints
  structure_grid = read_structure_file(harmonic_path//'/structure_grid.dat')
  qpoints_ibz = read_qpoints_file(harmonic_path//'/qpoints_ibz.dat')
  
  ! Read supercell structures
  allocate(supercells(no_supercells))
  do i=1,no_supercells
    sdir = wd//'/Supercell_'//i
    filename = harmonic_path//'/'//sdir//'/structure.dat' 
    supercells(i) = read_structure_file(filename)
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
  allocate(frequencies(structure%no_modes,size(qpoints_ibz)))
  allocate(energies(2*no_samples+1,structure%no_modes,size(qpoints_ibz)))
  do i=1,size(qpoints_ibz)
    if (.not. sc_acoustic(qpoints_ibz(i)%sc_id)) then
      ! Read frequencies
      frequencies_file =  harmonic_path//'/qpoint_'//i//'/frequencies.dat'
      do j=1,size(frequencies_file)
        frequencies(j,i) = dble(frequencies_file%line(j))
      enddo
      
      do j=1,structure%no_modes
        ! read energies
        do k=-no_samples,no_samples
          ddir = wd//'/qpoint_'//i//'/mode_'//j//'/amplitude_'//k
          
          dft_output_filename = make_dft_output_filename(dft_code,seedname)
          dft_output_filename = ddir//'/'//dft_output_filename
          
          if (file_exists(dft_output_filename)) then
            dft_output_file = read_dft_output_file(dft_code,dft_output_filename)
            energies(k,j,i) = dft_output_file%energy
          else
            energies(k,j,i) = static_energies(qpoints_ibz(i)%sc_id)
          endif
        enddo
      enddo
    endif
  enddo
  
  ! ----------------------------------------
  ! Process data.
  ! ----------------------------------------
  
  ! Calculate anharmonic 1-dimensional correction
  allocate(harmonic(basis_size,structure%no_modes,size(qpoints_ibz)))
  allocate(eigenvals(basis_size,structure%no_modes,size(qpoints_ibz)))
  do i=1,size(qpoints_ibz)
    if (i/=1 .or. .not. any(sc_acoustic)) then
      do j=1,structure%no_modes
        
        ! generate amplitudes, {(x,V(x))}
        ! generate potential at {u}
        allocate(amplitudes(2,2*no_samples+1))
        do k=1,2*no_samples+1
          amplitudes(1,k) = displacement*(k-no_samples-1)/real(no_samples,dp)
          amplitudes(2,k) = (energies(k,j,i)-energies(no_samples+1,j,i)) &
                        & / supercells(qpoints_ibz(i)%sc_id)%sc_size
        enddo
        
        ! fit splines
        ! interpolate potential onto integration_points points
        ! displacement is unchanged
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
      if (result_code/=0) then
        call err()
      endif
    endif
  enddo
  
  ! calculate free energy, F(T), for harmonic and anharmonic cases
  ! write output to anharmonic_correction.dat
  result_file = wd//'/anharmonic/anharmonic_correction.dat'
  call calculate_anharmonic_correction( structure,      &
                                      & structure_grid, &
                                      & qpoints_ibz,    &
                                      & basis_size,     &
                                      & harmonic,       &
                                      & eigenvals,      &
                                      & result_file)
end subroutine

end module
