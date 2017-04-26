! ======================================================================
! Program to calculate quadratic band gap correction.
! ======================================================================
module bs_quadratic_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function bs_quadratic_keywords() result(keywords)
  use help_module
  implicit none
  
  type(KeywordData) :: keywords(3)
  
  keywords = [ &
  & make_keyword('q-point', NO_ARGUMENT, 'q-point is the id of the q-point of &
     &interest.'),                                                            &
  & make_keyword('band', NO_ARGUMENT, 'band is the id of the band of interest.&
     &If the band is degenerate, the highest band is required.'),             &
  & make_keyword('degeneracy', '1', 'degeneracy is the number of degenerate &
     &bands.')                                                                ]
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine bs_quadratic(arguments)
  use constants_module, only : kb_in_au
  use utils_module,     only : mkdir
  use mapping_module
  use structure_module
  use bands_module
  use supercell_module
  use dictionary_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Working directory.
  type(String) :: wd
  
  ! Parameters.
  real(dp), parameter :: dtemperature = 50.0_dp
  
  ! User input variables.
  integer :: qpoint
  integer :: band
  integer :: degeneracy
  
  ! Previous user inputs.
  type(String) :: dft_code      ! The dft code name (castep,vasp,qe)
  type(String) :: seedname      ! The dft input file seedname
  type(String) :: harmonic_path ! The path to the harmonic directory
  
  ! Starting data.
  type(Dictionary)  :: setup_quadratic_arguments
  integer              :: no_sc
  integer              :: no_qpoints
  type(MappingData)                 :: mapping
  type(StructureData)               :: structure
  type(StructureData), allocatable  :: structure_scs(:)
  
  ! Band data.
  type(BandsData)           :: bands
  real(dp)                  :: band_energy
  type(String), allocatable :: frequencies_file(:)
  real(dp), allocatable     :: frequencies(:,:)
  real(dp), allocatable     :: bs(:,:)
  
  ! Qpoint data.
  integer, allocatable :: sc_ids(:)
  integer, allocatable :: gvectors(:)
  integer, allocatable :: multiplicity(:)
  integer, allocatable :: band_refs(:)
  integer              :: ref
  
  ! Temporary variables.
  integer        :: i,j,k
  type(String)   :: sdir
  type(String)   :: mdir
  type(String)   :: filename
  integer        :: amplitudes(2)
  integer        :: amplitude
  
  ! calculate_bs variables.
  real(dp)              :: renormalised_band
  real(dp)              :: renormalised_band_qpoint
  real(dp), allocatable :: deformation(:,:)
  real(dp)              :: temperature
  
  ! File contents.
  type(String), allocatable :: no_sc_file(:)
  type(String), allocatable :: ibz_file(:)
  type(String), allocatable :: line(:)
  
  ! File units.
  integer :: bgc_file
  integer :: bck_file
  
  ! --------------------------------------------------
  ! Get user inputs
  ! --------------------------------------------------
  wd = item(arguments, 'working_directory')
  qpoint = int(item(arguments, 'q-point'))
  band = int(item(arguments, 'band'))
  degeneracy = int(item(arguments, 'degeneracy'))
  
  ! --------------------------------------------------
  ! Read basic data
  ! --------------------------------------------------
  setup_quadratic_arguments = read_dictionary_file( &
     & wd//'/setup_quadratic.used_settings')
  dft_code = item(setup_quadratic_arguments, 'dft_code')
  seedname = item(setup_quadratic_arguments, 'seedname')
  harmonic_path = item(setup_quadratic_arguments, 'harmonic_path')
  
  no_sc_file = read_lines(wd//'/no_sc.dat')
  no_sc = int(no_sc_file(1))
  
  mapping = read_mapping_file(wd//'/mapping.dat')
  
  structure = read_structure_file(harmonic_path//'/structure.dat')
  
  ibz_file = read_lines(harmonic_path//'/qpoints_ibz.dat')
  no_qpoints = size(ibz_file)
  allocate(multiplicity(no_qpoints))
  allocate(sc_ids(no_qpoints))
  allocate(gvectors(no_qpoints))
  do i=1,no_qpoints
    line = split(ibz_file(i))
    multiplicity(i) = int(line(4))
    sc_ids(i) = int(line(5))
    gvectors(i) = int(line(6))
  enddo
  
  call mkdir(wd//'/bs')
  
  ! Obtain relevant band energy for each supercell
  filename = wd//'/Supercell_1/qpoint.'//qpoint//'.dat'
  bands = read_castep_bands_file(filename)
  band_energy = bands%bands(1,1)
  
  allocate(band_refs(no_sc))
  do i=1,no_sc
    sdir = wd//'/Supercell_'//i
    ! Obtain relevant band for each supercell
    filename = sdir//'/qpoint.'//qpoint//'.dat'
    bands = read_castep_bands_file(filename)
    band_refs(i) = minloc(abs(bands%bands(:,1)-band_energy),dim=1)
  enddo
  
  ! Read in supercell data
  allocate(structure_scs(no_sc))
  do i=1,size(structure_scs)
    sdir = wd//'/Supercell_'//i
    structure_scs(i) = read_structure_file(sdir//'/structure.dat')
  enddo
  
  ! Loop over q-points
  allocate(frequencies(structure%no_modes,no_qpoints))
  allocate(bs(structure%no_modes,no_qpoints))
  do i=1,no_qpoints
    ! Read frequencies
    frequencies_file = read_lines( &
       & harmonic_path//'/qpoint_'//i//'/frequencies.dat')
    frequencies(:,i) = dble(frequencies_file)
    
    ! Read bands
    ref = band_refs(sc_ids(i))
    do j=1,structure%no_modes
      amplitudes = [ mapping%first, mapping%last ]
      do k=1,2
        amplitude = amplitudes(k)
        if (amplitude/=0) then
          mdir = wd//'/qpoint_'//qpoint//'/mode_'//j//'/amplitude_'//amplitude
          if (file_exists(mdir//'/'//seedname//'.castep')) then
            bands = read_castep_bands_file(mdir//'/'//seedname//'.bands')
            bs(j,i) = bs(j,i)                                  &
                  & + sum(bands%bands(ref-degeneracy+1:ref,1)) &
                  & / (2.d0*degeneracy)
          endif
        endif
      enddo
      bs(j,i) = bs(j,i) - band_energy
    enddo
  enddo
  
  ! Calculate deformation potential
  allocate(deformation(structure%no_modes,no_qpoints))
  deformation = bs/mapping%max**2
  
  ! Calculate quadratic vibrational correction
  bgc_file = open_write_file(wd//'/bs/band_gap_correction.dat')
  bck_file = open_write_file(wd//'/bs/bg_correction_kp.dat')
  do k=0,20  ! loop over temperature
    renormalised_band=0.0
    temperature = k*dtemperature
    if(temperature<1.d-5)then
      do i=1,no_qpoints
        call print_line(bck_file,'q-point '//i)
        renormalised_band_qpoint=0.0
        do j=1,structure%no_modes
          renormalised_band = renormalised_band &
                          & + deformation(i,j)*multiplicity(i)/no_qpoints
          renormalised_band_qpoint = renormalised_band_qpoint+deformation(i,j)
          call print_line(bck_file,i//' '//j//' '//deformation(i,j))
        enddo
        call print_line(bck_file,i//' '//renormalised_band_qpoint)
      enddo
    else
      do i=1,no_qpoints
        do j=1,structure%no_modes
          renormalised_band = renormalised_band                               &
                          & + deformation(i,j)                                &
                          & * ( 1                                             &
                          &   + 2                                             &
                          &   / (exp(frequencies(i,j)/(temperature*kb_in_au)) &
                          &     - 1))                                         &
                          & * multiplicity(i)                                 &
                          & / no_qpoints
        enddo
      enddo
    endif ! temperature
    call print_line(bgc_file,temperature//' '//renormalised_band)
  enddo
  close(bgc_file)
  close(bck_file)
end subroutine
end module
