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
  use keyword_module
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [                                                                &
  & make_keyword( 'q-point',                                                  &
  &               'q-point is the id of the q-point of interest.'),           &
  & make_keyword( 'band',                                                     &
  &               'band is the id of the band of interest. If the band is &
  &degenerate, the highest band is required.'),                               &
  & make_keyword( 'degeneracy',                                               &
  &               'degeneracy is the number of degenerate bands.',            &
  &               default_value='1') ]
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine bs_quadratic(arguments)
  use constants_module, only : kb_in_au
  use utils_module,     only : mkdir
  use ifile_module
  use ofile_module
  use structure_module
  use bands_module
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
  integer      :: no_samples
  real(dp)     :: displacement
  
  ! Starting data.
  type(Dictionary)  :: setup_quadratic_arguments
  integer              :: no_supercells
  integer              :: no_qpoints
  type(StructureData)               :: structure
  type(StructureData), allocatable  :: supercells(:)
  
  ! Band data.
  type(BandsData)       :: bands
  real(dp)              :: band_energy
  type(IFile)           :: frequencies_file
  real(dp), allocatable :: frequencies(:,:)
  real(dp), allocatable :: bs(:,:)
  
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
  type(IFile)               :: no_supercells_file
  type(IFile)               :: ibz_file
  type(String), allocatable :: line(:)
  
  ! File units.
  type(OFile) :: bgc_file
  type(OFile) :: bck_file
  
  ! --------------------------------------------------
  ! Get user inputs
  ! --------------------------------------------------
  wd = arguments%value('working_directory')
  qpoint = int(arguments%value('q-point'))
  band = int(arguments%value('band'))
  degeneracy = int(arguments%value('degeneracy'))
  
  ! --------------------------------------------------
  ! Read basic data
  ! --------------------------------------------------
  call setup_quadratic_arguments%read_file( &
     & wd//'/setup_quadratic.used_settings')
  dft_code = setup_quadratic_arguments%value('dft_code')
  seedname = setup_quadratic_arguments%value('seedname')
  harmonic_path = setup_quadratic_arguments%value('harmonic_path')
  no_samples = int(setup_quadratic_arguments%value('no_samples'))
  displacement = dble(setup_quadratic_arguments%value('displacement'))
  
  no_supercells_file = wd//'/no_supercells.dat'
  no_supercells = int(no_supercells_file%line(1))
  
  structure = read_structure_file(harmonic_path//'/structure.dat')
  
  ibz_file = harmonic_path//'/qpoints_ibz.dat'
  no_qpoints = size(ibz_file)
  allocate(multiplicity(no_qpoints))
  allocate(sc_ids(no_qpoints))
  allocate(gvectors(no_qpoints))
  do i=1,no_qpoints
    line = split(ibz_file%line(i))
    multiplicity(i) = int(line(4))
    sc_ids(i) = int(line(5))
    gvectors(i) = int(line(6))
  enddo
  
  call mkdir(wd//'/bs')
  
  ! Obtain relevant band energy for each supercell
  filename = wd//'/Supercell_1/qpoint.'//qpoint//'.dat'
  bands = read_castep_bands_file(filename)
  band_energy = bands%bands(1,1)
  
  allocate(band_refs(no_supercells))
  do i=1,no_supercells
    sdir = wd//'/Supercell_'//i
    ! Obtain relevant band for each supercell
    filename = sdir//'/qpoint.'//qpoint//'.dat'
    bands = read_castep_bands_file(filename)
    band_refs(i) = minloc(abs(bands%bands(:,1)-band_energy),dim=1)
  enddo
  
  ! Read in supercell data
  allocate(supercells(no_supercells))
  do i=1,size(supercells)
    sdir = wd//'/Supercell_'//i
    supercells(i) = read_structure_file(sdir//'/structure.dat')
  enddo
  
  ! Loop over q-points
  allocate(frequencies(structure%no_modes,no_qpoints))
  allocate(bs(structure%no_modes,no_qpoints))
  do i=1,no_qpoints
    ! Read frequencies
    frequencies_file = harmonic_path//'/qpoint_'//i//'/frequencies.dat'
    do j=1,size(frequencies_file)
      frequencies(j,i) = dble(frequencies_file%line(j))
    enddo
    
    ! Read bands
    ref = band_refs(sc_ids(i))
    do j=1,structure%no_modes
      amplitudes = [ -no_samples, no_samples ]
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
  deformation = bs/displacement**2
  
  ! Calculate quadratic vibrational correction
  bgc_file = wd//'/bs/band_gap_correction.dat'
  bck_file = wd//'/bs/bg_correction_kp.dat'
  do k=0,20  ! loop over temperature
    renormalised_band=0.0
    temperature = k*dtemperature
    if(temperature<1.d-5)then
      do i=1,no_qpoints
        call bck_file%print_line('q-point '//i)
        renormalised_band_qpoint=0.0
        do j=1,structure%no_modes
          renormalised_band = renormalised_band &
                          & + deformation(i,j)*multiplicity(i)/no_qpoints
          renormalised_band_qpoint = renormalised_band_qpoint+deformation(i,j)
          call bck_file%print_line(i//' '//j//' '//deformation(i,j))
        enddo
        call bck_file%print_line(i//' '//renormalised_band_qpoint)
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
    call bgc_file%print_line(temperature//' '//renormalised_band)
  enddo
end subroutine
end module
