module bs_quadratic_module

contains

! ----------------------------------------------------------------------
! Program to calculate quadratic band gap correction
! ----------------------------------------------------------------------
subroutine bs_quadratic()
  use constants, only : dp, kB
  use utils,     only : mkdir
  use mapping_module
  use string_module
  use structure_module
  use file_module
  use bands_module
  use displacement_patterns_module
  use supercell_module
  implicit none
  
  ! Parameters
  real(dp), parameter :: dtemperature = 50.0d0
  
  ! User input variables
  integer :: kpoint
  integer :: band
  integer :: degeneracy
  
  ! Previous user inputs
  type(String) :: dft_code      ! The dft code name (castep,vasp,qe)
  type(String) :: seedname      ! The dft input file seedname
  type(String) :: harmonic_path ! The path to the harmonic directory
  
  ! Starting data
  integer              :: no_sc
  integer              :: no_kpoints
  type(MappingData)                 :: mapping
  type(StructureData)               :: structure
  type(StructureData), allocatable  :: structure_scs(:)
  type(DispPatterns)                :: disp_patterns
  
  ! Band data
  type(BandsData)       :: bands
  real(dp)              :: band_energy
  real(dp), allocatable :: frequencies(:,:)
  real(dp), allocatable :: bs(:,:)
  
  ! Kpoint data
  integer, allocatable :: sc_ids(:)
  integer, allocatable :: gvectors(:)
  
  integer, allocatable :: multiplicity(:)
  
  integer, allocatable :: band_refs(:)
  integer              :: ref
  
  ! Temporary variables
  integer        :: i,j,k
  type(String)   :: sdir
  type(String)   :: mdir
  type(String)   :: filename
  integer        :: configs(2)
  integer        :: config
  character(100) :: dump
  
  ! calculate_bs variables
  real(dp)              :: renormalised_band
  real(dp)              :: renormalised_band_kpoint
  real(dp), allocatable :: deformation(:,:)
  real(dp)              :: temperature
  
  ! File contents
  type(String), allocatable :: user_inputs(:)
  
  ! File units
  integer :: no_sc_file
  integer :: list_file
  integer :: bgc_file
  integer :: bck_file
  
  ! --------------------------------------------------
  ! Get user inputs
  ! --------------------------------------------------
  call print_line("What is the k-point of interest?")
  kpoint = int(read_line_from_user())
  
  call print_line("What is the band number of interest &
     &(for the primitive cell)?")
  call print_line("(In case of band degeneracy, the highest band is required)")
  band = int(read_line_from_user())
  
  call print_line("What is the band degeneracy?")
  degeneracy = int(read_line_from_user())
  
  ! --------------------------------------------------
  ! Read basic data
  ! --------------------------------------------------
  user_inputs = read_lines('user_input.txt')
  dft_code = user_inputs(1)
  seedname = user_inputs(2)
  harmonic_path = user_inputs(3)
  
  no_sc_file = open_read_file('no_sc.dat')
  read(no_sc_file,*) no_sc
  close(no_sc_file)
  
  mapping = read_mapping_file('mapping.dat')
  
  structure = read_structure_file(harmonic_path//'/structure.dat')
  
  filename = harmonic_path//'/list.dat'
  no_kpoints = count_lines(filename)
  allocate(multiplicity(no_kpoints))
  allocate(sc_ids(no_kpoints))
  allocate(gvectors(no_kpoints))
  list_file = open_read_file(filename)
  do i=1,no_kpoints
    read(list_file,*) dump,dump,dump,multiplicity(i),sc_ids(i),gvectors(i)
  enddo
  close(list_file)
  
  call mkdir('bs')
  
  ! Obtain relevant band energy for each supercell
  filename = 'Supercell_1/static/kpoint.'//kpoint//'.dat'
  bands = read_castep_bands_file(filename)
  band_energy = bands%bands(1,1)
  
  allocate(band_refs(no_sc))
  do i=1,no_sc
    sdir = 'Supercell_'//i
    ! Obtain relevant band for each supercell
    filename = sdir//'/static/kpoint.'//kpoint//'.dat'
    bands = read_castep_bands_file(filename)
    band_refs(i) = minloc(abs(bands%bands(:,1)-band_energy),dim=1)
  enddo
  
  ! Read in supercell data
  allocate(structure_scs(no_sc))
  do i=1,size(structure_scs)
    sdir = 'Supercell_'//i
    structure_scs(i) = read_structure_file(sdir//'/structure.dat')
  enddo
    
  ! Loop over kpoints
  allocate(frequencies(structure%no_modes,no_kpoints))
  allocate(bs(structure%no_modes,no_kpoints))
  do i=1,no_kpoints
    ! Read frequencies
    filename = sdir//'/lte/disp_patterns.dat'
    disp_patterns = read_disp_patterns_file( filename, &
                                           & structure_scs(sc_ids(i))%no_modes)
    frequencies(:,i) = disp_patterns%frequencies(:,gvectors(i))
    
    ! Read bands
    ref = band_refs(sc_ids(i))
    do j=1,structure%no_modes
      configs = (/mapping%first,mapping%last/)
      do k=1,2
        config = configs(k)
        if (config/=0) then
          mdir='kpoint.'//kpoint//'/mode.'//j//'.'//config
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
  allocate(deformation(structure%no_modes,no_kpoints))
  deformation = bs/mapping%max**2
  
  ! Calculate quadratic vibrational correction
  bgc_file = open_write_file('bs/band_gap_correction.dat')
  bck_file = open_write_file('bs/bg_correction_kp.dat')
  do k=0,20  ! loop over temperature
    renormalised_band=0.0
    temperature = k*dtemperature
    if(temperature<1.d-5)then
      do i=1,no_kpoints
        call print_line(bck_file,'k-point '//i)
        renormalised_band_kpoint=0.0
        do j=1,structure%no_modes
          renormalised_band = renormalised_band &
                          & + deformation(i,j)*multiplicity(i)/no_kpoints
          renormalised_band_kpoint = renormalised_band_kpoint+deformation(i,j)
          call print_line(bck_file,i//' '//j//' '//deformation(i,j))
        enddo
        call print_line(bck_file,i//' '//renormalised_band_kpoint)
      enddo
    else
      do i=1,no_kpoints
        do j=1,structure%no_modes
          renormalised_band = renormalised_band                              &
                          & + deformation(i,j)                               &
                          & * ( 1.0                                          &
                          &   + 2.0                                          &
                          &   / (dexp(frequencies(i,j)/(temperature*kB))-1)) &
                          & * multiplicity(i)                                &
                          & / no_kpoints
        enddo
      enddo
    endif ! temperature
    call print_line(bgc_file,temperature//' '//renormalised_band)
  enddo
  close(bgc_file)
  close(bck_file)
end subroutine
end module
