module bs_quadratic_module

contains

! ----------------------------------------------------------------------
! Program to calculate quadratic band gap correction
! ----------------------------------------------------------------------
subroutine bs_quadratic()
  use constants, only : dp,kB
  use mapping_module
  use string_module
  use structure_module
  use file_module
  use bands_module
  use displacement_patterns_module
  implicit none
  
  ! Parameters
  real(dp), parameter :: dtemperature = 50.0d0
  
  ! User input variables
  integer :: kpoint
  integer :: band
  integer :: degeneracy
  
  ! Starting data
  type(String)        :: harmonic_path
  type(String)        :: seedname
  integer             :: no_sc
  integer             :: no_kpoints
  type(MappingData)   :: mapping
  type(StructureData) :: structure
  type(StructureData) :: structure_sc
  type(DispPatterns)  :: disp_patterns
  
  ! Band data
  type(BandsData)       :: bands
  real(dp)              :: band_energy
  real(dp), allocatable :: frequencies(:,:)
  real(dp), allocatable :: bs(:,:)
  
  ! Kpoint data
  integer, allocatable :: kpoints(:)
  integer, allocatable :: gvectors(:)
  integer, allocatable :: sc_ids(:)
  
  integer, allocatable :: multiplicity(:)
  
  integer, allocatable :: band_refs(:)
  integer              :: ref
  
  ! File units
  integer :: harmonic_path_file
  integer :: no_sc_file
  integer :: seedname_file
  integer :: list_file
  integer :: ibz_file
  integer :: bgc_file
  integer :: bck_file
  
  ! Temporary variables
  character(100) :: temp_char
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
  
  ! --------------------------------------------------
  ! Get user inputs
  ! --------------------------------------------------
  write(*,"(a)") "What is the k-point of interest?"
  read(*,*) kpoint
  
  write(*,"(a)") "What is the band number of interest &
     &(for the primitive cell)?"
  write(*,"(a)") "(In case of band degeneracy, the highest band is required)"
  read(*,*) band
  
  write(*,"(a)") "What is the band degeneracy?"
  read(*,*) degeneracy
  
  ! --------------------------------------------------
  ! Read basic data
  ! --------------------------------------------------
  harmonic_path_file = open_read_file('harmonic_path.dat')
  read(harmonic_path_file,*) temp_char
  close(harmonic_path_file)
  harmonic_path = temp_char
  
  seedname_file = open_read_file('seedname.txt')
  read(seedname_file,*) temp_char
  close(seedname_file)
  seedname = temp_char
  
  no_sc_file = open_read_file('no_sc.dat')
  read(no_sc_file,*) no_sc
  close(no_sc_file)
  
  no_kpoints = count_lines(harmonic_path//'/ibz.dat')
  
  mapping = read_mapping_file('mapping.dat')
  
  structure = read_structure_file(harmonic_path//'/structure.dat')
  
  filename = harmonic_path//'/list.dat'
  no_kpoints = count_lines(filename)
  allocate(kpoints(no_kpoints))
  allocate(gvectors(no_kpoints))
  allocate(sc_ids(no_kpoints))
  list_file = open_read_file(filename)
  do i=1,no_kpoints
    read(list_file,*) kpoints(i),gvectors(i),sc_ids(i)
  enddo
  close(list_file)
  
  call system('mkdir bs')
  
  ! Obtain relevant band energy for each supercell
  filename = str('Supercell_1/static/kpoint.')//kpoint//'.dat'
  bands = read_castep_bands_file(filename)
  band_energy = bands%bands(1,1)
  
  allocate(band_refs(no_sc))
  do i=1,no_sc
    ! Obtain relevant band for each supercell
    filename = str('Supercell_')//i//'/static/kpoint.'//kpoint//'.dat'
    bands = read_castep_bands_file(filename)
    band_refs(i) = minloc(abs(bands%bands(:,1)-band_energy),dim=1)
  enddo
    
  ! Loop over kpoints
  allocate(frequencies(structure%no_modes,no_kpoints))
  allocate(bs(structure%no_modes,no_kpoints))
  do i=1,no_kpoints
    ! Read in supercell structure
    structure_sc = read_structure_file(sdir//'/structure.dat')
    
    ! Read frequencies
    filename = sdir//'/lte/disp_patterns.dat'
    disp_patterns = read_disp_patterns_file(filename,structure_sc%no_modes)
    frequencies(:,i) = disp_patterns%frequencies(:,gvectors(i))
    
    ! Read bands
    ref = band_refs(sc_ids(i))
    do j=1,structure%no_modes
      configs = (/mapping%first,mapping%last/)
      do k=1,2
        config = configs(k)
        if (config/=0) then
          mdir=str('kpoint.')//kpoint//'/configurations/mode.'//j//'.'//config
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
  
  ! read ibz file
  allocate(multiplicity(no_kpoints))
  ibz_file = open_read_file(harmonic_path//'/ibz.dat')
  do i=1,no_kpoints 
    read(ibz_file,*)dump,dump,dump,multiplicity(i)
  enddo
  close(ibz_file)
  
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
        write(bck_file,*)'k-point',i
        renormalised_band_kpoint=0.0
        do j=1,structure%no_modes
          renormalised_band = renormalised_band &
                          & + deformation(i,j)*multiplicity(i)/no_kpoints
          renormalised_band_kpoint = renormalised_band_kpoint+deformation(i,j)
          write(bck_file,*)i,j,deformation(i,j)
        enddo
        write(bck_file,*)i,renormalised_band_kpoint
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
    write(bgc_file,*)temperature,renormalised_band 
  enddo
  close(bgc_file)
  close(bck_file)
end subroutine
end module
