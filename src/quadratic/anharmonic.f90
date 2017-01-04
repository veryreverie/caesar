! ------------------------------------------------------------
! program to calculate anharmonic 1-dimensional correction
! ------------------------------------------------------------

module anharmonic_module
contains

subroutine anharmonic()
  use constants,                   only : dp
  use utils,                       only : i2s, file_exists, count_lines
  use process,                     only : ProcessResult, system_process
  use mapping_module,              only : Mapping, read_mapping
  use generate_amplitudes_module,  only : generate_amplitudes
  use calculate_anharmonic_module, only : calculate_anharmonic
  use quadratic_spline_module,     only : quadratic_spline
  use vscf_1d_module,              only : VscfReturn, vscf_1d, drop
  use file_io,                     only : open_read_file, open_write_file
  implicit none
  
  ! ----------------------------------------
  ! Parameters
  ! ----------------------------------------
  integer, parameter :: integration_points = 5000
  integer, parameter :: Nbasis = 20
  
  ! ----------------------------------------
  ! Working variables
  ! ----------------------------------------
  ! TODO: multiplicity and sc_n_kpoints might be the same thing.
  integer               :: no_supercells   ! no. of supercells
  integer               :: no_atoms        ! no. atoms in unit cell
  integer               :: no_modes        ! no_atoms*3
  integer, allocatable  :: no_atoms_sc(:)  ! no. atoms in supercell
  integer, allocatable  :: no_cells(:)     ! no_atoms_sc/no_atoms
  character(100)        :: castep          ! seedname.castep
  type(Mapping)         :: map             ! mapping.dat
  integer, allocatable  :: kpoints(:)      ! first column of all list.dats
  integer               :: kpoint
  integer, allocatable  :: sc_n_kpoints(:) ! no. lines in Supercell_*/list.dat
  integer, allocatable  :: sizes(:)
  integer, allocatable  :: multiplicity(:) ! fourth column of ibz.dat
  real(dp), allocatable :: energies(:,:,:)
  real(dp)              :: static_energy
  real(dp), allocatable :: frequencies(:,:) ! harmonic frequencies
  real(dp), allocatable :: amplitudes(:,:)
  real(dp), allocatable :: spline(:,:)
  type(VscfReturn)      :: vscf
  logical, allocatable  :: sc_acoustic(:)  ! if Supercell_i/acoustic.dat exists
  
  character(100)        :: k_str
  character(100)        :: l_str
  
  real(dp), allocatable :: eigenvals(:,:,:)
  real(dp), allocatable :: harmonic(:,:,:)
  
  character(100)        :: sdir        ! Supercell_*/ directory name
  
  ! ----------------------------------------
  ! Temporary variables
  ! ----------------------------------------
  integer             :: i, j, k, l  ! loop variables
  type(ProcessResult) :: proc        ! temporary process result
  real(dp)            :: temp_real
  character(100)      :: filename
  
  ! ----------------------------------------
  ! File units
  ! ----------------------------------------
  integer :: equilibrium_file
  integer :: seedname_file
  integer :: mapping_file
  integer :: list_file
  integer :: super_file
  integer :: energy_file
  integer :: frequency_file
  integer :: ibz_file
  integer :: result_file
  integer :: static_energy_file
  
  ! ----------------------------------------
  ! Read in data
  ! ----------------------------------------
  
  ! read the number of Supercell_* directories into no_supercells
  proc = system_process('ls -1d Supercell_* | wc -l')
  read(proc%stdout,*) no_supercells
  
  ! allocate arrays of size no_supercells
  allocate(sc_acoustic(no_supercells))
  allocate(no_atoms_sc(no_supercells))
  allocate(no_cells(no_supercells))
  
  ! read the number of atoms and modes
  equilibrium_file = open_read_file('Supercell_1/equilibrium.dat')
  read(equilibrium_file,*) no_atoms
  close(equilibrium_file)
  no_modes = no_atoms*3
  
  ! read the castep seedname into castep variable
  seedname_file = open_read_file('seedname.txt')
  read(seedname_file,*) castep
  close(seedname_file)
  castep = trim(castep)//'.castep'
  
  ! read sampling data from mapping.dat
  mapping_file = open_read_file('mapping.dat')
  map = read_mapping(mapping_file)
  close(mapping_file)
  
  ! check for Supercell_*/acoustic.dat
  do i=1,no_supercells
    sc_acoustic(i) = file_exists('Supercell_'//trim(i2s(i))//'/acoustic.dat')
  enddo
  
  ! count kpoints from Supercell_*/list.dat
  allocate(sc_n_kpoints(no_supercells))
  do i=1,no_supercells
    if (sc_acoustic(i)) then
      sc_n_kpoints(i) = 0
    else
      list_file = open_read_file('Supercell_'//trim(i2s(i))//'/list.dat')
      sc_n_kpoints(i) = count_lines(list_file)
      close(list_file)
    endif
  enddo
  
  ! allocate arrays with size sc_n_kpoints
  allocate(kpoints(sum(sc_n_kpoints)))
  allocate(sizes(size(kpoints)))
  
  ! read kpoints from Supercell_*/list.dat
  j = 1
  do i=1,no_supercells
    if (.not. sc_acoustic(i)) then
      list_file = open_read_file('Supercell_'//trim(i2s(i))//'/list.dat')
      do k=1,sc_n_kpoints(i)
        read(list_file,*) kpoints(j)
        j = j + 1
      enddo
      close(list_file)
    endif
  enddo
  
  ! read no_atoms_sc and no_cells from Supercell_*/super_equilibrium.dat
  do i=1,no_supercells
    if (.not. sc_acoustic(i)) then
      sdir = 'Supercell_'//trim(i2s(i))
      super_file = open_read_file(trim(sdir)//'/super_equilibrium.dat')
      read(super_file,*) no_atoms_sc(i)
      close(super_file)
      no_cells(i) = no_atoms_sc(i)/no_atoms
    endif
  enddo
  
  ! allocate arguments for calculate_anharmonic
  allocate(energies(size(kpoints),no_modes,map%count))
  allocate(frequencies(size(kpoints),no_modes))
  allocate(harmonic(size(kpoints),no_modes,Nbasis))
  allocate(eigenvals(size(kpoints),no_modes,Nbasis))
  
  ! read multiplicity from ibz.dat
  ibz_file = open_read_file('ibz.dat')
  allocate(multiplicity(count_lines(ibz_file)))
  do i=1,size(multiplicity)
    read(ibz_file,*) temp_real, temp_real, temp_real, multiplicity(i)
  enddo
  close(ibz_file)
  
  ! read data from supercells
  j=1
  do i=1,no_supercells
    if (.not. sc_acoustic(i)) then
      sdir = 'Supercell_'//trim(i2s(i))
      static_energy_file = open_read_file(trim(sdir)//'/static/energy.dat')
      read(static_energy_file,*) static_energy
      close(static_energy_file)
      
      ! j is a cumulative sum of kpoints
      do j=j,j+sc_n_kpoints(i)-1
        kpoint = kpoints(j)
        
        ! set sizes
        sizes(j) = no_cells(i)
        
        ! read energies
        do k=1,no_modes
          k_str = trim(sdir)//'/&
            &kpoint.'//trim(i2s(kpoint))//'/&
            &configurations/&
            &mode.'//trim(i2s(k))//'.'
          filename = trim(k_str)//trim(i2s(map%first))//'/'//trim(castep)
          if (file_exists(filename)) then
            do l=map%first,map%last
              l_str = trim(k_str)//trim(i2s(l))
              if (file_exists(trim(l_str)//'/'//trim(castep))) then
                energy_file = open_read_file(trim(l_str)//'/energy.dat')
                read(energy_file,*) energies(j,k,l)
                close(energy_file)
              else
                energies(j,k,l) = static_energy
              endif
            enddo ! loop over sampling points per mode
          endif
          
          ! read frequencies
          k_str = trim(sdir)//'/kpoint.'//trim(i2s(kpoint))
          frequency_file = open_read_file(k_str//'/&
            &frequency.'//trim(i2s(kpoint))//'.'//trim(i2s(k))//'.dat')
          read(frequency_file,*) frequencies(j,k)
          close(frequency_file)
          
        enddo
      enddo ! loop over big points
    endif
  enddo ! loop over supercells
  
  ! ----------------------------------------
  ! Process data.
  ! ----------------------------------------
  
  ! Calculate anharmonic 1-dimensional correction
  do i=1,size(kpoints)
    if (kpoints(i)/=1 .or. .not. any(sc_acoustic)) then
      do j=1,no_modes
        
        ! generate amplitudes
        ! generate potential at {q} defined by map
        amplitudes = generate_amplitudes(map, energies(i,j,:),&
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
  call calculate_anharmonic(multiplicity,no_modes,Nbasis,harmonic,&
    &eigenvals,result_file)
  close(result_file)
  
  deallocate(eigenvals)
  deallocate(harmonic)
  deallocate(frequencies)
  deallocate(energies)
  deallocate(multiplicity)
end subroutine

end module
