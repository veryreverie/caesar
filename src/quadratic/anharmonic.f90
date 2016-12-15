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
  use generate_amplitudes_module,  only : generate_amplitudes2
  use calculate_anharmonic_module, only : calculate_anharmonic2
  use quadratic_spline_module,     only : quadratic_spline2
  use vscf_1d_module,              only : VscfReturn, vscf_1d_2
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
  ! TODO: multiplicity and no_bigpoints might be the same thing.
  integer               :: no_supercells   ! no. of supercells
  integer               :: no_atoms        ! no. atoms in unit cell
  integer               :: no_modes        ! no_atoms*3
  integer, allocatable  :: no_atoms_sc(:)  ! no. atoms in supercell
  integer, allocatable  :: no_cells(:)     ! no_atoms_sc/no_atoms
  character(32)         :: castep          ! seedname.castep
  type(Mapping)         :: map             ! mapping.dat
  integer, allocatable  :: big_points(:)   ! first column of all list.dats
  integer               :: big_point
  integer, allocatable  :: no_bigpoints(:) ! no. lines in Supercell_*/list.dat
  integer, allocatable  :: sizes(:)
  integer, allocatable  :: multiplicity(:) ! fourth column of ibz.dat
  real(dp), allocatable :: energies(:,:,:)
  real(dp)              :: static_energy
  real(dp), allocatable :: frequencies(:,:) ! harmonic frequencies
  real(dp), allocatable :: amplitudes(:,:)
  real(dp), allocatable :: spline(:,:)
  type(VscfReturn)      :: vscf
  logical, allocatable  :: acoustics(:)    ! if Supercell_i/acoustic.dat exists
  logical               :: acoustic        ! if any */acoustic.dat exists
  
  character(32)        :: k_str
  character(32)        :: l_str
  
  real(dp), allocatable :: eigenvals(:,:,:)
  real(dp), allocatable :: harmonic(:,:,:)
  
  character(80)       :: sdir        ! Supercell_*/ directory name
  
  ! ----------------------------------------
  ! Temporary variables
  ! ----------------------------------------
  integer             :: i, j, k, l  ! loop variables
  type(ProcessResult) :: proc        ! temporary process result
  real(dp)            :: temp_real
  character(80)       :: filename
  
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
  ! Read initial data
  ! ----------------------------------------
  ! read the number of Supercell_* directories into no_supercells
  proc = system_process('ls -1d Supercell_* | wc -l')
  read(proc%stdout,*) no_supercells
  
  ! allocate arrays of size no_supercells
  allocate(acoustics(no_supercells))
  allocate(no_atoms_sc(no_supercells))
  allocate(no_cells(no_supercells))
  
  ! read the number of atoms and modes
  equilibrium_file = open_read_file('Supercell_1/equilibrium.dat')
  read(equilibrium_file,*) no_atoms
  close(equilibrium_file)
  no_modes = no_atoms*3
  
  ! read the castep seedname into castep variable
  seedname_file = open_read_file('Supercell_1/seedname.txt')
  read(seedname_file,*) castep
  close(seedname_file)
  castep = trim(castep)//'.castep'
  
  ! read sampling data from mapping.dat
  mapping_file = open_read_file('mapping.dat')
  map = read_mapping(mapping_file)
  close(mapping_file)
  
  ! check for Supercell_*/acoustic.dat
  acoustic=.false.
  do i=1,no_supercells
    sdir = 'Supercell_'//trim(i2s(i))
    if (file_exists(trim(sdir)//'/acoustic.dat')) then
      acoustics(i) = .true.
      acoustic = .true.
      call system('cp '//trim(sdir)//'/acoustic.dat anharmonic')
    else
      acoustics(i) = .false.
    endif
  enddo
  
  ! count big_points from Supercell_*/list.dat
  allocate(no_bigpoints(no_supercells))
  do i=1,no_supercells
    if (acoustics(i)) then
      no_bigpoints(i) = 0
    else
      sdir = 'Supercell_'//trim(i2s(i))
      list_file = open_read_file(trim(sdir)//'/list.dat')
      no_bigpoints(i) = count_lines(list_file)
      close(list_file)
    endif
  enddo
  
  ! allocate arrays with size no_bigpoints
  allocate(big_points(sum(no_bigpoints)))
  allocate(sizes(size(big_points)))
  
  ! read big_points from Supercell_*/list.dat
  j = 1
  do i=1,no_supercells
    if (.not. acoustics(i)) then
      sdir = 'Supercell_'//trim(i2s(i))
      list_file = open_read_file(trim(sdir)//'/list.dat')
      do k=1,no_bigpoints(i)
        read(list_file,*) big_points(j)
        j = j + 1
      enddo
      close(list_file)
    endif
  enddo
  
  ! read no_atoms_sc and no_cells from Supercell_*/super_equilibrium.dat
  do i=1,no_supercells
    if (.not. acoustics(i)) then
      sdir = 'Supercell_'//trim(i2s(i))
      super_file = open_read_file(trim(sdir)//'/super_equilibrium.dat')
      read(super_file,*) no_atoms_sc(i)
      close(super_file)
      no_cells(i) = no_atoms_sc(i)/no_atoms
    endif
  enddo
  
  ! allocate arguments for calculate_anharmonic
  allocate(energies(size(big_points),no_modes,map%count))
  allocate(frequencies(size(big_points),no_modes))
  allocate(harmonic(size(big_points),no_modes,Nbasis))
  allocate(eigenvals(size(big_points),no_modes,Nbasis))
  
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
    if (.not. acoustics(i)) then
      sdir = 'Supercell_'//trim(i2s(i))
      static_energy_file = open_read_file(trim(sdir)//'/static/energy.dat')
      read(static_energy_file,*) static_energy
      close(static_energy_file)
      
      ! j is a cumulative sum of big_points
      do j=j,j+no_bigpoints(i)-1
        big_point = big_points(j)
        
        ! set sizes
        sizes(j) = no_cells(i)
        
        ! read energies
        do k=1,no_modes
          k_str = trim(sdir)//'/&
            &kpoint.'//trim(i2s(big_point))//'/&
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
          k_str = trim(sdir)//'/kpoint.'//trim(i2s(big_point))
          frequency_file = open_read_file(k_str//'/&
            &frequency.'//trim(i2s(big_point))//'.'//trim(i2s(k))//'.dat')
          read(frequency_file,*) frequencies(j,k)
          close(frequency_file)
          
        enddo
      enddo ! loop over big points
    endif
  enddo ! loop over supercells
  
  ! ----------------------------------------
  ! Calculate anharmonic 1-dimensional correction
  ! ----------------------------------------
  
  do i=1,size(big_points)
    if (big_points(i)/=1 .or. .not. acoustic) then
      do j=1,no_modes
        
        ! generate amplitudes
        ! generate potential at {q} defined by map
        amplitudes = generate_amplitudes2(map, energies(i,j,:),&
          &frequencies(i,j),sizes(i))
        
        ! fit splines
        ! interpolate potential onto integration_points points
        ! min(q) and max(q) are unchanged
        spline = quadratic_spline2(integration_points, amplitudes)
        
        ! calculate 1-d anharmonic energy
        vscf = vscf_1d_2(frequencies(i,j), spline, Nbasis)
        
        harmonic(i,j,:) = vscf%harmonic
        eigenvals(i,j,:) = vscf%eigenvals
        
        deallocate(amplitudes)
        deallocate(spline)
        deallocate(vscf%anh_pot)
        deallocate(vscf%hamiltonian)
        deallocate(vscf%harmonic)
        deallocate(vscf%eigenvals)
        deallocate(vscf%eigenvecs)
        
      enddo ! loop over modes
    endif
  enddo ! loop over list.dat
  
  ! ----------------------------------------
  ! calculate free energy, F(T), for harmonic and anharmonic cases
  ! write output to anharmonic_correction.dat
  ! ----------------------------------------
  call system('mkdir anharmonic')
  result_file = open_write_file('anharmonic/anharmonic_correction.dat')
  call calculate_anharmonic2(multiplicity,no_modes,Nbasis,harmonic,&
    &eigenvals,result_file)
  close(result_file)
  
  deallocate(eigenvals)
  deallocate(harmonic)
  deallocate(frequencies)
  deallocate(energies)
  deallocate(multiplicity)
end subroutine

end module
