! ------------------------------------------------------------
! program to calculate anharmonic 1-dimensional correction
! ------------------------------------------------------------

module anharmonic_module
contains

subroutine anharmonic()
  use constants,                   only : dp
  use utils,                       only : i2s, file_exists, count_lines
  use process,                     only : ProcessResult, system_process
  use supercell_module,            only : Supercell, read_supercell
  use mapping_module,              only : Mapping, read_mapping
  use generate_amplitudes_module,  only : generate_amplitudes2
  use calculate_anharmonic_module, only : calculate_anharmonic2
  use quadratic_spline_module,     only : quadratic_spline2
  use vscf_1d_module,              only : VscfReturn, vscf_1d_2
  implicit none
  
  ! --------------------
  ! parameters
  ! --------------------
  integer, parameter :: integration_points = 5000
  integer, parameter :: Nbasis = 20
  
  ! --------------------
  ! working variables
  ! --------------------
  integer                      :: no_sc           ! no. of supercells
  character(32)                :: castep          ! seedname.castep
  type(Mapping)                :: map             ! mapping.dat
  type(Supercell), allocatable :: supercells(:)   ! Supercell_* info
  integer                      :: big_point
  integer                      :: no_kpoints      ! number of lines in ibz.dat
  integer, allocatable         :: multiplicity(:) ! fourth column of ibz.dat
  real(dp), allocatable        :: energy(:)
  real(dp)                     :: frequency
  integer                      :: cell_size
  real(dp), allocatable        :: amplitudes(:,:)
  real(dp), allocatable        :: spline(:,:)
  type(VscfReturn)             :: vscf
  logical                      :: acoustic
  
  character(32)        :: k_str
  character(32)        :: l_str
  
  ! --------------------
  ! temporary variables
  ! --------------------
  integer             :: i, j, k, l  ! loop variables
  type(ProcessResult) :: proc        ! temporary process result
  logical             :: temp_bool
  integer             :: no_lines    ! length of file
  real(dp)            :: temp_real
  character(80)       :: temp_string
  character(80)       :: filename
  character(80)       :: efile       ! energy file name
  character(80)       :: sdir        ! Supercell_*/ directory name
  
  ! read the number of Supercell_* directories into no_sc
  proc = system_process('ls -1d Supercell_* | wc -l')
  read(proc%stdout,*) no_sc
  allocate(supercells(no_sc))
  
  ! read the castep seedname into castep variable
  open(100, file='Supercell_1/seedname.txt', status='old', action='read')
  read(100,*) castep
  close(100)
  castep = trim(castep)//'.castep'
  
  ! read sampling data from mapping.dat
  map = read_mapping('mapping.dat')
  
  ! make anharmonic directory, and copy mapping.dat into it
  call system('mkdir anharmonic')
  call system('cp mapping.dat anharmonic')
  
  ! loop over supercells
  do i=1,no_sc
    sdir = 'Supercell_'//trim(i2s(i))
    call system('cat '//trim(sdir)//'/list.dat >> anharmonic/list.dat')
    
    write(*,*)'Working with Supercell ',trim(i2s(i))
    
    if (file_exists(trim(sdir)//'/acoustic.dat')) then
      acoustic = .true.
      call system('cp '//trim(sdir)//'/acoustic.dat anharmonic')
    else
      acoustic = .false.
      ! read supercell info from equilibrium.dat and super_equilibrium.dat
      filename = trim(sdir)//'/equilibrium.dat'
      open(unit=100, file=filename, status='old', action='read')
      filename = trim(sdir)//'/super_equilibrium.dat'
      open(unit=101, file=filename, status='old', action='read')
      supercells(i) = read_supercell(100, 101)
      close(101)
      close(100)
      
      ! read list.dat line by line
      filename = trim(sdir)//'/list.dat'
      open(unit=100, file=filename, status='old', action='read')
      no_lines = count_lines(100)
      do j=1,no_lines
        read(100,*) big_point
        open(unit=101,&
            &file='anharmonic/size.'//trim(i2s(big_point))//'.dat',&
            &action='write')
        write(101,*) supercells(i)%no_cells
        close(101)
        
        do k=1,supercells(i)%no_modes
          k_str = trim(sdir)//'/&
            &kpoint.'//trim(i2s(big_point))//'/&
            &configurations/&
            &mode.'//trim(i2s(j))//'.'
          filename = trim(k_str)//trim(i2s(map%first))//'/'//trim(castep)
          if (file_exists(filename)) then
            efile = 'anharmonic/&
              &energy.'//trim(i2s(big_point))//'.'//trim(i2s(j))//'.dat'
            do l=map%first,map%last
              l_str = trim(k_str)//trim(i2s(k))
              if (file_exists(trim(l_str)//'/'//trim(castep))) then
                call system('cat '//trim(l_str)//'/energy.dat >>'//efile)
              elseif (l==0) then ! TODO: bug? l==map%first?
                call system('cat '//trim(sdir)//'/static/energy.dat >> &
                  &'//efile)
              endif
            enddo ! loop over sampling points per mode
          endif
        enddo ! loop over modes
        
        call system('cp '//&
          &trim(sdir)//'/kpoint.'//trim(i2s(big_point))//'/frequency.*.dat &
          &anharmonic/')
      enddo ! loop over list.dat
      close(100)
      
    endif
  enddo
  
  ! Calculate anharmonic 1-dimensional correction
  call system('cp ibz.dat anharmonic')
  call system('cd anharmonic')
  
  ! read list.dat line by line
  open(unit=100, file='list.dat', status='old', action='read')
  no_lines = count_lines(100)
  do i=1,no_lines
    read(100,*) big_point
      
    ! write big_point to kpoint.dat for bs_quadratic.sh
    ! TODO: bug? this will be overwritten many times
    open(unit=101, file='kpoint.dat', action='write')
    write(101,*) big_point
    close(101)
    
    if (big_point/=1 .or. .not. acoustic) then
      do j=1,supercells(no_sc)%no_modes ! TODO: bug? (also lower down)
        ! --------------------
        ! generate amplitudes
        ! --------------------
        temp_string = '.'//trim(i2s(big_point))//'.'//trim(i2s(j))//'.dat'
        if (file_exists('energy'//trim(temp_string))) then
          allocate(energy(map%count))
          open(unit=101, file='energy'//trim(temp_string), status='old',&
            &action='read')
          do k=1,map%count
            read(101,*) energy(k)
          enddo
          close(101)
          
          open(unit=101, file='frequency'//trim(temp_string), status='old',&
            &action='read')
          read(101,*) frequency
          close(101)
          
          open(unit=101, file='size'//trim(temp_string), status='old',&
            &action='read')
          read(101,*) cell_size
          close(101)
          
          amplitudes = generate_amplitudes2(map,energy,frequency,cell_size)
          deallocate(energy)
          
          ! --------------------
          ! fit splines
          ! --------------------
          spline = quadratic_spline2(integration_points, amplitudes)
          deallocate(amplitudes)
          
          ! --------------------
          ! calculate 1-d anharmonic energy
          ! --------------------
          vscf = vscf_1d_2(frequency, spline(1,1), spline(2,:), Nbasis)
          
          deallocate(spline)
          deallocate(vscf%anh_pot)
          deallocate(vscf%hamiltonian)
          deallocate(vscf%eigenvals)
          deallocate(vscf%eigenvecs)
          
        endif
      enddo ! loop over modes
    endif
  enddo ! loop over list.dat
  
  ! read multiplicity and no_kpoints from ibz.dat
  open(100, file='ibz.dat', status='old', action='read')
  no_kpoints = count_lines(100)
  allocate(multiplicity(no_kpoints))
  do i=1,no_kpoints
    read(100,*) temp_real, temp_real, temp_real, multiplicity(no_kpoints)
  enddo
  close(100)
  
  ! call calculate_anharmonic
  call calculate_anharmonic2(multiplicity, supercells(no_sc)%no_modes, Nbasis)
  
  deallocate(multiplicity)
  deallocate(supercells)
end subroutine

end module
