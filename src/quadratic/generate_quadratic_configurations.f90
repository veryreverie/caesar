! Program to generate configurations along normal modes
module generate_quadratic_configurations_module
contains

subroutine generate_quadratic_configurations(args)
  use constants, only : dp,thermal
  use file_io,   only : open_read_file, open_write_file
  use string_module
  implicit none
  
  type(String), intent(in) :: args(:)
  
  ! file units
  integer :: super_eqm_file
  integer :: disp_patt_file
  integer :: positions_file
  
  ! filenames
  type(String) :: super_equilibrium_filename
  type(String) :: disp_patterns_filename
  type(String) :: positions_filename

  ! Input variables
  integer :: sampling_point,no_sampling_points,no_atoms
  real(dp) :: frequency,max_amplitude
  integer :: frequency_line ! the line in disp_patt.dat where freqency appears
  real(dp),allocatable :: atoms(:,:),mass(:),disp_patt(:,:)
  character(2),allocatable :: species(:)
  ! Working variables
  real(dp) :: amplitude,positions(3),quad_amplitude
  real(dp), parameter :: temperature=0.d0,tol=1.d-5
  real(dp), parameter :: thermal_energy=temperature/thermal
  
  ! temporary variables
  integer       :: i
  
  ! read input arguments
  max_amplitude = dble(args(1))
  sampling_point = int(args(2))
  no_sampling_points = int(args(3))
  frequency = dble(args(4))
  frequency_line = int(args(5))
  super_equilibrium_filename = args(6)
  disp_patterns_filename = args(7)
  positions_filename = args(8)
  
  ! Number of sampling points as read in is last_point-first_point, but we only
  ! want number of sampling points on each side of 0, as sampling_point goes
  ! from negative to positive
  no_sampling_points=no_sampling_points/2 
  
  ! take |frequency|
  frequency = dabs(frequency)
  
  ! Read in atoms
  super_eqm_file = open_read_file(super_equilibrium_filename)
  read(super_eqm_file,*) no_atoms
  allocate(atoms(no_atoms,3),mass(no_atoms),species(no_atoms))
  allocate(disp_patt(no_atoms,4))
  do i=1,no_atoms
    read(super_eqm_file,*) species(i),mass(i), atoms(i,:)
  enddo ! i
  close(super_eqm_file)
  
  ! Read in displacement pattern
  disp_patt_file = open_read_file(disp_patterns_filename)
  do i=1,frequency_line+2
    read(disp_patt_file,*) ! scroll forwards to the relevant lines
  enddo
  do i=1,no_atoms
    read(disp_patt_file,*) disp_patt(i,:)
  enddo ! i
  close(disp_patt_file)
  
  ! skip acoustic modes and equilibrium configuration
  if (frequency<=tol .or. sampling_point==0) stop
  
  ! calculate quad_amplitude
  ! The normal mode amplitudes to sample are calculated 
  IF(temperature<1.d-6)THEN
    ! This is really the normal mode amplitude squared,
    ! with units of 1/(E_h) as required.
    quad_amplitude = dsqrt(0.5d0/frequency)
  ELSE
    quad_amplitude = dsqrt( (1.d0/(dexp(frequency/thermal_energy)-1.d0)+0.5d0)&
                          & /frequency)
  END IF

  ! Calculate amplitude
  amplitude = max_amplitude                              &
          & * (sampling_point/(1.d0*no_sampling_points)) &
          & * quad_amplitude
  
  ! Write out displacement pattern
  positions_file = open_write_file(positions_filename)
  do i=1,no_atoms
    positions = atoms(i,:) + amplitude*disp_patt(i,1:3)*disp_patt(i,4:6)
    write(positions_file,*) species(i),mass(i),positions(:)
  enddo ! i
  close(positions_file)
end subroutine
end module
