! Program to generate configurations along normal modes

module generate_quadratic_configurations_module
contains

subroutine generate_quadratic_configurations(args)
  use constants, only : dp,thermal
  use file_io,   only : open_read_file, open_write_file
  implicit none
  
  character(*), intent(in) :: args(:)
  
  ! file units
  integer :: super_eqm_file     ! super_equilibrium.dat
  integer :: disp_patt_file     ! disp_patterns_temp.dat
  integer :: positions_file     ! positions.dat

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
  
  ! read non-filename input arguments
  read(args(1),*) max_amplitude
  read(args(2),*) sampling_point
  read(args(3),*) no_sampling_points
  read(args(4),*) frequency
  read(args(5),*) frequency_line
  
  ! Number of sampling points as read in is last_point-first_point, but we only
  ! want number of sampling points on each side of 0, as sampling_point goes
  ! from negative to positive
  no_sampling_points=no_sampling_points/2 
  
  ! take |frequency|
  frequency = dabs(frequency)
  
  ! Read in atoms
  super_eqm_file = open_read_file(args(6))
  read(super_eqm_file,*) no_atoms
  allocate(atoms(no_atoms,3),mass(no_atoms),species(no_atoms))
  allocate(disp_patt(no_atoms,4))
  do i=1,no_atoms
    read(super_eqm_file,*) species(i),mass(i), atoms(i,:)
  enddo ! i
  close(super_eqm_file)
  
  ! Read in displacement pattern
  disp_patt_file = open_read_file(args(7))
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
  positions_file = open_write_file(args(8))
  do i=1,no_atoms
    positions(1)=atoms(i,1)+amplitude*disp_patt(i,1)*disp_patt(i,4)
    positions(2)=atoms(i,2)+amplitude*disp_patt(i,2)*disp_patt(i,4)
    positions(3)=atoms(i,3)+amplitude*disp_patt(i,3)*disp_patt(i,4)
    write(positions_file,*) species(i),mass(i),positions(:)
  enddo ! i
  close(positions_file)
end subroutine
end module
