! Program to generate configurations along normal modes

module generate_quadratic_configurations_module
  use constants, only : dp
  implicit none
contains

subroutine generate_quadratic_configurations()
  implicit none
  ! Input variables
  integer :: sampling_point,no_sampling_points,no_atoms
  real(dp) :: frequency,max_amplitude
  real(dp),allocatable :: atoms(:,:),mass(:),disp_patt(:,:)
  character(2),allocatable :: species(:)
  ! Working variables
  integer :: i
  real(dp) :: amplitude,positions(3),quad_amplitude
  character(80) :: ch_dump
  real(dp) :: temperature=0.d0,tol=1.d-5
  real(dp) :: thermal_energy

  ! Read in configuration
  open(1,file='configuration.dat')
  read(1,*)max_amplitude,sampling_point,no_sampling_points
  close(1)
  ! Number of samplings points as read in is last_point-first_point, but we only
  ! want number of sampling points on each side of 0, as sampling_point goes
  ! from negative to positive
  no_sampling_points=no_sampling_points/2 
  ! Read in atoms
  open(1,file='super_equilibrium.dat')
  read(1,*)no_atoms
  allocate(atoms(no_atoms,3),mass(no_atoms),species(no_atoms))
  allocate(disp_patt(no_atoms,4))
  do i=1,no_atoms
    read(1,*)species(i),mass(i),atoms(i,:)
  enddo ! i
  close(1)
  ! Read in displacement pattern
  open(1,file='disp_patterns_temp.dat')
  read(1,*)ch_dump,ch_dump,frequency
  read(1,*)
  read(1,*)
  do i=1,no_atoms
    read(1,*)disp_patt(i,:)
  enddo ! i
  close(1)
  open(1,file='frequency.dat')
  write(1,*)frequency*27.21138602d0
  close(1)
  
  ! calculate quad_amplitude
  ! The normal mode amplitudes to sample are calculated 
  IF(temperature<1.d-6)THEN
    quad_amplitude=sqrt(0.5/abs(frequency)) ! This is really the normal mode amplitude squared, with units of 1/(E_h) as required.
  ELSE
    thermal_energy=temperature/3.1577464E5 ! thermal energy in a.u.
    quad_amplitude=sqrt((1/(EXP(abs(frequency)/thermal_energy)-1)+0.5)/abs(frequency))
  END IF

  ! Calculate amplitude
  amplitude=max_amplitude*(sampling_point/(1.d0*no_sampling_points))*quad_amplitude
  
  ! Write out displacement pattern
  if(abs(frequency)>tol)then
    if(sampling_point/=0)then
      open(1,file='positions.dat')
      do i=1,no_atoms
        positions(1)=atoms(i,1)+amplitude*disp_patt(i,1)*disp_patt(i,4)
        positions(2)=atoms(i,2)+amplitude*disp_patt(i,2)*disp_patt(i,4)
        positions(3)=atoms(i,3)+amplitude*disp_patt(i,3)*disp_patt(i,4)
        write(1,*)species(i),mass(i),positions(:)
      enddo ! i
      close(1)
    endif ! skip equilibrium configuration
  endif ! skip acoustic modes
end subroutine
end module
