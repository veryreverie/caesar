! ======================================================================
! Program to calculate the positions of the atoms along a normal mode
! Author: B. Monserrat
! Created: 30 November 2011
!
! Modifications:
! 30/11/2011 - Included estimation of normal mode amplitude at temperature T
! 09/10/2015 - Modified by JCAP to work with QUIP+GAP. equilibrium.dat is now in angstroms, not bohr, and in extended xyz format
! ======================================================================
module positions_module
  use string_module
  use io_module
  use constants_module, only : dp
contains

! ----------------------------------------------------------------------
! The normal mode amplitudes to sample are calculated 
! ----------------------------------------------------------------------
function max_amplitude(frequency,temperature) result(output)
  use constants_module, only : kb_in_au
  implicit none
  
  real(dp), intent(in) :: frequency
  real(dp), intent(in) :: temperature
  real(dp)             :: output
  
  real(dp),parameter :: tolerance=1.0e-6_dp
  real(dp) :: normal_mode_amplitude, thermal_energy
  
  if(temperature<tolerance)then
    ! n_m_a is really the normal mode amplitude squared,
    !    with units of 1/(e_h) as required.
    ! The maximum amplitude has units of sqrt(1/e_h) as required
    normal_mode_amplitude = 0.5_dp/frequency
    output = 5*sqrt(normal_mode_amplitude)
  else
    thermal_energy = 0.5_dp*temperature*kb_in_au
    normal_mode_amplitude = (1/(exp(frequency/thermal_energy)-1)+0.5_dp) &
                        & / frequency
    output = 5*sqrt(normal_mode_amplitude)
    call print_line('maximum normal mode amplitude is (a.u.) '//output)
  end if
end function max_amplitude

subroutine positions(structure)
  use structure_module
  use linear_algebra_module
  implicit none
  
  type(StructureData), intent(in) :: structure

  real(dp), parameter :: tol=1.0e-6_dp
  real(dp) :: scaling1, scaling2
  real(dp) :: frequency1, frequency2, temperature
  integer :: number_fit_points, current_number1, current_number2
  real(dp) :: amplitude
  
  type(RealVector) :: disp1
  type(RealVector) :: disp2
  type(RealVector) :: dpos1
  type(RealVector) :: dpos2
  type(RealVector) :: pos
  
  type(String), allocatable :: disp1_file(:)
  type(String), allocatable :: disp2_file(:)
  
  integer :: positions_file
  
  integer :: i
  type(String), allocatable :: line(:)
  
  disp1_file = read_lines('disp1.dat')
  line = split(disp1_file(1))
  frequency1 = dble(line(1))
  temperature = dble(line(2))
  number_fit_points = int(line(3))
  current_number1 = int(line(4))
  
  disp2_file = read_lines('disp2.dat')
  line = split(disp1_file(1))
  frequency2 = dble(line(1))
  temperature = dble(line(2))
  number_fit_points = int(line(3))
  current_number2 = int(line(4))
  
  positions_file = open_write_file('positions.dat')
  
  if (abs(frequency1)<tol .or. abs(frequency2)<tol) then
    scaling1 = 0
    scaling2 = 0
  else
    ! Calculate points on each side of zero.
    amplitude = (number_fit_points-1)/2.0_dp
    
    ! Calculate current scaling.
    scaling1 = -1 + (current_number1-1)/amplitude
    scaling2 = -1 + (current_number2-1)/amplitude
    
    scaling1 = scaling1 * max_amplitude(abs(frequency1),temperature)
    scaling2 = scaling2 * max_amplitude(abs(frequency2),temperature)
  endif
  call print_line(positions_file, scaling1//' '//scaling2//' '//1)

  call print_line('amplitude of first mode is '//scaling1)
  call print_line('amplitude of second mode is '//scaling2)

  do i=1,structure%no_atoms
    disp1 = dble(split(disp1_file(i+1)))
    disp2 = dble(split(disp2_file(i+1)))
    
    dpos1 = scaling1*disp1
    dpos2 = scaling2*disp2
    
    pos = structure%atoms(i) + dpos1 + dpos2
    
    call print_ine(positions_file, structure%species(i)//' '//pos)
  end do
  
  close(positions_file)
end subroutine
end module
