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
  use linear_algebra_module
  
  type :: PositionData
    logical                       :: scaling_check
    real(dp)                      :: scaling1
    real(dp)                      :: scaling2
    type(RealVector), allocatable :: positions(:)
  end type
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

function positions(structure, temperature, number_fit_points, frequency1, &
   &frequency2, current_number1, current_number2, disp1, disp2) result(output)
  use structure_module
  use linear_algebra_module
  implicit none
  
  type(StructureData), intent(in) :: structure
  real(dp),            intent(in) :: temperature
  integer,             intent(in) :: number_fit_points
  real(dp),            intent(in) :: frequency1
  real(dp),            intent(in) :: frequency2
  integer,             intent(in) :: current_number1
  integer,             intent(in) :: current_number2
  type(RealVector),    intent(in) :: disp1(:)
  type(RealVector),    intent(in) :: disp2(:)
  type(PositionData)              :: output

  real(dp), parameter :: tol=1.0e-6_dp
  
  integer :: i,ialloc
  
  if (abs(frequency1)<tol .or. abs(frequency2)<tol) then
    output%scaling_check = .false.
    output%scaling1 = 0
    output%scaling2 = 0
    allocate( output%positions(size(structure%atoms)), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(structure%atoms)
      output%positions(i) = structure%atoms(i)%cartesian_position()
    enddo
  else
    output%scaling_check = .true.
    output%scaling1 = (-1 + (current_number1-1)*2.0_dp/(number_fit_points-1)) &
                  & * max_amplitude(abs(frequency1),temperature)
    output%scaling2 = (-1 + (current_number2-1)*2.0_dp/(number_fit_points-1)) &
                  & * max_amplitude(abs(frequency2),temperature)
    allocate( output%positions(size(structure%atoms)), &
            & stat=ialloc); call err(ialloc)
    do i=1,structure%no_atoms
      output%positions(i) = structure%atoms(i)%cartesian_position() &
                        & + output%scaling1*disp1(i)                &
                        & + output%scaling2*disp2(i)
    enddo
  endif
end function
end module
