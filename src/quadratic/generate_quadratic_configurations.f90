! Program to generate configurations along normal modes
module generate_quadratic_configurations_module
contains

subroutine generate_quadratic_configurations(args)
  use constants, only : dp,thermal
  use file_module
  use string_module
  use structure_module
  use displacement_patterns_module
  implicit none
  
  type(String), intent(in) :: args(:)
  
  ! Parameters
  real(dp), parameter :: temperature=0.d0
  real(dp), parameter :: tol=1.d-5
  real(dp), parameter :: thermal_energy=temperature/thermal
  
  ! filenames
  type(String) :: structure_filename
  type(String) :: structure_sc_filename
  type(String) :: disp_patterns_filename
  type(String) :: mode_structure_filename

  ! Input variables
  real(dp) :: frequency,max_amplitude
  real(dp),allocatable :: disp_patt(:,:)
  
  integer             :: gvector
  integer             :: mode
  integer             :: sampling_point
  integer             :: no_sampling_points
  
  ! Structure data
  type(StructureData) :: structure
  type(StructureData) :: structure_sc
  
  ! Displacement patterns data
  type(DispPatterns) :: disp_patterns
  
  ! Working variables
  real(dp)       :: amplitude
  real(dp)       :: quad_amplitude
  
  ! temporary variables
  integer :: i
  
  ! read input arguments
  max_amplitude = dble(args(1))
  gvector = int(args(2))
  mode = int(args(3))
  sampling_point = int(args(4))
  no_sampling_points = int(args(5))
  structure_filename = args(9)
  structure_sc_filename = args(10)
  disp_patterns_filename = args(11)
  mode_structure_filename = args(12)
  
  ! Number of sampling points as read in is last_point-first_point, but we only
  ! want number of sampling points on each side of 0, as sampling_point goes
  ! from negative to positive
  no_sampling_points=no_sampling_points/2 
  
  ! Read in input structures
  structure = read_structure_file(structure_filename)
  structure_sc = read_structure_file(structure_sc_filename)
  
  ! Allocate disp_patt
  allocate(disp_patt(6,structure_sc%no_atoms))
  
  disp_patterns = read_disp_patterns_file( disp_patterns_filename, &
                                         & structure_sc%no_modes)
  frequency = disp_patterns%frequencies(mode,gvector)
  disp_patt = disp_patterns%disp_patterns(:,:,mode,gvector)
  
  ! skip acoustic modes and equilibrium configuration
  if (frequency<=tol .or. sampling_point==0) stop
  
  ! calculate quad_amplitude
  ! The normal mode amplitudes to sample are calculated 
  if(temperature<1.d-6)then
    ! This is really the normal mode amplitude squared,
    ! with units of 1/(E_h) as required.
    quad_amplitude = dsqrt(0.5d0/frequency)
  else
    quad_amplitude = dsqrt( (1.d0/(dexp(frequency/thermal_energy)-1.d0)+0.5d0)&
                          & /frequency)
  end if

  ! Calculate amplitude
  amplitude = max_amplitude                              &
          & * (sampling_point/(1.d0*no_sampling_points)) &
          & * quad_amplitude
  
  ! Calculate new positions
  do i=1,structure_sc%no_atoms
    structure_sc%atoms(:,i) = structure_sc%atoms(:,i) &
                          & + amplitude*disp_patt(1:3,i)*disp_patt(4:6,i)
  enddo
  
  call write_structure_file(structure_sc,mode_structure_filename)
end subroutine
end module
