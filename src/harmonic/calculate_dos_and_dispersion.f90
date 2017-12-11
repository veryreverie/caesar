! ======================================================================
! Calculates:
!    - The density of states under the harmonic approximation.
!    - The phonon dispersion relation under the harmonic approximation.
! ======================================================================
! Should be run after calculate_harmonic.
! Not required for further calculations.
module calculate_dos_and_dispersion_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function calculate_dos_and_dispersion_keywords() result(keywords)
  use keyword_module
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [                                                                &
  & KeywordData( 'temperature',                                               &
  &              'temperature is the temperature used to broaden the phonon &
  &density of states. Temperature should be given in Kelvin.',                &
  &              default_value='0'),                                          &
  & KeywordData( 'path',                                                      &
  &              'path is the path through fractional reciprocal space which &
  &will be mapped by the phonon dispersion curve. The path should be &
  &specified as labels and q-points, separated by commas. The Gamma-point &
  &should be labelled G.', &
  &              default_value='G 0.0 0.0 0.0, R 0.5 0.5 0.5, M 0.0 0.5 0.5, &
  &G 0.0 0.0 0.0, X 0.0 0.0 0.5') ]
end function

function calculate_dos_and_dispersion_mode() result(output)
  use caesar_modes_module
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'calculate_dos_and_dispersion'
  output%description = 'Calculates the density of states and phonon &
     &dispersion under the harmonic approximation. Should be run after &
     &calculate_harmonic.'
  output%keywords = calculate_dos_and_dispersion_keywords()
  output%main_subroutine => calculate_dos_and_dispersion
end function

! ----------------------------------------------------------------------
! The main subroutine.
! ----------------------------------------------------------------------
subroutine calculate_dos_and_dispersion(arguments)
  use dictionary_module
  use structure_module
  use qpoints_module
  use group_module
  use dynamical_matrix_module
  use lte_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Inputs.
  type(String)                  :: wd
  real(dp)                      :: temperature
  type(String),     allocatable :: path_string(:)
  type(String),     allocatable :: path_labels(:)
  type(RealVector), allocatable :: path_qpoints(:)
  
  ! Previously calculated data.
  type(StructureData)                :: structure
  type(StructureData)                :: large_supercell
  type(QpointData),      allocatable :: qpoints(:)
  type(DynamicalMatrix), allocatable :: dynamical_matrices(:)
  
  ! Temporary variables.
  type(String), allocatable :: path_point(:)
  integer                   :: i,ialloc
  
  ! --------------------------------------------------
  ! Read in arguments from user.
  ! --------------------------------------------------
  wd = arguments%value('working_directory')
  temperature = dble(arguments%value('temperature'))
  path_string = split(arguments%value('path'), ',')
  
  ! --------------------------------------------------
  ! Generate path for dispersion calculation.
  ! --------------------------------------------------
  allocate( path_qpoints(size(path_string)), &
          & path_labels(size(path_string)),  &
          & stat=ialloc); call err(ialloc)
  do i=1,size(path_string)
    path_point = split(path_string(i))
    path_labels(i) = path_point(1)
    path_qpoints(i) = dble(path_point(2:4))
  enddo
  
  ! --------------------------------------------------
  ! Read in previously calculated data.
  ! --------------------------------------------------
  structure = read_structure_file(wd//'/structure.dat')
  
  large_supercell = read_structure_file(wd//'/large_supercell.dat')
  
  qpoints = read_qpoints_file(wd//'/qpoints.dat')
  
  allocate( dynamical_matrices(size(qpoints)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(qpoints)
    dynamical_matrices(i) = read_dynamical_matrix_file(  &
       & wd//'/qpoint_'//left_pad(i,str(size(qpoints)))// &
       & '/dynamical_matrix.dat')
  enddo
  
  ! --------------------------------------------------
  ! Calculate harmonic DOS and phonon dispersion.
  ! --------------------------------------------------
  call fourier_interpolation(              &
     & qpoints,                            &
     & dynamical_matrices,                 &
     & structure,                          &
     & temperature,                        &
     & large_supercell,                    &
     & path_labels,                        &
     & path_qpoints,                       &
     & wd//'/phonon_dispersion_curve.dat', &
     & wd//'/high_symmetry_points.dat',    &
     & wd//'/free_energy.dat',             &
     & wd//'/freq_dos.dat')
end subroutine
end module