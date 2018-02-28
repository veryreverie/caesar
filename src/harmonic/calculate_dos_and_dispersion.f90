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
  &should be labelled G.',                                                    &
  &              default_value='G 0.0 0.0 0.0, R 0.5 0.5 0.5, M 0.0 0.5 0.5, &
  &G 0.0 0.0 0.0, X 0.0 0.0 0.5'),                                            &
  & KeywordData( 'no_dos_samples',                                            &
  &              'no_dos_samples is the number of points in reciprocal space &
  &at which the normal modes are calculated when calculating the vibrational &
  &density of states.',                                                       &
  &              default_value='100000') ]
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
  use force_constants_module
  use dynamical_matrix_module
  use min_images_module
  use lte_module
  use ofile_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Inputs.
  type(String)                  :: wd
  real(dp)                      :: temperature
  type(String),     allocatable :: path_string(:)
  type(String),     allocatable :: path_labels(:)
  type(RealVector), allocatable :: path_qpoints(:)
  integer                       :: no_dos_samples
  
  ! Previously calculated data.
  type(StructureData)                :: structure
  type(StructureData)                :: large_supercell
  type(QpointData),      allocatable :: qpoints(:)
  type(DynamicalMatrix), allocatable :: dynamical_matrices(:)
  
  ! Working variables.
  type(ForceConstants)          :: force_constants
  type(MinImages),  allocatable :: min_images(:,:)
  
  ! Logfile.
  type(OFile) :: logfile
  type(DynamicalMatrix) :: dyn_mat
  
  ! Temporary variables.
  type(String), allocatable :: path_point(:)
  integer                   :: i,ialloc
  
  ! --------------------------------------------------
  ! Read in arguments from user.
  ! --------------------------------------------------
  wd = arguments%value('working_directory')
  temperature = dble(arguments%value('temperature'))
  path_string = split(arguments%value('path'), ',')
  no_dos_samples = int(arguments%value('no_dos_samples'))
  
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
  ! Run calculations.
  ! --------------------------------------------------
  
  logfile = wd//'/dos_and_dispersion_log.dat'
  
  ! Construct the matrix of force constants from dynamical matrices.
  force_constants = reconstruct_force_constants( large_supercell,    &
                                               & qpoints,            &
                                               & dynamical_matrices, &
                                               & logfile)
  
  ! Calculate minimum image distances.
  min_images = calculate_min_images(large_supercell)
  
  ! Reconstruct calculated dynamical matrices, to check for consistency.
  do i=1,size(qpoints)
    dyn_mat = DynamicalMatrix( dblevec(qpoints(i)%qpoint), &
                             & large_supercell,            &
                             & force_constants,            &
                             & min_images)
    call logfile%print_line('Comparing dynamical matrices before and after &
       &reconstruction of force constants.')
    call compare_dynamical_matrices(dynamical_matrices(i),dyn_mat,logfile)
  enddo
  
  ! Generate harmonic phonon dispersion curve by interpolating between
  !    calculated q-points using Fourier interpolation.
  call generate_dispersion( large_supercell,                    &
                          & min_images,                         &
                          & force_constants,                    &
                          & path_labels,                        &
                          & path_qpoints,                       &
                          & wd//'/phonon_dispersion_curve.dat', &
                          & wd//'/high_symmetry_points.dat',    &
                          & logfile)
  
  ! Generate harmonic phonon density of states, interpolating as above.
  call generate_dos( large_supercell,        &
                   & min_images,             &
                   & force_constants,        &
                   & temperature,            &
                   & no_dos_samples,         &
                   & wd//'/free_energy.dat', &
                   & wd//'/freq_dos.dat',    &
                   & logfile)
end subroutine
end module
