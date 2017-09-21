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
  & make_keyword( 'temperature',                                              &
  &               'temperature is the temperature in Kelvin, used when &
  &calculating the density of states and phonon dispersion curve.',           &
  &               default_value='0'),                                         &
  & make_keyword( 'path',                                                     &
  &               'path is the path through fractional reciprocal space which &
  &will be mapped by the phonon dispersion curve. The path should be &
  &specified as vectors in reciprocal space separated by commas.',            &
  &               default_value='0.0 0.0 0.0, 0.5 0.5 0.5, 0.0 0.5 0.5, &
  &0.0 0.0 0.0, 0.0 0.5 0.0') ]
end function

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
  type(RealVector), allocatable :: disp_qpoints(:)
  
  ! Previously calculated data.
  type(StructureData)                :: structure
  type(StructureData)                :: large_supercell
  type(QpointData),      allocatable :: qpoints_ibz(:)
  type(Group),           allocatable :: atom_symmetry_group(:)
  type(DynamicalMatrix), allocatable :: ibz_dynamical_matrices(:)
  
  ! Temporary variables.
  integer :: i,ialloc
  
  ! --------------------------------------------------
  ! Read in arguments from user.
  ! --------------------------------------------------
  wd = arguments%value('working_directory')
  temperature = dble(arguments%value('temperature'))
  path_string = split(arguments%value('path'), ',')
  
  ! --------------------------------------------------
  ! Generate path for dispersion calculation.
  ! --------------------------------------------------
  allocate(disp_qpoints(size(path_string)), stat=ialloc); call err(ialloc)
  do i=1,size(path_string)
    disp_qpoints(i) = dble(split(path_string(i)))
  enddo
  
  ! --------------------------------------------------
  ! Read in previously calculated data.
  ! --------------------------------------------------
  structure = read_structure_file(wd//'/structure.dat')
  
  large_supercell = read_structure_file(wd//'/large_supercell.dat')
  
  qpoints_ibz = read_qpoints_file(wd//'/qpoints_ibz.dat')
  
  atom_symmetry_group = read_group_file( &
     & wd//'/Supercell_1/atom_symmetry_group.dat')
  
  allocate( ibz_dynamical_matrices(size(qpoints_ibz)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(qpoints_ibz)
    ibz_dynamical_matrices(i) = read_dynamical_matrix_file( &
       & wd//'/qpoint_'//i//'/dynamical_matrix.dat')
  enddo
  
  ! --------------------------------------------------
  ! Calculate harmonic DOS and phonon dispersion.
  ! --------------------------------------------------
  call fourier_interpolation(              &
     & ibz_dynamical_matrices,             &
     & structure,                          &
     & temperature,                        &
     & large_supercell,                    &
     & qpoints_ibz,                        &
     & disp_qpoints,                       &
     & atom_symmetry_group,                &
     & wd//'/phonon_dispersion_curve.dat', &
     & wd//'/high_symmetry_points.dat',    &
     & wd//'/free_energy.dat',             &
     & wd//'/freq_dos.dat')
end subroutine
end module
