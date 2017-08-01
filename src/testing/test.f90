! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module test_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function test_keywords() result(keywords)
  use help_module
  implicit none
  
  type(KeywordData) :: keywords(15)
  
  keywords = [                                                                &
  & make_keyword( 'seed_name',                                                &
  &               'seed_name is the DFT seedname from which file names are &
  &constructed.'),                                                            &
  & make_keyword( 'num_indep_data',                                           &
  &               'num_indep_data is the number of data points per mode. It &
  &should be an odd integer.',                                                &
  &               default_value='11'),                                        &
  & make_keyword( 'temperature', &
  &               'temperature is the temperature in Kelvin at which &
  &thermodynamic quantities are calculated.'),                                &
  & make_keyword( 'first_mode',                                               &
  &               'first_mode is the first mode to be considered.',           &
  &               default_value='4'),                                         &
  & make_keyword( 'last_mode',                                                &
  &               'last_mode is the last mode to be considered.',             &
  &               default_value='4'),                                         &
  & make_keyword( 'first_amp',                                                &
  &               'first_amp is the first amplitude to be considered.',       &
  &               default_value='1'),                                         &
  & make_keyword( 'last_amp',                                                 &
  &               'last_amp is the last amplitude to be considered.',         &
  &               default_keyword='num_indep_data'),                          &
  & make_keyword( 'mc_sampling',                                              &
  &               'mc_sampling should be specified to turn on Monte-Carlo &
  &sampling.',                                                                &
  &               is_boolean=.true.),                                         &
  & make_keyword( 'mc_data_points',                                           &
  &               'mc_data_points is the number of Monte-Carlo data points. &
  &Should only be set if mc_sampling is set.',                                &
  &               default_value='20'),                                        &
  & make_keyword( 'mc_continuation',                                          &
  &               'mc_continuation is the Monte-Carlo continuation. Should &
  &only be set if mc_sampling is set.',                                       &
  &               default_value='20'),                                        &
  & make_keyword( 'coupled_sampling',                                         &
  &               'coupled_sampling should be specified to turn on coupled &
  &sampling.',                                                                &
  &                is_boolean=.true.),                                        &
  & make_keyword( 'num_2body_data',                                           &
  &               'num_2body_data is the number of two-body data points. &
  &Should only be set if coupled_sampling is set.',                           &
  &               default_value='11'),                                        &
  & make_keyword( 'first_amp_2body',                                          &
  &               'first_amp_2body is the first two-body amplitude &
  &considered. Must be <= num_2body_data.',                                   &
  &               default_value='1'),                                         &
  & make_keyword( 'last_amp_2body',                                           &
  &               'last_amp_2body is the last two-body amplitude considered. &
  &Must be >= first_amp_2body and <= num_2body data.',                        &
  &               default_keyword='num_2body_data'),                          &
  & make_keyword( 'magres',                                                   &
  &               'magres specifies whether or not the DFT calculation is &
  &magres.',                                                                  &
  &               is_boolean=.true.) ]
end function

subroutine test(arguments)
  use dictionary_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: wd
  type(String) :: filename
  
  complex(dp)               :: a
  complex(dp)               :: b
  integer                   :: output_file
  type(String), allocatable :: input_file(:)
  
  wd = arguments%value('working_directory')
  filename = wd//'/test.test'
  
  output_file = open_write_file(filename)
  
  a = cmplx( 1.0, 2.0,dp)
  call print_line(a)
  call print_line(output_file, a)
  a = cmplx(-1.0, 2.0,dp)
  call print_line(a)
  call print_line(output_file, a)
  a = cmplx( 1.0,-2.0,dp)
  call print_line(a)
  call print_line(output_file, a)
  a = cmplx(-1.0,-2.0,dp)
  call print_line(a)
  call print_line(output_file, a)
  
  close(output_file)
  
  input_file = read_lines(filename)
  b = cmplx(input_file(1))
  call print_line(b)
  b = cmplx(input_file(2))
  call print_line(b)
  b = cmplx(input_file(3))
  call print_line(b)
  b = cmplx(input_file(4))
  call print_line(b)
  
end subroutine
end module
