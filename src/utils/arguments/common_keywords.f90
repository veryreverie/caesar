! ======================================================================
! Provides keywords common to all modes.
! ======================================================================
module common_keywords_submodule
  use random_module
  use io_module
  
  use keyword_submodule
  implicit none
  
  private
  
  public :: common_keywords
contains

function common_keywords() result(output)
  implicit none
  
  type(KeywordData), allocatable :: output(:)
  
  type(RandomReal) :: random_generator
  type(String)     :: seed
  
  ! Initialise a random number generator to get a default value of random_seed.
  random_generator = RandomReal()
  seed = str(random_generator%get_seed())
  
  output = [                                                                  &
  & KeywordData( 'interactive',                                               &
  &              'interactive specifies whether or not keywords can be &
  &specified interactively.',                                                 &
  &              default_value='true',                                        &
  &              allowed_in_file=.false.,                                     &
  &              can_be_interactive=.false.,                                  &
  &              flag_without_arguments='i'),                                 &
  & KeywordData( 'help',                                                      &
  &              'help requests helptext rather than running calculation.',   &
  &              is_optional=.true.,                                          &
  &              allowed_in_file=.false.,                                     &
  &              can_be_interactive=.false.,                                  &
  &              flag_with_arguments='h'),                                    &
  & KeywordData( 'input_file',                                                &
  &              'input_file specifies a file from which further settings &
  &will be read.',                                                            &
  &              is_optional=.true.,                                          &
  &              allowed_in_file=.false.,                                     &
  &              can_be_interactive=.false.,                                  &
  &              flag_with_arguments='false'),                                &
  & KeywordData( 'working_directory',                                         &
  &              'working_directory specifies the directory where all files &
  &and subsequent directories will be made.',                                 &
  &              default_value='.',                                           &
  &              is_path=.true.,                                              &
  &              allowed_in_file=.false.,                                     &
  &              can_be_interactive=.false.,                                  &
  &              flag_with_arguments='d'),                                    &
  & KeywordData( 'output_file',                                               &
  &              'output_file specifies a file to which terminal output will &
  &be written. This also disables terminal formatting, so should be favoured &
  &over piping to file. If unset, terminal output will go to the terminal.',  &
  &              is_optional=.true.,                                          &
  &              is_path=.true.,                                              &
  &              flag_with_arguments='o'),                                    &
  & KeywordData( 'random_seed',                                               &
  &              'random_seed specifies the seed which will be used to &
  &initialise any random number generation, allowing computations to be &
  &repeated exactly. If unset, the seed will be set to the current time in &
  &milliseconds.',                                                            &
  &              default_value=char(seed)) ]
end function
end module
