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
  use keyword_module
  implicit none
  
  type(KeywordData) :: keywords(6)
  
  keywords = [                                                                &
  & make_keyword( 'harmonic_states_cutoff',                                   &
  &               'harmonic_states_cutoff is the number of harmonic &
  &eigenstates in the direction of each normal mode.'),                       &
  & make_keyword( 'potential_basis_cutoff',                                &
  &               'potential_basis_cutoff is the order up to which the &
  &potential is expanded. e.g. a cubic expansion would be order 3.'),         &
  & make_keyword( 'scf_convergence_threshold',                                &
  &               'scf_convergence_threshold is the energy to within which &
  &the VSCF calculation will be converged.'),                                 &
  & make_keyword( 'max_scf_cycles',                                           &
  &               'max_scf_cycles is the maximum number of SCF cycles which &
  &will be carried out as part of the VSCF calculation.'),                    &
  & make_keyword( 'perturbation_order',                                       &
  &               'perturbation_order is the order up to which perturbation &
  &theory will be run',                                                       &
  &               is_optional=.true.),                                        &
  & make_keyword( 'perturb_states_to_same_order',                             &
  &               'perturb_states_to_same_order specifies whether or not to &
  &calculate state correction at the same order as energy corrections.',      &
  &               default_value='f') ]
end function

subroutine test(arguments)
  use dictionary_module
  use linear_algebra_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  call print_line('END TEST')
end subroutine
end module
