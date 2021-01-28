! ======================================================================
! An array of type SampleResult.
! ======================================================================
! Exists to allow heterogeneous storage of sampling results.
module caesar_sample_results_module
  use caesar_common_module
  
  use caesar_sample_result_module
  implicit none
  
  private
  
  public :: SampleResults
  
  type, extends(NoDefaultConstructor) :: SampleResults
    type(SampleResult), allocatable :: results(:)
  end type
contains
end module
