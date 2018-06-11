! ======================================================================
! An array of type SampleResult.
! ======================================================================
! Exists to allow heterogeneous storage of sampling results.
module sample_results_module
  use common_module
  
  use sample_result_module
  implicit none
  
  private
  
  public :: SampleResults
  
  type, extends(NoDefaultConstructor) :: SampleResults
    type(SampleResult), allocatable :: results(:)
  end type
contains
end module
