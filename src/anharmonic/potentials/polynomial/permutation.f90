! ======================================================================
! Generates and stores the permutation of one list against another.
! ======================================================================
! Takes two lists of integers, a and b, which may or may not contain
!    duplicate elements.
! Generates all unique permutations of a against b.
! See permutation_example for an example of how to use this class.
module caesar_permutation_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: PermutationData
  public :: permutation_example
  
  type, extends(NoDefaultConstructor) :: PermutationData
    integer, allocatable, private :: a_(:)
    integer, allocatable, private :: b_(:)
    integer, allocatable, private :: bins_(:,:)
    integer, allocatable, private :: a_permutation_(:)
    integer, allocatable, private :: b_permutation_(:)
    logical,              private :: all_permutations_done_
  contains
    procedure, public :: next_permutation => &
                       & next_permutation_PermutationData
    procedure, public :: a => &
                       & a_PermutationData
    procedure, public :: b => &
                       & b_PermutationData
    procedure, public :: log_no_permutations => &
                       & log_no_permutations_PermutationData
    procedure, public :: all_permutations_done => &
                       & all_permutations_done_PermutationData
  end type
  
  interface PermutationData
    module function new_PermutationData(a,b) result(this) 
      integer, intent(in)   :: a(:)
      integer, intent(in)   :: b(:)
      type(PermutationData) :: this
    end function
  end interface
  
  interface
    ! Generates the next permutation of a against b.
    ! Keeps each bin of a sorted in ascending order.
    ! Does not permute b.
    module subroutine next_permutation_PermutationData(this) 
      class(PermutationData), intent(inout) :: this
    end subroutine
  end interface
  
  interface
    module function a_PermutationData(this) result(output) 
      class(PermutationData), intent(in) :: this
      integer, allocatable               :: output(:)
    end function
  end interface
  
  interface
    module function b_PermutationData(this) result(output) 
      class(PermutationData), intent(in) :: this
      integer, allocatable               :: output(:)
    end function
  end interface
  
  interface
    module function log_no_permutations_PermutationData(this) result(output) 
      class(PermutationData), intent(in) :: this
      real(dp)                           :: output
    end function
  end interface
  
  interface
    module function all_permutations_done_PermutationData(this) result(output) 
      class(PermutationData), intent(in) :: this
      logical                            :: output
    end function
  end interface
  
  interface
    ! --------------------------------------------------
    ! Usage example.
    ! --------------------------------------------------
    module subroutine permutation_example() 
    end subroutine
  end interface
end module
