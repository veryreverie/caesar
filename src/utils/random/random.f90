!> Provides the [[RandomReal(type)]] random number generator.
module caesar_random_module
  use caesar_foundations_module
  use caesar_abstract_module
  use caesar_io_module
  implicit none
  
  private
  
  public :: RandomReal
  
  !> A random number generator.
  type, extends(NoDefaultConstructor) :: RandomReal
    logical, private :: initialised_ = .false.
    integer, private :: random_seed_
  contains
    procedure, public :: get_seed
    procedure, public :: random_number  => random_number_RandomReal
    procedure, public :: random_numbers => random_numbers_RandomReal
  end type

  interface RandomReal
    !> Constructor. If `seed` is given, this is used to seed the generator.
    !>    otherwise, the current time is used to seed the generator.
    module function new_RandomReal(seed) result(this) 
      integer, intent(in), optional :: seed
      type(RandomReal)              :: this
    end function
  end interface

  interface
    !> Return the seed used to initialise the generator.
    module function get_seed(this) result(output) 
      class(RandomReal), intent(in) :: this
      integer                       :: output
    end function
  end interface

  interface
    !> Return a random number, between `minimum` and `maximum`
    !>    which default to `0` and `1`.
    !> If one of `minimum` and `maximum` is specified, the other must also be
    !>    specified.
    module function random_number_RandomReal(this,minimum,maximum, &
       & log_distributed) result(output) 
      class(RandomReal), intent(in)           :: this
      real(dp),          intent(in), optional :: minimum
      real(dp),          intent(in), optional :: maximum
      !> Specifies that the random number should be drawn from a
      !>    log distribution, such that e.g. 0.1->1 is as likely as 0.01->0.1.
      logical,           intent(in), optional :: log_distributed
      real(dp)                                :: output
    end function
  end interface

  interface
    !> Returns an array of random numbers, of size `no_numbers`.
    module function random_numbers_RandomReal(this,no_numbers) result(output) 
      class(RandomReal), intent(in) :: this
      integer,           intent(in) :: no_numbers
      real(dp), allocatable         :: output(:)
    end function
  end interface
end module
