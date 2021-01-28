! ======================================================================
! The harmonic approximation to the stress, in cartesian co-ordinates.
! ======================================================================
module caesar_stress_hessian_module
  use caesar_utils_module
  
  use caesar_structure_module
  use caesar_normal_mode_module
  implicit none
  
  private
  
  public :: StressHessian
  
  type, extends(Stringsable) :: StressHessian
    type(CartesianHessian), allocatable :: elements(:,:)
  contains
    ! I/O.
    procedure, public :: read  => read_StressHessian
    procedure, public :: write => write_StressHessian
  end type
  
  interface StressHessian
    module procedure new_StressHessian
    module procedure new_StressHessian_Strings
    module procedure new_StressHessian_StringArray
  end interface
contains

function new_StressHessian(elements) result(this)
  implicit none
  
  type(CartesianHessian), intent(in) :: elements(:,:)
  type(StressHessian)                :: this
  
  this%elements = elements
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_StressHessian(this,input)
  implicit none
  
  class(StressHessian), intent(out) :: this
  type(String),         intent(in)  :: input(:)
  
  type(StringArray), allocatable :: sections(:)
  
  type(CartesianHessian), allocatable :: elements(:,:)
  
  integer :: i
  
  select type(this); type is(StressHessian)
    sections = split_into_sections( input,                         &
                                  & separating_line=repeat('-',50) )
    do i=1,size(sections)
      sections(i)%strings = sections(i)%strings(2:)
    enddo
    elements = transpose(reshape(CartesianHessian(sections), [3,3]))
    this = StressHessian(elements)
  class default
    call err()
  end select
end subroutine

function write_StressHessian(this) result(output)
  implicit none
  
  class(StressHessian), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  character(1), parameter :: directions(3) = ['x', 'y', 'z']
  
  integer :: i,j
  
  select type(this); type is(StressHessian)
    output = [str(repeat('-',50))]
    do i=1,3
      do j=1,3
        output = [                                                        &
           & output,                                                      &
           & str('Stress component '//directions(i)//directions(j)//':'), &
           & str(this%elements(i,j)),                                     &
           & str(repeat('-',50))                                          ]
      enddo
    enddo
  class default
    call err()
  end select
end function

function new_StressHessian_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(StressHessian)      :: this
  
  call this%read(input)
end function

impure elemental function new_StressHessian_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(StressHessian)           :: this
  
  this = StressHessian(str(input))
end function
end module
