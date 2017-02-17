module algebra_module
  use utils, only : dp
  use string_module
  implicit none
  
  type Vector3
    real(dp), private :: contents(3)
  end type
  
  type Matrix33
    real(dp), private :: contents(3,3)
  end type
  
  interface assignment(=)
    module procedure assign_Vector3_reals
    module procedure assign_Matrix33_reals_1d
    module procedure assign_Matrix33_reals_2d
  end interface
  
  interface str
    module procedure str_Vector3
  end interface
  
  interface strs
    module procedure strs_Matrix33
  end interface
  
  interface vec3
    module procedure vec3_String
  end interface
  
  interface mat33
    module procedure mat33_Strings
  end interface

contains

pure subroutine assign_Vector3_reals(output,input)
  implicit none
  
  real(dp),      intent(in)  :: input(3)
  type(Vector3), intent(out) :: output
  
  output%contents = input
end subroutine

pure subroutine assign_Matrix33_reals_1d(output,input)
  implicit none
  
  real(dp),       intent(in)  :: input(9)
  type(Matrix33), intent(out) :: output
  
  output%contents(1,:) = input(1:3)
  output%contents(2,:) = input(4:6)
  output%contents(3,:) = input(7:9)
end subroutine

pure subroutine assign_Matrix33_reals_2d(output,input)
  implicit none
  
  real(dp),       intent(in)  :: input(3,3)
  type(Matrix33), intent(out) :: output
  
  output%contents = input
end subroutine

pure function str_Vector3(this) result(output)
  implicit none
  
  type(Vector3), intent(in) :: this
  type(String)              :: output
  
  output = str(this%contents(1))//' '//this%contents(2)//' '//this%contents(3)
end function

pure function strs_Matrix33(this) result(output)
  implicit none
  
  type(Matrix33), intent(in) :: this
  type(String)               :: output(3)
  
  output(1) = str(this%contents(1,1))//' '// &
            &     this%contents(1,2) //' '// &
            &     this%contents(1,3)
  output(2) = str(this%contents(2,1))//' '// &
            &     this%contents(2,2) //' '// &
            &     this%contents(2,3)
  output(3) = str(this%contents(3,1))//' '// &
            &     this%contents(3,2) //' '// &
            &     this%contents(3,3)
end function

pure function vec3_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  type(Vector3)            :: output
  
  output%contents = dble(split(this))
end function

pure function mat33_Strings(this) result(output)
  implicit none
  
  type(String), intent(in) :: this(3)
  type(Matrix33)           :: output
  
  integer        :: i
  
  do i=1,3
    output%contents(i,:) = dble(split(this(i)))
  enddo
end function
end module
