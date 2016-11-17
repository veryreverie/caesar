! ------------------------------------------------------------
! various utilities
! ------------------------------------------------------------
module utils
  use constants,      only : dp
  use linear_algebra, only : inv_33
  implicit none
  
!  interface determinant33
!    module procedure determinant33_integer, determinant33_real
!  end interface

contains
  
  ! ----------------------------------------
  ! reports an error and stops
  ! ----------------------------------------
  subroutine errstop(subroutine_name, message)
    implicit none
    
    character(*), intent(in) :: subroutine_name ! where errstop is called
    character(*), intent(in) :: message         ! error message
    
    write(*,*) ''
    write(*,*) 'error in subroutine '//trim(adjustl(subroutine_name))//'.'
    write(*,*) ''
    call wordwrap(trim(adjustl(message)))
    write(*,*) ''
    stop
  end subroutine
  
  ! ----------------------------------------
  ! reports an allocation error and stops
  ! ----------------------------------------
  subroutine erralloc(arg)
    implicit none
    
    character(*), intent(in) :: arg
    
    write(*,*) ''
    write(*,*) 'Problem allocating '//trim(adjustl(arg))//' array.'
    write(*,*)
  end subroutine
  
  ! ----------------------------------------
  ! prints out the contents of the character string 'text',
  ! ensuring that line breaks only occur at space characters. The output
  ! is written to unit unit_in if this parameter is supplies; otherwise the
  ! output is written to unit o. The maximum length of each line is given
  ! by linelength_in if this is supplied; otherwise the default line length
  ! is 79 characters.
  ! ----------------------------------------
  subroutine wordwrap(text, unit_in, linelength_in)
    implicit none
    
    integer, intent(in), optional :: unit_in
    integer, intent(in), optional :: linelength_in
    character(*), intent(in)      :: text
    character(260)                :: temp
    integer                       :: i
    integer                       :: unit
    integer                       :: lentext
    integer                       :: startpoint
    integer                       :: stoppoint
    integer                       :: lastpos
    integer                       :: linelength
    
    ! check if unit_in is supplied
    if (present(unit_in)) then
      unit = unit_in
    else
      unit=6
    endif
    
    lentext = len(trim(text))
    if (lentext<1) then ! text is empty
      write(unit,*) ""
      return
    endif
    
    ! check if linelength_in is supplied
    if (present(linelength_in)) then
      if (linelength_in>=2) then
        linelength=linelength_in
      else
        linelength=2
      endif
    else
      linelength=79
    endif
    
    startpoint=1
    do i=1, huge(1) ! loop over lines
      stoppoint = startpoint+linelength-1
      if (stoppoint<=lentext) then
        lastpos = index(trim(text(startpoint:stoppoint)),' ',.true.)
        if (lastpos>0) stoppoint = startpoint+lastpos-1
      else
        stoppoint = lentext
      endif
      
      if (i==1) then
        ! allow the user to indent the first line, if they wish
        temp = text(startpoint:stoppoint) ! or pathscale.f90 fails to compile
        write(unit,*) trim(temp)
      else
        temp = text(startpoint:stoppoint) ! or pathscale.f90 fails to compile
        write(unit,*) trim(adjustl(temp))
      endif
      
      if (stoppoint==lentext) then
        exit
      else
        startpoint = stoppoint+1
      endif ! finished text?
    enddo
  end subroutine
  
  ! ----------------------------------------
  ! converts integers to left justified strings that can be printed in the
  ! middle of a sentence without introducing large amounts of white space.
  ! ----------------------------------------
  ! TODO: There are much better ways of doing this.
  character(12) function i2s(n)
    implicit none
    
    integer, intent(in) :: n                 ! input integer
    integer             :: i                 ! characters left to process
    integer             :: j                 ! loop counter
    integer, parameter  :: ichar0=ichar('0') ! the character '0'
    
    i2s = ''
    i=abs(n)
    do j=len(i2s),1,-1
      i2s(j:j)=achar(ichar0+mod(i,10))
      i=i/10
      if (i==0) exit
    enddo
    if (n<0) then
      i2s='-'//adjustl(i2s)
    else
      i2s=adjustl(i2s)
    endif
  end function
  
!  ! ----------------------------------------
!  ! given a 3x3 matrix A, returns det(A)
!  ! ----------------------------------------
!  function determinant33_integer(A) result(determinant)
!    implicit none
!    
!    integer, intent(in) :: A(3,3)
!    integer             :: determinant
!    
!    determinant = A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))&
!               &+ A(1,2)*(A(3,1)*A(2,3)-A(2,1)*A(3,3))&
!               &+ A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
!  end function
!  
!  function determinant33_real(A) result(determinant)
!    implicit none
!    
!    real(dp), intent(in) :: A(3,3)
!    real(dp)             :: determinant
!    
!    determinant = A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))&
!               &+ A(1,2)*(A(3,1)*A(2,3)-A(2,1)*A(3,3))&
!               &+ A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
!  end function
!  
!  ! ----------------------------------------
!  ! calculates the inverse, B, of matrix A
!  ! A and B are real, 3x3 matrices
!  ! ----------------------------------------
!  ! TODO: |d|<epsilon would be a more stable check than d==0
!  subroutine inv_33(A,B)
!    implicit none
!    
!    real(dp), intent(in)  :: A(3,3)
!    real(dp), intent(out) :: B(3,3)
!    real(dp)              :: d      ! det(A)
!    
!    d = A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))&
!     &+ A(1,2)*(A(3,1)*A(2,3)-A(2,1)*A(3,3))&
!     &+ A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
!    
!    if (d==0.d0) then
!      write(*,*) 'Error in inv_33: singular matrix.'
!      stop
!    endif
!    
!    d = 1.d0/d
!    
!    B(1,1) = (A(2,2)*A(3,3)-A(2,3)*A(3,2))*d
!    B(1,2) = (A(3,2)*A(1,3)-A(1,3)*A(3,2))*d
!    B(1,3) = (A(1,2)*A(2,3)-A(1,3)*A(3,2))*d
!    B(2,1) = (A(3,1)*A(2,3)-A(2,3)*A(3,2))*d
!    B(2,2) = (A(1,1)*A(3,3)-A(2,3)*A(3,2))*d
!    B(2,3) = (A(2,1)*A(1,3)-A(2,3)*A(3,2))*d
!    B(3,1) = (A(2,1)*A(3,2)-A(2,3)*A(3,2))*d
!    B(3,2) = (A(3,1)*A(1,2)-A(2,3)*A(3,2))*d
!    B(3,3) = (A(1,1)*A(2,2)-A(2,3)*A(3,2))*d
!  end subroutine
  
  ! ----------------------------------------
  ! checks if pos is in the supercell
  ! ----------------------------------------
  logical function is_in_supercell(pos, super_lattice)
    implicit none
    
    real(dp), intent(in) :: pos(3)
    real(dp), intent(in) :: super_lattice(3,3)
    real(dp)             :: trans_super_lattice(3,3)
    real(dp)             :: inv_super_lattice(3,3)
    real(dp)             :: frac_pos(3)
    real(dp)             :: t
    real(dp)             :: a
    real(dp)             :: b
    real(dp)             :: c
    real(dp)             :: f1
    real(dp)             :: f2
    real(dp)             :: f3
    real(dp), parameter  :: tol=1.d-2
    
    trans_super_lattice = transpose(super_lattice)
    call inv_33(trans_super_lattice, inv_super_lattice)
    frac_pos(1:3) = pos(1)*inv_super_lattice(1:3,1)&
                 &+ pos(2)*inv_super_lattice(1:3,2)&
                 &+ pos(3)*inv_super_lattice(1:3,3)
    
    if (frac_pos(1)>-tol .and. frac_pos(1)<(1.d0-tol)) then
      if (frac_pos(2)>-tol .and. frac_pos(2)<(1.d0-tol)) then
        is_in_supercell = (frac_pos(3)>-tol .and. frac_pos(3)<(1.d0-tol))
      else
        is_in_supercell = .false.
      endif
    else
      is_in_supercell = .false.
    endif
  end function
  
  ! ----------------------------------------
  ! returns an array containing the command line arguments
  ! ----------------------------------------
  function command_line_args() result(args)
    implicit none

    integer                        :: i         ! loop index
    integer                        :: arg_count ! no. of command line args
    character(len=32), allocatable :: args(:)   ! return value

    ! read in the number of arguments
    arg_count = iargc()
    ! allocate the return array
    allocate (args(arg_count))
    ! read the arguments into the array
    do i=1,arg_count
      call getarg(i, args(i))
    enddo

    return
  end function
  
end module
