! ======================================================================
! Defines type-specific array operations.
! ======================================================================
! N.B. MACRO_TYPE_NAME must be defined before this file is included.
! N.B. this file should be included twice for each type:
!    - once before the 'contains' line. MACRO_BODY must not be defined.
!    - once after the 'contains' line. MACRO_BODY must be defined.
! N.B. This file can be included for multiple types,
!    with a different MACRO_TYPE_NAME each time.
! N.B. this file may only be included in files which use caesar_io_module.

! ----------------------------------------------------------------------
! Functionality defined in this file:
! ----------------------------------------------------------------------
! append:
!    If x is allocated then append(x,y) is equivalent to x = [x,y].
!    If x is not allocated then append(x,y) is equivalent to x = [y].
!    In both cases, using append avoids a memory leak in gfortran.
! ----------------------------------------------------------------------

! Include useful macros.
#include "macros.fpp"

! Define MACRO_TYPE_DECLARATION from MACRO_TYPE_NAME.
#include "type_declaration.fpp"

! Check whether this inclusion is in the header (before 'contains'),
!    or the body (after 'contains').
#ifndef MACRO_BODY
! ----------------------------------------------------------------------
! The header:
!    - public declarations.
!    - interfaces.
! ----------------------------------------------------------------------

public :: append

interface append
  module procedure MACRO_CAT(append_scalar_,MACRO_TYPE_NAME)
  module procedure MACRO_CAT(append_array_,MACRO_TYPE_NAME)
end interface

! ----------------------------------------------------------------------
#else
! ----------------------------------------------------------------------
! The body:
!    - subroutine and function definitions.
! ----------------------------------------------------------------------

! append_scalar_TYPE_NAME
subroutine MACRO_CAT(append_scalar_,MACRO_TYPE_NAME)(lhs,rhs)
  implicit none
  
  MACRO_TYPE_DECLARATION, intent(inout), allocatable :: lhs(:)
  MACRO_TYPE_DECLARATION, intent(in)                 :: rhs
  
  MACRO_TYPE_DECLARATION, allocatable :: temp(:)
  
  integer :: ialloc
  
  if (.not. allocated(lhs)) then
    allocate(lhs(1), stat=ialloc); call err(ialloc)
    lhs(1) = rhs
  else
    allocate(temp(size(lhs)+1), stat=ialloc); call err(ialloc)
    temp(:size(lhs)) = lhs(:)
    temp(size(lhs)+1) = rhs
    deallocate(lhs, stat=ialloc); call err(ialloc)
    allocate(lhs(size(temp)), stat=ialloc); call err(ialloc)
    lhs(:) = temp(:)
    deallocate(temp, stat=ialloc); call err(ialloc)
  endif
end subroutine

! append_array_TYPE_NAME
subroutine MACRO_CAT(append_array_,MACRO_TYPE_NAME)(lhs,rhs)
  implicit none
  
  MACRO_TYPE_DECLARATION, intent(inout), allocatable :: lhs(:)
  MACRO_TYPE_DECLARATION, intent(in)                 :: rhs(:)
  
  MACRO_TYPE_DECLARATION, allocatable :: temp(:)
  
  integer :: ialloc
  
  if (.not. allocated(lhs)) then
    allocate(lhs(size(rhs)), stat=ialloc); call err(ialloc)
    lhs(:) = rhs(:)
  else
    allocate(temp(size(lhs)+size(rhs)), stat=ialloc); call err(ialloc)
    temp(:size(lhs)) = lhs(:)
    temp(size(lhs)+1:) = rhs(:)
    deallocate(lhs, stat=ialloc); call err(ialloc)
    allocate(lhs(size(temp)), stat=ialloc); call err(ialloc)
    lhs(:) = temp(:)
    deallocate(temp, stat=ialloc); call err(ialloc)
  endif
end subroutine

! ----------------------------------------------------------------------
#endif
