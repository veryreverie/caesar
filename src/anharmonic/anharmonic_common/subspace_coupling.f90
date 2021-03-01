! ======================================================================
! A "subspace" is a set of degenerate modes.
! "Coupled subspaces" are a set of subspaces which are coupled to one another.
! Coupling here means that there are terms in the Hamiltonian which contain
!    products of one subspace with another.
! ======================================================================
module caesar_subspace_coupling_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: SubspaceCoupling
  public :: generate_coupled_subspaces
  public :: size
  public :: operator(==)
  public :: operator(/=)
  
  type, extends(Stringable) :: SubspaceCoupling
    ! The ids of the degenerate subspaces which are coupled together.
    ! ids is assumed to always be sorted in ascending order.
    integer, allocatable :: ids(:)
  contains
    ! Returns the subspaces in this coupling.
    procedure, public :: coupled_subspaces
    ! I/O.
    procedure, public :: read  => read_SubspaceCoupling
    procedure, public :: write => write_SubspaceCoupling
  end type
  
  interface SubspaceCoupling
    ! ----------------------------------------------------------------------
    ! Basic functionality:
    !    - Constructor.
    !    - Concatenation with a DegenerateSubspace subspace.
    !    - size() module function.
    !    - == and /= operators.
    ! ----------------------------------------------------------------------
    module function new_SubspaceCoupling(ids) result(output) 
      integer, intent(in), optional :: ids(:)
      type(SubspaceCoupling)        :: output
    end function
  end interface
  
  interface operator(//)
    module function concatenate_SubspaceCoupling_DegenerateSubspace(input, &
       & subspace) result(output) 
      type(SubspaceCoupling),   intent(in) :: input
      type(DegenerateSubspace), intent(in) :: subspace
      type(SubspaceCoupling)               :: output
    end function
  end interface
  
  interface size
    module function size_SubspaceCoupling(this) result(output) 
      type(SubspaceCoupling), intent(in) :: this
      integer                            :: output
    end function
  end interface
  
  interface operator(==)
    impure elemental module function equality_SubspaceCoupling(this,that) &
       & result(output) 
      type(SubspaceCoupling), intent(in) :: this
      type(SubspaceCoupling), intent(in) :: that
      logical                            :: output
    end function
  end interface
  
  interface operator(/=)
    impure elemental module function non_equality_SubspaceCoupling(this,that) &
       & result(output) 
      type(SubspaceCoupling), intent(in) :: this
      type(SubspaceCoupling), intent(in) :: that
      logical                            :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Generates all sets of coupled subspaces up to a given order.
    ! ----------------------------------------------------------------------
    ! Calls its helper module function,which recursively calls itself.
    module function generate_coupled_subspaces(subspaces, &
       & maximum_coupling_order) result(output) 
      type(DegenerateSubspace), intent(in) :: subspaces(:)
      integer,                  intent(in) :: maximum_coupling_order
      type(SubspaceCoupling), allocatable  :: output(:)
    end function
  end interface
  
  interface
    ! Recursive helper module function for generate_coupled_subspaces.
    recursive module function generate_coupled_subspaces_helper(coupling_in, &
       & subspaces,coupling_order) result(output) 
      type(SubspaceCoupling),   intent(in) :: coupling_in
      type(DegenerateSubspace), intent(in) :: subspaces(:)
      integer,                  intent(in) :: coupling_order
      type(SubspaceCoupling), allocatable  :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns the coupled subspaces in this subspace coupling.
    ! ----------------------------------------------------------------------
    module function coupled_subspaces(this,subspaces) result(output) 
      class(SubspaceCoupling),  intent(in)  :: this
      type(DegenerateSubspace), intent(in)  :: subspaces(:)
      type(DegenerateSubspace), allocatable :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_SubspaceCoupling(this,input) 
      class(SubspaceCoupling), intent(out) :: this
      type(String),            intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_SubspaceCoupling(this) result(output) 
      class(SubspaceCoupling), intent(in) :: this
      type(String)                        :: output
    end function
  end interface
  
  interface SubspaceCoupling
    impure elemental module function new_SubspaceCoupling_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(SubspaceCoupling)   :: this
    end function
  end interface
end module
