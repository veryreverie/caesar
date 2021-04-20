! ======================================================================
! A "subspace" is a set of degenerate modes.
! A subspace monomial is a product of subspaces.
! ======================================================================
module caesar_subspace_monomial_module
  use caesar_common_module
  
  use caesar_subspace_coupling_module
  implicit none
  
  private
  
  public :: SubspaceMonomial
  public :: size
  public :: generate_subspace_monomials
  public :: operator(//)
  public :: operator(==)
  public :: operator(/=)
  public :: generate_complex_monomials
  public :: generate_paired_monomials
  
  type, extends(Stringable) :: SubspaceMonomial
    ! The ids and powers of the degenerate subspaces in the monomial.
    ! e.g. if ids=[1,2,5] and powers=[1,2,1] then the monomial is
    !    (s1).(s2)^2.(s5), where s1 is the subspace with id 1.
    integer, allocatable :: ids(:)
    integer, allocatable :: powers(:)
  contains
    ! Returns the coupled subspaces.
    procedure, public :: coupled_subspaces
    
    ! Check if this is subsidiary to the given coupling, e.g. the couplings
    !    [],[1],[2],[3],[1,2],[1,3] and [2,3] are subsidiaries of [1,2,3].
    procedure, public :: is_subsidiary_of
    
    ! I/O.
    procedure, public :: read  => read_SubspaceMonomial
    procedure, public :: write => write_SubspaceMonomial
  end type
  
  interface SubspaceMonomial
    ! ----------------------------------------------------------------------
    ! Basic functionality.
    !    - Constructors.
    !    - Concatenation with DegenerateSubspace.
    !    - size() module function.
    !    - equality and non-equality with other SubspaceMonomials.
    ! ----------------------------------------------------------------------
    module function new_SubspaceMonomial() result(this) 
      type(SubspaceMonomial) :: this
    end function
  
    module function new_SubspaceMonomial_ids_powers(ids,powers) result(this) 
      integer, intent(in)    :: ids(:)
      integer, intent(in)    :: powers(:)
      type(SubspaceMonomial) :: this
    end function
  
    module function new_SubspaceMonomial_DegenerateSubspaces(subspaces) &
       & result(this) 
      type(DegenerateSubspace), intent(in) :: subspaces(:)
      type(SubspaceMonomial)               :: this
    end function
  end interface
  
  interface operator(//)
    module function concatenate_SubspaceMonomial_DegenerateSubspace(this, &
       & subspace) result(output) 
      type(SubspaceMonomial),   intent(in) :: this
      type(DegenerateSubspace), intent(in) :: subspace
      type(SubspaceMonomial)               :: output
    end function
  end interface
  
  interface size
    module function size_SubspaceMonomial(this) result(output) 
      type(SubspaceMonomial), intent(in) :: this
      integer                            :: output
    end function
  end interface
  
  interface operator(==)
    impure elemental module function equality_SubspaceMonomial_SubspaceMonomial(this,that) result(output) 
      type(SubspaceMonomial), intent(in) :: this
      type(SubspaceMonomial), intent(in) :: that
      logical                            :: output
    end function
  end interface
  
  interface operator(/=)
    impure elemental module function non_equality_SubspaceMonomial_SubspaceMonomial(   this,that) result(output) 
      type(SubspaceMonomial), intent(in) :: this
      type(SubspaceMonomial), intent(in) :: that
      logical                            :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Use stored IDs to return the actual objects the IDs represent.
    ! ----------------------------------------------------------------------
    module function coupled_subspaces(this,subspaces) result(output) 
      class(SubspaceMonomial),  intent(in)  :: this
      type(DegenerateSubspace), intent(in)  :: subspaces(:)
      type(DegenerateSubspace), allocatable :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Check if a coupling is a subsidiary of another coupling.
    ! ----------------------------------------------------------------------
    module function is_subsidiary_of(this,that) result(output) 
      Class(SubspaceMonomial), intent(in) :: this
      type(SubspaceMonomial),  intent(in) :: that
      logical                             :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Generates all coupling monomials corresponding to a given subspace coupling,
    !    up to the given potential expansion order.
    ! ----------------------------------------------------------------------
    ! e.g. given the subspace coupling with ids [3,5] and expansion order 3,
    !    will return coupling monomials [3,3,5] and [3,5,5].
    module function generate_subspace_monomials(subspace_coupling,subspaces, &
       & minimum_expansion_order,maximum_expansion_order) result(output) 
      type(SubspaceCoupling),   intent(in) :: subspace_coupling
      type(DegenerateSubspace), intent(in) :: subspaces(:)
      integer,                  intent(in) :: minimum_expansion_order
      integer,                  intent(in) :: maximum_expansion_order
      type(SubspaceMonomial), allocatable  :: output(:)
    end function
  end interface
  
  interface
    ! Helper module function for generate_subspace_monomials.
    ! This module function appends a number of copies of the first subspace in
    !    coupled_subspaces, then calls itself to append the next subspace etc.
    ! The optional argument monomial_in is a recursive argument which should only
    !    be provided by recursive calls.
    recursive module function generate_subspace_monomials_helper(coupled_subspaces,minimum_expansion_order,maximum_expansion_order,monomial_in) result(output) 
      type(DegenerateSubspace), intent(in)           :: coupled_subspaces(:)
      integer,                  intent(in)           :: minimum_expansion_order
      integer,                  intent(in)           :: maximum_expansion_order
      type(SubspaceMonomial),   intent(in), optional :: monomial_in
      type(SubspaceMonomial), allocatable            :: output(:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Generate the ComplexMonomials corresponding to a given SubspaceMonomial.
    ! ----------------------------------------------------------------------
    ! The coefficients are chosen such that symmetry operations are unitary
    !    in the basis of monomials.
    module function generate_complex_monomials(this,maximum_coupling_order, &
       & subspaces,modes,qpoints,conserve_momentum,                         &
       & conserve_subspace_momentum) result(output) 
      type(SubspaceMonomial),   intent(in) :: this
      integer,                  intent(in) :: maximum_coupling_order
      type(DegenerateSubspace), intent(in) :: subspaces(:)
      type(ComplexMode),        intent(in) :: modes(:)
      type(QpointData),         intent(in) :: qpoints(:)
      logical,                  intent(in) :: conserve_momentum
      logical,                  intent(in) :: conserve_subspace_momentum
      type(ComplexMonomial), allocatable   :: output(:)
    end function
  end interface
  
  interface
    module function generate_paired_monomials(this,maximum_coupling_order, &
       & subspaces,modes,qpoints,conserve_momentum,                        &
       & conserve_subspace_momentum) result(output) 
      type(SubspaceMonomial),   intent(in) :: this
      integer,                  intent(in) :: maximum_coupling_order
      type(DegenerateSubspace), intent(in) :: subspaces(:)
      type(ComplexMode),        intent(in) :: modes(:)
      type(QpointData),         intent(in) :: qpoints(:)
      logical,                  intent(in) :: conserve_momentum
      logical,                  intent(in) :: conserve_subspace_momentum
      type(PairedMonomial), allocatable    :: output(:)
    end function
  end interface
  
  interface
    ! Helper functions for the above. These generate monomials for a single
    !    subspace.
    recursive module function generate_subspace_complex_monomials(modes, &
       & power,root) result(output) 
      type(ComplexMode),     intent(in)  :: modes(:)
      integer,               intent(in)  :: power
      type(ComplexMonomial), intent(in)  :: root
      type(ComplexMonomial), allocatable :: output(:)
    end function
  end interface
  
  interface
    module function no_permutations(input) result(output) 
      type(ComplexMonomial), intent(in) :: input
      real(dp)                          :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_SubspaceMonomial(this,input) 
      class(SubspaceMonomial), intent(out) :: this
      type(String),            intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_SubspaceMonomial(this) result(output) 
      class(SubspaceMonomial), intent(in) :: this
      type(String)                        :: output
    end function
  end interface
  
  interface SubspaceMonomial
    impure elemental module function new_SubspaceMonomial_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(SubspaceMonomial)   :: this
    end function
  end interface
end module
