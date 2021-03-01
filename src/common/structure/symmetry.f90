! ======================================================================
! A symmetry operator, in various representations.
! ======================================================================
module caesar_symmetry_module
  use caesar_utils_module
  
  use caesar_atom_module
  
  use caesar_basic_symmetry_module
  use caesar_qpoint_module
  implicit none
  
  private
  
  public :: SymmetryOperator
  public :: operator(*)
  public :: BasicSymmetry
  public :: operators_commute
  public :: superposed_operators_commute
  public :: operators_anticommute
  public :: loto_breaks_symmetry
  
  ! ----------------------------------------------------------------------
  ! A symmetry operation.
  ! ----------------------------------------------------------------------
  type :: SymmetryOperator
    ! The ID of the symmetry.
    integer :: id
    
    ! A symmetry is defined by the tensor, T, and translation, t,
    !    both defined in fractional co-ordinates.
    ! The symmetry S acts on the fractional co-ordinate r as:
    !    r -> T.r + t
    type(IntMatrix)  :: tensor
    type(RealVector) :: translation
    
    ! The tensor and translation in cartesian co-ordinates.
    ! (L^T)T(L^-T) and (L^T)t.
    type(RealMatrix) :: cartesian_tensor
    type(RealVector) :: cartesian_translation
    
    ! The tensor in fractional reciprocal co-ordinates.
    ! Used for transforming q-points, and other reciprocal space objects.
    type(FractionMatrix), private :: recip_tensor_
    
    ! The mapping from atoms to other atoms.
    ! r_i and r_j are the equilibrium positions of atom i and j.
    ! T * r_i + t = r_j + R_i
    !    - atom_group*i = j,
    !    - rvectors(i) = R_i in fractional supercell co-ordinates.
    type(Group)                  :: atom_group
    type(IntVector), allocatable :: rvectors(:)
  contains
    procedure, public :: inverse_transform
    
    procedure, public :: symmetry_order
  end type
  
  interface SymmetryOperator
    ! ----------------------------------------------------------------------
    ! Constructs a SymmetryOperator.
    ! ----------------------------------------------------------------------
    impure elemental module function new_SymmetryOperator(symmetry,lattice, &
       & recip_lattice) result(this) 
      type(BasicSymmetry), intent(in) :: symmetry
      type(RealMatrix),    intent(in) :: lattice
      type(RealMatrix),    intent(in) :: recip_lattice
      type(SymmetryOperator)          :: this
    end function
  end interface
  
  interface BasicSymmetry
    ! ----------------------------------------------------------------------
    ! Constructs a BasicSymmetry from a SymmetryOperator.
    ! ----------------------------------------------------------------------
    impure elemental module function new_BasicSymmetry_SymmetryOperator(this) &
       & result(output) 
      type(SymmetryOperator), intent(in) :: this
      type(BasicSymmetry)                :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns whether or not two symmetry operators commute or anti-commute.
    ! ----------------------------------------------------------------------
    ! Each operator is a tensor product of the cartesian tensor and an atom tensor.
    ! Two operators commute if:
    !    - Their cartesian and atom tensors both commute.
    !    - Their cartesian and atom tensors both anti-commute.
    ! Two operators anti-commute if:
    !    - One tensor commutes and the other anticommutes.
    impure elemental module function operators_commute(a,b,qpoint) &
       & result(output) 
      type(SymmetryOperator), intent(in) :: a
      type(SymmetryOperator), intent(in) :: b
      type(QpointData),       intent(in) :: qpoint
      logical                            :: output
    end function
  end interface
  
  interface
    impure elemental module function operators_anticommute(a,b,qpoint) &
       & result(output) 
      type(SymmetryOperator), intent(in) :: a
      type(SymmetryOperator), intent(in) :: b
      type(QpointData),       intent(in) :: qpoint
      logical                            :: output
    end function
  end interface
  
  interface
    ! For two symmetries A and B: A1=A,A2=A*,B1=B,B2=B*.
    !    if positive_superposition=.true.  : returns whether [A1+A2,B1+B2]=0.
    !    if positive_superposition=.false. : returns whether [A1-A2,B1-B2]=0.
    ! Each symmetry is the tensor product of the cartesian tensor T and the
    !    atom tensor G: A1=TA1^GA1, A2=TA2^GA2, B1=TB1^GB1 and B2=TB2^GB2.
    ! [A1+A2,B1+B2] = A1B1+A1B2+A2B1+A2B2-B1A1-B1A2-B2A1-B2A2
    ! [A1-A2,B1-B2] = A1B1-A1B2-A2B1+A2B2-B1A1+B1A2+B2A1-B2A2
    impure elemental module function superposed_operators_commute(a,b, &
       & qpoint,positive_superposition) result(output) 
      type(SymmetryOperator), intent(in) :: a
      type(SymmetryOperator), intent(in) :: b
      type(QpointData),       intent(in) :: qpoint
      logical,                intent(in) :: positive_superposition
      logical                            :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Helper functions for operators_commute and operators_anticommute,
    !    returning whether or not the operators' R-vectors (anti-)commute.
    ! ----------------------------------------------------------------------
    ! If the primitive cell has N atoms then the atom tensor is an N*N matrix,
    !    with N non-zero elements corresponding to each r_i -> r_j + R_i,
    !    so the ij element is e^(2*pi*i*q*R_i).
    impure elemental module function atom_tensors_commute(a,b,qpoint) &
       & result(output) 
      type(SymmetryOperator), intent(in) :: a
      type(SymmetryOperator), intent(in) :: b
      type(QpointData),       intent(in) :: qpoint
      logical                            :: output
    end function
  end interface
  
  interface
    impure elemental module function atom_tensors_anticommute(a,b,qpoint) &
       & result(output) 
      type(SymmetryOperator), intent(in) :: a
      type(SymmetryOperator), intent(in) :: b
      type(QpointData),       intent(in) :: qpoint
      logical                            :: output
    end function
  end interface
  
  interface operator(*)
    ! ----------------------------------------------------------------------
    ! Transforms a q-point, and translates it back to
    !    the primitive reciprocal cell.
    ! ----------------------------------------------------------------------
    impure elemental module function transform_qpoint(this,qpoint) &
       & result(output) 
      class(SymmetryOperator), intent(in) :: this
      type(QpointData),        intent(in) :: qpoint
      type(QpointData)                    :: output
    end function
  end interface
  
  interface
    impure elemental module function inverse_transform(this,qpoint) &
       & result(output) 
      class(SymmetryOperator), intent(in) :: this
      type(QpointData),        intent(in) :: qpoint
      type(QpointData)                    :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Calculates the order of a symmetry operation.
    ! The order is the smallest integer n>0 s.t. S^n=I, where I is the identity.
    ! If qpoint is given, then the S^n is only considered to be the identity when
    !    the associated translation, R, is such that q.R is a multiple of two*pi.
    ! ----------------------------------------------------------------------
    impure elemental module function symmetry_order(this,qpoint) &
       & result(output) 
      class(SymmetryOperator), intent(in)           :: this
      type(QpointData),        intent(in), optional :: qpoint
      integer                                       :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns whether or not the LO/TO direction breaks a given symmetry.
    ! ----------------------------------------------------------------------
    ! A symmetry is broken if its tensor transforms the LO/TO direction
    !    onto a different axis.
    ! Symmetries which map the LO/TO direction onto minus itself are allowed.
    impure elemental module function loto_breaks_symmetry(tensor, &
       & loto_direction) result(output) 
      type(IntMatrix),      intent(in) :: tensor
      type(FractionVector), intent(in) :: loto_direction
      logical                          :: output
    end function
  end interface
end module
