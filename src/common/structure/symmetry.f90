! ======================================================================
! A symmetry operator, in various representations.
! ======================================================================
module symmetry_submodule
  use utils_module
  
  use basic_structure_submodule
  use calculate_symmetry_submodule
  use qpoint_submodule
  implicit none
  
  private
  
  public :: SymmetryOperator
  public :: operators_commute
  
  ! ----------------------------------------------------------------------
  ! A symmetry operation.
  ! ----------------------------------------------------------------------
  type :: SymmetryOperator
    ! The rotation and translation in fractional co-ordinates.
    ! R and T.
    type(IntMatrix)  :: rotation
    type(RealVector) :: translation
    
    ! The rotation and translation in cartesian co-ordinates.
    ! (L^T)R(L^-T) and (L^T)T.
    type(RealMatrix) :: cartesian_rotation
    type(RealVector) :: cartesian_translation
    
    ! The rotation in fractional reciprocal co-ordinates.
    ! Used for rotating q-points, and other reciprocal space objects.
    type(FractionMatrix), private :: recip_rotation_
    
    ! The mapping from atoms to other atoms.
    ! rho_i and rho_j are the equilibrium positions of atom i and j.
    ! R * rho_i + T = rho_j + rvector_i
    !    - atom_group*i = j,
    !    - rvector(i) = rvector_i in fractional supercell co-ordinates.
    type(Group)                  :: atom_group
    type(IntVector), allocatable :: rvector(:)
    
    ! The mapping from atoms in the primitive cell to other atoms in the
    !    primitive cell.
    ! The entries in prim_rvector are given in fractional primitive cell
    !    co-ordinates.
    type(Group)                  :: prim_atom_group
    type(IntVector), allocatable :: prim_rvector(:)
  contains
    ! Rotating a q-point.
    generic,   public  :: operator(*) => rotate_qpoint
    procedure, private ::                rotate_qpoint
    
    procedure, public :: symmetry_order
  end type
  
  interface SymmetryOperator
    module procedure new_SymmetryOperator
  end interface
contains

! ----------------------------------------------------------------------
! Constructs a SymmetryOperator.
! ----------------------------------------------------------------------
function new_SymmetryOperator(symmetry,lattice,recip_lattice,atom_group, &
   & rvector,prim_atom_group,prim_rvector) result(output)
  implicit none
  
  type(BasicSymmetry), intent(in) :: symmetry
  type(RealMatrix),    intent(in) :: lattice
  type(RealMatrix),    intent(in) :: recip_lattice
  type(Group),         intent(in) :: atom_group
  type(IntVector),     intent(in) :: rvector(:)
  type(Group),         intent(in) :: prim_atom_group
  type(IntVector),     intent(in) :: prim_rvector(:)
  type(SymmetryOperator)          :: output
  
  output%rotation = symmetry%rotation
  output%translation = symmetry%translation
  
  output%cartesian_rotation = transpose(lattice) &
                          & * symmetry%rotation  &
                          & * recip_lattice
  output%cartesian_translation = transpose(lattice) &
                             & * symmetry%translation
  
  output%recip_rotation_ = transpose(invert(symmetry%rotation))
  
  output%atom_group      = atom_group
  output%rvector         = rvector
  output%prim_atom_group = prim_atom_group
  output%prim_rvector    = prim_rvector
end function

! ----------------------------------------------------------------------
! Returns whether or not two symmetry operators commute.
! ----------------------------------------------------------------------
function operators_commute(this,that) result(output)
  implicit none
  
  type(SymmetryOperator), intent(in) :: this
  type(SymmetryOperator), intent(in) :: that
  logical                            :: output
  
  ! Check that the symmetries' rotations commute.
  if (.not. matrices_commute(this%rotation,that%rotation)) then
    output = .false.
    return
  endif
  
  ! Check that the atom to atom transformations commute.
  if (this%atom_group*that%atom_group /= that%atom_group*this%atom_group) then
    output = .false.
    return
  endif
  
  output = .true.
end function

! ----------------------------------------------------------------------
! Rotates a q-point, and translates it back to the primitive reciprocal cell.
! ----------------------------------------------------------------------
impure elemental function rotate_qpoint(this,qpoint) result(output)
  implicit none
  
  class(SymmetryOperator), intent(in) :: this
  type(QpointData),        intent(in) :: qpoint
  type(QpointData)                    :: output
  
  ! Transfer across all metadata.
  output = qpoint
  
  ! Rotate q-point.
  output%qpoint = this%recip_rotation_ * qpoint%qpoint
  
  ! Translate q-point back into primitive reciprocal cell if necessary.
  call output%translate_to_primitive()
end function

! ----------------------------------------------------------------------
! Calculates the order of a symmetry operation.
! The order is the smallest integer n>0 s.t. S^n=I, where I is the identity.
! If qpoint is given, then the S^n is only considered to be the identity when
!    the associated translation, R, is such that q.R is a multiple of two*pi.
! ----------------------------------------------------------------------
function symmetry_order(this,qpoint) result(output)
  implicit none
  
  class(SymmetryOperator), intent(in)           :: this
  type(QpointData),        intent(in), optional :: qpoint
  integer                                       :: output
  
  type(IntMatrix)   :: identity ! The identity matrix, I.
  type(IntMatrix)   :: rotation ! The rotation after n operations.
  type(IntVector)   :: rvector  ! The translation after n operations.
  integer           :: atom_1   ! atom 1 maps to this after n operations.
  type(IntFraction) :: qr       ! qpoint*rvector.
  
  integer :: i
  
  identity = make_identity_matrix(3)
  
  ! After (S^i) : r -> rotation * r + rvector.
  ! After (S^0) : r -> r = I*r+0.
  rotation = identity
  rvector  = zeroes(3)
  atom_1   = 1
  
  ! Calculate n s.t. R^n=I. At this point, S^n=e^{2piiq.R}.
  output = 0
  do i=1,6
    ! After S^i, r -> rotation * r + rvector.
    rotation = this%rotation * rotation
    rvector = this%rotation * rvector + this%prim_rvector(atom_1)
    atom_1 = this%atom_group * atom_1
    if (rotation==identity) then
      if (atom_1/=1) then
        call print_line(CODE_ERROR//': Unable to find order of symmetry.')
        call err()
      endif
      output = i
      exit
    endif
  enddo
  
  if (output==0) then
    call print_line(CODE_ERROR//': Unable to find order of symmetry.')
    call err()
  endif
  
  ! Include the effect of the phase change.
  ! If S^n=e^{2piiq.R}, and q.R=a/b then S^(n*b)=I.
  if (present(qpoint)) then
    qr = qpoint%qpoint * rvector
    output = output * qr%denominator()
  endif
end function
end module
