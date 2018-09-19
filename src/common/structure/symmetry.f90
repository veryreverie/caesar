! ======================================================================
! A symmetry operator, in various representations.
! ======================================================================
module symmetry_submodule
  use utils_module
  
  use basic_symmetry_submodule
  use basic_structure_submodule
  use qpoint_submodule
  implicit none
  
  private
  
  public :: SymmetryOperator
  public :: operator(*)
  public :: BasicSymmetry
  public :: operators_commute
  
  ! ----------------------------------------------------------------------
  ! A symmetry operation.
  ! ----------------------------------------------------------------------
  type :: SymmetryOperator
    ! The ID of the symmetry.
    integer :: id
    
    ! The ID of the symmetry S s.t. this symmetry * S is the identity,
    !    up to R-vector translation.
    integer :: inverse_symmetry_id
    
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
    
    ! How the symmetry acts on other symmetries.
    ! If this symmetry * symmetry_j = symmetry_k then symmetry_group*j=k.
    ! N.B. these relations are true only up to R-vector translations.
    type(Group) :: symmetry_group
  contains
    procedure, public :: inverse_transform
    
    procedure, public :: symmetry_order
  end type
  
  interface SymmetryOperator
    module procedure new_SymmetryOperator
  end interface
  
  ! Transform a q-point.
  interface operator(*)
    module procedure transform_qpoint
  end interface
  
  interface BasicSymmetry
    module procedure new_BasicSymmetry_SymmetryOperator
  end interface
contains

! ----------------------------------------------------------------------
! Constructs a SymmetryOperator.
! ----------------------------------------------------------------------
function new_SymmetryOperator(symmetry,lattice,recip_lattice,symmetry_group, &
   & inverse_symmetry_id) result(this)
  implicit none
  
  type(BasicSymmetry), intent(in) :: symmetry
  type(RealMatrix),    intent(in) :: lattice
  type(RealMatrix),    intent(in) :: recip_lattice
  type(Group),         intent(in) :: symmetry_group
  integer,             intent(in) :: inverse_symmetry_id
  type(SymmetryOperator)          :: this
  
  this%id          = symmetry%id
  this%tensor      = symmetry%tensor
  this%translation = symmetry%translation
  this%atom_group  = symmetry%atom_group
  this%rvectors    = symmetry%rvectors
  
  this%cartesian_tensor = transpose(lattice) &
                      & * symmetry%tensor    &
                      & * recip_lattice
  this%cartesian_translation = transpose(lattice) &
                           & * symmetry%translation
  
  this%recip_tensor_ = transpose(invert(symmetry%tensor))
  
  this%symmetry_group = symmetry_group
  this%inverse_symmetry_id = inverse_symmetry_id
end function

! ----------------------------------------------------------------------
! Constructs a BasicSymmetry from a SymmetryOperator.
! ----------------------------------------------------------------------
impure elemental function new_BasicSymmetry_SymmetryOperator(this) &
   & result(output)
  implicit none
  
  type(SymmetryOperator), intent(in) :: this
  type(BasicSymmetry)                :: output
  
  output = BasicSymmetry( id          = this%id,          &
                        & tensor      = this%tensor,      &
                        & translation = this%translation, &
                        & atom_group  = this%atom_group,  &
                        & rvectors    = this%rvectors     )
end function

! ----------------------------------------------------------------------
! Returns whether or not two symmetry operators commute.
! ----------------------------------------------------------------------
function operators_commute(this,that,qpoint) result(output)
  implicit none
  
  type(SymmetryOperator), intent(in) :: this
  type(SymmetryOperator), intent(in) :: that
  type(QpointData),       intent(in) :: qpoint
  logical                            :: output
  
  type(IntVector) :: commutator_rvector
  
  integer         :: i
  
  ! Check that the symmetries' tensors commute.
  if (.not. matrices_commute(this%tensor,that%tensor)) then
    output = .false.
    return
  endif
  
  ! Check that the atom to atom transformations commute
  !    within a single unit cell.
  if (this%atom_group*that%atom_group /= that%atom_group*this%atom_group) then
    output = .false.
    return
  endif
  
  ! Check that R-vector translations commute.
  ! The condition on for two operators to commute (in addition to those above)
  !    is that the R-vector translation of their commutator is s.t. q.R=0.
  do i=1,size(this%atom_group)
    commutator_rvector = this%tensor*that%rvectors(i)     &
                     & + that%rvectors(this%atom_group*i) &
                     & - that%tensor*this%rvectors(i)     &
                     & - this%rvectors(that%atom_group*i)
    if (.not. is_int(qpoint%qpoint*commutator_rvector)) then
      output = .false.
      return
    endif
  enddo
  
  output = .true.
end function

! ----------------------------------------------------------------------
! Transforms a q-point, and translates it back to
!    the primitive reciprocal cell.
! ----------------------------------------------------------------------
impure elemental function transform_qpoint(this,qpoint) result(output)
  implicit none
  
  class(SymmetryOperator), intent(in) :: this
  type(QpointData),        intent(in) :: qpoint
  type(QpointData)                    :: output
  
  ! Transfer across all metadata.
  output = qpoint
  
  ! Transform q-point.
  output%qpoint = this%recip_tensor_ * qpoint%qpoint
  
  ! Translate q-point back into primitive reciprocal cell if necessary.
  call output%translate_to_primitive()
end function

impure elemental function inverse_transform(this,qpoint) result(output)
  implicit none
  
  class(SymmetryOperator), intent(in) :: this
  type(QpointData),        intent(in) :: qpoint
  type(QpointData)                    :: output
  
  ! Transfer across all metadata.
  output = qpoint
  
  ! Transform q-point.
  output%qpoint = transpose(this%tensor) * qpoint%qpoint
  
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
  type(IntMatrix)   :: tensor   ! The tensor after n operations.
  type(IntVector)   :: rvector  ! The translation after n operations.
  integer           :: atom_1   ! atom 1 maps to this after n operations.
  type(IntFraction) :: qr       ! qpoint*rvector.
  
  integer :: i
  
  identity = make_identity_matrix(3)
  
  ! After (S^i) : r -> tensor * r + rvector.
  ! After (S^0) : r -> r = I*r+0.
  tensor  = identity
  rvector = zeroes(3)
  atom_1  = 1
  
  ! Calculate n s.t. R^n=I. At this point, S^n=e^{2piiq.R}.
  output = 0
  do i=1,6
    ! After S^i, r -> tensor * r + rvector.
    tensor = this%tensor * tensor
    rvector = this%tensor * rvector + this%rvectors(atom_1)
    atom_1 = this%atom_group * atom_1
    if (tensor==identity .and. atom_1==1) then
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
