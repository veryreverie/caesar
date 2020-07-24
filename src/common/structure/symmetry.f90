! ======================================================================
! A symmetry operator, in various representations.
! ======================================================================
module symmetry_module
  use utils_module
  
  use atom_module
  
  use basic_symmetry_module
  use qpoint_module
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
impure elemental function new_SymmetryOperator(symmetry,lattice, &
   & recip_lattice) result(this)
  implicit none
  
  type(BasicSymmetry), intent(in) :: symmetry
  type(RealMatrix),    intent(in) :: lattice
  type(RealMatrix),    intent(in) :: recip_lattice
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
! Returns whether or not two symmetry operators commute or anti-commute.
! ----------------------------------------------------------------------
! Each operator is a tensor product of the cartesian tensor and an atom tensor.
! Two operators commute if:
!    - Their cartesian and atom tensors both commute.
!    - Their cartesian and atom tensors both anti-commute.
! Two operators anti-commute if:
!    - One tensor commutes and the other anticommutes.
impure elemental function operators_commute(a,b,qpoint) result(output)
  implicit none
  
  type(SymmetryOperator), intent(in) :: a
  type(SymmetryOperator), intent(in) :: b
  type(QpointData),       intent(in) :: qpoint
  logical                            :: output
  
  output = .false.
  
  ! Check that the atom groups commute.
  ! This is a necessary condition for the atom tensors commuting
  !    or anti-commuting.
  if (a%atom_group*b%atom_group /= b%atom_group*a%atom_group) then
    return
  endif
  
  ! Check that either the cartesian and atom tensors both commute
  !    or both anticommute.
  if (matrices_commute(a%tensor,b%tensor)) then
    if (atom_tensors_commute(a,b,qpoint)) then
      output = .true.
    endif
  elseif (matrices_anticommute(a%tensor,b%tensor)) then
    if (atom_tensors_anticommute(a,b,qpoint)) then
      output = .true.
    endif
  endif
end function

impure elemental function operators_anticommute(a,b,qpoint) &
   & result(output)
  implicit none
  
  type(SymmetryOperator), intent(in) :: a
  type(SymmetryOperator), intent(in) :: b
  type(QpointData),       intent(in) :: qpoint
  logical                            :: output
  
  output = .false.
  
  ! Check that the atom groups commute.
  ! This is a necessary condition for the atom tensors commuting
  !    or anti-commuting.
  if (a%atom_group*b%atom_group /= b%atom_group*a%atom_group) then
    return
  endif
  
  ! Check that one of the cartesian and atom tensors commute
  !    and the other anticommutes.
  if (matrices_commute(a%tensor,b%tensor)) then
    if (atom_tensors_anticommute(a,b,qpoint)) then
      output = .true.
    endif
  elseif (matrices_anticommute(a%tensor,b%tensor)) then
    if (atom_tensors_commute(a,b,qpoint)) then
      output = .true.
    endif
  endif
end function

! For two symmetries A and B: A1=A,A2=A*,B1=B,B2=B*.
!    if positive_superposition=.true.  : returns whether [A1+A2,B1+B2]=0.
!    if positive_superposition=.false. : returns whether [A1-A2,B1-B2]=0.
! Each symmetry is the tensor product of the cartesian tensor T and the
!    atom tensor G: A1=TA1^GA1, A2=TA2^GA2, B1=TB1^GB1 and B2=TB2^GB2.
! [A1+A2,B1+B2] = A1B1+A1B2+A2B1+A2B2-B1A1-B1A2-B2A1-B2A2
! [A1-A2,B1-B2] = A1B1-A1B2-A2B1+A2B2-B1A1+B1A2+B2A1-B2A2
impure elemental function superposed_operators_commute(a,b,qpoint, &
   & positive_superposition) result(output)
  implicit none
  
  type(SymmetryOperator), intent(in) :: a
  type(SymmetryOperator), intent(in) :: b
  type(QpointData),       intent(in) :: qpoint
  logical,                intent(in) :: positive_superposition
  logical                            :: output
  
  type(IntMatrix)              :: a2_tensor
  type(Group)                  :: a2_atom_group
  type(IntVector), allocatable :: a2_rvectors(:)
  
  type(IntMatrix)              :: b2_tensor
  type(Group)                  :: b2_atom_group
  type(IntVector), allocatable :: b2_rvectors(:)
  
  type(Group)          :: groups(8)
  type(IntMatrix)      :: tensors(8)
  integer              :: signatures(8)
  type(IntFraction)    :: qrs(8)
  logical              :: matched(8)
  integer, allocatable :: matches(:)
  
  type(IntMatrix) :: tensor
  
  integer :: i,j,k,ialloc
  
  output = .true.
  
  ! Construct A2=A* from A, and B2=B* from B.
  a2_tensor = intmat(invert(a%tensor))
  a2_atom_group = a%atom_group%inverse()
  allocate(a2_rvectors(size(a%rvectors)), stat=ialloc); call err(ialloc)
  do i=1,size(a%rvectors)
    a2_rvectors(a%atom_group*i) = -a%rvectors(i)
  enddo
  
  b2_tensor = intmat(invert(a%tensor))
  b2_atom_group = a%atom_group%inverse()
  allocate(b2_rvectors(size(a%rvectors)), stat=ialloc); call err(ialloc)
  do i=1,size(a%rvectors)
    b2_rvectors(a%atom_group*i) = -a%rvectors(i)
  enddo
  
  ! Construct the atom groups for the four AB operators,
  !    and for the four BA operators.
  groups = [ a%atom_group*b%atom_group,   &
           & a%atom_group*b2_atom_group,  &
           & a2_atom_group*b%atom_group,  &
           & a2_atom_group*b2_atom_group, &
           & b%atom_group*a%atom_group,   &
           & b%atom_group*a2_atom_group,  &
           & b2_atom_group*a%atom_group,  &
           & b2_atom_group*a2_atom_group  ]
  
  ! Do the same for the cartesian tensors.
  tensors = [ a%tensor*b%tensor,   &
            & a%tensor*b2_tensor,  &
            & a2_tensor*b%tensor,  &
            & a2_tensor*b2_tensor, &
            & b%tensor*a%tensor,   &
            & b%tensor*a2_tensor,  &
            & b2_tensor*a%tensor,  &
            & b2_tensor*a2_tensor  ]
  
  ! Construct the sign signatures.
  if (positive_superposition) then
    signatures = [ 1,  1,  1,  1, -1, -1, -1, -1]
  else
    signatures = [ 1, -1, -1,  1, -1,  1,  1, -1]
  endif
  
  ! Loop over atoms.
  ! For each atom, identify which operator combinations add together.
  ! For two operator combinations to add, they must have matching atom groups
  !    and a relative phase which is either +1 or -1.
  ! Check that the tensor from each matching set is zero.
  do i=1,size(a%atom_group)
    qrs = qpoint%qpoint * [                                       &
       & a%tensor*b%rvectors(i) + a%rvectors(b%atom_group*i),     &
       & a%tensor*b2_rvectors(i) + a%rvectors(b2_atom_group*i),   &
       & a2_tensor*b%rvectors(i) + a2_rvectors(b%atom_group*i),   &
       & a2_tensor*b2_rvectors(i) + a2_rvectors(b2_atom_group*i), &
       & b%tensor*a%rvectors(i) + b%rvectors(a%atom_group*i),     &
       & b%tensor*a2_rvectors(i) + b%rvectors(a2_atom_group*i),   &
       & b2_tensor*a%rvectors(i) + b2_rvectors(a%atom_group*i),   &
       & b2_tensor*a2_rvectors(i) + b2_rvectors(a2_atom_group*i)  ]
    
    matched = .false.
    do j=1,8
      if (matched(j)) then
        cycle
      endif
      
      matches = filter(       groups*i==groups(j)*i  &
                      & .and. is_int(2*(qrs-qrs(j))) )
      
      tensor = zeroes(3,3)
      
      do k=1,size(matches)
        if (is_int(qrs(matches(k))-qrs(j))) then
          ! e^2pii(q.r(matches(k))-q.r(j)) is 1.
          tensor = tensor+tensors(matches(k))*signatures(matches(k))
        else
          ! e^2pii(q.r(matches(k))-q.r(j)) is -1.
          tensor = tensor-tensors(matches(k))*signatures(matches(k))
        endif
      enddo
      
      if (tensor/=zeroes(3,3)) then
        output = .false.
        return
      endif
      
      matched(matches) = .true.
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Helper functions for operators_commute and operators_anticommute,
!    returning whether or not the operators' R-vectors (anti-)commute.
! ----------------------------------------------------------------------
! If the primitive cell has N atoms then the atom tensor is an N*N matrix,
!    with N non-zero elements corresponding to each r_i -> r_j + R_i,
!    so the ij element is e^(2*pi*i*q*R_i).
impure elemental function atom_tensors_commute(a,b,qpoint) &
   & result(output)
  implicit none
  
  type(SymmetryOperator), intent(in) :: a
  type(SymmetryOperator), intent(in) :: b
  type(QpointData),       intent(in) :: qpoint
  logical                            :: output
  
  type(IntVector)   :: r_ab
  type(IntVector)   :: r_ba
  type(IntFraction) :: qr
  integer           :: i
  
  output = .true.
  
  do i=1,size(a%atom_group)
    r_ab = a%tensor*b%rvectors(i) + a%rvectors(b%atom_group*i)
    r_ba = b%tensor*a%rvectors(i) + b%rvectors(a%atom_group*i)
    ! For e^(2pii.q.r_ab) to equal e^(2pii.q.r_ba),
    !    q.(r_ab-r_ba) must be an integer.
    qr = qpoint%qpoint*(r_ab-r_ba)
    if (.not. is_int(qr)) then
      output = .false.
      return
    endif
  enddo
end function

impure elemental function atom_tensors_anticommute(a,b,qpoint) &
   & result(output)
  implicit none
  
  type(SymmetryOperator), intent(in) :: a
  type(SymmetryOperator), intent(in) :: b
  type(QpointData),       intent(in) :: qpoint
  logical                            :: output
  
  type(IntVector)   :: r_ab
  type(IntVector)   :: r_ba
  type(IntFraction) :: qr
  integer           :: i
  
  output = .true.
  
  do i=1,size(a%atom_group)
    r_ab = a%tensor*b%rvectors(i) + a%rvectors(b%atom_group*i)
    r_ba = b%tensor*a%rvectors(i) + b%rvectors(a%atom_group*i)
    ! For e^(2pii.q.r_ab) to equal -e^(2pii.q.r_ba),
    !    q.(r_ab-r_ba) must be an integer plus 1/2.
    qr = qpoint%qpoint*(r_ab-r_ba)
    if (.not. (is_int(2*qr) .and. .not. is_int(qr))) then
      output = .false.
      return
    endif
  enddo
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
impure elemental function symmetry_order(this,qpoint) result(output)
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

! ----------------------------------------------------------------------
! Returns whether or not the LO/TO direction breaks a given symmetry.
! ----------------------------------------------------------------------
! A symmetry is broken if its tensor transforms the LO/TO direction
!    onto a different axis.
! Symmetries which map the LO/TO direction onto minus itself are allowed.
impure elemental function loto_breaks_symmetry(tensor,loto_direction) &
   & result(output)
  implicit none
  
  type(IntMatrix),      intent(in) :: tensor
  type(FractionVector), intent(in) :: loto_direction
  logical                          :: output
  
  type(FractionVector) :: transformed_direction
  
  ! A symmetry with tensor S acts on the LO/TO direction q as
  !    S: q -> S^(-T).q
  ! So the LO/TO direction breaks the symmetry if  q /= +/- S^(-T).q,
  !    or equivalently if q.S /= +/- q
  transformed_direction = loto_direction*tensor
  output = transformed_direction /=  loto_direction .and. &
         & transformed_direction /= -loto_direction
end function
end module
