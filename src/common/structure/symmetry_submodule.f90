submodule (caesar_symmetry_module) caesar_symmetry_submodule
  use caesar_structure_module
contains

module procedure new_SymmetryOperator
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
end procedure

module procedure new_BasicSymmetry_SymmetryOperator
  output = BasicSymmetry( id          = this%id,          &
                        & tensor      = this%tensor,      &
                        & translation = this%translation, &
                        & atom_group  = this%atom_group,  &
                        & rvectors    = this%rvectors     )
end procedure

module procedure operators_commute
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
end procedure

module procedure operators_anticommute
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
end procedure

module procedure superposed_operators_commute
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
end procedure

module procedure atom_tensors_commute
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
end procedure

module procedure atom_tensors_anticommute
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
end procedure

module procedure transform_qpoint
  ! Transfer across all metadata.
  output = qpoint
  
  ! Transform q-point.
  output%qpoint = this%recip_tensor_ * qpoint%qpoint
  
  ! Translate q-point back into primitive reciprocal cell if necessary.
  call output%translate_to_primitive()
end procedure

module procedure inverse_transform
  ! Transfer across all metadata.
  output = qpoint
  
  ! Transform q-point.
  output%qpoint = transpose(this%tensor) * qpoint%qpoint
  
  ! Translate q-point back into primitive reciprocal cell if necessary.
  call output%translate_to_primitive()
end procedure

module procedure symmetry_order
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
end procedure

module procedure loto_breaks_symmetry
  type(FractionVector) :: transformed_direction
  
  ! A symmetry with tensor S acts on the LO/TO direction q as
  !    S: q -> S^(-T).q
  ! So the LO/TO direction breaks the symmetry if  q /= +/- S^(-T).q,
  !    or equivalently if q.S /= +/- q
  transformed_direction = loto_direction*tensor
  output = transformed_direction /=  loto_direction .and. &
         & transformed_direction /= -loto_direction
end procedure
end submodule
