! ======================================================================
! A symmetry operator, in various representations.
! ======================================================================
module symmetry_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  use linear_algebra_module
  use fraction_algebra_module
  use group_module
  implicit none
  
  private
  
  public :: SymmetryOperator
  public :: calculate_symmetries
  
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
    type(FractionMatrix) :: recip_rotation
    
    ! The mapping from atoms to other atoms.
    ! rho_i and rho_j are the equilibrium positions of atom i and j.
    ! R * rho_i + T = rho_j + rvector_j
    !    - atom_group*i = j,
    !    - rvector(j) = rvector_j in fractional co-ordinates.
    type(Group)                  :: atom_group
    type(IntVector), allocatable :: rvector(:)
    
    ! The mapping from symmetry operators to other symmetry operators.
    ! If S * S_i = S_j then symmetry_group*i = j,
    !    where S is this operator, and S_i and S_j are other operators.
    type(Group) :: operator_group
    
    ! The id of the operator S_j, s.t. S*S_j = S_j*S = I.
    integer :: inverse
  end type
contains

! ----------------------------------------------------------------------
! Calculates symmetry operations, in various representations.
! ----------------------------------------------------------------------
function calculate_symmetries(basic_symmetries,lattice,recip_lattice,atoms) &
   & result(output)
  use atom_module
  use basic_symmetry_module
  implicit none
  
  type(BasicSymmetry), intent(in)     :: basic_symmetries(:)
  type(RealMatrix),    intent(in)     :: lattice
  type(RealMatrix),    intent(in)     :: recip_lattice
  type(AtomData),      intent(in)     :: atoms(:)
  type(SymmetryOperator), allocatable :: output(:)
  
  ! Atom group variables.
  type(RealVector)             :: transformed_position
  type(RealVector)             :: distance
  real(dp),        allocatable :: offsets(:)
  integer,         allocatable :: operations(:,:)
  type(IntVector), allocatable :: rvectors(:,:)
  
  ! Operator group variables.
  integer, allocatable :: operator_group(:)
  type(IntMatrix)      :: rotation_k
  type(Group)          :: operation_k
  
  ! Invers variables.
  integer :: identity
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  allocate(output(size(basic_symmetries)), stat=ialloc); call err(ialloc)
  
  ! Returns if the symmetry is only a dummy placeholder.
  if (size(output)==0) then
    return
  endif
  
  ! --------------------------------------------------
  ! Calculates basic symmetry properties.
  ! --------------------------------------------------
  do i=1,size(output)
    output(i)%rotation = basic_symmetries(i)%rotation
    output(i)%translation = basic_symmetries(i)%translation
    
    output(i)%cartesian_rotation = transpose(lattice) &
                               & * output(i)%rotation &
                               & * recip_lattice
    output(i)%cartesian_translation = transpose(lattice) &
                                  & * output(i)%translation
    
    output(i)%recip_rotation = transpose(invert(output(i)%rotation))
  enddo

  ! --------------------------------------------------
  ! Calculates how symmetries map atoms onto other atoms.
  ! --------------------------------------------------
  ! If symmetry i maps atoms(j)%frac_pos to atoms(k)%frac_pos + R then
  !    - output(i)%atom_group * j = k.
  !    - output(i)%rvector(j)     = R.
  
  ! Work out which atoms map to which atoms under each symmetry operation.
  allocate( offsets(size(atoms)),                  &
          & operations(size(atoms), size(output)), &
          & rvectors(size(atoms), size(output)),   &
          & stat=ialloc); call err(ialloc)
  do i=1,size(output)
    do j=1,size(atoms)
      ! Calculate the position of the transformed atom.
      transformed_position = output(i)%rotation             &
                         & * atoms(j)%fractional_position() &
                         & + output(i)%translation
      
      ! Identify which atom is closest to the transformed position,
      !    modulo supercell lattice vectors.
      do k=1,size(atoms)
        distance = transformed_position - atoms(k)%fractional_position()
        offsets(k) = l2_norm(distance - vec(nint(distance)))
      enddo
      operations(j,i) = minloc(offsets,1)
      
      ! R * position(i) + t = position(j) + an R-vector.
      ! Identify this R-vector, and record it.
      rvectors(j,i) = nint( transformed_position &
                        & - atoms(j)%fractional_position())
      
      ! Check that the transformed atom is acceptably close to its image.
      if (offsets(operations(j,i))>1.0e-10_dp) then
        call print_line(CODE_ERROR//': Error mapping atoms under symmetry.')
        call err()
      endif
    enddo
  enddo
  
  ! Check that each symmetry is one-to-one, and that mapped atoms are of the
  !    same species.
  do i=1,size(output)
    do j=1,size(atoms)
      if (count(operations(:,i)==j)/=1) then
        call print_line('Error: symmetry operation not one-to-one.')
        call err()
      endif
      
      if (atoms(operations(j,i))%species()/=atoms(j)%species()) then
        call print_line('Error: symmetry operation between different species.')
        call err()
      endif
    enddo
  enddo
  
  do i=1,size(output)
    output(i)%atom_group = operations(:,i)
    output(i)%rvector = rvectors(:,i)
  enddo

  ! --------------------------------------------------
  ! Calculates how symmetries map symmetries onto other symmetries.
  ! --------------------------------------------------
  ! output(i)%operator_group * j = k if symmetry i * symmetry j = symmetry k,
  !    modulo translations by lattice vectors.
  
  allocate(operator_group(size(output)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    do_j : do j=1,size(output)
      rotation_k = output(i)%rotation &
               & * output(j)%rotation
      operation_k = output(i)%atom_group &
                & * output(j)%atom_group
      do k=1,size(output)
        if ( rotation_k==output(k)%rotation .and. &
           & operation_k==output(k)%atom_group) then
          operator_group(j) = k
          cycle do_j
        endif
      enddo
      
      call print_line('Error: symmetry '//i//' times symmetry '//j//' is not &
         &itself a symmetry.')
      call err()
    enddo do_j
    
    output(i)%operator_group = operator_group
  enddo

  ! --------------------------------------------------
  ! Calculates the inverse of each symmetry.
  ! --------------------------------------------------
  ! If symmetry i * symmetry j = I then output(i)=j.
  
  ! Locate the identity operator. S*S=S iff S=I.
  identity = 0
  do i=1,size(output)
    if (output(i)%operator_group*i == i) then
      identity = i
      exit
    endif
  enddo
  
  if (identity==0) then
    call print_line('Error: The identity symmetry has not been found.')
    call err()
  endif
  
  ! Locate the inverse of each operator.
  do_i : do i=1,size(output)
    do j=1,size(output)
      if (output(i)%operator_group*j==identity) then
        output(i)%inverse = j
        
        if (j<i) then
          if (output(j)%inverse/=i) then
            call print_line('Error: operator inverses do not match.')
            call err()
          endif
        endif
        
        cycle do_i
      endif
    enddo
    
    call print_line('Error: operator inverse not found.')
    call err()
  enddo do_i
end function
end module
