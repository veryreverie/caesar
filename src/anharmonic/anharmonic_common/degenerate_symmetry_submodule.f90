submodule (caesar_degenerate_symmetry_module) caesar_degenerate_symmetry_submodule
  use caesar_anharmonic_common_module
contains

module procedure new_SingleModeSymmetry
  if (size(symmetric_mode_ids)/=size(symmetric_mode_coefficients)) then
    call err()
  endif
  
  this%mode_id = mode_id
  this%symmetric_mode_ids = symmetric_mode_ids
  this%symmetric_mode_coefficients = symmetric_mode_coefficients
end procedure

module procedure new_DegenerateSymmetry
  type(ComplexMatrix), allocatable :: subspace_symmetries(:)
  
  type(ComplexMode), allocatable :: degenerate_modes(:)
  type(QpointData),  allocatable :: degenerate_qpoints(:)
  
  type(ComplexMode)                     :: mode
  type(QpointData)                      :: qpoint
  type(QpointData)                      :: transformed_qpoint
  integer                               :: subspace_position
  type(DegenerateSubspace)              :: subspace
  type(ComplexMode),        allocatable :: symmetric_modes(:)
  integer                               :: mode_i_position
  integer,                  allocatable :: symmetric_mode_positions(:)
  complex(dp),              allocatable :: symmetric_mode_coefficients(:)
  complex(dp),              allocatable :: subspace_symmetry(:,:)
  
  integer :: i,j,ialloc
  
  this%symmetry_id = symmetry%id
  
  ! Calculate symmetries subspace by subspace.
  allocate(subspace_symmetries(size(subspaces)), stat=ialloc); call err(ialloc)
  do i=1,size(subspaces)
    degenerate_modes = subspaces(i)%modes(modes)
    degenerate_qpoints = subspaces(i)%qpoints(modes,qpoints)
    subspace_symmetries(i) =                                           &
       & calculate_symmetry_in_normal_coordinates( degenerate_modes,   &
       &                                           degenerate_qpoints, &
       &                                           symmetry            )
  enddo
  
  allocate(this%symmetries_(size(modes)), stat=ialloc); call err(ialloc)
  do i=1,size(this%symmetries_)
    mode = modes(i)
    
    ! Identify the q-point of mode i, and the q-point of the set of modes
    !    which the symmetry transfoms mode i to.
    qpoint = qpoints(first(qpoints%id==mode%qpoint_id))
    transformed_qpoint = symmetry*qpoint
    j = first(qpoints==transformed_qpoint, default=0)
    if (j/=0) then
      transformed_qpoint = qpoints(j)
    else
      call print_line(ERROR//': Unable to find symmetry-transformed q-point. &
         &This may be caused by the q-point grid not obeying the symmetry &
         &of the system.')
      call err()
    endif
    
    ! Identify the subspace to which mode i belongs, and its position in
    !    the array subspaces.
    subspace_position = first(subspaces%id==modes(i)%subspace_id)
    subspace = subspaces(subspace_position)
    
    ! Identify the set of modes which the symmetry transforms mode i to.
    ! These modes must be degenerate with mode i, and at the q-point which
    !    the symmetry transforms mode i's q-point to.
    symmetric_modes = subspace%modes(modes)
    symmetric_modes = symmetric_modes(                            &
       & filter(symmetric_modes%qpoint_id==transformed_qpoint%id) )
    
    ! Identify the position of mode i within the subspace's mode_ids,
    !    and the positions of the symmetric modes within the same.
    mode_i_position = first(subspace%mode_ids==mode%id)
    symmetric_mode_positions =                               &
       & [( first(subspace%mode_ids==symmetric_modes(j)%id), &
       & j=1,                                                &
       & size(symmetric_modes)                               )]
    
    ! Identify the coefficients of the symmetry.
    ! If mode i is ui, and symmetric_modes are {uj}, then these are the
    !    coefficients of S.ui in terms of uj.
    subspace_symmetry = cmplx(subspace_symmetries(subspace_position))
    symmetric_mode_coefficients = subspace_symmetry( &
                        & symmetric_mode_positions,  &
                        & mode_i_position            )
    
    ! Construct the symmetry.
    this%symmetries_(i) = SingleModeSymmetry( &
                & mode%id,                    &
                & symmetric_modes%id,         &
                & symmetric_mode_coefficients )
    
    if (abs(l2_norm(symmetric_mode_coefficients)-1)>1e-6_dp) then
      call print_line(CODE_ERROR//': Error constructing symmetry '// &
         & symmetry%id//' in normal mode co-ordinates.')
      call print_line(this%symmetries_(i))
      call err()
    endif
  enddo
end procedure

module procedure calculate_symmetry
  complex(dp), allocatable :: symmetry(:,:)
  
  integer :: i,j,k,ialloc
  
  type(ComplexPolynomial) :: transformed_input
  
  if (size(input)==0) then
    output = cmplxmat(zeroes(0,0))
    return
  endif
  
  allocate(symmetry(size(input), size(input)), stat=ialloc); call err(ialloc)
  symmetry = 0
  do i=1,size(input)
    transformed_input = this%transform_monomial(input(i), modes)
    
    do j=1,size(transformed_input%terms)
      k = first_equivalent( input,                      &
                          & transformed_input%terms(j), &
                          & compare_complex_monomials   )
      symmetry(k,i) = transformed_input%terms(j)%coefficient
      
      if (include_coefficients) then
        symmetry(k,i) = symmetry(k,i) &
                    & * input(i)%coefficient &
                    & / input(k)%coefficient
      endif
    enddo
  enddo
  
  output = mat(symmetry)
end procedure

module procedure transform_monomial
  integer,                  allocatable :: mode_ids(:)
  type(SingleModeSymmetry), allocatable :: symmetry
  type(ComplexMode)                     :: symmetric_mode
  type(ComplexMonomial),    allocatable :: monomials(:)
  
  integer :: i,j,k,ialloc
  
  ! Convert the monomial into a list of modes.
  allocate(mode_ids(input%total_power()), stat=ialloc); call err(ialloc)
  k = 0
  do i=1,size(input)
    do j=1,input%power(i)
      k = k+1
      mode_ids(k) = input%id(i)
    enddo
    
    if (input%id(i)/=input%paired_id(i)) then
      do j=1,input%paired_power(i)
        k = k+1
        mode_ids(k) = input%paired_id(i)
      enddo
    endif
  enddo
  
  do i=1,size(mode_ids)
    symmetry = this%symmetries_(first(this%symmetries_%mode_id==mode_ids(i)))
    allocate( monomials(size(symmetry%symmetric_mode_ids)), &
            & stat=ialloc); call err(ialloc)
    do j=1,size(symmetry%symmetric_mode_ids)
      symmetric_mode = modes(first(modes%id==symmetry%symmetric_mode_ids(j)))
      monomials(j) = ComplexMonomial(                                 &
         & coefficient = symmetry%symmetric_mode_coefficients(j),     &
         & modes       = [ComplexUnivariate(symmetric_mode, power=1)] )
    enddo
    
    if (i==1) then
      output = ComplexPolynomial(monomials)
    else
      output%terms = [(                                          &
         & (output%terms(i)*monomials(j), j=1, size(monomials)), &
         & i=1,                                                  &
         & size(output)                                          )]
    endif
    deallocate(monomials, stat=ialloc); call err(ialloc)
  enddo
  
  call output%simplify()
end procedure

module procedure read_SingleModeSymmetry
  type(String), allocatable :: line(:)
  type(String), allocatable :: symmetric_mode(:)
  
  integer                  :: mode_id
  integer,     allocatable :: symmetric_mode_ids(:)
  complex(dp), allocatable :: symmetric_mode_coefficients(:)
  
  integer :: i,ialloc
  
  select type(this); type is(SingleModeSymmetry)
    ! If the symmetry takes u1 to a*u2+b*u3+... then
    ! Input = 'u1 -> au2 + bu3 + ...'
    ! line = ['u1', '->', 'au2', '+', 'bu3', ...]
    line = split_line(input)
    
    ! Split off the ID of the input mode. (trim the 'u' off the front.)
    mode_id = int(slice(line(1),2,len(line(1))))
    
    ! Split each 'au2', 'bu3' etc. term by the 'u', to leave the coefficient
    !    and mode id of each.
    allocate( symmetric_mode_ids((size(line)-1)/2),          &
            & symmetric_mode_coefficients((size(line)-1)/2), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(symmetric_mode_ids)
      symmetric_mode = split_line(line(2*i+1), delimiter='u')
      symmetric_mode_ids(i) = int(symmetric_mode(2))
      symmetric_mode_coefficients(i) = cmplx(symmetric_mode(1))
    enddo
    
    this = SingleModeSymmetry( mode_id,                    &
                             & symmetric_mode_ids,         &
                             & symmetric_mode_coefficients )
  class default
    call err()
  end select
end procedure

module procedure write_SingleModeSymmetry
  type(String), allocatable :: symmetric_modes(:)
  
  integer :: i,ialloc
  
  select type(this); type is(SingleModeSymmetry)
    allocate( symmetric_modes(size(this%symmetric_mode_ids)), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(symmetric_modes)
      symmetric_modes(i) = this%symmetric_mode_coefficients(i)// &
                         & 'u'//this%symmetric_mode_ids(i)
    enddo
    output = 'u'//this%mode_id//' -> '//join(symmetric_modes,delimiter=' + ')
  class default
    call err()
  end select
end procedure

module procedure new_SingleModeSymmetry_String
  call this%read(input)
end procedure

module procedure read_DegenerateSymmetry
  type(String), allocatable :: line(:)
  
  select type(this); type is(DegenerateSymmetry)
    line = split_line(input(1))
    this%symmetry_id = int(line(3))
    this%symmetries_ = SingleModeSymmetry(input(2:))
  class default
    call err()
  end select
end procedure

module procedure write_DegenerateSymmetry
  select type(this); type is(DegenerateSymmetry)
    output = [ 'Symmetry : '//this%symmetry_id, &
             & str(this%symmetries_)            ]
  class default
    call err()
  end select
end procedure

module procedure new_DegenerateSymmetry_Strings
  call this%read(input)
end procedure

module procedure new_DegenerateSymmetry_StringArray
  this = DegenerateSymmetry(str(input))
end procedure
end submodule
