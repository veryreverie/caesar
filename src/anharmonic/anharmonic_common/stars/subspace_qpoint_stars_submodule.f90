submodule (caesar_subspace_qpoint_stars_module) &
   & caesar_subspace_qpoint_stars_submodule
  use caesar_stars_module
contains

module procedure new_SubspaceQpointStars
  this%subspace_id = subspace_id
  this%powers      = powers
end procedure

module procedure read_SubspaceQpointStars
  integer                        :: subspace_id
  type(QpointStars), allocatable :: powers(:)
  
  select type(this); type is(SubspaceQpointStars)
    subspace_id = int(token(input(1), 5))
    powers = QpointStars(split_into_sections( input(2:),                     &
                                            & separating_line=repeat('-',50) ))
    this = SubspaceQpointStars(subspace_id, powers)
  class default
    call err()
  end select
end procedure

module procedure write_SubspaceQpointStars
  select type(this); type is(SubspaceQpointStars)
    output = [ 'q-point stars in subspace '//this%subspace_id//' :', &
             & str(this%powers, separating_line=repeat('-',50))      ]
  class default
    call err()
  end select
end procedure

module procedure new_SubspaceQpointStars_Strings
  call this%read(input)
end procedure

module procedure new_SubspaceQpointStars_StringArray
  call this%read(str(input))
end procedure

module procedure generate_subspace_qpoint_stars
  type(QpointData), allocatable :: subspace_qpoints(:)
  
  integer, allocatable :: first_qpoint_ids(:)
  
  integer :: i,j,ialloc
  
  allocate( first_qpoint_ids(size(subspaces)), &
          & output(size(subspaces)),           &
          & stat=ialloc); call err(ialloc)
  output%subspace_id = subspaces%id
  do i=1,size(subspaces)
    ! Get the set of q-points corresponding to subspace i.
    subspace_qpoints = subspaces(i)%qpoints(modes, qpoints)
    subspace_qpoints = subspace_qpoints(set(subspace_qpoints%id))
    
    ! Check if subspace i involves the same q-points as a previously
    !    processed subspace.
    ! If so, simply copy the previous star.
    first_qpoint_ids(i) = subspace_qpoints(1)%id
    j = first(first_qpoint_ids(:i)==first_qpoint_ids(i))
    if (j<i) then
      output(i)%powers = output(j)%powers
      cycle
    endif
    
    ! Construct the q-point stars for subspace i.
    output(i)%powers = generate_qpoint_stars( subspace_qpoints,       &
                                            & qpoint_symmetry_groups, &
                                            & max_power,              &
                                            & max_qpoint_coupling,    &
                                            & conserve_momentum       )
  enddo
end procedure
end submodule
