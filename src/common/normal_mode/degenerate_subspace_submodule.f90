submodule (caesar_degenerate_subspace_module) caesar_degenerate_subspace_submodule
  use caesar_normal_mode_module
contains

module procedure new_DegenerateSubspace
  if (size(mode_ids)/=size(paired_ids)) then
    call print_line(CODE_ERROR//': mode IDs and paired IDs do not match.')
    call err()
  endif
  
  this%id         = id
  this%frequency  = frequency
  this%mode_ids   = mode_ids
  this%paired_ids = paired_ids
end procedure

module procedure size_DegenerateSubspace
  output = size(input%mode_ids)
end procedure

module procedure process_degeneracies
  integer,           allocatable :: subspace_ids(:)
  type(ComplexMode), allocatable :: subspace_modes(:)
  
  integer :: i,ialloc
  
  ! Make a list of degeneracy ids.
  subspace_ids = modes%subspace_id
  ! Remove the purely translational modes.
  subspace_ids = subspace_ids(filter(.not.modes%translational_mode))
  ! De-duplicate the list, so that each id appears exactly once.
  subspace_ids = subspace_ids(set(subspace_ids))
  
  ! Generate subspaces.
  allocate(output(size(subspace_ids)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    subspace_modes = modes(filter(modes%subspace_id==subspace_ids(i)))
    
    output(i) = DegenerateSubspace( id         = subspace_ids(i),             &
                                  & frequency  = subspace_modes(1)%frequency, &
                                  & mode_ids   = subspace_modes%id,           &
                                  & paired_ids = subspace_modes%paired_id     )
  enddo
end procedure

module procedure modes_DegenerateSubspace_ComplexModes
  integer :: i,j,ialloc
  
  !output = [( modes(first(modes%id==this%mode_ids(i))), i=1, size(this) )]
  ! WORKAROUND: To avoid a memory leak in ifort 19.1.0.166
  allocate(output(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    j = first(modes%id==this%mode_ids(i))
    output(i) = modes(j)
  enddo
end procedure

module procedure modes_DegenerateSubspace_RealModes
  integer :: i,j,ialloc
  
  !output = [( modes(first(modes%id==this%mode_ids(i))), i=1, size(this) )]
  ! WORKAROUND: To avoid a memory leak in ifort 19.1.0.166
  allocate(output(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    j = first(modes%id==this%mode_ids(i))
    output(i) = modes(j)
  enddo
end procedure

module procedure qpoints_DegenerateSubspace_ComplexModes
  type(ComplexMode), allocatable :: subspace_modes(:)
  
  integer :: i
  
  subspace_modes = this%modes(modes)
  output = [( qpoints(first(qpoints%id==subspace_modes(i)%qpoint_id)), &
            & i=1,                                                     &
            & size(subspace_modes)                                     )]
end procedure

module procedure qpoints_DegenerateSubspace_RealModes
  type(RealMode), allocatable :: subspace_modes(:)
  
  integer :: i
  
  subspace_modes = this%modes(modes)
  output = [( qpoints(first(qpoints%id==subspace_modes(i)%qpoint_id_plus)), &
            & i=1,                                                          &
            & size(subspace_modes)                                          )]
end procedure

module procedure read_DegenerateSubspace
  integer              :: id
  real(dp)             :: frequency
  integer, allocatable :: mode_ids(:)
  integer, allocatable :: paired_ids(:)
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(DegenerateSubspace)
    line = split_line(input(1))
    id = int(line(5))
    
    line = split_line(input(2))
    frequency = dble(line(5))
    
    line = split_line(input(3))
    mode_ids = int(line(5:))
    
    line = split_line(input(4))
    paired_ids = int(line(5:))
    
    this = DegenerateSubspace(id, frequency, mode_ids, paired_ids)
  class default
    call err()
  end select
end procedure

module procedure write_DegenerateSubspace
  select type(this); type is(DegenerateSubspace)
    output = [ 'Degenerate subspace ID : '//this%id,        &
             & 'Frequency of modes     : '//this%frequency, &
             & 'Degenerate mode IDs    : '//this%mode_ids,  &
             & 'Paired mode IDS        : '//this%paired_ids ]
  class default
    call err()
  end select
end procedure

module procedure new_DegenerateSubspace_Strings
  call this%read(input)
end procedure

module procedure new_DegenerateSubspace_StringArray
  this = DegenerateSubspace(str(input))
end procedure
end submodule
