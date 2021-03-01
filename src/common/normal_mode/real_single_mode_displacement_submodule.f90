submodule (caesar_real_single_mode_displacement_module) caesar_real_single_mode_displacement_submodule
  use caesar_normal_mode_module
contains

module procedure new_RealSingleDisplacement
  this%id        = id
  this%magnitude = magnitude
end procedure

module procedure new_RealSingleDisplacement_RealMode
  this = RealSingleDisplacement(mode%id, magnitude)
end procedure

module procedure multiply_real_RealSingleDisplacement
  output = RealSingleDisplacement( id        = that%id,            &
                                 & magnitude = this*that%magnitude )
end procedure

module procedure multiply_RealSingleDisplacement_real
  output = RealSingleDisplacement( id        = this%id,            &
                                 & magnitude = this%magnitude*that )
end procedure

module procedure divide_RealSingleDisplacement_real
  output = RealSingleDisplacement( id        = this%id,            &
                                 & magnitude = this%magnitude/that )
end procedure

module procedure add_RealSingleDisplacement_RealSingleDisplacement
  if (this%id/=that%id) then
    call print_line(ERROR//': Trying to add vectors along different modes.')
    call err()
  endif
  
  output = RealSingleDisplacement( id        = this%id,                      &
                                 & magnitude = this%magnitude+that%magnitude )
end procedure

module procedure negative_RealSingleDisplacement
  output = RealSingleDisplacement(id=this%id, magnitude=-this%magnitude)
end procedure

module procedure subtract_RealSingleDisplacement_RealSingleDisplacement
  if (this%id/=that%id) then
    call print_line(ERROR//': Trying to subtract vectors along different &
       &modes.')
    call err()
  endif
  
  output = RealSingleDisplacement( id        = this%id,                      &
                                 & magnitude = this%magnitude-that%magnitude )
end procedure

module procedure new_MassWeightedDisplacement_RealSingleDisplacement
  if (mode%id/=this%id) then
    call print_line(CODE_ERROR//': Mode and vector incompatible.')
    call err()
  elseif ( qpoint%id/=mode%qpoint_id_plus .and. &
         & qpoint%id/=mode%qpoint_id_minus      ) then
    call print_line(CODE_ERROR//': Mode and q-point incompatible.')
    call err()
  endif
  
  output = this%magnitude * MassWeightedDisplacement(mode,structure,qpoint)
end procedure

module procedure new_CartesianDisplacement_RealSingleDisplacement
  if (mode%id/=this%id) then
    call print_line(CODE_ERROR//': Mode and vector incompatible.')
    call err()
  elseif ( qpoint%id/=mode%qpoint_id_plus .and. &
         & qpoint%id/=mode%qpoint_id_minus      ) then
    call print_line(CODE_ERROR//': Mode and q-point incompatible.')
    call err()
  endif
  
  output = CartesianDisplacement(                            &
     & MassWeightedDisplacement(this,mode,structure,qpoint), &
     & structure                                             )
end procedure

module procedure new_RealSingleDisplacement_MassWeightedDisplacement
  type(MassWeightedDisplacement) :: mode_vector
  
  real(dp) :: magnitude
  
  ! Modes are orthonormal in mass-reduced co-ordinates.
  ! If M is the mass-weighting matrix, M_ab = 1/sqrt(m_a*m_b),
  !    where m_a and m_b are the masses of atoms a and b respectively, then
  !
  ! u_i and u_j are the cartesian representations of modes i and j.
  ! r is the cartesian vector.
  !
  ! u_i.M.u_j = 0 if i/=j
  !           = n if i=j, where n is the number of primitive cells.
  !
  ! r = sum_j[ a_j*u_j ]
  ! => u_i.M.r = sum_j[ a_j*u_i.M.u_j ] = a_i*u_i.M.u_i = a_i*n
  ! => a_j = u_i.M.r / n
  
  mode_vector = MassWeightedDisplacement(mode,structure,qpoint)
  
  magnitude = sum( vector%vectors       &
          &      * mode_vector%vectors) &
          & / structure%sc_size
  
  this = RealSingleDisplacement(id=mode%id, magnitude=magnitude)
end procedure

module procedure new_RealSingleDisplacement_CartesianDisplacement
  this = RealSingleDisplacement( mode,                                       &
                               & MassWeightedDisplacement(vector,structure), &
                               & structure,                                  &
                               & qpoint                                      )
end procedure

module procedure select_mode_RealSingleDisplacement
  output = modes(first(modes%id==displacement%id))
end procedure

module procedure select_modes_RealSingleDisplacements
  integer :: i,ialloc
  
  allocate(output(size(displacements)), stat=ialloc); call err(ialloc)
  do i=1,size(displacements)
    output(i) = select_mode(displacements(i), modes)
  enddo
end procedure

module procedure read_RealSingleDisplacement
  type(String), allocatable :: split_string(:)
  integer                   :: id
  real(dp)                  :: magnitude
  
  select type(this); type is(RealSingleDisplacement)
    split_string = split_line(input)
    if (size(split_string)/=3) then
      call print_line(ERROR//': unable to parse real single mode vector &
         &from string: '//input)
      call err()
    endif
    
    ! If e.g. id=3 and power=2.1 then split_string = ["u3","=","2.1"]
    ! The 'u' needs stripping off the first element to give the id.
    id = int(slice(split_string(1),2,len(split_string(1))))
    magnitude = dble(split_string(3))
    
    this = RealSingleDisplacement(id,magnitude)
  end select
end procedure

module procedure write_RealSingleDisplacement
  select type(this); type is(RealSingleDisplacement)
    output = 'u'//this%id//' = '//this%magnitude
  end select
end procedure

module procedure new_RealSingleDisplacement_String
  call this%read(input)
end procedure
end submodule
