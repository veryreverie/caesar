submodule (caesar_real_mode_module) caesar_real_mode_submodule
  use caesar_normal_mode_module
contains

module procedure new_RealMode
  this%id                 = id
  this%paired_id          = paired_id
  this%frequency          = frequency
  this%spring_constant    = spring_constant
  this%soft_mode          = soft_mode
  this%translational_mode = translational_mode
  this%cos_vector         = cos_vector
  this%sin_vector         = sin_vector
  this%qpoint_id_plus     = qpoint_id_plus
  this%qpoint_id_minus    = qpoint_id_minus
  this%subspace_id        = subspace_id
end procedure

module procedure new_MassWeightedDisplacement_RealMode
  output = MassWeightedDisplacement(this%construct_vector(structure,qpoint))
end procedure

module procedure new_MassWeightedForce_RealMode
  output = MassWeightedForce(this%construct_vector(structure,qpoint))
end procedure

module procedure new_CartesianDisplacement_RealMode
  output = CartesianDisplacement( MassWeightedDisplacement( this,         &
                                &                           structure,    &
                                &                           qpoint     ), &
                                & structure                               )
end procedure

module procedure new_CartesianForce_RealMode
  output = CartesianForce( MassWeightedForce( this,         &
                         &                    structure,    &
                         &                    qpoint     ), &
                         & structure                        )
end procedure

module procedure construct_vector
  type(AtomData)    :: atom
  type(IntFraction) :: qr
  
  integer :: i,ialloc
  
  if (qpoint%id/=this%qpoint_id_plus.and.qpoint%id/=this%qpoint_id_minus) then
    call print_line(CODE_ERROR//': Mode and q-point incompatible.')
    call err()
  endif
  
  allocate(output(structure%no_atoms), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    atom = structure%atoms(i)
    
    if (qpoint%id==this%qpoint_id_plus) then
      qr =  qpoint%qpoint*structure%rvectors(atom%rvec_id())
    else
      qr = -qpoint%qpoint*structure%rvectors(atom%rvec_id())
    endif
    
    output(i) = this%cos_vector(atom%prim_id()) * cos_2pi(qr) &
            & + this%sin_vector(atom%prim_id()) * sin_2pi(qr)
  enddo
end procedure

module procedure select_qpoint_RealMode
  output = qpoints(first(qpoints%id==mode%qpoint_id_plus))
end procedure

module procedure select_qpoints_RealModes
  integer :: i,ialloc
  
  allocate(output(size(modes)), stat=ialloc); call err(ialloc)
  do i=1,size(modes)
    output(i) = select_qpoint(modes(i), qpoints)
  enddo
end procedure

module procedure read_RealMode
  integer                       :: id
  integer                       :: paired_id
  real(dp)                      :: frequency
  real(dp)                      :: spring_constant
  logical                       :: soft_mode
  logical                       :: translational_mode
  type(RealVector), allocatable :: cos_vector(:)
  type(RealVector), allocatable :: sin_vector(:)
  integer                       :: qpoint_id_plus
  integer                       :: qpoint_id_minus
  integer                       :: subspace_id
  
  type(String), allocatable :: line(:)
  
  integer :: no_atoms
  
  select type(this); type is(RealMode)
    ! Read the id of this mode.
    line = split_line(input(1))
    id = int(line(4))
    
    ! Read the id of this mode's pair.
    line = split_line(input(2))
    paired_id = int(line(6))
    
    ! Read frequency and spring constant.
    line = split_line(input(3))
    frequency = dble(line(4))
    
    line = split_line(input(4))
    spring_constant = dble(line(4))
    
    ! Read whether or not this mode is soft.
    line = split_line(input(5))
    soft_mode = lgcl(line(5))
    
    ! Read whether or not this mode is purely translational.
    line = split_line(input(6))
    translational_mode = lgcl(line(5))
    
    ! Read the q-qpoint id of the plus mode.
    line = split_line(input(7))
    qpoint_id_plus = int(line(5))
    
    ! Read the q-qpoint id of the minus mode.
    line = split_line(input(8))
    qpoint_id_minus = int(line(5))
    
    ! Read the degeneracy id of this mode.
    line = split_line(input(9))
    subspace_id = int(line(4))
    
    ! Read in the vector associated with the mode.
    no_atoms = (size(input)-12)/2
    cos_vector = RealVector(input(12:11+no_atoms))
    sin_vector = RealVector(input(13+no_atoms:12+2*no_atoms))
    
    this = RealMode( id,                 &
                   & paired_id,          &
                   & frequency,          &
                   & spring_constant,    &
                   & soft_mode,          &
                   & translational_mode, &
                   & cos_vector,         &
                   & sin_vector,         &
                   & qpoint_id_plus,     &
                   & qpoint_id_minus,    &
                   & subspace_id         )
  class default
    call err()
  end select
end procedure

module procedure write_RealMode
  select type(this); type is(RealMode)
    output = [ 'Mode ID                   : '//this%id,                 &
             & 'ID of paired mode         : '//this%paired_id,          &
             & 'Mode frequency            : '//this%frequency,          &
             & 'Spring constant           : '//this%spring_constant,    &
             & 'Mode is soft              : '//this%soft_mode,          &
             & 'Mode purely translational : '//this%translational_mode, &
             & 'Positive q-point id       : '//this%qpoint_id_plus,     &
             & 'Negative q-point id       : '//this%qpoint_id_minus,    &
             & 'Subspace id               : '//this%subspace_id,        &
             & str('Mass-weighted unit vector in primitive cell:'),     &
             & str('cos component:'),                                   &
             & str(this%cos_vector),                                    &
             & str('sin component:'),                                   &
             & str(this%sin_vector)                                     ]
  class default
    call err()
  end select
end procedure

module procedure new_RealMode_Strings
  call this%read(input)
end procedure

module procedure new_RealMode_StringArray
  this = RealMode(str(input))
end procedure
end submodule
