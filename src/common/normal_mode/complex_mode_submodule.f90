submodule (caesar_complex_mode_module) caesar_complex_mode_submodule
  use caesar_normal_mode_module
contains

module procedure new_ComplexMode
  this%id                 = id
  this%paired_id          = paired_id
  this%frequency          = frequency
  this%spring_constant    = spring_constant
  this%soft_mode          = soft_mode
  this%translational_mode = translational_mode
  this%unit_vector        = unit_vector
  this%qpoint_id          = qpoint_id
  this%paired_qpoint_id   = paired_qpoint_id
  this%subspace_id        = subspace_id
end procedure

module procedure new_ComplexMode_unprocessed
  real(dp) :: spring_constant
  
  if (frequency>=0.0_dp) then
    spring_constant = frequency**2
  else
    spring_constant = -frequency**2
  endif
  
  this = ComplexMode( id                 = 0,                    &
                    & paired_id          = 0,                    &
                    & frequency          = frequency,            &
                    & spring_constant    = spring_constant,      &
                    & soft_mode          = frequency<-1.0e-6_dp, &
                    & translational_mode = .false.,              &
                    & unit_vector        = unit_vector,          &
                    & qpoint_id          = 0,                    &
                    & paired_qpoint_id   = 0,                    &
                    & subspace_id        = 0                     )
end procedure

module procedure new_ComplexMode_HermitianEigenstuff
  real(dp)                         :: frequency
  type(ComplexVector), allocatable :: unit_vector(:)
  
  integer :: i,ialloc
  
  if (size(eigenstuff%evec)/=structure%no_modes_prim) then
    call print_line(ERROR//': incompatible eigenvector and structure.')
    call err()
  endif
  
  ! Calculate frequency.
  ! V = sum[ 0.5*w^2*u^2 ] where w is the frequency.
  ! F = -2V = sum[ -w^2*u^2 ]
  ! -> the evals of F are -w^2.
  if (eigenstuff%eval>=0.0_dp) then
    ! Unstable mode.
    frequency = - sqrt(eigenstuff%eval)
  else
    ! Stable mode.
    frequency = sqrt(- eigenstuff%eval)
  endif
  
  ! Convert eigenvector with dimension no_modes into an array of
  !    3D vectors.
  allocate( unit_vector(structure%no_atoms_prim), &
          & stat=ialloc); call err(ialloc)
  do i=1,structure%no_atoms_prim
    unit_vector(i) = vec(eigenstuff%evec(3*i-2:3*i))
  enddo
  
  ! Call private constructor.
  this = ComplexMode( frequency   = frequency,  &
                    & unit_vector = unit_vector )
end procedure

module procedure conjg_ComplexMode
  output = ComplexMode(                               &
     & id                 = input%paired_id,          &
     & paired_id          = input%id,                 &
     & frequency          = input%frequency,          &
     & spring_constant    = input%spring_constant,    &
     & soft_mode          = input%soft_mode,          &
     & translational_mode = input%translational_mode, &
     & unit_vector        = conjg(input%unit_vector), &
     & qpoint_id          = input%paired_qpoint_id,   &
     & paired_qpoint_id   = input%qpoint_id,          &
     & subspace_id        = input%subspace_id         )
end procedure

module procedure transform_ComplexMode
  integer :: no_atoms
  integer :: no_modes
  integer :: atom_from
  integer :: atom_to
  
  type(IntVector) :: r
  
  type(ComplexVector), allocatable :: unit_vector(:)
  
  integer :: id
  integer :: paired_id
  
  integer :: ialloc
  
  ! Check that inputs are consistent.
  if (size(symmetry%atom_group)/=size(input%unit_vector)) then
    call print_line(CODE_ERROR//': Symmetry and Mode do not have the same &
       & number of atoms. Modes must only be transformed using symmetries of &
       & the primitive structure.')
  elseif (input%qpoint_id/=qpoint_from%id) then
    call print_line(CODE_ERROR//': mode and q-point not compatible.')
    call print_line(input%qpoint_id//' '//qpoint_from%id)
    call err()
  elseif (input%paired_qpoint_id/=qpoint_from%paired_qpoint_id) then
    call print_line(CODE_ERROR//': mode and q-point not compatible.')
    call err()
  elseif (symmetry * qpoint_from /= qpoint_to) then
    call print_line(CODE_ERROR//': Symmetry does not transform q-points as &
       &expected.')
    call err()
  endif
  
  no_atoms = size(input%unit_vector)
  no_modes = 3*no_atoms
  
  if (input%id==0) then
    id = 0
    paired_id = 0
  else
    id = input%id + no_modes*(qpoint_to%id-qpoint_from%id)
    paired_id = input%paired_id                         &
            & + no_modes*( qpoint_to%paired_qpoint_id   &
            &            - qpoint_from%paired_qpoint_id )
  endif
  
  allocate( unit_vector(no_atoms), &
          & stat=ialloc); call err(ialloc)
  output = input
  do atom_from=1,no_atoms
    atom_to = symmetry%atom_group * atom_from
    r = symmetry%rvectors(atom_from)
    
    unit_vector(atom_to) = symmetry%cartesian_tensor    &
                       & * input%unit_vector(atom_from) &
                       & * exp_2pii(-qpoint_to%qpoint*r)
  enddo
  
  output = ComplexMode( id                 = id,                         &
                      & paired_id          = paired_id,                  &
                      & frequency          = input%frequency,            &
                      & spring_constant    = input%spring_constant,      &
                      & soft_mode          = input%soft_mode,            &
                      & translational_mode = input%translational_mode,   &
                      & unit_vector        = unit_vector,                &
                      & qpoint_id          = qpoint_to%id,               &
                      & paired_qpoint_id   = qpoint_to%paired_qpoint_id, &
                      & subspace_id        = input%subspace_id)
end procedure

module procedure select_qpoint_ComplexMode
  output = qpoints(first(qpoints%id==mode%qpoint_id))
end procedure

module procedure select_qpoints_ComplexModes
  integer :: i,ialloc
  
  allocate(output(size(modes)), stat=ialloc); call err(ialloc)
  do i=1,size(modes)
    output(i) = select_qpoint(modes(i), qpoints)
  enddo
end procedure

module procedure generate_translational_modes
  type(QpointData) :: gamma_qpoint
  
  integer                          :: id
  complex(dp)                      :: vector(3)
  type(ComplexVector), allocatable :: unit_vector(:)
  
  integer :: i,j,ialloc
  
  allocate(output(3), stat=ialloc); call err(ialloc)
  
  gamma_qpoint = qpoints(first(qpoints%is_gvector()))
  
  do i=1,3
    id = (gamma_qpoint%id-1)*structure%no_modes+i
    vector = [((0.0_dp,0.0_dp),j=1,3)]
    vector(i) = cmplx(sqrt(1.0_dp/structure%no_modes),0.0_dp,dp)
    unit_vector = [(vec(vector),j=1,structure%no_modes)]
    output(i) = ComplexMode( id                 = id,              &
                           & paired_id          = id,              &
                           & frequency          = 0.0_dp,          &
                           & spring_constant    = 0.0_dp,          &
                           & soft_mode          = .false.,         &
                           & translational_mode = .true.,          &
                           & unit_vector        = unit_vector,     &
                           & qpoint_id          = gamma_qpoint%id, &
                           & paired_qpoint_id   = gamma_qpoint%id, &
                           & subspace_id        = 0                )
  enddo
end procedure

module procedure stress_prefactor
  integer :: i
  
  output = dblemat(zeroes(3,3))
  do i=1,size(this%unit_vector)
    if (present(that)) then
      if (this%qpoint_id==that%qpoint_id) then
        output = output                                         &
             & + real(outer_product( this%unit_vector(i),       &
             &                       conjg(that%unit_vector(i)) ))
      elseif (this%qpoint_id==that%paired_qpoint_id) then
        output = output                                   &
             & + real(outer_product( this%unit_vector(i), &
             &                       that%unit_vector(i)  ))
      else
        call print_line(ERROR//': Trying to construct the stress prefactor &
           &between modes at incompatible q-points.')
        call err()
      endif
    else
      output = output                                         &
           & + real(outer_product( this%unit_vector(i),       &
           &                       conjg(this%unit_vector(i)) ))
    endif
  enddo
end procedure

module procedure read_ComplexMode
  integer                          :: id
  integer                          :: paired_id
  real(dp)                         :: frequency
  real(dp)                         :: spring_constant
  logical                          :: soft_mode
  logical                          :: translational_mode
  type(ComplexVector), allocatable :: unit_vector(:)
  integer                          :: qpoint_id
  integer                          :: paired_qpoint_id
  integer                          :: subspace_id
  
  type(String), allocatable :: line(:)
  
  integer :: no_atoms
  
  select type(this); type is(ComplexMode)
    ! Read the id of this mode.
    line = split_line(input(1))
    id = int(line(4))
    
    ! Read the id of this mode's pair.
    line = split_line(input(2))
    paired_id = int(line(6))
    
    ! Read frequency and spring_constant.
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
    
    ! Read the q-point id of this mode.
    line = split_line(input(7))
    qpoint_id = int(line(4))
    
    ! Read the q-point id of this mode's pair.
    line = split_line(input(8))
    paired_qpoint_id = int(line(5))
    
    ! Read the degeneracy id of this mode.
    line = split_line(input(9))
    subspace_id = int(line(4))
    
    ! Read in the vectors along the direction of the mode.
    no_atoms = size(input)-10
    unit_vector = ComplexVector(input(11:10+no_atoms))
    
    this = ComplexMode( id,                 &
                      & paired_id,          &
                      & frequency,          &
                      & spring_constant,    &
                      & soft_mode,          &
                      & translational_mode, &
                      & unit_vector,        &
                      & qpoint_id,          &
                      & paired_qpoint_id,   &
                      & subspace_id)
  class default
    call err()
  end select
end procedure

module procedure write_ComplexMode
  select type(this); type is(ComplexMode)
    output = [ 'Mode ID                   : '//this%id,                 &
             & 'ID of paired mode         : '//this%paired_id,          &
             & 'Mode frequency            : '//this%frequency,          &
             & 'Spring constant           : '//this%spring_constant,    &
             & 'Mode is soft              : '//this%soft_mode,          &
             & 'Mode purely translational : '//this%translational_mode, &
             & 'q-point id                : '//this%qpoint_id,          &
             & 'Paired q-point id         : '//this%paired_qpoint_id,   &
             & 'Subspace id               : '//this%subspace_id,        &
             & str('Mass-weighted displacements in primitive cell:'),   &
             & str(this%unit_vector)                                    ]
  class default
    call err()
  end select
end procedure

module procedure new_ComplexMode_Strings
  call this%read(input)
end procedure

module procedure new_ComplexMode_StringArray
  this = ComplexMode(str(input))
end procedure
end submodule
