! ======================================================================
! Convert a ModeVector to a set of atomic displacements,
!    or convert a set of forces to a ModeVector.
! ======================================================================
module mode_vector_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
  
  public :: DisplacementData
  public :: ForceData
  
  public :: real_mode_to_displacement
  public :: force_to_real_mode
  
  ! --------------------------------------------------
  ! Atomic displacements and forces in cartesian co-ordinates.
  ! --------------------------------------------------
  type DisplacementData
    type(RealVector), allocatable :: displacements(:)
  end type
  
  type ForceData
    type(RealVector), allocatable :: forces(:)
  end type
contains

! ----------------------------------------------------------------------
! Conversions between co-ordinate systems.
! ----------------------------------------------------------------------
! Converts a vector in normal mode co-ordinates to cartesian co-ordinates.
function real_mode_to_displacement(input,modes,qpoint,supercell) result(output)
  use constants_module, only : pi
  use qpoints_module
  use structure_module
  use linear_algebra_module
  use atom_module
  implicit none
  
  type(RealModeVector), intent(in) :: input
  type(RealMode),       intent(in) :: modes(:)
  type(QpointData),     intent(in) :: qpoint
  type(StructureData),  intent(in) :: supercell
  type(DisplacementData)           :: output
  
  ! The q-point q, q.R, cos(q.R) and sin(q.R).
  type(RealVector) :: q
  real(dp)         :: qr
  real(dp)         :: cos_qr
  real(dp)         :: sin_qr
  
  ! Atom data.
  type(AtomData) :: atom
  integer        :: prim
  
  ! Temporary variables
  integer :: i,j,ialloc
  
  ! Check that inputs are as expected.
  if (size(input%vector)/=size(modes)) then
    call print_line(CODE_ERROR//': The modes and mode vector do not match.')
    call err()
  elseif (supercell%no_modes_prim/=size(modes)) then
    call print_line(CODE_ERROR//': The number of normal modes do not match up.')
  endif
  
  ! Perform conversion.
  allocate( output%displacements(size(supercell%atoms)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(supercell%atoms)
    output%displacements(i) = dble(int(zeroes(3)))
    
    atom = supercell%atoms(i)
    prim = atom%prim_id()
    
    ! Calculate 2*pi*q.R, cos(2*pi*i*q.R) and sin(2*pi*i*q.R).
    q = dblevec(qpoint%qpoint)
    qr = 2 * pi * q * supercell%rvectors(prim)
    cos_qr = cos(qr)
    sin_qr = sin(qr)
    
    ! Calculate displacements in cartesian co-ordinates.
    do j=1,size(modes)
      if (modes(j)%at_paired_qpoint) then
        output%displacements(i) = output%displacements(i)                  &
                              & + input%vector(j,1)                        &
                              & * modes(j)%primitive_displacements(prim,1) &
                              & * cos_qr                                   &
                              & / supercell%atoms(i)%mass()
      else
        output%displacements(i) = output%displacements(i)                  &
                              & + input%vector(j,1)                        &
                              & * modes(j)%primitive_displacements(prim,1) &
                              & * cos_qr                                   &
                              & * sqrt(2.0_dp)                             &
                              & / supercell%atoms(i)%mass()                &
                              & + input%vector(j,2)                        &
                              & * modes(j)%primitive_displacements(prim,2) &
                              & * sin_qr                                   &
                              & * sqrt(2.0_dp)                             &
                              & / supercell%atoms(i)%mass()
      endif
    enddo
  enddo
end function

! Converts a force in cartesian co-ordinates into real normal mode co-ordinates.
function force_to_real_mode(input,modes,qpoint,supercell) result(output)
  use constants_module, only : pi
  use qpoints_module
  use structure_module
  use linear_algebra_module
  use atom_module
  implicit none
  
  type(ForceData),     intent(in) :: input
  type(RealMode),      intent(in) :: modes(:)
  type(QpointData),    intent(in) :: qpoint
  type(StructureData), intent(in) :: supercell
  type(RealModeVector)            :: output
  
  ! The q-point q, q.R, cos(q.R) and sin(q.R).
  type(RealVector) :: q
  real(dp)         :: qr
  real(dp)         :: cos_qr
  real(dp)         :: sin_qr
  
  ! Atom data.
  type(AtomData) :: atom
  integer        :: prim
  
  ! Temporary variables
  integer :: i,j,ialloc
  
  ! Allocate output.
  output%at_paired_qpoint = modes(1)%at_paired_qpoint
  if (output%at_paired_qpoint) then
    allocate( output%vector(size(modes),1), &
            & stat=ialloc); call err(ialloc)
  else
    allocate( output%vector(size(modes),2), &
            & stat=ialloc); call err(ialloc)
  endif
  output%vector = 0
  
  ! Convert force into normal mode co-ordinates.
  do i=1,supercell%no_atoms
    atom = supercell%atoms(i)
    prim = atom%prim_id()
    
    ! Calculate 2*pi*q.R, sin(2*pi*i*q.R) and cos(2*pi*i*q.R).
    q = dblevec(qpoint%qpoint)
    qr = 2 * pi * q * supercell%rvectors(atom%rvec_id())
    cos_qr = cos(qr)
    sin_qr = sin(qr)
    
    ! Calculate displacements in normal-mode co-ordinates.
    do j=1,size(modes)
      ! Calculate the dot product of the input vector with the normal mode.
      if (output%at_paired_qpoint) then
        output%vector(j,1) = output%vector(j,1)                       &
                         & + input%forces(i)                          &
                         & * modes(j)%primitive_displacements(prim,1) &
                         & * cos_qr                                   &
                         & / ( supercell%atoms(i)%mass()              &
                         &   * supercell%sc_size )
      else
        output%vector(j,1) = output%vector(j,1)                       &
                         & + input%forces(i)                          &
                         & * modes(j)%primitive_displacements(prim,1) &
                         & * cos_qr                                   &
                         & * sqrt(2.0_dp)                             &
                         & / ( supercell%atoms(i)%mass()              &
                         &   * supercell%sc_size )
        output%vector(j,2) = output%vector(j,2)                       &
                         & + input%forces(i)                          &
                         & * modes(j)%primitive_displacements(prim,2) &
                         & * sin_qr                                   &
                         & * sqrt(2.0_dp)                             &
                         & / ( supercell%atoms(i)%mass()              &
                         &   * supercell%sc_size )
      endif
    enddo
  enddo
end function
end module
