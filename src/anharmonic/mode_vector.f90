! ======================================================================
! Vectors of normal modes.
! ======================================================================
module mode_vector_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
  
  private
  
  public :: DisplacementData
  public :: ForceData
  
  public :: ComplexModeVector
  public :: RealModeVector
  
  public :: complex_to_real
  public :: real_to_complex
  
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
  
  ! --------------------------------------------------
  ! The coefficients of a point in terms of complex normal mode co-ordinates.
  ! --------------------------------------------------
  type ComplexModeVector
    ! Whether or not 2q=G. If true, there is only one mode, and it is real.
    !    If not, there are two modes, and they are conjugates of one another.
    logical :: at_paired_qpoint
    
    complex(dp), allocatable :: vector(:,:)
  contains
    generic,   public  :: operator(+) => add_ComplexModeVectors
    procedure, private ::                add_ComplexModeVectors
    generic,   public  :: operator(-) => subtract_ComplexModeVectors
    procedure, private ::                subtract_ComplexModeVectors
  end type
  
  ! --------------------------------------------------
  ! The coefficients of a point in terms of real normal mode co-ordinates.
  ! --------------------------------------------------
  type RealModeVector
    ! Whether or not 2q=G. If true, there is only the cosine mode.
    !    If not, there is a cosine mode and a sine mode.
    logical :: at_paired_qpoint
    
    real(dp), allocatable :: vector(:,:)
  contains
    generic,   public  :: operator(+) => add_RealModeVectors
    procedure, private ::                add_RealModeVectors
    generic,   public  :: operator(-) => subtract_RealModeVectors
    procedure, private ::                subtract_RealModeVectors
  end type
  
  ! --------------------------------------------------
  ! Conversions between complex and real co-ordinates.
  ! --------------------------------------------------
  interface complex_to_real
    module procedure complex_to_real_ModeVector
  end interface
  
  interface real_to_complex
    module procedure real_to_complex_ModeVector
  end interface
contains

! ----------------------------------------------------------------------
! ComplexModeVector procedures.
! ----------------------------------------------------------------------
function add_ComplexModeVectors(this,that) result(output)
  implicit none
  
  class(ComplexModeVector), intent(in) :: this
  class(ComplexModeVector), intent(in) :: that
  type(ComplexModeVector)              :: output
  
  ! Check that inputs match and are as expected.
  if (this%at_paired_qpoint .neqv. that%at_paired_qpoint) then
    call print_line(CODE_ERROR//': Attempted to add mode vectors at different &
       &q-points.')
    call err()
  elseif (any(shape(this%vector)/=shape(that%vector))) then
    call print_line(CODE_ERROR//': Attempted to add mode vectors of different &
       &shapes.')
    call err()
  elseif (this%at_paired_qpoint .and. size(this%vector,2)/=1) then
    call print_line(CODE_ERROR//': Mode vector of unexpected shape.')
    call err()
  elseif ((.not. this%at_paired_qpoint) .and. size(this%vector,2)/=2) then
    call print_line(CODE_ERROR//': Mode vector of unexpected shape.')
    call err()
  endif
  
  ! Add vectors together.
  output%at_paired_qpoint = this%at_paired_qpoint
  output%vector = this%vector + that%vector
end function

function subtract_ComplexModeVectors(this,that) result(output)
  implicit none
  
  class(ComplexModeVector), intent(in) :: this
  class(ComplexModeVector), intent(in) :: that
  type(ComplexModeVector)              :: output
  
  ! Check that inputs match and are as expected.
  if (this%at_paired_qpoint .neqv. that%at_paired_qpoint) then
    call print_line(CODE_ERROR//': Attempted to subtract mode vectors at &
       &different q-points.')
    call err()
  elseif (any(shape(this%vector)/=shape(that%vector))) then
    call print_line(CODE_ERROR//': Attempted to subtract mode vectors of &
       &different shapes.')
    call err()
  elseif (this%at_paired_qpoint .and. size(this%vector,2)/=1) then
    call print_line(CODE_ERROR//': Mode vector of unexpected shape.')
    call err()
  elseif ((.not. this%at_paired_qpoint) .and. size(this%vector,2)/=2) then
    call print_line(CODE_ERROR//': Mode vector of unexpected shape.')
    call err()
  endif
  
  ! Subtract vectors.
  output%at_paired_qpoint = this%at_paired_qpoint
  output%vector = this%vector - that%vector
end function

! ----------------------------------------------------------------------
! RealModeVector procedures.
! ----------------------------------------------------------------------
function add_RealModeVectors(this,that) result(output)
  implicit none
  
  class(RealModeVector), intent(in) :: this
  class(RealModeVector), intent(in) :: that
  type(RealModeVector)              :: output
  
  ! Check that inputs match and are as expected.
  if (this%at_paired_qpoint .neqv. that%at_paired_qpoint) then
    call print_line(CODE_ERROR//': Attempted to add mode vectors at different &
       &q-points.')
    call err()
  elseif (any(shape(this%vector)/=shape(that%vector))) then
    call print_line(CODE_ERROR//': Attempted to add mode vectors of different &
       &shapes.')
    call err()
  elseif (this%at_paired_qpoint .and. size(this%vector,2)/=1) then
    call print_line(CODE_ERROR//': Mode vector of unexpected shape.')
    call err()
  elseif ((.not. this%at_paired_qpoint) .and. size(this%vector,2)/=2) then
    call print_line(CODE_ERROR//': Mode vector of unexpected shape.')
    call err()
  endif
  
  ! Add vectors together.
  output%at_paired_qpoint = this%at_paired_qpoint
  output%vector = this%vector + that%vector
end function

function subtract_RealModeVectors(this,that) result(output)
  implicit none
  
  class(RealModeVector), intent(in) :: this
  class(RealModeVector), intent(in) :: that
  type(RealModeVector)              :: output
  
  ! Check that inputs match and are as expected.
  if (this%at_paired_qpoint .neqv. that%at_paired_qpoint) then
    call print_line(CODE_ERROR//': Attempted to subtract mode vectors at &
       &different q-points.')
    call err()
  elseif (any(shape(this%vector)/=shape(that%vector))) then
    call print_line(CODE_ERROR//': Attempted to subtract mode vectors of &
       &different shapes.')
    call err()
  elseif (this%at_paired_qpoint .and. size(this%vector,2)/=1) then
    call print_line(CODE_ERROR//': Mode vector of unexpected shape.')
    call err()
  elseif ((.not. this%at_paired_qpoint) .and. size(this%vector,2)/=2) then
    call print_line(CODE_ERROR//': Mode vector of unexpected shape.')
    call err()
  endif
  
  ! Subtract vectors.
  output%at_paired_qpoint = this%at_paired_qpoint
  output%vector = this%vector - that%vector
end function

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
  if (input%at_paired_qpoint .neqv. modes(1)%at_paired_qpoint) then
    call print_line(CODE_ERROR//': The modes and mode vector do not match.')
    call err()
  elseif (size(input%vector,2)/=size(modes)) then
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
      if (input%at_paired_qpoint) then
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

function complex_to_real_ModeVector(input) result(output)
  implicit none
  
  type(ComplexModeVector), intent(in) :: input
  type(RealModeVector)                :: output
  
  integer :: ialloc
  
  if (size(input%vector,1)==1) then
    ! Out*x = In*x
    ! => Out = In
    output%vector = real(input%vector)
  elseif (size(input%vector,1)==2) then
    ! Out1*(x+ + x-)/sqrt(2) + Out2*(x+ - x-)/(sqrt(2)i) = In1*x+ + In2*x-
    ! => Out1 =  (In1+In2)/sqrt(2)    = Real(In1)*sqrt(2)
    !    Out2 = -(In1-In2)/(sqrt(2)i) = -Imag(In1)*sqrt(2)
    allocate( output%vector(size(input%vector,1),size(input%vector,2)), &
            & stat=ialloc); call err(ialloc)
    output%vector(1,:) = real(input%vector(1,:)) * sqrt(2.0_dp)
    output%vector(2,:) = -aimag(input%vector(1,:)) * sqrt(2.0_dp)
  else
    call print_line(CODE_ERROR//': Normal mode vector of &
       &unexpected shape.')
    call err()
  endif
end function

function real_to_complex_ModeVector(input) result(output)
  implicit none
  
  type(RealModeVector), intent(in) :: input
  type(ComplexModeVector)          :: output
  
  integer :: ialloc
  
  if (size(input%vector,1)==1) then
    ! O1*x(q) = I1*x(q)
    ! => O1 = I1
    output%vector = input%vector
  elseif (size(input%vector,1)==2) then
    ! O1*x(q)+O2*x(G-q) = I1*(x(q)+x(G-q))/sqrt(2) + I2*(x(q)-x(G-q))/(sqrt(2)i)
    ! => 01 = (I1-iI2)/sqrt(2)
    !    O2 = (I1+iI2)/sqrt(2)
    allocate( output%vector(size(input%vector,1),size(input%vector,2)), &
            & stat=ialloc); call err(ialloc)
    output%vector(1,:) = cmplx(input%vector(1,:),-input%vector(1,:),dp) &
                     & / sqrt(2.0_dp)
    output%vector(2,:) = cmplx(input%vector(1,:), input%vector(1,:),dp) &
                     & / sqrt(2.0_dp)
  else
    call print_line(CODE_ERROR//': Normal mode vector of &
       &unexpected shape.')
    call err()
  endif
end function
end module
