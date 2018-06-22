! ======================================================================
! Harmonic normal modes, in complex co-ordinates.
! ======================================================================
module complex_mode_submodule
  use utils_module
  
  use structure_module
  implicit none
  
  private
  
  public :: ComplexMode
  public :: conjg
  public :: transform
  
  ! A normal mode in complex co-ordinates.
  type, extends(Stringsable) :: ComplexMode
    ! An id which is unique to each mode, and the unique id of the equivalent
    !    mode at the q-point q' s.t. q+q' is a G-vector.
    ! If 2q is a G-vector, then paired_id=id.
    integer :: id
    integer :: paired_id
    
    ! The frequency, and frequency-relevant information.
    real(dp) :: frequency
    logical  :: soft_mode          ! True if frequency < -epsilon.
    logical  :: translational_mode ! True if frequency=0 and at gamma.
    
    ! Unit vectors along the normal mode, in cartesian and mass-weighted
    !    co-ordinates.
    ! Normal modes are orthogonal in mass-weighted co-ordinates, so in general
    !    they are not orthogonal in cartesian co-ordinates.
    ! In both cases, the unit vectors are stored as displacements from
    !    equilibrium for each atom in the primitive cell.
    type(ComplexVector), allocatable :: mass_weighted_vector(:)
    type(ComplexVector), allocatable :: cartesian_vector(:)
    
    ! The ID of the q-point at which this mode exists,
    !    and that at which the pair to this mode exists.
    integer :: qpoint_id
    integer :: paired_qpoint_id
    
    ! The ID of the subspace modes which are degenerate with this mode.
    integer :: subspace_id
  contains
    procedure, public :: read  => read_ComplexMode
    procedure, public :: write => write_ComplexMode
  end type
  
  interface ComplexMode
    module procedure new_ComplexMode
    module procedure new_ComplexMode_HermitianEigenstuff
    module procedure new_ComplexMode_StringArray
  end interface
  
  interface conjg
    module procedure conjg_ComplexMode
  end interface
  
  interface transform
    module procedure transform_ComplexMode
  end interface
contains

! ----------------------------------------------------------------------
! Basic constructor.
! ----------------------------------------------------------------------
function new_ComplexMode(id,paired_id,frequency,soft_mode,translational_mode, &
   & mass_weighted_vector,cartesian_vector,qpoint_id,paired_qpoint_id,        &
   & subspace_id) result(this)
  implicit none
  
  integer,             intent(in) :: id
  integer,             intent(in) :: paired_id
  real(dp),            intent(in) :: frequency
  logical,             intent(in) :: soft_mode
  logical,             intent(in) :: translational_mode
  type(ComplexVector), intent(in) :: mass_weighted_vector(:)
  type(ComplexVector), intent(in) :: cartesian_vector(:)
  integer,             intent(in) :: qpoint_id
  integer,             intent(in) :: paired_qpoint_id
  integer,             intent(in) :: subspace_id
  type(ComplexMode)               :: this
  
  this%id                   = id
  this%paired_id            = paired_id
  this%frequency            = frequency
  this%soft_mode            = soft_mode
  this%translational_mode   = translational_mode
  this%mass_weighted_vector = mass_weighted_vector
  this%cartesian_vector     = cartesian_vector
  this%qpoint_id            = qpoint_id
  this%paired_qpoint_id     = paired_qpoint_id
  this%subspace_id          = subspace_id
end function

! ----------------------------------------------------------------------
! Construct a complex mode from the eigenstuff of the dynamical matrix.
! ----------------------------------------------------------------------
impure elemental function new_ComplexMode_HermitianEigenstuff(eigenstuff, &
   & structure) result(this)
  implicit none
  
  type(HermitianEigenstuff), intent(in) :: eigenstuff
  type(StructureData),       intent(in) :: structure
  type(ComplexMode)                     :: this
  
  real(dp) :: frequency
  
  type(ComplexVector), allocatable :: mass_weighted_vector(:)
  type(ComplexVector), allocatable :: cartesian_vector(:)
  
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
  allocate( mass_weighted_vector(structure%no_atoms_prim), &
          & cartesian_vector(structure%no_atoms_prim),     &
          & stat=ialloc); call err(ialloc)
  do i=1,structure%no_atoms_prim
    mass_weighted_vector(i) = eigenstuff%evec(3*i-2:3*i)
    cartesian_vector(i) = mass_weighted_vector(i) &
                      & * sqrt(structure%atoms(i)%mass())
  enddo
  
  ! Normalise cartesian vector. (Mass-weighted vector is already normalised).
  cartesian_vector = cartesian_vector &
                 & / sqrt(sum(real(cartesian_vector*conjg(cartesian_vector))))
  
  ! Call private constructor.
  this = ComplexMode( id                   = 0,                    &
                    & paired_id            = 0,                    &
                    & frequency            = frequency,            &
                    & soft_mode            = frequency<-1.0e-6_dp, &
                    & translational_mode   = .false.,              &
                    & mass_weighted_vector = mass_weighted_vector,  &
                    & cartesian_vector     = cartesian_vector,     &
                    & qpoint_id            = 0,                    &
                    & paired_qpoint_id     = 0,                    &
                    & subspace_id          = 0)
end function

! ----------------------------------------------------------------------
! Take the conjugate of the mode.
! This is equivalent to the transformation q -> -q.
! ----------------------------------------------------------------------
impure elemental function conjg_ComplexMode(input) result(output)
  implicit none
  
  type(ComplexMode), intent(in) :: input
  type(ComplexMode)             :: output
  
  output = ComplexMode(                                         &
     & id                   = input%paired_id,                  &
     & paired_id            = input%id,                         &
     & frequency            = input%frequency,                  &
     & soft_mode            = input%soft_mode,                  &
     & translational_mode   = input%translational_mode,         &
     & mass_weighted_vector = conjg(input%mass_weighted_vector), &
     & cartesian_vector     = conjg(input%cartesian_vector),    &
     & qpoint_id            = input%paired_qpoint_id,           &
     & paired_qpoint_id     = input%qpoint_id,                  &
     & subspace_id          = input%subspace_id)
end function

! ----------------------------------------------------------------------
! Transform a mode by a symmetry operator.
! ----------------------------------------------------------------------
impure elemental function transform_ComplexMode(input,symmetry,qpoint_from, &
   & qpoint_to) result(output)
  implicit none
  
  type(ComplexMode),      intent(in) :: input
  type(SymmetryOperator), intent(in) :: symmetry
  type(QpointData),       intent(in) :: qpoint_from
  type(QpointData),       intent(in) :: qpoint_to
  type(ComplexMode)                  :: output
  
  integer :: no_atoms
  integer :: atom_from
  integer :: atom_to
  
  type(IntVector) :: r
  
  type(ComplexVector), allocatable :: mass_weighted_vector(:)
  type(ComplexVector), allocatable :: cartesian_vector(:)
  
  integer :: ialloc
  
  ! Check that inputs are consistent.
  if (input%qpoint_id/=qpoint_from%id) then
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
  
  no_atoms = size(input%mass_weighted_vector)
  
  ! Allocate output, and transfer across all data.
  ! (Displacements need rotating, but everything else stays the same.)
  allocate( mass_weighted_vector(no_atoms), &
          & cartesian_vector(no_atoms),     &
          & stat=ialloc); call err(ialloc)
  output = input
  do atom_from=1,no_atoms
    atom_to = symmetry%prim_atom_group * atom_from
    r = symmetry%prim_rvector(atom_from)
    
    mass_weighted_vector(atom_to) =              &
       &   symmetry%cartesian_rotation           &
       & * input%mass_weighted_vector(atom_from) &
       & * exp_2pii(qpoint_to%qpoint*r)
    
    cartesian_vector(atom_to) =              &
       &   symmetry%cartesian_rotation       &
       & * input%cartesian_vector(atom_from) &
       & * exp_2pii(qpoint_to%qpoint*r)
  enddo
  
  output = ComplexMode( id                   = 0,                          &
                      & paired_id            = 0,                          &
                      & frequency            = input%frequency,            &
                      & soft_mode            = input%soft_mode,            &
                      & translational_mode   = input%translational_mode,   &
                      & mass_weighted_vector = mass_weighted_vector,       &
                      & cartesian_vector     = cartesian_vector,           &
                      & qpoint_id            = qpoint_to%id,               &
                      & paired_qpoint_id     = qpoint_to%paired_qpoint_id, &
                      & subspace_id          = input%subspace_id)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_ComplexMode(this,input)
  implicit none
  
  class(ComplexMode), intent(out) :: this
  type(String),       intent(in)  :: input(:)
  
  integer                          :: id
  integer                          :: paired_id
  real(dp)                         :: frequency
  logical                          :: soft_mode
  logical                          :: translational_mode
  type(ComplexVector), allocatable :: mass_weighted_vector(:)
  type(ComplexVector), allocatable :: cartesian_vector(:)
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
    
    ! Read frequency.
    line = split_line(input(3))
    frequency = dble(line(4))
    
    ! Read whether or not this mode is soft.
    line = split_line(input(4))
    soft_mode = lgcl(line(5))
    
    ! Read whether or not this mode is purely translational.
    line = split_line(input(5))
    translational_mode = lgcl(line(5))
    
    ! Read the q-point id of this mode.
    line = split_line(input(6))
    qpoint_id = int(line(4))
    
    ! Read the q-point id of this mode's pair.
    line = split_line(input(7))
    paired_qpoint_id = int(line(5))
    
    ! Read the degeneracy id of this mode.
    line = split_line(input(8))
    subspace_id = int(line(4))
    
    ! Read in the vectors along the direction of the mode.
    no_atoms = (size(input)-10)/2
    mass_weighted_vector = ComplexVector(input(10:9+no_atoms))
    cartesian_vector = ComplexVector(input(size(input)-no_atoms+1:))
    
    this = ComplexMode( id,                   &
                      & paired_id,            &
                      & frequency,            &
                      & soft_mode,            &
                      & translational_mode,   &
                      & mass_weighted_vector, &
                      & cartesian_vector,     &
                      & qpoint_id,            &
                      & paired_qpoint_id,     &
                      & subspace_id)
  end select
end subroutine

function write_ComplexMode(this) result(output)
  implicit none
  
  class(ComplexMode), intent(in) :: this
  type(String), allocatable      :: output(:)
  
  select type(this); type is(ComplexMode)
    output = [ 'Mode ID                   : '//this%id,                 &
             & 'ID of paired mode         : '//this%paired_id,          &
             & 'Mode frequency            : '//this%frequency,          &
             & 'Mode is soft              : '//this%soft_mode,          &
             & 'Mode purely translational : '//this%translational_mode, &
             & 'q-point id                : '//this%qpoint_id,          &
             & 'Paired q-point id         : '//this%paired_qpoint_id,   &
             & 'Subspace id               : '//this%subspace_id,        &
             & str('Mass-weighted displacements in primitive cell:'),   &
             & str(this%mass_weighted_vector),                          &
             & str('Cartesian displacements in primitive cell:'),       &
             & str(this%cartesian_vector)                               ]
  end select
end function

impure elemental function new_ComplexMode_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(ComplexMode)             :: this
  
  this = input
end function
end module
