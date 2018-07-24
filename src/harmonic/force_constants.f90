! ======================================================================
! Reads harmonic forces, and generates the matrix of force constants.
! ======================================================================
module force_constants_module
  use common_module
  
  use unique_directions_module
  implicit none
  
  private
  
  public :: ForceConstants
  
  ! The matrix of force constants, F, such that F.x=f, where x and f are the
  !    displacement and force respectively.
  type, extends(NoDefaultConstructor) :: ForceConstants
    type(RealMatrix), allocatable, private :: constants_(:,:)
  contains
    procedure, public :: constants
    procedure, public :: check
  end type
  
  interface ForceConstants
    module procedure new_ForceConstants_forces
    module procedure new_ForceConstants_constants
  end interface
contains

! ----------------------------------------------------------------------
! Uses symmetry operations to construct force constants.
! ----------------------------------------------------------------------
! x is the collective vector of displacements, {xi}.
! f is the collective vector of forces, {fi}.
! Both are supercell%no_modes long.
!
! Under the harmonic approximation, the force constants, F, are defined as:
!    U = - x.F.x / 2
!
!    f = -dU/dx = F.x
!
! Under the symmetry s:
!    x -> x' = R.x
!    f -> f' = R.f
!
! F can be found by minimising L:
!    L = sum(x,s)[ (f'-F.x')^2 ]
!      = sum(x,s)[ f'.f' - 2*f'.F.x' + x'.F.F.x' ]
!
! => 0 = dL/dF = -2 * sum(x,s)[ f'^x' - F.x'^x' ]
! 
! where ^ is the outer product.
!
! => F = sum(x,s)[f'^x'] . inverse( sum(x,s)[x'^x'] )
!
! sum(x,s)[x'^x'] is block diagonal, so can be inverted in 3x3 blocks.
function new_ForceConstants_forces(supercell,unique_directions,sdir, &
   & acoustic_sum_rule_forces,calculation_reader,logfile)  result(output)
  implicit none
  
  type(StructureData),     intent(in)    :: supercell
  type(UniqueDirection),   intent(in)    :: unique_directions(:)
  type(String),            intent(in)    :: sdir
  logical,                 intent(in)    :: acoustic_sum_rule_forces
  type(CalculationReader), intent(inout) :: calculation_reader
  type(OFile),             intent(inout) :: logfile
  type(ForceConstants)                   :: output
  
  
  ! Forces (mass reduced).
  type(RealVector), allocatable :: forces(:,:)
  
  ! sum(s)[ f'^x' ].
  type(RealMatrix), allocatable :: fx(:,:)
  
  ! sum(s)[ x'^x' ] (diagonal blocks only).
  type(RealMatrix), allocatable :: xx(:)
  
  ! Read in forces (mass reduced).
  forces = read_forces( supercell,                &
                      & unique_directions,        &
                      & sdir,                     &
                      & acoustic_sum_rule_forces, &
                      & calculation_reader        )
  
  ! Construct sum[x'^x'].
  xx = construct_xx(unique_directions, supercell, logfile)
  
  ! Construct sum[f'^x'].
  fx = construct_fx(unique_directions, forces, supercell, logfile)
  
  ! Construct F = sum[f'^x'] . inverse(sum[x'^x']).
  output = construct_f(xx,fx,supercell,logfile)
  
  ! Check output.
  call output%check( forces,            &
                   & supercell,         &
                   & unique_directions, &
                   & logfile)
end function

! ----------------------------------------------------------------------
! Constructs and checks force constants from given matrices.
! ----------------------------------------------------------------------
function new_ForceConstants_constants(structure,force_constants,logfile) &
   & result(this)
  implicit none
  
  type(StructureData), intent(in)    :: structure
  type(RealMatrix),    intent(in)    :: force_constants(:,:)
  type(OFile),         intent(inout) :: logfile
  type(ForceConstants)               :: this
  
  ! Variables for checking force constants.
  type(AtomData)   :: atom_i
  type(AtomData)   :: atom_j
  type(AtomData)   :: atom_k
  integer          :: rvec_i
  integer          :: rvec_j
  integer          :: rvec_k
  integer          :: rvec_ij
  type(RealMatrix) :: matrix
  type(RealMatrix) :: copy
  real(dp)         :: average
  real(dp)         :: difference
  
  integer :: i,j,k
  
  ! Copy force constants into ForceConstants.
  this%constants_ = force_constants
  
  ! Check that the force constants between atom i and atom j are the same as
  !    those between atom i + R and atom j + R.
  average = 0.0_dp
  difference = 0.0_dp
  do i=1,structure%no_atoms
    atom_i = structure%atoms(i)
    rvec_i = atom_i%rvec_id()
    do j=1,structure%no_atoms
      atom_j = structure%atoms(j)
      rvec_j = atom_j%rvec_id()
      rvec_ij = structure%paired_rvector_group(rvec_i) * rvec_j
      
      do k=1,structure%no_atoms
        atom_k = structure%atoms(k)
        rvec_k = atom_k%rvec_id()
        if (atom_k%prim_id()/=atom_j%prim_id()) then
          cycle
        elseif (rvec_k/=rvec_ij) then
          cycle
        endif
        matrix = this%constants_(atom_i%id(),atom_j%id())
        copy = this%constants_(atom_i%prim_id(),atom_k%id())
        average = average + sum_squares((matrix+copy)/2)
        difference = difference + sum_squares(matrix-copy)
      enddo
    enddo
  enddo
  call logfile%print_line('Fractional L2 difference in force constants at &
     &different R-vectors: '//sqrt(difference/average))
  if (sqrt(difference/average)>1.0e-10_dp) then
    call print_line(WARNING//': Reconstructed force constants do not obey &
       &R-vector symmetries. Please check log files.')
    call err()
  endif
end function

! ----------------------------------------------------------------------
! Return the force constants between two atoms.
! ----------------------------------------------------------------------
function constants(this,a,b) result(output)
  implicit none
  
  class(ForceConstants), intent(in) :: this
  type(AtomData),        intent(in) :: a
  type(AtomData),        intent(in) :: b
  type(RealMatrix)                  :: output
  
  output = this%constants_(a%id(),b%id())
end function
  
! ----------------------------------------------------------------------
! Read in forces, and mass-weight them.
! ----------------------------------------------------------------------
function read_forces(supercell,unique_directions,sdir, &
   & acoustic_sum_rule_forces,calculation_reader) result(output)
  implicit none
  
  type(StructureData),     intent(in)    :: supercell
  type(UniqueDirection),   intent(in)    :: unique_directions(:)
  type(String),            intent(in)    :: sdir
  logical,                 intent(in)    :: acoustic_sum_rule_forces
  type(CalculationReader), intent(inout) :: calculation_reader
  type(RealVector), allocatable          :: output(:,:)
  
  ! DFT output data.
  type(ElectronicStructure) :: electronic_structure
  
  ! Direction information.
  type(String)   :: directory
  type(AtomData) :: atom
  type(String)   :: direction
  type(String)   :: atom_string
  
  ! Temporary variables.
  integer          :: i,j,ialloc
  type(RealVector) :: total
  
  allocate( output(supercell%no_atoms, size(unique_directions)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(unique_directions)
    atom = supercell%atoms(unique_directions(i)%atom_id)
    direction = unique_directions(i)%direction
    atom_string = left_pad(atom%id(),str(supercell%no_atoms_prim))
    directory = sdir//'/atom.'//atom_string//'.'//direction
    
    electronic_structure = calculation_reader%read_calculation(directory)
    output(:,i) = electronic_structure%forces%vectors
    
    ! Enforce continuous translational symmetry,
    !    i.e. ensure sum(output,1)=0.
    ! This is the 'Accoustic sum rule'.
    if (acoustic_sum_rule_forces) then
      total = dblevec(zeroes(3))
      do j=1,supercell%no_atoms
        total = total+output(j,i)
      enddo
      do j=1,supercell%no_atoms
        output(j,i) = output(j,i)-total/supercell%no_atoms
      enddo
    endif
    
    ! Mass-reduce forces.
    do j=1,supercell%no_atoms
      output(j,i) = output(j,i) / sqrt(atom%mass()*supercell%atoms(j)%mass())
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Construct and check sum[x'^x'].
! ----------------------------------------------------------------------
function construct_xx(unique_directions,supercell,logfile) result(output)
  implicit none
  
  type(UniqueDirection), intent(in)    :: unique_directions(:)
  type(StructureData),   intent(in)    :: supercell
  type(OFile),           intent(inout) :: logfile
  type(RealMatrix), allocatable        :: output(:)
  
  ! Atoms.
  type(AtomData) :: atom_1
  type(AtomData) :: atom_1p
  
  ! A displacement.
  type(RealVector) :: x
  
  ! Variables for checking the output.
  type(RealMatrix) :: matrix
  type(RealMatrix) :: symmetric
  real(dp)         :: average
  real(dp)         :: difference
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  ! Construct sum[x'^x'].
  allocate(output(supercell%no_atoms), stat=ialloc); call err(ialloc)
  output = dblemat(zeroes(3,3))
  do i=1,size(supercell%symmetries)
    do j=1,size(unique_directions)
      atom_1 = supercell%atoms(unique_directions(j)%atom_id)
      atom_1p = supercell%atoms( supercell%symmetries(i)%atom_group &
                             & * atom_1%id())
      
      x = supercell%symmetries(i)%cartesian_tensor &
      & * unique_directions(j)%atomic_displacement
      
      output(atom_1p%id()) = output(atom_1p%id()) + outer_product(x,x)
    enddo
  enddo
  
  ! Check sum[x'^x'] obeys correct symmetries.
  average = 0
  difference = 0
  do i=1,size(supercell%symmetries)
    do j=1,size(unique_directions)
      atom_1 = supercell%atoms(unique_directions(j)%atom_id)
      atom_1p = supercell%atoms( supercell%symmetries(i)%atom_group &
                             & * atom_1%id())
      matrix = output(atom_1p%id())
      symmetric = supercell%symmetries(i)%cartesian_tensor &
              & * output(atom_1%id())                      &
              & * transpose(supercell%symmetries(i)%cartesian_tensor)
      average = average + sum_squares((matrix+symmetric)/2)
      difference = difference + sum_squares(matrix-symmetric)
    enddo
  enddo
  call logfile%print_line(                                  &
     & 'Fractional L2 error in symmetry of x^x        : '// &
     & sqrt(difference/average))
  if (sqrt(difference/average) > 1e-10_dp) then
    call print_line(WARNING//': x^x is not as symmetric as expected. Please &
       &check log files.')
  endif
end function

! ----------------------------------------------------------------------
! Construct and check sum[f'^x'].
! ----------------------------------------------------------------------
function construct_fx(unique_directions,forces,supercell,logfile) &
   & result(output)
  implicit none
  
  type(UniqueDirection), intent(in)    :: unique_directions(:)
  type(RealVector),      intent(in)    :: forces(:,:)
  type(StructureData),   intent(in)    :: supercell
  type(OFile),           intent(inout) :: logfile
  type(RealMatrix), allocatable        :: output(:,:)
  
  ! Atoms.
  type(AtomData) :: atom_1
  type(AtomData) :: atom_1p
  type(AtomData) :: atom_2
  type(AtomData) :: atom_2p
  
  ! A displacement and a force.
  type(RealVector) :: x
  type(RealVector) :: f
  
  ! Variables for checking the output.
  logical, allocatable :: updated(:,:)
  type(RealMatrix)     :: matrix
  type(RealMatrix)     :: symmetric
  real(dp)             :: average
  real(dp)             :: difference
  
  ! Temporary variables.
  integer :: i,j,k,ialloc

  ! Construct sum[f'^x'].
  allocate( output(supercell%no_atoms,supercell%no_atoms),  &
          & updated(supercell%no_atoms,supercell%no_atoms), &
          & stat=ialloc); call err(ialloc)
  output = dblemat(zeroes(3,3))
  updated = .false.
  do i=1,size(supercell%symmetries)
    do j=1,size(unique_directions)
      atom_1 = supercell%atoms(unique_directions(j)%atom_id)
      atom_1p = supercell%atoms( supercell%symmetries(i)%atom_group &
                             & * atom_1%id())
      
      x = supercell%symmetries(i)%cartesian_tensor &
      & * unique_directions(j)%atomic_displacement
      
      do k=1,supercell%no_atoms
        atom_2 = supercell%atoms(k)
        atom_2p = supercell%atoms( supercell%symmetries(i)%atom_group &
                               & * atom_2%id())
        
        f = supercell%symmetries(i)%cartesian_tensor * forces(atom_2%id(),j)
        
        output(atom_2p%id(),atom_1p%id()) = output(atom_2p%id(),atom_1p%id()) &
                                        & + outer_product(f,x)
        updated(atom_2p%id(),atom_1p%id()) = .true.
      enddo
    enddo
  enddo
  
  ! Check that all elements of fx have been updated.
  if (any(.not. updated)) then
    call print_line(CODE_ERROR//': Some elements of f^x were not updated.')
    call err()
  endif
  
  ! Check sum[f'^x'] obeys correct symmetries.
  average = 0
  difference = 0
  do i=1,size(supercell%symmetries)
    do j=1,size(unique_directions)
      atom_1 = supercell%atoms(unique_directions(j)%atom_id)
      atom_1p = supercell%atoms( supercell%symmetries(i)%atom_group &
                             & * atom_1%id())
      do k=1,supercell%no_atoms
        atom_2 = supercell%atoms(k)
        atom_2p = supercell%atoms( supercell%symmetries(i)%atom_group &
                               & * atom_2%id())
        matrix = output(atom_1p%id(),atom_2p%id())
        symmetric = supercell%symmetries(i)%cartesian_tensor &
                & * output(atom_1%id(),atom_2%id())          &
                & * transpose(supercell%symmetries(i)%cartesian_tensor)
        average = average + sum_squares((matrix+symmetric)/2)
        difference = difference + sum_squares(matrix-symmetric)
      enddo
    enddo
  enddo
  call logfile%print_line(                                  &
     & 'Fractional L2 error in symmetry of f^x        : '// &
     & sqrt(difference/average))
  if (sqrt(difference/average) > 1e-10_dp) then
    call print_line(WARNING//': f^x is not as symmetric as expected. Please &
       &check log files.')
  endif
end function

! ----------------------------------------------------------------------
! Construct and check F = s * sum[f'^x'] . inverse(sum[x'^x']).
! ----------------------------------------------------------------------
! sum[x'^x'] is block diagonal, so only the 3x3 diagonal blocks need inverting.
! N.B. the factor of s (sc_size) is to correct for the fact that only one
!    atom in every s primitive cells is displaced, so the magnitude of x is
!    effectively sqrt(s) times shorter than the equivalent displacement in
!    the 1x1x1 supercell.
function construct_f(xx,fx,supercell,logfile) result(output)
  implicit none
  
  type(RealMatrix),    intent(in)    :: xx(:)
  type(RealMatrix),    intent(in)    :: fx(:,:)
  type(StructureData), intent(in)    :: supercell
  type(OFile),         intent(inout) :: logfile
  type(ForceConstants)               :: output
  
  type(RealMatrix) :: xx_inverse
  
  ! Variables for checking the output.
  type(AtomData)   :: atom_1
  type(AtomData)   :: atom_1p
  type(AtomData)   :: atom_2
  type(AtomData)   :: atom_2p
  type(RealMatrix) :: matrix
  type(RealMatrix) :: symmetric
  real(dp)         :: average
  real(dp)         :: difference
  
  integer :: i,j,k,ialloc
  
  ! Construct F(i1,i2).
  allocate(output%constants_(supercell%no_atoms,supercell%no_atoms), &
     & stat=ialloc); call err(ialloc)
  do i=1,supercell%no_atoms
    xx_inverse = invert(xx(i))
    do j=1,supercell%no_atoms
      output%constants_(j,i) = fx(j,i) * xx_inverse * supercell%sc_size
    enddo
  enddo
  
  ! Symmetrise F(i1,i2) under exchange of co-ordinates.
  do i=1,supercell%no_atoms
    do j=1,i
      output%constants_(j,i) = ( output%constants_(j,i)              &
                           &   + transpose(output%constants_(i,j)) ) &
                           & / 2.0_dp
      output%constants_(i,j) = transpose(output%constants_(j,i))
    enddo
  enddo
  
  ! Check F(i1,i2) transforms correctly under symmetry operators.
  average = 0
  difference = 0
  do i=1,size(supercell%symmetries)
    do j=1,supercell%no_atoms_prim
      atom_1 = supercell%atoms(j)
      atom_1p = supercell%atoms( supercell%symmetries(i)%atom_group &
                             & * atom_1%id())
      do k=1,supercell%no_atoms
        atom_2 = supercell%atoms(k)
        atom_2p = supercell%atoms( supercell%symmetries(i)%atom_group &
                               & * atom_2%id())
        matrix = output%constants(atom_1p,atom_2p)
        symmetric = supercell%symmetries(i)%cartesian_tensor &
                & * output%constants(atom_1,atom_2)          &
                & * transpose(supercell%symmetries(i)%cartesian_tensor)
        average = average + sum_squares((matrix+symmetric)/2)
        difference = difference + sum_squares(matrix-symmetric)
      enddo
    enddo
  enddo
  call logfile%print_line(                                  &
     & 'Fractional L2 error in symmetry of F(i1,i2)   : '// &
     & sqrt(difference/average))
  if (sqrt(difference/average) > 1e-10_dp) then
    call print_line(WARNING//': F(i1,i2) is not as symmetric as expected. &
       &Please check log files.')
  endif
  
  ! Check F(i1,i2)=F(i2,i1).
  average = 0
  difference = 0
  do i=1,supercell%no_atoms
    atom_1 = supercell%atoms(i)
    do j=1,supercell%no_atoms
      atom_2 = supercell%atoms(j)
      matrix = output%constants(atom_2,atom_1)
      symmetric = transpose(output%constants(atom_1,atom_2))
      average = average + sum_squares((matrix+symmetric)/2)
      difference = difference + sum_squares(matrix-symmetric)
    enddo
  enddo
  call logfile%print_line(                                  &
     & 'Fractional L2 error in F(i1,i2)=F(i2,i1)      : '// &
     & sqrt(difference/average))
  if (sqrt(difference/average) > 1e-10_dp) then
    call print_line(WARNING//': F(i1,i2)/=F(i2,i1). Please check log files.')
  endif
end function

! ----------------------------------------------------------------------
! Checks the force constants against the calculated forces and symmetries.
! ----------------------------------------------------------------------
! Checks that the elements corresponding to calculated forces have not been
!    changed too much by symmetrisation.
! Checks that the force constants have the correct symmetry properties.
subroutine check(this,forces,supercell,unique_directions,logfile)
  implicit none
  
  class(ForceConstants), intent(in)    :: this
  type(RealVector),      intent(in)    :: forces(:,:)
  type(StructureData),   intent(in)    :: supercell
  type(UniqueDirection), intent(in)    :: unique_directions(:)
  type(OFile),           intent(inout) :: logfile
  
  ! Atoms.
  type(AtomData) :: atom_1
  type(AtomData) :: atom_2
  
  real(dp) :: average
  real(dp) :: difference
  
  type(RealVector) :: calculated
  type(RealVector) :: fitted
  
  integer :: i,j
  
  ! Check force constants against raw forces.
  average = 0
  difference = 0
  do i=1,size(unique_directions)
    atom_1 = supercell%atoms(unique_directions(i)%atom_id)
    
    do j=1,supercell%no_atoms
      atom_2 = supercell%atoms(j)
      
      calculated = forces(j,i)
      fitted = this%constants(atom_2,atom_1)            &
           & * unique_directions(i)%atomic_displacement &
           & / supercell%sc_size
      
      average = average + sum_squares((calculated+fitted)/2)
      difference = difference + sum_squares(calculated-fitted)
    enddo
  enddo
  call logfile%print_line('Fractional L2 difference between forces &
     &before and after symmetrisation (this may be large): '// &
     &sqrt(difference/average))
end subroutine
end module
