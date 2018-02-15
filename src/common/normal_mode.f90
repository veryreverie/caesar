! ======================================================================
! Harmonic normal modes.
! ======================================================================
!
! Most calculation happens in normal mode co-ordinates, but setting up
!    DFT calculations requires conversion to cartesian, and the results of
!    DFT have to be converted from cartesian forces.
!
! - Co-ordinates can be converted between real and normal modes via the
!      relevant class constructors.
! - Displacements in real mode co-ordinates can be converted into cartesian
!      co-ordinates using real_mode_to_displacement.
! - Forces can be converted into real mode co-ordinates using
!      force_to_real_mode.
!
! d(R,j) = sqrt(m_j)(r(R,j)-r0(R,j))
! f(R,j) = -dV/dd(R,j)
!
! x(q,j) = 1/N * sum_R [ d(R,j) * exp( i q.R) ]
! d(R,j) =       sum_q [ x(q,j) * exp(-i q.R) ]
!
! N.B. if 2q=G then q.R is a half-integer.
!    => exp(i q.R) = cos(q.R) = +/- 1
!    => x(q,j) is real for all j.
!
! x(q,j) = 1/N       * sum_R [ d(R,j) * cos(q.R) ]           (for 2q=G only)
! c(q,j) = sqrt(2)/N * sum_R [ d(R,j) * cos(q.R) ]           (for 2q/=G)
! s(q,j) = sqrt(2)/N * sum_R [ d(R,j) * sin(q.R) ]           (for 2q/=G)
! d(R,j) =           sum_q [x_q * cos(q,R)]                  (sum over 2q=G)
!        + sqrt(2) * sum_q [c_q * cos(q.R) + s_q * sin(q.R)] (sum over 2q/=G)
!
! ux(q,k) = sum_j [ U(q,j,k) * x(q,j) ]        (for 2q=G only)
! uc(q,k) = sum_j [ U(q,j,k) * c(q,j) ]        (for 2q/=G)
! us(q,k) = sum_j [ U(q,j,k) * s(q,j) ]        (for 2q/=G)
module normal_mode_module
  use constants_module, only : dp
  use string_module
  use io_module
  use linear_algebra_module
  implicit none
  
  private
  
  public :: ComplexMode
  public :: RealMode
  
  public :: calculate_modes
  public :: rotate_complex_modes
  
  public :: complex_to_real
  public :: real_to_complex
  
  ! --------------------------------------------------
  ! A normal mode in complex co-ordinates.
  ! --------------------------------------------------
  type ComplexMode
    ! Whether or not 2q=G. If true, there is only one mode, and it is real.
    !    If not, there are two modes, and they are conjugates of one another.
    logical :: at_paired_qpoint
    
    ! The frequency, and frequency-relevant information.
    real(dp) :: frequency
    logical  :: soft_mode          ! True if frequency < -epsilon.
    logical  :: translational_mode ! True if frequency=0 and at gamma.
    
    ! The displacements of atoms in the primitive cell.
    ! An n x 1 array of vectors if at_paired_qpoint, or an n x 2 array if not.
    ! The first column corresponds to e^(iq.R), the second to e^(-iq.R).
    type(ComplexVector), allocatable :: primitive_displacements(:,:)
    
    ! An id which is shared between degenerate states, and different otherwise.
    integer :: degeneracy_id
  contains
    procedure, public :: write_file => write_file_ComplexMode
  end type
  
  interface ComplexMode
    module procedure read_file_ComplexMode
  end interface
  
  interface operator(*)
    module procedure dot_ComplexModes
  end interface
  
  interface l2_norm
    module procedure l2_norm_ComplexMode
  end interface
  
  ! --------------------------------------------------
  ! A normal mode in real co-ordinates.
  ! --------------------------------------------------
  type RealMode
    ! Whether or not 2q=G. If true, there is only the cosine mode.
    !    If not, there is a cosine mode and a sine mode.
    logical :: at_paired_qpoint
    
    ! The frequency, and frequency-relevant information.
    real(dp) :: frequency
    logical  :: soft_mode          ! True if frequency < -epsilon.
    logical  :: translational_mode ! True if frequency=0 and at gamma.
    
    ! The displacements of atoms in the primitive cell.
    ! An n x 1 array of vectors if at_paired_qpoint, or an n x 2 array if not.
    ! The first column corresponds to cos(q.R), the second to sin(q.R).
    type(RealVector), allocatable :: primitive_displacements(:,:)
  contains
    procedure, public :: write_file => write_file_RealMode
  end type
  
  interface RealMode
    module procedure read_file_RealMode
  end interface
  
  interface calculate_modes
    module procedure calculate_modes_interpolated
    module procedure calculate_modes_calculated
  end interface
  
  ! --------------------------------------------------
  ! Conversions between complex and real co-ordinates.
  ! --------------------------------------------------
  interface complex_to_real
    module procedure complex_to_real_Mode
  end interface
  
  interface real_to_complex
    module procedure real_to_complex_Mode
  end interface
contains

! ----------------------------------------------------------------------
! ComplexMode procedures.
! ----------------------------------------------------------------------
subroutine write_file_ComplexMode(this,filename)
  use ofile_module
  implicit none
  
  class(ComplexMode), intent(in) :: this
  type(String),       intent(in) :: filename
  
  type(OFile) :: mode_file
  
  integer :: i
  
  mode_file = filename
  call mode_file%print_line('Mode at paired q-point: '//this%at_paired_qpoint)
  call mode_file%print_line('Frequency:              '//this%frequency)
  call mode_file%print_line('Soft mode:              '//this%soft_mode)
  call mode_file%print_line('Translational mode:     '//this%translational_mode)
  call mode_file%print_line('Degeneracy id:          '//this%degeneracy_id)
  call mode_file%print_line('Primitive Displacements of e^(iq.R) mode:')
  do i=1,size(this%primitive_displacements,1)
    call mode_file%print_line(this%primitive_displacements(i,1))
  enddo
  if (.not. this%at_paired_qpoint) then
    call mode_file%print_line('Primitive Displacements of e^(-iq.R) mode:')
    do i=1,size(this%primitive_displacements,1)
      call mode_file%print_line(this%primitive_displacements(i,2))
    enddo
  endif
end subroutine

function read_file_ComplexMode(filename) result(this)
  use ifile_module
  implicit none
  
  type(String), intent(in) :: filename
  type(ComplexMode)        :: this
  
  type(IFile)               :: mode_file
  type(String), allocatable :: line(:)
  integer                   :: no_atoms
  
  integer :: i,ialloc
  
  ! Read in mode file.
  mode_file = filename
  
  ! Read whether or not mode is at a q-point such that 2q=G.
  line = split(mode_file%line(1))
  this%at_paired_qpoint = lgcl(line(5))
  
  ! Read frequency.
  line = split(mode_file%line(2))
  this%frequency = dble(line(2))
  
  ! Read whether or not this mode is soft.
  line = split(mode_file%line(3))
  this%soft_mode = lgcl(line(3))
  
  ! Read whether or not this mode is purely translational.
  line = split(mode_file%line(4))
  this%translational_mode = lgcl(line(3))
  
  ! Read the degeneracy id of this mode.
  line = split(mode_file%line(5))
  this%degeneracy_id = int(line(3))
  
  ! Read in the displacement associated with the mode.
  if (this%at_paired_qpoint) then
    no_atoms = size(mode_file)-6
    allocate( this%primitive_displacements(no_atoms,1), &
            & stat=ialloc); call err(ialloc)
  else
    no_atoms = (size(mode_file)-6)/2
    allocate( this%primitive_displacements(no_atoms,2), &
            & stat=ialloc); call err(ialloc)
  endif
  
  do i=1,no_atoms
    line = split(mode_file%line(6+i))
    this%primitive_displacements(i,1) = cmplx(line)
    
    if (.not. this%at_paired_qpoint) then
      line = split(mode_file%line(7+no_atoms+i))
      this%primitive_displacements(i,2) = cmplx(line)
    endif
  enddo
end function

! ----------------------------------------------------------------------
! RealMode procedures.
! ----------------------------------------------------------------------
subroutine write_file_RealMode(this,filename)
  use ofile_module
  implicit none
  
  class(RealMode), intent(in) :: this
  type(String),    intent(in) :: filename
  
  type(OFile) :: mode_file
  
  integer :: i
  
  mode_file = filename
  call mode_file%print_line('Mode at paired q-point: '//this%at_paired_qpoint)
  call mode_file%print_line('Frequency:              '//this%frequency)
  call mode_file%print_line('Soft mode:              '//this%soft_mode)
  call mode_file%print_line('Translational mode:     '//this%translational_mode)
  call mode_file%print_line('Primitive Displacements of cos(q.R) mode:')
  do i=1,size(this%primitive_displacements)
    call mode_file%print_line(this%primitive_displacements(i,1))
  enddo
  if (.not. this%at_paired_qpoint) then
    call mode_file%print_line('Primitive Displacements of sin(q.R) mode:')
    do i=1,size(this%primitive_displacements)
      call mode_file%print_line(this%primitive_displacements(i,2))
    enddo
  endif
end subroutine

function read_file_RealMode(filename) result(this)
  use ifile_module
  implicit none
  
  type(String), intent(in) :: filename
  type(RealMode)           :: this
  
  type(IFile)               :: mode_file
  type(String), allocatable :: line(:)
  integer                   :: no_atoms
  
  integer :: i,ialloc
  
  ! Read in mode file.
  mode_file = filename
  
  ! Read kind.
  line = split(mode_file%line(1))
  this%at_paired_qpoint = lgcl(line(5))
  
  ! Read frequency.
  line = split(mode_file%line(2))
  this%frequency = dble(line(2))
  
  ! Read whether or not this mode is soft.
  line = split(mode_file%line(3))
  this%soft_mode = lgcl(line(3))
  
  ! Read whether or not this mode is purely translational.
  line = split(mode_file%line(4))
  this%translational_mode = lgcl(line(3))
  
  ! Read in the displacement associated with the mode.
  if (this%at_paired_qpoint) then
    no_atoms = size(mode_file)-5
    allocate( this%primitive_displacements(no_atoms,1), &
            & stat=ialloc); call err(ialloc)
  else
    no_atoms = (size(mode_file)-5)/2
    allocate( this%primitive_displacements(no_atoms,2), &
            & stat=ialloc); call err(ialloc)
  endif
  do i=1,no_atoms
    line = split(mode_file%line(5+i))
    this%primitive_displacements(i,1) = dble(line)
    if (.not. this%at_paired_qpoint) then
      line = split(mode_file%line(6+no_atoms+i))
      this%primitive_displacements(i,2) = dble(line)
    endif
  enddo
end function

! ----------------------------------------------------------------------
! Calculates complex modes by diagonalising a dynamical matrix.
! ----------------------------------------------------------------------
! N.B. Structure may be any supercell.

! Calculate modes for a q-point other than one of the calculated q-points.
function calculate_modes_interpolated(matrices,structure) result(output)
  use structure_module
  use atom_module
  use eigenstuff_module
  implicit none
  
  type(ComplexMatrix), intent(in) :: matrices(:,:)
  type(StructureData), intent(in) :: structure
  type(ComplexMode), allocatable  :: output(:)
  
  complex(dp), allocatable :: dyn_mat(:,:)
  
  type(HermitianEigenstuff), allocatable :: estuff(:)
  
  type(ComplexVector) :: displacement
  
  integer :: i,j,k,ialloc
  
  ! Convert (3x3Matrix) x no_atoms x no_atoms to no_modes x no_modes
  allocate( dyn_mat(structure%no_modes_prim,structure%no_modes_prim), &
          & stat=ialloc); call err(ialloc)
  do i=1,structure%no_atoms_prim
    do j=1,structure%no_atoms_prim
      dyn_mat(3*j-2:3*j, 3*i-2:3*i) = cmplx(matrices(j,i))
    enddo
  enddo
  
  ! Diagonalise dynamical matrix.
  estuff = diagonalise_hermitian(dyn_mat)
  
  ! Calculate normal mode frequencies and displacements.
  !          V = sum_i[ 0.5 * freq_i**2 * u_i**2]
  ! -> F = -2V = sum_i[ - freq_i**2 * u_i**2 ]
  allocate( output(structure%no_modes_prim),   &
          & stat=ialloc); call err(ialloc)
  do i=1,structure%no_modes_prim
    
    ! The eigenvalues are in descending order, but the normal modes should
    !    be in ascending order of frequency. i->k reverses the order.
    k = structure%no_modes_prim - i + 1
    
    if (estuff(k)%eval>=0.0_dp) then
      ! Unstable mode.
      output(i)%frequency = - sqrt(estuff(k)%eval)
    else
      ! Stable mode.
      output(i)%frequency = sqrt(- estuff(k)%eval)
    endif
    
    output(i)%soft_mode = output(i)%frequency < -1.0e-6_dp
    output(i)%translational_mode = .false.
    output(i)%at_paired_qpoint = .false.
    
    ! Calculate displacements in the primitive cell,
    !    which are the non-mass-reduced eigenvectors of the dynamical matrix.
    allocate( output(i)%primitive_displacements(structure%no_atoms_prim,2), &
            & stat=ialloc); call err(ialloc)
    do j=1,structure%no_atoms_prim
      displacement = estuff(k)%evec(3*j-2:3*j) &
                 & * sqrt(structure%atoms(j)%mass())
      
      output(i)%primitive_displacements(j,:) = [ displacement, &
                                                 & conjg(displacement) ]
    enddo
    
    ! Re-normalise modes, now in non-mass-reduced co-ordinates.
    output(i)%primitive_displacements = output(i)%primitive_displacements &
                                    & / l2_norm(output(i))
  enddo
end function

! Calculate modes for one of the calculated q-points.
! Will lift degeneracies using symmetries.
function calculate_modes_calculated(matrices,structure,qpoint, &
   &degenerate_energy,degeneracy_id,logfile) result(output)
  use utils_module, only : sum_squares
  use structure_module
  use qpoints_module
  use atom_module
  use eigenstuff_module
  use ofile_module
  use symmetry_module
  use logic_module
  implicit none
  
  type(ComplexMatrix), intent(in)    :: matrices(:,:)
  type(StructureData), intent(in)    :: structure
  type(QpointData),    intent(in)    :: qpoint
  real(dp),            intent(in)    :: degenerate_energy
  integer,             intent(in)    :: degeneracy_id
  type(OFile),         intent(inout) :: logfile
  type(ComplexMode), allocatable     :: output(:)
  
  real(dp), allocatable :: frequencies(:)
  logical,  allocatable :: translational(:)
  
  ! Symmetry data.
  type(SymmetryOperator), allocatable :: symmetries(:)
  
  integer, allocatable :: states(:)
  
  ! Error checking variables.
  type(ComplexMatrix) :: symmetry
  real(dp)            :: check
  
  integer :: i,j,k,ialloc
  
  ! Calculate normal modes as if at an arbitrary q-point.
  output = calculate_modes(matrices,structure)
  
  ! Identify purely translational modes (at the gamma-point only).
  if (is_int(qpoint%qpoint)) then
    allocate( translational(structure%no_modes_prim), &
            & frequencies(structure%no_modes_prim),   &
            & stat=ialloc); call err(ialloc)
    translational = .false.
    do i=1,structure%no_modes_prim
      frequencies(i) = output(i)%frequency
    enddo
    do i=1,3
      j = minloc(abs(frequencies),dim=1,mask=.not.translational)
      translational(j) = .true.
      output(j)%translational_mode = .true.
    enddo
  endif
  
  ! Pair down complex modes if 2q=G, i.e. if the modes are real.
  if (is_int(2*qpoint%qpoint)) then
    do i=1,structure%no_modes_prim
      output(i)%at_paired_qpoint = .true.
      
      output(i)%primitive_displacements = &
         & output(i)%primitive_displacements(:,1:1)
    enddo
  endif
  
  ! Identify the symmetries which map the q-point to itself.
  symmetries = structure%symmetries( filter( structure%symmetries, &
                                           & leaves_q_invariant) )
  
  ! Assign degeneracy ids, which are equal if two states are degenerate,
  !    and different if they are not.
  output(1)%degeneracy_id = degeneracy_id
  do i=2,size(output)
    if (abs(output(i)%frequency-output(i-1)%frequency)<degenerate_energy) then
      output(i)%degeneracy_id = output(i-1)%degeneracy_id
    else
      output(i)%degeneracy_id = output(i-1)%degeneracy_id + 1
    endif
  enddo
  
  call print_line('')
  call print_line('q-point        : '//qpoint%qpoint)
  call print_line('Degeneracy ids : '//output%degeneracy_id)
  
  ! Loop over degeneracy ids, checking each degenerate subspace, and lifting
  !    degeneracy using symmetry operators.
  do i=degeneracy_id,output(size(output))%degeneracy_id
    ! Find the set of states with degeneracy id i.
    states = filter(output%degeneracy_id==i)
    
    call print_line('States: '//states)
    
    ! Check that degenerate states are consistent, i.e. that if two states
    !    are both degenerate with a third state that they are also degenerate
    !    with one another.
    if ( maxval(output(states)%frequency) - minval(output(states)%frequency) &
     & > degenerate_energy ) then
      call print_line(ERROR//': Modes inconsistently degenerate. Please try &
         &adjusting degenerate_energy.')
      call err()
    endif
    
    ! Set the frequencies of each degenerate state to
    !    the average of their frequencies.
    output(states)%frequency = sum(output(states)%frequency) / size(states)
    
    ! Check that all symmetries map the degenerate subspace onto itself.
    do j=1,size(symmetries)
      symmetry = calculate_symmetry_in_normal_coordinates( qpoint,         &
                                                         & output(states), &
                                                         & qpoint,         &
                                                         & output(states), &
                                                         & symmetries(j),  &
                                                         & logfile)
      check = sqrt(sum_squares( symmetry*hermitian(symmetry) &
                            & - cmplxmat(make_identity_matrix(size(states)))))
      if (check>1e-2_dp) then
        call print_line(ERROR//': Rotation between degenerate modes at &
           &q-point '//qpoint%qpoint//' not unitary. Please adjust &
           &degenerate_energy so that it accurately captures degeneracies.')
        stop
      endif
    enddo
    
    ! Lift each degeneracy using symmetry operators.
    if (size(states)>1) then
      output(states) = lift_degeneracies( output(states), &
                                        & symmetries,     &
                                        & qpoint,         &
                                        & logfile)
    endif
  enddo
contains
  ! Lambda for identifying if a symmetry leaves the q-point invariant.
  ! Captures qpoint.
  function leaves_q_invariant(input) result(output)
    implicit none
    
    class(*), intent(in) :: input
    logical              :: output
    
    select type(input); type is(SymmetryOperator)
      output = input%recip_rotation * qpoint%qpoint == qpoint%qpoint
    end select
  end function
end function

! --------------------------------------------------
! Helper function for lift_degeneracies.
! Recursively lifts degeneracies.
! --------------------------------------------------
! Input must be a list of degenerate modes.
! Symmetries must all take the q-point to itself.
recursive function lift_degeneracies(input,symmetries,qpoint,logfile) &
   & result(output)
  use symmetry_module
  use qpoints_module
  use eigenstuff_module
  use phase_module
  use logic_module
  use ofile_module
  implicit none
  
  type(ComplexMode),      intent(in)    :: input(:)
  type(SymmetryOperator), intent(in)    :: symmetries(:)
  type(QpointData),       intent(in)    :: qpoint
  type(OFile),            intent(inout) :: logfile
  type(ComplexMode), allocatable        :: output(:)
  
  integer :: no_modes
  
  ! Symmetry information.
  integer :: sym_id
  integer :: order
  
  type(SymmetryOperator) :: first_symmetry
  type(ComplexMatrix) :: symmetry
  
  type(UnitaryEigenstuff), allocatable :: estuff(:)
  
  type(SymmetryOperator), allocatable :: commuting_symmetries(:)
  
  ! Error checking.
  real(dp) :: check
  
  integer :: i,j,ialloc
  
  ! All q-point data and eigenvalues will be unchanged.
  ! Copy over all data, and only change that which changes.
  output = input
  
  ! Count the number of modes.
  no_modes = size(input)
  
  if (no_modes==1) then
    call print_line(CODE_ERROR//': Trying to lift the degeneracy of only one &
       &state.')
    call err()
  endif
  
  if (size(symmetries)==0) then
    call print_line(ERROR//': Unable to lift degeneracies using symmetry. &
       &Please try reducing degenerate_energy.')
    stop
  endif

  first_symmetry = symmetries(1)
  
  ! Construct the first symmetry in normal mode co-ordinates.
  symmetry = calculate_symmetry_in_normal_coordinates( qpoint,         &
                                                     & input,          &
                                                     & qpoint,         &
                                                     & input,          &
                                                     & first_symmetry, &
                                                     & logfile)
  
  ! Calculate the order of the first symmetry, n s.t. S^n=I.
  order = calculate_symmetry_order(first_symmetry, qpoint)
  
  ! Diagonalise the first symmetry, and construct diagonalised displacements.
  ! Only transform displacements if this symmetry lifts degeneracy.
  estuff = diagonalise_unitary(symmetry,order)
  if (estuff(1)%eval/=estuff(no_modes)%eval) then
    do i=1,no_modes
      output(i)%primitive_displacements = cmplxvec(zeroes(3))
      do j=1,no_modes
        output(i)%primitive_displacements = output(i)%primitive_displacements &
                                        & + estuff(i)%evec(j)                 &
                                        & * input(j)%primitive_displacements
      enddo
      
      check = abs(l2_norm(output(i))-1)
      call logfile%print_line('Error in mode normalisation: '// &
         & check)
      if (check>1e-10_dp) then
        call print_line(WARNING//': Error in mode normalisation. Please check &
           &log files.')
      endif
    enddo
  endif
  
  ! Lift remaining degeneracies using remaining symmetries.
  i = 1
  do while(i<=no_modes)
    ! The range i:j is the set of modes degenerate with mode i under the
    !    first symmetry.
    j = last(estuff%eval==estuff(i)%eval)
    
    if (j<i) then
      call err()
    endif
    
    if (j>i) then
      ! Select only the symmetry operators which commute with the
      !    first symmetry.
      commuting_symmetries = symmetries(filter(symmetries,commutes_with_first))
      
      ! Lift further degeneracies using symmetries which commute with the
      !    first symmetry, not including the first symmetry.
      output(i:j) = lift_degeneracies( output(i:j),              &
                                     & commuting_symmetries(2:), &
                                     & qpoint,                   &
                                     & logfile)
    endif
    i = j+1
  enddo
contains
  ! Lambda for determining whether or not a symmetry commutes with the first
  !    symmetry.
  ! Captures first_symmetry.
  function commutes_with_first(input) result(output)
    implicit none
    
    class(*), intent(in) :: input
    logical              :: output
    
    select type(input); type is(SymmetryOperator)
      output = operators_commute(input,first_symmetry)
    end select
  end function
end function

! --------------------------------------------------
! Calculates a symmetry in normal mode co-ordinates.
! --------------------------------------------------
! Takes q1, {u1}, q2, {u2} and S. Outputs {u2.S.u1}.
function calculate_symmetry_in_normal_coordinates(qpoint_from,modes_from, &
   & qpoint_to,modes_to,symmetry,logfile) result(output)
  use utils_module, only : sum_squares
  use qpoints_module
  use symmetry_module
  use ofile_module
  implicit none
  
  type(QpointData),       intent(in)    :: qpoint_from
  type(ComplexMode),      intent(in)    :: modes_from(:)
  type(QpointData),       intent(in)    :: qpoint_to
  type(ComplexMode),      intent(in)    :: modes_to(:)
  type(SymmetryOperator), intent(in)    :: symmetry
  type(OFile),            intent(inout) :: logfile
  type(ComplexMatrix)                   :: output
  
  integer :: no_modes
  
  type(ComplexMode), allocatable :: rotated_modes(:)
  complex(dp),       allocatable :: dot_products(:,:)
  
  real(dp) :: check
  
  integer :: i,j,ialloc
  
  no_modes = size(modes_from)
  if (no_modes/=size(modes_to)) then
    call print_line(CODE_ERROR//': Inconsistent number of modes.')
    call err()
  endif
  
  ! Construct the symmetry in normal mode co-ordinates.
  rotated_modes = rotate_complex_modes( modes_from,         &
                                      & symmetry,           &
                                      & qpoint_from%qpoint, &
                                      & qpoint_to%qpoint)
  
  ! Construct the overlap matrix, u2.S.u1.
  allocate( dot_products(no_modes,no_modes), &
          & stat=ialloc); call err(ialloc)
  do i=1,no_modes
    do j=1,no_modes
      dot_products(j,i) = modes_to(j) * rotated_modes(i)
    enddo
  enddo
  
  output = dot_products
  
  ! Check that the symmetry is unitary.
  check = sqrt(sum_squares( output*hermitian(output) &
                        & - cmplxmat(make_identity_matrix(no_modes))))
  call logfile%print_line('Error in unitarity of rotation: ' &
     & //check)
  if (check>1e-10_dp) then
    call print_line(WARNING//': Rotation between degenerate modes not &
       &unitary. Please try adjusting degenerate_energy. Please check log &
       &files.')
  endif
end function

! ----------------------------------------------------------------------
! Calculates the order of a symmetry operation at a give q-point.
! The order is the smallest integer n>0 s.t. S^n=I, where I is the identity.
! ----------------------------------------------------------------------
! N.B. This calculation assumes that the symmetry changes relative phases.
!    If this is not the case, then n will be too large.
function calculate_symmetry_order(symmetry,qpoint) result(output)
  use utils_module, only : lcm
  use symmetry_module
  use qpoints_module
  implicit none
  
  type(SymmetryOperator), intent(in) :: symmetry
  type(QpointData),       intent(in) :: qpoint
  integer                            :: output
  
  type(IntMatrix) :: identity
  type(IntMatrix) :: rotation ! R^n.
  
  integer :: i
  
  identity = make_identity_matrix(3)
  rotation = identity
  
  ! Calculate n s.t. R^n=I.
  output = 0
  do i=1,6
    rotation = symmetry%rotation * rotation
    if (rotation==identity) then
      output = i
      exit
    endif
  enddo
  
  if (output==0) then
    call print_line(CODE_ERROR//': Unable to find order of symmetry.')
    call err()
  endif
  
  ! Assume that the symmetry changes relative phases.
  output = lcm(output,qpoint%min_sc_size())
end function

! ----------------------------------------------------------------------
! Find the dot product between two complex modes.
! ----------------------------------------------------------------------
function dot_ComplexModes(this,that) result(output)
  implicit none
  
  type(ComplexMode), intent(in) :: this
  type(ComplexMode), intent(in) :: that
  complex(dp)                   :: output
  
  output = sum( this%primitive_displacements(:,1) &
            & * conjg(that%primitive_displacements(:,1)) )
end function

function l2_norm_ComplexMode(this) result(output)
  implicit none
  
  type(ComplexMode), intent(in) :: this
  real(dp)                      :: output
  
  output = sqrt(real(this*this))
end function

! ----------------------------------------------------------------------
! Rotate complex modes from one q-point to another.
! ----------------------------------------------------------------------
function rotate_complex_modes(input,symmetry,qpoint_from,qpoint_to) &
   & result(output)
  use utils_module, only : exp_2pii
  use symmetry_module
  use qpoints_module
  implicit none
  
  type(ComplexMode),      intent(in) :: input(:)
  type(SymmetryOperator), intent(in) :: symmetry
  type(FractionVector),   intent(in) :: qpoint_from
  type(FractionVector),   intent(in) :: qpoint_to
  type(ComplexMode), allocatable     :: output(:)
  
  type(FractionVector) :: q
  type(IntVector)      :: r
  
  integer :: no_atoms
  
  integer :: mode
  integer :: atom_1
  integer :: atom_1p
  
  no_atoms = size(input(1)%primitive_displacements,1)
  q = qpoint_to
  
  ! Check that the symmetry rotates the q-point as expected.
  if (symmetry%recip_rotation * qpoint_from /= qpoint_to) then
    call print_line(CODE_ERROR//': Symmetry does not transform q-points as &
       &expected.')
    call err()
  endif
  
  ! Allocate output, and transfer across all data.
  ! (Displacements need rotating, but everything else stays the same.)
  output = input
  
  ! Rotate displacements.
  do mode=1,size(input)
    do atom_1=1,no_atoms
      atom_1p = symmetry%atom_group * atom_1
      r = symmetry%rvector(atom_1)
      if (input(1)%at_paired_qpoint) then
        output(mode)%primitive_displacements(atom_1p,1) =    &
           & symmetry%cartesian_rotation                     &
           & * input(mode)%primitive_displacements(atom_1,1) &
           & * exp_2pii(q*r)
      else
        output(mode)%primitive_displacements(atom_1p,1) =    &
           & symmetry%cartesian_rotation                     &
           & * input(mode)%primitive_displacements(atom_1,1) &
           & * exp_2pii(q*r)
        output(mode)%primitive_displacements(atom_1p,2) =    &
           & symmetry%cartesian_rotation                     &
           & * input(mode)%primitive_displacements(atom_1,2) &
           & * exp_2pii(q*r)
      endif
    enddo
  enddo
end function

function complex_to_real_Mode(input) result(output)
  implicit none
  
  type(ComplexMode), intent(in) :: input
  type(RealMode)                :: output
  
  integer :: no_atoms_prim
  
  integer :: i,ialloc
  
  ! Copy over common information.
  output%at_paired_qpoint   = input%at_paired_qpoint
  output%frequency          = input%frequency
  output%soft_mode          = input%soft_mode
  output%translational_mode = input%translational_mode
  
  no_atoms_prim = size(input%primitive_displacements,1)
  
  ! Convert displacements.
  if (input%at_paired_qpoint) then
    if (size(input%primitive_displacements,2)/=1) then
      call print_line(CODE_ERROR//': a normal mode at a qpoint &
         &where 2q=G has a paired mode.')
      call err()
    endif
    allocate( output%primitive_displacements(no_atoms_prim,1), &
            & stat=ialloc); call err(ialloc)
    ! x = x
    do i=1,no_atoms_prim
      output%primitive_displacements(i,1) = &
         & real(input%primitive_displacements(i,1))
    enddo
  else
    if (size(input%primitive_displacements,2)/=2) then
      call print_line(CODE_ERROR//': a normal mode at a qpoint &
         &where 2q/=G does not have a paired mode.')
      call err()
    endif
    allocate( output%primitive_displacements(no_atoms_prim,2), &
            & stat=ialloc); call err(ialloc)
    do i=1,no_atoms_prim
      ! c = (x+ + x-)/sqrt(2) = real(x+) * sqrt(2)
      output%primitive_displacements(i,1) = &
         & real(input%primitive_displacements(i,1)) * sqrt(2.0_dp)
      ! s = (x+ - x-)/(sqrt(2)i) = imag(x+) * sqrt(2)
      output%primitive_displacements(i,2) = &
         & aimag(input%primitive_displacements(i,1)) * sqrt(2.0_dp)
    enddo
  endif
end function

function real_to_complex_Mode(input) result(output)
  implicit none
  
  type(RealMode), intent(in) :: input
  type(ComplexMode)          :: output
  
  integer :: no_atoms_prim
  
  integer :: i,ialloc
  
  ! Copy over common information.
  output%at_paired_qpoint   = input%at_paired_qpoint
  output%frequency          = input%frequency
  output%soft_mode          = input%soft_mode
  output%translational_mode = input%translational_mode
  
  no_atoms_prim = size(input%primitive_displacements,1)
  
  ! Convert displacements.
  if (input%at_paired_qpoint) then
    if (size(input%primitive_displacements,2)/=1) then
      call print_line(CODE_ERROR//': a normal mode at a qpoint &
         &where 2q=G has a paired mode.')
      call err()
    endif
    allocate( output%primitive_displacements(no_atoms_prim,1), &
            & stat=ialloc); call err(ialloc)
    do i=1,no_atoms_prim
      ! x = x
      output%primitive_displacements(i,1) = &
         & cmplxvec(input%primitive_displacements(i,1))
    enddo
  else
    if (size(input%primitive_displacements,2)/=2) then
      call print_line(CODE_ERROR//': a normal mode at a qpoint &
         &where 2q/=G does not have a paired mode.')
      call err()
    endif
    allocate( output%primitive_displacements(no_atoms_prim,2), &
            & stat=ialloc); call err(ialloc)
    do i=1,no_atoms_prim
      ! x+ = (c+is)/sqrt(2)
      output%primitive_displacements(i,1) =               &
         & cmplxvec(  input%primitive_displacements(i,1), &
         &            input%primitive_displacements(i,2)) &
         & / sqrt(2.0_dp)
      ! x- = (c-is)/sqrt(2)
      output%primitive_displacements(i,2) =               &
         & cmplxvec(  input%primitive_displacements(i,1), &
         &           -input%primitive_displacements(i,2)) &
         & / sqrt(2.0_dp)
    enddo
  endif
end function
end module
