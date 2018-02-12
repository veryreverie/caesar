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
    this%primitive_displacements(i,1) = cmplx(line)
    
    if (.not. this%at_paired_qpoint) then
      line = split(mode_file%line(6+no_atoms+i))
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

! Will lift degeneracies iff degenerate_energy is present.
! Calculate modes for one of the calculated q-points.
function calculate_modes_calculated(matrices,structure,qpoint, &
   &degenerate_energy,logfile) result(output)
  use structure_module
  use qpoints_module
  use atom_module
  use eigenstuff_module
  use ofile_module
  implicit none
  
  type(ComplexMatrix), intent(in)    :: matrices(:,:)
  type(StructureData), intent(in)    :: structure
  type(QpointData),    intent(in)    :: qpoint
  real(dp),            intent(in)    :: degenerate_energy
  type(OFile),         intent(inout) :: logfile
  type(ComplexMode), allocatable     :: output(:)
  
  complex(dp), allocatable :: dyn_mat(:,:)
  
  type(HermitianEigenstuff) :: estuff
  
  real(dp), allocatable :: frequencies(:)
  logical,  allocatable :: translational(:)
  
  type(ComplexVector)              :: displacement
  
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
  
  ! Lift degeneracies, expressing degenerate states in terms of
  !    the eigenvectors of symmetry operators.
  output = lift_degeneracies( output,            &
                            & structure,         &
                            & qpoint,            &
                            & degenerate_energy, &
                            & logfile)
end function

function lift_degeneracies(input,structure,qpoint,degenerate_energy,logfile) &
   & result(output)
  use structure_module
  use linear_algebra_module
  use group_module
  use symmetry_module
  use qpoints_module
  use logic_module
  use ofile_module
  implicit none
  
  type(ComplexMode),   intent(in)    :: input(:)
  type(StructureData), intent(in)    :: structure
  type(QpointData),    intent(in)    :: qpoint
  real(dp),            intent(in)    :: degenerate_energy
  type(OFile),         intent(inout) :: logfile
  type(ComplexMode), allocatable     :: output(:)
  
  ! The id of the first mode with which mode i is degenerate.
  ! e.g. if there are four modes, and the middle two are degenerate, then
  !    degeneracy_ids = [ 1, 2, 2, 4 ].
  integer, allocatable :: degeneracy_ids(:)
  
  logical, allocatable :: degenerate(:)
  
  real(dp) :: frequency
  real(dp) :: prev_frequency
  real(dp) :: prev_degenerate_frequency
  
  integer :: i,j,ialloc
  
  output = input
  
  ! Identify which modes are degenerate.
  allocate(degeneracy_ids(size(output)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    if (i==1) then
      degeneracy_ids(1) = 1
    else
      ! Find the frequency of this mode, the previous mode, and the first mode
      !    which is degenerate with the previous mode.
      frequency = output(i)%frequency
      prev_frequency = output(i-1)%frequency
      prev_degenerate_frequency = output(degeneracy_ids(i-1))%frequency
      if ( abs(frequency-prev_frequency)<degenerate_energy .and. &
         & abs(frequency-prev_degenerate_frequency)<degenerate_energy) then
        ! This mode is degenerage with the previous set.
        degeneracy_ids(i) = degeneracy_ids(i-1)
      elseif (abs(frequency-prev_frequency)<degenerate_energy) then
        ! This mode is degenerate with the previous mode, but not with the
        !    modes which that mode is degenerate with. This is an error.
        call print_line(ERROR//': Modes inconsistently degenerate. Try &
           &adjusting degenerate_energy.')
        call err()
      else
        ! This mode is not degenerate with the previous set.
        degeneracy_ids(i) = i
      endif
    endif
  enddo
  
  ! Lift the degeneracies.
  allocate(degenerate(size(output)), stat=ialloc); call err(ialloc)
  degenerate = .true.
  do i=1,size(output)
    ! Ignore this mode if it has already been dealt with.
    if (.not. degenerate(i)) then
      cycle
    endif
    
    ! Find the id of the last mode which is degenerate with mode i.
    j = last(degeneracy_ids==i)
    
    if (j==i) then
      degenerate(i) = .false.
    elseif (j>i) then
      output(i:j) = lift_degeneracies_2( output(i:j),          &
                                       & structure%symmetries, &
                                       & qpoint,               &
                                       & logfile)
      degenerate(i:j) = .false.
    else
      call err()
    endif
  enddo
end function

! --------------------------------------------------
! Helper function for lift_degeneracies.
! Recursively lifts degeneracies.
! --------------------------------------------------
! Input must be a list of degenerate modes.
recursive function lift_degeneracies_2(input,symmetries,qpoint,logfile) &
   &result(output)
  use utils_module, only : sum_squares
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
  
  type(ComplexMode), allocatable :: rotated_modes(:)
  complex(dp),       allocatable :: dot_products(:,:)
  
  type(UnitaryEigenstuff), allocatable :: estuff(:)
  
  integer :: sym_id
  
  integer :: order
  
  integer, allocatable :: eval_ids(:)
  integer, allocatable :: sort_ids(:)
  logical, allocatable :: done(:)
  
  real(dp) :: check
  
  integer :: no_modes
  
  integer :: i,j,ialloc
  
  ! All q-point data and eigenvalues will be unchanged.
  ! Copy over all data, and only change that which changes.
  output = input
  
  ! Count the number of modes.
  no_modes = size(input)
  if (no_modes==1) then
    return
  endif
  
  ! Find a symmetry which maps q onto itself.
  sym_id = 0
  do i=1,size(symmetries)
    if (symmetries(i)%recip_rotation*qpoint%qpoint == qpoint%qpoint) then
      sym_id = i
      exit
    endif
  enddo
  if (sym_id==0) then
    call print_line(ERROR//': Unable to lift degeneracies using symmetry.')
    call err()
  endif
  
  ! Calculate the order of the symmetry, n s.t. S^n=I.
  order = calculate_symmetry_order(symmetries(sym_id), qpoint)
  
  ! Construct the symmetry in normal mode co-ordinates.
  rotated_modes = rotate_complex_modes( input,              &
                                      & symmetries(sym_id), &
                                      & qpoint%qpoint,      &
                                      & qpoint%qpoint)
  allocate( dot_products(size(input),size(input)), &
          & stat=ialloc); call err(ialloc)
  do i=1,no_modes
    do j=1,no_modes
      dot_products(j,i) = input(j) * rotated_modes(i)
    enddo
  enddo
  
  ! Check that the symmetry only maps the degenerate modes onto each another.
  check = 0
  do i=1,no_modes
    check = check &
        & + abs(dot_product(dot_products(:,i),dot_products(:,i))-1)
  enddo
  call logfile%print_line('Fractional L2 error in rotation between modes: ' &
     & //sqrt(check))
  if (sqrt(check)>1e-6_dp) then
    call print_line(WARNING//': Rotations and degeneracies inconsistent. &
       &Please check log files.')
    call print_line('check: '//sqrt(check))
  endif
  
  ! Check that the symmetry is unitary.
  check = sqrt(sum_squares(mat(dot_products)*hermitian(mat(dot_products))-cmplxmat(make_identity_matrix(3))))
  
  ! Diagonalise the symmetry, and construct diagonalised displacements.
  ! Only transform displacements if this symmetry lifts degeneracy.
  estuff = diagonalise_unitary(dot_products,order)
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
      if (check>1e-8_dp) then
        call print_line(WARNING//': Error in mode normalisation. Please check &
           &log files.')
      endif
    enddo
  endif
  
  ! Lift remaining degeneracies recursively.
  i = 1
  do while(i<=no_modes)
    j = 0
    do j=no_modes,i,-1
      if (estuff(j)%eval==estuff(i)%eval) then
        exit
      endif
    enddo
    
    if (j<i) then
      call err()
    endif
    
    if (j>i) then
      output(i:j) = lift_degeneracies_2( output(i:j),           &
                                       & symmetries(sym_id+1:), &
                                       & qpoint,                &
                                       & logfile)
    endif
    
    i = j+1
  enddo
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
  
  integer :: i
  
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
  
  integer :: ialloc
  
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
