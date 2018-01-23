! ======================================================================
! The forces between atoms at a given q-point.
! ======================================================================
module dynamical_matrix_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  use linear_algebra_module
  use normal_mode_module
  implicit none
  
  public :: DynamicalMatrix
  public :: write_dynamical_matrix_file
  public :: read_dynamical_matrix_file
  public :: compare_dynamical_matrices
  
  type :: DynamicalMatrix
    type(ComplexMatrix), allocatable :: matrices(:,:)
    type(ComplexMode),   allocatable :: complex_modes(:)
  contains
    procedure, public :: check
    procedure, public :: frequencies
  end type
  
  interface DynamicalMatrix
    module procedure new_DynamicalMatrix
    module procedure new_DynamicalMatrix_calculated
    module procedure new_DynamicalMatrix_interpolated
  end interface
  
  interface conjg
    module procedure conjg_DynamicalMatrix
  end interface
contains

! ----------------------------------------------------------------------
! Constructs the matrix of force constants in q-point co-ordinates,
!    given the matrix of force constants in supercell co-ordinates.
! ----------------------------------------------------------------------
function new_DynamicalMatrix(q,at_gamma,paired_qpoint,supercell, &
   & force_constants,min_images) result(this)
  use constants_module, only : pi
  use structure_module
  use min_images_module
  use force_constants_module
  use atom_module
  use ofile_module
  implicit none
  
  type(RealVector),     intent(in)    :: q
  logical,              intent(in)    :: at_gamma
  logical,              intent(in)    :: paired_qpoint
  type(StructureData),  intent(in)    :: supercell
  type(ForceConstants), intent(in)    :: force_constants
  type(MinImages),      intent(in)    :: min_images(:,:)
  type(DynamicalMatrix)               :: this
  
  type(AtomData) :: atom_1
  type(AtomData) :: atom_2
  
  
  type(IntVector), allocatable :: rvectors(:)
  real(dp)                     :: qr
  complex(dp)                  :: exp_minus_iqr
  
  integer :: i,j,k,ialloc
  
  ! Evaluate the dynamical matrix.
  allocate( this%matrices( supercell%no_atoms_prim,  &
          &                supercell%no_atoms_prim), &
          & stat=ialloc); call err(ialloc)
  this%matrices = mat(cmplx(zeroes(3,3)))
  do i=1,supercell%no_atoms_prim
    atom_1 = supercell%atoms(i)
    do j=1,supercell%no_atoms
      atom_2 = supercell%atoms(j)
      rvectors = min_images(atom_2%id(),atom_1%id())%image_rvectors
      exp_minus_iqr = 0.0_dp
      do k=1,size(rvectors)
        qr = 2*pi*q*rvectors(k)
        exp_minus_iqr = exp_minus_iqr              &
                    & + cmplx(cos(qr),-sin(qr),dp) &
                    & / size(rvectors)
      enddo
      this%matrices(atom_2%prim_id(),atom_1%prim_id()) =      &
         & this%matrices(atom_2%prim_id(),atom_1%prim_id())   &
         & + force_constants%constants(atom_2,atom_1)         &
         & * exp_minus_iqr                                    &
         & / supercell%sc_size
    enddo
  enddo
  
  ! Diagonalise the dynamical matrix, to obtain the normal mode
  !    co-ordinates (eigenvectors) and harmonic frequencies (eigenvalues).
  this%complex_modes = calculate_modes( this%matrices, &
                                      & at_gamma,      &
                                      & paired_qpoint, &
                                      & supercell)
end function

! ----------------------------------------------------------------------
! Construct a dynamical matrix at a q-point commensurate with the supercell,
!    i.e. S.q is a vector of integers.
! ----------------------------------------------------------------------
function new_DynamicalMatrix_calculated(qpoint,at_gamma,paired_qpoint, &
   & supercell,force_constants) result(this)
  use structure_module
  use min_images_module
  use force_constants_module
  use atom_module
  use qpoints_module
  use ofile_module
  implicit none
  
  type(QpointData),     intent(in) :: qpoint
  logical,              intent(in) :: at_gamma
  logical,              intent(in) :: paired_qpoint
  type(StructureData),  intent(in) :: supercell
  type(ForceConstants), intent(in) :: force_constants
  type(DynamicalMatrix)            :: this
  
  type(RealVector)             :: q
  type(IntVector)              :: rvector_i
  type(IntVector)              :: rvector_j
  type(IntVector)              :: rvector_ji
  type(MinImages), allocatable :: min_images(:,:)
  
  integer :: i,j,ialloc
  
  ! Generate a dummy min_images.
  ! Since q is commensurate with the supercell, it doesn't matter which
  !    supercell R-vector, Rs, is chosen, since q.Rs=0.
  allocate( min_images( supercell%no_atoms,  &
          &             supercell%no_atoms), &
          & stat=ialloc); call err(ialloc)
  do i=1,supercell%no_atoms
    rvector_i = supercell%rvectors(supercell%atoms(i)%rvec_id())
    do j=1,supercell%no_atoms
      rvector_j = supercell%rvectors(supercell%atoms(j)%rvec_id())
      rvector_ji = rvector_i - rvector_j
      min_images(j,i) = MinImages([rvector_ji])
    enddo
  enddo
  
  q = dble(qpoint%qpoint)
  this = DynamicalMatrix( q,               &
                        & at_gamma,        &
                        & paired_qpoint,   &
                        & supercell,       &
                        & force_constants, &
                        & min_images)
end function

! ----------------------------------------------------------------------
! Construct a dynamical matrix at a q-point which is not commensurate
!    with the supercell.
! This is only an approximation, using a minimum-image convention.
! ----------------------------------------------------------------------
function new_DynamicalMatrix_interpolated(q,supercell,force_constants, &
   & min_images) result(this)
  use structure_module
  use min_images_module
  use force_constants_module
  use atom_module
  use qpoints_module
  use ofile_module
  implicit none
  
  type(RealVector),     intent(in) :: q
  type(StructureData),  intent(in) :: supercell
  type(ForceConstants), intent(in) :: force_constants
  type(MinImages),      intent(in) :: min_images(:,:)
  type(DynamicalMatrix)            :: this
  
  this = DynamicalMatrix( q,               &
                        & .false.,         &
                        & .false.,         &
                        & supercell,       &
                        & force_constants, &
                        & min_images)
end function

! ----------------------------------------------------------------------
! Diagonalise the dynamical matrix, to obtain the normal mode
!    co-ordinates (eigenvectors) and harmonic frequencies (eigenvalues).
! ----------------------------------------------------------------------
! Structure may be any supercell.
function calculate_modes(matrices,at_gamma,paired_qpoint,structure) &
   & result(output)
  use structure_module
  use min_images_module
  use force_constants_module
  use atom_module
  implicit none
  
  type(ComplexMatrix), intent(in) :: matrices(:,:)
  logical,             intent(in) :: at_gamma
  logical,             intent(in) :: paired_qpoint
  type(StructureData), intent(in) :: structure
  type(ComplexMode), allocatable  :: output(:)
  
  complex(dp), allocatable :: dyn_mat(:,:)
  
  type(ComplexEigenstuff) :: estuff
  
  logical,             allocatable :: translational(:)
  type(ComplexVector)              :: displacement
  
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
  estuff = calculate_eigenstuff(dyn_mat)
  
  ! Identify purely translational modes (at the gamma-point only).
  allocate( translational(structure%no_modes_prim), &
          & stat=ialloc); call err(ialloc)
  translational = .false.
  if (at_gamma) then
    do i=1,3
      j = minloc(abs(estuff%evals),dim=1,mask=.not.translational)
      translational(j) = .true.
    enddo
  endif
  
  ! Calculate normal mode frequencies and displacements.
  !          V = sum_i[ 0.5 * freq_i**2 * u_i**2]
  ! -> F = -2V = sum_i[ - freq_i**2 * u_i**2 ]
  allocate( output(structure%no_modes_prim),   &
          & stat=ialloc); call err(ialloc)
  do i=1,structure%no_modes_prim
    
    ! The eigenvalues are in descending order, but the normal modes should
    !    be in ascending order of frequency. i->k reverses the order.
    k = structure%no_modes_prim - i + 1
    
    if (estuff%evals(k)>=0.0_dp) then
      ! Unstable mode.
      output(i)%frequency = - sqrt(estuff%evals(k))
    else
      ! Stable mode.
      output(i)%frequency = sqrt(- estuff%evals(k))
    endif
    
    output(i)%soft_mode = output(i)%frequency < -1.0e-6_dp
    output(i)%translational_mode = translational(k)
    
    ! Calculate displacements in the primitive cell,
    !    which are the non-mass-reduced eigenvectors of the dynamical matrix.
    if (paired_qpoint) then
      output(i)%at_paired_qpoint = .true.
      allocate( output(i)%primitive_displacements(structure%no_atoms_prim,1), &
              & stat=ialloc); call err(ialloc)
    else
      output(i)%at_paired_qpoint = .false.
      allocate( output(i)%primitive_displacements(structure%no_atoms_prim,2), &
              & stat=ialloc); call err(ialloc)
    endif
    
    do j=1,structure%no_atoms_prim
      displacement = estuff%evecs(3*j-2:3*j, k) &
                 & / sqrt(structure%atoms(j)%mass())
      
      if (paired_qpoint) then
        output(i)%primitive_displacements(j,:) = [displacement]
      else
        output(i)%primitive_displacements(j,:) = [ displacement, &
                                                 & conjg(displacement) ]
      endif
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Construct the dynamical matrix and normal modes at the q-point q'=G-q,
!    given the dynamical matrix and normal modes at the q-point q.
! ----------------------------------------------------------------------
function conjg_DynamicalMatrix(input) result(output)
  use structure_module
  use ofile_module
  implicit none
  
  type(DynamicalMatrix), intent(in) :: input
  type(DynamicalMatrix)             :: output
  
  ! Array sizes.
  integer :: no_atoms
  integer :: no_modes
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  no_atoms = size(input%matrices,1)
  no_modes = size(input%complex_modes)
  
  ! The dynamical matrix at G-q is the complex conjugate of that at q.
  allocate( output%matrices(no_atoms,no_atoms), &
          & stat=ialloc); call err(ialloc)
  do i=1,no_atoms
    do j=1,no_atoms
      output%matrices(j,i) = conjg(input%matrices(j,i))
    enddo
  enddo
  
  ! The frequencies of the normal modes at G-q are the same as those at q.
  allocate(output%complex_modes(no_modes), stat=ialloc); call err(ialloc)
  do i=1,no_modes
    output%complex_modes(i)%frequency = input%complex_modes(i)%frequency
    output%complex_modes(i)%soft_mode = input%complex_modes(i)%soft_mode
    output%complex_modes(i)%translational_mode = &
       & input%complex_modes(i)%translational_mode
  enddo
  
  ! The displacements of the normal modes at G-q are the complex conjugates
  !    of those at q.
  do i=1,no_modes
    allocate( output%complex_modes(i)%primitive_displacements(no_atoms,2), &
            & stat=ialloc); call err(ialloc)
    do j=1,no_atoms
      output%complex_modes(i)%primitive_displacements(j,1) = &
         & input%complex_modes(i)%primitive_displacements(j,2)
      output%complex_modes(i)%primitive_displacements(j,2) = &
         & input%complex_modes(i)%primitive_displacements(j,1)
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Check a dynamical matrix.
! ----------------------------------------------------------------------
! Always checks that the matrix is Hermitian.
! If check_eigenstuff is .true., also checks that the normal modes match
!    the dynamical matrix.
! check_eigenstuff defaults to .true..
! Structure may be any supercell.
subroutine check(this,structure,logfile,check_eigenstuff)
  use utils_module, only : sum_squares
  use structure_module
  use ofile_module
  implicit none
  
  class(DynamicalMatrix), intent(in)           :: this
  type(StructureData),    intent(in)           :: structure
  type(OFile),            intent(inout)        :: logfile
  logical,                intent(in), optional :: check_eigenstuff
  
  type(ComplexMatrix)            :: matrix
  type(ComplexMatrix)            :: hermitian_matrix
  type(ComplexMode), allocatable :: modes(:)
  real(dp)                       :: freq_1
  real(dp)                       :: freq_2
  type(ComplexVector)            :: prim_disp_1
  type(ComplexVector)            :: prim_disp_2
  real(dp)                       :: average
  real(dp)                       :: difference
  logical                        :: check_estuff
  
  integer :: no_atoms
  integer :: i,j
  
  if (present(check_eigenstuff)) then
    check_estuff = check_eigenstuff
  else
    check_estuff = .true.
  endif
  
  no_atoms = size(this%matrices,1)
  
  ! Check the dynamical matrix is Hermitian.
  average = 0.0_dp
  difference = 0.0_dp
  do i=1,no_atoms
    do j=1,i
      matrix = this%matrices(j,i)
      hermitian_matrix = hermitian(this%matrices(i,j))
      
      average = average + sum_squares((matrix+hermitian_matrix)/2.0_dp)
      difference = difference + sum_squares(matrix-hermitian_matrix)
    enddo
  enddo
  call logfile%print_line(                                        &
     & 'Fractional L2 error in hermicity of dynamical matrix :'// &
     & sqrt(difference/average))
  if (sqrt(difference/average)>1.0e-10_dp) then
    call print_line(WARNING//': Dynamical matrix is not hermitian. Please &
       &check log files.')
  endif
  
  ! Check that dynamical matrix and normal modes match.
  if (check_estuff) then
    modes = calculate_modes( this%matrices,                          &
                           & .false.,                                &
                           & this%complex_modes(1)%at_paired_qpoint, &
                           & structure)
    
    ! Check that eigenfrequencies match.
    average = 0.0_dp
    difference = 0.0_dp
    do i=1,structure%no_modes_prim
      if (this%complex_modes(i)%translational_mode) then
        cycle
      endif
      freq_1 = this%complex_modes(i)%frequency
      freq_2 = modes(i)%frequency
      average = average + ((freq_1+freq_2)/2)**2
      difference = difference + (freq_1-freq_2)**2
    enddo
    call logfile%print_line(                           &
       & 'Fractional L2 error in eigenfrequencies: '// &
       & sqrt(difference/average))
    if (sqrt(difference/average) > 1e-10_dp) then
      call print_line(WARNING//': Eigenfrequencies do not match. Please &
         &check log files.')
    endif
    
    ! Check that the primitive displacements match.
    ! N.B. the global phase of the displacements is not necessarily consistent.
    average = 0.0_dp
    difference = 0.0_dp
    do i=1,structure%no_modes_prim
      ! Ignore degenerate modes.
      if (this%complex_modes(i)%translational_mode) then
        cycle
      endif
      if (i>1) then
        freq_1 = modes(i)%frequency
        freq_2 = modes(i-1)%frequency
        if (abs((freq_1-freq_2)/freq_1)<0.01) then
          cycle
        endif
      endif
      if (i<structure%no_modes_prim) then
        freq_1 = modes(i)%frequency
        freq_2 = modes(i+1)%frequency
        if (abs((freq_1-freq_2)/freq_1)<0.01) then
          cycle
        endif
      endif
      
      do j=1,structure%no_atoms_prim
        prim_disp_1 = this%complex_modes(i)%primitive_displacements(j,1)
        prim_disp_2 = modes(i)%primitive_displacements(j,1)
        ! Ignore phases.
        prim_disp_1 = cmplx(vec(abs(cmplx(prim_disp_1))))
        prim_disp_2 = cmplx(vec(abs(cmplx(prim_disp_2))))
        average = average + sum_squares((prim_disp_1+prim_disp_2)/2)
        difference = difference + sum_squares(prim_disp_1-prim_disp_2)
      enddo
    enddo
    call logfile%print_line(                                          &
       & 'Fractional L2 error in rotated primitive displacements: '// &
       & sqrt(difference/average))
    if (sqrt(difference/average) > 1e-10_dp) then
      call print_line(WARNING//': Error in primitive displacements. &
         &Please check log files.')
      call print_line(difference//' / '//average)
    endif
  endif
end subroutine

! ----------------------------------------------------------------------
! Construct the matrix of force constants for the large supercell, given the
!    dynamical matrices at each q-point.
! ----------------------------------------------------------------------
function reconstruct_force_constants(large_supercell,qpoints, &
   & dynamical_matrices,logfile) result(output)
  use constants_module, only : pi
  use structure_module
  use qpoints_module
  use force_constants_module
  use atom_module
  use ofile_module
  implicit none
  
  type(StructureData),   intent(in)    :: large_supercell
  type(QpointData),      intent(in)    :: qpoints(:)
  type(DynamicalMatrix), intent(in)    :: dynamical_matrices(:)
  type(OFile),           intent(inout) :: logfile
  type(ForceConstants)                 :: output
  
  type(RealMatrix), allocatable :: force_constants(:,:)
  
  type(AtomData) :: atom_1
  type(AtomData) :: atom_2
  
  type(IntVector)  :: r
  type(RealVector) :: q
  real(dp)         :: qr
  complex(dp)      :: exp_iqr
  
  integer :: i,j,k,ialloc
  
  allocate( force_constants( large_supercell%no_atoms,  &
          &                  large_supercell%no_atoms), &
          & stat=ialloc); call err(ialloc)
  force_constants = mat(dble(zeroes(3,3)))
  
  ! Loop across q-points, summing up the contribution from the dynamical matrix
  !    at each.
  do i=1,large_supercell%sc_size
    q = dble(qpoints(i)%qpoint)
    do j=1,large_supercell%no_atoms
      atom_1 = large_supercell%atoms(j)
      do k=1,large_supercell%no_atoms
        atom_2 = large_supercell%atoms(k)
        
        ! Calculate exp(2*pi*i * q.(R2-R1))
        r = large_supercell%rvectors(atom_2%rvec_id()) &
        & - large_supercell%rvectors(atom_1%rvec_id())
        qr = 2 * pi * q * r
        exp_iqr = cmplx(cos(qr),sin(qr),dp)
        
        ! Add in the contribution to the force constant matrix.
        force_constants(atom_1%id(),atom_2%id()) =                     &
           &   force_constants(atom_1%id(),atom_2%id())                &
           & + real( dynamical_matrices(i)%matrices( atom_1%prim_id(), &
           &                                         atom_2%prim_id()) &
           &       * exp_iqr)
      enddo
    enddo
  enddo
  
  output = ForceConstants(large_supercell, force_constants, logfile)
end function

! ----------------------------------------------------------------------
! Writes a dynamical matrix to file.
! ----------------------------------------------------------------------
subroutine write_dynamical_matrix_file(dynamical_matrix,filename)
  use ofile_module
  implicit none
  
  type(DynamicalMatrix), intent(in) :: dynamical_matrix
  type(String),          intent(in) :: filename
  
  type(OFile) :: matrix_file
  
  integer :: no_atoms
  
  integer :: i,j
  
  no_atoms = size(dynamical_matrix%matrices,1)
  if (size(dynamical_matrix%matrices,2)/=no_atoms) then
    call err()
  endif
  
  matrix_file = filename
  do i=1,no_atoms
    do j=1,no_atoms
      call matrix_file%print_line('Atoms: '//i//' and '//j//'.')
      call matrix_file%print_line(dynamical_matrix%matrices(j,i))
      call matrix_file%print_line('')
    enddo
  enddo
end subroutine

! ----------------------------------------------------------------------
! Reads a dynamical matrix from file.
! ----------------------------------------------------------------------
function read_dynamical_matrix_file(filename) result(dynamical_matrix)
  use utils_module, only : int_sqrt
  use ifile_module
  implicit none
  
  type(String), intent(in) :: filename
  type(DynamicalMatrix)    :: dynamical_matrix
  
  type(IFile) :: matrix_file
  
  integer :: no_atoms
  
  ! Temporary variables
  integer                   :: i,j,k,ialloc
  type(String), allocatable :: line(:)
  complex(dp)               :: matrix(3,3)
  
  matrix_file = filename
  
  no_atoms = int_sqrt(size(matrix_file)/5)
  
  allocate( dynamical_matrix%matrices(no_atoms,no_atoms), &
          & stat=ialloc); call err(ialloc)
  
  do i=1,no_atoms
    do j=1,no_atoms
      do k=1,3
        line = split(matrix_file%line(5*(no_atoms*(i-1)+(j-1))+1+k))
        matrix(k,:) = cmplx(line)
      enddo
      dynamical_matrix%matrices(j,i) = matrix
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Rotate a dynamical matrix and set of normal modes onto a new q-point.
! ----------------------------------------------------------------------
! Construct data at q_new from data at q_old, where
!    R . q_old = q_new.
! N.B. the provided q-point should be q_new not q_old.
function rotate_modes(input,symmetry,qpoint) result(output)
  use constants_module, only : pi
  use utils_module, only : sum_squares
  use qpoints_module
  use symmetry_module
  implicit none
  
  type(DynamicalMatrix),  intent(in)    :: input
  type(SymmetryOperator), intent(in)    :: symmetry
  type(QpointData),       intent(in)    :: qpoint
  type(DynamicalMatrix)                 :: output
  
  ! q.Rj an exp(i q.Rj), where Rj is symmetry%rvector(j)
  type(RealVector)         :: q
  real(dp)                 :: qr
  complex(dp), allocatable :: exp_iqr(:)
  
  ! Array sizes.
  integer :: no_atoms
  integer :: no_modes
  
  ! Atom labels.
  integer :: atom_1
  integer :: atom_1p
  integer :: atom_2
  integer :: atom_2p
  integer :: mode
  
  ! Temporary variables.
  integer :: i,ialloc
  
  no_atoms = size(input%matrices,1)
  no_modes = 3*no_atoms
  
  ! Construct phases.
  ! The symmetry, S, maps equilibrium position ri to rj+R,
  !    and q-point q to q'.
  ! The phase change is then (q').(rj+R)-q.ri.
  
  ! S : r -> r' + R
  ! S : q -> q'
  ! q'.r'-q.r = q'.R
  
  allocate(exp_iqr(no_atoms), stat=ialloc); call err(ialloc)
  do i=1,no_atoms
    q = dble(qpoint%qpoint)
    qr = 2*pi*q*symmetry%rvector(i)
    exp_iqr(i) = cmplx(cos(qr),sin(qr),dp)
  enddo
  
  ! Rotate dynamical matrix.
  allocate(output%matrices(no_atoms,no_atoms), stat=ialloc); call err(ialloc)
  do atom_1=1,no_atoms
    atom_1p = symmetry%atom_group * atom_1
    do atom_2=1,no_atoms
      atom_2p = symmetry%atom_group * atom_2
      
      output%matrices(atom_2p,atom_1p) = symmetry%cartesian_rotation   &
                                     & * exp_iqr(atom_2)               &
                                     & * input%matrices(atom_2,atom_1) &
                                     & * conjg(exp_iqr(atom_1))        &
                                     & * transpose(symmetry%cartesian_rotation)
    enddo
  enddo
  
  ! Rotate normal modes.
  if ( input%complex_modes(1)%at_paired_qpoint .neqv. &
     & qpoint%is_paired_qpoint) then
    call print_line(CODE_ERROR//': Error matching rotated q-points.')
    call err()
  endif
  
  allocate(output%complex_modes(no_atoms), stat=ialloc); call err(ialloc)
  output%complex_modes = input%complex_modes
  do mode=1,no_modes
    do atom_1=1,no_atoms
      atom_1p = symmetry%atom_group * atom_1
      if (qpoint%is_paired_qpoint) then
        output%complex_modes(mode)%primitive_displacements(atom_1p,1) =    &
           & symmetry%cartesian_rotation                                   &
           & * input%complex_modes(mode)%primitive_displacements(atom_1,1) &
           & * exp_iqr(atom_1)
      else
        output%complex_modes(mode)%primitive_displacements(atom_1p,1) =    &
           & symmetry%cartesian_rotation                                   &
           & * input%complex_modes(mode)%primitive_displacements(atom_1,1) &
           & * exp_iqr(atom_1)
        output%complex_modes(mode)%primitive_displacements(atom_1p,2) =    &
           & symmetry%cartesian_rotation                                   &
           & * input%complex_modes(mode)%primitive_displacements(atom_1,2) &
           & * exp_iqr(atom_1)
      endif
    enddo
  enddo
end function

subroutine print_dyn_mat(input)
  implicit none
  
  type(DynamicalMatrix), intent(in) :: input
  
  integer :: no_atoms
  integer :: i,j,k,l
  
  complex(dp) :: thingy(3,3)
  character(50) :: thing
  
  type(String) :: line
  
  no_atoms = size(input%matrices,1)
  
  do i=1,no_atoms
    do j=1,3
      line = ''
      do k=1,no_atoms
        do l=1,3
          thingy = cmplx(input%matrices(i,k))
          write(thing,fmt='(es8.0,sp,es7.0,"i")') thingy(j,l)
          line = line//trim(thing)
        enddo
      enddo
      call print_line(line)
    enddo
  enddo
end subroutine

! ----------------------------------------------------------------------
! Returns mode frequencies as a single array.
! ----------------------------------------------------------------------
function frequencies(this) result(output)
  implicit none
  
  class(DynamicalMatrix), intent(in) :: this
  real(dp), allocatable              :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(this%complex_modes)), stat=ialloc); call err(ialloc)
  do i=1,size(this%complex_modes)
    output(i) = this%complex_modes(i)%frequency
  enddo
end function

! ----------------------------------------------------------------------
! Compares two dynamical matrices.
! ----------------------------------------------------------------------
subroutine compare_dynamical_matrices(a,b,logfile)
  use utils_module, only : sum_squares
  use ofile_module
  implicit none
  
  type(DynamicalMatrix), intent(in)    :: a
  type(DynamicalMatrix), intent(in)    :: b
  type(OFile),           intent(inout) :: logfile
  
  integer :: no_atoms
  
  type(ComplexMatrix) :: mat_a
  type(ComplexMatrix) :: mat_b
  
  real(dp) :: average
  real(dp) :: difference
  
  integer :: i,j
  
  no_atoms = size(a%matrices,1)
  if (size(b%matrices,1)/=no_atoms) then
    call print_line(CODE_ERROR//': dynamical matrices a and b have different &
       &sizes.')
    call err()
  endif
  
  average = 0.0_dp
  difference = 0.0_dp
  do i=1,no_atoms
    do j=1,no_atoms
      mat_a = a%matrices(j,i)
      mat_b = b%matrices(j,i)
      average = average + sum_squares((mat_a+mat_b)/2)
      difference = difference + sum_squares(mat_a-mat_b)
    enddo
  enddo
  call logfile%print_line('Fractional L2 difference between dynamical &
     &matrices: '//sqrt(difference/average))
  if (sqrt(difference/average)>1e-10_dp) then
    call print_line(WARNING//': Dynamical matrices differ. Please check &
       &log files.')
  endif
end subroutine
end module
