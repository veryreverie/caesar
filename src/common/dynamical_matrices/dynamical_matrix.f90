! ======================================================================
! The forces between atoms at a given q-point.
! ======================================================================
module dynamical_matrix_module
  use utils_module
  
  use structure_module
  use normal_mode_module
  
  use min_images_module
  use cartesian_hessian_module
  implicit none
  
  private
  
  public :: DynamicalMatrix
  public :: reconstruct_hessian
  public :: ComplexMode
  
  type, extends(Stringsable) :: DynamicalMatrix
    type(ComplexMatrix), allocatable, private :: matrices_(:,:)
  contains
    procedure, public :: check
    
    procedure, public :: matrices => matrices_DynamicalMatrix
    
    ! I/O.
    procedure, public :: read  => read_DynamicalMatrix
    procedure, public :: write => write_DynamicalMatrix
  end type
  
  interface DynamicalMatrix
    module procedure new_DynamicalMatrix
    module procedure new_DynamicalMatrix_interpolated
    module procedure new_DynamicalMatrix_ComplexModes
    module procedure new_DynamicalMatrix_Strings
    module procedure new_DynamicalMatrix_StringArray
  end interface
  
  interface ComplexMode
    module procedure new_ComplexMode_DynamicalMatrix
  end interface
  
  interface calculate_modes
    module procedure calculate_modes_complex
    module procedure calculate_modes_real
  end interface
contains

! Constructors.
function new_DynamicalMatrix(matrices) result(this)
  implicit none
  
  type(ComplexMatrix), intent(in) :: matrices(:,:)
  type(DynamicalMatrix)           :: this
  
  this%matrices_ = matrices
end function

! Getter for matrices.
function matrices_DynamicalMatrix(this) result(output)
  implicit none
  
  class(DynamicalMatrix), intent(in) :: this
  type(ComplexMatrix), allocatable   :: output(:,:)
  
  output = this%matrices_
end function

! ----------------------------------------------------------------------
! Constructs the dynamical matrix, which is the matrix of force constants in
!    q-point co-ordinates,
!    given the Hessian, which is the matrix of force constants in
!    cartesian supercell co-ordinates.
! ----------------------------------------------------------------------

! --------------------------------------------------
! Construct a dynamical matrix at a q-point which is not commensurate
!    with the supercell.
! This is only an approximation, using a minimum-image convention.
! --------------------------------------------------
function new_DynamicalMatrix_interpolated(q,supercell,hessian,min_images) &
   & result(this)
  implicit none
  
  type(RealVector),       intent(in)           :: q
  type(StructureData),    intent(in)           :: supercell
  type(CartesianHessian), intent(in)           :: hessian
  type(MinImages),        intent(in), optional :: min_images(:,:)
  type(DynamicalMatrix)                        :: this
  
  ! Evaluate the dynamical matrix.
  this%matrices_ = calculate_dynamical_matrix( q,         &
                                             & supercell, &
                                             & hessian,   &
                                             & min_images )
end function

function calculate_dynamical_matrix(q,supercell,hessian,min_images) &
   & result(output)
  implicit none
  
  type(RealVector),       intent(in)           :: q
  type(StructureData),    intent(in)           :: supercell
  type(CartesianHessian), intent(in)           :: hessian
  type(MinImages),        intent(in), optional :: min_images(:,:)
  type(ComplexMatrix), allocatable             :: output(:,:)
  
  type(AtomData)               :: atom_1
  type(AtomData)               :: atom_2
  type(IntVector)              :: rvector
  type(IntVector), allocatable :: rvectors(:)
  
  integer :: i,j,k,ialloc
  
  ! Allocate and zero output.
  allocate( output( supercell%no_atoms_prim,  &
          &         supercell%no_atoms_prim), &
          & stat=ialloc); call err(ialloc)
  output = cmplxmat(zeroes(3,3))
  
  ! Add up contributions to the dynamical matrix from
  !    the force constants between each pair of atoms.
  do i=1,supercell%no_atoms
    atom_1 = supercell%atoms(i)
    
    ! Atom 2 will always be at R=0, so the R-vector from atom 2 to atom 1
    !    is simply that of atom 1.
    rvector = supercell%rvectors(atom_1%rvec_id())
    do j=1,supercell%no_atoms_prim
      atom_2 = supercell%atoms(j)
      
      if (present(min_images)) then
        rvectors = min_images(atom_2%id(),atom_1%id())%image_rvectors
        ! Check that all min image R-vectors actually point from atom 2 to
        !    a copy of atom 1.
        do k=1,size(rvectors)
          if (.not. is_int( supercell%recip_supercell &
                        & * (rvectors(k)-rvector))) then
            call print_line(CODE_ERROR// &
               & ': Problem with minimum image R-vectors.')
            call err()
          endif
        enddo
      else
        rvectors = [rvector]
      endif
      
      output(atom_2%prim_id(),atom_1%prim_id()) =    &
         & output(atom_2%prim_id(),atom_1%prim_id()) &
         & + hessian%elements(atom_2,atom_1)         &
         & * sum(exp_2pii(q*rvectors))               &
         & / (supercell%sc_size*size(rvectors))
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
subroutine check(this,structure,logfile)
  implicit none
  
  class(DynamicalMatrix), intent(in)    :: this
  type(StructureData),    intent(in)    :: structure
  type(OFile),            intent(inout) :: logfile
  
  type(ComplexMatrix)            :: matrix
  type(ComplexMatrix)            :: hermitian_matrix
  type(ComplexMode), allocatable :: modes(:)
  real(dp)                       :: freq_1
  real(dp)                       :: freq_2
  type(ComplexVector)            :: prim_vec_1
  type(ComplexVector)            :: prim_vec_2
  real(dp)                       :: average
  real(dp)                       :: difference
  
  integer :: no_atoms
  integer :: i,j
  
  no_atoms = size(this%matrices_,1)
  
  ! Check the dynamical matrix is Hermitian.
  average = 0.0_dp
  difference = 0.0_dp
  do i=1,no_atoms
    do j=1,i
      matrix = this%matrices_(j,i)
      hermitian_matrix = hermitian(this%matrices_(i,j))
      
      average = average + sum_squares((matrix+hermitian_matrix)/2.0_dp)
      difference = difference + sum_squares(matrix-hermitian_matrix)
    enddo
  enddo
  call logfile%print_line(                                        &
     & 'Fractional L2 error in hermicity of dynamical matrix :'// &
     & sqrt(difference/average))
  if (average>1e-20_dp) then
    if (sqrt(difference/average)>1.0e-10_dp) then
      call print_line(WARNING//': Dynamical matrix is not Hermitian. Please &
         &check log files.')
      call print_line('Fractional L2 error in hermicity: '// &
         & sqrt(difference/average))
    endif
  endif
end subroutine

! ----------------------------------------------------------------------
! Calculates complex modes by diagonalising a dynamical matrix.
! ----------------------------------------------------------------------
! N.B. Structure may be any supercell.

function new_ComplexMode_DynamicalMatrix(dynamical_matrix,structure, &
   & modes_real) result(output)
  implicit none
  
  type(DynamicalMatrix), intent(in)           :: dynamical_matrix
  type(StructureData),   intent(in)           :: structure
  logical,               intent(in), optional :: modes_real
  type(ComplexMode), allocatable              :: output(:)
  
  if (present(modes_real)) then
    if (modes_real) then
      output = calculate_modes(real(dynamical_matrix%matrices()), structure)
      return
    endif
  endif
  
  output = calculate_modes(dynamical_matrix%matrices(), structure)
end function

! Calculate modes at a q-point where 2q/=G.
function calculate_modes_complex(matrices,structure) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: matrices(:,:)
  type(StructureData), intent(in) :: structure
  type(ComplexMode), allocatable  :: output(:)
  
  complex(dp),               allocatable :: dyn_mat(:,:)
  type(HermitianEigenstuff), allocatable :: estuff(:)
  
  integer :: i,j,ialloc
  
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
  
  ! Calculate normal modes.
  ! Eigenvalues are reversed so that the complex modes are in ascending order.
  output = ComplexMode(estuff(size(estuff):1:-1), structure)
end function

! Calculate modes at a q-point where 2q=G.
function calculate_modes_real(matrices,structure) result(output)
  implicit none
  
  type(RealMatrix),    intent(in) :: matrices(:,:)
  type(StructureData), intent(in) :: structure
  type(ComplexMode), allocatable  :: output(:)
  
  real(dp),                  allocatable :: dyn_mat(:,:)
  type(SymmetricEigenstuff), allocatable :: real_estuff(:)
  type(HermitianEigenstuff), allocatable :: estuff(:)
  
  integer :: i,j,ialloc
  
  ! Convert (3x3Matrix) x no_atoms x no_atoms to no_modes x no_modes
  allocate( dyn_mat(structure%no_modes_prim,structure%no_modes_prim), &
          & stat=ialloc); call err(ialloc)
  do i=1,structure%no_atoms_prim
    do j=1,structure%no_atoms_prim
      dyn_mat(3*j-2:3*j, 3*i-2:3*i) = dble(matrices(j,i))
    enddo
  enddo
  
  ! Diagonalise dynamical matrix.
  real_estuff = diagonalise_symmetric(dyn_mat)
  estuff = [(                                                         &
     & HermitianEigenstuff( real_estuff(i)%eval,                      &
     &                      [cmplx(real_estuff(i)%evec,0.0_dp,dp)] ), &
     & i=1,                                                           &
     & size(real_estuff)                                              )]
  
  ! Calculate normal modes.
  ! Eigenvalues are reversed so that the complex modes are in ascending order.
  output = ComplexMode(estuff(size(estuff):1:-1), structure)
end function

! ----------------------------------------------------------------------
! Construct a DynamicalMatrix from the ComplexModes at a given q-point.
! ----------------------------------------------------------------------
function new_DynamicalMatrix_ComplexModes(modes,frequencies) result(this)
  implicit none
  
  type(ComplexMode), intent(in)           :: modes(:)
  real(dp),          intent(in), optional :: frequencies(:)
  type(DynamicalMatrix)                   :: this
  
  type(ComplexMode), allocatable :: new_modes(:)
  
  integer :: no_atoms
  
  integer :: i,j,k,ialloc
  
  if (present(frequencies)) then
    if (size(modes)/=size(frequencies)) then
      call print_line(CODE_ERROR//': modes and frequencies do not match.')
      call err()
    endif
  endif
  
  if (size(modes)==0) then
    no_atoms = 1
  else
    no_atoms = size(modes(1)%unit_vector)
  endif
  
  if (size(modes)/=3*no_atoms .and. size(modes)/=3*(no_atoms-1)) then
    call print_line(CODE_ERROR//': modes and no_atoms do not match.')
  endif
  
  do i=1,size(modes)
    if (size(modes(i)%unit_vector)/=no_atoms) then
      call print_line(CODE_ERROR//': inconsistent no_atoms.')
      call err()
    endif
  enddo
  
  new_modes = modes
  if (present(frequencies)) then
    do i=1,size(frequencies)
      new_modes(i)%frequency = frequencies(i)
      if (frequencies(i)>=0) then
        new_modes(i)%spring_constant = frequencies(i)**2
      else
        new_modes(i)%spring_constant = -frequencies(i)**2
      endif
    enddo
  endif
  
  ! A dynamical matrix is D = sum_i e_i u_i^(u_i)*.
  ! The eigenvalue e_i is the negative of the mode's spring constant.
  allocate(this%matrices_(no_atoms,no_atoms), stat=ialloc); call err(ialloc)
  this%matrices_ = cmplxmat(zeroes(3,3))
  do i=1,size(modes)
    do j=1,no_atoms
      do k=1,no_atoms
        this%matrices_(k,j) =                                    &
           &   this%matrices_(k,j)                               &
           & - new_modes(i)%spring_constant                      &
           & * outer_product( new_modes(i)%unit_vector(k),       &
           &                  conjg(new_modes(i)%unit_vector(j)) )
      enddo
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Construct the Hessian for the large supercell, given the
!    dynamical matrices at each q-point.
! ----------------------------------------------------------------------
function reconstruct_hessian(large_supercell,qpoints,dynamical_matrices, &
   & logfile) result(output)
  implicit none
  
  type(StructureData),   intent(in)    :: large_supercell
  type(QpointData),      intent(in)    :: qpoints(:)
  type(DynamicalMatrix), intent(in)    :: dynamical_matrices(:)
  type(OFile),           intent(inout) :: logfile
  type(CartesianHessian)               :: output
  
  type(RealMatrix), allocatable :: hessian(:,:)
  
  type(AtomData) :: atom_1
  type(AtomData) :: atom_2
  
  type(IntVector)  :: r
  type(RealVector) :: q
  
  integer :: i,j,k,ialloc
  
  allocate( hessian( large_supercell%no_atoms,  &
          &          large_supercell%no_atoms), &
          & stat=ialloc); call err(ialloc)
  hessian = dblemat(zeroes(3,3))
  
  ! Loop across q-points, summing up the contribution from
  !    the dynamical matrix at each.
  do i=1,large_supercell%sc_size
    q = dblevec(qpoints(i)%qpoint)
    do j=1,large_supercell%no_atoms
      atom_1 = large_supercell%atoms(j)
      do k=1,large_supercell%no_atoms
        atom_2 = large_supercell%atoms(k)
        
        r = large_supercell%rvectors(atom_2%rvec_id()) &
        & - large_supercell%rvectors(atom_1%rvec_id())
        
        ! Add in the contribution to the Hessian.
        hessian(atom_1%id(),atom_2%id()) =                               &
           &   hessian(atom_1%id(),atom_2%id())                          &
           & + real( dynamical_matrices(i)%matrices_( atom_1%prim_id(),  &
           &                                          atom_2%prim_id() ) &
           &       * exp_2pii(-q*r)                                      )
      enddo
    enddo
  enddo
  
  output = CartesianHessian(large_supercell, hessian, logfile)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_DynamicalMatrix(this,input)
  implicit none
  
  class(DynamicalMatrix), intent(out) :: this
  type(String),           intent(in)  :: input(:)
  
  type(StringArray), allocatable :: elements(:)
  integer                        :: no_atoms
  
  integer :: i,j,k,ialloc
  
  select type(this); type is(DynamicalMatrix)
    elements = split_into_sections(input)
    no_atoms = int_sqrt(size(elements))
    allocate(this%matrices_(no_atoms,no_atoms), stat=ialloc); call err(ialloc)
    k = 0
    do i=1,no_atoms
      do j=1,no_atoms
        k = k+1
        this%matrices_(j,i) = ComplexMatrix(elements(k)%strings(2:4))
      enddo
    enddo
  class default
    call err()
  end select
end subroutine

function write_DynamicalMatrix(this) result(output)
  implicit none
  
  class(DynamicalMatrix), intent(in) :: this
  type(String), allocatable          :: output(:)
  
  integer :: no_atoms
  
  type(String) :: matrix_strings(3)
  
  integer :: i,j,k,ialloc
  
  select type(this); type is(DynamicalMatrix)
    no_atoms = size(this%matrices_,1)
    if (size(this%matrices_,2)/=no_atoms) then
      call err()
    endif
    
    allocate(output(5*no_atoms*no_atoms), stat=ialloc); call err(ialloc)
    k = 0
    do i=1,no_atoms
      do j=1,no_atoms
        k = k+1
        matrix_strings = str(this%matrices_(j,i))
        output(5*k-4) = 'Atoms: ('//j//' '//i//')'
        output(5*k-3) = matrix_strings(1)
        output(5*k-2) = matrix_strings(2)
        output(5*k-1) = matrix_strings(3)
        output(5*k)   = ''
      enddo
    enddo
  class default
    call err()
  end select
end function

function new_DynamicalMatrix_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(DynamicalMatrix)    :: this
  
  call this%read(input)
end function

impure elemental function new_DynamicalMatrix_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(DynamicalMatrix)         :: this
  
  this = DynamicalMatrix(str(input))
end function
end module
