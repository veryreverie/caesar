! ======================================================================
! The forces between atoms at a given q-point.
! ======================================================================
module dynamical_matrix_module
  use utils_module
  
  use structure_module
  use normal_mode_module
  
  use min_images_module
  implicit none
  
  private
  
  public :: DynamicalMatrix
  public :: reconstruct_hessian
  public :: ComplexMode
  public :: conjg
  
  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: operator(/)
  
  type, extends(Stringsable) :: DynamicalMatrix
    type(ComplexMatrix), allocatable, private :: elements_(:,:)
  contains
    procedure, public :: check
    
    procedure, public :: elements => elements_DynamicalMatrix
    
    procedure, public :: expectation => expectation_DynamicalMatrix
    
    ! I/O.
    procedure, public :: read  => read_DynamicalMatrix
    procedure, public :: write => write_DynamicalMatrix
  end type
  
  interface DynamicalMatrix
    module procedure new_DynamicalMatrix
    module procedure new_DynamicalMatrix_zeroes
    module procedure new_DynamicalMatrix_interpolated
    module procedure new_DynamicalMatrix_ComplexModes
    module procedure new_DynamicalMatrix_ComplexMode_ComplexMode
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
  
  interface operator(+)
    module procedure add_DynamicalMatrix_DynamicalMatrix
  end interface
  
  interface operator(-)
    module procedure negative_DynamicalMatrix
    module procedure subtract_DynamicalMatrix_DynamicalMatrix
  end interface
  
  interface operator(*)
    module procedure multiply_DynamicalMatrix_real
    module procedure multiply_real_DynamicalMatrix
    module procedure multiply_DynamicalMatrix_complex
    module procedure multiply_complex_DynamicalMatrix
  end interface
  
  interface operator(/)
    module procedure divide_DynamicalMatrix_real
    module procedure divide_DynamicalMatrix_complex
  end interface
  
  interface conjg
    module procedure conjg_DynamicalMatrix
  end interface
contains

! Constructors.
function new_DynamicalMatrix(elements) result(this)
  implicit none
  
  type(ComplexMatrix), intent(in) :: elements(:,:)
  type(DynamicalMatrix)           :: this
  
  this%elements_ = elements
end function

function new_DynamicalMatrix_zeroes(no_atoms) result(this)
  implicit none
  
  integer, intent(in)   :: no_atoms
  type(DynamicalMatrix) :: this
  
  integer :: ialloc
  
  allocate(this%elements_(no_atoms,no_atoms), stat=ialloc); call err(ialloc)
  this%elements_ = cmplxmat(zeroes(3,3))
end function

! Getter for elements.
function elements_DynamicalMatrix(this) result(output)
  implicit none
  
  class(DynamicalMatrix), intent(in) :: this
  type(ComplexMatrix), allocatable   :: output(:,:)
  
  output = this%elements_
end function

! The expectation of the dynamical matrix w/r/t a given mode.
impure elemental function expectation_DynamicalMatrix(this,mode) result(output)
  implicit none
  
  class(DynamicalMatrix), intent(in) :: this
  type(ComplexMode),      intent(in) :: mode
  real(dp)                           :: output
  
  integer :: i,j
  
  output = 0
  do i=1,size(mode%unit_vector)
    do j=1,size(mode%unit_vector)
      output = output + conjg(mode%unit_vector(i)) &
                    & * this%elements_(i,j)        &
                    & * mode%unit_vector(j)
    enddo
  enddo
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
  
  type(ComplexMatrix), allocatable :: elements(:,:)
  
  type(AtomData)               :: atom_1
  type(AtomData)               :: atom_2
  type(IntVector)              :: rvector
  type(IntVector), allocatable :: rvectors(:)
  
  integer :: i,j,k,ialloc
  
  ! Allocate and zero elements.
  allocate( elements( supercell%no_atoms_prim,  &
          &           supercell%no_atoms_prim), &
          & stat=ialloc); call err(ialloc)
  elements = cmplxmat(zeroes(3,3))
  
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
      
      elements(atom_2%prim_id(),atom_1%prim_id()) =    &
         & elements(atom_2%prim_id(),atom_1%prim_id()) &
         & + hessian%elements(atom_2,atom_1)           &
         & * sum(exp_2pii(q*rvectors))                 &
         & / (supercell%sc_size*size(rvectors))
    enddo
  enddo
  
  ! Construct output.
  this = DynamicalMatrix(elements)
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
  
  no_atoms = size(this%elements_,1)
  
  ! Check the dynamical matrix is Hermitian.
  average = 0.0_dp
  difference = 0.0_dp
  do i=1,no_atoms
    do j=1,i
      matrix = this%elements_(j,i)
      hermitian_matrix = hermitian(this%elements_(i,j))
      
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
      output = calculate_modes(real(dynamical_matrix%elements()), structure)
      return
    endif
  endif
  
  output = calculate_modes(dynamical_matrix%elements(), structure)
end function

! Calculate modes at a q-point where 2q/=G.
function calculate_modes_complex(elements,structure) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in) :: elements(:,:)
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
      dyn_mat(3*j-2:3*j, 3*i-2:3*i) = cmplx(elements(j,i))
    enddo
  enddo
  
  ! Diagonalise dynamical matrix.
  estuff = diagonalise_hermitian(dyn_mat)
  
  ! Calculate normal modes.
  ! Eigenvalues are reversed so that the complex modes are in ascending order.
  output = ComplexMode(estuff(size(estuff):1:-1), structure)
end function

! Calculate modes at a q-point where 2q=G.
function calculate_modes_real(elements,structure) result(output)
  implicit none
  
  type(RealMatrix),    intent(in) :: elements(:,:)
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
      dyn_mat(3*j-2:3*j, 3*i-2:3*i) = dble(elements(j,i))
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
  
  integer :: i
  
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
  
  ! Construct the dynamical matrix as the sum of the outer products of each
  !    mode times that mode's eigenvalue, e_i.
  ! The eigenvalue e_i is the negative of the mode's spring constant.
  this = DynamicalMatrix(no_atoms)
  do i=1,size(modes)
    this = this                                           &
       & + DynamicalMatrix( new_modes(i),                 &
       &                    new_modes(i),                 &
       &                    -new_modes(i)%spring_constant )
  enddo
end function

! ----------------------------------------------------------------------
! Construct the contribution to a DynamicalMatrix two ComplexModes.
! If they are at the same q-point, the conjugate of the second mode is taken.
! If they are at paired q-points, no conjugate is taken.
! If they are neither at the same q-point nor at paired q-points,
!    an error is thrown.
! ----------------------------------------------------------------------
function new_DynamicalMatrix_ComplexMode_ComplexMode(mode1,mode2,coefficient) &
   & result(this)
  implicit none
  
  type(ComplexMode), intent(in)           :: mode1
  type(ComplexMode), intent(in)           :: mode2
  real(dp),          intent(in), optional :: coefficient
  type(DynamicalMatrix)                   :: this
  
  real(dp) :: coefficient_
  
  integer :: no_atoms
  
  integer :: i,j
  
  if (present(coefficient)) then
    coefficient_ = coefficient
  else
    coefficient_ = 1
  endif
  
  no_atoms = size(mode1%unit_vector)
  if (size(mode2%unit_vector)/=no_atoms) then
    call print_line(CODE_ERROR//': Trying to construct a dynamical matrix &
       &from two modes with different numbers of atoms.')
    call err()
  endif
  
  this = DynamicalMatrix(no_atoms)
  if (mode1%qpoint_id==mode2%qpoint_id) then
    ! If q_1 = q_2 then the ij element of D is the outer product of the i
    !    component of mode1 with the j component of (mode2)*.
    do i=1,no_atoms
      do j=1,no_atoms
        this%elements_(j,i) = this%elements_(j,i)                        &
                          & + coefficient                                &
                          & * outer_product( mode1%unit_vector(j),       &
                          &                  conjg(mode2%unit_vector(i)) )
      enddo
    enddo
  elseif (mode1%qpoint_id==mode2%paired_qpoint_id) then
    ! If q_1 = -q_2 then the ij element of D is the outer product of the i
    !    component of mode1 with the j component of mode2.
    do i=1,no_atoms
      do j=1,no_atoms
        this%elements_(j,i) = this%elements_(j,i)                  &
                          & + coefficient                          &
                          & * outer_product( mode1%unit_vector(j), &
                          &                  mode2%unit_vector(i)  )
      enddo
    enddo
  else
    call print_line(CODE_ERROR//': Trying to construct a dynamical matrix &
       &from two modes at q-points which are neither q1=q2 nor q1=-q2.')
    call err()
  endif
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
  
  allocate( hessian( large_supercell%no_atoms_prim, &
          &          large_supercell%no_atoms),     &
          & stat=ialloc); call err(ialloc)
  hessian = dblemat(zeroes(3,3))
  
  ! Loop across q-points, summing up the contribution from
  !    the dynamical matrix at each.
  do i=1,large_supercell%sc_size
    q = dblevec(qpoints(i)%qpoint)
    do j=1,large_supercell%no_atoms_prim
      atom_1 = large_supercell%atoms(j)
      do k=1,large_supercell%no_atoms
        atom_2 = large_supercell%atoms(k)
        
        r = large_supercell%rvectors(atom_2%rvec_id())
        
        ! Add in the contribution to the Hessian.
        hessian(atom_1%id(),atom_2%id()) =                               &
           &   hessian(atom_1%id(),atom_2%id())                          &
           & + real( dynamical_matrices(i)%elements_( atom_1%prim_id(),  &
           &                                          atom_2%prim_id() ) &
           &       * exp_2pii(-q*r)                                      )
      enddo
    enddo
  enddo
  
  output = CartesianHessian( supercell      = large_supercell, &
                           & elements       = hessian,         &
                           & check_symmetry = .true.,          &
                           & logfile        = logfile          )
end function

! ----------------------------------------------------------------------
! Algebra involving dynamical matrices.
! ----------------------------------------------------------------------
impure elemental function add_DynamicalMatrix_DynamicalMatrix(this,that) &
   & result(output)
  implicit none
  
  type(DynamicalMatrix), intent(in) :: this
  type(DynamicalMatrix), intent(in) :: that
  type(DynamicalMatrix)             :: output
  
  output = DynamicalMatrix(this%elements_ + that%elements_)
end function

impure elemental function negative_DynamicalMatrix(this) result(output)
  implicit none
  
  type(DynamicalMatrix), intent(in) :: this
  type(DynamicalMatrix)             :: output
  
  output = DynamicalMatrix(-this%elements_)
end function

impure elemental function subtract_DynamicalMatrix_DynamicalMatrix(this,that) &
   & result(output)
  implicit none
  
  type(DynamicalMatrix), intent(in) :: this
  type(DynamicalMatrix), intent(in) :: that
  type(DynamicalMatrix)             :: output
  
  output = DynamicalMatrix(this%elements_ - that%elements_)
end function

impure elemental function multiply_DynamicalMatrix_real(this,that) &
   & result(output) 
  implicit none
  
  type(DynamicalMatrix), intent(in) :: this
  real(dp),              intent(in) :: that
  type(DynamicalMatrix)             :: output
  
  output = DynamicalMatrix(this%elements_*that)
end function

impure elemental function multiply_real_DynamicalMatrix(this,that) &
   & result(output) 
  implicit none
  
  real(dp),              intent(in) :: this
  type(DynamicalMatrix), intent(in) :: that
  type(DynamicalMatrix)             :: output
  
  output = DynamicalMatrix(this*that%elements_)
end function

impure elemental function multiply_DynamicalMatrix_complex(this,that) &
   & result(output) 
  implicit none
  
  type(DynamicalMatrix), intent(in) :: this
  complex(dp),           intent(in) :: that
  type(DynamicalMatrix)             :: output
  
  output = DynamicalMatrix(this%elements_*that)
end function

impure elemental function multiply_complex_DynamicalMatrix(this,that) &
   & result(output) 
  implicit none
  
  complex(dp),           intent(in) :: this
  type(DynamicalMatrix), intent(in) :: that
  type(DynamicalMatrix)             :: output
  
  output = DynamicalMatrix(this*that%elements_)
end function

impure elemental function divide_DynamicalMatrix_real(this,that) &
   & result(output) 
  implicit none
  
  type(DynamicalMatrix), intent(in) :: this
  real(dp),              intent(in) :: that
  type(DynamicalMatrix)             :: output
  
  output = DynamicalMatrix(this%elements_/that)
end function

impure elemental function divide_DynamicalMatrix_complex(this,that) &
   & result(output) 
  implicit none
  
  type(DynamicalMatrix), intent(in) :: this
  complex(dp),           intent(in) :: that
  type(DynamicalMatrix)             :: output
  
  output = DynamicalMatrix(this%elements_/that)
end function

impure elemental function conjg_DynamicalMatrix(input) result(output)
  implicit none
  
  type(DynamicalMatrix), intent(in) :: input
  type(DynamicalMatrix)             :: output
  
  output = DynamicalMatrix(conjg(input%elements_))
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
    allocate(this%elements_(no_atoms,no_atoms), stat=ialloc); call err(ialloc)
    k = 0
    do i=1,no_atoms
      do j=1,no_atoms
        k = k+1
        this%elements_(j,i) = ComplexMatrix(elements(k)%strings(2:4))
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
    no_atoms = size(this%elements_,1)
    if (size(this%elements_,2)/=no_atoms) then
      call err()
    endif
    
    allocate(output(5*no_atoms*no_atoms), stat=ialloc); call err(ialloc)
    k = 0
    do i=1,no_atoms
      do j=1,no_atoms
        k = k+1
        matrix_strings = str(this%elements_(j,i))
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
