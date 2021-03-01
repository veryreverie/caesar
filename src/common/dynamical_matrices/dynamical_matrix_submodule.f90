submodule (caesar_dynamical_matrix_module) caesar_dynamical_matrix_submodule
  use caesar_dynamical_matrices_module
contains

module procedure new_DynamicalMatrix
  this%elements_ = elements
end procedure

module procedure new_DynamicalMatrix_zeroes
  type(ComplexMatrix) :: zero_matrix
  
  integer :: i,j,ialloc
  
  allocate(this%elements_(no_atoms,no_atoms), stat=ialloc); call err(ialloc)
  ! this%elements_ = cmplxmat(zeroes(3,3))
  ! WORKAROUND for a memory leak in ifort 19.1.0.166.
  zero_matrix = cmplxmat(zeroes(3,3))
  do i=1,size(this%elements_,2)
    do j=1,size(this%elements_,1)
      this%elements_(j,i) = zero_matrix
    enddo
  enddo
end procedure

module procedure elements_DynamicalMatrix
  output = this%elements_
end procedure

module procedure expectation_DynamicalMatrix
  integer :: i,j
  
  output = 0
  do i=1,size(mode%unit_vector)
    do j=1,size(mode%unit_vector)
      output = output + real( conjg(mode%unit_vector(i)) &
                          & * this%elements_(i,j)        &
                          & * mode%unit_vector(j)        )
    enddo
  enddo
end procedure

module procedure new_DynamicalMatrix_interpolated
  type(ComplexMatrix), allocatable :: elements(:,:)
  
  type(AtomData)               :: atom_1
  type(AtomData)               :: atom_2
  type(IntVector)              :: rvector
  type(IntVector), allocatable :: rvectors(:)
  
  integer :: i,j,ialloc
  
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
end procedure

module procedure new_DynamicalMatrix_qpoint
  this = DynamicalMatrix( dblevec(qpoint%qpoint), &
                        & supercell,              &
                        & hessian,                &
                        & min_images              )
end procedure

module procedure check
  type(ComplexMatrix) :: matrix
  type(ComplexMatrix) :: hermitian_matrix
  real(dp)            :: average
  real(dp)            :: difference
  
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
    if (sqrt(difference/average)>1.0e-6_dp) then
      call print_line(WARNING//': Dynamical matrix is not Hermitian. Please &
         &check log files.')
      call print_line('Fractional L2 error in hermicity: '// &
         & sqrt(difference/average))
    endif
  endif
end procedure

module procedure new_ComplexMode_DynamicalMatrix
  if (present(modes_real)) then
    if (modes_real) then
      output = calculate_modes(real(dynamical_matrix%elements()), structure)
      return
    endif
  endif
  
  output = calculate_modes(dynamical_matrix%elements(), structure)
end procedure

module procedure calculate_modes_complex
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
  
  ! Reverse the eigenvalues so that the complex modes are in ascending order
  !    of eigenvalue.
  estuff = estuff(size(estuff):1:-1)
  
  ! Calculate normal modes.
  !output = ComplexMode(estuff, structure)
  ! WORKAROUND for a memory leak in ifort 19.1.0.166.
  allocate(output(size(estuff)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    output(i) = ComplexMode(estuff(i), structure)
  enddo
end procedure

module procedure calculate_modes_real
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
end procedure

module procedure new_DynamicalMatrix_ComplexModes
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
end procedure

module procedure new_DynamicalMatrix_ComplexMode_ComplexMode
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
end procedure

module procedure reconstruct_hessian
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
end procedure

module procedure add_DynamicalMatrix_DynamicalMatrix
  type(ComplexMatrix), allocatable :: elements(:,:)
  
  integer :: i,j
  
  ! output = DynamicalMatrix(this%elements+that%elements)
  ! WORKAROUND to avoid memory leak in ifort 19.1.0.166.
  elements = this%elements_
  do i=1,size(elements,2)
    do j=1,size(elements,1)
      elements(j,i) = elements(j,i) + that%elements_(j,i)
    enddo
  enddo
  output = DynamicalMatrix(elements)
end procedure

module procedure negative_DynamicalMatrix
  type(ComplexMatrix), allocatable :: elements(:,:)
  
  ! output = DynamicalMatrix(-this%elements_)
  ! WORKAROUND to avoid memory leak in ifort 19.1.0.166.
  elements = -this%elements_
  output = DynamicalMatrix(elements)
end procedure

module procedure subtract_DynamicalMatrix_DynamicalMatrix
  output = this + (-that)
end procedure

module procedure multiply_DynamicalMatrix_real
  type(ComplexMatrix), allocatable :: elements(:,:)
  
  integer :: i,j
  
  ! output = DynamicalMatrix(this%elements_*that)
  ! WORKAROUND to avoid memory leak in ifort 19.1.0.166.
  elements = this%elements_
  do i=1,size(elements,2)
    do j=1,size(elements,1)
      elements(j,i) = elements(j,i) * that
    enddo
  enddo
  output = DynamicalMatrix(elements)
end procedure

module procedure multiply_real_DynamicalMatrix
  type(ComplexMatrix), allocatable :: elements(:,:)
  
  integer :: i,j
  
  ! output = DynamicalMatrix(this*that%elements_)
  ! WORKAROUND to avoid memory leak in ifort 19.1.0.166.
  elements = that%elements_
  do i=1,size(elements,2)
    do j=1,size(elements,1)
      elements(j,i) = elements(j,i) * this
    enddo
  enddo
  output = DynamicalMatrix(elements)
end procedure

module procedure multiply_DynamicalMatrix_complex
  type(ComplexMatrix), allocatable :: elements(:,:)
  
  integer :: i,j
  
  ! output = DynamicalMatrix(this%elements_*that)
  ! WORKAROUND to avoid memory leak in ifort 19.1.0.166.
  elements = this%elements_
  do i=1,size(elements,2)
    do j=1,size(elements,1)
      elements(j,i) = elements(j,i) * that
    enddo
  enddo
  output = DynamicalMatrix(elements)
end procedure

module procedure multiply_complex_DynamicalMatrix
  type(ComplexMatrix), allocatable :: elements(:,:)
  
  integer :: i,j
  
  ! output = DynamicalMatrix(this*that%elements_)
  ! WORKAROUND to avoid memory leak in ifort 19.1.0.166.
  elements = that%elements_
  do i=1,size(elements,2)
    do j=1,size(elements,1)
      elements(j,i) = elements(j,i) * this
    enddo
  enddo
  output = DynamicalMatrix(elements)
end procedure

module procedure divide_DynamicalMatrix_real
  type(ComplexMatrix), allocatable :: elements(:,:)
  
  integer :: i,j
  
  ! output = DynamicalMatrix(this%elements_*that)
  ! WORKAROUND to avoid memory leak in ifort 19.1.0.166.
  elements = this%elements_
  do i=1,size(elements,2)
    do j=1,size(elements,1)
      elements(j,i) = elements(j,i) / that
    enddo
  enddo
  output = DynamicalMatrix(elements)
end procedure

module procedure divide_DynamicalMatrix_complex
  type(ComplexMatrix), allocatable :: elements(:,:)
  
  integer :: i,j
  
  ! output = DynamicalMatrix(this%elements_*that)
  ! WORKAROUND to avoid memory leak in ifort 19.1.0.166.
  elements = this%elements_
  do i=1,size(elements,2)
    do j=1,size(elements,1)
      elements(j,i) = elements(j,i) / that
    enddo
  enddo
  output = DynamicalMatrix(elements)
end procedure

module procedure conjg_DynamicalMatrix
  type(ComplexMatrix), allocatable :: elements(:,:)
  
  ! elements = DynamicalMatrix(conjg(input%elements_))
  ! WORKAROUND to avoid memory leak in ifort 19.1.0.166.
  elements = input%elements_
  elements = conjg(elements)
  output = DynamicalMatrix(elements)
end procedure

module procedure read_DynamicalMatrix
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
end procedure

module procedure write_DynamicalMatrix
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
end procedure

module procedure new_DynamicalMatrix_Strings
  call this%read(input)
end procedure

module procedure new_DynamicalMatrix_StringArray
  this = DynamicalMatrix(str(input))
end procedure
end submodule
