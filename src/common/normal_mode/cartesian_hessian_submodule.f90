submodule (caesar_cartesian_hessian_module) caesar_cartesian_hessian_submodule
  use caesar_normal_mode_module
contains

module procedure new_CartesianHessian_elements
  ! Variables for checking Hessian.
  type(AtomData)   :: atom_i
  type(AtomData)   :: atom_j
  type(AtomData)   :: atom_i2
  type(AtomData)   :: atom_j2
  type(RealMatrix) :: matrix
  type(RealMatrix) :: symmetric
  real(dp)         :: average
  real(dp)         :: difference
  
  integer :: i,j,k,k2
  
  ! Check elements array has the right shape.
  if (any(shape(elements)/=[supercell%no_atoms_prim,supercell%no_atoms])) then
    call print_line(CODE_ERROR//': Hessian elements must have dimension &
       &no_atoms_prim*no_atoms.')
    call err()
  endif
  
  ! Copy elements into CartesianHessian.
  this%elements_ = elements
  
  ! Check Hessian is symmetric, i.e. F(i,j+R)=F(j,i-R)^T
  average = 0
  difference = 0
  do i=1,supercell%no_atoms_prim
    do j=1,supercell%no_atoms_prim
      do k=1,supercell%sc_size
        k2 = supercell%paired_rvector_id(k)
        matrix = this%elements(                               &
           & supercell%atoms(i),                              &
           & supercell%atoms(j+(k-1)*supercell%no_atoms_prim) )
        symmetric = transpose(this%elements(                   &
           & supercell%atoms(j),                               &
           & supercell%atoms(i+(k2-1)*supercell%no_atoms_prim) ))
        average = average + sum_squares((matrix+symmetric)/2)
        difference = difference + sum_squares(matrix-symmetric)
      enddo
    enddo
  enddo
  if (present(logfile)) then
    call logfile%print_line(                                  &
       & 'Fractional L2 error in F(i1,i2)=F(i2,i1)      : '// &
       & sqrt(difference/average))
  endif
  if (sqrt(difference/average) > 1e-6_dp) then
    call print_line(WARNING//': F(i1,i2)/=F(i2,i1). Please check log files.')
    call print_line('Fractional L2 error in F(i1,i2)=F(i2,i1) : '// &
       & sqrt(difference/average))
  endif
  
  ! Check F(i1,i2) transforms correctly under symmetry operators.
  if (check_symmetry) then
    average = 0
    difference = 0
    do i=1,size(supercell%symmetries)
      do j=1,supercell%no_atoms_prim
        atom_i = supercell%atoms(j)
        atom_i2 = supercell%atoms( supercell%symmetries(i)%atom_group &
                               & * atom_i%id())
        if (atom_i2%id()/=atom_i2%prim_id()) then
          cycle
        endif
        do k=1,supercell%no_atoms
          atom_j = supercell%atoms(k)
          atom_j2 = supercell%atoms( supercell%symmetries(i)%atom_group &
                                 & * atom_j%id())
          matrix = this%elements(atom_i2,atom_j2)
          symmetric = supercell%symmetries(i)%cartesian_tensor &
                  & * this%elements(atom_i,atom_j)             &
                  & * transpose(supercell%symmetries(i)%cartesian_tensor)
          average = average + sum_squares((matrix+symmetric)/2)
          difference = difference + sum_squares(matrix-symmetric)
        enddo
      enddo
    enddo
    if (present(logfile)) then
      call logfile%print_line(                                  &
         & 'Fractional L2 error in symmetry of F(i1,i2)   : '// &
         & sqrt(difference/average))
    endif
    if (sqrt(difference/average) > 1e-6_dp) then
      call print_line(WARNING//': F(i1,i2) is not as symmetric as expected. &
         &Please check log files.')
      call print_line('Fractional L2 error in symmetry: '//&
         &sqrt(difference/average))
    endif
  endif
end procedure

module procedure elements
  if (a%id()/=a%prim_id()) then
    call print_line(CODE_ERROR//": The first atom in a call to a Hessian's &
       &elements() module procedure must be in the primitive cell.")
    call err()
  endif
  
  output = this%elements_(a%id(),b%id())
end procedure

module procedure read_CartesianHessian
  type(StringArray), allocatable :: elements(:)
  integer,           allocatable :: is(:)
  integer,           allocatable :: js(:)
  integer                        :: no_atoms_prim
  integer                        :: no_atoms
  
  integer :: i,j,k,ialloc
  
  select type(this); type is(CartesianHessian)
    elements = split_into_sections(input)
    js = [( int(token(elements(i)%strings(1), 3)), &
          & i=1,                                   &
          & size(elements)                         )]
    is = [( int(token(elements(i)%strings(1), 4)), &
          & i=1, &
          & size(elements))]
    no_atoms_prim = maxval(js)
    no_atoms = maxval(is)
    no_atoms = int_sqrt(size(elements))
    allocate( this%elements_(no_atoms_prim,no_atoms), &
            & stat=ialloc); call err(ialloc)
    k = 0
    do i=1,no_atoms
      do j=1,no_atoms_prim
        k = k+1
        this%elements_(j,i) = RealMatrix(elements(k)%strings(2:4))
      enddo
    enddo
  class default
    call err()
  end select
end procedure

module procedure write_CartesianHessian
  integer :: no_atoms_prim
  integer :: no_atoms
  
  type(String) :: matrix_strings(3)
  
  integer :: i,j,k,ialloc
  
  select type(this); type is(CartesianHessian)
    no_atoms_prim = size(this%elements_,1)
    no_atoms = size(this%elements_,2)
    
    allocate(output(5*no_atoms_prim*no_atoms), stat=ialloc); call err(ialloc)
    k = 0
    do i=1,no_atoms
      do j=1,no_atoms_prim
        k = k+1
        matrix_strings = str(this%elements_(j,i))
        output(5*k-4) = 'Atoms: ( '//j//' '//i//' )'
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

module procedure new_CartesianHessian_Strings
  call this%read(input)
end procedure

module procedure new_CartesianHessian_StringArray
  this = CartesianHessian(str(input))
end procedure
end submodule
