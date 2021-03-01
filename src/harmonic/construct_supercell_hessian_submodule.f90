submodule (caesar_construct_supercell_hessian_module) caesar_construct_supercell_hessian_submodule
  use caesar_harmonic_module
contains

module procedure construct_supercell_hessian
  ! Electronic structure information.
  type(UniqueDirection),     allocatable :: unique_directions(:)
  type(ElectronicStructure), allocatable :: electronic_structure(:)
  
  ! Files and directories.
  type(IFile)  :: unique_directions_file
  type(String) :: calculation_directory
  
  ! Forces (mass reduced).
  type(RealVector), allocatable :: forces(:,:)
  
  ! sum(s)[ f'^x' ].
  type(RealMatrix), allocatable :: fx(:,:)
  
  ! sum(s)[ x'^x' ] (diagonal blocks only).
  type(RealMatrix), allocatable :: xx(:)
  
  ! Temporary variables.
  integer :: i,ialloc
  
  ! Read in symmetry group and unique atoms.
  unique_directions_file = IFile(directory//'/unique_directions.dat')
  unique_directions = UniqueDirection(unique_directions_file%sections())
  
  ! Read in electronic structure.
  allocate( electronic_structure(size(unique_directions)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(unique_directions)
    calculation_directory =                                    &
       & directory//'/atom.'                                // &
       & left_pad( unique_directions(i)%atom_id,               &
       &           str(maxval(unique_directions%atom_id)) ) // &
       & '.'//unique_directions(i)%direction
    
    electronic_structure(i) = calculation_reader%read_calculation( &
           & calculation_directory,                                &
           & CartesianDisplacement(unique_directions(i),supercell) )
    
    if (size(electronic_structure(i)%forces()) /= supercell%no_atoms) then
      call print_line( ERROR//': Wrong number of forces in '//            &
                     & calculation_directory//'/electronic_structure.dat' )
      call err()
    endif
  enddo
  
  ! Parse forces from electronic structure data.
  forces = parse_forces( supercell,               &
                       & unique_directions,       &
                       & electronic_structure,    &
                       & acoustic_sum_rule_forces )
  
  ! Construct sum[x'^x'].
  xx = construct_xx(unique_directions, supercell, logfile)
  
  ! Construct sum[f'^x'].
  fx = construct_fx(unique_directions, forces, supercell, logfile)
  
  ! Construct F = sum[f'^x'] . inverse(sum[x'^x']).
  output = construct_f(xx,fx,supercell,logfile)
  
  ! Check output.
  call check_hessian( output,            &
                    & forces,            &
                    & supercell,         &
                    & unique_directions, &
                    & logfile            )
end procedure

module procedure parse_forces
  type(CartesianForce), allocatable :: forces
  
  type(RealVector) :: sum_forces
  
  type(AtomData) :: atom_i
  
  integer :: i,j,ialloc
  
  allocate( output(supercell%no_atoms, size(unique_directions)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(unique_directions)
    forces = electronic_structure(i)%forces()
    output(:,i) = forces%vectors
    
    ! Enforce continuous translational symmetry,
    !    i.e. ensure sum(output,1)=0.
    ! This is the 'Accoustic sum rule'.
    if (acoustic_sum_rule_forces) then
      sum_forces = sum(output(:,i))
      do j=1,supercell%no_atoms
        output(j,i) = output(j,i) - sum_forces/supercell%no_atoms
      enddo
    endif
    
    ! Mass-reduce forces.
    atom_i = supercell%atoms(unique_directions(i)%atom_id)
    output(:,i) = output(:,i) / sqrt(supercell%atoms%mass()*atom_i%mass())
  enddo
end procedure

module procedure construct_xx
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
    call print_line('Fractional L2 error: '//sqrt(difference/average))
  endif
end procedure

module procedure construct_fx
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
    call print_line('Fractional L2 error: '//sqrt(difference/average))
  endif
end procedure

module procedure construct_f
  type(RealMatrix), allocatable :: elements(:,:)
  type(RealMatrix)              :: xx_inverse
  
  ! Variables for checking the output.
  type(AtomData) :: atom_1
  type(AtomData) :: atom_2
  
  integer :: i,j,k,k2,ialloc
  
  ! Construct F(i1,i2).
  allocate(elements(supercell%no_atoms_prim,supercell%no_atoms), &
     & stat=ialloc); call err(ialloc)
  do i=1,supercell%no_atoms
    xx_inverse = invert(xx(i))
    do j=1,supercell%no_atoms_prim
      elements(j,i) = fx(j,i) * xx_inverse * supercell%sc_size
    enddo
  enddo
  
  ! Symmetrise F(i1,i2) under exchange of co-ordinates.
  do i=1,supercell%no_atoms_prim
    do j=1,i
      do k=1,supercell%sc_size
        k2 = supercell%paired_rvector_id(k)
        atom_1 = supercell%atoms(i+(k-1)*supercell%no_atoms_prim)
        atom_2 = supercell%atoms(j+(k2-1)*supercell%no_atoms_prim)
        ! F(i,j+R) = F(j,i-R)^T
        elements(atom_1%prim_id(),atom_2%id()) =                   &
           & ( elements(atom_1%prim_id(),atom_2%id())              &
           & + transpose(elements(atom_2%prim_id(),atom_1%id())) ) &
           & / 2
        elements(atom_2%prim_id(),atom_1%id()) = &
           & transpose(elements(atom_1%prim_id(),atom_2%id()))
      enddo
    enddo
  enddo
  
  ! Construct output.
  output = CartesianHessian( supercell      = supercell, &
                           & elements       = elements,  &
                           & check_symmetry = .true.,    &
                           & logfile        = logfile    )
end procedure

module procedure check_hessian
  ! Atoms.
  type(AtomData) :: atom_1
  type(AtomData) :: atom_2
  
  real(dp) :: average
  real(dp) :: difference
  
  type(RealVector) :: calculated
  type(RealVector) :: fitted
  
  integer :: i,j
  
  ! Check Hessian against raw forces.
  average = 0
  difference = 0
  do i=1,size(unique_directions)
    atom_1 = supercell%atoms(unique_directions(i)%atom_id)
    
    do j=1,supercell%no_atoms
      atom_2 = supercell%atoms(j)
      
      calculated = forces(j,i)
      fitted = transpose(hessian%elements(atom_1,atom_2)) &
           & * unique_directions(i)%atomic_displacement   &
           & / supercell%sc_size
      
      average = average + sum_squares((calculated+fitted)/2)
      difference = difference + sum_squares(calculated-fitted)
    enddo
  enddo
  call logfile%print_line('Fractional L2 difference between forces &
     &before and after symmetrisation (this may be large): '// &
     &sqrt(difference/average))
end procedure
end submodule
