! ======================================================================
! Reads harmonic forces, and generates the cartesian Hessian matrix.
! ======================================================================
module construct_supercell_hessian_module
  use common_module
  implicit none
  
  private
  
  public :: construct_supercell_hessian
contains

! ----------------------------------------------------------------------
! Read in the calculated forces, and use them with symmetry operations
!    to construct the Hessian.
! ----------------------------------------------------------------------
! x is the collective vector of displacements, {xi}.
! f is the collective vector of forces, {fi}.
! Both are supercell%no_modes long.
!
! Under the harmonic approximation, the hessian, F, is defined as:
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
impure elemental function construct_supercell_hessian(supercell,directory, &
   & calculation_reader,acoustic_sum_rule_forces,logfile) result(output)
  implicit none
  
  type(StructureData),     intent(in)    :: supercell
  type(String),            intent(in)    :: directory
  type(CalculationReader), intent(inout) :: calculation_reader
  logical,                 intent(in)    :: acoustic_sum_rule_forces
  type(OFile),             intent(inout) :: logfile
  type(CartesianHessian)                 :: output
  
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
end function

! ----------------------------------------------------------------------
! Parse forces from electronic structure data,
!    enforce acoustic sum rules if requires,
!    and mass-reduce forces.
! ----------------------------------------------------------------------
function parse_forces(supercell,unique_directions,electronic_structure, &
   & acoustic_sum_rule_forces) result(output)
  implicit none
  
  type(StructureData),       intent(in) :: supercell
  type(UniqueDirection),     intent(in) :: unique_directions(:)
  type(ElectronicStructure), intent(in) :: electronic_structure(:)
  logical,                   intent(in) :: acoustic_sum_rule_forces
  type(RealVector), allocatable         :: output(:,:)
  
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
    call print_line('Fractional L2 error: '//sqrt(difference/average))
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
    call print_line('Fractional L2 error: '//sqrt(difference/average))
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
  type(CartesianHessian)             :: output
  
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
end function

! ----------------------------------------------------------------------
! Checks the Hessian against the calculated forces and symmetries.
! ----------------------------------------------------------------------
! Checks that the elements corresponding to calculated forces have not been
!    changed too much by symmetrisation.
! Checks that the Hessian has the correct symmetry properties.
subroutine check_hessian(hessian,forces,supercell,unique_directions,logfile)
  implicit none
  
  type(CartesianHessian), intent(in)    :: hessian
  type(RealVector),       intent(in)    :: forces(:,:)
  type(StructureData),    intent(in)    :: supercell
  type(UniqueDirection),  intent(in)    :: unique_directions(:)
  type(OFile),            intent(inout) :: logfile
  
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
end subroutine
end module
