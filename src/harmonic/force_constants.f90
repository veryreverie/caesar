! ======================================================================
! Reads harmonic forces, and generates the matrix of force constants.
! ======================================================================
module force_constants_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
  
  private
  
  public :: construct_force_constants
  
  interface sum_squares
    module procedure sum_squares_RealVector
    module procedure sum_squares_RealMatrix
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
function construct_force_constants(supercell,unique_directions,sdir, &
   & file_type,seedname,log_filename) result(output)
  use linear_algebra_module
  use structure_module
  use unique_directions_module
  use group_module
  use atom_module
  use ofile_module
  implicit none
  
  type(StructureData),   intent(in) :: supercell
  type(UniqueDirection), intent(in) :: unique_directions(:)
  type(String),          intent(in) :: sdir
  type(String),          intent(in) :: file_type
  type(String),          intent(in) :: seedname
  type(String),          intent(in) :: log_filename
  type(RealMatrix), allocatable     :: output(:,:,:)
  
  ! Forces (mass reduced).
  type(RealVector), allocatable :: forces(:,:)
  
  ! sum(s)[ f'^x' ].
  type(RealMatrix), allocatable :: fx(:,:)
  
  ! sum(s)[ x'^x' ] (diagonal blocks only).
  type(RealMatrix), allocatable :: xx(:)
  
  ! Force constants, F.
  type(RealMatrix), allocatable :: force_constants(:,:)
  
  ! Log file.
  type(OFile) :: logfile
  
  ! Open logfile.
  logfile = log_filename
  
  ! Read in forces (mass reduced).
  forces = read_forces( supercell,         &
                      & unique_directions, &
                      & sdir,              &
                      & file_type,         &
                      & seedname)
  
  ! Construct sum[x'^x'].
  xx = construct_xx(unique_directions, supercell, logfile)
  
  ! Construct sum[f'^x'].
  fx = construct_fx(unique_directions, forces, supercell, logfile)
  
  ! Construct F = sum[f'^x'] . inverse(sum[x'^x'])
  force_constants = construct_f(xx,fx,supercell,logfile)
  
  ! Average F across R-vectors.
  ! F((R1,i1),(R2,i2)) -> F(R,i1,i2) where R=R2-R1.
  output = average_f(force_constants,supercell,logfile)
  
  call check_force_constants( forces,                &
                            & output,                &
                            & supercell,             &
                            & unique_directions,     &
                            & logfile)
end function
  
! ----------------------------------------------------------------------
! Read in forces, and mass-weight them.
! ----------------------------------------------------------------------
function read_forces(supercell,unique_directions,sdir,file_type,seedname) &
   & result(output)
  use structure_module
  use unique_directions_module
  use output_file_module
  use linear_algebra_module
  use atom_module
  implicit none
  
  type(StructureData),   intent(in) :: supercell
  type(UniqueDirection), intent(in) :: unique_directions(:)
  type(String),          intent(in) :: sdir
  type(String),          intent(in) :: file_type
  type(String),          intent(in) :: seedname
  type(RealVector), allocatable     :: output(:,:)
  
  ! DFT output data.
  type(String)     :: output_filename
  type(OutputFile) :: output_file
  
  ! Direction information.
  type(AtomData) :: atom
  type(String)   :: direction
  type(String)   :: atom_string
  
  ! Temporary variables.
  integer          :: i,j,ialloc
  type(RealVector) :: total
  
  output_filename = make_output_filename(file_type,seedname)
  
  allocate( output(supercell%no_atoms, size(unique_directions)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(unique_directions)
    atom = supercell%atoms(unique_directions(i)%atom_id)
    direction = unique_directions(i)%direction
    atom_string = left_pad(atom%id(),str(supercell%no_atoms_prim))
    
    output_file = read_output_file(                                           &
       & file_type,                                                           &
       & sdir//'/atom.'//atom_string//'.'//direction//'/'//output_filename, &
       & supercell)
    
    output(:,i) = output_file%forces(:)
    
    ! Enforce continuous translational symmetry,
    !    i.e. ensure sum(output,1)=0.
    ! This is the 'Accoustic sum rule'.
    total = dble(zeroes(3))
    do j=1,supercell%no_atoms
      total = total+output(j,i)
    enddo
    do j=1,supercell%no_atoms
      output(j,i) = output(j,i)-total/supercell%no_atoms
    enddo
    
    ! Mass-reduce forces.
    do j=1,supercell%no_atoms
      output(j,i) = output(j,i) / sqrt(atom%mass()*supercell%atoms(j)%mass())
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Calculates the sum of squares of the elements of a matrix.
! ----------------------------------------------------------------------
function sum_squares_RealVector(input) result(output)
  use linear_algebra_module
  implicit none
  
  type(RealVector), intent(in) :: input
  real(dp)                     :: output
  
  output = sum(dble(input)*dble(input))
end function

function sum_squares_RealMatrix(input) result(output)
  use linear_algebra_module
  implicit none
  
  type(RealMatrix), intent(in) :: input
  real(dp)                     :: output
  
  output = sum(dble(input)*dble(input))
end function

! ----------------------------------------------------------------------
! Construct and check sum[x'^x'].
! ----------------------------------------------------------------------
function construct_xx(unique_directions,supercell,logfile) result(output)
  use unique_directions_module
  use structure_module
  use linear_algebra_module
  use atom_module
  use ofile_module
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
  output = mat(dble(zeroes(3,3)))
  do i=1,size(supercell%symmetries)
    do j=1,size(unique_directions)
      atom_1 = supercell%atoms(unique_directions(j)%atom_id)
      atom_1p = supercell%atoms( supercell%symmetries(i)%atom_group &
                             & * atom_1%id())
      
      x = supercell%symmetries(i)%cartesian_rotation &
      & * unique_directions(j)%displacement
      
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
      symmetric = supercell%symmetries(i)%cartesian_rotation &
              & * output(atom_1%id())                        &
              & * transpose(supercell%symmetries(i)%cartesian_rotation)
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
  use unique_directions_module
  use structure_module
  use linear_algebra_module
  use atom_module
  use ofile_module
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
  type(RealMatrix) :: matrix
  type(RealMatrix) :: symmetric
  real(dp)         :: average
  real(dp)         :: difference
  
  ! Temporary variables.
  integer :: i,j,k,ialloc

  ! Construct sum[f'^x'].
  allocate( output(supercell%no_atoms,supercell%no_atoms), &
          & stat=ialloc); call err(ialloc)
  output = mat(dble(zeroes(3,3)))
  do i=1,size(supercell%symmetries)
    do j=1,size(unique_directions)
      atom_1 = supercell%atoms(unique_directions(j)%atom_id)
      atom_1p = supercell%atoms( supercell%symmetries(i)%atom_group &
                             & * atom_1%id())
      
      x = supercell%symmetries(i)%cartesian_rotation &
      & * unique_directions(j)%displacement
      
      do k=1,supercell%no_atoms
        atom_2 = supercell%atoms(k)
        atom_2p = supercell%atoms( supercell%symmetries(i)%atom_group &
                               & * atom_2%id())
        
        f = supercell%symmetries(i)%cartesian_rotation * forces(atom_2%id(),j)
        
        output(atom_1p%id(),atom_2p%id()) = output(atom_1p%id(),atom_2p%id()) &
                                        & + outer_product(f,x)
      enddo
    enddo
  enddo
  
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
        symmetric = supercell%symmetries(i)%cartesian_rotation &
                & * output(atom_1%id(),atom_2%id())            &
                & * transpose(supercell%symmetries(i)%cartesian_rotation)
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
! Construct and check F = sum[f'^x'] . inverse(sum[x'^x']).
! ----------------------------------------------------------------------
! sum[x'^x'] is block diagonal, so only the 3x3 diagonal blocks need inverting.
function construct_f(xx,fx,supercell,logfile) result(output)
  use linear_algebra_module
  use structure_module
  use ofile_module
  use atom_module
  implicit none
  
  type(RealMatrix),    intent(in)    :: xx(:)
  type(RealMatrix),    intent(in)    :: fx(:,:)
  type(StructureData), intent(in)    :: supercell
  type(OFile),         intent(inout) :: logfile
  type(RealMatrix), allocatable      :: output(:,:)
  
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
  allocate(output(supercell%no_atoms,supercell%no_atoms), &
     & stat=ialloc); call err(ialloc)
  do i=1,supercell%no_atoms
    call print_line(xx(i))
    xx_inverse = invert(xx(i))
    do j=1,supercell%no_atoms
      output(j,i) = fx(j,i) * xx_inverse
    enddo
  enddo
  
  ! Symmetrise F(i1,i2) under exchange of co-ordinates.
  do i=1,supercell%no_atoms
    do j=1,i
      output(j,i) = (output(j,i) + transpose(output(i,j))) / 2.0_dp
      output(i,j) = transpose(output(j,i))
    enddo
  enddo
  
  ! Check F(i1,i2) obeys correct symmetries.
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
        matrix = output(atom_1p%id(),atom_2p%id())
        symmetric = supercell%symmetries(i)%cartesian_rotation &
                & * output(atom_1%id(),atom_2%id())            &
                & * transpose(supercell%symmetries(i)%cartesian_rotation)
        average = average + sum_squares((matrix+symmetric)/2)
        difference = difference + sum_squares(matrix-symmetric)
      enddo
    enddo
  enddo
  call logfile%print_line( &
     & 'Fractional L2 error in symmetry of F(i1,i2)   : '// &
     & sqrt(difference/average))
  if (sqrt(difference/average) > 1e-10_dp) then
    call print_line(WARNING//': F(i1,i2) is not as symmetric as expected. &
       &Please check log files.')
  endif
end function

! ----------------------------------------------------------------------
! Convert F from atom-atom representation to rvec-prim-prim representation.
! ----------------------------------------------------------------------
! Due to translational symmetry, F((R1,i1),(R2,i2)) = F((0,i1),(R2-R1,i2)),
!    where R1 and i1 are the R-vector and prim atom of atom 1 respectively.
! This function averages across R1, to produce F(R,i1,i2), where R=R2-R1.
function average_f(f,supercell,logfile) result(output)
  use linear_algebra_module
  use structure_module
  use ofile_module
  use atom_module
  use group_module
  implicit none
  
  type(RealMatrix),    intent(in)    :: f(:,:)
  type(StructureData), intent(in)    :: supercell
  type(OFile),         intent(inout) :: logfile
  type(RealMatrix), allocatable      :: output(:,:,:)
  
  ! Atoms.
  type(AtomData) :: atom_1
  type(AtomData) :: atom_2
  
  ! The group s.t. group(i) * j = k if Ri+Rj=Rk.
  type(Group), allocatable :: rvector_group(:)
  
  ! The id of R-vector R = R2-R1.
  integer :: rvector
  
  ! Variables for checking the output.
  type(AtomData)   :: atom_1p
  type(AtomData)   :: atom_2p
  integer          :: rvector_p
  type(RealMatrix) :: matrix
  type(RealMatrix) :: symmetric
  real(dp)         :: average
  real(dp)         :: difference
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  ! Construct F(R,i1,i2).
  allocate( output( supercell%sc_size,        &
          &         supercell%no_atoms_prim,  &
          &         supercell%no_atoms_prim), &
          & stat=ialloc); call err(ialloc)
  output = mat(dble(zeroes(3,3)))
  rvector_group = supercell%calculate_rvector_group()
  do i=1,supercell%no_atoms
    atom_1 = supercell%atoms(i)
    do j=1,supercell%no_atoms
      atom_2 = supercell%atoms(j)
      
      ! Find the R-vector from atom 1 to atom 2.
      ! R = R2 - R1 = R2 + (G-R1) = R2 + R1p.
      rvector = rvector_group(supercell%paired_rvec(atom_1%rvec_id())) &
            & * atom_2%rvec_id()
      
      output(rvector,atom_1%prim_id(),atom_2%prim_id()) =    &
         & output(rvector,atom_1%prim_id(),atom_2%prim_id()) &
         & + f(atom_1%id(),atom_2%id())                      &
         & / supercell%sc_size
    enddo
  enddo
  
  ! Check F(R,i1,i2) obeys correct symmetries.
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
        rvector = rvector_group(supercell%paired_rvec(atom_1%rvec_id())) &
              & * atom_2%rvec_id()
        rvector_p = rvector_group(supercell%paired_rvec(atom_1p%rvec_id())) &
                & * atom_2p%rvec_id()
        matrix = output(rvector_p,atom_1p%prim_id(),atom_2p%prim_id())
        symmetric = supercell%symmetries(i)%cartesian_rotation        &
                & * output(rvector,atom_1%prim_id(),atom_2%prim_id()) &
                & * transpose(supercell%symmetries(i)%cartesian_rotation)
        average = average + sum_squares((matrix+symmetric)/2)
        difference = difference + sum_squares(matrix-symmetric)
      enddo
    enddo
  enddo
  call logfile%print_line(                                  &
     & 'Fractional L2 error in symmetry of F(R,i1,i2) : '// &
     & sqrt(difference/average))
  if (sqrt(difference/average) > 1e-10_dp) then
    call print_line(WARNING//': F(R,i1,i2) is not as symmetric as expected. &
       &Please check log files.')
  endif
end function

! ----------------------------------------------------------------------
! Checks the force constants against the calculated forces and symmetries.
! ----------------------------------------------------------------------
! Checks that the elements corresponding to calculated forces have not been
!    changed too much by symmetrisation.
! Checks that the force constants have the correct symmetry properties.
subroutine check_force_constants(forces,force_constants,supercell, &
   & unique_directions,logfile)
  use linear_algebra_module
  use structure_module
  use unique_directions_module
  use group_module
  use ofile_module
  use atom_module
  implicit none
  
  type(RealVector),      intent(in)    :: forces(:,:)
  type(RealMatrix),      intent(in)    :: force_constants(:,:,:)
  type(StructureData),   intent(in)    :: supercell
  type(UniqueDirection), intent(in)    :: unique_directions(:)
  type(OFile),           intent(inout) :: logfile
  
  ! Atoms.
  type(AtomData) :: atom_1
  type(AtomData) :: atom_2
  
  integer :: rvector
  
  real(dp) :: average
  real(dp) :: difference
  
  type(RealVector) :: calculated
  type(RealVector) :: fitted
  
  type(Group), allocatable :: rvector_group(:)
  
  integer :: i,j
  
  rvector_group = supercell%calculate_rvector_group()
  
  ! Check force constants against calculated forces.
  average = 0
  difference = 0
  do i=1,size(unique_directions)
    atom_1 = supercell%atoms(unique_directions(i)%atom_id)
    
    do j=1,supercell%no_atoms
      atom_2 = supercell%atoms(j)
      
      rvector = rvector_group(supercell%paired_rvec(atom_1%rvec_id())) &
            & * atom_2%rvec_id()
      
      calculated = forces(j,i)
      fitted = force_constants(rvector,atom_1%prim_id(),atom_2%prim_id()) &
           & * unique_directions(i)%displacement
      
      average = average + sum_squares((calculated+fitted)/2)
      difference = difference + sum_squares(calculated-fitted)
    enddo
  enddo
  call logfile%print_line('Fractional L2 difference between forces &
     &before and after symmetrisation: '//sqrt(difference/average))
end subroutine
end module
