! ======================================================================
! The third stage of Caesar.
! Uses the forces calculated previously to calculate harmonic normal modes.
! ======================================================================
module calculate_harmonic_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function calculate_harmonic_keywords() result(keywords)
  use keyword_module
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [KeywordData::]
end function

  
! |x> is the collective vector of displacements, {xi}.
! |f> is the collective vector of forces, {fi}.
! Both are supercell%no_modes long.
!
! Under the harmonic approximation, the force constants, F, are defined as:
!    U = - 1/2 <x|F|x>
!
!    |f> = -dU/d|x> = F|x>
!
! Under the symmetry s:
!    |x> -> |xs> = Rs|x>
!    |f> -> |fs> = Rs|f>
!
! F can be found by minimising L:
!    L = sum(s)[ (|fs>-F|xs>)^2 ]
!      = sum(s)[ <fs|fs> - 2<fs|F|xs> + <xs|FF|xs> ]
!
! => 0 = dL/dF = -2 * sum(s)[ |fs><xs| - F|xs><xs| ]
! => F = sum(s)[|fs><xs|] . (sum(s)[|xs><xs|])^-1
!
! sum(s)[|xs><xs|] is block diagonal, so can be inverted in 3x3 blocks.

! ----------------------------------------------------------------------
! Read in |f>, average over +/- |x>, and divide by | |x> |.
! ----------------------------------------------------------------------
function read_forces(supercell,unique_directions,sdir,dft_code,seedname) &
   & result(output)
  use structure_module
  use unique_directions_module
  use dft_output_file_module
  use linear_algebra_module
  implicit none
  
  type(StructureData),    intent(in) :: supercell
  type(UniqueDirections), intent(in) :: unique_directions
  type(String),           intent(in) :: sdir
  type(String),           intent(in) :: dft_code
  type(String),           intent(in) :: seedname
  type(RealVector), allocatable      :: output(:,:)
  
  ! DFT output data.
  type(String)        :: dft_output_filename
  type(DftOutputFile) :: positive
  type(DftOutputFile) :: negative
  
  ! Direction information.
  integer      :: atom
  character(1) :: direction
  
  ! Temporary variables.
  integer          :: i,j,ialloc
  type(RealVector) :: total
  
  dft_output_filename = make_dft_output_filename(dft_code,seedname)
  
  allocate( output(supercell%no_atoms, size(unique_directions)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(unique_directions)
    atom = unique_directions%atoms(i)
    direction = unique_directions%directions_char(i)
    
    positive = read_dft_output_file(dft_code,      &
       & sdir//'/atom.'//atom//'.+d'//direction//'/'//dft_output_filename)
    negative = read_dft_output_file(dft_code,      &
       & sdir//'/atom.'//atom//'.-d'//direction//'/'//dft_output_filename)
    
    do j=1,supercell%no_atoms
      output(j,i) = (positive%forces(j)-negative%forces(j)) / 0.02_dp
    enddo
    
    ! Enforce continuous translational symmetry,
    !    i.e. ensure sum(output,1)=0.
    total = vec([0.0_dp,0.0_dp,0.0_dp])
    do j=1,supercell%no_atoms
      total = total+output(j,i)
    enddo
    do j=1,supercell%no_atoms
      output(j,i) = output(j,i)-total/supercell%no_atoms
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Uses symmetry operations to construct force constants.
! ----------------------------------------------------------------------
function construct_force_constants(forces,supercell,unique_directions, &
   & atom_symmetry_group) result(output)
  use linear_algebra_module
  use structure_module
  use unique_directions_module
  use group_module
  implicit none
  
  type(RealVector),       intent(in) :: forces(:,:)
  type(StructureData),    intent(in) :: supercell
  type(UniqueDirections), intent(in) :: unique_directions
  type(Group),            intent(in) :: atom_symmetry_group(:)
  type(RealMatrix), allocatable      :: output(:,:,:)
  
  ! Atom ids.
  integer :: atom_1,atom_1p,atom_2,atom_2p
  
  ! Rotations in cartesian co-ordinates.
  type(RealMatrix), allocatable :: rotations_cart(:)
  
  ! Parts of |x> and |f>.
  type(RealVector) :: x
  type(RealVector) :: f
  
  ! sum(s)[ |fs><xs| ].
  type(RealMatrix), allocatable :: fx(:,:)
  
  ! sum(s)[ |xs><xs| ] (diagonal blocks only).
  type(RealMatrix), allocatable :: xx(:)
  type(RealMatrix), allocatable :: xx_inverse(:)
  
  ! Force constants, F.
  type(RealMatrix), allocatable :: force_constants(:,:)
  
  ! R-vector information.
  type(Group), allocatable :: rvector_group(:)
  integer                  :: rvector_1
  integer                  :: rvector_2
  integer                  :: rvector
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  ! --------------------------------------------------
  ! Get symmetries in cartesian co-ordinates.
  ! --------------------------------------------------
  rotations_cart = calculate_cartesian_rotations(supercell)
  
  ! --------------------------------------------------
  ! Construct xx and fx.
  ! --------------------------------------------------
  allocate( xx(supercell%no_atoms),                    &
            fx(supercell%no_atoms,supercell%no_atoms), &
          & stat=ialloc); call err(ialloc)
  xx = mat([ 0.0_dp,0.0_dp,0.0_dp, &
           & 0.0_dp,0.0_dp,0.0_dp, &
           & 0.0_dp,0.0_dp,0.0_dp], 3,3)
  fx = mat([ 0.0_dp,0.0_dp,0.0_dp, &
           & 0.0_dp,0.0_dp,0.0_dp, &
           & 0.0_dp,0.0_dp,0.0_dp], 3,3)
  do i=1,supercell%no_symmetries
    do j=1,size(unique_directions)
      atom_1 = unique_directions%atoms(j)
      atom_1p = atom_symmetry_group(i) * atom_1
      
      if (unique_directions%directions_char(j)=='x') then
        x = [ 1.0_dp, 0.0_dp, 0.0_dp ]
      elseif (unique_directions%directions_char(j)=='y') then
        x = [ 0.0_dp, 1.0_dp, 0.0_dp ]
      else
        x = [ 0.0_dp, 0.0_dp, 1.0_dp ]
      endif
      x = rotations_cart(i) * x
      
      xx(atom_1p) = xx(atom_1p) + outer_product(x,x)
      
      do atom_2=1,supercell%no_atoms
        atom_2p = atom_symmetry_group(i) * atom_2
        
        f = rotations_cart(i) * forces(atom_2,j)
        
        fx(atom_1p,atom_2p) = fx(atom_1p,atom_2p) + outer_product(f,x)
      enddo
    enddo
  enddo
  
  ! --------------------------------------------------
  ! Construct xx_inverse.
  ! --------------------------------------------------
  allocate(xx_inverse(supercell%no_atoms),stat=ialloc); call err(ialloc)
  do i=1,supercell%no_atoms
    xx_inverse(i) = invert(xx(i))
  enddo
  
  ! --------------------------------------------------
  ! Construct F.
  ! --------------------------------------------------
  allocate(force_constants(supercell%no_atoms,supercell%no_atoms), &
     & stat=ialloc); call err(ialloc)
  do i=1,supercell%no_atoms
    do j=1,supercell%no_atoms
      force_constants(j,i) = fx(j,i) * xx_inverse(i)
    enddo
  enddo
  
  ! --------------------------------------------------
  ! Symmetrise F.
  ! --------------------------------------------------
  do i=1,supercell%no_atoms
    do j=1,supercell%no_atoms
      force_constants(j,i) = ( force_constants(j,i)            &
                         & + transpose(force_constants(i,j)) &
                         & ) / 2.0_dp
      force_constants(i,j) = transpose(force_constants(i,j))
    enddo
  enddo
  
  ! --------------------------------------------------
  ! Average F across primitive lattice R-vectors, and
  !    convert to a mode-mode-Rvector representation.
  ! --------------------------------------------------
  allocate( output( supercell%sc_size,        &
          &         supercell%no_atoms_prim,  &
          &         supercell%no_atoms_prim), &
          & stat=ialloc); call err(ialloc)
  output = mat([ 0.0_dp,0.0_dp,0.0_dp, &
               & 0.0_dp,0.0_dp,0.0_dp, &
               & 0.0_dp,0.0_dp,0.0_dp  ], 3,3)
  
  rvector_group = calculate_rvector_group(supercell)
  
  do atom_1=1,supercell%no_atoms
    atom_1p = supercell%atom_to_prim(atom_1)
    rvector_1 = supercell%atom_to_rvec(atom_1)
    do atom_2=1,supercell%no_atoms
      atom_2p = supercell%atom_to_prim(atom_2)
      rvector_2 = supercell%atom_to_rvec(atom_2)
      
      rvector = rvector_group(supercell%paired_rvec(rvector_1)) * rvector_2
      
      output(rvector,atom_1p,atom_2p) = output(rvector,atom_1p,atom_2p) &
                                    & + force_constants(atom_1,atom_2)  &
                                    & / supercell%sc_size
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine calculate_harmonic(arguments)
  use utils_module,          only : mkdir
  use linear_algebra_module
  use setup_harmonic_module
  use structure_module
  use dft_output_file_module
  use lte_module
  use unique_directions_module
  use group_module
  use qpoints_module
  use dictionary_module
  use normal_mode_module
  use dynamical_matrix_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Working directory.
  type(String) :: wd
  
  ! Files.
  type(Dictionary)          :: setup_harmonic_arguments
  type(String), allocatable :: no_supercells_file(:)
  
  ! Setup data.
  integer             :: no_supercells
  type(String)        :: dft_code
  type(String)        :: seedname
  type(StructureData) :: structure
  type(StructureData) :: supercell
  
  ! Force constant data.
  type(Group),           allocatable :: atom_symmetry_group(:)
  type(UniqueDirections)             :: unique_directions
  type(RealVector),      allocatable :: forces(:,:)
  type(RealMatrix),      allocatable :: force_constants(:,:,:)
  
  ! q-point data.
  type(StructureData)           :: large_supercell
  type(QpointData), allocatable :: qpoints_ibz(:)
  
  ! Lte output data.
  type(LteReturn)                  :: lte_result
  integer                          :: mode
  integer                          :: atom
  integer                          :: gvector
  
  ! Temporary variables.
  integer                   :: i,j,k
  type(String)              :: sdir,qdir
  
  ! --------------------------------------------------
  ! Read in arguments from user.
  ! --------------------------------------------------
  wd = arguments%value('working_directory')
  
  ! --------------------------------------------------
  ! Read in previous arguments.
  ! --------------------------------------------------
  setup_harmonic_arguments = setup_harmonic_keywords()
  call setup_harmonic_arguments%read_file(wd//'/setup_harmonic.used_settings')
  dft_code = setup_harmonic_arguments%value('dft_code')
  seedname = setup_harmonic_arguments%value('seedname')
  
  no_supercells_file = read_lines(wd//'/no_supercells.dat')
  no_supercells = int(no_supercells_file(1))
  
  structure = read_structure_file(wd//'/structure.dat')
  
  large_supercell = read_structure_file(wd//'/large_supercell.dat')
  
  qpoints_ibz = read_qpoints_file(wd//'/qpoints_ibz.dat')
  
  ! --------------------------------------------------
  ! Loop over supercells
  ! --------------------------------------------------
  do i=1,no_supercells
    sdir = wd//'/Supercell_'//i
    
    ! Read in supercell structure data.
    supercell = read_structure_file(sdir//'/structure.dat')
    
    ! Read in symmetry group and unique atoms.
    atom_symmetry_group = read_group_file(sdir//'/atom_symmetry_group.dat')
    unique_directions = read_unique_directions_file( &
       & sdir//'/unique_directions.dat')
    
    ! Calculate force constants.
    forces = read_forces(supercell,unique_directions,sdir,dft_code, &
       & seedname)
    ! Mass-reduce forces.
    do j=1,size(unique_directions)
      atom = unique_directions%atoms(j)
      do k=1,supercell%no_atoms
        forces(k,j) = forces(k,j) &
                    & / sqrt(supercell%mass(atom)*supercell%mass(k))
      enddo
    enddo
    force_constants = construct_force_constants(forces,supercell, &
       & unique_directions,atom_symmetry_group)
    
    ! Run normal mode analysis.
    lte_result = evaluate_freqs_on_grid(supercell, force_constants)
    
    deallocate(force_constants)
    
    do j=1,size(qpoints_ibz)
      if (qpoints_ibz(j)%sc_id/=i) then
        cycle
      endif
      
      gvector = qpoints_ibz(j)%gvector_id
      
      ! Write out dynamical matrix.
      call write_dynamical_matrix_file(            &
         & lte_result%dynamical_matrices(gvector), &
         & qdir//'/dynamical_matrix.dat')
      
      ! Write out normal modes.
      qdir = wd//'/qpoint_'//j
      call mkdir(qdir)
      do mode=1,structure%no_modes
        call write_normal_mode_file( lte_result%normal_modes(mode,gvector), &
                                   & qdir//'/mode_'//mode//'.dat')
      enddo
    enddo
  enddo
end subroutine
end module
