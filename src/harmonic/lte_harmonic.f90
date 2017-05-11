! ======================================================================
! Reads dft output files in sdir, and generates matrix of force constants.
! ======================================================================
module lte_harmonic_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function lte_harmonic_keywords() result(keywords)
  use help_module
  implicit none
  
  type(KeywordData) :: keywords(1)
  
  keywords = [ &
  & make_keyword('temperature', '0', 'temperature is the temperature in &
      &Kelvin, used when calculating the density of states and phonon &
      &dispersion curve.') ]
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
  use utils_module, only : make_dft_output_filename
  use structure_module
  use unique_directions_module
  use dft_output_file_module
  implicit none
  
  type(StructureData),    intent(in) :: supercell
  type(UniqueDirections), intent(in) :: unique_directions
  type(String),           intent(in) :: sdir
  type(String),           intent(in) :: dft_code
  type(String),           intent(in) :: seedname
  real(dp), allocatable              :: output(:,:,:)
  
  ! DFT output data.
  type(String)        :: dft_output_filename
  type(DftOutputFile) :: positive
  type(DftOutputFile) :: negative
  
  ! Direction information.
  integer      :: atom
  character(1) :: direction
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  dft_output_filename = make_dft_output_filename(dft_code,seedname)
  
  allocate( output(3,supercell%no_atoms,size(unique_directions)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(unique_directions)
    atom = unique_directions%atoms(i)
    direction = unique_directions%directions_char(i)
    
    positive = read_dft_output_file(dft_code,      &
       & sdir//'/atom.'//atom//'.+d'//direction//'/'//dft_output_filename)
    negative = read_dft_output_file(dft_code,      &
       & sdir//'/atom.'//atom//'.-d'//direction//'/'//dft_output_filename)
    
    output(:,:,i) = (positive%forces-negative%forces) / 0.02_dp
    
    ! Enforce continuous translational symmetry,
    !    i.e. ensure sum(output,1)=0.
    do j=1,3
      output(j,:,i) = output(j,:,i) - sum(output(j,:,i))/size(output,2)
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Uses symmetry operations to construct force constants.
! ----------------------------------------------------------------------
function construct_force_constants(forces,supercell,unique_directions, &
   & symmetry_group) result(output)
  use linear_algebra_module
  use structure_module
  use unique_directions_module
  use group_module
  implicit none
  
  real(dp),               intent(in) :: forces(:,:,:)
  type(StructureData),    intent(in) :: supercell
  type(UniqueDirections), intent(in) :: unique_directions
  type(Group),            intent(in) :: symmetry_group(:)
  real(dp), allocatable              :: output(:,:,:)
  
  ! Atom ids.
  integer :: atom_1,atom_1p,atom_2,atom_2p
  integer :: mode_1, mode_2
  
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
      atom_1p = operate(symmetry_group(i), atom_1)
      
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
        atom_2p = operate(symmetry_group(i), atom_2)
        
        f = rotations_cart(i) * vec(forces(:,atom_2,j))
        
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
  allocate( output( supercell%sc_size,                        &
          &         supercell%no_modes/supercell%sc_size,  &
          &         supercell%no_modes/supercell%sc_size), &
          & stat=ialloc); call err(ialloc)
  output = 0.0_dp
  
  rvector_group = calculate_rvector_group(supercell)
  
  do atom_1=1,supercell%no_atoms
    atom_1p = supercell%atom_to_prim(atom_1)
    mode_1 = (atom_1p-1)*3 + 1
    rvector_1 = supercell%atom_to_rvec(atom_1)
    do atom_2=1,supercell%no_atoms
      atom_2p = supercell%atom_to_prim(atom_2)
      mode_2 = (atom_2p-1)*3 + 1
      rvector_2 = supercell%atom_to_rvec(atom_2)
      
      rvector = operate( rvector_group(supercell%paired_rvec(rvector_1)), &
                       & rvector_2)
      
      output(rvector,mode_1:mode_1+2,mode_2:mode_2+2) = &
         &   output(rvector,mode_1:mode_1+2,mode_2:mode_2+2) &
         & + dble(force_constants(atom_1,atom_2))
    enddo
  enddo
  
  output = output / supercell%sc_size
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine lte_harmonic(arguments)
  use utils_module,          only : mkdir
  use linear_algebra_module, only : invert
  use structure_module
  use dft_output_file_module
  use lte_module
  use unique_directions_module
  use group_module
  use qpoints_module
  use dictionary_module
  use normal_mode_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Working directory.
  type(String) :: wd
  
  ! User-input temperature
  real(dp) :: temperature
  
  ! File contents
  type(Dictionary)          :: setup_harmonic_arguments
  type(String), allocatable :: no_supercells_file(:)
  
  ! Setup data.
  integer             :: no_supercells
  type(String)        :: dft_code
  type(String)        :: seedname
  type(StructureData) :: structure
  type(StructureData) :: supercell
  
  ! Force constant data.
  type(Group),           allocatable :: symmetry_group(:)
  type(UniqueDirections)             :: unique_directions
  real(dp),              allocatable :: forces(:,:,:)
  real(dp),              allocatable :: force_constants(:,:,:)
  integer                            :: forces_file
  integer                            :: force_constants_file
  
  ! q-point data.
  type(StructureData)           :: structure_grid
  type(QpointData), allocatable :: qpoints_ibz(:)
  
  ! lte input data.
  integer               :: no_kspace_lines
  real(dp), allocatable :: disp_qpoints(:,:)
  
  ! Lte output data.
  type(LteReturn)          :: lte_result
  complex(dp), allocatable :: ibz_dynamical_matrices(:,:,:)
  integer                  :: mode
  integer                  :: atom
  
  ! Normal mode data.
  integer          :: gvector
  type(NormalMode) :: normal_mode
  
  ! Temporary variables.
  integer                   :: i,j,k,l,ialloc
  type(String)              :: sdir,qdir
  
  ! --------------------------------------------------
  ! Read in arguments from user.
  ! --------------------------------------------------
  wd = item(arguments, 'working_directory')
  temperature = dble(item(arguments, 'temperature'))
  
  ! --------------------------------------------------
  ! Read in previous arguments.
  ! --------------------------------------------------
  setup_harmonic_arguments = read_dictionary_file( &
     & wd//'/setup_harmonic.used_settings')
  dft_code = item(setup_harmonic_arguments, 'dft_code')
  seedname = item(setup_harmonic_arguments, 'seedname')
  
  no_supercells_file = read_lines(wd//'/no_sc.dat')
  no_supercells = int(no_supercells_file(1))
  
  structure = read_structure_file(wd//'/structure.dat')
  
  ! Read q-points and G-vectors.
  structure_grid = read_structure_file(wd//'/structure_grid.dat')
  qpoints_ibz = read_qpoints_file(wd//'/qpoints_ibz.dat')
  
  allocate( ibz_dynamical_matrices( structure%no_modes, &
          &                         structure%no_modes, &
          &                         size(qpoints_ibz)), &
          & stat=ialloc); call err(ialloc)
  
  ! --------------------------------------------------
  ! Loop over supercells
  ! --------------------------------------------------
  do i=1,no_supercells
    sdir = wd//'/Supercell_'//i
    
    ! Read in supercell structure data.
    supercell = read_structure_file(sdir//'/structure.dat')
    
    ! Read in symmetry group and unique atoms.
    symmetry_group = read_group_file(sdir//'/symmetry_group.dat')
    unique_directions = read_unique_directions_file( &
       & sdir//'/unique_directions.dat')
    
    ! Calculate force constants.
    forces = read_forces(supercell,unique_directions,sdir,dft_code, &
       & seedname)
    ! Mass-reduce forces.
    do j=1,size(forces,3)
      atom = unique_directions%atoms(j)
      do k=1,size(forces,2)
        forces(:,k,j) = forces(:,k,j) &
                    & / sqrt(supercell%mass(atom)*supercell%mass(k))
      enddo
    enddo
    force_constants = construct_force_constants(forces,supercell, &
       & unique_directions,symmetry_group)
    
    ! Write out forces and force constants for debugging purposes.
    forces_file = open_write_file(sdir//'/forces.dat')
    do j=1,size(forces,3)
      do k=1,structure%no_modes
        call print_line(forces_file, &
           & 'Force constant on mode '//k//' as a result of mode '// &
           & unique_directions%modes(j)//' being displaced.')
        do l=1,supercell%sc_size
          atom = supercell%rvec_and_prim_to_atom((k-1)/3+1,l)
          call print_line(forces_file, 'R-vector '//k//': &
             &acceleration/displacement = '//forces(modulo(k-1,3)+1,atom,j))
        enddo
        call print_line(forces_file,'')
      enddo
    enddo
    close(forces_file)
    
    force_constants_file = open_write_file(sdir//'/force_constants.dat')
    do j=1,size(force_constants,3)
      do k=1,size(force_constants,2)
        call print_line(force_constants_file, 'Symmetrised force constant on &
           & mode '//j//' as a result of mode '//k//' being displaced.')
        do l=1,size(force_constants,1)
          call print_line(force_constants_file, &
             & 'R-vector '//l//': acceleration/displacement = '// &
             & force_constants(l,k,j))
        enddo
        call print_line(force_constants_file, '')
      enddo
    enddo
    
    ! Run normal mode analysis.
    lte_result = evaluate_freqs_on_grid(supercell, force_constants)
    
    deallocate(force_constants)
    
    do j=1,size(qpoints_ibz)
      if (qpoints_ibz(j)%sc_id/=i) then
        cycle
      endif
      
      gvector = qpoints_ibz(j)%gvector_id
      
      ! Move dynamical matrices into ibz_dynamical matrices.
      ibz_dynamical_matrices(:,:,j) = &
         & lte_result%dynamical_matrices(:,:,gvector)
      
      ! Write out normal modes.
      qdir = wd//'/qpoint_'//j
      call mkdir(qdir)
      do mode=1,structure%no_modes
        call new(normal_mode, supercell%no_atoms)
        normal_mode%frequency = lte_result%frequencies(mode, gvector)
        normal_mode%soft_mode = lte_result%soft_modes(mode, gvector)
        normal_mode%displacements = lte_result%displacements(:,:,mode, gvector)
        call write_normal_mode_file(normal_mode, qdir//'/mode_'//mode//'.dat')
      enddo
    enddo
  enddo
  
!  ! Write path for fourier interpolation
!  no_kspace_lines = 4
!  allocate(disp_qpoints(3,no_kspace_lines+1))
!  disp_qpoints(:,1) = [ 0.0_dp, 0.0_dp, 0.0_dp ] ! GM
!  disp_qpoints(:,2) = [ 0.5_dp, 0.5_dp, 0.5_dp ] ! T
!  disp_qpoints(:,3) = [ 0.0_dp, 0.5_dp, 0.5_dp ] ! FB
!  disp_qpoints(:,4) = [ 0.0_dp, 0.0_dp, 0.0_dp ] ! GM
!  disp_qpoints(:,5) = [ 0.0_dp, 0.5_dp, 0.0_dp ] ! L
!  
!  ! Read in primitive symmetry group.
!  symmetry_group = read_group_file(wd//'/Supercell_1/symmetry_group.dat')
!  
!  call print_line('')
!  call print_line('Running fourier interpolation (this may take some time).')
!  call fourier_interpolation(              &
!     & ibz_dynamical_matrices,             &
!     & structure,                          &
!     & temperature,                        &
!     & structure_grid,                     &
!     & qpoints_ibz,                        &
!     & disp_qpoints,                       &
!     & symmetry_group,                     &
!     & wd//'/phonon_dispersion_curve.dat', &
!     & wd//'/high_symmetry_points.dat',    &
!     & wd//'/free_energy.dat',             &
!     & wd//'/freq_dos.dat')
end subroutine
end module
