! ======================================================================
! Runs unit tests on lte_harmonic.
! ======================================================================
module lte_harmonic_test_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function lte_harmonic_test_keywords() result(keywords)
  use help_module
  use lte_harmonic_module
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  type(KeywordData), allocatable :: keywords_lte_harmonic(:)
  type(KeywordData)              :: keywords_test(0)
  
  integer :: ialloc
  
  keywords_lte_harmonic = lte_harmonic_keywords()
  
  allocate( keywords(size(keywords_lte_harmonic)+size(keywords_test)), &
          & stat=ialloc); call err(ialloc)
  
  keywords(:size(keywords_lte_harmonic)) = keywords_lte_harmonic
  keywords(size(keywords_lte_harmonic)+1:) = keywords_test
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine lte_harmonic_test(arguments)
  use constants_module, only : pi
  use utils_module,     only : mkdir
  use structure_module
  use group_module
  use atom_mapping_module
  use lte_module
  use unique_directions_module
  use displacement_patterns_module
  use lte_harmonic_module
  use qpoints_module
  use dictionary_module
  use linear_algebra_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Previous user inputs.
  type(Dictionary) :: setup_harmonic_arguments
  
  ! Directories and files.
  type(String) :: wd
  type(String) :: sdir
  type(String) :: qdir
  
  ! Setup data
  type(String), allocatable :: no_supercells_file(:)
  integer                   :: no_supercells
  type(String)              :: dft_code
  type(String)              :: seedname
  
  ! Structure data.
  type(StructureData)              :: structure
  type(StructureData), allocatable :: supercells(:)
  type(StructureData)              :: supercell
  type(UniqueDirections)           :: unique_directions
  
  ! Symmetry data.
  type(Group),      allocatable :: atom_symmetry_group(:)
  type(RealMatrix), allocatable :: rotations_cart(:)
  
  ! R-vector data.
  type(Group), allocatable  :: rvector_group(:)
  
  ! Force constant data.
  type(RealVector), allocatable :: forces(:,:)
  type(RealMatrix), allocatable :: force_constants(:,:,:)
  
  ! q-point data.
  type(StructureData)           :: structure_grid
  type(QpointData), allocatable :: qpoints_ibz(:)
  type(RealVector)              :: qpoint
  
  ! Dynamical matrix data.
  integer                          :: gvector
  type(ComplexMatrix), allocatable :: dynamical_matrix(:,:)
  real(dp)                         :: exponent
  
  ! Displacement pattern variables.
  type(LteReturn), allocatable :: lte_results(:)
  integer :: mode,atom,atom_sc
  
  ! Fourier interpolation data.
  integer  :: grid(3)
  real(dp) :: temperature
  
  ! Atom, mode and R-vector IDs.
  integer :: atom_1,atom_2
  integer :: atom_1p,atom_2p
  integer :: atom_1_sc,atom_2_sc
  integer :: atom_1p_sc,atom_2p_sc
  integer :: rvector,rvector_p,rvector_1p,rvector_2p
  
  ! Errors.
  real(dp) :: error, force
  type(RealMatrix) :: matrix, matrix_p
  type(RealVector) :: x, f, fx
  real(dp) :: total
  
  ! Temporary variables.
  integer                   :: i,j,k,l,ialloc
  
  ! --------------------------------------------------
  ! Read in arguments from user.
  ! --------------------------------------------------
  wd = arguments%value('working_directory')
  temperature = dble(arguments%value('temperature'))
  
  ! --------------------------------------------------
  ! Copy out settings to lte_harmonic.used_settings
  ! --------------------------------------------------
  call arguments%write_file(wd//'/lte_harmonic.used_settings')
  
  ! --------------------------------------------------
  ! Run new lte_harmonic.
  ! --------------------------------------------------
  call lte_harmonic(arguments)
  
  ! --------------------------------------------------
  ! Read in previous arguments.
  ! --------------------------------------------------
  call setup_harmonic_arguments%read_file(wd//'/setup_harmonic.used_settings')
  dft_code = setup_harmonic_arguments%value('dft_code')
  seedname = setup_harmonic_arguments%value('seedname')
  grid = int(split(setup_harmonic_arguments%value('q-point_grid')))
  
  ! Read in setup data.
  no_supercells_file = read_lines(wd//'/no_supercells.dat')
  no_supercells = int(no_supercells_file(1))
  
  ! Read in structure.
  structure = read_structure_file(wd//'/structure.dat')
  
  ! Read in supercells.
  allocate( supercells(no_supercells), &
          & lte_results(no_supercells), &
          & stat=ialloc); call err(ialloc)
  do i=1,no_supercells
    sdir = wd//'/Supercell_'//i
    supercells(i) = read_structure_file(sdir//'/structure.dat')
  enddo
  
  ! --------------------------------------------------
  ! Loop across supercells, testing each in turn.
  ! --------------------------------------------------
  
  do i=1,no_supercells
    call print_line('')
    call print_line('Checking supercell '//i)
    sdir = wd//'/Supercell_'//i
    
    supercell = supercells(i)
    
    rotations_cart = calculate_cartesian_rotations(supercell)
    
    ! Read in symmetry group and unique atoms.
    atom_symmetry_group = read_group_file(sdir//'/atom_symmetry_group.dat')
    unique_directions = read_unique_directions_file( &
       & sdir//'/unique_directions.dat')
    
    ! Read in forces.
    forces = read_forces(supercell, unique_directions, sdir, dft_code, &
       & seedname)
    
    ! --------------------------------------------------
    ! Check force constants.
    ! --------------------------------------------------
    ! Mass reduce forces.
    do j=1,size(unique_directions)
      do k=1,supercell%no_atoms
        forces(k,j) = forces(k,j)            &
                  & / sqrt ( supercell%mass(k) &
                  &        * supercell%mass(unique_directions%atoms(j)))
      enddo
    enddo
    
    ! Calculate force constants.
    force_constants = construct_force_constants(forces,supercell,&
       & unique_directions,atom_symmetry_group)
    
    ! Calculate L2 error between forces and force constants.
    error = 0.0_dp
    force = 0.0_dp
    do j=1,supercell%no_symmetries
      do k=1,size(unique_directions)
        
        if (unique_directions%directions_char(k)=='x') then
          x = [1.0_dp, 0.0_dp, 0.0_dp]
        elseif (unique_directions%directions_char(k)=='y') then
          x = [0.0_dp, 1.0_dp, 0.0_dp]
        else
          x = [0.0_dp, 0.0_dp, 1.0_dp]
        endif
        x = rotations_cart(j) * x
        
        do l=1,supercell%no_atoms
          atom_1 = atom_symmetry_group(j) * unique_directions%atoms(k)
          if (supercell%atom_to_rvec(atom_1)/=1) then
            cycle
          endif
          atom_2 = atom_symmetry_group(j) * l
          rvector = supercell%atom_to_rvec(atom_2)
          f = rotations_cart(j) * forces(l,k)
          fx = force_constants(rvector, supercelL%atom_to_prim(atom_2), supercell%atom_to_prim(atom_1)) * x
          error = error + (fx-f)*(fx-f)
          force = force + f*f
        enddo
      enddo
    enddo
    error = sqrt(error/force)
    call print_line( &
       & 'Fractional L2 difference between symmetrised and raw forces:')
    call print_line('New code: '//error)
    
    rvector_group = calculate_rvector_group(supercell)
    
    ! Calculate L2 error between symmetries.
    error = 0.0_dp
    total = 0.0_dp
    do j=1,supercell%no_symmetries
      do atom_1=1,structure%no_atoms
        do atom_2=1,structure%no_atoms
          do rvector=1,supercell%sc_size
            atom_1_sc = supercell%rvec_and_prim_to_atom(atom_1,1)
            atom_2_sc = supercell%rvec_and_prim_to_atom(atom_2,rvector)
            atom_1p_sc = atom_symmetry_group(j) * atom_1_sc
            atom_2p_sc = atom_symmetry_group(j) * atom_2_sc
            atom_1p = supercell%atom_to_prim(atom_1p_sc)
            atom_2p = supercell%atom_to_prim(atom_2p_sc)
            rvector_1p = supercell%atom_to_rvec(atom_1p_sc)
            rvector_2p = supercell%atom_to_rvec(atom_2p_sc)
            rvector_p = rvector_group(supercell%paired_rvec(rvector_1p)) &
                    & * rvector_2p
            matrix = force_constants(rvector,atom_2,atom_1)
            matrix_p = force_constants(rvector_p,atom_2p,atom_1p)
            matrix_p =  transpose(rotations_cart(j)) &
                       & * matrix_p                  &
                       & * rotations_cart(j)
            total = total + sum(dble(matrix)**2)
            error = error + sum(dble(matrix-matrix_p)**2)
          enddo
        enddo
      enddo
    enddo
    
    if (total < 1.0e-10_dp .and. error > 1.0e-10_dp) then
      call print_line('Error in symmetrisation of force constants.')
      call err()
    endif
    
    error = sqrt(error/total)
    
    if (total > 1.0e-10_dp .and. error > 1.0e-10_dp) then
      call print_line('Error in symmetrisation of force constants.')
      call err()
    endif
    
    call print_line( &
       & 'Fractional L2 error in symmetrisation of force constants:')
    call print_line('New code: '//error)
    
    ! Run lte on each supercell.
    lte_results(i) = evaluate_freqs_on_grid(supercell, force_constants)
    
    ! Check dynamical matrices.
    do gvector=1,supercell%sc_size
      qpoint = transpose(supercell%recip_supercell) &
           & * supercell%gvectors(gvector) / dble(supercell%sc_size)
      allocate( dynamical_matrix( supercell%no_atoms_prim,  &
              &                   supercell%no_atoms_prim), &
              & stat=ialloc); call err(ialloc)
      dynamical_matrix = mat([ &
   & cmplx(0.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp), &
   & cmplx(0.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp), &
   & cmplx(0.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp)],&
   & 3,3)
      do rvector=1,supercell%sc_size
        exponent = -2*pi*qpoint*supercell%rvectors(rvector)
        do j=1,supercell%no_atoms_prim
          do k=1,supercell%no_atoms_prim
            dynamical_matrix(k,j) = dynamical_matrix(k,j)        &
                                & + force_constants(rvector,k,j) &
                                & * cmplx(cos(exponent),sin(exponent),dp)
          enddo
        enddo
      enddo
      
      error = 0.0_dp
      total = 0.0_dp
      do j=1,supercell%no_atoms_prim
        do k=1,supercell%no_atoms_prim
          error = error+sum(abs(cmplx(dynamical_matrix(k,j)-lte_results(i)%dynamical_matrices(k,j,gvector)))**2)
          total = total+sum(abs(cmplx(dynamical_matrix(k,j)))**2)
        enddo
      enddo
      if (error > 1.0e-10_dp) then
        error = sqrt(error / total)
        if (error > 1.0e-10_dp) then
          call print_line('G-vector '//gvector//' dynamical matrix incorrect.')
          call err()
        endif
      endif
      call print_line('G-vector '//gvector//' dynamical matrix correct.')
      
      ! Check polarisation vectors and frequencies.
      dynamical_matrix = mat([ &
   & cmplx(0.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp), &
   & cmplx(0.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp), &
   & cmplx(0.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp)],&
   & 3,3)
      do mode=1,supercell%no_modes_prim
        do atom=1,supercell%no_atoms_prim
          do rvector=1,supercell%sc_size
            atom_sc = supercell%rvec_and_prim_to_atom(atom,rvector)
            ! TODO: fill this out.
          enddo
        enddo
      enddo
      
      deallocate(dynamical_matrix, stat=ialloc); call err(ialloc)
    enddo
  enddo
  
  call print_line('')
  call print_line('All force constants and dynamical matrices constructed.')
  
  ! --------------------------------------------------
  ! Loop across q-points, testing each in turn.
  ! --------------------------------------------------
  
  ! Read qpoint data.
  structure_grid = read_structure_file(wd//'/structure_grid.dat')
  qpoints_ibz = read_qpoints_file(wd//'/qpoints_ibz.dat')
  
  do i=1,size(qpoints_ibz)
    qdir = wd//'qpoint_'//i
    
    ! Read in frequencies.
    
    ! Read in displacements.
    
    ! Reconstruct real part of 
  enddo
end subroutine
end module
