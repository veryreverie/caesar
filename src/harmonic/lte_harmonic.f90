module lte_harmonic_module

contains

! Program to construct and execute LTE
subroutine lte_harmonic()
  use constants,      only : dp, eV_per_A_to_au
  use linear_algebra, only : invert
  use file_module
  use string_module
  use structure_module
  use dft_output_file_module
  use lte_module
  use fourier_interpolation_module
  use supercell_module
  use unique_directions_module
  use group_module
  implicit none
  
  ! User-input temperature
  real(dp) :: temperature
  
  ! Setup data
  integer             :: no_sc
  type(String)        :: dft_code
  type(String)        :: seedname
  type(StructureData) :: structure
  integer             :: grid(3)
  
  ! Supercell-specific setup data
  type(SupercellData), allocatable :: supercells(:)
  type(StructureData), allocatable :: structure_scs(:)
  
  ! Force constant data
  character(1)           :: direction
  type(Group)            :: symmetry_group
  type(UniqueDirections) :: unique_directions
  type(DftOutputFile)    :: dft
  real(dp),  allocatable :: force_constants(:,:,:,:)
  integer                :: xy,xz,yz
  integer                :: atom_1,atom_2,atom_1p,atom_2p
  real(dp)               :: force_consts_symm(3,3)
  real(dp)               :: transformation(3,3)
  logical                :: forces_calculated(3)
  logical,   allocatable :: atom_calculated(:)
  real(dp),  allocatable :: row_avgs(:,:,:)
  real(dp),  allocatable :: col_avgs(:,:,:)
  real(dp)               :: avg(3,3)
  integer                :: sc_size
  
  ! kpoint data
  integer               :: no_kpoints
  integer,  allocatable :: kpoints(:,:)
  integer,  allocatable :: multiplicity(:)
  integer,  allocatable :: sc_ids(:)
  integer, allocatable  :: gvector_ids(:)
  
  ! lte input data
  integer               :: no_kspace_lines
  real(dp), allocatable :: disp_kpoints(:,:)
  
  ! Temporary variables
  integer        :: i,j,k
  type(String)   :: sdir
  
  ! File contents
  type(String), allocatable :: user_inputs(:)
  
  ! File units
  integer :: no_sc_file
  integer :: ibz_file
  integer :: grid_file
  
  ! ----------------------------------------------------------------------
  ! Get temperature from user
  ! ----------------------------------------------------------------------
  write(*,"(a)") "What temperature (K)?"
  read(*,*) temperature
  
  ! ----------------------------------------------------------------------
  ! Read in initial data
  ! ----------------------------------------------------------------------
  user_inputs = read_lines('user_input.txt')
  dft_code = user_inputs(1)
  seedname = user_inputs(2)
  
  no_sc_file = open_read_file('no_sc.dat')
  read(no_sc_file,*) no_sc
  close(no_sc_file)
  
  structure = read_structure_file('structure.dat',identity_supercell())
  
  ! Read grid file
  grid_file = open_read_file('grid.dat')
  read(grid_file,*) grid
  close(grid_file)
  
  ! Read kpoints from ibz.dat
  no_kpoints = count_lines('ibz.dat')
  allocate(kpoints(3,no_kpoints))
  allocate(multiplicity(no_kpoints))
  allocate(sc_ids(no_kpoints))
  allocate(gvector_ids(no_kpoints))
  ibz_file = open_read_file('ibz.dat')
  do i=1,no_kpoints
    read(ibz_file,*) kpoints(:,i), multiplicity(i), sc_ids(i), gvector_ids(i)
  enddo
  
  ! Read in supercells
  supercells = read_supercells_file(str('supercells.dat'))
  
  ! Read in supercell structures
  allocate(structure_scs(no_sc))
  do i=1,no_sc
    sdir = str('Supercell_')//i
    structure_scs(i) = read_structure_file( sdir//'/structure.dat', &
                                          & supercells(i))
  enddo
  
  ! ----------------------------------------------------------------------
  ! Make directories
  ! ----------------------------------------------------------------------
  call system('mkdir lte')
  do i=1,no_sc
    sdir = str('Supercell_')//i
    call system('mkdir '//sdir//'/lte')
  enddo
  
  ! ----------------------------------------------------------------------
  ! Loop over supercells
  ! ----------------------------------------------------------------------
  do i=1,no_sc
    sdir = str('Supercell_')//i
    sc_size = structure_scs(i)%supercell%sc_size
    
    ! Read in symmetry group and unique atoms.
    symmetry_group = read_group_file(sdir//'/symmetry_group.dat')
    unique_directions = read_unique_directions_file( &
       & sdir//'/unique_directions.dat')
    
    ! Read forces from DFT output files.
    allocate(force_constants(3,3,structure_scs(i)%no_atoms, &
                                &structure_scs(i)%no_atoms))
    force_constants = 1e100_dp ! To make uninitialised values obvious.
    do j=1,size(unique_directions)
      atom_1 = unique_directions%unique_atoms(j)
      forces_calculated = (/ .true.,.true.,.true. /)
      if (unique_directions%xy_symmetry(j)/=0) then
        forces_calculated(2) = .false.
      endif
      if ( unique_directions%xz_symmetry(j)/=0 .or. &
         & unique_directions%yz_symmetry(j)/=0) then
        forces_calculated(3) = .false.
      endif
      
      do k=1,3
        if (.not. forces_calculated(k)) then
          cycle
        endif
        direction = directions(j)
        
        ! Calculate second derivatives of energy, using finite differences
        !    of force data.
        dft = read_dft_output_file(dft_code, &
           &                       sdir//'/atom.'//atom_1//'.+d'//direction, &
           &                       seedname)
        force_constants(k,:,:,atom_1) = dft%forces
        dft = read_dft_output_file(dft_code, &
           &                       sdir//'/atom.'//atom_1//'.-d'//direction, &
           &                       seedname)
        force_constants(k,:,:,atom_1) = force_constants(k,:,:,atom_1) &
                                    & - dft%forces
        force_constants(k,:,:,atom_1) = force_constants(k,:,:,atom_1) &
                                    & *eV_per_A_to_au / 0.02_dp
        
        ! Mass-reduce force constants.
        do atom_2=1,structure_scs(i)%no_atoms
          force_constants(:,:,atom_2,atom_1) =        &
             &   force_constants(:,:,atom_2,atom_1)   &
             & / dsqrt( structure_scs(i)%mass(atom_1) &
             &        * structure_scs(i)%mass(atom_2))
        enddo
      enddo
    enddo
    
    ! Reconstruct missing directions from symmetry.
    do j=1,size(unique_directions)
      atom_1 = unique_directions%unique_atoms(j)
      
      xy = unique_directions%xy_symmetry(j)
      xz = unique_directions%xz_symmetry(j)
      yz = unique_directions%yz_symmetry(j)
      
      do atom_2=1,structure_scs(i)%no_atoms
        ! Construct force constants in the space of independent vectors
        !    produced by the symmetry operations, and the transformation from
        !    cartesian co-ordinates to this basis.
        force_consts_symm(1,:) = force_constants(1,:,atom_2,atom_1)
        transformation(:,1) = (/ 1.0_dp,0.0_dp,0.0_dp /)
        
        if (xy==0) then
          force_consts_symm(2,:) = force_constants(2,:,atom_2,atom_1)
          transformation(:,2) = (/ 0.0_dp,1.0_dp,0.0_dp /)
        else
          force_consts_symm(2,:) = force_consts_symm(1,:)
          transformation(:,2) = structure_scs(i)%rotation_matrices(:,1,xy)
        endif
        
        if (xz==0 .and. yz==0) then
          force_consts_symm(3,:) = force_constants(3,:,atom_2,atom_1)
          transformation(:,3) = (/ 0.0_dp,0.0_dp,1.0_dp /)
        elseif (xz/=0) then
          force_consts_symm(3,:) = force_consts_symm(1,:)
          transformation(:,3) = structure_scs(i)%rotation_matrices(:,1,xz)
        else
          force_consts_symm(3,:) = force_consts_symm(2,:)
          transformation(:,3) = structure_scs(i)%rotation_matrices(:,2,yz)
        endif
        
        ! Transform force constants into cartesian co-ordinates.
        force_constants(:,:,atom_2,atom_1) = matmul( force_consts_symm, &
                                                   & invert(transformation))
      enddo
    enddo
    
    ! Copy force constants to all positions related by symmetry.
    atom_calculated = .false.
    do j=1,size(unique_directions)
      atom_1 = unique_directions%unique_atoms(j)
      atom_calculated(atom_1) = .true.
      do k=1,size(symmetry_group)
        atom_1p = operate(atom_1,symmetry_group,k)
        if (atom_calculated(atom_1)) then
          cycle
        endif
        do atom_2=1,structure_scs(i)%no_atoms
          atom_2p = operate(atom_2,symmetry_group,k)
          force_constants(:,:,atom_2p,atom_1p) = matmul(matmul( &
             & structure_scs(i)%rotation_matrices(:,:,k),       &
             & force_constants(:,:,atom_2,atom_1)),             &
             & transpose(structure_scs(i)%rotation_matrices(:,:,k)))
        enddo
      enddo
    enddo
    
    ! Impose order-of-differentiation symmetry.
    do atom_1=1,structure_scs(i)%no_atoms
      do atom_2=1,structure_scs(i)%no_atoms
        force_constants(:,:,atom_2,atom_1) = (               &
           & force_constants(:,:,atom_2,atom_1)              &
           & + transpose(force_constants(:,:,atom_1,atom_2)) &
           & ) / 2
        force_constants(:,:,atom_1,atom_2) = &
           & transpose(force_constants(:,:,atom_2,atom_1))
      enddo
    enddo
    
    ! Impose Newton III symmetry.
    allocate(row_avgs(3,3,structure_scs(i)%no_atoms))
    allocate(col_avgs(3,3,structure_scs(i)%no_atoms))
    row_avgs = sum(force_constants,3)/structure_scs(i)%no_atoms
    col_avgs = sum(force_constants,4)/structure_scs(i)%no_atoms
    avg = sum(row_avgs,3)
    do atom_1=1,structure_scs(i)%no_atoms
      do atom_2=1,structure_scs(i)%no_atoms
        force_constants(:,:,atom_2,atom_1) =    &
           & force_constants(:,:,atom_2,atom_1) &
           & - row_avgs(:,:,atom_1)             &
           & - col_avgs(:,:,atom_2)             &
           & + avg
      enddo
    enddo
    deallocate(row_avgs)
    deallocate(col_avgs)
        
    call lte_4( structure,                         &
              & structure_scs(i),                  &
              & force_constants,                   &
              & 0.0_dp,                            &
              & sdir//'/lte/kpairs.dat',           &
              & sdir//'/lte/freq_grids.dat',       &
              & sdir//'/lte/disp_patterns.dat',    &
              & sdir//'/lte/kdisp_patterns.dat',   &
              & sdir//'/lte/pol_vec.dat',          &
              & sdir//'/lte/error.txt',            &
              & sdir//'/lte/dyn_mat.')
    
    if (file_exists(sdir//'/lte/error.txt')) then
      write(*,"(a)") "There is an error in lte: check error.txt file."
      stop
    endif
    
    deallocate(force_constants)
  enddo
  
  ! Write path for fourier interpolation
  no_kspace_lines = 4
  allocate(disp_kpoints(3,0:no_kspace_lines))
  disp_kpoints(:,0) = (/ 0.0_dp, 0.0_dp, 0.0_dp /) ! GM
  disp_kpoints(:,1) = (/ 0.5_dp, 0.5_dp, 0.5_dp /) ! T
  disp_kpoints(:,2) = (/ 0.0_dp, 0.5_dp, 0.5_dp /) ! FB
  disp_kpoints(:,3) = (/ 0.0_dp, 0.0_dp, 0.0_dp /) ! GM
  disp_kpoints(:,4) = (/ 0.0_dp, 0.5_dp, 0.0_dp /) ! L
    
  call fourier_interpolation(                  &
     & structure,                              &
     & grid,                                   &
     & temperature,                            &
     & kpoints, sc_ids, gvector_ids,           &
     & str('lte/atoms_in_primitive_cell.'),    &
     & str('lte/dyn_mat.'),                    &! Supercell_*/lte/dyn_mat.*.dat
     & disp_kpoints,                           &
     & str('lte/phonon_dispersion_curve.dat'), &
     & str('lte/high_symmetry_points.dat'),    &
     & str('lte/free_energy.dat'),             &
     & str('lte/freq_dos.dat'))
  
end subroutine
end module
