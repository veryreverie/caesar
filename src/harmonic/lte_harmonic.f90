module lte_harmonic_module

contains

! ----------------------------------------------------------------------
! Reads dft output files in sdir, and generates matrix of force constants.
! ----------------------------------------------------------------------
function calculate_force_constants(structure,structure_sc,symmetry_group, &
   & unique_directions,sdir,dft_code,seedname) result(output)
  use constants,      only : dp, eV_per_A_to_au, directions
  use utils,          only : mkdir
  use linear_algebra, only : invert
  use file_module
  use string_module
  use structure_module
  use dft_output_file_module
  use lte_module
  use supercell_module
  use unique_directions_module
  use group_module
  implicit none
  
  type(StructureData),    intent(in) :: structure
  type(StructureData),    intent(in) :: structure_sc
  type(Group),            intent(in) :: symmetry_group(:)
  type(UniqueDirections), intent(in) :: unique_directions
  type(String),           intent(in) :: sdir
  type(String),           intent(in) :: dft_code
  type(String),           intent(in) :: seedname
  real(dp), allocatable :: output(:,:,:)
  
  ! Setup data
  
  ! Force constant data
  character(1)                       :: direction
  type(DftOutputFile)                :: positive
  type(DftOutputFile)                :: negative
  real(dp),              allocatable :: force_constants(:,:,:,:)
  integer                            :: xy,xz,yz
  logical                            :: forces_calculated(3)
  logical,               allocatable :: atom_calculated(:)
  real(dp),              allocatable :: totals(:,:,:)
  real(dp)                           :: total(3,3)
  
  ! Atom ids in various representations.
  !    _sc is in supercell ordering.
  !    _prim is in primitive cell ordering.
  integer                            :: atom_1_prim,atom_2_prim
  integer                            :: atom_1_sc,atom_2_sc
  integer                            :: atom_1p_prim
  integer                            :: atom_1p_sc,atom_2p_sc
  integer                            :: mode_1,mode_2
  
  ! Force constant construction by symmetry.
  integer  :: a
  integer  :: b
  integer  :: c
  integer  :: d
  real(dp) :: Fab(3,3)
  real(dp) :: Fac(3,3)
  real(dp) :: Fad(3,3)
  real(dp) :: Rcb(3,3)
  real(dp) :: Rdb(3,3)
  
  ! Temporary variables
  integer        :: j,k
  
  allocate(atom_calculated(structure%no_atoms))
  atom_calculated = .false.
  
  ! Read forces from DFT output files.
  allocate(force_constants(3,3,structure_sc%no_atoms,structure%no_atoms))
  
  force_constants = 0.0_dp
  do j=1,size(unique_directions)
    atom_1_sc = unique_directions%unique_atoms(j)
    atom_1_prim = structure_sc%atom_to_prim(atom_1_sc)
    atom_calculated(atom_1_prim) = .true.
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
      
      direction = directions(k)
      
      ! Calculate second derivatives of energy, using finite differences
      !    of force data.
      positive = read_dft_output_file(dft_code,      &
         & sdir//'/atom.'//atom_1_sc//'.+d'//direction, &
         & seedname)
      negative = read_dft_output_file(dft_code,      &
         & sdir//'/atom.'//atom_1_sc//'.-d'//direction, &
         & seedname)
      
      force_constants(k,:,:,atom_1_prim) = (positive%forces-negative%forces) &
                                       & * eV_per_A_to_au / 0.02_dp
    enddo
  enddo
  
  ! Reconstruct missing directions from symmetry.
  !
  ! A symmetry operation takes atom a to atom a, and atom b to atom c
  ! The rotation part of this symmetry is Rcb.
  !
  ! Rotation matrix = R = sum_i(|R(i,:)><i|)
  ! Force constants = F = sum_i(|F(i,:)><i|)
  !
  ! => Fac.Rcb = Rcb.Fab
  ! => <Fac(i,:)|Rcb = <Rcb(i,:)|Fab
  do j=1,size(unique_directions)
    atom_1_sc = unique_directions%unique_atoms(j)
    a = structure_sc%atom_to_prim(atom_1_sc)
    
    xy = unique_directions%xy_symmetry(j)
    xz = unique_directions%xz_symmetry(j)
    yz = unique_directions%yz_symmetry(j)
    
    do b=1,structure_sc%no_atoms
      Fab = force_constants(:,:,b,a)
      if (xy/=0 .and. xz==0) then
        ! Construct Fab(2,:).
        
        ! <Fac(1,:)|Rcb = <Rcb(1,:)|Fab
        !               = Rcb(1,1)Fab(1,:)+Rcb(1,2)Fab(2,:)+Rcb(1,3)Fab(3,:)
        ! => Fab(2,:) = [Fac(1,:)Rcb-Rcb(1,:)Fab'] / Rcb(1,2)
        ! Where Fab'(1,:)=Fab(1,:), Fab'(3,:)=Fab(3,:) and Fab'(2,:)=0
        c = operate(symmetry_group(xy),b)
        Rcb = structure_sc%rotation_matrices(:,:,xy)
        Fac = force_constants(:,:,c,a)
        Fab(2,:) = (matmul(Fac(1,:),Rcb)-matmul(Rcb(1,:),Fab)) / Rcb(1,2)
      elseif (xy==0 .and. xz/=0) then
        ! Construct Fab(3,:).
        
        ! As above, but with Fab'(3,:)=0 rather than Fab'(2,:)=0
        ! => Fab(3,:) = [Fac(1,:)Rcb-Rcb(1,:)F'] / R(1,3)
        c = operate(symmetry_group(xz),b)
        Rcb = structure_sc%rotation_matrices(:,:,xz)
        Fac = force_constants(:,:,c,a)
        Fab(3,:) = (matmul(Fac(1,:),Rcb)-matmul(Rcb(1,:),Fab)) / Rcb(1,3)
      elseif (yz/=0) then
        ! Construct Fab(3,:).
        
        ! Similarly,
        ! => Fab(3,:) = [Fac(2,:)Rcb-Rcb(2,:)Fab'] / Rcb(2,3)
        c = operate(symmetry_group(yz),b)
        Rcb = structure_sc%rotation_matrices(:,:,yz)
        Fac = force_constants(:,:,c,a)
        Fab(3,:) = (matmul(Fac(2,:),Rcb)-matmul(Rcb(2,:),Fab)) / Rcb(2,3)
      elseif (xy/=0 .and. xz/=0) then
        ! Construct Fab(2,:).
        
        ! Fac(1,:)Rcb = Rcb(1,1)Fab(1,:)+Rcb(1,2)Fab(2,:)+Rdb(1,3)Fab(3,:)
        ! Fad(1,:)Rdb = Rdb(1,1)Fab(1,:)+Rdb(1,2)Fab(2,:)+Rdb(1,3)Fab(3,:)
        ! =>
        ! /Rcb(1,2) Rcb(1,3)\ /Fab(2,:)\ = /Fac(1,:)Rcb-Rcb(1,1)Fab(1,:)\ 
        ! \Rdb(1,2) Rdb(1,3)/ \Fab(3,:)/   \Fad(1,:)Rdb-Rdb(1,1)Fab(1,:)/
        ! =>
        ! /Fab(2,:)\ = / Rdb(1,3) -Rcb(1,3)\ /Fac(1,:)Rcb-Rcb(1,1)Fab(1,:)\ 
        ! \Fab(3,:)/   \-Rdb(1,2)  Rcb(1,2)/ \Fad(1,:)Rdb-Rdb(1,1)Fab(1,:)/
        !              ----------------------------------------------------
        !                        Rcb(1,2)Rdb(1,3)-Rcb(1,3)Rdb(1,2)
        c = operate(symmetry_group(xy),b)
        d = operate(symmetry_group(xz),b)
        Rcb = structure_sc%rotation_matrices(:,:,xy)
        Rdb = structure_sc%rotation_matrices(:,:,xz)
        Fac = force_constants(:,:,c,a)
        Fad = force_constants(:,:,d,a)
        Fab(2,:) = ( Rdb(1,3)*(matmul(Fac(1,:),Rcb)-Rcb(1,1)*Fab(1,:))  &
               &   - Rcb(1,3)*(matmul(Fad(2,:),Rdb)-Rdb(1,1)*Fab(1,:))) &
               & / (Rcb(1,2)*Rdb(1,3)-Rcb(1,3)*Rdb(1,2))
        Fab(2,:) = (-Rdb(1,2)*(matmul(Fac(1,:),Rcb)-Rcb(1,1)*Fab(1,:))  &
               &   + Rcb(1,2)*(matmul(Fad(2,:),Rdb)-Rdb(1,1)*Fab(1,:))) &
               & / (Rcb(1,2)*Rdb(1,3)-Rcb(1,3)*Rdb(1,2))
      endif
      force_constants(:,:,b,a) = Fab
    enddo
  enddo
  
  ! Average across symmetries.
  do j=1,size(unique_directions)
    atom_1_sc = unique_directions%unique_atoms(j)
    atom_1_prim = structure_sc%atom_to_prim(atom_1_sc)
    do k=1,size(symmetry_group)
      atom_1p_sc = operate(symmetry_group(k),atom_1_sc)
      atom_1p_prim = structure_sc%atom_to_prim(atom_1p_sc)
      
      if (structure_sc%atom_to_gvec(atom_1p_sc)/=1) then
        cycle
      endif
      
      if (.not. atom_calculated(atom_1p_prim)) then
        cycle
      endif
      
      do atom_2_sc=1,structure_sc%no_atoms
        atom_2p_sc = operate(symmetry_group(k),atom_2_sc)
        
        force_constants(:,:,atom_2p_sc,atom_1p_prim) = (         &
           &   matmul(matmul(                                    &
           &   structure_sc%rotation_matrices(:,:,k),            &
           &   force_constants(:,:,atom_2_sc,atom_1_prim)),      &
           &   transpose(structure_sc%rotation_matrices(:,:,k))) &
           & +                                                   &
           &   force_constants(:,:,atom_2p_sc,atom_1p_prim)      &
           & ) / 2
        
        force_constants(:,:,atom_2_sc,atom_1_prim) =             &
           &   matmul(matmul(                                    &
           &   transpose(structure_sc%rotation_matrices(:,:,k)), &
           &   force_constants(:,:,atom_2p_sc,atom_1p_prim)),    &
           &   structure_sc%rotation_matrices(:,:,k))
      enddo
    enddo
  enddo
  
  ! Copy force constants to all positions related by symmetry.
  do j=1,size(unique_directions)
    atom_1_sc = unique_directions%unique_atoms(j)
    atom_1_prim = structure_sc%atom_to_prim(atom_1_sc)
    do k=1,size(symmetry_group)
      atom_1p_sc = operate(symmetry_group(k),atom_1_sc)
      atom_1p_prim = structure_sc%atom_to_prim(atom_1p_sc)
      
      if (structure_sc%atom_to_gvec(atom_1p_sc)/=1) then
        cycle
      endif
      
      if (atom_calculated(atom_1p_prim)) then
        cycle
      endif
      
      atom_calculated(atom_1p_prim) = .true.
      
      do atom_2_sc=1,structure_sc%no_atoms
        atom_2p_sc = operate(symmetry_group(k),atom_2_sc)
        force_constants(:,:,atom_2p_sc,atom_1p_prim) = matmul(matmul( &
           & structure_sc%rotation_matrices(:,:,k),                   &
           & force_constants(:,:,atom_2_sc,atom_1_prim)),             &
           & transpose(structure_sc%rotation_matrices(:,:,k)))
      enddo
    enddo
  enddo
  
  deallocate(atom_calculated)
  
  ! Impose order-of-differentiation symmetry.
  do atom_1_prim=1,structure%no_atoms
    do atom_2_prim=1,atom_1_prim
      do j=1,structure_sc%sc_size
        k = structure_sc%paired_gvec(j)
        atom_1_sc = structure_sc%gvec_and_prim_to_atom(atom_1_prim,j)
        atom_2_sc = structure_sc%gvec_and_prim_to_atom(atom_2_prim,k)
        
        force_constants(:,:,atom_2_sc,atom_1_prim) = (               &
           &   force_constants(:,:,atom_2_sc,atom_1_prim)            &
           & + transpose(force_constants(:,:,atom_1_sc,atom_2_prim)) &
           & ) / 2
        
        force_constants(:,:,atom_1_sc,atom_2_prim) = &
           & transpose(force_constants(:,:,atom_2_sc,atom_1_prim))
      enddo
    enddo
  enddo
  
  ! Impose translational symmetry.
  allocate(totals(3,3,structure%no_atoms))
  totals = 0.0_dp
  do j=1,structure_sc%sc_size
    do atom_1_prim=1,structure%no_atoms
      atom_1_sc = structure_sc%gvec_and_prim_to_atom(atom_1_prim,1)
      do atom_2_prim=1,structure%no_atoms
        atom_2_sc = structure_sc%gvec_and_prim_to_atom(atom_2_prim,j)
        
        if (atom_1_sc==atom_2_sc) then
          cycle
        endif
        
        totals(:,:,atom_2_prim) = totals(:,:,atom_2_prim) &
                              & + force_constants(:,:,atom_2_sc,atom_1_prim)
      enddo
    enddo
  enddo
  total = sum(totals,3)
  
  do j=1,structure_sc%sc_size
    do atom_1_prim=1,structure%no_atoms
      atom_1_sc = structure_sc%gvec_and_prim_to_atom(atom_1_prim,1)
      do atom_2_prim=1,structure%no_atoms
        atom_2_sc = structure_sc%gvec_and_prim_to_atom(atom_2_prim,j)
        
        if (atom_1_sc==atom_2_sc) then
          cycle
        endif
        
!        force_constants(:,:,atom_2_sc,atom_1_prim) =                   &
!           & force_constants(:,:,atom_2_sc,atom_1_prim)                &
!           & - (totals(:,:,atom_1_prim)+totals(:,:,atom_2_prim)-total) &
!           & / (structure_sc%no_atoms-1)
      enddo
    enddo
  enddo
  deallocate(totals)
  
  ! Mass-reduce force constants.
  do j=1,structure_sc%sc_size
    do atom_1_prim=1,structure%no_atoms
      do atom_2_prim=1,structure%no_atoms
        atom_2_sc = structure_sc%gvec_and_prim_to_atom(atom_2_prim,j)
        force_constants(:,:,atom_2_sc,atom_1_prim) =        &
           &   force_constants(:,:,atom_2_sc,atom_1_prim)   &
           & / dsqrt(structure%mass(atom_1_prim)*structure_sc%mass(atom_2_sc))
      enddo
    enddo
  enddo
  
  ! Convert into mode-mode-Gvector representation.
  allocate(output(structure%no_modes,structure%no_modes,structure_sc%sc_size))
  do atom_1_prim=1,structure%no_atoms
    mode_1 = (atom_1_prim-1)*3+1
    do atom_2_prim=1,structure%no_atoms
      mode_2 = (atom_2_prim-1)*3+1
      do j=1,structure_sc%sc_size
        atom_2_sc = structure_sc%gvec_and_prim_to_atom(atom_2_prim,j)
        output(mode_2:mode_2+2,mode_1:mode_1+2,j) = &
           & transpose(force_constants(:,:,atom_2_sc,atom_1_prim))
      enddo
    enddo
  enddo
end function

! Program to construct and execute LTE
subroutine lte_harmonic()
  use constants,      only : dp
  use utils,          only : mkdir
  use linear_algebra, only : invert
  use file_module
  use string_module
  use structure_module
  use dft_output_file_module
  use lte_module
  use supercell_module
  use unique_directions_module
  use group_module
  use err_module
  implicit none
  
  ! User-input temperature
  real(dp) :: temperature
  
  ! File contents
  type(String), allocatable :: user_inputs(:)
  type(String), allocatable :: grid_file(:)
  type(String), allocatable :: no_sc_file(:)
  type(String), allocatable :: ibz_file(:)
  
  ! Setup data
  integer             :: no_sc
  type(String)        :: dft_code
  type(String)        :: seedname
  type(StructureData) :: structure
  integer             :: grid(3)
  
  ! Supercell-specific setup data
  type(StructureData) :: structure_sc
  
  ! Force constant data
  type(Group),           allocatable :: symmetry_group(:)
  type(UniqueDirections)             :: unique_directions
  real(dp),              allocatable :: force_constants(:,:,:)
  
  ! kpoint data
  integer              :: no_kpoints
  integer, allocatable :: kpoints(:,:)
  integer, allocatable :: multiplicity(:)
  integer, allocatable :: sc_ids(:)
  integer, allocatable :: gvector_ids(:)
  
  ! lte input data
  integer               :: no_kspace_lines
  real(dp), allocatable :: disp_kpoints(:,:)
  
  ! Lte output data
  type(LteReturn)          :: lte_result
  complex(dp), allocatable :: ibz_dynamical_matrices(:,:,:)
  integer                  :: frequencies_file
  integer                  :: prefactors_file
  integer                  :: displacement_pattern_file
  integer                  :: mode
  integer                  :: atom
  
  ! Temporary variables
  integer                   :: i,j
  type(String)              :: sdir
  type(String), allocatable :: line(:)
  
  ! ----------------------------------------------------------------------
  ! Get temperature from user
  ! ----------------------------------------------------------------------
  call print_line('What temperature (K)?')
  temperature = dble(read_line_from_user())
  
  ! ----------------------------------------------------------------------
  ! Read in initial data
  ! ----------------------------------------------------------------------
  user_inputs = read_lines('user_input.txt')
  dft_code = user_inputs(1)
  seedname = user_inputs(2)
  
  no_sc_file = read_lines('no_sc.dat')
  no_sc = int(no_sc_file(1))
  
  structure = read_structure_file('structure.dat')
  
  ! Read grid file
  grid_file = read_lines('grid.dat')
  grid = int(split(grid_file(1)))
  
  ! Read kpoints from ibz.dat
  ibz_file = read_lines('ibz.dat')
  no_kpoints = size(ibz_file)
  allocate(kpoints(3,no_kpoints))
  allocate(multiplicity(no_kpoints))
  allocate(sc_ids(no_kpoints))
  allocate(gvector_ids(no_kpoints))
  do i=1,no_kpoints
    line = split(ibz_file(i))
    kpoints(:,i) = int(line(1:3))
    multiplicity(i) = int(line(4))
    sc_ids(i) = int(line(5))
    gvector_ids(i) = int(line(6))
  enddo
  
  ! ----------------------------------------------------------------------
  ! Make directories
  ! ----------------------------------------------------------------------
  allocate(ibz_dynamical_matrices( structure%no_modes, &
                                 & structure%no_modes, &
                                 & no_kpoints))
  
  ! ----------------------------------------------------------------------
  ! Loop over supercells
  ! ----------------------------------------------------------------------
  do i=1,no_sc
    sdir = 'Supercell_'//i
    
    ! Read in supercell structure data.
    structure_sc = read_structure_file(sdir//'/structure.dat')
    
    ! Read in symmetry group and unique atoms.
    symmetry_group = read_group_file(sdir//'/symmetry_group.dat')
    unique_directions = read_unique_directions_file( &
       & sdir//'/unique_directions.dat')
    
    ! Calculate force constants.
    force_constants = calculate_force_constants(structure,structure_sc, &
       & symmetry_group,unique_directions,sdir,dft_code,seedname)
    
    lte_result = evaluate_freqs_on_grid( structure,                       &
                                       & structure_sc,                    &
                                       & force_constants)
    
    deallocate(force_constants)
    
    do j=1,no_kpoints
      if (sc_ids(j)/=i) then
        cycle
      endif
      
      ! Move dynamical matrices into ibz_dynamical matrices.
      ibz_dynamical_matrices(:,:,j) = &
         & lte_result%dynamical_matrices(:,:,gvector_ids(j))
    
      ! Write out frequencies.
      frequencies_file = open_write_file( &
         & 'frequencies.kpoint_'//j//'.dat')
      do mode=1,structure%no_modes
        call print_line(frequencies_file, &
           & lte_result%frequencies(mode,gvector_ids(j)))
      enddo
      
      ! Write out prefactors.
      prefactors_file = open_write_file( &
         & 'prefactors.kpoint_'//j//'.dat')
      do mode=1,structure%no_modes
        call print_line(prefactors_file,'Mode : '//mode)
        do atom=1,structure_sc%no_atoms
          call print_line(prefactors_file, &
             & lte_result%prefactors(atom,mode,gvector_ids(j)))
        enddo
        call print_line(prefactors_file,'')
      enddo
      
      ! Write out displacement patterns.
      displacement_pattern_file = open_write_file( &
         & 'displacements.kpoint_'//j//'.dat')
      do mode=1,structure%no_modes
        call print_line(displacement_pattern_file,'Mode : '//mode)
        do atom=1,structure_sc%no_atoms
          call print_line(displacement_pattern_file, &
             & lte_result%displacements(:,atom,mode,gvector_ids(j)))
        enddo
        call print_line(displacement_pattern_file,'')
      enddo
      close(displacement_pattern_file)
    enddo
  enddo
  
  ! Write path for fourier interpolation
  no_kspace_lines = 4
  allocate(disp_kpoints(3,no_kspace_lines+1))
  disp_kpoints(:,1) = (/ 0.0_dp, 0.0_dp, 0.0_dp /) ! GM
  disp_kpoints(:,2) = (/ 0.5_dp, 0.5_dp, 0.5_dp /) ! T
  disp_kpoints(:,3) = (/ 0.0_dp, 0.5_dp, 0.5_dp /) ! FB
  disp_kpoints(:,4) = (/ 0.0_dp, 0.0_dp, 0.0_dp /) ! GM
  disp_kpoints(:,5) = (/ 0.0_dp, 0.5_dp, 0.0_dp /) ! L
  
  ! Read in primitive symmetry group.
  symmetry_group = read_group_file(str('Supercell_1/symmetry_group.dat'))
  
!  call print_line('Running fourier interpolation (this may take some time).')
!  call fourier_interpolation(              &
!     & ibz_dynamical_matrices,             &
!     & structure,                          &
!     & grid,                               &
!     & temperature,                        &
!     & kpoints,                            &
!     & disp_kpoints,                       &
!     & symmetry_group,                     &
!     & str('phonon_dispersion_curve.dat'), &
!     & str('high_symmetry_points.dat'),    &
!     & str('free_energy.dat'),             &
!     & str('freq_dos.dat'))
end subroutine
end module
