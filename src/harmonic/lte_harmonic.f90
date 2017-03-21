module lte_harmonic_module

contains

! Program to construct and execute LTE
subroutine lte_harmonic()
  use constants,      only : dp, eV_per_A_to_au
  use utils,          only : mkdir
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
  type(StructureData) :: structure_sc
  
  ! Force constant data
  character(1)                       :: direction
  type(Group),           allocatable :: symmetry_group(:)
  type(UniqueDirections)             :: unique_directions
  type(DftOutputFile)                :: positive
  type(DftOutputFile)                :: negative
  real(dp),              allocatable :: force_constants(:,:,:,:)
  real(dp),              allocatable :: force_constants_2(:,:,:)
  integer                            :: xy,xz,yz
  integer                            :: atom_1,atom_2,atom_1p,atom_2p
  integer                            :: mode_1,mode_2
  logical                            :: forces_calculated(3)
  logical,               allocatable :: atom_calculated(:)
  real(dp),              allocatable :: atom_1_sums(:,:,:)
  real(dp),              allocatable :: atom_2_sums(:,:,:)
  real(dp)                           :: sums(3,3)
  
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
  
  ! kpoint data
  integer              :: no_kpoints
  integer, allocatable :: kpoints(:,:)
  integer, allocatable :: multiplicity(:)
  integer, allocatable :: sc_ids(:)
  integer, allocatable :: gvector_ids(:)
  
  ! lte input data
  integer               :: no_kspace_lines
  real(dp), allocatable :: disp_kpoints(:,:)
  
  ! Temporary variables
  integer        :: i,j,k
  type(String)   :: sdir
  
  ! File contents
  type(String), allocatable :: user_inputs(:)
  type(String), allocatable :: grid_file(:)
  
  ! File units
  integer :: no_sc_file
  integer :: ibz_file
  integer :: force_constants_file
  
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
  
  no_sc_file = open_read_file('no_sc.dat')
  read(no_sc_file,*) no_sc
  close(no_sc_file)
  
  structure = read_structure_file('structure.dat')
  
  ! Read grid file
  grid_file = read_lines('grid.dat')
  grid = int(split(grid_file(1)))
  
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
  
  ! ----------------------------------------------------------------------
  ! Make directories
  ! ----------------------------------------------------------------------
  call mkdir('lte')
  do i=1,no_sc
    sdir = 'Supercell_'//i
    call mkdir(sdir//'/lte')
  enddo
  
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
    
    allocate(atom_calculated(structure%no_atoms))
    atom_calculated = .false.
    
    ! Read forces from DFT output files.
    allocate(force_constants(3,3,structure_sc%no_atoms,structure%no_atoms))
    force_constants = 0.0_dp
    do j=1,size(unique_directions)
      atom_1p = unique_directions%unique_atoms(j)
      atom_1 = structure_sc%atom_to_prim(atom_1p)
      atom_calculated(atom_1) = .true.
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
           & sdir//'/atom.'//atom_1p//'.+d'//direction, &
           & seedname)
        negative = read_dft_output_file(dft_code,      &
           & sdir//'/atom.'//atom_1p//'.-d'//direction, &
           & seedname)
        
        force_constants(k,:,:,atom_1) = (positive%forces-negative%forces) &
                                    & * eV_per_A_to_au / 0.02_dp
      enddo
    enddo
    
    ! Reconstruct missing directions from symmetry.
    
    ! A symmetry operation takes atom a to atom a, and atom b to atom c
    ! The rotation part of this symmetry is Rcb.
    
    ! Rotation matrix = R = sum_i(|R(i,:)><i|)
    ! Force constants = F = sum_i(|F(i,:)><i|)
    
    ! => Fac.Rcb = Rcb.Fab
    ! => <Fac(i,:)|Rcb = <Rcb(i,:)|Fab
    do j=1,size(unique_directions)
      atom_1p = unique_directions%unique_atoms(j)
      a = structure_sc%atom_to_prim(atom_1p)
      
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
      atom_1 = unique_directions%unique_atoms(j)
      do k=1,size(symmetry_group)
        a = operate(symmetry_group(k),atom_1)
        atom_1p = structure_sc%atom_to_prim(a)
        
        if (structure_sc%atom_to_gvec(a)/=1) then
          cycle
        endif
        
        if (.not. atom_calculated(atom_1p)) then
          cycle
        endif
        
        do atom_2=1,structure_sc%no_atoms
          atom_2p = operate(symmetry_group(k),atom_2)
          
          force_constants(:,:,atom_2p,atom_1p) = (                 &
             &   matmul(matmul(                                    &
             &   structure_sc%rotation_matrices(:,:,k),            &
             &   force_constants(:,:,atom_2,atom_1)),              &
             &   transpose(structure_sc%rotation_matrices(:,:,k))) &
             & +                                                   &
             &   force_constants(:,:,atom_2p,atom_1p)              &
             & ) / 2
          
          force_constants(:,:,atom_2,atom_1) =                     &
             &   matmul(matmul(                                    &
             &   transpose(structure_sc%rotation_matrices(:,:,k)), &
             &   force_constants(:,:,atom_2,atom_1)),              &
             &   structure_sc%rotation_matrices(:,:,k))
        enddo
      enddo
    enddo
    
    ! Copy force constants to all positions related by symmetry.
    do j=1,size(unique_directions)
      atom_1 = unique_directions%unique_atoms(j)
      do k=1,size(symmetry_group)
        a = operate(symmetry_group(k),atom_1)
        atom_1p = structure_sc%atom_to_prim(a)
        
        if (structure_sc%atom_to_gvec(a)/=1) then
          cycle
        endif
        
        if (atom_calculated(atom_1p)) then
          cycle
        endif
        
        atom_calculated(atom_1p) = .true.
        
        do atom_2=1,structure_sc%no_atoms
          atom_2p = operate(symmetry_group(k),atom_2)
          force_constants(:,:,atom_2p,atom_1p) = matmul(matmul( &
             & structure_sc%rotation_matrices(:,:,k),           &
             & force_constants(:,:,atom_2,atom_1)),             &
             & transpose(structure_sc%rotation_matrices(:,:,k)))
        enddo
      enddo
    enddo
    
    deallocate(atom_calculated)
    
    ! Impose order-of-differentiation symmetry.
    do atom_1=1,structure%no_atoms
      do atom_2=1,structure%no_atoms
        do j=1,structure_sc%sc_size
          k = structure_sc%paired_gvec(j)
          atom_1p = structure_sc%gvec_and_prim_to_atom(atom_1,j)
          atom_2p = structure_sc%gvec_and_prim_to_atom(atom_2,k)
          
          force_constants(:,:,atom_2p,atom_1) = (               &
             &   force_constants(:,:,atom_2p,atom_1)            &
             & + transpose(force_constants(:,:,atom_1p,atom_2)) &
             & ) / 2
          
          force_constants(:,:,atom_1p,atom_1) = &
             & transpose(force_constants(:,:,atom_2p,atom_1))
        enddo
      enddo
    enddo
    
    ! Impose translational symmetry.
    allocate(atom_1_sums(3,3,structure%no_atoms))
    allocate(atom_2_sums(3,3,structure%no_atoms))
    do j=1,structure_sc%sc_size
      atom_1_sums = 0.0_dp
      atom_2_sums = 0.0_dp
      do atom_1=1,structure%no_atoms
        do atom_2=1,structure%no_atoms
          atom_2p = structure_sc%gvec_and_prim_to_atom(atom_2,j)
          atom_1_sums(:,:,atom_2) = atom_1_sums(:,:,atom_2) &
                                & + force_constants(:,:,atom_2p,atom_1)
          atom_2_sums(:,:,atom_1) = atom_2_sums(:,:,atom_1)  &
                                & + force_constants(:,:,atom_2p,atom_1)
        enddo
      enddo
      sums = sum(atom_1_sums,3) ! also = sum(atom_2_sums,3)
      
      do atom_1=1,structure%no_atoms
        do atom_2=1,structure%no_atoms
          atom_2p = structure_sc%gvec_and_prim_to_atom(atom_2,j)
          
          force_constants(:,:,atom_2p,atom_1) =             &
             &   force_constants(:,:,atom_2p,atom_1)        &
             & - atom_1_sums(:,:,atom_2)/structure%no_atoms &
             & - atom_2_sums(:,:,atom_1)/structure%no_atoms &
             & + sums/structure%no_atoms
        enddo
      enddo
    enddo
    
    deallocate(atom_1_sums)
    deallocate(atom_2_sums)
    
    ! Mass-reduce force constants.
    do atom_1=1,structure%no_atoms
      do atom_2=1,structure_sc%no_atoms
        force_constants(:,:,atom_2,atom_1) =        &
           &   force_constants(:,:,atom_2,atom_1)   &
           & / dsqrt(structure%mass(atom_1)*structure_sc%mass(atom_2))
      enddo
    enddo
    
    ! Write out force constants for debugging purposes.
    force_constants_file = open_write_file(sdir//'/force_constants.dat')
    do atom_1=1,structure%no_atoms
      do atom_2=1,structure_sc%no_atoms
        call print_line(force_constants_file,atom_1//' '//atom_2)
        do j=1,3
          call print_line( force_constants_file, &
                         & force_constants(j,:,atom_2,atom_1))
        enddo
        call print_line(force_constants_file,'')
      enddo
    enddo
    
    ! Convert into mode-mode-Gvector representation.
    allocate(force_constants_2( structure%no_modes, &
                              & structure%no_modes, &
                              & structure_sc%sc_size))
    do atom_1=1,structure%no_atoms
      mode_1 = (atom_1-1)*3+1
      do atom_2=1,structure%no_atoms
        do j=1,structure_sc%sc_size
          mode_2 = (atom_2-1)*3+1
          atom_2p = structure_sc%gvec_and_prim_to_atom(atom_2,j)
          force_constants_2(mode_2:mode_2+2,mode_1:mode_1+2,j) = &
             & transpose(force_constants(:,:,atom_2p,atom_1))
        enddo
      enddo
    enddo
        
    call lte_4( structure,                       &
              & structure_sc,                    &
              & force_constants_2,               &
              & 0.0_dp,                          &
              & sdir//'/lte/freq_grids.dat',     &
              & sdir//'/lte/disp_patterns.dat',  &
              & sdir//'/lte/kdisp_patterns.dat', &
              & sdir//'/lte/pol_vec.dat',        &
              & sdir//'/lte/error.txt',          &
              & sdir//'/lte/dyn_mat.')
    
    if (file_exists(sdir//'/lte/error.txt')) then
      call print_line('There is an error in lte: check error.txt file.')
      stop
    endif
    
    deallocate(force_constants)
    deallocate(force_constants_2)
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
  call fourier_interpolation(                  &
     & structure,                              &
     & grid,                                   &
     & temperature,                            &
     & kpoints, sc_ids, gvector_ids,           &
     & str('lte/dyn_mat.'),                    &! Supercell_*/lte/dyn_mat.*.dat
     & disp_kpoints,                           &
     & symmetry_group,                         &
     & str('lte/phonon_dispersion_curve.dat'), &
     & str('lte/high_symmetry_points.dat'),    &
     & str('lte/free_energy.dat'),             &
     & str('lte/freq_dos.dat'))
end subroutine
end module
