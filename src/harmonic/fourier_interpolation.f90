module fourier_interpolation_module
  use constants
  use utils
  implicit none
contains

! ----------------------------------------------------------------------
! Determine the mapping under each symmetry operation for each k-point on the
! grid.
! ----------------------------------------------------------------------
function calculate_kpoint_symmetry_group(structure,structure_sc,grid) &
   & result(output)
  use structure_module
  use group_module
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: structure_sc
  integer,             intent(in) :: grid(3)
  type(Group), allocatable        :: output(:)
  
  integer :: i_grid,i_symm,j_grid
  
  real(dp) :: rotation(3,3)
  integer  :: rvec(3)
  
  integer, allocatable :: operations(:,:)
  
  allocate(operations(structure_sc%sc_size,structure%no_symmetries))
  operations = 0
  do i_symm=1,structure%no_symmetries
    do i_grid=1,structure_sc%sc_size
      ! Work out rotation in fractional lattice coordinates.
      rotation = matmul(matmul( structure%lattice, &
                              & structure%rotation_matrices(:,:,i_symm)), &
                              & transpose(structure%recip_lattice))
      ! Rotate the G-vector, and map it back to the IBZ.
      rvec = reduce_to_ibz( nint(matmul( rotation,                          &
                                       & structure_sc%gvectors(:,i_grid))), &
                          & grid)
      do j_grid=1,structure_sc%sc_size
        if (all(rvec == structure_sc%gvectors(:,j_grid))) then
          operations(i_grid,i_symm)=j_grid
        endif
      enddo
    enddo
  enddo
  
  allocate(output(structure%no_symmetries))
  do i_symm=1,structure%no_symmetries
    call new(output(i_symm),structure_sc%sc_size)
    output(i_symm) = operations(:,i_symm)
  enddo
end function

! ----------------------------------------------------------------------
! Read in dynamical matrices at each k-point in the IBZ.
! ----------------------------------------------------------------------
function read_dyn_mats(structure,sc_ids,gvector_ids,dyn_mat_fileroot) &
   & result(dyn_mats_ibz)
  use file_module
  use string_module
  use structure_module
  implicit none
  
  type(StructureData), intent(in) :: structure
  integer,             intent(in) :: sc_ids(:)
  integer,             intent(in) :: gvector_ids(:)
  type(String),        intent(in) :: dyn_mat_fileroot
  complex(dp), allocatable        :: dyn_mats_ibz(:,:,:,:,:)
  
  ! File information.
  type(String)              :: filename
  type(String), allocatable :: contents(:)
  type(String), allocatable :: line(:)
  
  ! Temporary variables
  integer :: i
  integer :: no_kpoints
  integer :: line_no
  
  integer :: atom_1,atom_2,j,k
  
  no_kpoints = size(sc_ids)
  
  allocate(dyn_mats_ibz(3,3,structure%no_atoms,structure%no_atoms,no_kpoints))
  dyn_mats_ibz = cmplx(0.0_dp,0.0_dp,dp)
  
  ! Loop across k-points
  do i=1,no_kpoints
    ! Construct the relevant filename, and read the file.
    filename = 'Supercell_'//sc_ids(i)//'/'// &
             & dyn_mat_fileroot//gvector_ids(i)//'.dat'
    contents = read_lines(filename)
    
    ! Read elements into dyn_mats_ibz
    line_no = 1
    do atom_1=1,structure%no_atoms
      do j=1,3
        do atom_2=1,structure%no_atoms
          do k=1,3
            line = split(contents(line_no))
            dyn_mats_ibz(j,k,atom_2,atom_1,i) = &
               & cmplx( dble(line(5)), dble(line(6)), dp)
            line_no = line_no + 1
          enddo
        enddo
      enddo
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Construct the dynamical matrix at an arbitrary wave vector, taking into
! account all minimum image primitive cells in the supercell. 
! ----------------------------------------------------------------------
function construct_dyn_mat(structure,structure_sc,kpoint, &
   & min_im_cell_pos,force_consts) result(dyn_mat)
  use structure_module
  use min_images_module
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: structure_sc
  real(dp),            intent(in) :: kpoint(3)
  type(MinImages),     intent(in) :: min_im_cell_pos(:)
  real(dp),            intent(in) :: force_consts(:,:,:)
  complex(dp), allocatable        :: dyn_mat(:,:)
  
  integer :: i_cell,i_im
  real(dp) :: k_dot_r
  
  complex(dp), allocatable :: phases(:)
  
  ! Calculate phases.
  allocate(phases(structure_sc%sc_size))
  phases = cmplx(0.0_dp,0.0_dp,dp)
  do i_cell=1,structure_sc%sc_size
    do i_im=1,size(min_im_cell_pos(i_cell))
      k_dot_r = dot_product(kpoint,min_im_cell_pos(i_cell)%images(:,i_im))
      phases(i_cell) = phases(i_cell) + cmplx(cos(k_dot_r), -sin(k_dot_r), dp)
    enddo
    phases(i_cell) = phases(i_cell)/size(min_im_cell_pos(i_cell))
  enddo
  
  ! Calculate dynamic matrix.
  allocate(dyn_mat(structure%no_modes,structure%no_modes))
  dyn_mat=cmplx(0.d0,0.d0,dp)
  do i_cell=1,structure_sc%sc_size
    dyn_mat = dyn_mat + force_consts(:,:,i_cell)*phases(i_cell)
  enddo
end function

! ----------------------------------------------------------------------
! Diagonalise the dynamical matrix and calculate its eigenvalues.
! ----------------------------------------------------------------------
function calculate_frequencies(dyn_mat) result(freqs)
  use linear_algebra
  implicit none
  
  complex(dp), intent(in) :: dyn_mat(:,:)
  real(dp), allocatable   :: freqs(:)
  
  integer                  :: no_modes
  type(ComplexEigenstuff)  :: estuff
  integer                  :: i,j
  
  no_modes = size(dyn_mat,1)
  
  ! Calculate eigenvalues and eigenvectors.
  estuff = calculate_eigenstuff(dyn_mat)
  
  ! Calculate frequencies.
  allocate(freqs(no_modes))
  i = no_modes ! i=no_modes-j+1
  do j=1,no_modes
    if(estuff%evals(i)>=0.d0)then
      freqs(j)=-dsqrt(estuff%evals(i))
    else
      freqs(j)=dsqrt(-estuff%evals(i))
    endif
    i = i-1
  enddo
end function

subroutine generate_dispersion(structure,structure_sc,&
   & min_im_cell_pos,force_consts,path,phonon_dispersion_curve_filename, &
   & high_symmetry_points_filename)
  use file_module
  use string_module
  use structure_module
  use min_images_module
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: structure_sc
  type(MinImages),     intent(in) :: min_im_cell_pos(:)
  real(dp),            intent(in) :: force_consts(:,:,:)
  real(dp),            intent(in) :: path(:,:)
  type(String),        intent(in) :: phonon_dispersion_curve_filename
  type(String),        intent(in) :: high_symmetry_points_filename
  
  integer                  :: no_paths
  real(dp),    allocatable :: kpoints(:,:)
  real(dp),    allocatable :: path_length(:)
  real(dp),    allocatable :: cumulative_length(:)
  integer,     allocatable :: no_points(:)
  real(dp)                 :: kpoint(3)
  complex(dp), allocatable :: dyn_mat(:,:)
  real(dp),    allocatable :: omega(:)
  
  ! File units.
  integer :: phonon_dispersion_curve_file
  integer :: high_symmetry_points_file
  
  ! Temporary variables.
  integer :: i,j
  
  no_paths = size(path,2)-1
  
  ! Transform k-points into reciprocal space.
  allocate(kpoints(3,no_paths+1))
  kpoints = 2*pi*matmul(structure%recip_lattice, path)
  
  ! Work out distances in k-space.
  allocate(path_length(no_paths))
  do i=1,no_paths
    path_length(i) = norm2(kpoints(:,i+1)-kpoints(:,i))
  enddo
  
  allocate(cumulative_length(no_paths+1))
  cumulative_length(1) = 0.0_dp
  do i=2,no_paths+1
    cumulative_length(i) = cumulative_length(i-1)+path_length(i-1)
  enddo
  
  ! Space sampling points across the path, in proportion with path length.
  allocate(no_points(no_paths))
  do i=1,no_paths
    no_points(i) = nint(1000*path_length(i)/cumulative_length(no_paths+1))
  enddo
  
  ! Write path lengths to file.
  high_symmetry_points_file = open_write_file(high_symmetry_points_filename)
  do i=1,no_paths+1
    call print_line(high_symmetry_points_file,i//' '//cumulative_length(i))
  enddo
  close(high_symmetry_points_file)
  
  ! Travel along k-space paths, calculating frequencies at each point.
  phonon_dispersion_curve_file = &
     & open_write_file(phonon_dispersion_curve_filename)
  do i=1,no_paths
    do j=0,no_points(i)-1
      kpoint = ((no_points(i)-j)*kpoints(:,i)+j*kpoints(:,i+1))/no_points(i)
      dyn_mat = construct_dyn_mat(structure,structure_sc,kpoint, &
         & min_im_cell_pos,force_consts)
      omega = calculate_frequencies(dyn_mat)
      call print_line(phonon_dispersion_curve_file, &
         & (cumulative_length(i)+j*path_length(i))//' '//omega)
    enddo
  enddo
  
  ! Calculate frequencies at final k-space point.
  kpoint = kpoints(:,no_paths+1)
  dyn_mat = construct_dyn_mat(structure,structure_sc,kpoint,min_im_cell_pos, &
     & force_consts)
  omega = calculate_frequencies(dyn_mat)
  call print_line(phonon_dispersion_curve_file, &
     & (cumulative_length(no_paths+1))//' '//omega)
  
  close(phonon_dispersion_curve_file)
end subroutine

subroutine generate_dos(structure,structure_sc,min_im_cell_pos, &
   & force_consts,temperature,free_energy_filename,freq_dos_filename)
  use string_module
  use file_module
  use structure_module
  use min_images_module
  use err_module
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: structure_sc
  type(MinImages),     intent(in) :: min_im_cell_pos(:)
  real(dp),            intent(in) :: force_consts(:,:,:)
  real(dp),            intent(in) :: temperature
  type(String),        intent(in) :: free_energy_filename
  type(String),        intent(in) :: freq_dos_filename
  
  integer,parameter :: no_bins=1000,no_prelims=10000,no_samples=1000000
  real(dp),parameter :: freq_tol=1.d-8,safety_factor=1.1d0
  
  integer :: i_sample,i_freq,i_bin
  real(dp) :: max_freq,min_freq,frac(3),kpoint(3),bin_width,&
    &rec_bin_width,freq_dos(no_bins),free_energy,omega
  LOGICAL :: soft_modes
  
  complex(dp), allocatable :: dyn_mat(:,:)
  real(dp),    allocatable :: freqs(:)
  
  ! file units
  integer :: free_energy_file
  integer :: freq_dos_file
  
  ! Initialise the random number generator
  call random_seed()
  
  max_freq=-1.d0
  min_freq=huge(1.d0)
  
  do i_sample=1,no_prelims
    call random_number(frac)
    kpoint = 2*pi*matmul(structure%recip_lattice,frac)
    dyn_mat = construct_dyn_mat(structure,structure_sc,kpoint, &
       & min_im_cell_pos,force_consts)
    freqs = calculate_frequencies(dyn_mat)
    
    if (freqs(1)<min_freq) then
      min_freq=freqs(1)
    endif
    
    if (freqs(3*structure%no_atoms)>max_freq) then
      max_freq = freqs(3*structure%no_atoms)
    endif
  enddo
  
  soft_modes=(min_freq<-freq_tol)
  if (max_freq<=0.d0) then
    call print_line('The system is pathologically unstable.')
    call err()
  endif
  
  bin_width=safety_factor*max_freq/dble(no_bins)
  rec_bin_width=1.d0/bin_width
  freq_dos=0.d0
  
  do i_sample=1,no_samples
    call random_number(frac)
    kpoint = 2*pi*matmul(structure%recip_lattice,frac)
    dyn_mat = construct_dyn_mat(structure,structure_sc,kpoint, &
       & min_im_cell_pos,force_consts)
    freqs = calculate_frequencies(dyn_mat)
    
    soft_modes=(freqs(1)<-freq_tol)
    do i_freq=1,3*structure%no_atoms
      if(freqs(i_freq)>-freq_tol)then
        i_bin = max(1,ceiling(rec_bin_width*freqs(i_freq)))
        if (i_bin>no_bins) then
          call print_line('Frequency too high to be binned.')
          call err()
        endif
        freq_dos(i_bin) = freq_dos(i_bin)+1.d0
      endif
    enddo
  enddo
  
  free_energy=0.d0
  do i_bin=1,no_bins
    omega=bin_width*(dble(i_bin)-0.5d0)
    free_energy = freq_dos(i_bin)                         &
              & * harmonic_free_energy(temperature,omega) &
              & / dble(no_samples)                        &
              & + free_energy
  enddo
  
  free_energy_file = open_write_file(free_energy_filename)
  write(free_energy_file,*) free_energy
  close(free_energy_file)
  
  freq_dos=freq_dos*rec_bin_width/dble(no_samples)
  
  freq_dos_file = open_write_file(freq_dos_filename)
  do i_bin=1,no_bins
    write(freq_dos_file,*)bin_width*(dble(i_bin)-0.5d0),freq_dos(i_bin)
  enddo
  close(freq_dos_file)
  
  if(soft_modes)write(*,*)'Soft modes present.'
end subroutine

real(dp) function harmonic_free_energy(temperature,omega)
  implicit none
  
  REAL(dp),PARAMETER :: tol=1.d-8
  REAL(dp),PARAMETER :: kB_au_per_K=3.16679002948702D-006
  REAL(dp),INTENT(in) :: temperature,omega
  REAL(dp) :: difference,kT
  
  if(temperature<tol)then
    harmonic_free_energy=0.5d0*omega
  else
    kT=kB_au_per_K*temperature
    difference=1.d0-exp(-omega/kT)
    if(difference>0.d0)then
      harmonic_free_energy=0.5d0*omega+kT*log(difference)
    else
      harmonic_free_energy=-huge(0.d0)
    endif
  endif
end function

subroutine fourier_interpolation(dyn_mats_ibz,structure,grid,temperature,kpoints, &
   & path,atom_symmetry_group,   &
   & phonon_dispersion_curve_filename,high_symmetry_points_filename,        &
   & free_energy_filename,freq_dos_filename)
  use constants, only : dp
  use utils,     only : reduce_interval
  use linear_algebra
  use file_module
  use structure_module
  use string_module
  use supercell_module
  use group_module
  use min_images_module
  use construct_supercell_module
  implicit none
  
  ! filenames
  complex(dp),         intent(in) :: dyn_mats_ibz(:,:,:)
  type(StructureData), intent(in) :: structure
  integer,             intent(in) :: grid(3)
  real(dp),            intent(in) :: temperature
  integer,             intent(in) :: kpoints(:,:)
  real(dp),            intent(in) :: path(:,:)
  type(Group),         intent(in) :: atom_symmetry_group(:)
  type(String),        intent(in) :: phonon_dispersion_curve_filename
  type(String),        intent(in) :: high_symmetry_points_filename
  type(String),        intent(in) :: free_energy_filename
  type(String),        intent(in) :: freq_dos_filename
  
  ! variables
  type(MinImages), allocatable :: min_im_cell_pos(:)
  real(dp),        allocatable :: force_consts(:,:,:)
  
  complex(dp), allocatable :: dyn_mats_grid(:,:,:)
  complex(dp), allocatable :: phase(:,:,:)
  
  INTEGER :: i_cell
  INTEGER :: i_grid
  
  REAL(dp) :: k_dot_r
  real(dp) :: exponent
  real(dp) :: gvec_cart(3)
  
  integer :: i,j,k,l
  integer :: kpoint_id
  integer :: atom_1,atom_2,atom_1p,atom_2p
  integer :: mode_1,mode_2,mode_1p,mode_2p
  
  ! Supercell data
  type(SupercellData) :: grid_supercell
  type(StructureData) :: structure_sc
  
  ! Symmetry group data.
  type(Group), allocatable :: kpoint_symmetry_group(:)
  
  ! --------------------------------------------------
  ! Construct a supercell across the entire grid.
  ! --------------------------------------------------
  call new(grid_supercell,product(grid))
  grid_supercell%supercell(1,:) = (/ grid(1), 0      , 0       /)
  grid_supercell%supercell(2,:) = (/ 0      , grid(2), 0       /)
  grid_supercell%supercell(3,:) = (/ 0      , 0      , grid(3) /)
  
  grid_supercell%recip_supercell = invert_int(grid_supercell%supercell)
  
  l = 1
  do i=0,grid(1)-1
    do j=0,grid(2)-1
      do k=0,grid(3)-1
        grid_supercell%gvectors(:,l) = (/ i,j,k /)
      enddo
    enddo
  enddo
  
  structure_sc = construct_supercell(structure,grid_supercell)
  
  ! --------------------------------------------------
  ! Calculate the symmetry group for the kpoints.
  ! --------------------------------------------------
  kpoint_symmetry_group = calculate_kpoint_symmetry_group( structure,    &
                                                         & structure_sc, &
                                                         & grid)
  
  ! --------------------------------------------------
  ! Read in the dynamical matrix at each k-point in the IBZ
  ! --------------------------------------------------
  !dyn_mats_ibz = read_dyn_mats(structure,sc_ids,gvector_ids,dyn_mat_fileroot)
  
  ! --------------------------------------------------
  ! Use symmetries to construct all dynamical matrices
  !   from the dynamical matrices of the IBZ.
  ! --------------------------------------------------
  
  ! Calculate the relative phases between atoms at each k-point.
  allocate(phase(structure%no_atoms,structure%no_atoms,structure_sc%sc_size))
  do i=1,structure_sc%sc_size
    do atom_1=1,structure%no_atoms
      do atom_2=1,structure%no_atoms
        ! Calculate k.dx
        exponent = dot_product( kpoints(:,i), &
                              &   matmul( structure%recip_lattice, &
                                         &   structure%atoms(:,atom_2) &
                                         & - structure%atoms(:,atom_1)) &
                              & / grid)
        ! Calculate exp(i.k.dx)
        phase(atom_2,atom_1,i) = cmplx(cos(exponent),sin(exponent),dp)
      enddo
    enddo
  enddo
  
  ! Construct dynamical matrices.
  allocate(dyn_mats_grid(structure%no_modes,structure%no_modes,structure_sc%sc_size))
  do_i : do i=1,structure_sc%sc_size
    do j=1,structure%no_symmetries
      kpoint_id = operate(kpoint_symmetry_group(j),i)
      if (kpoint_id==0) then
        cycle
      endif
      
      ! Find the rotated kpoint in the IBZ.
      do k=1,size(kpoints,2)
        if (.not. all(structure_sc%gvectors(:,kpoint_id)==kpoints(:,k))) then
          cycle
        endif
        
        ! Construct the element of the dynamical matrix from that in the IBZ.
        do atom_1=1,structure%no_atoms
          mode_1 = (atom_1-1)*3+1
          atom_1p = operate(atom_symmetry_group(j),atom_1)
          mode_1p = (atom_1p-1)*3+1
          do atom_2=1,structure%no_atoms
            mode_2 = (atom_2-1)*3+1
            atom_2p = operate(atom_symmetry_group(j),atom_2)
            mode_2p = (atom_2p-1)*3+1
            
            dyn_mats_grid(mode_2p:mode_2p+2,mode_1p:mode_1p+2,i) =    &
               & matmul(matmul(                                       &
               &    structure%rotation_matrices(:,:,j),               &
               &    dyn_mats_ibz(mode_2:mode_2+2,mode_1:mode_1+2,k)), &
               &    transpose(structure%rotation_matrices(:,:,j)))    &
               & * phase(atom_2p,atom_2,i)*conjg(phase(atom_1p,atom_1,i))
          enddo
        enddo
        
        ! Apply time reversal symmetry if required.
        if (all(modulo(structure_sc%gvectors(:,i)*2,grid)==0)) then
          dyn_mats_grid(:,:,i) = cmplx(real(dyn_mats_grid(:,:,i)), 0.0_dp, dp)
        endif
        
        cycle do_i
      enddo
    enddo
  enddo do_i
  
  ! --------------------------------------------------
  ! Construct the matrix of force constants
  ! --------------------------------------------------
  allocate(force_consts( structure%no_modes, &
                       & structure%no_modes, &
                       & structure_sc%sc_size))
  force_consts=0.0_dp
  
  do i_cell=1,structure_sc%sc_size
    do i_grid=1,structure_sc%sc_size
      k_dot_r = dot_product(kpoints(:,i_grid),structure_sc%gvectors(:,i_cell))
      force_consts(:,:,i_cell) = force_consts(:,:,i_cell) &
                             & + real( dyn_mats_grid(:,:,i_grid) &
                             &       * cmplx(cos(k_dot_r),sin(k_dot_r),dp))
    enddo
  enddo
  
  do atom_1=1,structure%no_atoms
    mode_1 = (atom_1-1)*3+1
    do atom_2=1,structure%no_atoms
      mode_2 = (atom_2-1)*3+1
      force_consts(mode_2:mode_2+2,mode_1:mode_1+2,:) = &
         &   force_consts(mode_2:mode_2+2,mode_1:mode_1+2,:) &
         & / structure_sc%sc_size
    enddo
  enddo
  
  ! --------------------------------------------------
  ! Calculate minimum image distances.
  ! --------------------------------------------------
  allocate(min_im_cell_pos(structure_sc%sc_size))
  do i=1,structure_sc%sc_size
    gvec_cart = matmul(transpose(structure%lattice),structure_sc%gvectors(:,i))
    min_im_cell_pos(i) = min_images_brute_force(gvec_cart,structure_sc)
  enddo
  
  ! --------------------------------------------------
  ! Generate dispersion and density of states.
  ! --------------------------------------------------
  call generate_dispersion(structure,structure_sc,&
     & min_im_cell_pos,force_consts,path,phonon_dispersion_curve_filename, &
     & high_symmetry_points_filename)
  
  call generate_dos(structure,structure_sc,min_im_cell_pos, &
     & force_consts,temperature,free_energy_filename,freq_dos_filename)
end subroutine
end module
