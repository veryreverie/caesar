module phonon
  use constants
  use utils
  use linear_algebra
  use symmetry
  implicit none

contains

! ----------------------------------------------------------------------
! Construct the dynamical matrix at an arbitrary wave vector. We take into
! account all minimum image primitive cells in the supercell. 
! ----------------------------------------------------------------------
subroutine construct_dyn_mat(kpoint,basis,mass,no_cells,no_ims,cell_vecs,&
   &force_consts,dyn_mat)
  implicit none
  
  INTEGER,INTENT(in) :: basis,no_cells,no_ims(no_cells)
  REAL(dp),INTENT(in) :: kpoint(3),mass(basis),cell_vecs(3,8,no_cells),&
    &force_consts(basis,3,basis,3,no_cells)
  COMPLEX(dp),INTENT(out) :: dyn_mat(basis,3,basis,3)
  INTEGER :: i_atom,j_atom,i_cart,j_cart,i_cell,i_im
  REAL(dp) :: k_dot_r,prefactor
  COMPLEX(dp) :: phase
  
  dyn_mat=cmplx(0.d0,0.d0,dp)
  
  do i_atom=1,basis
    do i_cart=1,3
      do j_atom=1,basis
        do j_cart=1,3
          prefactor=1.d0/sqrt(mass(i_atom)*mass(j_atom))
          do i_cell=1,no_cells
            phase=cmplx(0.d0,0.d0,dp)
            do i_im=1,no_ims(i_cell)
              k_dot_r=dot_product(kpoint(1:3),cell_vecs(1:3,i_im,i_cell))
              phase=phase+cmplx(cos(k_dot_r),-sin(k_dot_r),dp)
            enddo ! i_im
            phase=phase/dble(no_ims(i_cell))
            dyn_mat(i_atom,i_cart,j_atom,j_cart)=&
              &dyn_mat(i_atom,i_cart,j_atom,j_cart)+&
              &force_consts(i_atom,i_cart,j_atom,j_cart,i_cell)*phase
          enddo ! i_cell
          dyn_mat(i_atom,i_cart,j_atom,j_cart)=prefactor*&
            &dyn_mat(i_atom,i_cart,j_atom,j_cart)
        enddo ! j_cart
      enddo ! j_atom
    enddo ! i_cart
  enddo ! i_atom
end subroutine

! ----------------------------------------------------------------------
! Diagonalise the dynamical matrix and calculate its eigenvalues.
! ----------------------------------------------------------------------
subroutine calculate_frequencies(basis,dyn_mat,freqs)
  use linear_algebra
  implicit none
  
  INTEGER,INTENT(in) :: basis
  COMPLEX(dp),INTENT(in) :: dyn_mat(basis,3,basis,3)
  REAL(dp),INTENT(out) :: freqs(3*basis)
  INTEGER :: i_atom,j_atom,i_cart,j_cart,i_index,j_index
  REAL(dp) :: minus_freqs_sq(3*basis)
  COMPLEX(dp) :: temp,temp_mat(3*basis,3*basis)
  
  type(ComplexEigenstuff) :: estuff

  i_index=0
  do i_atom=1,basis
    do i_cart=1,3
      i_index=i_index+1
      j_index=0
      do j_atom=1,basis
        do j_cart=1,3
          j_index=j_index+1
          temp_mat(i_index,j_index)=dyn_mat(i_atom,i_cart,j_atom,j_cart)
        enddo ! j_cart
      enddo ! j_atom
    enddo ! i_cart
  enddo ! i_atom

  do i_index=1,3*basis
    temp_mat(i_index,i_index)=cmplx(real(temp_mat(i_index,i_index)),0.d0,dp)
    do j_index=i_index+1,3*basis
      temp=0.5d0*(temp_mat(i_index,j_index)+conjg(temp_mat(j_index,i_index)))
      temp_mat(i_index,j_index)=temp
      temp_mat(j_index,i_index)=conjg(temp)
    enddo ! j_index
  enddo ! i_index
  
  estuff = calculate_eigenstuff(temp_mat)
  minus_freqs_sq = estuff%evals

  i_index=3*basis
  do j_index=1,3*basis
    if(minus_freqs_sq(i_index)>=0.d0)then
      freqs(j_index)=-sqrt(minus_freqs_sq(i_index))
    else
      freqs(j_index)=sqrt(-minus_freqs_sq(i_index))
    endif ! minus_freqs_sq
    i_index=i_index-1
  enddo ! j_index
end subroutine

subroutine generate_dispersion(rec_vecs,basis,mass,no_cells,no_ims,          &
   & cell_vecs,force_consts,path,phonon_dispersion_curve_filename, &
   & high_symmetry_points_filename)
  use file_module
  use string_module
  implicit none
  
  ! inputs
  real(dp),     intent(in) :: rec_vecs(3,3)
  integer,      intent(in) :: basis
  real(dp),     intent(in) :: mass(:)
  integer,      intent(in) :: no_cells
  integer,      intent(in) :: no_ims(:)
  real(dp),     intent(in) :: cell_vecs(:,:,:)
  real(dp),     intent(in) :: force_consts(:,:,:,:,:)
  real(dp),     intent(in) :: path(:,:)
  type(String), intent(in) :: phonon_dispersion_curve_filename
  type(String), intent(in) :: high_symmetry_points_filename
  
  integer :: i_path,i_cart,path_length,i_point
  real(dp) :: k_start(3),k_stop(3),k_diff(3),k_dist,total_k_dist,delta_k,&
    &kpoint(3),omega(3*basis)
  complex(dp) :: dyn_mat(basis,3,basis,3)
  
  integer :: no_points
  
  ! file units
  integer :: phonon_dispersion_curve_file
  integer :: high_symmetry_points_file
  
  no_points = size(path,2)
  
  total_k_dist=0.d0
  do i_point=1,no_points-1
    do i_cart=1,3
      k_start(i_cart)=dot_product(path(1:3,i_point),rec_vecs(1:3,i_cart))
      k_stop(i_cart)=dot_product(path(1:3,i_point+1),rec_vecs(1:3,i_cart))
    enddo ! i_cart
    k_diff=k_stop-k_start
    total_k_dist=total_k_dist+sqrt(dot_product(k_diff,k_diff))
  enddo ! i_point

  delta_k=total_k_dist/1000.d0

  phonon_dispersion_curve_file = open_write_file(phonon_dispersion_curve_filename)
  high_symmetry_points_file = open_write_file(high_symmetry_points_filename)

  total_k_dist=0.d0
  do i_point=1,no_points-1
    write(high_symmetry_points_file,*)i_point,total_k_dist
    do i_cart=1,3
      k_start(i_cart)=dot_product(path(1:3,i_point),rec_vecs(1:3,i_cart))
      k_stop(i_cart)=dot_product(path(1:3,i_point+1),rec_vecs(1:3,i_cart))
    enddo ! i_cart
    k_diff=k_stop-k_start
    k_dist=sqrt(dot_product(k_diff,k_diff))
    path_length=int(k_dist/delta_k)
    do i_path=0,path_length-1
      kpoint=k_start+dble(i_path)*k_diff/dble(path_length)
      call construct_dyn_mat(kpoint,basis,mass,no_cells,no_ims,cell_vecs,&
        &force_consts,dyn_mat)
      call calculate_frequencies(basis,dyn_mat,omega)
      write(phonon_dispersion_curve_file,*)total_k_dist,omega
      total_k_dist=total_k_dist+k_dist/dble(path_length)
    enddo ! i_path
  enddo ! i_point
  call construct_dyn_mat(k_stop,basis,mass,no_cells,no_ims,cell_vecs,&
    &force_consts,dyn_mat)
  call calculate_frequencies(basis,dyn_mat,omega)
  write(phonon_dispersion_curve_file,*)total_k_dist,omega
  write(high_symmetry_points_file,*)no_points,total_k_dist
  
  close(phonon_dispersion_curve_file)
  close(high_symmetry_points_file)
end subroutine

subroutine generate_dos(rec_vecs,basis,mass,no_cells,no_ims,cell_vecs, &
   &force_consts,temperature,free_energy_filename,freq_dos_filename)
  use string_module
  use file_module
  implicit none
  
  real(dp),     intent(in) :: temperature
  type(String), intent(in) :: free_energy_filename
  type(String), intent(in) :: freq_dos_filename
  
  INTEGER,PARAMETER :: no_bins=1000,no_prelims=10000,no_samples=1000000
  REAL(dp),PARAMETER :: freq_tol=1.d-8,safety_factor=1.1d0
  REAL(dp) :: T
  INTEGER,INTENT(in) :: basis,no_cells,no_ims(no_cells)
  REAL(dp),INTENT(in) :: rec_vecs(3,3),mass(basis),cell_vecs(3,8,no_cells),&
    &force_consts(basis,3,basis,3,no_cells)
  INTEGER :: i_sample,i_cart,i_freq,i_bin
  REAL(dp) :: max_freq,min_freq,frac(3),kpoint(3),freqs(3*basis),bin_width,&
    &rec_bin_width,freq_dos(no_bins),free_energy,omega
  COMPLEX(dp) :: dyn_mat(basis,3,basis,3)
  LOGICAL :: soft_modes
  
  ! file units
  integer :: free_energy_file
  integer :: freq_dos_file
  
  ! Initialise the random number generator
  call random_seed()
  
  ! Read in temperature
  T = temperature
  
  max_freq=-1.d0
  min_freq=huge(1.d0)
  
  do i_sample=1,no_prelims
    call random_number(frac)
    do i_cart=1,3
      kpoint(i_cart)=dot_product(frac(1:3),rec_vecs(1:3,i_cart))
    enddo ! i_cart
    call construct_dyn_mat(kpoint,basis,mass,no_cells,no_ims,cell_vecs,&
        &force_consts,dyn_mat)
    call calculate_frequencies(basis,dyn_mat,freqs)
    if(freqs(1)<min_freq)min_freq=freqs(1)
    if(freqs(3*basis)>max_freq)max_freq=freqs(3*basis)
  enddo ! i_sample
  soft_modes=(min_freq<-freq_tol)
  if(max_freq<=0.d0)call errstop('GENERATE_DOS','The system is pathologically &
    &unstable.')
  
  bin_width=safety_factor*max_freq/dble(no_bins)
  rec_bin_width=1.d0/bin_width
  freq_dos=0.d0
  
  do i_sample=1,no_samples
    call random_number(frac)
    do i_cart=1,3
      kpoint(i_cart)=dot_product(frac(1:3),rec_vecs(1:3,i_cart))
    enddo ! i_cart
    call construct_dyn_mat(kpoint,basis,mass,no_cells,no_ims,cell_vecs,&
        &force_consts,dyn_mat)
    call calculate_frequencies(basis,dyn_mat,freqs)
    soft_modes=(freqs(1)<-freq_tol)
    do i_freq=1,3*basis
      if(freqs(i_freq)>-freq_tol)then
        i_bin=max(1,ceiling(rec_bin_width*freqs(i_freq)))
        if(i_bin>no_bins)call errstop('GENERATE_DOS','Frequency too high to be &
          &binned.')
        freq_dos(i_bin)=freq_dos(i_bin)+1.d0
      endif ! freqs
    enddo ! i_freq
  enddo ! i_sample
  
  free_energy=0.d0
  do i_bin=1,no_bins
    omega=bin_width*(dble(i_bin)-0.5d0)
    free_energy=freq_dos(i_bin)*harmonic_free_energy(T,omega)/dble(no_samples)+&
      &free_energy
  enddo ! i_bin
  
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

real(dp) function harmonic_free_energy(T,omega)
  implicit none
  
  REAL(dp),PARAMETER :: tol=1.d-8
  REAL(dp),PARAMETER :: kB_au_per_K=3.16679002948702D-006
  REAL(dp),INTENT(in) :: T,omega
  REAL(dp) :: difference,kT
  
  if(T<tol)then
    harmonic_free_energy=0.5d0*omega
  else
    kT=kB_au_per_K*T
    difference=1.d0-exp(-omega/kT)
    if(difference>0.d0)then
      harmonic_free_energy=0.5d0*omega+kT*log(difference)
    else
      harmonic_free_energy=-huge(0.d0)
    endif
  endif
end function
end module

module fourier_interpolation_module
  use constants
  use utils
  implicit none
contains

! ----------------------------------------------------------------------
! Read in dynamical matrices at each k-point in the IBZ.
! ----------------------------------------------------------------------
subroutine read_dyn_mats(basis,mass,atom_prim_frac,no_kpoints,dyn_mats, &
   & kpoint_to_supercell,atoms_in_primitive_cell_fileroot,dyn_mat_fileroot,gvector_ids)
  use file_module
  use string_module
  implicit none
  
  type(String), intent(in) :: atoms_in_primitive_cell_fileroot
  type(String), intent(in) :: dyn_mat_fileroot
  integer,      intent(in) :: gvector_ids(:)
  
  type(String) :: filename
  
  ! file units
  integer :: atoms_in_primitive_cell_file
  integer :: dyn_mat_file
  
  INTEGER,INTENT(in) :: basis,no_kpoints,kpoint_to_supercell(no_kpoints)
  REAL(dp),INTENT(out) :: mass(basis),atom_prim_frac(3,basis)
  COMPLEX(dp),INTENT(out) :: dyn_mats(basis,3,basis,3,no_kpoints)
  REAL(dp),PARAMETER :: mass_tol=1.d-4,frac_tol=1.d-8
  INTEGER :: ierr,i_atom,i_cart,j_atom,j_cart,atom1,cart1,atom2,cart2,ibz_point,&
    &supercell,atom_map(basis)
  REAL(dp) :: real_part,imag_part,temp_mass,temp_frac(3)
  LOGICAL :: found_atom(basis)
  
  atoms_in_primitive_cell_file = open_read_file( atoms_in_primitive_cell_fileroot&
                                            & // '1.dat')
  do i_atom=1,basis
    read(atoms_in_primitive_cell_file,*,iostat=ierr) mass(i_atom), &
                                                   & atom_prim_frac(1:3,i_atom)
    if(ierr/=0)call errstop('READ_DYN_MATS','Problem reading &
    &atoms_in_primitive_cell.1.dat file.')
  enddo ! i_atom
  if(any(atom_prim_frac(1:3,1:basis)<-1.d-4.or.&
    &atom_prim_frac(1:3,1:basis)>=1.d0))call errstop('READ_DYN_MATS',&
      &'Fractional atomic coordinates are not in range [0.0,1.0)')
  close(atoms_in_primitive_cell_file)
  
  filename = str('Supercell_')//kpoint_to_supercell(1)//'/' &
     & //dyn_mat_fileroot//gvector_ids(1)//'/dat'
  dyn_mat_file = open_read_file(filename)
  do i_atom=1,basis
    do i_cart=1,3
      do j_atom=1,basis
        do j_cart=1,3
          read(dyn_mat_file,*,iostat=ierr) atom1,     &
                                         & cart1,     &
                                         & atom2,     &
                                         & cart2,     &
                                         & real_part, &
                                         & imag_part
          if(ierr/=0)call errstop('READ_DYN_MATS','Problem reading '// &
             & char(filename)//' file.')
          if(atom1/=i_atom.or.cart1/=i_cart.or.atom2/=j_atom.or.cart2/=j_cart)call &
            errstop('READ_DYN_MATS',char(filename)// &
               & ' file does not seem to be in the expected order.')
            dyn_mats(atom1,cart1,atom2,cart2,1)=cmplx(real_part,imag_part,dp)
        enddo ! j_cart
      enddo ! j_atom
    enddo ! i_cart
  enddo ! i_atom
  close(dyn_mat_file)
  
  do ibz_point=2,no_kpoints
    supercell=kpoint_to_supercell(ibz_point)
    atoms_in_primitive_cell_file = open_read_file( atoms_in_primitive_cell_fileroot&
                                              & // supercell//'.dat')
    found_atom(1:basis)=.false.
    atom_map(1:basis)=0
    do i_atom=1,basis
      read(atoms_in_primitive_cell_file,*,iostat=ierr)temp_mass,temp_frac(1:3)
      if(ierr/=0)call errstop('READ_DYN_MATS','Problem reading &
        &atoms_in_primitive_cell.'//trim(i2s(supercell))//'.dat file.')
      do j_atom=1,basis
        if(abs(temp_mass-mass(j_atom))<mass_tol)then
          if(all(abs(temp_frac(1:3)-atom_prim_frac(1:3,j_atom))<frac_tol))then
            found_atom(j_atom)=.true.
            atom_map(i_atom)=j_atom
          endif ! frac_tol
        endif ! mass_tol
      enddo ! j_atom
    enddo ! i_atom
    if(.not.any(found_atom(1:basis)))call errstop('READ_DYN_MATS','Unable to &
      &find all atoms in supercell '//trim(i2s(supercell))//'.')
    close(atoms_in_primitive_cell_file)
    
    filename = str('Supercell_')//kpoint_to_supercell(ibz_point)//'/'// &
       & dyn_mat_fileroot//gvector_ids(ibz_point)//'.dat'
    dyn_mat_file = open_read_file(filename)
    do i_atom=1,basis
      do i_cart=1,3
        do j_atom=1,basis
          do j_cart=1,3
            read(dyn_mat_file,*,iostat=ierr) atom1,     &
                                           & cart1,     &
                                           & atom2,     &
                                           & cart2,     &
                                           & real_part, &
                                           & imag_part
            if(ierr/=0)call errstop('READ_DYN_MATS','Problem reading '// &
               & char(filename)//' file')
            if(atom1/=i_atom.or.cart1/=i_cart.or.atom2/=j_atom.or.cart2/=j_cart)call &
              errstop('READ_DYN_MATS','dyn_mat.'//trim(i2s(ibz_point))//'.dat file &
                &does not seem to be in the expected order.')
              dyn_mats(atom_map(atom1),cart1,atom_map(atom2),cart2,ibz_point)=&
                &cmplx(real_part,imag_part,dp)
          enddo ! j_cart
        enddo ! j_atom
      enddo ! i_cart
    enddo ! i_atom
    close(dyn_mat_file)
  enddo ! ibz_point
end subroutine

subroutine fourier_interpolation(structure,grid,temperature,kpoints,sc_ids, &
   & gvector_ids, atoms_in_primitive_cell_fileroot,dyn_mat_fileroot,path,   &
   & phonon_dispersion_curve_filename,high_symmetry_points_filename,        &
   & free_energy_filename,freq_dos_filename)
  use constants, only : dp, pi
  use utils,     only : reduce_interval
  use linear_algebra
  use min_images
  use symmetry
  use phonon
  use file_module
  use structure_module
  use string_module
  implicit none
  
  ! filenames
  type(StructureData), intent(in) :: structure
  integer,             intent(in) :: grid(3)
  real(dp),            intent(in) :: temperature
  integer,             intent(in) :: kpoints(:,:)
  integer,             intent(in) :: sc_ids(:)
  integer,             intent(in) :: gvector_ids(:)
  type(String),        intent(in) :: atoms_in_primitive_cell_fileroot! append *.dat
  type(String),        intent(in) :: dyn_mat_fileroot                ! append *.dat
  real(dp),            intent(in) :: path(:,:)
  type(String),        intent(in) :: phonon_dispersion_curve_filename
  type(String),        intent(in) :: high_symmetry_points_filename
  type(String),        intent(in) :: free_energy_filename
  type(String),        intent(in) :: freq_dos_filename
  
  ! parameter
  REAL(dp),PARAMETER :: tol=1.d-8 
  
  ! variables
  integer :: i
  
  INTEGER,ALLOCATABLE :: atom_map_symm_backwards(:,:)
  INTEGER,ALLOCATABLE :: atom_map_symm_forwards(:,:)
  INTEGER,ALLOCATABLE :: grid_map_symm_backwards(:,:)
  INTEGER,ALLOCATABLE :: grid_map_symm_forwards(:,:)
  INTEGER,ALLOCATABLE :: grid_to_ibz_map(:)
  INTEGER,ALLOCATABLE :: ibz_to_grid_symm(:)
  INTEGER,ALLOCATABLE :: ibz_to_supercell_map(:)
  INTEGER,ALLOCATABLE :: identity_map(:)
  INTEGER,ALLOCATABLE :: no_im_cells(:)
  
  REAL(dp),ALLOCATABLE :: atom_pos_cart(:,:)
  REAL(dp),ALLOCATABLE :: atom_pos_frac(:,:)
  REAL(dp),ALLOCATABLE :: cell_pos_cart(:,:)
  REAL(dp),ALLOCATABLE :: min_im_cell_pos(:,:,:)
  REAL(dp),ALLOCATABLE :: force_consts(:,:,:,:,:)
  REAL(dp),ALLOCATABLE :: grid_points_cart(:,:)
  REAL(dp),ALLOCATABLE :: grid_points_frac(:,:)
  REAL(dp),ALLOCATABLE :: ibz_points_cart(:,:)
  REAL(dp),ALLOCATABLE :: ibz_points_frac(:,:)
  REAL(dp),ALLOCATABLE :: mass(:)
  
  COMPLEX(dp),ALLOCATABLE :: dyn_mats_grid(:,:,:,:,:)
  COMPLEX(dp),ALLOCATABLE :: dyn_mats_ibz(:,:,:,:,:)
  COMPLEX(dp),ALLOCATABLE :: phase(:,:)
  COMPLEX(dp),ALLOCATABLE :: temp_dyn_mat(:,:,:,:)
  COMPLEX(dp),ALLOCATABLE :: dyn_mats_symm(:,:,:,:,:)
  
  LOGICAL,ALLOCATABLE :: time_reversal(:)
  
  INTEGER :: counter
  INTEGER :: i_atom,j_atom
  INTEGER :: i_back
  INTEGER :: i_cart,j_cart
  INTEGER :: i_cell
  INTEGER :: i_grid
  INTEGER :: i_point
  INTEGER :: i_symm
  INTEGER :: m1,m2,m3
  
  REAL(dp) :: identity(3,3)
  REAL(dp) :: k_dot_r
  REAL(dp) :: kpoint(3)
  REAL(dp) :: prefactor
  REAL(dp) :: prim_rec_vecs(3,3)
  REAL(dp) :: super_latt_vecs(3,3)
  REAL(dp) :: super_rec_vecs(3,3)
  
  ! Input array lengths
  integer :: no_kpoints_path
  integer :: no_grid_points
  integer :: no_kpoints
  
  ! Read basic input files and allocate corresponding arrays
  
  no_grid_points=product(grid(:))
  
  allocate(identity_map(structure%no_atoms))
  allocate(grid_to_ibz_map(no_grid_points))
  allocate(ibz_to_grid_symm(no_grid_points))
  allocate(time_reversal(no_grid_points))
  allocate(grid_map_symm_backwards(no_grid_points,structure%no_symmetries))
  allocate(grid_map_symm_forwards(no_grid_points,structure%no_symmetries))
  allocate(grid_points_cart(3,no_grid_points))
  allocate(grid_points_frac(3,no_grid_points))
  allocate(dyn_mats_grid(structure%no_atoms,3,structure%no_atoms,3,no_grid_points))
  allocate(temp_dyn_mat(structure%no_atoms,3,structure%no_atoms,3))
  allocate(dyn_mats_symm(structure%no_atoms,3,structure%no_atoms,3,no_grid_points))
  allocate(force_consts(structure%no_atoms,3,structure%no_atoms,3,no_grid_points))
  
  identity(1,:) = (/ 1.0_dp, 0.0_dp, 0.0_dp /)
  identity(2,:) = (/ 0.0_dp, 1.0_dp, 0.0_dp /)
  identity(3,:) = (/ 0.0_dp, 0.0_dp, 1.0_dp /)
  
  do i_atom=1,structure%no_atoms
    identity_map(i_atom)=i_atom
  enddo
  
  prim_rec_vecs = 2*pi*structure%recip_lattice
  
  do i=1,3
    super_latt_vecs(i,:) = grid(i)*structure%lattice(i,:)
  enddo
  
  super_rec_vecs = 2*pi*transpose(inv_33(super_latt_vecs))
  
  i_grid=0
  do m1=0,grid(1)-1
    do m2=0,grid(2)-1
      do m3=0,grid(3)-1
        i_grid=i_grid+1
        grid_points_frac(:,i_grid) = reduce_interval( (/m1,m2,m3/)/dble(grid),&
                                                    & tol)
      enddo
    enddo
  enddo
  
  ! Convert grid points from fractional to Cartesian coodinates
  grid_points_cart = matmul(transpose(prim_rec_vecs),grid_points_frac)
  
  ! Determine which symmetry operations are in the point group and inverse
  ! group for each wave vector
  call kpoint_symmetry_maps(prim_rec_vecs,no_grid_points,grid_points_cart,&
    &structure%no_symmetries,structure%rotation_matrices,grid_map_symm_forwards,grid_map_symm_backwards)
  
  allocate(atom_pos_frac(3,structure%no_atoms))
  allocate(atom_pos_cart(3,structure%no_atoms))
  allocate(mass(structure%no_atoms))
  allocate(atom_map_symm_forwards(structure%no_atoms,structure%no_symmetries))
  allocate(atom_map_symm_backwards(structure%no_atoms,structure%no_symmetries))
  allocate(phase(structure%no_atoms,structure%no_atoms))
  allocate(cell_pos_cart(3,no_grid_points))
  allocate(no_im_cells(no_grid_points))
  allocate(min_im_cell_pos(3,8,no_grid_points))
  
  i_cell=0
  do m1=0,grid(1)-1
    do m2=0,grid(2)-1
      do m3=0,grid(3)-1
        i_cell=i_cell+1
        cell_pos_cart(:,i_cell) = matmul((/m1,m2,m3/),structure%lattice)
        call min_images_brute_force(cell_pos_cart(:,i_cell),super_latt_vecs,&
          &min_im_cell_pos(:,:,i_cell),no_im_cells(i_cell))
      enddo ! m3
    enddo ! m2
  enddo ! m1
  
  ! Get the number of k-points in the IBZ and allocate corresponding arrays
  no_kpoints = size(sc_ids)
  
  allocate(ibz_points_cart(3,no_kpoints))
  allocate(ibz_points_frac(3,no_kpoints))
  allocate(ibz_to_supercell_map(no_kpoints))
  allocate(dyn_mats_ibz(structure%no_atoms,3,structure%no_atoms,3,no_kpoints))
  
  ! Read input files related to k-points in the IBZ
  do i=1,size(kpoints)
    ibz_points_frac(:,i) = kpoints(:,i)/grid
  enddo
  ibz_to_supercell_map = sc_ids
  
  ! Convert IBZ points from fractional to Cartesian coodinates
  ibz_points_cart = matmul(transpose(prim_rec_vecs),ibz_points_frac)
  
  ! Read in the dynamical matrix at each k-point in the IBZ
  call read_dyn_mats(structure%no_atoms,mass,atom_pos_frac,no_kpoints,dyn_mats_ibz, &
      & ibz_to_supercell_map,atoms_in_primitive_cell_fileroot,dyn_mat_fileroot,gvector_ids)
  
  atom_pos_cart = matmul(transpose(structure%lattice),atom_pos_frac)
  
  ! Determine the mapping of the atoms in the primitive cell under each symmetry
  ! operation
    call atom_symmetry_maps(structure%lattice,structure%no_atoms,atom_pos_cart,structure%no_symmetries,&
    &structure%rotation_matrices,structure%offsets,atom_map_symm_forwards,atom_map_symm_backwards)

  ! Map each k-point on the grid to one in the IBZ
  call match_kpoints(prim_rec_vecs,no_grid_points,grid_points_cart,&
    &no_kpoints,ibz_points_cart,structure%no_symmetries,structure%rotation_matrices,grid_to_ibz_map,&
    &ibz_to_grid_symm,time_reversal)

  ! Determine the dynamical matrix at each point on the k-point grid by applying 
  ! the appropriate unitary trasformation to the dynamical matrix at a point in 
  ! the IBZ
  do i_grid=1,no_grid_points
    i_point=grid_to_ibz_map(i_grid)
    kpoint=ibz_points_cart(1:3,i_point)
    i_symm=ibz_to_grid_symm(i_grid)
    call g_matrix_phases(kpoint,structure%rotation_matrices(:,:,i_symm),structure%offsets(:,i_symm),&
      &structure%no_atoms,atom_pos_cart,phase)
    call apply_symmetry_to_dyn_mat(structure%no_atoms,atom_map_symm_forwards(:,i_symm),phase,&
      &dyn_mats_ibz(:,:,:,:,i_point),structure%rotation_matrices(:,:,i_symm),&
      &dyn_mats_grid(:,:,:,:,i_grid))
    if(time_reversal(i_grid))then
      dyn_mats_grid(:,:,:,:,i_grid)=&
        &cmplx(real(dyn_mats_grid(:,:,:,:,i_grid)),0.d0,dp)
    endif ! time_reversal
  enddo ! i_grid

  ! Symmetrize the dynamical matrix at each wave vector with respect to the
  ! translation operations of the supercell
  dyn_mats_symm=cmplx(0.d0,0.d0,dp)
  do i_grid=1,no_grid_points
    kpoint=grid_points_cart(1:3,i_grid)
    do i_cell=1,no_grid_points
      call g_matrix_phases(kpoint,identity,cell_pos_cart(:,i_cell),structure%no_atoms,&
        &atom_pos_cart,phase)
      call apply_symmetry_to_dyn_mat(structure%no_atoms,identity_map,phase,&
        &dyn_mats_grid(:,:,:,:,i_grid),identity,temp_dyn_mat)
      dyn_mats_symm(:,:,:,:,i_grid)=temp_dyn_mat+dyn_mats_symm(:,:,:,:,i_grid)
    enddo ! i_cell
    dyn_mats_symm(:,:,:,:,i_grid)=1.d0/dble(no_grid_points)*&
      &dyn_mats_symm(:,:,:,:,i_grid)
    if(time_reversal(i_grid))then
      dyn_mats_symm(:,:,:,:,i_grid)=&
        &cmplx(real(dyn_mats_symm(:,:,:,:,i_grid)),0.d0,dp)
    endif ! time_reversal
  enddo ! i_grid
  dyn_mats_grid=dyn_mats_symm
  
  ! Symmetrize the dynamical matrix at each wave vector with respect to its
  ! point and inverse groups
  dyn_mats_symm=cmplx(0.d0,0.d0,dp)
  do i_grid=1,no_grid_points
    counter=0
    do i_symm=1,structure%no_symmetries
      if(grid_map_symm_backwards(i_grid,i_symm)>0)then
        counter=counter+1
        i_back=grid_map_symm_backwards(i_grid,i_symm)
        kpoint=grid_points_cart(1:3,i_back)
        call g_matrix_phases(kpoint,structure%rotation_matrices(:,:,i_symm),structure%offsets(:,i_symm),&
          &structure%no_atoms,atom_pos_cart,phase)
        call apply_symmetry_to_dyn_mat(structure%no_atoms,atom_map_symm_forwards(:,i_symm),&
          &phase,dyn_mats_grid(:,:,:,:,i_back),structure%rotation_matrices(:,:,i_symm),temp_dyn_mat)
        dyn_mats_symm(:,:,:,:,i_grid)=temp_dyn_mat+dyn_mats_symm(:,:,:,:,i_grid)
      endif ! grid_map_symm_backwards
    enddo ! i_symm
    dyn_mats_symm(:,:,:,:,i_grid)=1.d0/dble(counter)*dyn_mats_symm(:,:,:,:,i_grid)
    if(time_reversal(i_grid))then
      dyn_mats_symm(:,:,:,:,i_grid)=&
        &cmplx(real(dyn_mats_symm(:,:,:,:,i_grid)),0.d0,dp)
    endif ! time_reversal
  enddo ! i_grid
  dyn_mats_grid=dyn_mats_symm
  
  ! Construct the matrix of force constants
  force_consts=0.d0
  do i_cell=1,no_grid_points
    do i_atom=1,structure%no_atoms
      do i_cart=1,3
        do j_atom=1,structure%no_atoms
          do j_cart=1,3
            prefactor=sqrt(mass(i_atom)*mass(j_atom))/&
              &dble(no_grid_points)
            do i_grid=1,no_grid_points
              kpoint=grid_points_cart(:,i_grid)
              k_dot_r=dot_product(kpoint(:),cell_pos_cart(:,i_cell))
              force_consts(i_atom,i_cart,j_atom,j_cart,i_cell)=&
                &force_consts(i_atom,i_cart,j_atom,j_cart,i_cell)+&
                &real(dyn_mats_grid(i_atom,i_cart,j_atom,j_cart,i_grid)*&
                &cmplx(cos(k_dot_r),sin(k_dot_r),dp))
            enddo ! i_grid
            force_consts(i_atom,i_cart,j_atom,j_cart,i_cell)=prefactor*&
              &force_consts(i_atom,i_cart,j_atom,j_cart,i_cell)
          enddo ! j_cart
        enddo ! j_atom
      enddo ! i_cart
    enddo ! i_atom
  enddo ! i_cell
  
  no_kpoints_path = size(path,2)
  
  call generate_dispersion(prim_rec_vecs,structure%no_atoms,mass,no_grid_points,    &
      & no_im_cells,min_im_cell_pos,force_consts,path, &
      & phonon_dispersion_curve_filename,high_symmetry_points_filename)

  call generate_dos(prim_rec_vecs,structure%no_atoms,mass,no_grid_points,no_im_cells,      &
    & min_im_cell_pos,force_consts,temperature,free_energy_filename, &
    & freq_dos_filename)
end subroutine
end module
