! ======================================================================
! Routines to calculate a number of properties of the crystal given a
!    matrix of force constants.
! ======================================================================

! --------------------------------------------------
! Original file comments.
! --------------------------------------------------

! Lattice Thermal Energy, Neil Drummond, 12/2004-1/2005.
! Modified by B. Monserrat

! Calculates the frequencies and patterns of atomic displacement at the 
!    supercell G vectors. This is useful if you are interested in following
!    a mode.

! CHANGES TO CODE
! ===============
! 20/07/05 Added calculation of speed of sound for monatomic crystals.
!          (Needs eigenvectors of dynamical matrix.)
! 21/07/05 Bugfix: take min image of all prim-cell vectors w.r.t supercell.
! 12/11/05 Fix for nonzero translations: r->b+Rr instead of r->b+R(r-b)
! 30/01/08 Divided freq_dos into separate sets, to allow error estimates.
! 02/02/08 Bugfix: include multiple images where necessary.
! 29/02/08 Applied fix r->b=Rr instead of r->b+R(r-b) everywhere.
! 17/04/09 Eliminated binning of imaginary frequencies, to avoid large free
!          energy due to lowest-freq bin at T>0.  Tidied, introduced i2s, etc.
!          Inserted better min-image routine.
! 08/05/09 Sped up initialisation by eliminating min_image calls.
! 13/05/09 Added evaluation of frequencies on supercell G vectors.
! 18/05/09 Fixed min-image bug introduced on 17/04/09.  Removed more min-images.
! 20/04/10 Enabled calculation of pattern of atomic displacement at each G.
! 17/06/11 Introduced BLAS & LAPACK.  Fixed bug in randomisation of theta for
!          speed-of-sound calculation.
! 2016     Integrated into Caesar. See Caesar git history.

! --------------------------------------------------

module lte_module
  use constants_module, only : dp
  use string_module
  use io_module
  use linear_algebra_module
  use normal_mode_module
  use dynamical_matrix_module
  implicit none
  
  ! Contains all the information generated by evaluate_freqs_on_grid.
  type LteReturn
    type(DynamicalMatrix), allocatable :: dynamical_matrices(:)
    type(NormalMode),      allocatable :: normal_modes(:,:)
  end type
  
  ! Contains the result of calculate_frequencies_and_polarisations
  type FreqsPolVecs
    real(dp),            allocatable :: frequencies(:)
    type(ComplexVector), allocatable :: polarisation_vectors(:,:)
  end type
  
contains

! ----------------------------------------------------------------------
! Returns delta_prim(k,b,a), the set of G-vectors which are equivalent to 
!    supercell G-vector k, modulo the supercell lattice.
! The only G-vectors given are such that the real-space distance between
!    atom a in the primitive (gamma-point) cell and atom b at G-vector k
!    is minimum, within a tolerance.
! Result is given in fractional supercell co-ordinates.
! ----------------------------------------------------------------------
function calculate_delta_prim(supercell) result(delta_prim)
  use structure_module
  use min_images_module
  implicit none
  
  type(StructureData), intent(in) :: supercell
  type(MinImages), allocatable    :: delta_prim(:,:,:)
  
  type(RealVector) :: delta_r
  type(RealVector) :: delta_r_corr
  integer  :: p,im
  
  integer :: atom_1_prim,atom_2_prim
  integer :: atom_1_sc_1,atom_2_sc_1,atom_2_sc_p
  
  integer :: ialloc

  ! Work out number of equivalent images and Delta Prim. Lattice Vectors
  ! for pairs of atoms (used in evaluation of force-constant matrix).
  allocate( delta_prim( supercell%sc_size,        &
          &             supercell%no_atoms_prim,  &
          &             supercell%no_atoms_prim), &
          & stat=ialloc); call err(ialloc)
  do atom_1_prim=1,supercell%no_atoms_prim
    atom_1_sc_1 = supercell%rvec_and_prim_to_atom(atom_1_prim,1)
    do atom_2_prim=1,supercell%no_atoms_prim
      atom_2_sc_1 = supercell%rvec_and_prim_to_atom(atom_2_prim,1)
      delta_r_corr = supercell%recip_lattice        &
                 & * ( supercell%atoms(atom_2_sc_1) &
                 &   - supercell%atoms(atom_1_sc_1))
      do p=1,supercell%sc_size
        ! Work out minimum distance(s) between atom_1 at gamma and
        !    atom_2 at G-vector p.
        atom_2_sc_p = supercell%rvec_and_prim_to_atom(atom_2_prim,p)
        delta_r = supercell%recip_lattice        &
              & * ( supercell%atoms(atom_2_sc_p) &
              &   - supercell%atoms(atom_1_sc_1))
        delta_prim(p,atom_2_prim,atom_1_prim) = &
           & min_images_brute_force(delta_r,supercell)
        
        ! Turn this into the corresponding difference(s) of lattice vectors.
        do im=1,size(delta_prim(p,atom_2_prim,atom_1_prim))
          delta_prim(p,atom_2_prim,atom_1_prim)%images(im) =    &
             & delta_prim(p,atom_2_prim,atom_1_prim)%images(im) &
                                       & - delta_r_corr
        enddo
      enddo
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Constructs the matrix of force constants in q-point co-ordinates,
!    given the matrix of force constants in supercell co-ordinates.
! ----------------------------------------------------------------------
function construct_dynamical_matrix(q,supercell, &
   & force_constants,delta_prim) result(dynamical_matrix)
  use constants_module, only : pi
  use structure_module
  use min_images_module
  implicit none
  
  type(RealVector),    intent(in)  :: q
  type(StructureData), intent(in)  :: supercell
  type(RealMatrix),    intent(in)  :: force_constants(:,:,:)
  type(MinImages),     intent(in)  :: delta_prim(:,:,:)
  type(DynamicalMatrix)            :: dynamical_matrix
  
  integer :: p,n,m,im
  complex(dp) :: tempc
  complex(dp), allocatable :: exp_iqr(:,:,:)
  real(dp) :: qr
  
  integer :: ialloc
  
  ! Precompute exp(-iq.(R-R')) to go in the dynamical matrix.
  allocate( exp_iqr( supercell%sc_size, &
          &          supercell%no_atoms_prim,   &
          &          supercell%no_atoms_prim),  &
          & stat=ialloc); call err(ialloc)
  do n=1,supercell%no_atoms_prim
    do m=1,supercell%no_atoms_prim
      do p=1,supercell%sc_size
        tempc = cmplx(0.0_dp,0.0_dp,dp)
        do im=1,size(delta_prim(p,m,n))
          qr = -2*pi*q*delta_prim(p,m,n)%images(im)
          tempc = tempc + cmplx(cos(qr),sin(qr),dp)
        enddo
        exp_iqr(p,m,n) = tempc / size(delta_prim(p,m,n))
      enddo
    enddo
  enddo
  
  ! Evaluate the dynamical matrix.
  allocate( dynamical_matrix%matrices( supercell%no_atoms_prim,  &
          &                            supercell%no_atoms_prim), &
          & stat=ialloc); call err(ialloc)
  dynamical_matrix%matrices = mat(cmplx(zeroes(3,3)))
  do n=1,supercell%no_atoms_prim
    do m=1,supercell%no_atoms_prim
      do p=1,supercell%sc_size
        dynamical_matrix%matrices(m,n) = dynamical_matrix%matrices(m,n)  &
                                     & + force_constants(p,m,n)          &
                                     & * exp_iqr(p,m,n)
      enddo
      
      ! Enforce Hermiticity on the dynamical matrix.
      if (n>=m) then
        dynamical_matrix%matrices(m,n) =                  &
           & (           dynamical_matrix%matrices(m,n)   &
           & + hermitian(dynamical_matrix%matrices(n,m))) &
           & / 2.0_dp
      endif
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Diagonalise the matrix of force constants, to obtain the normal mode
!    co-ordinates (eigenvectors) and harmonic frequencies (eigenvalues).
! ----------------------------------------------------------------------
function calculate_frequencies_and_polarisations(dynamical_matrix) &
   & result(output)
  use linear_algebra_module
  use structure_module
  implicit none
  
  type(DynamicalMatrix), intent(in) :: dynamical_matrix
  type(FreqsPolVecs)                :: output
  
  integer :: no_atoms,no_modes
  integer :: i,j,k,ialloc
  
  complex(dp), allocatable :: dyn_mat(:,:)
  
  type(ComplexEigenstuff) :: estuff
  
  no_atoms = size(dynamical_matrix%matrices,1)
  no_modes = 3*no_atoms
  
  ! Convert (3x3Matrix) x no_atoms x no_atoms to no_modes x no_modes
  allocate(dyn_mat(no_modes,no_modes),stat=ialloc); call err(ialloc)
  do i=1,no_atoms
    do j=1,no_atoms
      dyn_mat(3*j-2:3*j, 3*i-2:3*i) = cmplx(dynamical_matrix%matrices(j,i))
    enddo
  enddo
  
  ! Diagonalise dynamical matrix.
  estuff = calculate_eigenstuff(dyn_mat)
  
  ! Calculate frequencies and polarisation vectors.
  allocate( output%frequencies(3*no_atoms),                   &
          & output%polarisation_vectors(no_atoms,3*no_atoms), &
          & stat=ialloc); call err(ialloc)
  do i=1,no_modes
    k = no_modes - i + 1
    
    if (estuff%evals(k)>=0.0_dp) then
      ! Unstable mode.
      output%frequencies(i) = - sqrt(estuff%evals(k))
    else
      ! Stable mode.
      output%frequencies(i) = sqrt(- estuff%evals(k))
    endif
    
    do j=1,no_atoms
      output%polarisation_vectors(j,i) = estuff%evecs(3*j-2:3*j,k)
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Calculates the mean thermal energy of an isolated harmonic
!    oscillator of a given frequency.
! ----------------------------------------------------------------------
function harmonic_energy(temperature,frequency) result(output)
  use constants_module, only : kb_in_au
  implicit none
  
  real(dp), intent(in) :: temperature
  real(dp), intent(in) :: frequency
  real(dp)             :: output
  
  real(dp) :: denominator ! Bose factor = 1/denominator.
  
  if (temperature<=0) then
    ! Zero-point energy.
    output = frequency/2
  else
    denominator = exp(frequency/(kb_in_au*temperature))-1
    if (denominator>0) then
      ! General case.
      output = (1/denominator+0.5_dp)*frequency
    else
      ! High-temperature limit.
      output = kb_in_au*temperature
    endif
  endif
end function

! ----------------------------------------------------------------------
! Calculates the mean free energy of an isolated harmonic
!    oscillator of a given frequency.
! ----------------------------------------------------------------------
function harmonic_free_energy(temperature,frequency) result(output)
  use constants_module, only : kb_in_au
  implicit none
  
  real(dp), intent(in) :: temperature
  real(dp), intent(in) :: frequency
  real(dp)             :: output
  
  real(dp) :: difference
  real(dp) :: thermal_energy
  
  if (temperature<=0) then
    ! Zero-point energy.
    output = frequency/2
  else
    thermal_energy = kb_in_au*temperature
    difference = 1-exp(-frequency/thermal_energy)
    if (difference>0) then
      output = frequency/2 + thermal_energy*log(difference)
    else
      ! High-temperature limit.
      output = -huge(0.0_dp)
    endif
  endif
end function

! ----------------------------------------------------------------------
! Calculate the frequency density-of-states by Monte Carlo sampling of
!    the Brillouin zone.
! ----------------------------------------------------------------------
subroutine generate_dos(supercell,delta_prim, &
   & force_consts,temperature,free_energy_filename,freq_dos_filename)
  use ofile_module
  use structure_module
  use min_images_module
  implicit none
  
  type(StructureData), intent(in) :: supercell
  type(MinImages),     intent(in) :: delta_prim(:,:,:)
  type(RealMatrix),    intent(in) :: force_consts(:,:,:)
  real(dp),            intent(in) :: temperature
  type(String),        intent(in) :: free_energy_filename
  type(String),        intent(in) :: freq_dos_filename
  
  integer,parameter :: no_bins=1000,no_prelims=10000,no_samples=1000000
  real(dp),parameter :: freq_tol=1.0e-8_dp,safety_factor=1.15_dp
  
  integer :: i_sample,i_freq,i_bin
  real(dp) :: max_freq,min_freq,frac(3),bin_width,&
    &freq_dos(no_bins),free_energy,omega
  type(RealVector) :: qpoint
  
  type(DynamicalMatrix) :: dyn_mat
  type(FreqsPolVecs)    :: frequencies
  
  ! files.
  type(OFile) :: free_energy_file
  type(OFile) :: freq_dos_file
  
  ! Initialise the random number generator
  call random_seed()
  
  max_freq=-1.0_dp
  min_freq=huge(1.0_dp)
  
  ! Establish (approximate) maximum and minimum frequencies and hence
  ! choose the bin width.
  do i_sample=1,no_prelims
    call random_number(frac)
    qpoint = transpose(supercell%recip_supercell)*vec(frac)/supercell%sc_size
    dyn_mat = construct_dynamical_matrix(qpoint,supercell, &
       & force_consts,delta_prim)
    frequencies = calculate_frequencies_and_polarisations(dyn_mat)
    
    min_freq = min(min_freq,frequencies%frequencies(1))
    max_freq = max(max_freq,frequencies%frequencies(supercell%no_modes_prim))
  enddo
  
  if (max_freq<=0.0_dp) then
    call print_line('The system is pathologically unstable.')
    call err()
  endif
  
  bin_width=safety_factor*max_freq/dble(no_bins)
  freq_dos=0.0_dp
  
  do i_sample=1,no_samples
    call random_number(frac)
    qpoint = transpose(supercell%recip_supercell)*vec(frac)/supercell%sc_size
    dyn_mat = construct_dynamical_matrix(qpoint,supercell, &
       & force_consts,delta_prim)
    frequencies = calculate_frequencies_and_polarisations(dyn_mat)
    
    do i_freq=1,supercell%no_modes_prim
      ! Only bin positive frequencies.
      if (frequencies%frequencies(i_freq) > -freq_tol) then
        i_bin = max(1,ceiling(frequencies%frequencies(i_freq)/bin_width))
        if (i_bin>no_bins) then
          call print_line('Frequency too high to be binned.')
          call err()
        endif
        freq_dos(i_bin) = freq_dos(i_bin)+1.0_dp
      endif
    enddo
  enddo
  
  free_energy = 0.0_dp
  do i_bin=1,no_bins
    omega = bin_width*(dble(i_bin)-0.5_dp)
    free_energy = free_energy                               &
              & +   freq_dos(i_bin)                         &
              &   * harmonic_free_energy(temperature,omega) &
              &   / no_samples
  enddo
  
  free_energy_file = free_energy_filename
  call free_energy_file%print_line(free_energy)
  
  ! Normalise frequency DoS so that its integral is the number of
  !    degrees of freedom in the primitive cell. Note that the total
  !    number of frequencies sampled is no_samples*supercell%no_modes_prim.
  ! (Imaginary frequencies are ignored, however.)
  freq_dos = freq_dos/(no_samples*bin_width)
  
  ! Write out the frequency DoS.
  freq_dos_file = freq_dos_filename
  do i_bin=1,no_bins
    call freq_dos_file%print_line( bin_width*(i_bin-0.5_dp)//' '// &
                                 & freq_dos(i_bin))
  enddo
end subroutine

! ----------------------------------------------------------------------
! Use the frequency density-of-states to evaluate the lattice thermal
!    energy of the crystal as a function of the temperature.
! Repeat this for each set of frequency DoS data, to estimate the error
!    in the LTE.
! ----------------------------------------------------------------------
subroutine calc_lte(bin_width,temperature,freq_dos,tdependence1_filename)
  use constants_module, only : max_bin, no_fdos_sets
  use ofile_module
  implicit none
  
  real(dp), intent(in) :: bin_width
  real(dp), intent(in) :: temperature
  real(dp), intent(in) :: freq_dos(:,:)
  
  ! File name.
  type(String), intent(in) :: tdependence1_filename
  
  ! File.
  type(OFile) :: tdependence1_file
  
  integer :: bin,j
  real(dp) :: omega,lte_val,lte_sq,E_H(0:max_bin),lte,lte_err
  
  do bin=0,max_bin
    ! omega is the frequency in the middle of the corresponding bin.
    omega=(DBLE(bin)+0.5_dp)*bin_width
    ! Array of harmonic energies at each frequency.
    E_H(bin)=harmonic_energy(temperature,omega) 
  enddo ! bin
  lte=0.0_dp
  lte_sq=0.0_dp
  do j=1,no_fdos_sets
    ! LAPACK commented out because it isn't working. 9/1/2017
    ! lte_val=ddot(max_bin+1,freq_dos(0,j),1,E_H(0),1)
    lte_val = dot_product(freq_dos(:,j),E_H(:))
    lte=lte+lte_val ; lte_sq=lte_sq+lte_val**2
  enddo ! j
  lte=bin_width*lte/DBLE(no_fdos_sets)
  lte_sq=bin_width**2*lte_sq/DBLE(no_fdos_sets)
  lte_err=SQRT((lte_sq-lte**2)/DBLE(no_fdos_sets-1))
  call print_line('Done. LTE per primitive cell: '//lte//' +/- '//lte_err)
   
  tdependence1_file = tdependence1_filename
  call tdependence1_file%print_line(lte)
end subroutine

! ----------------------------------------------------------------------
! Use the frequency density-of-states to evaluate the lattice thermal
!    free energy of the crystal as a function of the temperature.
! Repeat this for each set of frequency DoS data, to estimate the error
!    in the LTFE.
! ----------------------------------------------------------------------
subroutine calc_ltfe(bin_width,temperature,freq_dos,tdependence2_filename)
  use constants_module, only : max_bin, no_fdos_sets
  use ofile_module
  implicit none
  
  real(dp), intent(in) :: bin_width
  real(dp), intent(in) :: temperature
  real(dp), intent(in) :: freq_dos(:,:)
  
  ! File name.
  type(String), intent(in) :: tdependence2_filename
  
  ! File.
  type(OFile) :: tdependence2_file
  
  integer :: bin,j
  real(dp) :: omega,ltfe_sq,ltfe_val,FE_H(0:max_bin),ltfe,ltfe_err
  
  do bin=0,max_bin
    ! omega is the frequency in the middle of the corresponding bin.
    omega=(DBLE(bin)+0.5_dp)*bin_width
    ! Array of harmonic energies at each frequency.
    FE_H(bin)=harmonic_free_energy(temperature,omega)
  enddo ! bin
  ltfe=0.0_dp
  ltfe_sq=0.0_dp
  do j=1,no_fdos_sets
    ltfe_val=DOT_PRODUCT(freq_dos(:,j),FE_H(:))
    ltfe=ltfe+ltfe_val ; ltfe_sq=ltfe_sq+ltfe_val**2
  enddo ! j
  ltfe=bin_width*ltfe/DBLE(no_fdos_sets)
  ltfe_sq=bin_width**2*ltfe_sq/DBLE(no_fdos_sets)
  ltfe_err=SQRT((ltfe_sq-ltfe**2)/DBLE(no_fdos_sets-1))
  call print_line('and LTFE per primitive cell   : '//ltfe//' +/- '//ltfe_err)
  
  tdependence2_file = tdependence2_filename
  call tdependence2_file%print_line(ltfe)
end subroutine

! ----------------------------------------------------------------------
! Generates a dispersion curve file, which contains all the branches of the 
!    phonon dispersion curve in a format that xmgrace can read.
! The branches of the dispersion curve are plotted against the total distance 
!    travelled along the specified paths in q-space.
! ----------------------------------------------------------------------
subroutine generate_dispersion(structure,supercell,&
   & delta_prim,force_consts,path,phonon_dispersion_curve_filename, &
   & high_symmetry_points_filename)
  use constants_module, only : pi
  use ofile_module
  use structure_module
  use min_images_module
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: supercell
  type(MinImages),     intent(in) :: delta_prim(:,:,:)
  type(RealMatrix),    intent(in) :: force_consts(:,:,:)
  type(RealVector),    intent(in) :: path(:)
  type(String),        intent(in) :: phonon_dispersion_curve_filename
  type(String),        intent(in) :: high_symmetry_points_filename
  
  integer                  :: no_paths
  type(RealVector), allocatable :: qpoints(:)
  real(dp),    allocatable :: path_length(:)
  real(dp),    allocatable :: cumulative_length(:)
  integer,     allocatable :: no_points(:)
  type(RealVector)         :: qpoint
  type(DynamicalMatrix)    :: dyn_mat
  type(FreqsPolVecs)       :: frequencies
  
  ! File units.
  type(OFile) :: phonon_dispersion_curve_file
  type(OFile) :: high_symmetry_points_file
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  no_paths = size(path)-1
  
  ! Transform q-points into reciprocal space (from fractional co-ords.)
  allocate(qpoints(no_paths+1), stat=ialloc); call err(ialloc)
  do i=1,no_paths+1
    qpoints(i) = 2*pi*transpose(structure%recip_lattice)*path(i)
  enddo
  
  ! Work out distances in reciprocal space.
  allocate(path_length(no_paths), stat=ialloc); call err(ialloc)
  do i=1,no_paths
    path_length(i) = l2_norm(qpoints(i+1)-qpoints(i))
  enddo
  
  allocate(cumulative_length(no_paths+1), stat=ialloc); call err(ialloc)
  cumulative_length(1) = 0.0_dp
  do i=2,no_paths+1
    cumulative_length(i) = cumulative_length(i-1)+path_length(i-1)
  enddo
  
  ! Space sampling points across the path, in proportion with path length.
  allocate(no_points(no_paths), stat=ialloc); call err(ialloc)
  do i=1,no_paths
    no_points(i) = nint(1000*path_length(i)/cumulative_length(no_paths+1))
  enddo
  
  ! Write path lengths to file.
  high_symmetry_points_file = high_symmetry_points_filename
  do i=1,no_paths+1
    call high_symmetry_points_file%print_line(i//' '//cumulative_length(i))
  enddo
  
  ! Travel along k-space paths, calculating frequencies at each point.
  phonon_dispersion_curve_file = phonon_dispersion_curve_filename
  do i=1,no_paths
    do j=0,no_points(i)-1
      qpoint = ((no_points(i)-j)*qpoints(i)+j*qpoints(i+1))/no_points(i)
      dyn_mat = construct_dynamical_matrix(qpoint,supercell, &
         & force_consts,delta_prim)
      frequencies = calculate_frequencies_and_polarisations(dyn_mat)
      call phonon_dispersion_curve_file%print_line(     &
         & cumulative_length(i)+j*path_length(i)//' '// &
         & frequencies%frequencies)
    enddo
  enddo
  
  ! Calculate frequencies at final k-space point.
  qpoint = qpoints(no_paths+1)
  dyn_mat = construct_dynamical_matrix(qpoint,supercell, &
     & force_consts,delta_prim)
  frequencies = calculate_frequencies_and_polarisations(dyn_mat)
  call phonon_dispersion_curve_file%print_line( &
       & (cumulative_length(no_paths+1))//' '// &
       & frequencies%frequencies)
end subroutine

! ----------------------------------------------------------------------
! Calculates the set of phonon normal modes at each supercell G-vector.
! Returns the real part of the non-mass-reduced polarisation vector, which
!    is the pattern of displacement corresponding to the normal mode.
! ----------------------------------------------------------------------
function evaluate_freqs_on_grid(supercell,force_constants) &
   & result(output)
  use structure_module
  use min_images_module
  implicit none
  
  type(StructureData), intent(in) :: supercell
  type(RealMatrix),    intent(in) :: force_constants(:,:,:)
  type(LteReturn)                 :: output
  
  ! Working variables.
  type(MinImages),     allocatable :: delta_prim(:,:,:)
  type(RealVector)                 :: qpoint
  type(FreqsPolVecs)               :: frequencies_polarisations
  type(ComplexVector), allocatable :: pol_vec(:,:)
  real(dp),            allocatable :: frequencies(:)
  logical,             allocatable :: translational(:)
  
  ! Indices.
  integer :: atom
  integer :: mode
  integer :: gvector,gvector_p
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  ! Allocate output and translational modes label.
  allocate( output%dynamical_matrices(supercell%sc_size),                   &
          & output%normal_modes(supercell%no_modes_prim,supercell%sc_size), &
          & stat=ialloc); call err(ialloc)
  do i=1,supercell%sc_size
    do j=1,supercell%no_modes_prim
      allocate( output%normal_modes(j,i)%displacements( &
              &               supercell%no_atoms_prim), &
              & stat=ialloc); call err(ialloc)
    enddo
  enddo
  
  ! Calculate minimum R-vector separations between atoms.
  delta_prim = calculate_delta_prim(supercell)
  
  ! Calculate dynamical matrices, and their eigenstuff.
  do gvector=1,supercell%sc_size
    gvector_p = supercell%paired_gvec(gvector)
    if (gvector_p<gvector) then
      cycle
    endif
    
    qpoint = dble(supercell%gvectors(gvector))
    output%dynamical_matrices(gvector) = construct_dynamical_matrix( &
                                                  & qpoint,          &
                                                  & supercell,       &
                                                  & force_constants, &
                                                  & delta_prim)
    if (gvector_p/=gvector) then
      output%dynamical_matrices(gvector_p) = output%dynamical_matrices(gvector)
      do i=1,supercell%no_atoms_prim
        do j=1,supercell%no_atoms_prim
          output%dynamical_matrices(gvector_p)%matrices(j,i) = &
             & hermitian(output%dynamical_matrices(gvector)%matrices(i,j))
        enddo
      enddo
    endif
    
    ! Calculate normal modes.
    frequencies_polarisations = calculate_frequencies_and_polarisations( &
       & output%dynamical_matrices(gvector))
    frequencies = frequencies_polarisations%frequencies
    pol_vec = frequencies_polarisations%polarisation_vectors
    
    ! Identify purely translational modes (at the gamma-point only).
    allocate(translational(size(frequencies)), stat=ialloc); call err(ialloc)
    translational = .false.
    if (all(int(supercell%gvectors(gvector))==0)) then
      do i=1,3
        j = minloc(frequencies,dim=1,mask=.not.translational)
        translational(j) = .true.
      enddo
    endif
  
    ! Copy frequency data to each normal mode.
    do i=1,supercell%no_modes_prim
      output%normal_modes(i,gvector)%frequency = frequencies(i)
      output%normal_modes(i,gvector)%soft_mode = frequencies(i) < -1.0e-6_dp
      output%normal_modes(i,gvector)%translational_mode = translational(i)
    enddo
    deallocate(translational, stat=ialloc); call err(ialloc)
    
    ! Calculate normal mode displacements from polarisation vectors.
    do mode=1,supercell%no_modes_prim
      do atom=1,supercell%no_atoms_prim
        output%normal_modes(mode,gvector)%displacements(atom) = &
                &   pol_vec(supercell%atom_to_prim(atom), mode) &
                & / sqrt(supercell%mass(atom))
      enddo
    enddo
    
    ! Calculate the normal mode at the paired g-vector.
    if (gvector_p/=gvector) then
      do i=1,supercell%no_modes_prim
        output%normal_modes(i,gvector_p) = output%normal_modes(i,gvector)
        do j=1,supercell%no_atoms_prim
          output%normal_modes(i,gvector_p)%displacements(j) = conjg( &
             & cmplx(output%normal_modes(i,gvector)%displacements(j)))
        enddo
      enddo
    endif
  enddo
end function

! ----------------------------------------------------------------------
! Calculates energy and free energy for a single supercell.
! ----------------------------------------------------------------------
subroutine calculate_lte_and_ltfe(supercell,force_constants, &
   & temperature,free_energy_filename,freq_dos_filename, &
   & tdependence1_filename,tdependence2_filename)
  use structure_module
  use min_images_module
  implicit none
  
  ! ----------------------------------------
  ! Inputs
  ! ----------------------------------------
  type(StructureData), intent(in) :: supercell
  type(RealMatrix),    intent(in) :: force_constants(:,:,:)
  real(dp),            intent(in) :: temperature
  
  ! ----------------------------------------
  ! filenames
  ! ----------------------------------------
  type(String), intent(in) :: free_energy_filename
  type(String), intent(in) :: freq_dos_filename
  type(String), intent(in) :: tdependence1_filename
  type(String), intent(in) :: tdependence2_filename
  
  ! ----------------------------------------
  ! previously global variables
  ! ----------------------------------------
  real(dp) :: bin_width
  real(dp), allocatable :: freq_dos(:,:)
  
  type(MinImages), allocatable :: delta_prim(:,:,:)
  
  delta_prim = calculate_delta_prim(supercell)

  call print_line('Temperature (K)                    :'//temperature)
  if (temperature<0.0_dp) then
    call print_line('Temperature should be non-negative.')
    call err()
  endif
  
  if (temperature<=0.0_dp) then
    call print_line('(i.e. the zero-point energy is to be calculated.)')
  endif
  
  call print_line('')
  call print_line('The mean thermal energy and the free energy will &
    &be calculated.')
  call print_line('Calculating the frequency density-of-states function...')
  call generate_dos(supercell,delta_prim, &
   & force_constants,temperature,free_energy_filename,freq_dos_filename)
  call print_line('Done.  Frequency density-of-states function written to &
    &freq_dos.dat.  (Please view this file using XMGrace.)')
  call print_line('')

  call print_line('Calculating the lattice thermal energy (LTE) and free energy &
    &(LTFE)...')
  call calc_lte(bin_width,temperature,freq_dos,tdependence1_filename)
  call calc_ltfe(bin_width,temperature,freq_dos,tdependence2_filename)
  call print_line('')
end subroutine

! ----------------------------------------------------------------------
! Calculates a dispersion curve for a single supercell.
! ----------------------------------------------------------------------
subroutine calculate_dispersion_curve(structure,supercell,force_constants, &
   & no_kspace_lines,path, &
   & phonon_dispersion_curve_filename,high_symmetry_points_filename)
  use structure_module
  use min_images_module
  implicit none
  
  ! ----------------------------------------
  ! Inputs
  ! ----------------------------------------
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: supercell
  type(RealMatrix),    intent(in) :: force_constants(:,:,:)
  integer,             intent(in) :: no_kspace_lines
  type(RealVector),    intent(in) :: path(:)
  
  ! ----------------------------------------
  ! filenames
  ! ----------------------------------------
  type(String), intent(in) :: phonon_dispersion_curve_filename
  type(String), intent(in) :: high_symmetry_points_filename
  
  ! ----------------------------------------
  ! previously global variables
  ! ----------------------------------------
  type(MinImages), allocatable :: delta_prim(:,:,:)
  
  integer :: i
  
  delta_prim = calculate_delta_prim(supercell)

  call print_line('Number of lines in k-space to plot     : '// &
     & no_kspace_lines)
  
  if (no_kspace_lines<1) then
    call print_line('Need to supply more lines in k-space!')
    call err()
  endif
  
  call print_line('Points along walk in reciprocal space &
    &(Cartesian components in a.u.):')
  do i=0,no_kspace_lines
    call print_line(path(i))
  enddo
  call print_line('Have read in points for dispersion curve.')
  
  call print_line('A dispersion curve will be calculated.')
  call print_line('Calculating the requested dispersion curve.')
  call generate_dispersion(structure,supercell,&
     & delta_prim,force_constants,path,phonon_dispersion_curve_filename, &
     & high_symmetry_points_filename)
  call print_line('Done.  dispersion_curve.dat has been generated.  (Please &
    &view this file using XMGrace.)')
  call print_line('')
end subroutine

! ----------------------------------------------------------------------
! Calculates density of states and phonon dispersion at the harmonic
!    approximation.
! Interpolates between phonons at each q-point.
! ----------------------------------------------------------------------
subroutine fourier_interpolation(dyn_mats_ibz,structure,temperature, &
   & structure_grid,qpoints_ibz, &
   & path,atom_symmetry_group,   &
   & phonon_dispersion_curve_filename,high_symmetry_points_filename,        &
   & free_energy_filename,freq_dos_filename)
  use linear_algebra_module
  use structure_module
  use group_module
  use min_images_module
  use construct_supercell_module
  use qpoints_module
  implicit none
  
  ! filenames
  type(DynamicalMatrix), intent(in) :: dyn_mats_ibz(:)
  type(StructureData),   intent(in) :: structure
  real(dp),              intent(in) :: temperature
  type(StructureData),   intent(in) :: structure_grid
  type(QpointData),      intent(in) :: qpoints_ibz(:)
  type(RealVector),      intent(in) :: path(:)
  type(Group),           intent(in) :: atom_symmetry_group(:)
  type(String),          intent(in) :: phonon_dispersion_curve_filename
  type(String),          intent(in) :: high_symmetry_points_filename
  type(String),          intent(in) :: free_energy_filename
  type(String),          intent(in) :: freq_dos_filename
  
  ! variables
  type(RealMatrix), allocatable :: rotations_cart(:)
  
  type(MinImages), allocatable :: delta_prim(:,:,:)
  type(RealMatrix), allocatable :: force_consts(:,:,:)
  
  type(DynamicalMatrix), allocatable :: dyn_mats_grid(:)
  complex(dp), allocatable :: phase(:,:,:)
  
  real(dp) :: exponent
  real(dp) :: qr
  
  integer :: i,j,k,l
  integer :: atom_1,atom_2,atom_1p,atom_2p
  integer :: rvec,gvec
  
  integer :: ialloc
  
  ! --------------------------------------------------
  ! Use symmetries to construct all dynamical matrices
  !   from the dynamical matrices of the IBZ.
  ! --------------------------------------------------
  
  ! Calculate the relative phases between atoms at each q-point.
  allocate( phase( structure%no_atoms,      &
          &        structure%no_atoms,      &
          &        structure_grid%sc_size), &
          & stat=ialloc); call err(ialloc)
  do i=1,structure_grid%sc_size
    do atom_1=1,structure%no_atoms
      do atom_2=1,structure%no_atoms
        ! Calculate k.dx
        exponent = ( structure_grid%gvectors(i)                        &
                 & * structure%recip_lattice                           &
                 & * (structure%atoms(atom_2)-structure%atoms(atom_1)) &
                 & ) / structure_grid%sc_size
        ! Calculate exp(i.k.dx)
        phase(atom_2,atom_1,i) = cmplx(cos(exponent),sin(exponent),dp)
      enddo
    enddo
  enddo
  
  ! Construct dynamical matrices.
  allocate( dyn_mats_grid(structure_grid%sc_size), &
          & stat=ialloc); call err(ialloc)
  rotations_cart = calculate_cartesian_rotations(structure)
  do i=1,size(qpoints_ibz)
    do j=1,size(qpoints_ibz(i)%gvectors)
      k = qpoints_ibz(i)%gvectors(j)
      l = qpoints_ibz(i)%rotations(j)
    
      ! Construct the element of the dynamical matrix from that in the IBZ.
      do atom_1=1,structure%no_atoms
        atom_1p = atom_symmetry_group(j) * atom_1
        do atom_2=1,structure%no_atoms
          atom_2p = atom_symmetry_group(j) * atom_2
          
          dyn_mats_grid(k)%matrices(atom_2p,atom_1p) =                        &
                                  &   rotations_cart(l)                       &
                                  & * dyn_mats_ibz(i)%matrices(atom_2,atom_1) &
                                  & * transpose(rotations_cart(l))            &
                                  & * phase(atom_2p,atom_2,k)                 &
                                  & * conjg(phase(atom_1p,atom_1,k))
          
          ! Apply time reversal symmetry if required.
          ! (imag(dyn_mats_grid) = 0)
          if (structure_grid%paired_gvec(i) == i) then
            dyn_mats_grid(i)%matrices(atom_2p,atom_1p) =                 &
               & cmplx(                                                  &
               & dble(real(dyn_mats_grid(i)%matrices(atom_2p,atom_1p))), &
               & 0.0_dp, dp)
          endif
        enddo
      enddo
    enddo
  enddo
  
  ! --------------------------------------------------
  ! Construct the matrix of force constants
  ! --------------------------------------------------
  allocate( force_consts( structure_grid%sc_size, &
          &               structure%no_atoms,     &
          &               structure%no_atoms),    &
          & stat=ialloc); call err(ialloc)
  force_consts = mat([ 0.0_dp,0.0_dp,0.0_dp, &
                     & 0.0_dp,0.0_dp,0.0_dp, &
                     & 0.0_dp,0.0_dp,0.0_dp  ], 3,3)
  
  do rvec=1,structure_grid%sc_size
    do i=1,structure%no_atoms
      do j=1,structure%no_atoms
        do gvec=1,structure_grid%sc_size
          qr = structure_grid%gvectors(gvec)  &
           & * structure_grid%recip_supercell &
           & * structure_grid%rvectors(rvec)
          
          force_consts(rvec,j,i) = force_consts(rvec,j,i)                  &
                               & + real( dyn_mats_grid(gvec)%matrices(j,k) &
                               &      * cmplx(cos(qr),sin(qr),dp))         &
                               & / real(structure_grid%sc_size,dp)
        enddo
      enddo
    enddo
  enddo
  
  ! --------------------------------------------------
  ! Calculate minimum image distances.
  ! --------------------------------------------------
  allocate( delta_prim( structure%no_atoms,       &
          &             structure%no_atoms,       &
          &             structure_grid%no_modes), &
          & stat=ialloc); call err(ialloc)
  delta_prim = calculate_delta_prim(structure_grid)
  
  ! --------------------------------------------------
  ! Generate dispersion and density of states.
  ! --------------------------------------------------
  call generate_dispersion(structure,structure_grid,&
     & delta_prim,force_consts,path,phonon_dispersion_curve_filename, &
     & high_symmetry_points_filename)
  
  call generate_dos(structure_grid,delta_prim, &
     & force_consts,temperature,free_energy_filename,freq_dos_filename)
end subroutine
end module
