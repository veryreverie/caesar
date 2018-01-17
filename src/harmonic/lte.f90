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
contains

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
subroutine generate_dos(supercell,min_images,force_constants,temperature, &
   & free_energy_filename,dos_filename,log_filename)
  use ofile_module
  use structure_module
  use min_images_module
  use force_constants_module
  implicit none
  
  type(StructureData),  intent(in) :: supercell
  type(MinImages),      intent(in) :: min_images(:,:)
  type(ForceConstants), intent(in) :: force_constants
  real(dp),             intent(in) :: temperature
  type(String),         intent(in) :: free_energy_filename
  type(String),         intent(in) :: dos_filename
  type(String),         intent(in) :: log_filename
  
  integer,  parameter :: no_bins=1000,no_prelims=10000,no_samples=100000
  integer,  parameter :: print_every = 10000
  real(dp), parameter :: freq_tol=1.0e-8_dp,safety_margin=0.15_dp
  
  integer :: i_sample,i_freq,i_bin
  real(dp) :: max_freq,min_freq,freq_spread,frac(3),bin_width,&
    &freq_dos(no_bins),free_energy,omega
  type(RealVector) :: qpoint
  
  type(DynamicalMatrix) :: dyn_mat
  
  ! files.
  type(OFile) :: free_energy_file
  type(OFile) :: dos_file
  type(OFile) :: logfile
  
  ! Open logfile.
  logfile = log_filename
  
  ! Initialise the random number generator
  call random_seed()
  
  max_freq=-1.0_dp
  min_freq=huge(1.0_dp)
  
  ! Establish (approximate) maximum and minimum frequencies and hence
  ! choose the bin width.
  do i_sample=1,no_prelims
    call random_number(frac)
    qpoint = vec(frac)
    dyn_mat = DynamicalMatrix( qpoint,          &
                             & supercell,       &
                             & force_constants, &
                             & min_images,      &
                             & logfile)
    
    min_freq = min( min_freq, &
                  & dyn_mat%complex_modes(1)%frequency)
    max_freq = max( max_freq, &
                  & dyn_mat%complex_modes(supercell%no_modes_prim)%frequency)
  enddo
  
  if (max_freq<=0.0_dp) then
    call print_line(WARNING//': The system is pathologically unstable.')
  endif
  
  ! Spread out min and max frequencies to leave safety margin.
  freq_spread = max_freq-min_freq
  min_freq = max_freq-(1+safety_margin)*freq_spread
  max_freq = min_freq+(1+2*safety_margin)*freq_spread
  freq_spread = freq_spread*(1+2*safety_margin)
  bin_width=(max_freq-min_freq)/no_bins
  freq_dos=0.0_dp
  
  do i_sample=1,no_samples
    call random_number(frac)
    qpoint = vec(frac)
    dyn_mat = DynamicalMatrix( qpoint,          &
                             & supercell,       &
                             & force_constants, &
                             & min_images,      &
                             & logfile)
    
    do i_freq=1,supercell%no_modes_prim
      i_bin = ceiling( (dyn_mat%complex_modes(i_freq)%frequency-min_freq) &
                   & / bin_width)
      if (i_bin<1) then
        call print_line(ERROR//': Frequency too low to be binned.')
        call err()
      elseif (i_bin>no_bins) then
        call print_line(ERROR//': Frequency too high to be binned.')
        call err()
      endif
      freq_dos(i_bin) = freq_dos(i_bin)+1.0_dp
    enddo
    
    if (modulo(i_sample,print_every)==0) then
      call print_line('Calculating DOS: '//i_sample//' of '//no_samples// &
         & ' samples complete.')
    endif
  enddo
  
  !free_energy = 0.0_dp
  !do i_bin=1,no_bins
  !  omega = bin_width*(dble(i_bin)-0.5_dp)
  !  free_energy = free_energy                               &
  !            & +   freq_dos(i_bin)                         &
  !            &   * harmonic_free_energy(temperature,omega) &
  !            &   / no_samples
  !enddo
  !
  !free_energy_file = free_energy_filename
  !call free_energy_file%print_line(free_energy)
  
  ! Normalise frequency DoS so that its integral is the number of
  !    degrees of freedom in the primitive cell. Note that the total
  !    number of frequencies sampled is no_samples*supercell%no_modes_prim.
  ! (Imaginary frequencies are ignored, however.)
  freq_dos = freq_dos/(no_samples*bin_width)
  
  ! Write out the frequency DoS.
  dos_file = dos_filename
  do i_bin=1,no_bins
    call dos_file%print_line(                    &
       & 'Bin: '//min_freq+bin_width*(i_bin-1)// &
       & ' to '//min_freq+bin_width*i_bin)
    call dos_file%print_line( 'Bin DOS: '//freq_dos(i_bin))
    call dos_file%print_line( '')
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
subroutine generate_dispersion(structure,supercell,min_images,force_constants, &
   & path_labels,path_qpoints,dispersion_filename,                             &
   & high_symmetry_points_filename,log_filename)
  use constants_module, only : pi
  use ofile_module
  use structure_module
  use min_images_module
  use force_constants_module
  implicit none
  
  type(StructureData),  intent(in) :: structure
  type(StructureData),  intent(in) :: supercell
  type(MinImages),      intent(in) :: min_images(:,:)
  type(ForceConstants), intent(in) :: force_constants
  type(String),         intent(in) :: path_labels(:)
  type(RealVector),     intent(in) :: path_qpoints(:)
  type(String),         intent(in) :: dispersion_filename
  type(String),         intent(in) :: high_symmetry_points_filename
  type(String),         intent(in) :: log_filename
  
  ! Path variables.
  integer,  parameter :: total_no_points = 1000
  real(dp), parameter :: fractional_separation = 1.0_dp/total_no_points
  
  integer               :: no_segments
  integer               :: no_vertices
  real(dp), allocatable :: fractional_lengths(:)
  real(dp), allocatable :: fractional_distances(:)
  integer,  allocatable :: points_per_segment(:)
  
  ! q-point and dynamical matrix variables.
  type(RealVector)      :: qpoint
  type(DynamicalMatrix) :: dyn_mat
  
  ! File units.
  type(OFile) :: dispersion_file
  type(OFile) :: high_symmetry_points_file
  type(OFile) :: logfile
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  ! Open logfile.
  logfile = log_filename
  
  no_vertices = size(path_labels)
  no_segments = size(path_labels)-1
  if (size(path_qpoints) /= no_vertices) then
    call print_line(CODE_ERROR//': The number of q-points does not match the &
       &number of q-point labels.')
    call err()
  endif
  
  ! Work out distances in terms of fractions of the path.
  allocate( fractional_lengths(no_segments),   &
          & fractional_distances(no_vertices), &
          & points_per_segment(no_segments),   &
          & stat=ialloc); call err(ialloc)
  
  do i=1,no_segments
    fractional_lengths(i) = l2_norm(path_qpoints(i+1)-path_qpoints(i))
  enddo
  fractional_lengths = fractional_lengths / sum(fractional_lengths)
  
  fractional_distances(1) = 0.0_dp
  do i=1,no_segments
    fractional_distances(i+1) = fractional_distances(i) &
                            & + fractional_lengths(i)
  enddo
  
  ! Space sampling points along the path, in proportion with path length.
  do i=1,no_segments
    points_per_segment(i) = nint(total_no_points*fractional_lengths(i))
  enddo
  
  ! Write path to file.
  high_symmetry_points_file = high_symmetry_points_filename
  do i=1,no_vertices
    call high_symmetry_points_file%print_line( &
       & 'q-point: '//path_labels(i)//' '//path_qpoints(i))
    call high_symmetry_points_file%print_line( &
       & 'Fraction along path: '//fractional_distances(i))
    call high_symmetry_points_file%print_line('')
  enddo
  
  ! Travel along q-space paths, calculating frequencies at each point.
  dispersion_file = dispersion_filename
  do i=1,no_segments
    do j=0,points_per_segment(i)-1
      qpoint = ( (points_per_segment(i)-j)*path_qpoints(i)     &
           &   + j                        *path_qpoints(i+1) ) &
           & / points_per_segment(i)
      dyn_mat = DynamicalMatrix( qpoint,          &
                               & supercell,       &
                               & force_constants, &
                               & min_images,      &
                               & logfile)
      call dispersion_file%print_line( &
         & 'Fraction along path: '//                &
         & fractional_distances(i)+j*fractional_separation)
      call dispersion_file%print_line( &
         & 'Frequencies: '//dyn_mat%frequencies())
    enddo
  enddo
  
  ! Calculate frequencies at final k-space point.
  qpoint = path_qpoints(no_vertices)
  dyn_mat = DynamicalMatrix( qpoint,          &
                           & supercell,       &
                           & force_constants, &
                           & min_images,      &
                           & logfile)
  call dispersion_file%print_line( &
     & 'Fraction along path: '//1.0_dp)
  call dispersion_file%print_line( &
     & 'Frequencies: '//dyn_mat%frequencies())
end subroutine

! ----------------------------------------------------------------------
! Calculates energy and free energy for a single supercell.
! ----------------------------------------------------------------------
subroutine calculate_lte_and_ltfe(supercell,force_constants, &
   & temperature,free_energy_filename,freq_dos_filename, &
   & tdependence1_filename,tdependence2_filename,log_filename)
  use structure_module
  use min_images_module
  use force_constants_module
  implicit none
  
  ! ----------------------------------------
  ! Inputs
  ! ----------------------------------------
  type(StructureData),  intent(in) :: supercell
  type(ForceConstants), intent(in) :: force_constants
  real(dp),             intent(in) :: temperature
  
  ! ----------------------------------------
  ! filenames
  ! ----------------------------------------
  type(String), intent(in) :: free_energy_filename
  type(String), intent(in) :: freq_dos_filename
  type(String), intent(in) :: tdependence1_filename
  type(String), intent(in) :: tdependence2_filename
  type(String), intent(in) :: log_filename
  
  ! ----------------------------------------
  ! previously global variables
  ! ----------------------------------------
  real(dp) :: bin_width
  real(dp), allocatable :: freq_dos(:,:)
  
  type(MinImages), allocatable :: min_images(:,:)
  
  min_images = calculate_min_images(supercell)

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
  call generate_dos( supercell,            &
                   & min_images,           &
                   & force_constants,      &
                   & temperature,          &
                   & free_energy_filename, &
                   & freq_dos_filename,    &
                   & log_filename)
  call print_line('Done.  Frequency density-of-states function written to &
    &freq_dos.dat.  (Please view this file using XMGrace.)')
  call print_line('')

  call print_line('Calculating the lattice thermal energy (LTE) and free energy &
    &(LTFE)...')
  call calc_lte(bin_width,temperature,freq_dos,tdependence1_filename)
  call calc_ltfe(bin_width,temperature,freq_dos,tdependence2_filename)
  call print_line('')
end subroutine
end module
