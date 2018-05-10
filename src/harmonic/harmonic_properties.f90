! ======================================================================
! Routines to calculate a number of harmonic properties of the crystal given
!    the matrix of force constants.
! ======================================================================
module harmonic_properties_module
  use common_module
  
  use min_images_module
  use force_constants_module
  use dynamical_matrix_module
  use harmonic_thermodynamics_module
  implicit none
contains

! ----------------------------------------------------------------------
! Calculate the frequency density-of-states by random sampling of
!    the Brillouin zone.
! ----------------------------------------------------------------------
subroutine generate_dos(supercell,min_images,force_constants,       &
   & thermal_energies,min_frequency,no_dos_samples,energy_filename, &
   & dos_filename,logfile)
  implicit none
  
  type(StructureData),  intent(in)    :: supercell
  type(MinImages),      intent(in)    :: min_images(:,:)
  type(ForceConstants), intent(in)    :: force_constants
  real(dp),             intent(in)    :: thermal_energies(:)
  real(dp),             intent(in)    :: min_frequency
  integer,              intent(in)    :: no_dos_samples
  type(String),         intent(in)    :: energy_filename
  type(String),         intent(in)    :: dos_filename
  type(OFile),          intent(inout) :: logfile
  
  ! Calculation parameter.
  integer  :: no_bins
  integer  :: no_prelims
  integer  :: print_every
  real(dp) :: safety_margin
  
  ! Parameters calculated based on preliminary calculation.
  real(dp) :: min_freq
  real(dp) :: max_freq
  real(dp) :: freq_spread
  real(dp) :: bin_width
  
  ! Parameter for adding extra bins if needed.
  integer :: extra_bins
  
  ! Output variables.
  real(dp), allocatable :: freq_dos(:)
  real(dp), allocatable :: energy(:)
  real(dp), allocatable :: free_energy(:)
  real(dp), allocatable :: entropy(:)
  
  ! Working variables.
  real(dp)                                  :: qpoint(3)
  type(DynamicalMatrix)                     :: dyn_mat
  real(dp)                                  :: frequency
  integer                                   :: no_frequencies_ignored
  type(ThermodynamicVariables), allocatable :: thermodynamics(:)
  
  ! files.
  type(OFile) :: energy_file
  type(OFile) :: dos_file
  
  ! Temporary variables.
  integer :: bin
  integer :: i,j,ialloc
  
  ! Set parameters for calculation.
  no_bins       = no_dos_samples/100
  no_prelims    = no_dos_samples/10
  print_every   = no_dos_samples/10
  safety_margin = 0.15_dp
  
  ! Initialise the random number generator
  call random_seed()
  
  ! Establish (approximate) maximum and minimum frequencies and hence
  !    choose the bin width.
  max_freq = 0.0_dp
  min_freq = 0.0_dp
  do i=1,no_prelims
    call random_number(qpoint)
    dyn_mat = DynamicalMatrix( vec(qpoint),     &
                             & supercell,       &
                             & force_constants, &
                             & min_images)
    call dyn_mat%check( supercell, &
                      & logfile,   &
                      & check_eigenstuff=.false.)
    
    min_freq = min( min_freq, &
                  & dyn_mat%complex_modes(1)%frequency)
    max_freq = max( max_freq, &
                  & dyn_mat%complex_modes(supercell%no_modes_prim)%frequency)
  enddo
  
  if (max_freq<=min_frequency) then
    call print_line(ERROR//': The system is pathologically unstable; all &
       &frequencies are less than min_frequency.')
    stop
  endif
  
  ! Spread out min and max frequencies to leave safety margin.
  freq_spread =  max_freq - min_freq
  min_freq    =  min_freq - safety_margin*freq_spread
  max_freq    =  max_freq + safety_margin*freq_spread
  bin_width   = (max_freq - min_freq) / no_bins
  
  ! Calculate density of states.
  allocate( freq_dos(no_bins),                   &
          & energy(size(thermal_energies)),      &
          & free_energy(size(thermal_energies)), &
          & entropy(size(thermal_energies)),     &
          & stat=ialloc); call err(ialloc)
  freq_dos=0.0_dp
  energy = 0.0_dp
  free_energy = 0.0_dp
  entropy = 0.0_dp
  no_frequencies_ignored = 0
  do i=1,no_dos_samples
    call random_number(qpoint)
    dyn_mat = DynamicalMatrix( vec(qpoint),     &
                             & supercell,       &
                             & force_constants, &
                             & min_images)
    call dyn_mat%check( supercell,         &
                      & logfile,           &
                      & check_eigenstuff=.false.)
    
    do j=1,supercell%no_modes_prim
      frequency = dyn_mat%complex_modes(j)%frequency
      
      ! Bin frequency for density of states.
      bin = ceiling( (frequency-min_freq) / bin_width)
      if (bin<1) then
        extra_bins  = 1-bin
        min_freq    = min_freq    - extra_bins*bin_width
        freq_spread = freq_spread + extra_bins*bin_width
        freq_dos    = [dble(dblevec(zeroes(extra_bins))), freq_dos]
        no_bins     = no_bins + extra_bins
        bin         = bin     + extra_bins
      elseif (bin>no_bins) then
        extra_bins  = bin-no_bins
        max_freq    = max_freq    + extra_bins*bin_width
        freq_spread = freq_spread + extra_bins*bin_width
        freq_dos    = [freq_dos, dble(dblevec(zeroes(extra_bins)))]
        no_bins     = no_bins + extra_bins
      endif
      freq_dos(bin) = freq_dos(bin) + 1
      
      ! Calculate thermodynamic quantities.
      if (frequency<min_frequency) then
        no_frequencies_ignored = no_frequencies_ignored + 1
      else
        thermodynamics = ThermodynamicVariables(thermal_energies,frequency)
        energy = energy + thermodynamics%energy
        free_energy = free_energy + thermodynamics%free_energy
        entropy = entropy + thermodynamics%entropy
      endif
    enddo
    
    if (modulo(i,print_every)==0) then
      call print_line('Sampling q-points: '//i//' of '//no_dos_samples// &
         & ' samples complete.')
    endif
  enddo
  
  ! Normalise variables to be per unit cell.
  freq_dos    = freq_dos    / (no_dos_samples*bin_width)
  energy      = energy      / no_dos_samples
  free_energy = free_energy / no_dos_samples
  entropy     = entropy     / no_dos_samples
  
  call print_line(WARNING//': '//no_frequencies_ignored//' modes ignored for &
     &having frequencies less than min_frequency, out of '// &
     & no_dos_samples*supercell%no_modes_prim//' modes sampled.')
  
  ! Write out density of states.
  dos_file = OFile(dos_filename)
  do i=1,no_bins
    call dos_file%print_line(                &
       & 'Bin: '//min_freq+bin_width*(i-1)// &
       & ' to '//min_freq+bin_width*i)
    call dos_file%print_line( 'Bin DOS: '//freq_dos(i))
    call dos_file%print_line( '')
  enddo
  
  ! Write out thermodynamic variables.
  energy_file = OFile(energy_filename)
  call energy_file%print_line( 'kB * temperature (Hartree) | &
                              &U=<E> (Hartree) | &
                              &F=U-TS (Hartree) | &
                              &S (dimensionless)')
  do i=1,size(thermal_energies)
    call energy_file%print_line( thermal_energies(i) //' '// &
                               & energy(i)           //' '// &
                               & free_energy(i)      //' '// &
                               & entropy(i))
  enddo
end subroutine

! ----------------------------------------------------------------------
! Generates the phonon dispersion curve.
! ----------------------------------------------------------------------
subroutine generate_dispersion(large_supercell,min_images,force_constants, &
   & path_labels,path_qpoints,dispersion_filename,                         &
   & high_symmetry_points_filename,logfile)
  implicit none
  
  type(StructureData),  intent(in)    :: large_supercell
  type(MinImages),      intent(in)    :: min_images(:,:)
  type(ForceConstants), intent(in)    :: force_constants
  type(String),         intent(in)    :: path_labels(:)
  type(RealVector),     intent(in)    :: path_qpoints(:)
  type(String),         intent(in)    :: dispersion_filename
  type(String),         intent(in)    :: high_symmetry_points_filename
  type(OFile),          intent(inout) :: logfile
  
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
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
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
  high_symmetry_points_file = OFile(high_symmetry_points_filename)
  do i=1,no_vertices
    call high_symmetry_points_file%print_line( &
       & 'q-point: '//path_labels(i)//' '//path_qpoints(i))
    call high_symmetry_points_file%print_line( &
       & 'Fraction along path: '//fractional_distances(i))
    call high_symmetry_points_file%print_line('')
  enddo
  
  ! Travel along q-space paths, calculating frequencies at each point.
  dispersion_file = OFile(dispersion_filename)
  do i=1,no_segments
    do j=0,points_per_segment(i)-1
      qpoint = ( (points_per_segment(i)-j)*path_qpoints(i)     &
           &   + j                        *path_qpoints(i+1) ) &
           & / points_per_segment(i)
      dyn_mat = DynamicalMatrix( qpoint,          &
                               & large_supercell, &
                               & force_constants, &
                               & min_images)
      call dyn_mat%check( large_supercell, &
                        & logfile,         &
                        & check_eigenstuff=.false.)
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
                           & large_supercell, &
                           & force_constants, &
                           & min_images)
  call dyn_mat%check( large_supercell, &
                    & logfile,         &
                    & check_eigenstuff=.false.)
  call dispersion_file%print_line( &
     & 'Fraction along path: '//1.0_dp)
  call dispersion_file%print_line( &
     & 'Frequencies: '//dyn_mat%frequencies())
end subroutine
end module