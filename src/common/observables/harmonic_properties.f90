! ======================================================================
! Routines to calculate a number of harmonic properties of the crystal given
!    the matrix of force constants.
! ======================================================================
module harmonic_properties_module
  use utils_module
  
  use structure_module
  use dynamical_matrices_module
  
  use harmonic_thermodynamics_module
  implicit none
  
  ! The major points along the q-point path.
  type, extends(Stringsable) :: QpointPath
    type(String)     :: label
    type(RealVector) :: qpoint
    real(dp)         :: path_fraction
  contains
    procedure, public :: read  => read_QpointPath
    procedure, public :: write => write_QpointPath
  end type
  
  interface QpointPath
    module procedure new_QpointPath
    module procedure new_QpointPath_Strings
    module procedure new_QpointPath_StringArray
  end interface
  
  ! The minor points along the q-point path.
  type, extends(Stringsable) :: PathFrequencies
    real(dp)              :: path_fraction
    real(dp), allocatable :: frequencies(:)
  contains
    procedure, public :: read  => read_PathFrequencies
    procedure, public :: write => write_PathFrequencies
  end type
  
  interface PathFrequencies
    module procedure new_PathFrequencies
    module procedure new_PathFrequencies_Strings
    module procedure new_PathFrequencies_StringArray
  end interface
contains

! Constructors.
impure elemental function new_QpointPath(label,qpoint,path_fraction) &
   & result(this)
  implicit none
  
  type(String),     intent(in) :: label
  type(RealVector), intent(in) :: qpoint
  real(dp),         intent(in) :: path_fraction
  type(QpointPath)             :: this
  
  this%label = label
  this%qpoint = qpoint
  this%path_fraction = path_fraction
end function

function new_PathFrequencies(path_fraction,frequencies) result(this)
  implicit none
  
  real(dp), intent(in)  :: path_fraction
  real(dp), intent(in)  :: frequencies(:)
  type(PathFrequencies) :: this
  
  this%path_fraction = path_fraction
  this%frequencies = frequencies
end function

! ----------------------------------------------------------------------
! Calculate the frequency density-of-states by random sampling of
!    the Brillouin zone.
! ----------------------------------------------------------------------
subroutine generate_dos(supercell,min_images,force_constants,            &
   & thermal_energies,min_frequency,no_dos_samples,sampled_qpoints_file, &
   & thermodynamic_file,pdos_file,logfile,random_generator)
  implicit none
  
  type(StructureData),  intent(in)    :: supercell
  type(MinImages),      intent(in)    :: min_images(:,:)
  type(ForceConstants), intent(in)    :: force_constants
  real(dp),             intent(in)    :: thermal_energies(:)
  real(dp),             intent(in)    :: min_frequency
  integer,              intent(in)    :: no_dos_samples
  type(OFile),          intent(inout) :: sampled_qpoints_file
  type(OFile),          intent(inout) :: thermodynamic_file
  type(OFile),          intent(inout) :: pdos_file
  type(OFile),          intent(inout) :: logfile
  type(RandomReal),     intent(in)    :: random_generator
  
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
  type(RealVector)                          :: qpoint
  type(DynamicalMatrix)                     :: dyn_mat
  real(dp)                                  :: frequency
  integer,                      allocatable :: no_frequencies_ignored(:)
  type(ThermodynamicVariables), allocatable :: thermodynamics(:)
  
  ! Temporary variables.
  integer :: bin
  integer :: i,j,ialloc
  
  ! Set parameters for calculation.
  no_bins       = no_dos_samples/100
  no_prelims    = no_dos_samples/10
  print_every   = no_dos_samples/10
  safety_margin = 0.15_dp
  
  ! Establish (approximate) maximum and minimum frequencies and hence
  !    choose the bin width.
  max_freq = 0.0_dp
  min_freq = 0.0_dp
  do i=1,no_prelims
    qpoint = random_generator%random_numbers(3)
    dyn_mat = DynamicalMatrix( qpoint,          &
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
  allocate( freq_dos(no_bins),                      &
          & energy(size(thermal_energies)),         &
          & free_energy(size(thermal_energies)),    &
          & entropy(size(thermal_energies)),        &
          & no_frequencies_ignored(no_dos_samples), &
          & stat=ialloc); call err(ialloc)
  freq_dos=0.0_dp
  energy = 0.0_dp
  free_energy = 0.0_dp
  entropy = 0.0_dp
  no_frequencies_ignored = 0
  call sampled_qpoints_file%print_line('q-point (x,y,z) | &
                                       &number of frequencies ignored')
  do i=1,no_dos_samples
    qpoint = random_generator%random_numbers(3)
    dyn_mat = DynamicalMatrix( qpoint,          &
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
        no_frequencies_ignored(i) = no_frequencies_ignored(i) + 1
      else
        thermodynamics = ThermodynamicVariables(thermal_energies,frequency)
        energy = energy + thermodynamics%energy
        free_energy = free_energy + thermodynamics%free_energy
        entropy = entropy + thermodynamics%entropy
      endif
    enddo
    
    call sampled_qpoints_file%print_line( qpoint//' '// &
                                        & no_frequencies_ignored(i))
    
    if (modulo(i,print_every)==0) then
      call print_line('Sampling q-points: '//i//' of '//no_dos_samples// &
         & ' samples complete.')
    endif
  enddo
  
  ! Normalise variables to be per unit cell.
  ! N.B. the divisor is not corrected for ignored frequencies, since ignored
  !    frequencies are considered to contribute zero energy, F and S.
  ! (The contribution of a single low-frequency mode diverges, but assuming
  !    that such modes are localised around Gamma with a typical phonon
  !    dispersion then the integral across them approaches zero.)
  freq_dos    = freq_dos    / (no_dos_samples*bin_width)
  energy      = energy      / no_dos_samples
  free_energy = free_energy / no_dos_samples
  entropy     = entropy     / no_dos_samples
  
  if (any(no_frequencies_ignored/=0)) then
    call print_line(WARNING//': '//sum(no_frequencies_ignored)//' modes &
       &ignored for having frequencies less than min_frequency, out of '// &
       & no_dos_samples*supercell%no_modes_prim//' modes sampled.')
  endif
  
  ! Write out phonon density of states.
  do i=1,no_bins
    call pdos_file%print_line(                &
        & 'Bin: '//min_freq+bin_width*(i-1)// &
        & ' to '//min_freq+bin_width*i)
    call pdos_file%print_line( 'Bin DOS: '//freq_dos(i))
    call pdos_file%print_line( '')
  enddo
  
  ! Write out thermodynamic variables.
  call thermodynamic_file%print_line( &
     &'kB * temperature (Hartree per cell) | &
     &Vibrational Energy per cell, U=<E>, (Hartree) | &
     &Vibrational Free Energy per cell, F=U-TS, (Hartree) | &
     &Vibrational Shannon Entropy per cell, S/k_B, (arb. units)')
  do i=1,size(thermal_energies)
    call thermodynamic_file%print_line( thermal_energies(i) //' '// &
                                      & energy(i)           //' '// &
                                      & free_energy(i)      //' '// &
                                      & entropy(i))
  enddo
end subroutine

! ----------------------------------------------------------------------
! Generates the phonon dispersion curve.
! ----------------------------------------------------------------------
subroutine generate_dispersion(large_supercell,min_images,force_constants, &
   & path_labels,path_qpoints,dispersion_file,symmetry_points_file,logfile)
  implicit none
  
  type(StructureData),  intent(in)    :: large_supercell
  type(MinImages),      intent(in)    :: min_images(:,:)
  type(ForceConstants), intent(in)    :: force_constants
  type(String),         intent(in)    :: path_labels(:)
  type(RealVector),     intent(in)    :: path_qpoints(:)
  type(OFile),          intent(inout) :: dispersion_file
  type(OFile),          intent(inout) :: symmetry_points_file
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
  
  ! Output variables.
  type(QpointPath),      allocatable :: path(:)
  type(PathFrequencies), allocatable :: frequencies(:)
  
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
  
  ! Store q-point path.
  path = QpointPath(path_labels, path_qpoints, fractional_distances)
  
  ! Travel along q-space paths, calculating frequencies at each point.
  frequencies = [PathFrequencies::]
  do i=1,no_segments
    do j=0,points_per_segment(i)-1
      qpoint = ( (points_per_segment(i)-j)*path_qpoints(i)     &
           &   + j                        *path_qpoints(i+1) ) &
           & / points_per_segment(i)
      dyn_mat = DynamicalMatrix( qpoint,          &
                               & large_supercell, &
                               & force_constants, &
                               & min_images       )
      call dyn_mat%check( large_supercell,         &
                        & logfile,                 &
                        & check_eigenstuff=.false. )
      frequencies = [                                                         &
         & frequencies,                                                       &
         & PathFrequencies(                                                   &
         &    path_fraction=fractional_distances(i)+j*fractional_separation,  &
         &    frequencies=dyn_mat%complex_modes%frequency                    )]
    enddo
  enddo
  
  ! Calculate frequencies at final k-space point.
  qpoint = path_qpoints(no_vertices)
  dyn_mat = DynamicalMatrix( qpoint,          &
                           & large_supercell, &
                           & force_constants, &
                           & min_images       )
  call dyn_mat%check( large_supercell,         &
                    & logfile,                 &
                    & check_eigenstuff=.false. )
  frequencies = [                                     &
     & frequencies,                                   &
     & PathFrequencies(                               &
     &    path_fraction=1.0_dp,                       &
     &    frequencies=dyn_mat%complex_modes%frequency )]
  
  ! Write outputs to file.
  call symmetry_points_file%print_lines(path, separating_line='')
  call dispersion_file%print_lines(frequencies)
end subroutine

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_QpointPath(this,input)
  implicit none
  
  class(QpointPath), intent(out) :: this
  type(String),      intent(in)  :: input(:)
  
  type(String)     :: label
  type(RealVector) :: qpoint
  real(dp)         :: path_fraction
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(QpointPath)
    line = split_line(input(1))
    label = line(2)
    qpoint = vec(dble(line(3:)))
    
    line = split_line(input(2))
    path_fraction = dble(line(4))
    
    this = QpointPath(label,qpoint,path_fraction)
  class default
    call err()
  end select
end subroutine

function write_QpointPath(this) result(output)
  implicit none
  
  class(QpointPath), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  select type(this); type is(QpointPath)
    output = [ 'q-point: '//this%label//' '//this%qpoint,  &
             & 'Fraction along path: '//this%path_fraction ]
  class default
    call err()
  end select
end function

function new_QpointPath_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(QpointPath)         :: this
  
  call this%read(input)
end function

impure elemental function new_QpointPath_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(QpointPath)              :: this
  
  this = QpointPath(str(input))
end function

subroutine read_PathFrequencies(this,input)
  implicit none
  
  class(PathFrequencies), intent(out) :: this
  type(String),      intent(in)  :: input(:)
  
  real(dp)              :: path_fraction
  real(dp), allocatable :: frequencies(:)
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(PathFrequencies)
    line = split_line(input(1))
    path_fraction = dble(line(4))
    
    line = split_line(input(2))
    frequencies = dble(line(2:))
    
    this = PathFrequencies(path_fraction,frequencies)
  class default
    call err()
  end select
end subroutine

function write_PathFrequencies(this) result(output)
  implicit none
  
  class(PathFrequencies), intent(in) :: this
  type(String), allocatable          :: output(:)
  
  select type(this); type is(PathFrequencies)
    output = [ 'Fraction along path: '//this%path_fraction, &
             & 'Frequencies: '//this%frequencies            ]
  class default
    call err()
  end select
end function

function new_PathFrequencies_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(PathFrequencies)    :: this
  
  call this%read(input)
end function

impure elemental function new_PathFrequencies_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(PathFrequencies)         :: this
  
  this = PathFrequencies(str(input))
end function
end module
