submodule (caesar_harmonic_properties_module) caesar_harmonic_properties_submodule
  use caesar_observables_module
contains

module procedure new_PathFrequencies
  this%path_fraction = path_fraction
  this%qpoint = qpoint
  this%frequencies = frequencies
  this%dynamical_matrix = dynamical_matrix
end procedure

module procedure new_SampledQpoint
  this%qpoint = qpoint
  this%no_soft_frequencies = no_soft_frequencies
end procedure

module procedure new_PdosBin
  this%min_frequency = min_frequency
  this%max_frequency = max_frequency
  this%occupation = occupation
end procedure

module procedure new_PhononDispersion
  this%path_ = path
  this%frequencies_ = frequencies
end procedure

module procedure new_PhononDos
  this%qpoints = qpoints
  this%pdos = pdos
  this%thermodynamic_data = thermodynamic_data
end procedure

module procedure new_PhononDispersion_CartesianHessian
  ! q-point and dynamical matrix variables.
  type(QpointPath)                         :: path
  type(PathQpointAndFraction), allocatable :: path_qpoints(:)
  type(RealVector)                         :: qpoint
  type(DynamicalMatrix)                    :: dynamical_matrix
  type(ComplexMode),           allocatable :: complex_modes(:)
  
  ! Output variables.
  type(PathFrequencies), allocatable :: frequencies(:)
  
  ! Temporary variables.
  integer :: print_every
  integer :: i,ialloc
  
  ! Generate q-point path.
  path = QpointPath(path_string, no_path_points)
  
  ! Travel along q-space paths, calculating frequencies at each point.
  path_qpoints = path%path_qpoints()
  print_every = size(path_qpoints)/10
  allocate(frequencies(size(path_qpoints)), stat=ialloc); call err(ialloc)
  do i=1,size(path_qpoints)
    qpoint = path_qpoints(i)%qpoint
    dynamical_matrix = DynamicalMatrix( qpoint,    &
                                      & supercell, &
                                      & hessian,   &
                                      & min_images )
    call dynamical_matrix%check(supercell, logfile)
    complex_modes = ComplexMode(dynamical_matrix, supercell)
    frequencies(i) = PathFrequencies(                       &
       &  path_fraction    = path_qpoints(i)%path_fraction, &
       &  qpoint           = qpoint,                        &
       &  frequencies      = complex_modes%frequency,       &
       &  dynamical_matrix = dynamical_matrix               )
    
    if (modulo(i,print_every)==0) then
      call print_line('Dispersion curve: '//i//' of '//size(path_qpoints)// &
         & ' q-points sampled.')
    endif
  enddo
  
  ! Write outputs to file.
  this = PhononDispersion(path, frequencies)
end procedure

module procedure path_PhononDispersion
  output = str(this%path_%vertices, separating_line='')
end procedure

module procedure frequencies_PhononDispersion
  output = str(this%frequencies_)
end procedure

module procedure json_PhononDispersion
  type(RealVector) :: vector
  
  type(String)              :: lattice(3)
  type(String)              :: atom_types
  type(String)              :: atom_numbers
  type(String), allocatable :: atom_pos_car(:)
  type(String), allocatable :: atom_pos_red(:)
  type(String), allocatable :: highsym_qpts(:)
  type(String), allocatable :: qpoints(:)
  type(String), allocatable :: distances(:)
  type(String), allocatable :: eigenvalues(:)
  type(String), allocatable :: vectors(:)
  
  type(ComplexMode), allocatable :: modes(:)
  
  complex(dp) :: element
  
  integer :: i,j,k,l,n,ialloc
  
  ! Construct lattice strings.
  do i=1,3
    lattice(i) = '    ['
    do j=1,3
      lattice(i) = lattice(i)//structure%lattice%element(i,j)
      if (j/=3) then
        lattice(i) = lattice(i)//','
      endif
    enddo
    lattice(i) = lattice(i)//']'
    if (i/=3) then
      lattice(i) = lattice(i)//','
    endif
  enddo
  
  ! Construct atom_types string.
  atom_types = '['
  do i=1,size(structure%atoms)
    atom_types = atom_types//'"'//structure%atoms(i)%species()//'"'
    if (i/=size(structure%atoms)) then
      atom_types = atom_types//','
    endif
  enddo
  atom_types = atom_types//']'
  
  ! Construct atom_numbers string.
  atom_numbers = '['
  do i=1,size(structure%atoms)
    atom_numbers = atom_numbers//atomic_symbol_to_number( &
                           & structure%atoms(i)%species() )
    if (i/=size(structure%atoms)) then
      atom_numbers = atom_numbers//','
    endif
  enddo
  atom_numbers = atom_numbers//']'
  
  ! Construct atom_pos_car strings.
  allocate(atom_pos_car(size(structure%atoms)), stat=ialloc); call err(ialloc)
  do i=1,size(structure%atoms)
    atom_pos_car(i) = '    ['
    vector = structure%atoms(i)%cartesian_position()
    do j=1,3
      atom_pos_car(i) = atom_pos_car(i)//vector%element(j)
      if (j/=3) then
        atom_pos_car(i) = atom_pos_car(i)//','
      endif
    enddo
    atom_pos_car(i) = atom_pos_car(i)//']'
    if (i/=size(structure%atoms)) then
      atom_pos_car(i) = atom_pos_car(i)//','
    endif
  enddo
  
  ! Construct atom_pos_red strings.
  allocate(atom_pos_red(size(structure%atoms)), stat=ialloc); call err(ialloc)
  do i=1,size(structure%atoms)
    atom_pos_red(i) = '    ['
    vector = structure%atoms(i)%fractional_position()
    do j=1,3
      atom_pos_red(i) = atom_pos_red(i)//vector%element(j)
      if (j/=3) then
        atom_pos_red(i) = atom_pos_red(i)//','
      endif
    enddo
    atom_pos_red(i) = atom_pos_red(i)//']'
    if (i/=size(structure%atoms)) then
      atom_pos_red(i) = atom_pos_red(i)//','
    endif
  enddo
  
  ! Construct highsym_qpts.
  allocate( highsym_qpts(size(this%path_%vertices)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(this%path_%vertices)
    highsym_qpts(i) = '    ['                            // &
                    & this%path_%vertices(i)%point_index // &
                    & ',"'                               // &
                    & this%path_%vertices(i)%label       // &
                    & '"]'
    if (i/=size(this%path_%vertices)) then
      highsym_qpts(i) = highsym_qpts(i)//','
    endif
  enddo
  
  ! Construct qpoints.
  allocate(qpoints(size(this%frequencies_)), stat=ialloc); call err(ialloc)
  do i=1,size(this%frequencies_)
    qpoints(i) = '    ['
    do j=1,3
      qpoints(i) = qpoints(i)//this%frequencies_(i)%qpoint%element(j)
      if (j/=3) then
        qpoints(i) = qpoints(i)//','
      endif
    enddo
    qpoints(i) = qpoints(i)//']'
    if (i/=size(this%frequencies_)) then
      qpoints(i) = qpoints(i)//','
    endif
  enddo
  
  ! Construct distances.
  allocate(distances(size(this%frequencies_)), stat=ialloc); call err(ialloc)
  do i=1,size(this%frequencies_)
    distances(i) = '  '//this%frequencies_(i)%path_fraction
    if (i/=size(this%frequencies_)) then
      distances(i) = distances(i)//','
    endif
  enddo
  
  ! Construct eigenvalues strings.
  allocate(eigenvalues(size(this%frequencies_)), stat=ialloc); call err(ialloc)
  do i=1,size(this%frequencies_)
    eigenvalues(i) = '    ['
    do j=1,size(this%frequencies_(i)%frequencies)
      eigenvalues(i) = eigenvalues(i)//this%frequencies_(i)%frequencies(j)
      if (j/=size(this%frequencies_(i)%frequencies)) then
        eigenvalues(i) = eigenvalues(i)//','
      endif
    enddo
    eigenvalues(i) = eigenvalues(i)//']'
    if (i/=size(this%frequencies_)) then
      eigenvalues(i) = eigenvalues(i)//','
    endif
  enddo
  
  ! Construct eigenvectors strings.
  allocate( vectors(size(this%frequencies_)*(structure%no_modes+2)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(this%frequencies_)
    modes = ComplexMode(this%frequencies_(i)%dynamical_matrix, structure)
    n = (i-1)*(structure%no_modes+2) + 1
    vectors(n) = '    ['
    do j=1,structure%no_modes
      n = (i-1)*(structure%no_modes+2) + j + 1
      vectors(n) = '      ['
      do k=1,structure%no_atoms
        vectors(n) = vectors(n)//'['
        do l=1,3
          element = modes(j)%unit_vector(k)%element(l)
          vectors(n) = vectors(n)//'['//real(element)//','//aimag(element)//']'
          if (l/=3) then
            vectors(n) = vectors(n)//','
          endif
        enddo
        vectors(n) = vectors(n)//']'
        if (k/=structure%no_atoms) then
          vectors(n) = vectors(n)//','
        endif
      enddo
      vectors(n) = vectors(n)//']'
      if (j/=structure%no_modes) then
        vectors(n) = vectors(n)//','
      endif
    enddo
    n = i*(structure%no_modes+2)
    vectors(n) = '    ]'
    if (i/=size(this%frequencies_)) then
      vectors(n) = vectors(n)//','
    endif
  enddo
  
  ! Construct output.
  output = [ str('{'),                                   &
           & '  "name": "'//seedname//'",',              &
           & '  "natoms": '//size(structure%atoms)//',', &
           & str('  "lattice": ['),                      &
           & lattice,                                    &
           & str('  ],'),                                &
           & '  "atom_types": '//atom_types//',',        &
           & '  "atom_numbers": '//atom_numbers//',',    &
           & '  "formula": "'//seedname//'",',           &
           & str('  "repetitions": [1,1,1],'),           &
           & str('  "atom_pos_car": ['),                 &
           & atom_pos_car,                               &
           & str('  ],'),                                &
           & str('  "atom_pos_red": ['),                 &
           & atom_pos_red,                               &
           & str('  ],'),                                &
           & str('  "highsym_qpts": ['),                 &
           & highsym_qpts,                               &
           & str('  ],'),                                &
           & str('  "qpoints": ['),                      &
           & qpoints,                                    &
           & str('  ],'),                                &
           & str('  "distances": ['),                    &
           & distances,                                  &
           & str('  ],'),                                &
           & str('  "eigenvalues": ['),                  &
           & eigenvalues,                                &
           & str('  ],'),                                &
           & str('  "vectors": ['),                      &
           & vectors,                                    &
           & str('  ]'),                                 &
           & str('}')                                    ]
end procedure

module procedure new_PhononDos_CartesianHessian
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
  real(dp),         allocatable :: freq_dos(:)
  real(dp),         allocatable :: energy(:)
  real(dp),         allocatable :: free_energy(:)
  real(dp),         allocatable :: entropy(:)
  type(RealMatrix), allocatable :: stress(:)
  real(dp),         allocatable :: volume
  real(dp),         allocatable :: enthalpy(:)
  real(dp),         allocatable :: gibbs(:)
  
  type(SampledQpoint),     allocatable :: qpoints(:)
  type(ThermodynamicData), allocatable :: thermodynamic_data(:)
  type(PdosBin),           allocatable :: pdos(:)
  
  ! Working variables.
  type(RealVector)                     :: qpoint
  type(DynamicalMatrix)                :: dynamical_matrix
  real(dp)                             :: frequency
  type(ComplexMode),       allocatable :: complex_modes(:)
  integer,                 allocatable :: no_frequencies_ignored(:)
  type(ThermodynamicData), allocatable :: thermodynamics(:)
  
  ! Stress variables.
  type(StressDynamicalMatrix) :: stress_dynamical_matrix
  type(RealMatrix)            :: stress_prefactor
  
  ! Temporary variables.
  integer :: bin
  integer :: i,j,k,ialloc
  
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
    qpoint = vec(random_generator%random_numbers(3))
    dynamical_matrix = DynamicalMatrix( qpoint,    &
                                      & supercell, &
                                      & hessian,   &
                                      & min_images )
    call dynamical_matrix%check(supercell, logfile)
    
    complex_modes = ComplexMode(dynamical_matrix, supercell)
    min_freq = min(min_freq, complex_modes(1)%frequency)
    max_freq = max(max_freq, complex_modes(size(complex_modes))%frequency)
  enddo
  
  if (max_freq<=min_frequency) then
    call print_line(ERROR//': The system is pathologically unstable; all &
       &frequencies are less than min_frequency.')
    call quit()
  endif
  
  ! Spread out min and max frequencies to leave safety margin.
  freq_spread =  max_freq - min_freq
  min_freq    =  min_freq - safety_margin*freq_spread
  max_freq    =  max_freq + safety_margin*freq_spread
  bin_width   = (max_freq - min_freq) / no_bins
  
  ! Initialise containers.
  freq_dos = [(0.0_dp, i=1, no_bins)]
  energy = [(0.0_dp, i=1, size(thermal_energies))]
  free_energy = [(0.0_dp, i=1, size(thermal_energies))]
  entropy = [(0.0_dp, i=1, size(thermal_energies))]
  no_frequencies_ignored = [(0, i=1, no_dos_samples)]
  if (present(stress_hessian)) then
    stress = [(dblemat(zeroes(3,3)), i=1, size(thermal_energies))]
    volume = supercell%volume / supercell%sc_size
  endif
  allocate(qpoints(no_dos_samples), stat=ialloc); call err(ialloc)
  
  ! Calculate density of states.
  do i=1,no_dos_samples
    qpoint = vec(random_generator%random_numbers(3))
    dynamical_matrix = DynamicalMatrix( qpoint,    &
                                      & supercell, &
                                      & hessian,   &
                                      & min_images )
    call dynamical_matrix%check(supercell, logfile)
    complex_modes = ComplexMode(dynamical_matrix, supercell)
    
    if (present(stress_hessian)) then
      stress_dynamical_matrix = StressDynamicalMatrix( qpoint,           &
                                                     & stress_supercell, &
                                                     & stress_hessian,   &
                                                     & stress_min_images )
    endif
    
    do j=1,supercell%no_modes_prim
      frequency = complex_modes(j)%frequency
      
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
        thermodynamics = [( ThermodynamicData( thermal_energies(k),    &
                          &                    frequency            ), &
                          & k=1,                                       &
                          & size(thermal_energies)                     )]
        energy = energy + thermodynamics%energy
        free_energy = free_energy + thermodynamics%free_energy
        entropy = entropy + thermodynamics%entropy
        
        if (present(stress_hessian)) then
          ! Add in the kinetic stress from the mode.
          ! For harmonic states, the kinetic stress is equal to UI/V,
          !    where I is the stress prefactor
          !    and V is the primitive cell volume.
          stress_prefactor = complex_modes(j)%stress_prefactor()
          stress = stress + thermodynamics%energy * stress_prefactor/volume
          
          ! Add in the potential stress from the mode.
          ! For harmonic states, U=2<V> = w^2<uu*>.
          ! <stress> = -1/2 * u*.stress_dyn_mat.u * <uu*>.
          !          = -1/2 * u*.stress_dyn_mat.u * U/w^2.
          stress = stress                                                &
               & - stress_dynamical_matrix%expectation(complex_modes(j)) &
               & * thermodynamics%energy / (2*complex_modes(j)%frequency**2)
        endif
      endif
    enddo
    
    qpoints(i) = SampledQpoint(qpoint, no_frequencies_ignored(i))
    
    ! Print progress.
    if (modulo(i,print_every)==0) then
      call print_line('Density of states: '//i//' of '//no_dos_samples// &
         & ' q-points sampled.')
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
  if (present(stress_hessian)) then
    stress = stress / no_dos_samples
    enthalpy = energy + trace(stress) * volume / 3
    gibbs = free_energy + trace(stress) * volume / 3
  endif
  
  if (any(no_frequencies_ignored/=0)) then
    call print_line(WARNING//': '//sum(no_frequencies_ignored)//' modes &
       &ignored for having frequencies less than min_frequency, out of '// &
       & no_dos_samples*supercell%no_modes_prim//' modes sampled.')
  endif
  
  pdos = [( PdosBin( min_frequency = min_freq+(i-1)*bin_width,    &
          &          max_frequency = min_freq+i*bin_width,        &
          &          occupation    = freq_dos(i)               ), &
          & i=1,                                                  &
          & no_bins                                               )]
  
  if (present(stress_hessian)) then
    thermodynamic_data = ThermodynamicData( thermal_energies, &
                                          & energy,           &
                                          & free_energy,      &
                                          & entropy,          &
                                          & stress,           &
                                          & volume,           &
                                          & enthalpy,         &
                                          & gibbs             )
  else
    thermodynamic_data = ThermodynamicData( thermal_energies, &
                                          & energy,           &
                                          & free_energy,      &
                                          & entropy           )
  endif
  
  this = PhononDos(qpoints, pdos, thermodynamic_data)
end procedure

module procedure write_PathFrequencies
  select type(this); type is(PathFrequencies)
    output = [ 'Fraction along path: '//this%path_fraction, &
             & 'Frequencies: '//this%frequencies            ]
  class default
    call err()
  end select
end procedure

module procedure read_SampledQpoint
  type(RealVector) :: qpoint
  integer          :: no_soft_frequencies
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(SampledQpoint)
    line = split_line(input)
    qpoint = vec(dble(line(:3)))
    no_soft_frequencies = int(line(4))
    
    this = SampledQpoint(qpoint,no_soft_frequencies)
  class default
    call err()
  end select
end procedure

module procedure write_SampledQpoint
  select type(this); type is(SampledQpoint)
    output = this%qpoint//' '//this%no_soft_frequencies
  class default
    call err()
  end select
end procedure

module procedure new_SampledQpoint_String
  call this%read(input)
end procedure

module procedure read_PdosBin
  real(dp) :: min_frequency
  real(dp) :: max_frequency
  real(dp) :: occupation
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(PdosBin)
    line = split_line(input)
    min_frequency = dble(line(3))
    max_frequency = dble(line(5))
    occupation = dble(line(7))
    
    this = PdosBin(min_frequency,max_frequency,occupation)
  class default
    call err()
  end select
end procedure

module procedure write_PdosBin
  select type(this); type is(PdosBin)
    output = 'Bin frequencies: '//this%min_frequency//' to &
             &'//this%max_frequency//' PDOS: '//this%occupation
  class default
    call err()
  end select
end procedure

module procedure new_PdosBin_String
  call this%read(input)
end procedure

module procedure write_files_PhononDispersion
  type(OFile) :: symmetry_points_file
  type(OFile) :: dispersion_file
  type(OFile) :: json_file
  
  symmetry_points_file = OFile(directory//'/high_symmetry_points.dat')
  call symmetry_points_file%print_lines(this%path())
  
  dispersion_file = OFile(directory//'/phonon_dispersion_curve.dat')
  call dispersion_file%print_lines(this%frequencies())
  
  json_file = OFile(directory//'/'//seedname//'.json')
  call json_file%print_lines(this%json(seedname,structure))
end procedure

module procedure write_files_PhononDos
  type(OFile) :: sampled_qpoints_file
  type(OFile) :: pdos_file
  
  sampled_qpoints_file = OFile(directory//'/sampled_qpoints.dat')
  call sampled_qpoints_file%print_line('q-point (x,y,z) | &
                                       &number of frequencies ignored')
  call sampled_qpoints_file%print_lines(this%qpoints)
  
  pdos_file = OFile(directory//'/phonon_density_of_states.dat')
  call pdos_file%print_lines(this%pdos)
end procedure
end submodule
