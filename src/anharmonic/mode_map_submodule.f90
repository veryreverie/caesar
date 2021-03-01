submodule (caesar_mode_map_module) caesar_mode_map_submodule
  use caesar_anharmonic_module
contains

module procedure new_ModeMap
  this%mode_id                    = mode_id
  this%harmonic_frequency         = harmonic_frequency
  this%mode_displacements         = mode_displacements
  this%l2_cartesian_displacements = l2_cartesian_displacements
  this%harmonic_energies          = harmonic_energies
  this%harmonic_forces            = harmonic_forces
  
  if (present(anharmonic_energies)) then
    this%anharmonic_energies = anharmonic_energies
  endif
  if (present(anharmonic_forces)) then
    this%anharmonic_forces = anharmonic_forces
  endif
  if (present(anharmonic_pressures)) then
    this%anharmonic_pressures = anharmonic_pressures
  endif
  if (present(sampled_energies)) then
    this%sampled_energies = sampled_energies
  endif
  if (present(sampled_forces)) then
    this%sampled_forces = sampled_forces
  endif
  if (present(sampled_pressures)) then
    this%sampled_pressures = sampled_pressures
  endif
end procedure

module procedure new_ModeMap_potential
  ! Output variables.
  real(dp), allocatable :: harmonic_energies(:)
  real(dp), allocatable :: harmonic_forces(:)
  real(dp), allocatable :: anharmonic_energies(:)
  real(dp), allocatable :: anharmonic_forces(:)
  real(dp), allocatable :: anharmonic_pressures(:)
  
  ! Force in real mode co-ordinates.
  type(RealModeDisplacement) :: displacement
  type(RealModeForce)        :: anharmonic_force
  
  ! Temporary variables.
  integer :: i,ialloc
  
  if (present(potential) .and. .not. present(anharmonic_data)) then
    call print_line(CODE_ERROR//': if potential is present then &
       &anharmonic_data must also be present.')
    call err()
  elseif (present(stress) .and. .not. present(potential)) then
    call print_line(CODE_ERROR//': if stress is present then potential must &
       &also be present.')
    call err()
  endif
  
  allocate( harmonic_energies(size(mode_displacements)),   &
          & harmonic_forces(size(mode_displacements)),     &
          & stat=ialloc); call err(ialloc)
  if (present(potential)) then
   allocate( anharmonic_energies(size(mode_displacements)), &
           & anharmonic_forces(size(mode_displacements)),   &
           & stat=ialloc); call err(ialloc)
  endif
  if (present(stress)) then
    allocate( anharmonic_pressures(size(mode_displacements)), &
            & stat=ialloc); call err(ialloc)
  endif
  do i=1,size(mode_displacements)
    displacement = RealModeDisplacement([mode],[mode_displacements(i)])
    
    harmonic_energies(i) = 0.5_dp                &
                       & * mode%spring_constant  &
                       & * mode_displacements(i) &
                       & * mode_displacements(i)
    harmonic_forces(i) = - mode%spring_constant &
                     & * mode_displacements(i)
    
    if (present(potential)) then
      anharmonic_energies(i) = ( potential%energy(displacement)   &
                           &   - potential%undisplaced_energy() ) &
                           & / anharmonic_data%anharmonic_supercell%sc_size
      anharmonic_force = potential%force(displacement)
      anharmonic_forces(i) = anharmonic_force%force(mode) &
                         & / anharmonic_data%anharmonic_supercell%sc_size
    endif
    
    if (present(stress)) then
      anharmonic_pressures(i) = trace( stress%stress(displacement)   &
                            &        - stress%undisplaced_stress() ) &
                            & / 3
    endif
  enddo
  
  this = ModeMap( mode%id,                    &
                & mode%frequency,             &
                & mode_displacements,         &
                & l2_cartesian_displacements, &
                & harmonic_energies,          &
                & harmonic_forces,            &
                & anharmonic_energies,        &
                & anharmonic_forces,          &
                & anharmonic_pressures        )
end procedure

module procedure read_ModeMap
  integer               :: mode_id
  real(dp)              :: harmonic_frequency
  real(dp), allocatable :: mode_displacements(:)
  real(dp), allocatable :: l2_cartesian_displacements(:)
  real(dp), allocatable :: harmonic_energies(:)
  real(dp), allocatable :: harmonic_forces(:)
  real(dp), allocatable :: anharmonic_energies(:)
  real(dp), allocatable :: anharmonic_forces(:)
  real(dp), allocatable :: anharmonic_pressures(:)
  real(dp), allocatable :: sampled_energies(:)
  real(dp), allocatable :: sampled_forces(:)
  real(dp), allocatable :: sampled_pressures(:)
  
  type(String), allocatable :: line(:)
  
  integer :: no
  
  integer :: i,j,ialloc
  
  select type(this); type is(ModeMap)
    line = split_line(input(1))
    mode_id = int(line(3))
    
    line = split_line(input(2))
    harmonic_frequency = dble(line(3))
    
    no = size(input)-3
    
    line = split_line(input(3))
    allocate( mode_displacements(no),         &
            & l2_cartesian_displacements(no), &
            & harmonic_energies(no),          &
            & harmonic_forces(no),            &
            & stat=ialloc); call err(ialloc)
    if (any(line==str('Anharmonic'))) then
      allocate( anharmonic_energies(no), &
              & anharmonic_forces(no),   &
              & stat=ialloc); call err(ialloc)
      if (any(line==str('pressure'))) then
        allocate(anharmonic_pressures(no), stat=ialloc); call err(ialloc)
      endif
    endif
    if (any(line==str('Sampled'))) then
      allocate( sampled_energies(no), &
              & sampled_forces(no),   &
              & stat=ialloc); call err(ialloc)
      if (any(line==str('pressure'))) then
        allocate(sampled_pressures(no), stat=ialloc); call err(ialloc)
      endif
    endif
    
    do i=1,size(input)-4
      line = split_line(input(i+4))
      mode_displacements(i) = dble(line(1))
      l2_cartesian_displacements(i) = dble(line(2))
      harmonic_energies(i) = dble(line(3))
      harmonic_forces(i) = dble(line(4))
      j = 4
      if (allocated(anharmonic_energies)) then
        anharmonic_energies(i) = dble(line(j+1))
        anharmonic_forces(i) = dble(line(j+2))
        j = j+2
      endif
      if (allocated(anharmonic_pressures)) then
        anharmonic_pressures(i) = dble(line(j+1))
        j = j+1
      endif
      if (allocated(sampled_energies)) then
        sampled_energies(i) = dble(line(j+1))
        sampled_forces(i) = dble(line(j+2))
        j = j+2
      endif
      if (allocated(sampled_pressures)) then
        sampled_pressures(i) = dble(line(j+1))
      endif
    enddo
    
    this = ModeMap( mode_id,                    &
                  & harmonic_frequency,         &
                  & mode_displacements,         &
                  & l2_cartesian_displacements, &
                  & harmonic_energies,          &
                  & harmonic_forces,            &
                  & anharmonic_energies,        &
                  & anharmonic_forces,          &
                  & anharmonic_pressures,       &
                  & sampled_energies,           &
                  & sampled_forces,             &
                  & sampled_pressures           )
  class default
    call err()
  end select
end procedure

module procedure write_ModeMap
  type(String) :: line
  
  integer :: i
  
  select type(this); type is(ModeMap)
    output = [ 'Mode ID: '//this%mode_id,                      &
             & 'Harmonic frequency: '//this%harmonic_frequency ]
    
    line = 'Mode displacement | L2 Cartesian displacement | Harmonic energy | &
           &Harmonic force'
    if (allocated(this%anharmonic_energies)) then
      line = line//' | Anharmonic energy | Anharmonic force'
    endif
    if (allocated(this%anharmonic_pressures)) then
      line = line//' | Anharmonic pressure'
    endif
    if (allocated(this%sampled_energies)) then
      line = line//' | Sampled energy | Sampled force'
    endif
    if (allocated(this%sampled_pressures)) then
      line = line//' | Sampled pressure'
    endif
    output = [output, line]
    
    do i=1,size(this%mode_displacements)
      line = this%mode_displacements(i)         //' '// &
           & this%l2_cartesian_displacements(i) //' '// &
           & this%harmonic_energies(i)          //' '// &
           & this%harmonic_forces(i)
      if (allocated(this%anharmonic_energies)) then
        line = line                        //' '// &
             & this%anharmonic_energies(i) //' '// &
             & this%anharmonic_forces(i)
      endif
      if (allocated(this%anharmonic_pressures)) then
        line = line //' '// &
             & this%anharmonic_pressures(i)
      endif
      if (allocated(this%sampled_energies)) then
        line = line                     //' '// &
             & this%sampled_energies(i) //' '// &
             & this%sampled_forces(i)
      endif
      if (allocated(this%sampled_pressures)) then
        line = line //' '// &
             & this%sampled_pressures(i)
      endif
      
      output = [output, line]
    enddo
  class default
    call err()
  end select
end procedure

module procedure new_ModeMap_Strings
  call this%read(input)
end procedure

module procedure new_ModeMap_StringArray
  this = ModeMap(str(input))
end procedure
end submodule
