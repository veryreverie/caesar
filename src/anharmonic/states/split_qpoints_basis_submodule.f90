submodule (caesar_split_qpoints_basis_module) caesar_split_qpoints_basis_submodule
  use caesar_states_module
contains

module procedure new_QpointModeIDs
  this%mode_ids = mode_ids(sort(mode_ids))
  this%paired_mode_ids = paired_mode_ids(sort(paired_mode_ids))
end procedure

module procedure new_SplitQpointsBasis
  integer :: i
  
  this%supercell_size        = supercell_size
  this%maximum_power         = maximum_power
  this%expansion_order       = expansion_order
  this%subspace_id           = subspace_id
  this%frequency             = frequency
  this%wavevectors           = wavevectors
  this%qpoints               = qpoints
  this%qpoints_to_integrate_ = [(i, i=1, size(qpoints))]
end procedure

module procedure representation_SplitQpointsBasis
  output = 'split q-points basis'
end procedure

module procedure new_SplitQpointsBasis_subspace
  type(ComplexMode), allocatable :: subspace_modes(:)
  type(QpointData),  allocatable :: subspace_qpoints(:)
  
  type(ComplexMode),   allocatable :: qpoint_modes(:)
  type(QpointModeIDs), allocatable :: qpoint_mode_ids(:)
  
  type(WavevectorBasis), allocatable :: wavevectors(:)
  
  integer :: i,ialloc
  
  ! List the modes and corresponding q-points in the subspace.
  ! The q-points are de-duplicated, and then q-points with id>paired_id are
  !    removed.
  subspace_modes = subspace%modes(modes)
  subspace_qpoints = subspace%qpoints(modes, qpoints)
  
  subspace_qpoints = subspace_qpoints(set(subspace_qpoints%id))
  subspace_qpoints = subspace_qpoints(                                &
     & filter(subspace_qpoints%id<=subspace_qpoints%paired_qpoint_id) )
  
  allocate( qpoint_mode_ids(size(subspace_qpoints)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(qpoint_mode_ids)
    qpoint_modes = subspace_modes(                                &
       & filter(subspace_modes%qpoint_id==subspace_qpoints(i)%id) )
    qpoint_mode_ids(i) = QpointModeIDs( qpoint_modes%id,       &
                                      & qpoint_modes%paired_id )
  enddo
  
  ! Generate a basis at the first q-point only.
  wavevectors = WavevectorBasis( subspace,                  &
                               & frequency,                 &
                               & modes,                     &
                               & qpoints,                   &
                               & maximum_power,             &
                               & potential_expansion_order, &
                               & supercell%sc_size,         &
                               & symmetries,                &
                               & subspace_qpoints(1)        )
  
  output = SplitQpointsBasis( supercell%sc_size,         &
                            & maximum_power,             &
                            & potential_expansion_order, &
                            & subspace%id,               &
                            & frequency,                 &
                            & wavevectors,               &
                            & qpoint_mode_ids            )
end procedure

module procedure ground_state_wavefunction
  real(dp) :: mass
  
  real(dp)                  :: coefficient
  type(String), allocatable :: terms(:)
  
  integer :: i,ialloc
  
  ! Calculate the (geometric) average mass.
  mass = product(supercell%atoms%mass())
  mass = mass**(1.0_dp/size(supercell%atoms))
  
  ! Calculate the coefficient.
  coefficient = 1
  allocate(terms(0), stat=ialloc); call err(ialloc)
  do i=1,size(subspace)
    if (subspace%mode_ids(i)==subspace%paired_ids(i)) then
      ! |0_i> = sqrt(sqrt(m*w/pi)) exp(- 1/2 N w (u_i)^2 )
      coefficient = coefficient * (mass*this%frequency/PI)**0.25_dp
      terms = [terms, 'u'//subspace%mode_ids(i)//'^2']
    elseif (subspace%mode_ids(i)<subspace%paired_ids(i)) then
      ! |0_i,0_j> = sqrt(2*m*w/pi) exp(- N w |u_i|^2 )
      coefficient = coefficient * (2*mass*this%frequency/PI)**0.5_dp
      terms = [ terms,                                                    &
              & '2*u'//subspace%mode_ids(i)//'*u'//subspace%paired_ids(i) ]
    endif
  enddo
  
  output = coefficient//'*e^('               // &
     & (-supercell%sc_size*this%frequency/2) // &
     & '*('//join(terms,'+')//'))'
end procedure

module procedure initial_states_SplitQpointsBasis
  associate( gamma_basis =>                                               &
           & this%wavevectors(first(is_int(this%wavevectors%wavevector))) )
    output = gamma_basis%initial_states( subspace,       &
                                       & thermal_energy, &
                                       & anharmonic_data )
  end associate
end procedure

module procedure calculate_states_SplitQpointsBasis
  output = calculate_states( this%wavevectors,        &
                           & subspace,                &
                           & subspace_potential,      &
                           & thermal_energy,          &
                           & state_energy_cutoff,     &
                           & convergence_data,        &
                           & anharmonic_data          )
end procedure

module procedure process_subspace_potential_SplitQpointsBasis
  type(SplitQpointsBasis) :: basis
  
  type(PotentialBasePointer) :: integrated_potential
  real(dp)                   :: correction_energy
  
  integer :: i
  
  ! Calculate (prod_{j/=i}<j|)V(prod_{j/=i}|j>).
  if (size(this%qpoints)==1) then
    return
  endif
  
  basis = this
  
  basis%qpoints_to_integrate_ = [(i,i=2,size(this%qpoints))]
  
  call potential%braket( states,                           &
                       & subspace        = subspace,       &
                       & subspace_basis  = basis,          &
                       & whole_subspace  = .false.,        &
                       & anharmonic_data = anharmonic_data )
  
  ! Calculate (prod_i<i|)V(prod_i|i>), and use this to correct for
  !    over-counting. See main vscf method for details.
  integrated_potential = PotentialBasePointer(potential)
  
  basis%qpoints_to_integrate_ = [1]
  
  call integrated_potential%braket( states,                           &
                                  & subspace        = subspace,       &
                                  & subspace_basis  = basis,          &
                                  & anharmonic_data = anharmonic_data )
  
  correction_energy = integrated_potential%undisplaced_energy() &
                  & * (1.0_dp-size(this%qpoints))/size(this%qpoints)
  
  call potential%add_constant(correction_energy)
end procedure

module procedure process_subspace_stress_SplitQpointsBasis
  type(SplitQpointsBasis) :: basis
  
  type(StressBasePointer) :: integrated_stress
  type(RealMatrix)        :: correction_stress
  
  integer :: i
  
  ! Calculate (prod_{j/=i}<j|)V(prod_{j/=i}|j>).
  if (size(this%qpoints)==1) then
    return
  endif
  
  basis = this
  
  basis%qpoints_to_integrate_ = [(i,i=2,size(this%qpoints))]
  
  call stress%braket( states,                           &
                    & subspace        = subspace,       &
                    & subspace_basis  = basis,          &
                    & whole_subspace  = .false.,        &
                    & anharmonic_data = anharmonic_data )
  
  ! Calculate (prod_i<i|)V(prod_i|i>), and use this to correct for
  !    over-counting. See main vscf method for details.
  integrated_stress = StressBasePointer(stress)
  
  basis%qpoints_to_integrate_ = [1]
  
  call integrated_stress%braket( states,                           &
                               & subspace        = subspace,       &
                               & subspace_basis  = basis,          &
                               & anharmonic_data = anharmonic_data )
  
  correction_stress = integrated_stress%undisplaced_stress() &
                  & * (1.0_dp-size(this%qpoints))/size(this%qpoints)
  
  call stress%add_constant(correction_stress)
end procedure

module procedure wavefunction_SplitQpointsBasis
  ! TODO
  output = ''
end procedure

module procedure mode_ids_SplitQpointsBasis
  integer :: i
  
  output = [( this%qpoints(this%qpoints_to_integrate_(i))%mode_ids, &
            & i=1,                                                  &
            & size(this%qpoints_to_integrate_)                      )]
end procedure

module procedure paired_mode_ids_SplitQpointsBasis
  integer :: i
  
  output = [( this%qpoints(this%qpoints_to_integrate_(i))%paired_mode_ids, &
            & i=1,                                                         &
            & size(this%qpoints_to_integrate_)                             )]
end procedure

module procedure inner_product_SplitQpointsBasis
  type(WavevectorState)              :: split_bra
  type(WavevectorState), allocatable :: split_ket
  
  integer :: i
  
  split_bra = WavevectorState(bra)
  if (present(ket)) then
    split_ket = WavevectorState(ket)
  endif
  
  i = first(this%wavevectors%wavevector == split_bra%wavevector)
  
  output = this%wavevectors(i)%inner_product( split_bra,       &
                                            & split_ket,       &
                                            & subspace,        &
                                            & anharmonic_data  )
  
  output = output * size(this%qpoints_to_integrate_)
end procedure

module procedure integrate_BasisState_SplitQpointsBasis
  type(WavevectorState)              :: split_bra
  type(WavevectorState), allocatable :: split_ket
  
  type(SparseMonomial), allocatable :: monomials(:)
  
  integer :: i
  
  ! Convert the bra (and the ket if present) to type WavevectorState.
  split_bra = WavevectorState(bra)
  if (present(ket)) then
    split_ket = WavevectorState(ket)
  endif
  
  ! Find the basis at the correct wavevector for the states.
  i = first(this%wavevectors%wavevector == split_bra%wavevector)
  
  ! Split the monomial by q-point.
  monomials = split_monomial(this,monomial)
  
  ! Integrate across all q-points to be integrated.
  output = product(this%wavevectors(i)%integrate( split_bra,       &
                                                & monomials,       &
                                                & split_ket,       &
                                                & subspace,        &
                                                & anharmonic_data  ))
end procedure

module procedure split_monomial
  type(QpointModeIDs), allocatable :: qpoints(:)
  
  integer :: modes_per_qpoint
  
  integer :: i,ialloc
  
  qpoints = this%qpoints(this%qpoints_to_integrate_)
  
  modes_per_qpoint = size(monomial%modes)/size(qpoints)
  
  allocate(output(size(qpoints)), stat=ialloc); call err(ialloc)
  do i=1,size(qpoints)
    ! Separate off the modes at each q-point.
    output(i) = SparseMonomial(                           &
       & modes=monomial%modes( (i-1)*modes_per_qpoint+1   &
       &                     :  i   *modes_per_qpoint   ) )
    ! Transform each mode to be at the first q-point.
    output(i)%modes%id = this%qpoints(1)%mode_ids
    output(i)%modes%paired_id = this%qpoints(1)%paired_mode_ids
  enddo
end procedure

module procedure kinetic_energy_SplitQpointsBasis
  type(WavevectorState)              :: split_bra
  type(WavevectorState), allocatable :: split_ket
  
  integer :: i
  
  split_bra = WavevectorState(bra)
  if (present(ket)) then
    split_ket = WavevectorState(ket)
  endif
  
  i = first(this%wavevectors%wavevector == split_bra%wavevector)
  
  output = this%wavevectors(i)%kinetic_energy( split_bra,      &
                                             & split_ket,      &
                                             & subspace,       &
                                             & anharmonic_data )
  
  output = output * size(this%qpoints_to_integrate_)
end procedure

module procedure harmonic_potential_energy_SplitQpointsBasis
  type(WavevectorState)              :: split_bra
  type(WavevectorState), allocatable :: split_ket
  
  integer :: i
  
  split_bra = WavevectorState(bra)
  if (present(ket)) then
    split_ket = WavevectorState(ket)
  endif
  
  i = first(this%wavevectors%wavevector == split_bra%wavevector)
  
  output = this%wavevectors(i)%harmonic_potential_energy( split_bra,      &
                                                        & split_ket,      &
                                                        & subspace,       &
                                                        & anharmonic_data )
  
  output = output * size(this%qpoints_to_integrate_)
end procedure

module procedure kinetic_stress_SplitQpointsBasis
  type(WavevectorState)              :: split_bra
  type(WavevectorState), allocatable :: split_ket
  
  integer :: i
  
  split_bra = WavevectorState(bra)
  if (present(ket)) then
    split_ket = WavevectorState(ket)
  endif
  
  i = first(this%wavevectors%wavevector == split_bra%wavevector)
  
  output = this%wavevectors(i)%kinetic_stress( split_bra,         &
                                             & split_ket,         &
                                             & subspace,          &
                                             & stress_prefactors, &
                                             & anharmonic_data    )
  
  output = output * size(this%qpoints_to_integrate_)
end procedure

module procedure thermodynamic_data_SplitQpointsBasis
  type(RealMatrix) :: stress
  
  integer :: i
  
  ! Calculate the thermodynamic properties for one q-point in the system.
  output = core_shell_thermodynamics(     &
     & thermal_energy,                    &
     & this%frequency,                    &
     & this%supercell_size,               &
     & size(subspace)/size(this%qpoints), &
     & subspace,                          &
     & this%wavevectors,                  &
     & states,                            &
     & subspace_potential,                &
     & subspace_stress,                   &
     & stress_prefactors,                 &
     & anharmonic_data                    )
  
  ! Multiply by the number of q-points, and add in the constant term.
  output = output * size(this%qpoints)
  
  ! Average the stress across all symmetries.
  ! Each subspace obeys all symmetries, so the stress must too.
  ! This is needed because the properties are only calculated at one q-point,
  !    which will not in general obey the symmetries.
  if (allocated(output%stress)) then
    associate(symmetries=>anharmonic_data%structure%symmetries)
      stress = dblemat(zeroes(3,3))
      do i=1,size(symmetries)
        stress = stress                         &
             & + symmetries(i)%cartesian_tensor &
             & * output%stress                  &
             & * transpose(symmetries(i)%cartesian_tensor)
      enddo
      output%stress = stress / size(symmetries)
    end associate
  endif
end procedure

module procedure wavefunctions_SplitQpointsBasis
  type(WavevectorStates), pointer :: split_states
  
  type(String)                    :: ground_state
  type(String), allocatable       :: state_wavefunctions(:)
  type(SplitQpointsWavefunctions) :: wavefunctions
  
  split_states => wavevector_states_pointer(states)
  
  ! Construct the wavefunction of |0>.
  ground_state = this%ground_state_wavefunction( &
          & subspace,                            &
          & anharmonic_data%anharmonic_supercell )
  
  ! Construct the wavefunctions of each state, in terms of |0>.
  state_wavefunctions = this%wavefunction(  &
     & split_states%states,                 &
     & anharmonic_data%anharmonic_supercell )
  
  ! Pack the wavefunctions into the output,
  !    of concrete type SplitQpointsWavefunctions.
  wavefunctions = SplitQpointsWavefunctions( subspace%id,           &
                                           & subspace%mode_ids,     &
                                           & subspace%paired_ids,   &
                                           & ground_state,          &
                                           & split_states%energies, &
                                           & state_wavefunctions    )
  
  ! Convert the output to abstract class SubspaceWavefunctions.
  output = SubspaceWavefunctionsPointer(wavefunctions)
end procedure

module procedure integrate_BasisStates_SplitQpointsBasis
  type(SparseMonomial), allocatable :: monomials(:)
  
  complex(dp) :: term
  
  integer :: i,j
  
  ! Split the monomial by q-point.
  monomials = split_monomial(this,monomial)
  
  ! Integrate across all q-points to be integrated.
  output = 1.0_dp
  do i=1,size(this%qpoints_to_integrate_)
    term = 0
    do j=1,size(this%wavevectors)
      term = term + this%wavevectors(j)%integrate( states,         &
                                                 & monomials(i),   &
                                                 & subspace,       &
                                                 & anharmonic_data )
    enddo
    output = output * term
  enddo
end procedure

module procedure free_energy_gradient_SplitQpointsBasis
  output = free_energy_gradient( this%wavevectors,      &
       &                         subspace_potential,    &
       &                         basis_functions,       &
       &                         subspace,              &
       &                         states,                &
       &                         thermal_energy,        &
       &                         state_energy_cutoff,   &
       &                         anharmonic_data      ) &
       & * size(this%qpoints)
end procedure

module procedure select_symmetries_SplitQpointsBasis
  logical, allocatable :: selected(:)
  
  type(ComplexMode) :: mode
  type(ComplexMode) :: paired_mode
  
  type(QpointData) :: qpoint
  type(QpointData) :: paired_qpoint
  
  type(SymmetryOperator) :: symmetry
  
  type(QpointData) :: transformed_qpoint
  
  integer :: i,ialloc
  
  allocate(selected(size(symmetries)), stat=ialloc); call err(ialloc)
  
  mode = anharmonic_data%complex_modes(          &
     & first( anharmonic_data%complex_modes%id   &
     &     == this%qpoints(1)%mode_ids(1)      ) )
  paired_mode = anharmonic_data%complex_modes(     &
     & first( anharmonic_data%complex_modes%id     &
     &     == this%qpoints(1)%paired_mode_ids(1) ) )
  
  qpoint = select_qpoint(mode, anharmonic_data%qpoints)
  paired_qpoint = select_qpoint(paired_mode, anharmonic_data%qpoints)
  
  do i=1,size(symmetries)
    symmetry = anharmonic_data%structure%symmetries(first(                  &
       & anharmonic_data%structure%symmetries%id==symmetries(i)%symmetry_id ))
    transformed_qpoint = symmetry * qpoint
    selected(i) = ( transformed_qpoint==qpoint        &
             & .or. transformed_qpoint==paired_qpoint )
  enddo
  
  output = symmetries(filter(selected))
end procedure

module procedure read_QpointModeIDs
  integer, allocatable :: mode_ids(:)
  integer, allocatable :: paired_mode_ids(:)
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(QpointModeIDs)
    line = tokens(input, delimiter=',')
    
    mode_ids = int(tokens(line(1), first=2))
    paired_mode_ids = int(tokens(line(2), first=3))
    
    this = QpointModeIDs(mode_ids, paired_mode_ids)
  class default
    call err()
  end select
end procedure

module procedure write_QpointModeIDs
  select type(this); type is(QpointModeIDs)
    output = 'Modes: '//this%mode_ids// &
           & ', Paired modes: '//this%paired_mode_ids
  class default
    call err()
  end select
end procedure

module procedure new_QpointModeIDs_String
  call this%read(input)
end procedure

module procedure read_SplitQpointsBasis
  !integer                            :: supercell_size
  !integer                            :: maximum_power
  !integer                            :: expansion_order
  !integer                            :: subspace_id
  !real(dp)                           :: frequency
  !type(WavevectorBasis), allocatable :: wavevectors(:)
  !
  !integer, allocatable :: starting_lines(:)
  !integer, allocatable :: ending_lines(:)
  !
  !type(String), allocatable :: line(:)
  !type(String), allocatable :: wavevector_lines(:)
  !
  !integer :: i,ialloc
  
  ! TODO
  call err()
  
  !select type(this); type is(SplitQpointsBasis)
  !  line = split_line(input(1))
  !  maximum_power = int(line(4))
  !  
  !  line = split_line(input(2))
  !  expansion_order = int(line(4))
  !  
  !  line = split_line(input(3))
  !  subspace_id = int(line(3))
  !  
  !  line = split_line(input(4))
  !  frequency = dble(line(3))
  !  
  !  starting_lines = [integer::]
  !  do i=5,size(input)
  !    line = split_line(input(i))
  !    if (size(line)>0) then
  !      if (line(1)=='Wavevector') then
  !        starting_lines = [starting_lines, i]
  !      endif
  !    endif
  !  enddo
  !  
  !  ending_lines = [starting_lines(2:)-1, size(input)]
  !  
  !  allocate(wavevectors(size(starting_lines)), stat=ialloc); call err(ialloc)
  !  do i=1,size(starting_lines)
  !    wavevector_lines = [input(:4), input(starting_lines(i):ending_lines(i))]
  !    wavevectors(i) = WavevectorBasis(wavevector_lines)
  !  enddo
  !  
  !  this = SplitQpointsBasis( maximum_power,   &
  !                          & expansion_order, &
  !                          & subspace_id,     &
  !                          & frequency,       &
  !                          & wavevectors      )
  !class default
  !  call err()
  !end select
end procedure

module procedure write_SplitQpointsBasis
  type(String), allocatable :: wavevector_strings(:)
  
  integer :: i
  
  select type(this); type is(SplitQpointsBasis)
    output = [ 'Supercell size   : '//this%supercell_size,  &
             & 'Maximum power    : '//this%maximum_power,   &
             & 'Expansion order  : '//this%expansion_order, &
             & 'Subspace         : '//this%subspace_id,     &
             & 'Frequency        : '//this%frequency,       &
             & str('Wavevectors:'),                         &
             & str(repeat('-',50))                          ]
    do i=1,size(this%wavevectors)
      wavevector_strings = str(this%wavevectors(i))
      output = [ output, wavevector_strings(5:), str(repeat('-',50)) ]
    enddo
    output = [output, str('q-point mode ids:')]
    do i=1,size(this%qpoints)
      output = [output, str(this%qpoints(i))]
    enddo
  class default
    call err()
  end select
end procedure

module procedure new_SplitQpointsBasis_Strings
  call this%read(input)
end procedure

module procedure new_SplitQpointsBasis_StringArray
  this = SplitQpointsBasis(str(input))
end procedure
end submodule
