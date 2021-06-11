submodule (caesar_anharmonic_data_module) caesar_anharmonic_data_submodule
  use caesar_anharmonic_common_module
contains

module procedure new_AnharmonicData
  integer :: i
  
  this%structure                 = structure
  this%anharmonic_supercell      = anharmonic_supercell
  this%qpoints                   = qpoints
  this%complex_modes             = complex_modes
  this%real_modes                = real_modes
  this%degenerate_subspaces      = degenerate_subspaces
  this%degenerate_symmetries     = degenerate_symmetries
  this%subspace_couplings        = subspace_couplings
  this%max_subspace_coupling     = max_subspace_coupling
  this%max_qpoint_coupling       = max_qpoint_coupling
  this%potential_expansion_order = potential_expansion_order
  this%vscf_basis_functions_only = vscf_basis_functions_only
  this%max_displacement          = max_displacement
  
  this%qpoint_symmetry_groups = [(                                  &
     & structure%symmetries(i)%qpoint_symmetry_group(this%qpoints), &
     & i=1,                                                         &
     & size(structure%symmetries)                                   )]
end procedure

module procedure new_AnharmonicData_data
  ! Degeneracy data.
  type(DegenerateSubspace), allocatable :: degenerate_subspaces(:)
  type(DegenerateSymmetry), allocatable :: degenerate_symmetries(:)
  
  ! Coupling data.
  type(SubspaceCoupling), allocatable :: subspace_couplings(:)
  
  ! Temporary variables.
  integer :: i,ialloc
  
  ! Retrieve data on how normal modes are grouped into subspaces
  !    of degenerate modes.
  degenerate_subspaces = process_degeneracies( &
        & interpolated_supercell%complex_modes )
  
  ! Generate the symmetry operators in each degenerate subspace.
  allocate( degenerate_symmetries(size(structure%symmetries)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(structure%symmetries)
    degenerate_symmetries(i) = DegenerateSymmetry( &
           & structure%symmetries(i),              &
           & degenerate_subspaces,                 &
           & interpolated_supercell%complex_modes, &
           & interpolated_supercell%qpoints        )
  enddo
  
  ! Generate all sets of coupled subspaces, up to max_subspace_coupling.
  call print_line('Generating couplings between subspaces.')
  subspace_couplings = generate_coupled_subspaces( degenerate_subspaces, &
                                                 & max_subspace_coupling )
  
  ! Load anharmonic data into container.
  this = AnharmonicData( structure,                            &
                       & interpolated_supercell%supercell,     &
                       & interpolated_supercell%qpoints,       &
                       & interpolated_supercell%complex_modes, &
                       & interpolated_supercell%real_modes,    &
                       & degenerate_subspaces,                 &
                       & degenerate_symmetries,                &
                       & subspace_couplings,                   &
                       & max_subspace_coupling,                &
                       & max_qpoint_coupling,                  &
                       & potential_expansion_order,            &
                       & vscf_basis_functions_only,            &
                       & max_displacement                      )
end procedure

module procedure read_AnharmonicData
  type(StructureData)                   :: structure
  type(StructureData)                   :: anharmonic_supercell
  type(QpointData),         allocatable :: qpoints(:)
  type(ComplexMode),        allocatable :: complex_modes(:)
  type(RealMode),           allocatable :: real_modes(:)
  type(DegenerateSubspace), allocatable :: degenerate_subspaces(:)
  type(DegenerateSymmetry), allocatable :: degenerate_symmetries(:)
  type(SubspaceCoupling),   allocatable :: subspace_couplings(:)
  integer                               :: max_subspace_coupling
  integer                               :: max_qpoint_coupling
  integer                               :: potential_expansion_order
  logical                               :: vscf_basis_functions_only
  real(dp)                              :: maximum_weighted_displacement
  real(dp)                              :: frequency_of_max_displacement
  type(MaxDisplacement)                 :: max_displacement
  
  type(String),      allocatable :: headers(:)
  type(StringArray), allocatable :: sections(:)
  
  type(String) :: separator
  
  integer :: i,ialloc
  
  separator = repeat('=',70)
  
  select type(this); type is(AnharmonicData)
    ! Split the file into sections.
    sections = split_into_sections(input, separating_line=separator)
    
    ! Trim the header line from each section.
    allocate(headers(size(sections)), stat=ialloc); call err(ialloc)
    do i=1,size(sections)
      headers(i) = trim(lower_case(sections(i)%strings(1)))
      sections(i)%strings = sections(i)%strings(2:)
    enddo
    
    ! Parse each section into the relevant variable.
    i = first(headers=='structure', default=0)
    if (i==0) then
      call print_line(ERROR//': anharmonic data must contain "structure".')
      call err()
    else
      structure = StructureData(sections(i))
    endif
    
    i = first(headers=='anharmonic supercell', default=0)
    if (i==0) then
      call print_line(ERROR//': anharmonic data must contain "anharmonic &
         &supercell".')
      call err()
    else
      anharmonic_supercell = StructureData(sections(i))
    endif
    
    i = first(headers=='q-points', default=0)
    if (i==0) then
      call print_line(ERROR//': anharmonic data must contain "q-points".')
      call err()
    else
      qpoints = QpointData(split_into_sections(sections(i)))
    endif
    
    i = first(headers=='complex modes', default=0)
    if (i==0) then
      call print_line(ERROR//': anharmonic data must contain "complex modes".')
      call err()
    else
      complex_modes = ComplexMode(split_into_sections(sections(i)))
    endif
    
    i = first(headers=='real modes', default=0)
    if (i==0) then
      call print_line(ERROR//': anharmonic data must contain "real modes".')
      call err()
    else
      real_modes = RealMode(split_into_sections(sections(i)))
    endif
    
    i = first(headers=='degenerate subspaces', default=0)
    if (i==0) then
      call print_line(ERROR//': anharmonic data must contain "degenerate &
         &subspaces".')
      call err()
    else
      degenerate_subspaces = DegenerateSubspace( &
              & split_into_sections(sections(i)) )
    endif
    
    i = first(headers=='degenerate symmetries', default=0)
    if (i==0) then
      call print_line(ERROR//': anharmonic data must contain "degenerate &
         &symmetries".')
      call err()
    else
      degenerate_symmetries = DegenerateSymmetry( &
               & split_into_sections(sections(i)) )
    endif
    
    i = first(headers=='subspace couplings', default=0)
    if (i==0) then
      call print_line(ERROR//': anharmonic data must contain "subspace &
         &couplings".')
      call err()
    else
      subspace_couplings = SubspaceCoupling(sections(i)%strings)
    endif
    
    ! maximum_coupling_order included for legacy reasons.
    i = first( headers == 'max subspace coupling'   &
        & .or. headers == 'maximum_coupling_order', &
        &      default=0                            )
    if (i==0) then
      call print_line(WARNING//': anharmonic data does not contain "max &
         &subspace coupling". Defaulting to 1.')
      max_subspace_coupling = 1
    else
      max_subspace_coupling = int(sections(i)%strings(1))
    endif
    
    i = first(headers=='max q-point coupling', default=0)
    if (i==0) then
      call print_line(WARNING//': anharmonic data does not contain "max &
         &q-point coupling". Defaulting to 1.')
      max_qpoint_coupling = 1
    else
      max_qpoint_coupling = int(sections(i)%strings(1))
    endif
    
    i = first(headers=='potential expansion order', default=0)
    if (i==0) then
      call print_line(ERROR//': anharmonic data must contain "potential &
         &expansion order".')
      call err()
    else
      potential_expansion_order = int(sections(i)%strings(1))
    endif
    
    i = first(headers=='vscf basis functions only', default=0)
    if (i==0) then
      call print_line(ERROR//': anharmonic data must contain "vscf basis &
         &functions only".')
      call err()
    else
      vscf_basis_functions_only = lgcl(sections(i)%strings(1))
    endif
    
    i = first(headers=='maximum weighted displacement', default=0)
    if (i==0) then
      call print_line(ERROR//': anharmonic data must contain "maximum &
         &weighted displacement".')
      call err()
    else
      maximum_weighted_displacement = dble(sections(i)%strings(1))
    endif
    
    i = first(headers=='frequency of maximum displacement', default=0)
    if (i==0) then
      call print_line(ERROR//': anharmonic data must contain "frequency of &
         &maximum displacement".')
      call err()
    else
      frequency_of_max_displacement = dble(sections(i)%strings(1))
    endif
    
    max_displacement = MaxDisplacement(                                 &
       & maximum_weighted_displacement = maximum_weighted_displacement, &
       & frequency_of_max_displacement = frequency_of_max_displacement  )
    
    ! Construct the output.
    this = AnharmonicData( structure,                 &
                         & anharmonic_supercell,      &
                         & qpoints,                   &
                         & complex_modes,             &
                         & real_modes,                &
                         & degenerate_subspaces,      &
                         & degenerate_symmetries,     &
                         & subspace_couplings,        &
                         & max_subspace_coupling,     &
                         & max_qpoint_coupling,       &
                         & potential_expansion_order, &
                         & vscf_basis_functions_only, &
                         & max_displacement           )
  class default
    call err()
  end select
end procedure

module procedure write_AnharmonicData
  type(String) :: separator
  
  separator = repeat('=',70)
  
  select type(this); type is(AnharmonicData)
    output = [ str('Structure'),                                         &
             & str(this%structure),                                      &
             & separator,                                                &
             & str('Anharmonic Supercell'),                              &
             & str(this%anharmonic_supercell),                           &
             & separator,                                                &
             & str('q-points'),                                          &
             & str(this%qpoints, separating_line=''),                    &
             & separator,                                                &
             & str('Complex modes'),                                     &
             & str(this%complex_modes, separating_line=''),              &
             & separator,                                                &
             & str('Real modes'),                                        &
             & str(this%real_modes, separating_line=''),                 &
             & separator,                                                &
             & str('Degenerate subspaces'),                              &
             & str(this%degenerate_subspaces, separating_line=''),       &
             & separator,                                                &
             & str('Degenerate symmetries'),                             &
             & str(this%degenerate_symmetries, separating_line=''),      &
             & separator,                                                &
             & str('Subspace couplings'),                                &
             & str(this%subspace_couplings),                             &
             & separator,                                                &
             & str('max subspace coupling'),                             &
             & str(this%max_subspace_coupling),                          &
             & separator,                                                &
             & str('max q-point coupling'),                              &
             & str(this%max_subspace_coupling),                          &
             & separator,                                                &
             & str('Potential expansion order'),                         &
             & str(this%potential_expansion_order),                      &
             & separator,                                                &
             & str('VSCF basis functions only'),                         &
             & str(this%vscf_basis_functions_only),                      &
             & separator,                                                &
             & str('Maximum weighted displacement'),                     &
             & str(this%max_displacement%maximum_weighted_displacement), &
             & separator,                                                &
             & str('Frequency of maximum displacement'),                 &
             & str(this%max_displacement%frequency_of_max_displacement)  ]
       
  class default
    call err()
  end select
end procedure

module procedure new_AnharmonicData_Strings
  call this%read(input)
end procedure

module procedure new_AnharmonicData_StringArray
  this = AnharmonicData(str(input))
end procedure
end submodule
