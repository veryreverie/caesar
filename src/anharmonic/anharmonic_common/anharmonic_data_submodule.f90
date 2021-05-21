submodule (caesar_anharmonic_data_module) caesar_anharmonic_data_submodule
  use caesar_anharmonic_common_module
contains

module procedure new_AnharmonicData
  this%structure                     = structure
  this%anharmonic_supercell          = anharmonic_supercell
  this%qpoints                       = qpoints
  this%complex_modes                 = complex_modes
  this%real_modes                    = real_modes
  this%degenerate_subspaces          = degenerate_subspaces
  this%degenerate_symmetries         = degenerate_symmetries
  this%subspace_couplings            = subspace_couplings
  this%maximum_coupling_order        = maximum_coupling_order
  this%potential_expansion_order     = potential_expansion_order
  this%vscf_basis_functions_only     = vscf_basis_functions_only
  this%maximum_weighted_displacement = maximum_weighted_displacement
  this%frequency_of_max_displacement = frequency_of_max_displacement
  
  this%subspace_qpoint_stars = generate_subspace_qpoint_stars( &
              & subspaces         = degenerate_subspaces,      &
              & modes             = complex_modes,             &
              & qpoints           = qpoints,                   &
              & symmetries        = structure%symmetries,      &
              & max_power         = potential_expansion_order, &
              & conserve_momentum = maximum_coupling_order==1  )
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
  
  ! Generate all sets of coupled subspaces, up to maximum_coupling_order.
  call print_line('Generating couplings between subspaces.')
  subspace_couplings = generate_coupled_subspaces( degenerate_subspaces,  &
                                                 & maximum_coupling_order )
  
  ! Load anharmonic data into container.
  this = AnharmonicData( structure,                                      &
                       & interpolated_supercell%supercell,               &
                       & interpolated_supercell%qpoints,                 &
                       & interpolated_supercell%complex_modes,           &
                       & interpolated_supercell%real_modes,              &
                       & degenerate_subspaces,                           &
                       & degenerate_symmetries,                          &
                       & subspace_couplings,                             &
                       & maximum_coupling_order,                         &
                       & potential_expansion_order,                      &
                       & vscf_basis_functions_only,                      &
                       & max_displacement%maximum_weighted_displacement, &
                       & max_displacement%frequency_of_max_displacement  )
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
  integer                               :: maximum_coupling_order
  integer                               :: potential_expansion_order
  logical                               :: vscf_basis_functions_only
  real(dp)                              :: maximum_weighted_displacement
  real(dp)                              :: frequency_of_max_displacement
  
  type(StringArray), allocatable :: sections(:)
  
  type(String) :: separator
  
  integer :: i
  
  separator = repeat('=',70)
  
  select type(this); type is(AnharmonicData)
    ! Split the file into sections.
    sections = split_into_sections(input, separating_line=separator)
    
    ! Trim the header line from each section.
    do i=1,size(sections)
      sections(i)%strings = sections(i)%strings(2:)
    enddo
    
    ! Parse each section into the relevant variable.
    structure = StructureData(sections(1))
    anharmonic_supercell = StructureData(sections(2))
    qpoints = QpointData(split_into_sections(sections(3)))
    complex_modes = ComplexMode(split_into_sections(sections(4)))
    real_modes = RealMode(split_into_sections(sections(5)))
    degenerate_subspaces = DegenerateSubspace(split_into_sections(sections(6)))
    degenerate_symmetries = &
       & DegenerateSymmetry(split_into_sections(sections(7)))
    subspace_couplings = SubspaceCoupling(sections(8)%strings)
    maximum_coupling_order = int(sections(9)%strings(1))
    potential_expansion_order = int(sections(10)%strings(1))
    vscf_basis_functions_only = lgcl(sections(11)%strings(1))
    maximum_weighted_displacement = dble(sections(12)%strings(1))
    frequency_of_max_displacement = dble(sections(13)%strings(1))
    
    ! Construct the output.
    this = AnharmonicData( structure,                     &
                         & anharmonic_supercell,          &
                         & qpoints,                       &
                         & complex_modes,                 &
                         & real_modes,                    &
                         & degenerate_subspaces,          &
                         & degenerate_symmetries,         &
                         & subspace_couplings,            &
                         & maximum_coupling_order,        &
                         & potential_expansion_order,     &
                         & vscf_basis_functions_only,     &
                         & maximum_weighted_displacement, &
                         & frequency_of_max_displacement  )
  class default
    call err()
  end select
end procedure

module procedure write_AnharmonicData
  type(String) :: separator
  
  separator = repeat('=',70)
  
  select type(this); type is(AnharmonicData)
    output = [ str('Structure'),                                    &
             & str(this%structure),                                 &
             & separator,                                           &
             & str('Anharmonic Supercell'),                         &
             & str(this%anharmonic_supercell),                      &
             & separator,                                           &
             & str('q-points'),                                     &
             & str(this%qpoints, separating_line=''),               &
             & separator,                                           &
             & str('Complex modes'),                                &
             & str(this%complex_modes, separating_line=''),         &
             & separator,                                           &
             & str('Real modes'),                                   &
             & str(this%real_modes, separating_line=''),            &
             & separator,                                           &
             & str('Degenerate subspaces'),                         &
             & str(this%degenerate_subspaces, separating_line=''),  &
             & separator,                                           &
             & str('Degenerate symmetries'),                        &
             & str(this%degenerate_symmetries, separating_line=''), &
             & separator,                                           &
             & str('Subspace couplings'),                           &
             & str(this%subspace_couplings),                        &
             & separator,                                           &
             & str('maximum coupling order'),                       &
             & str(this%maximum_coupling_order),                    &
             & separator,                                           &
             & str('Potential expansion order'),                    &
             & str(this%potential_expansion_order),                 &
             & separator,                                           &
             & str('VSCF basis functions only'),                    &
             & str(this%vscf_basis_functions_only),                 &
             & separator,                                           &
             & str('Maximum weighted displacement'),                &
             & str(this%maximum_weighted_displacement),             &
             & separator,                                           &
             & str('Frequency of maximum displacement'),            &
             & str(this%frequency_of_max_displacement)              ]
       
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
