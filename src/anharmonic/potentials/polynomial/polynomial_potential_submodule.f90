submodule (caesar_polynomial_potential_module) caesar_polynomial_potential_submodule
  use caesar_polynomial_module
contains

module procedure new_PolynomialPotential
  this%potential_expansion_order_ = anharmonic_data%potential_expansion_order
end procedure

module procedure new_PolynomialPotential_PotentialData
  select type(input); type is(PolynomialPotential)
    this = input
  type is(PotentialPointer)
    this = new_PolynomialPotential_PotentialData(input%potential())
  class default
    call err()
  end select
end procedure

module procedure new_PolynomialPotential_BasisFunctions
  this%potential_expansion_order_ = potential_expansion_order
  this%reference_energy_          = reference_energy
  this%basis_functions_           = basis_functions
end procedure

module procedure representation_PolynomialPotential
  output = 'polynomial'
end procedure

module procedure generate_sampling_points_PolynomialPotential
  ! Variables used when generating sampling points.
  type(SubspaceQpointStars),    allocatable :: subspace_qpoint_stars(:)
  type(CouplingBasisFunctions), allocatable :: basis_functions(:)
  type(SamplingPoints),         allocatable :: sampling_points(:)
  
  ! Supercell variables.
  type(RealMode),   allocatable :: sampling_point_modes(:)
  type(QpointData), allocatable :: sampling_point_qpoints(:)
  type(IntMatrix)               :: supercell_matrix
  type(StructureData)           :: supercell
  
  ! Displacement data.
  type(RealModeDisplacement)      :: sampling_point
  type(RealModeDisplacement)      :: transformed_sampling_point
  type(VscfRvectors), allocatable :: vscf_rvectors(:)
  type(CartesianDisplacement)     :: displacement
  type(StructureData)             :: displaced_structure
  
  ! Directories and files.
  type(String) :: equilibrium_dir
  type(String) :: coupling_dir
  type(OFile)  :: basis_functions_file
  type(OFile)  :: sampling_points_file
  type(String) :: sampling_dir
  type(OFile)  :: supercell_file
  type(OFile)  :: vscf_rvectors_file
  type(String) :: vscf_rvectors_dir
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  if (.not. use_forces) then
    call print_line(ERROR//': The polynomial potential cannot currently not &
       &use forces.')
    call err()
  elseif (use_hessians) then
    call print_line(ERROR//': The polynomial potential cannot currently &
       &use hessians.')
    call err()
  endif
  
  ! --------------------------------------------------
  ! Generate basis functions and sampling points.
  ! --------------------------------------------------
  
  subspace_qpoint_stars = generate_subspace_qpoint_stars(                  &
     & subspaces              = anharmonic_data%degenerate_subspaces,      &
     & modes                  = anharmonic_data%complex_modes,             &
     & qpoints                = anharmonic_data%qpoints,                   &
     & qpoint_symmetry_groups = anharmonic_data%qpoint_symmetry_groups,    &
     & max_power              = anharmonic_data%potential_expansion_order, &
     & max_qpoint_coupling    = anharmonic_data%max_qpoint_coupling,       &
     & conserve_momentum      = anharmonic_data%max_subspace_coupling==1   &
     &                     .or. anharmonic_data%vscf_basis_functions_only  )
  
  ! Loop over subspace couplings, generating basis functions and sampling
  !    points for each.
  allocate( basis_functions(size(anharmonic_data%subspace_couplings)), &
          & sampling_points(size(anharmonic_data%subspace_couplings)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(anharmonic_data%subspace_couplings)
    ! Generate basis functions at each coupling,
    !    and the sampling points from which the basis function coefficients
    !    can be constructed.
    call print_line('Generating sampling points in subspace coupling '// &
       & i//' of '//size(anharmonic_data%subspace_couplings)//': '// &
       & anharmonic_data%subspace_couplings(i))
    basis_functions(i) = generate_basis_functions(  &
       & anharmonic_data%subspace_couplings(i),     &
       & this%potential_expansion_order_,           &
       & anharmonic_data%max_subspace_coupling,     &
       & anharmonic_data%max_qpoint_coupling,       &
       & anharmonic_data%structure,                 &
       & anharmonic_data%complex_modes,             &
       & anharmonic_data%real_modes,                &
       & anharmonic_data%qpoints,                   &
       & anharmonic_data%degenerate_subspaces,      &
       & anharmonic_data%degenerate_symmetries,     &
       & anharmonic_data%vscf_basis_functions_only, &
       & anharmonic_data%qpoint_symmetry_groups,    &
       & subspace_qpoint_stars,                     &
       & logfile                                    )
    call print_line( 'Coupling contains '//size(basis_functions(i))// &
                   & ' basis functions.' )
    
    sampling_points(i) = generate_sampling_points( &
       & anharmonic_data%subspace_couplings(i),    &
       & basis_functions(i)%basis_functions(),     &
       & this%potential_expansion_order_,          &
       & anharmonic_data%max_displacement,         &
       & anharmonic_data%real_modes,               &
       & energy_to_force_ratio                     )
  enddo
  
  ! --------------------------------------------------
  ! Write out basis functions and sampling points.
  ! --------------------------------------------------
  
  ! Write out sampling point at equilibrium.
  equilibrium_dir = sampling_points_dir//'/equilibrium'
  call calculation_writer%write_calculation( anharmonic_data%structure, &
                                           & equilibrium_dir            )
  
  ! Write out all other sampling points.
  do i=1,size(sampling_points)
    if (modulo(i,max(size(anharmonic_data%subspace_couplings)/100,1))==0) then
      call print_line('Generating calculations in subspace coupling '// &
         & i//' of '//size(anharmonic_data%subspace_couplings)//'.')
    endif
    ! Make a directory for each coupling.
    coupling_dir = sampling_points_dir// &
       & '/coupling_'//left_pad(i,str(size(sampling_points)))
    call mkdir(coupling_dir)
    
    ! Write basis functions to file.
    basis_functions_file = OFile(coupling_dir//'/basis_functions.dat')
    call basis_functions_file%print_lines(basis_functions(i))
    
    ! Write sampling points to file.
    sampling_points_file = OFile(coupling_dir//'/sampling_points.dat')
    call sampling_points_file%print_lines(sampling_points(i))
    
    do j=1,size(sampling_points(i))
      ! Select the sampling point for clarity.
      sampling_point = sampling_points(i)%points(j)
      
      ! Make a directory for each sampling point.
      sampling_dir = coupling_dir// &
         & '/sampling_point_'//left_pad(j,str(size(sampling_points(i))))
      call mkdir(sampling_dir)
      
      ! Construct a supercell for each sampling point.
      sampling_point_modes = select_modes( sampling_point%vectors,    &
                                         & anharmonic_data%real_modes )
      sampling_point_qpoints = select_qpoints( sampling_point_modes,   &
                                             & anharmonic_data%qpoints )
      supercell_matrix = construct_supercell_matrix( &
                         & sampling_point_qpoints,   &
                         & anharmonic_data%structure )
      supercell = construct_supercell( anharmonic_data%structure, &
                                     & supercell_matrix           )
      
      ! Write out the supercell.
      supercell_file = OFile(sampling_dir//'/structure.dat')
      call supercell_file%print_lines(supercell)
      
      ! Construct VSCF R-vectors.
      vscf_rvectors = construct_vscf_rvectors( sampling_point,             &
                                             & supercell,                  &
                                             & anharmonic_data%real_modes, &
                                             & anharmonic_data%qpoints     )
      vscf_rvectors_file = OFile(sampling_dir//'/vscf_rvectors.dat')
      call vscf_rvectors_file%print_lines(vscf_rvectors,separating_line='')
      
      do k=1,size(vscf_rvectors)
        ! Transform the sampling point by the VSCF R-vector.
        transformed_sampling_point = vscf_rvectors(k)%transform( &
                                   & sampling_point,             &
                                   & anharmonic_data%real_modes, &
                                   & anharmonic_data%qpoints     )
        
        ! Construct displaced structure.
        displacement = CartesianDisplacement( transformed_sampling_point, &
                                            & supercell,                  &
                                            & anharmonic_data%real_modes, &
                                            & anharmonic_data%qpoints     )
        displaced_structure = displace_structure(supercell,displacement)
        
        ! Create directory and structure files for displaced structure.
        vscf_rvectors_dir = sampling_dir// &
           & '/vscf_rvectors_'//left_pad(k,str(size(vscf_rvectors)))
        call calculation_writer%write_calculation( displaced_structure, &
                                                 & vscf_rvectors_dir    )
      enddo
    enddo
  enddo
end procedure

module procedure generate_potential_PolynomialPotential
  ! Basis functions and sampling points corresponding to the un-displaced
  !    structure and the constant basis function.
  type(BasisFunction)         :: constant_basis_function
  type(RealModeDisplacement)  :: equilibrium_sampling_point
  type(SampleResult)          :: equilibrium_sample_result
  
  ! Basis functions and sampling points at each coupling.
  type(SubspaceCoupling)       :: coupling
  type(String)                 :: coupling_directory
  type(CouplingBasisFunctions) :: basis_functions
  type(SamplingPoints)         :: sampling_points
  type(SampleResults)          :: sample_results
  
  ! Temporary variables.
  integer :: i
  
  ! Generate the potential at zero displacement.
  call print_line('Generating potential at zero displacement.')
  constant_basis_function = generate_constant_basis_function()
  equilibrium_sampling_point = generate_equilibrium_sampling_point()
  equilibrium_sample_result = read_equilibrium_sample_result( &
                                & equilibrium_sampling_point, &
                                & sampling_points_dir,        &
                                & anharmonic_data,            &
                                & calculation_reader          )
  this%reference_energy_ = equilibrium_sample_result%energy
  
  this%basis_functions_ = CouplingBasisFunctions( &
             & anharmonic_data%subspace_couplings )
  
  do i=1,size(anharmonic_data%subspace_couplings)
    call print_line('Fitting potential in subspace coupling '//i//' of ' // &
       &size(anharmonic_data%subspace_couplings)//', containing &
       &subspaces '//anharmonic_data%subspace_couplings(i)%ids()//'.')
    
    coupling = anharmonic_data%subspace_couplings(i)
    
    ! Generate directory names.
    coupling_directory = sampling_points_dir//'/coupling_'//left_pad( &
                      & i,                                            &
                      & str(size(anharmonic_data%subspace_couplings)) )
    
    ! Read in basis functions.
    basis_functions = read_basis_functions(coupling_directory)
    
    ! Read in sampling points.
    sampling_points = read_sampling_points(coupling_directory)
    
    ! Read in electronic structure data.
    sample_results = read_sample_results( sampling_points,    &
                                        & coupling_directory, &
                                        & anharmonic_data,    &
                                        & calculation_reader  )
    
    ! Fit basis function coefficients.
    call basis_functions%fit_coefficients( sampling_points%points,      &
                                         & sample_results%results,      &
                                         & anharmonic_data%real_modes,  &
                                         & weighted_energy_force_ratio, &
                                         & this                         )
    
    ! Add fitted basis functions to potential.
    this%basis_functions_(i) = basis_functions
  enddo
end procedure

module procedure generate_stress_PolynomialPotential
  type(SubspaceQpointStars),    allocatable :: subspace_qpoint_stars(:)
  
  ! Electronic structure results corresponding to the un-displaced structure.
  type(CouplingStressBasisFunctions) :: zero_basis_functions(0)
  type(RealModeDisplacement)         :: equilibrium_sampling_point
  type(SampleResult)                 :: equilibrium_sample_result
  
  ! Stress basis functions.
  type(CouplingStressBasisFunctions), allocatable :: basis_functions(:)
  
  ! Electronic structure results.
  type(String)         :: coupling_directory
  type(SamplingPoints) :: sampling_points
  type(SampleResults)  :: sample_results
  
  ! The stress.
  type(PolynomialStress) :: stress
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  subspace_qpoint_stars = generate_subspace_qpoint_stars(                  &
     & subspaces              = anharmonic_data%degenerate_subspaces,      &
     & modes                  = anharmonic_data%complex_modes,             &
     & qpoints                = anharmonic_data%qpoints,                   &
     & qpoint_symmetry_groups = anharmonic_data%qpoint_symmetry_groups,    &
     & max_power              = anharmonic_data%potential_expansion_order, &
     & max_qpoint_coupling    = anharmonic_data%max_qpoint_coupling,       &
     & conserve_momentum      = anharmonic_data%max_subspace_coupling==1   &
     &                     .or. anharmonic_data%vscf_basis_functions_only  )
  
  ! Generate the stress tensor at zero displacement.
  call print_line('Generating stress at zero displacement.')
  equilibrium_sampling_point = generate_equilibrium_sampling_point()
  equilibrium_sample_result = read_equilibrium_sample_result( &
                                & equilibrium_sampling_point, &
                                & sampling_points_dir,        &
                                & anharmonic_data,            &
                                & calculation_reader          )
  stress = PolynomialStress( stress_expansion_order,             &
                           & equilibrium_sample_result%stress(), &
                           & zero_basis_functions                )
  
  ! Generate basis functions.
  allocate( basis_functions(size(stress_subspace_coupling)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(basis_functions)
    call print_line('Fitting stress in subspace coupling '//i//' of ' // &
       &size(stress_subspace_coupling)//', containing &
       &subspaces '//stress_subspace_coupling(i)%ids()//'.')
    basis_functions(i) = generate_stress_basis_functions( &
                & stress_subspace_coupling(i),            &
                & stress_expansion_order,                 &
                & anharmonic_data%max_subspace_coupling,  &
                & anharmonic_data%max_qpoint_coupling,    &
                & anharmonic_data%structure,              &
                & anharmonic_data%complex_modes,          &
                & anharmonic_data%qpoints,                &
                & anharmonic_data%degenerate_subspaces,   &
                & anharmonic_data%degenerate_symmetries,  &
                & vscf_basis_functions_only,              &
                & anharmonic_data%qpoint_symmetry_groups, &
                & subspace_qpoint_stars,                  &
                & logfile                                 )
    
    j = first(anharmonic_data%subspace_couplings==stress_subspace_coupling(i))
    coupling_directory = sampling_points_dir//'/coupling_'//left_pad( &
                      & j,                                            &
                      & str(size(anharmonic_data%subspace_couplings)) )
    
    ! Read in all sampling points and sample results.
    sampling_points = read_sampling_points(coupling_directory)
    
    sample_results = read_sample_results( sampling_points,    &
                                        & coupling_directory, &
                                        & anharmonic_data,    &
                                        & calculation_reader  )
    
    ! Fit basis functions.
    j = first(stress_subspace_coupling==basis_functions(i)%coupling)
    
    call basis_functions(i)%fit_coefficients( &
                    & sampling_points%points, &
                    & sample_results%results, &
                    & stress                  )
  enddo
  
  ! Assemble output.
  stress = PolynomialStress( stress_expansion_order,             &
                           & equilibrium_sample_result%stress(), &
                           & basis_functions                     )
  output = StressPointer(stress)
end procedure

module procedure read_sampling_points
  type(IFile) :: sampling_points_file
  
  sampling_points_file = IFile(coupling_directory//'/sampling_points.dat')
  output = SamplingPoints(sampling_points_file%lines())
end procedure

module procedure generate_equilibrium_sampling_point
  type(RealSingleDisplacement) :: zero_displacement(0)
  
  output = RealModeDisplacement(zero_displacement)
end procedure

module procedure read_basis_functions
  type(IFile) :: basis_functions_file
  
  basis_functions_file = IFile(coupling_directory//'/basis_functions.dat')
  output = CouplingBasisFunctions(basis_functions_file%lines())
end procedure

module procedure generate_constant_basis_function
  type(ComplexMonomial) :: constant_complex_monomial
  
  type(ComplexUnivariate) :: zero_complex(0)
  
  constant_complex_monomial = ComplexMonomial( &
         & coefficient = (1.0_dp,0.0_dp),      &
         & modes       = zero_complex          )
  output = BasisFunction(                                                     &
     & complex_representation = ComplexPolynomial([constant_complex_monomial]))
end procedure

module procedure read_sample_results
  type(StructureData)                    :: supercell
  type(CartesianDisplacement)            :: displacement
  type(VscfRvectors),        allocatable :: vscf_rvectors(:)
  type(ElectronicStructure), allocatable :: calculations(:)
  
  ! Files and directories.
  type(String) :: sampling_directory
  type(IFile)  :: supercell_file
  type(IFile)  :: vscf_rvectors_file
  type(String) :: vscf_rvectors_directory
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  allocate( output%results(size(sampling_points)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(sampling_points)
    sampling_directory = coupling_directory// &
       & '/sampling_point_'//left_pad(i,str(size(sampling_points)))
    
    ! Read in supercell and VSCF R-vectors.
    supercell_file = IFile(sampling_directory//'/structure.dat')
    supercell = StructureData(supercell_file%lines())
    
    displacement = CartesianDisplacement( sampling_points%points(i),  &
                                        & supercell,                  &
                                        & anharmonic_data%real_modes, &
                                        & anharmonic_data%qpoints     )
    
    vscf_rvectors_file = IFile(sampling_directory//'/vscf_rvectors.dat')
    vscf_rvectors = VscfRvectors(vscf_rvectors_file%sections())
    
    ! Read in electronic structure calculations.
    allocate( calculations(size(vscf_rvectors)), &
            & stat=ialloc); call err(ialloc)
    do j=1,size(vscf_rvectors)
      vscf_rvectors_directory = sampling_directory// &
         & '/vscf_rvectors_'//left_pad(j,str(size(vscf_rvectors)))
      calculations(j) = calculation_reader%read_calculation( &
                                  & vscf_rvectors_directory, &
                                  & displacement             )
    enddo
    
    ! Average electronic structure across VSCF R-vectors, and convert
    !    to correct normalisation and real mode co-ordinates.
    output%results(i) = SampleResult( vscf_rvectors,              &
                                    & calculations,               &
                                    & supercell,                  &
                                    & anharmonic_data%real_modes, &
                                    & anharmonic_data%qpoints,    &
                                    & anharmonic_data             )
    
    deallocate(calculations, stat=ialloc); call err(ialloc)
  enddo
end procedure

module procedure read_equilibrium_sample_result
  type(CartesianDisplacement) :: displacement
  type(ElectronicStructure)   :: electronic_structure
  
  ! The vector of length 0 as a CartesianDisplacement.
  displacement = CartesianDisplacement( &
          & equilibrium_sampling_point, &
          & anharmonic_data%structure,  &
          & anharmonic_data%real_modes, &
          & anharmonic_data%qpoints     )
  
  electronic_structure = calculation_reader%read_calculation( &
                 & sampling_points_directory//'/equilibrium', &
                 & displacement                               )
  
  output = SampleResult( electronic_structure,       &
                       & anharmonic_data%structure,  &
                       & anharmonic_data%real_modes, &
                       & anharmonic_data%qpoints,    &
                       & anharmonic_data             )
end procedure

module procedure zero_energy_PolynomialPotential
  this%reference_energy_ = this%reference_energy_ - this%undisplaced_energy()
end procedure

module procedure add_constant_PolynomialPotential
  this%reference_energy_ = this%reference_energy_ + input
end procedure

module procedure optimise_subspace_potential_PolynomialPotential
  if (size(this%basis_functions_)/=1) then
    call print_line(CODE_ERROR//': Calling optimise_subspace_potential &
       &on a potential with more than one coupling.')
    call err()
  endif
  
  this%reference_energy_ = this%reference_energy_ &
                       & + this%basis_functions_(1)%undisplaced_energy()
  
  call this%basis_functions_(1)%optimise( subspace,                 &
                                        & subspace_basis,           &
                                        & old_subspace_potential,   &
                                        & anharmonic_data           )
end procedure

module procedure energy_RealModeDisplacement_PolynomialPotential
  output = this%reference_energy_ &
       & + sum(this%basis_functions_%energy(displacement))
end procedure

module procedure energy_ComplexModeDisplacement_PolynomialPotential
  output = this%reference_energy_ &
       & + sum(this%basis_functions_%energy(displacement))
end procedure

module procedure force_RealModeDisplacement_PolynomialPotential
  if (size(this%basis_functions_)>0) then
    output = sum(this%basis_functions_%force(displacement))
  else
    output = RealModeForce(RealSingleForce(displacement%vectors%id,0.0_dp))
  endif
end procedure

module procedure force_ComplexModeDisplacement_PolynomialPotential
  if (size(this%basis_functions_)>0) then
    output = sum(this%basis_functions_%force(displacement))
  else
    output = ComplexModeForce(ComplexSingleForce( displacement%vectors%id, &
                                                & (0.0_dp,0.0_dp)          ))
  endif
end procedure

module procedure braket_SubspaceBraKet_PolynomialPotential
  integer :: i
  
  ! Integrate the reference energy (N.B. <i|e|j> = e<i|j> if e is a scalar.).
  this%reference_energy_ = this%reference_energy_ &
                       & * braket%inner_product(anharmonic_data)
  
  ! Integrate each basis function between the bra and the ket.
  do i=1,size(this%basis_functions_)
    call this%basis_functions_(i)%braket( braket,         &
                                        & whole_subspace, &
                                        & anharmonic_data )
  enddo
  
  ! Simplify the potential.
  if (set_default(whole_subspace,.true.)) then
    call this%simplify()
  endif
end procedure

module procedure braket_BasisState_PolynomialPotential
  integer :: i
  
  ! Integrate the reference energy (N.B. <i|e|j> = e<i|j> if e is a scalar.).
  this%reference_energy_ = this%reference_energy_                         &
                       & * subspace_basis%inner_product( bra,            &
                       &                                 ket,            &
                       &                                 subspace,       &
                       &                                 anharmonic_data )
  
  ! Integrate each basis function between the bra and the ket.
  do i=1,size(this%basis_functions_)
    call this%basis_functions_(i)%braket( bra,            &
                                        & ket,            &
                                        & subspace,       &
                                        & subspace_basis, &
                                        & whole_subspace, &
                                        & anharmonic_data )
  enddo
  
  ! Simplify the potential.
  if (set_default(whole_subspace,.true.)) then
    call this%simplify()
  endif
end procedure

module procedure braket_BasisStates_PolynomialPotential
  integer :: i
  
  ! Integrate each basis function between the bra and the ket.
  do i=1,size(this%basis_functions_)
    call this%basis_functions_(i)%braket( states,         &
                                        & subspace,       &
                                        & subspace_basis, &
                                        & whole_subspace, &
                                        & anharmonic_data )
  enddo
  
  ! Simplify the potential.
  if (set_default(whole_subspace,.true.)) then
    call this%simplify()
  endif
end procedure

module procedure simplify_PolynomialPotential
  logical, allocatable :: to_remove(:)
  
  integer :: i,j
  
  to_remove = [(.false., i=1, size(this%basis_functions_))]
  
  do i=1,size(this%basis_functions_)
    if (size(this%basis_functions_(i)%coupling)==0) then
      this%reference_energy_ =      &
         &   this%reference_energy_ &
         & + this%basis_functions_(i)%undisplaced_energy()
      to_remove(i) = .true.
    else
      do j=1,i-1
        if (    this%basis_functions_(i)%coupling &
           & == this%basis_functions_(j)%coupling ) then
          if (to_remove(j)) then
            call err()
          endif
          call this%basis_functions_(j)%append( &
                     & this%basis_functions_(i) )
          to_remove(i) = .true.
          exit
        endif
      enddo
    endif
  enddo
  
  ! Remove constant and duplicate terms.
  this%basis_functions_ = this%basis_functions_(filter(.not.to_remove))
end procedure

module procedure harmonic_expectation_PolynomialPotential
  output = this%reference_energy_                                          &
       & + sum(this%basis_functions_%harmonic_expectation( frequency,      &
       &                                                   thermal_energy, &
       &                                                   supercell_size, &
       &                                                   anharmonic_data ))
end procedure

module procedure potential_energy_SubspaceBraKet_PolynomialPotential
  integer :: i
  
  output = this%reference_energy_                               &
       & * braket%inner_product(anharmonic_data)                &
       & + sum([( this%basis_functions_(i)%potential_energy(    &
       &                                     braket,            &
       &                                     anharmonic_data ), &
       &          i=1,                                          &
       &          size(this%basis_functions_)                   )])
end procedure

module procedure potential_energy_BasisState_PolynomialPotential
  integer :: i
  
  output = this%reference_energy_                                      &
       & * subspace_basis%inner_product( bra,                          &
       &                                 ket,                          &
       &                                 subspace,                     &
       &                                 anharmonic_data )             &
       & + sum([( this%basis_functions_(i)%potential_energy(           &
       &                                     bra,                      &
       &                                     ket,                      &
       &                                     subspace,                 &
       &                                     subspace_basis,           &
       &                                     anharmonic_data ),        &
       &          i=1,                                                 &
       &          size(this%basis_functions_)                   )])
end procedure

module procedure coefficients_PolynomialPotential
  integer :: i
  
  output = [( this%basis_functions_(i)%coefficients(            &
            &    anharmonic_data%potential_expansion_order-1 ), &
            & i=1,                                              &
            & size(this%basis_functions_)                       )]
end procedure

module procedure set_coefficients_PolynomialPotential
  integer :: no_coefficients
  
  integer :: i,j
  
  j = 0
  do i=1,size(this%basis_functions_)
    no_coefficients = this%basis_functions_(i)%no_coefficients( &
                  & anharmonic_data%potential_expansion_order-1 )
    call this%basis_functions_(i)%set_coefficients(  &
       & coefficients(j+1:j+no_coefficients),        &
       & anharmonic_data%potential_expansion_order-1 )
    j = j+no_coefficients
  enddo
end procedure

module procedure all_basis_functions_PolynomialPotential
  integer :: i
  
  output = [( this%basis_functions_(i)%all_basis_functions(), &
            & i=1,                                            &
            & size(this%basis_functions_)                     )]
end procedure

module procedure variable_basis_functions_PolynomialPotential
  integer :: i
  
  output = [( this%basis_functions_(i)%variable_basis_functions(    &
            &        anharmonic_data%potential_expansion_order-1 ), &
            & i=1,                                                  &
            & size(this%basis_functions_)                           )]
end procedure

module procedure can_be_interpolated_PolynomialPotential
  output = .true.
end procedure

module procedure interpolate_potential_PolynomialPotential
  type(SubspaceQpointStars), allocatable :: subspace_qpoint_stars(:)
  
  type(PolynomialPotential) :: potential_
  
  type(PolynomialInterpolator) :: interpolator
  
  integer :: i,j,ialloc
  
  subspace_qpoint_stars = generate_subspace_qpoint_stars(              &
     & subspaces              =                                        &
     &         interpolated_anharmonic_data%degenerate_subspaces,      &
     & modes                  =                                        &
     &         interpolated_anharmonic_data%complex_modes,             &
     & qpoints                =                                        &
     &         interpolated_anharmonic_data%qpoints,                   &
     & qpoint_symmetry_groups =                                        &
     &         interpolated_anharmonic_data%qpoint_symmetry_groups,    &
     & max_power              =                                        &
     &         interpolated_anharmonic_data%potential_expansion_order, &
     & max_qpoint_coupling    =                                        &
     &         interpolated_anharmonic_data%max_qpoint_coupling,       &
     & conserve_momentum      =                                        &
     &         interpolated_anharmonic_data%max_subspace_coupling==1   &
     &    .or. interpolated_anharmonic_data%vscf_basis_functions_only  )
  
  potential_ = PolynomialPotential(potential)
  
  this%reference_energy_ =                                           &
     & ( potential_%reference_energy_                                &
     & * interpolated_anharmonic_data%anharmonic_supercell%sc_size ) &
     & / anharmonic_data%anharmonic_supercell%sc_size
  
  allocate( this%basis_functions_(                                      &
          &    size(interpolated_anharmonic_data%subspace_couplings) ), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(this%basis_functions_)
    call print_line('Generating basis functions for interpolated potential &
       &in subspace coupling '//i//' of '//size(this%basis_functions_)//'.' )
    ! Generate basis functions at each coupling,
    !    and the sampling points from which the basis function coefficients
    !    can be constructed.
    this%basis_functions_(i) = generate_basis_functions(         &
       & interpolated_anharmonic_data%subspace_couplings(i),     &
       & this%potential_expansion_order_,                        &
       & interpolated_anharmonic_data%max_subspace_coupling,     &
       & interpolated_anharmonic_data%max_qpoint_coupling,       &
       & interpolated_anharmonic_data%structure,                 &
       & interpolated_anharmonic_data%complex_modes,             &
       & interpolated_anharmonic_data%real_modes,                &
       & interpolated_anharmonic_data%qpoints,                   &
       & interpolated_anharmonic_data%degenerate_subspaces,      &
       & interpolated_anharmonic_data%degenerate_symmetries,     &
       & interpolated_anharmonic_data%vscf_basis_functions_only, &
       & interpolated_anharmonic_data%qpoint_symmetry_groups,    &
       & subspace_qpoint_stars,                                  &
       & logfile                                                 )
    call this%basis_functions_(i)%zero_coefficients()
  enddo
  
  ! Construct the polynomial interpolator.
  interpolator = PolynomialInterpolator(                           &
     & min_images                   = anharmonic_min_images,       &
     & anharmonic_data              = anharmonic_data,             &
     & interpolated_anharmonic_data = interpolated_anharmonic_data )
  
  do i=1,size(this%basis_functions_)
    call print_line( 'Interpolating coefficients in subspace coupling '//i// &
                   & ' of '//size(this%basis_functions_)//'.'                )
    do j=1,size(potential_%basis_functions_)
      call this%basis_functions_(i)%add_interpolated_contribution( &
                                 & potential_%basis_functions_(j), &
                                 & interpolator                    )
    enddo
    
    call this%basis_functions_(i)%add_harmonic_contribution( &
                            & difference_dynamical_matrices, &
                            & interpolated_anharmonic_data   )
  enddo
end procedure

module procedure interpolate_coefficient_PolynomialPotential
  output = sum(this%basis_functions_%interpolate_coefficient( monomial,    &
                                                            & interpolator ))
end procedure

module procedure calculate_dynamical_matrices_PolynomialPotential
  integer :: i
  
  output = [( DynamicalMatrix(anharmonic_data%structure%no_atoms), &
            & i=1,                                                 &
            & size(qpoints)                                        )]
  
  do i=1,size(this%basis_functions_)
    output = output                                                 &
         & + this%basis_functions_(i)%calculate_dynamical_matrices( &
         &                                         qpoints,         &
         &                                         thermal_energy,  &
         &                                         subspaces,       &
         &                                         subspace_bases,  &
         &                                         subspace_states, &
         &                                         anharmonic_data  )
  enddo
end procedure

module procedure energy_correction_PolynomialPotential
  integer :: i
  
  output = 0
  do i=1,size(this%basis_functions_)
    output = output + this%basis_functions_(i)%energy_correction( &
                                               & subspaces,       &
                                               & subspace_bases,  &
                                               & subspace_states, &
                                               & anharmonic_data  )
  enddo
end procedure

module procedure expansion_order_PolynomialPotential
  output = this%potential_expansion_order_
end procedure

module procedure read_PolynomialPotential
  integer                                   :: potential_expansion_order
  real(dp)                                  :: reference_energy
  type(CouplingBasisFunctions), allocatable :: basis_functions(:)
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(PolynomialPotential)
    line = split_line(input(1))
    potential_expansion_order = int(line(3))
    
    line = split_line(input(2))
    reference_energy = dble(line(3))
    
    basis_functions = CouplingBasisFunctions(split_into_sections( &
                                 & input(5:),                     &
                                 & separating_line=repeat('=',50) ))
    
    this = PolynomialPotential( potential_expansion_order, &
                              & reference_energy,          &
                              & basis_functions            )
  class default
    call err()
  end select
end procedure

module procedure write_PolynomialPotential
  select type(this); type is(PolynomialPotential)
    output = [ 'Expansion order: '//this%potential_expansion_order_,      &
             & 'Reference energy: '//this%reference_energy_,              &
             & str('Basis functions:'),                                   &
             & str(''),                                                   &
             & str(this%basis_functions_, separating_line=repeat('=',50)) ]
  class default
    call err()
  end select
end procedure

module procedure new_PolynomialPotential_Strings
  call this%read(input)
end procedure

module procedure new_PolynomialPotential_StringArray
  this = PolynomialPotential(str(input))
end procedure
end submodule
