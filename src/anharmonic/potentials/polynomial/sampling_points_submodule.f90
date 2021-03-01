submodule (caesar_sampling_points_module) caesar_sampling_points_submodule
  use caesar_polynomial_module
contains

module procedure new_SamplingPoints
  this%points = points
end procedure

module procedure size_SamplingPoints
  output = size(this%points)
end procedure

module procedure construct_sample_matrix
  real(dp)              :: energy
  type(RealModeForce)   :: forces
  real(dp), allocatable :: weight
  
  integer :: dims
  
  integer :: i,j,ialloc
  
  dims = 1+size(modes)
  
  allocate( output( size(sampling_points)*dims,    &
          &         size(basis_functions)       ), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(basis_functions)
    do j=1,size(sampling_points)
      energy = basis_functions(i)%energy(sampling_points(j))
      forces = basis_functions(i)%force(sampling_points(j))
      if (present(sample_weights)) then
        weight = sample_weights(j)
      endif
      
      output((j-1)*dims+1:j*dims, i) = make_sample_vector( &
                                     & energy,             &
                                     & forces,             &
                                     & modes,              &
                                     & energy_force_ratio, &
                                     & weight              )
    enddo
  enddo
end procedure

module procedure make_sample_vector
  integer :: i,ialloc
  
  allocate(output(1+size(modes)), stat=ialloc); call err(ialloc)
  output(1) = energy / energy_force_ratio
  do i=1,size(modes)
    output(1+i) = force%force(modes(i))
  enddo
  
  if (present(weight)) then
    output = output * weight
  endif
end procedure

module procedure generate_sampling_points
  type(RealModeDisplacement), allocatable :: basis_points(:)
  logical,                    allocatable :: basis_points_used(:)
  real(dp),                   allocatable :: determinants(:)
  
  type(RealModeDisplacement), allocatable :: sampling_points(:)
  
  integer :: points_per_basis_function
  
  real(dp), allocatable :: sample_matrix(:,:)
  
  integer :: dims
  
  integer :: i,j,k,l,ialloc
  
  dims = 1+size(real_modes)
  
  points_per_basis_function = 2
  
  allocate( sampling_points( size(basis_functions)       &
          &                * points_per_basis_function), &
          & sample_matrix( size(sampling_points)         &
          &              * dims,                         &
          &                size(basis_functions)         &
          &              * points_per_basis_function),   &
          & stat=ialloc); call err(ialloc)
  l = 0
  do i=1,size(basis_functions)
    sample_matrix(:l*dims,i:i) = construct_sample_matrix( &
                                 & basis_functions(i:i),  &
                                 & sampling_points(:l),   &
                                 & real_modes,            &
                                 & energy_to_force_ratio  )
    
    ! Generate the possible sampling points for this basis function.
    basis_points = generate_basis_points( basis_functions(i),            &
                                        & potential_expansion_order,     &
                                        & maximum_weighted_displacement, &
                                        & frequency_of_max_displacement, &
                                        & real_modes                     )
    basis_points_used = [(.false., j=1, size(basis_points))]
    determinants = [(0.0_dp, j=1, size(basis_points))]
    
    ! Find the sampling point which most resolves the first i basis functions.
    ! Repeat until points_per_basis_function sampling points have been found.
    do j=1,points_per_basis_function
      l = l+1
      do k=1,size(basis_points)
        if (.not. basis_points_used(k)) then
          sampling_points(l) = basis_points(k)
          
          sample_matrix((l-1)*dims+1:l*dims,:i) = construct_sample_matrix( &
                                                  & basis_functions(:i),   &
                                                  & sampling_points(l:l),  &
                                                  & real_modes,            &
                                                  & energy_to_force_ratio  )
          determinants(k) = calculate_average_eigenvalue( &
                         & mat(sample_matrix(:l*dims,:i)) )
        endif
      enddo
      k = maxloc(abs(determinants), 1, mask=.not.basis_points_used)
      sampling_points(l) = basis_points(k)
      basis_points_used(k) = .true.
      
      sample_matrix((l-1)*dims+1:l*dims,:i) = construct_sample_matrix( &
                                               & basis_functions(:i),  &
                                               & basis_points(k:k),    &
                                               & real_modes,           &
                                               & energy_to_force_ratio )
      
      if (i==size(basis_functions) .and. j==points_per_basis_function) then
        call print_line( 'Average |eigenvalue| of fitting matrix: '// &
                       & determinants(k)*energy_to_force_ratio//' (Ha).')
      endif
    enddo
  enddo
  
  output = SamplingPoints(sampling_points)
end procedure

module procedure generate_basis_points
  type(ComplexMonomial), allocatable :: monomials(:)
  
  type(RealSingleDisplacement), allocatable :: displacement(:)
  
  type(RealMode) :: mode
  real(dp)       :: frequency
  real(dp)       :: magnitude
  
  integer :: i,j,ialloc
  
  monomials = basis_function%terms()
  allocate(output(2*size(monomials)), stat=ialloc); call err(ialloc)
  do i=1,size(monomials)
    allocate(displacement(0), stat=ialloc); call err(ialloc)
    do j=1,size(monomials(i))
      if (monomials(i)%power(j)>0) then
        mode = modes(first(modes%id==monomials(i)%id(j)))
        frequency = max(mode%frequency, frequency_of_max_displacement)
        magnitude = maximum_weighted_displacement         &
                & * sqrt( monomials(i)%power(j)           &
                &       * frequency_of_max_displacement   &
                &       / frequency                     ) &
                & / potential_expansion_order
        displacement = [displacement, RealSingleDisplacement(mode, magnitude)]
      endif
      if (monomials(i)%id(j)/=monomials(i)%paired_id(j)) then
        if (monomials(i)%paired_power(j)>0) then
          mode = modes(first(modes%id==monomials(i)%paired_id(j)))
          frequency = max(mode%frequency, frequency_of_max_displacement)
          magnitude = maximum_weighted_displacement         &
                  & * sqrt( monomials(i)%paired_power(j)    &
                  &       * frequency_of_max_displacement   &
                  &       / frequency                     ) &
                  & / potential_expansion_order
          displacement = [ displacement,                           &
                         & RealSingleDisplacement(mode, magnitude) ]
        endif
      endif
    enddo
    output(2*i-1:2*i) = [  RealModeDisplacement(displacement), &
                        & -RealModeDisplacement(displacement)  ]
    deallocate(displacement, stat=ialloc); call err(ialloc)
  enddo
end procedure

module procedure calculate_average_eigenvalue
  type(RealQRDecomposition) :: qr
  real(dp), allocatable     :: evals(:)
  
  integer :: i
  
  qr = qr_decomposition(transpose(matrix)*matrix)
  
  ! The eigenvalues are the diagonal entries of the R matrix
  !    in a QR decomposition.
  evals = [(qr%r(i,i), i=1, size(qr%r,1))]
  
  evals = abs(evals)
  
  if (any(evals<1e-200_dp)) then
    output = 0.0_dp
    return
  endif
  
  output = product(evals**(1.0_dp/size(qr%r,1)))
end procedure

module procedure read_SamplingPoints
  select type(this); type is(SamplingPoints)
    this = SamplingPoints(RealModeDisplacement(split_into_sections(input)))
  class default
    call err()
  end select
end procedure

module procedure write_SamplingPoints
  select type(this); type is(SamplingPoints)
    output = str(this%points,separating_line='')
  class default
    call err()
  end select
end procedure

module procedure new_SamplingPoints_Strings
  call this%read(input)
end procedure

module procedure new_SamplingPoints_StringArray
  this = SamplingPoints(str(input))
end procedure
end submodule
