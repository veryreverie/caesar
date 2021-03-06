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
  type(RealMode), allocatable :: coupling_modes(:)
  
  type(RealModeDisplacement), allocatable :: basis_points(:)
  logical,                    allocatable :: basis_points_used(:)
  real(dp),                   allocatable :: determinants(:)
  
  type(RealModeDisplacement), allocatable :: sampling_points(:)
  
  integer :: points_per_basis_function
  
  real(dp), allocatable :: sample_matrix(:,:)
  real(dp), allocatable :: old_fitting_matrix(:,:)
  real(dp), allocatable :: new_fitting_matrix(:,:)
  
  type(RealVector), allocatable :: reduced_sampling_vector
  type(RealMatrix), allocatable :: reduced_sampling_matrix
  
  type(RealQRDecomposition) :: qr
  real(dp), allocatable     :: evals(:)
  
  integer :: dims
  
  integer :: i,j,k,l,m,ialloc
  
  coupling_modes = real_modes(                                   &
     & filter(real_modes%subspace_id.in.subspace_coupling%ids()) )
  
  dims = 1+size(coupling_modes)
  
  points_per_basis_function = 2
  
  allocate( sampling_points( size(basis_functions)            &
          &                * points_per_basis_function ),     &
          & sample_matrix( size(sampling_points)              &
          &              * dims,                              &
          &                size(basis_functions)              &
          &              * points_per_basis_function ),       &
          & old_fitting_matrix( size(basis_functions)         &
          &                   * points_per_basis_function,    &
          &                     size(basis_functions)         &
          &                   * points_per_basis_function  ), &
          & new_fitting_matrix( size(basis_functions)         &
          &                   * points_per_basis_function,    &
          &                     size(basis_functions)         &
          &                   * points_per_basis_function  ), &
          & stat=ialloc); call err(ialloc)
  l = 0
  do i=1,size(basis_functions)
    sample_matrix(:l*dims,i:i) = construct_sample_matrix( &
                                 & basis_functions(i:i),  &
                                 & sampling_points(:l),   &
                                 & coupling_modes,        &
                                 & energy_to_force_ratio  )
    reduced_sampling_vector = vec(sample_matrix(:l*dims,i))
    reduced_sampling_matrix = mat(sample_matrix(:l*dims,:i))
    old_fitting_matrix(i,:i) = dble( reduced_sampling_vector &
                                 & * reduced_sampling_matrix )
    old_fitting_matrix(:i-1,i) = old_fitting_matrix(i,:i-1)
    
    ! Generate the possible sampling points for this basis function.
    basis_points = generate_basis_points( basis_functions(i),        &
                                        & potential_expansion_order, &
                                        & max_displacement,          &
                                        & coupling_modes             )
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
                                                  & coupling_modes,        &
                                                  & energy_to_force_ratio  )
          reduced_sampling_matrix = mat(sample_matrix((l-1)*dims+1:l*dims,:i))
          new_fitting_matrix(:i,:i) = dble(         &
             &   transpose(reduced_sampling_matrix) &
             & * reduced_sampling_matrix            )
  
          qr = qr_decomposition( old_fitting_matrix(:i,:i) &
                             & + new_fitting_matrix(:i,:i) )
          
          evals = [(abs(qr%r(m,m)), m=1, size(qr%r,1))]
  
          if (any(evals<1e-200_dp)) then
            determinants(k) = 0.0_dp
          else
            determinants(k) = product(evals**(1.0_dp/size(qr%r,1)))
          endif
        endif
      enddo
      k = maxloc(abs(determinants), 1, mask=.not.basis_points_used)
      sampling_points(l) = basis_points(k)
      basis_points_used(k) = .true.
      if (i==size(basis_functions) .and. j==points_per_basis_function) then
        call print_line( 'Average |eigenvalue| of fitting matrix: '// &
                       & determinants(k)*energy_to_force_ratio//' (Ha).')
      endif
      
      sample_matrix((l-1)*dims+1:l*dims,:i) = construct_sample_matrix( &
                                               & basis_functions(:i),  &
                                               & basis_points(k:k),    &
                                               & coupling_modes,       &
                                               & energy_to_force_ratio )
      reduced_sampling_matrix = mat(sample_matrix((l-1)*dims+1:l*dims,:i))
      old_fitting_matrix(:i,:i) = old_fitting_matrix(:i,:i)                &
                              & + dble( transpose(reduced_sampling_matrix) &
                              &       * reduced_sampling_matrix            )
    enddo
  enddo
  
  output = SamplingPoints(sampling_points)
end procedure

module procedure generate_basis_points
  type(ComplexMonomial), allocatable :: monomials(:)
  
  integer :: monomial_power
  
  type(RealSingleDisplacement), allocatable :: displacement(:)
  
  type(RealMode) :: mode
  integer        :: sign_power
  real(dp)       :: sign
  real(dp)       :: magnitude
  
  integer :: i,j,ialloc
  
  monomials = basis_function%terms()
  
  allocate(output(2*size(monomials)), stat=ialloc); call err(ialloc)
  do i=1,size(monomials)
    ! Each monomial `prod_j (u_j)^(n_j)` is converted to a displacement
    !    sum_j s_j*sqrt(n_j/n_total)*u_j,
    ! where n_total = sum_j n_j
    ! s_j is a sign, which is
    !    -> +1 if (j==1 and n_total even) or (j/=1 and n_j even).
    !    -> -1 otherwise.
    ! The choice of sign is to help reduce duplication if the system has an
    !    inversion symmetry.
    monomial_power = monomials(i)%total_power()
    allocate(displacement(0), stat=ialloc); call err(ialloc)
    do j=1,size(monomials(i))
      mode = modes(first(modes%id==monomials(i)%id(j)))
      
      if (monomials(i)%power(j)>0) then
        if (j==1) then
          sign_power = monomial_power
        else
          sign_power = monomials(i)%power(j)
        endif
        
        if (modulo(sign_power,2)==1) then
          sign = -1.0_dp
        else
          sign = 1.0_dp
        endif
        
        magnitude = sign                                              &
                & * max_displacement%max_displacement(mode%frequency) &
                & * sqrt(monomials(i)%power(j)*1.0_dp/monomial_power)
        displacement = [ displacement,                              &
                       & RealSingleDisplacement(mode%id, magnitude) ]
      endif
      
      if (monomials(i)%id(j)/=monomials(i)%paired_id(j)) then
        if (monomials(i)%paired_power(j)>0) then
          if (j==1 .and. monomials(i)%power(j)==0) then
            sign_power = monomial_power
          else
            sign_power = monomials(i)%power(j)
          endif
          
          if (modulo(sign_power,2)==1) then
            sign = -1.0_dp
          else
            sign = 1.0_dp
          endif
          
          magnitude = sign                                              &
                  & * max_displacement%max_displacement(mode%frequency) &
                  & * sqrt(monomials(i)%paired_power(j)*1.0_dp/monomial_power)
          displacement = [ displacement,                                     &
                         & RealSingleDisplacement(mode%paired_id, magnitude) ]
        endif
      endif
    enddo
    
    ! The displacements are included with two lengths.
    !    -> n_total and n_total+1 if n_total is odd.
    !    -> n_total-1 and n_total if n_total is even.
    if (modulo(monomial_power,2)==1) then
      output(2*i-1) = RealModeDisplacement(displacement) &
                  & * (monomial_power*1.0_dp/potential_expansion_order)
      output(2*i) = RealModeDisplacement(displacement) &
                & * ((monomial_power+1)*1.0_dp/potential_expansion_order)
    else
      output(2*i-1) = RealModeDisplacement(displacement) &
                  & * ((monomial_power-1)*1.0_dp/potential_expansion_order)
      output(2*i) = RealModeDisplacement(displacement) &
                & * (monomial_power*1.0_dp/potential_expansion_order)
    endif
    deallocate(displacement, stat=ialloc); call err(ialloc)
  enddo
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
