submodule (caesar_real_complex_conversion_module) caesar_real_complex_conversion_submodule
  use caesar_normal_mode_module
contains

module procedure complex_to_real_mode
  type(RealVector), allocatable :: a(:)
  type(RealVector), allocatable :: b(:)
  
  type(RealVector), allocatable :: cos_vector(:)
  type(RealVector), allocatable :: sin_vector(:)
  
  integer :: qpoint_id_plus
  integer :: qpoint_id_minus
  
  if (input%id==input%paired_id) then
    ! output = input is its own pair. It is already real.
    ! c = u+, s = u- = 0.
    qpoint_id_plus = input%qpoint_id
    qpoint_id_minus = qpoint_id_plus
    a = real(input%unit_vector)
    b = a * 0.0_dp
    cos_vector = a
    sin_vector = b
  elseif (input%id<input%paired_id) then
    ! output is c.
    ! input is u+.
    ! u+ = (a+ib)e^{ 2pi i q.R}
    ! c  = sqrt(2)(a*cos(2pi q.R) - b*sin(2pi q.R))
    qpoint_id_plus = input%qpoint_id
    qpoint_id_minus = input%paired_qpoint_id
    a = real( input%unit_vector)
    b = aimag(input%unit_vector)
    cos_vector =  sqrt(2.0_dp) * a
    sin_vector = -sqrt(2.0_dp) * b
  elseif (input%id>input%paired_id) then
    ! output is s.
    ! input is u-.
    ! u- = (a-ib)e^{-2pi i q.R}
    ! s  = sqrt(2)(a*sin(2pi q.R) + b*cos(2pi q.R))
    qpoint_id_plus = input%paired_qpoint_id
    qpoint_id_minus = input%qpoint_id
    a =  real( input%unit_vector)
    b = -aimag(input%unit_vector)
    cos_vector = sqrt(2.0_dp) * b
    sin_vector = sqrt(2.0_dp) * a
  else
    call err()
  endif
  
  output = RealMode(                                  &
     & id                 = input%id,                 &
     & paired_id          = input%paired_id,          &
     & frequency          = input%frequency,          &
     & spring_constant    = input%spring_constant,    &
     & soft_mode          = input%soft_mode,          &
     & translational_mode = input%translational_mode, &
     & cos_vector         = cos_vector,               &
     & sin_vector         = sin_vector,               &
     & qpoint_id_plus     = qpoint_id_plus,           &
     & qpoint_id_minus    = qpoint_id_minus,          &
     & subspace_id        = input%subspace_id         )
end procedure

module procedure real_to_complex_mode
  type(RealVector), allocatable :: a(:)
  type(RealVector), allocatable :: b(:)
  
  type(ComplexVector), allocatable :: unit_vector(:)
  
  integer :: qpoint_id
  integer :: paired_qpoint_id
  
  ! Construct vectors and q-point ids.
  if (input%id==input%paired_id) then
    ! output = input is its own pair. It is already real.
    ! c = u+, s = u- = 0.
    qpoint_id = input%qpoint_id_plus
    paired_qpoint_id = qpoint_id
    a = input%cos_vector
    unit_vector = cmplxvec(a)
  elseif (input%id<input%paired_id) then
    ! output is u+.
    ! input is c.
    ! c  = sqrt(2)(a*cos(2pi q.R) - b*sin(2pi q.R))
    ! u+ = (a+ib)e^{2pi i q.R}
    qpoint_id = input%qpoint_id_plus
    paired_qpoint_id = input%qpoint_id_minus
    a =  input%cos_vector / sqrt(2.0_dp)
    b = -input%sin_vector / sqrt(2.0_dp)
    unit_vector = cmplxvec(a,b)
  elseif (input%paired_id<input%id) then
    ! output is u-.
    ! input is s.
    ! s  = sqrt(2)(a*sin(2pi q.R) + b*cos(2pi q.R))
    ! u- = (a-ib)e^{-2pi i q.R}
    qpoint_id = input%qpoint_id_minus
    paired_qpoint_id = input%qpoint_id_plus
    a = input%sin_vector / sqrt(2.0_dp)
    b = input%cos_vector / sqrt(2.0_dp)
    unit_vector = cmplxvec(a,-b)
  else
    call err()
  endif
  
  output = ComplexMode(                               &
     & id                 = input%id,                 &
     & paired_id          = input%paired_id,          &
     & frequency          = input%frequency,          &
     & spring_constant    = input%spring_constant,    &
     & soft_mode          = input%soft_mode,          &
     & translational_mode = input%translational_mode, &
     & unit_vector        = unit_vector,              &
     & qpoint_id          = qpoint_id,                &
     & paired_qpoint_id   = paired_qpoint_id,         &
     & subspace_id        = input%subspace_id         )
end procedure

module procedure basis_conversion_matrix_ComplexMonomials_RealMonomials
  complex(dp), allocatable :: matrix(:,:)
  
  integer :: i,j,ialloc
  
  allocate(matrix(size(this),size(that)), stat=ialloc); call err(ialloc)
  do i=1,size(that)
    do j=1,size(this)
      matrix(j,i) = element(this(j),that(i),include_coefficients)
    enddo
  enddo
  
  output = mat(matrix)
end procedure

module procedure basis_conversion_matrix_RealMonomials_ComplexMonomials
  complex(dp), allocatable :: matrix(:,:)
  
  integer :: i,j,ialloc
  
  allocate(matrix(size(this),size(that)), stat=ialloc); call err(ialloc)
  do i=1,size(that)
    do j=1,size(this)
      matrix(j,i) = element(this(j),that(i),include_coefficients)
    enddo
  enddo
  
  output = mat(matrix)
end procedure

module procedure coefficient_conversion_matrix_ComplexMonomials_RealMonomials
  output = transpose(basis_conversion_matrix(that,this,include_coefficients))
end procedure

module procedure coefficient_conversion_matrix_RealMonomials_ComplexMonomials
  output = transpose(basis_conversion_matrix(that,this,include_coefficients))
end procedure

module procedure element_ComplexMonomial_RealMonomial
  type(ComplexUnivariate) :: this_mode
  type(RealUnivariate)    :: that_mode
  
  logical, allocatable :: that_mode_accounted(:)
  
  complex(dp) :: overlap
  
  integer :: i,j,p,q,l,n
  
  if (this%total_power()/=that%total_power()) then
    output = 0
    return
  endif
  
  output = 1.0_dp
  
  that_mode_accounted = [(.false., i=1, size(that))]
  
  ! Loop over modes, multiplying the output by the overlap of each mode pair
  !    in 'this' with the eqivalent mode pair in 'that'.
  do i=1,size(this)
    this_mode = this%mode(i)
    
    ! Find the location of this%modes(i) in 'that'.
    j = first_equivalent(that%ids(), this%id(i), default=0, sorted=.true.)
    
    if (j==0) then
      ! If the modes doesn't exist in that, and the mode in 'this' has
      !    non-zero power, then the overlap between 'this' and 'that' is zero.
      ! If the mode in 'this' has zero power, then the mode has zero power
      !    in both, and can be neglected.
      if (this_mode%power/=0 .or. this_mode%paired_power/=0) then
        output = 0
        return
      endif
    else
      that_mode = that%mode(j)
      ! If the mode exists in both, account for it in the overlap.
      if (this_mode%total_power()/=that_mode%total_power()) then
        ! There is no overlap between this and that.
        output = 0
        return
      elseif (this_mode%id==this_mode%paired_id) then
        ! The mode is real. The overlap is 1.
        continue
      else
        ! this_mode = (u+)^p * (u-)^{n-p}.
        ! that_mode = (c )^q * (s )^{n-q}.
        !
        ! (u+)^p * (u-)^{n-p} = sum_{q=0}^n M(p,q,n) c^q * s*{n-q}.
        !
        ! (u+)^p * (u-)^{n-p} = (c+is)^p * (c-is)^(n-p) / sqrt(2)^n
        !
        ! => this = sum_{l=0}^p sum_{m=0}^{n-p}
        !           [ bin(p,l) bin(n-p,m) i^(-n+2j-l+m) c^{l+m} s^{n-l-m} ]
        !           / sqrt(2)^n
        !
        ! q = l+m
        !
        ! => sum_{l=0}^n sum_{m=0}^{n-p} = sum_{q=0}^n sum_{l=max(0,p+q-n)}^{min(p,q)}
        !
        ! => that = sum_{q,l} [ bin(p,l)*bin(n-p,q-l) i^(-n+2j+q-2l) c^q s^{n-q} ]
        !         / sqrt(2)^n
        !
        ! => M(p,q,n) = i^{-n+2j+q} / sqrt(2)^n
        !             * sum_{l=max(0,p+q-n)}^{min(p,q)} bin(p,l)*bin(n-p,q-l)*(-1)^l
        p = this_mode%power
        q = that_mode%power
        n = this_mode%total_power()
        overlap = 0
        do l=max(0,p+q-n),min(p,q)
          overlap = overlap                                      &
                & + exp(log_binomial(p,l)+log_binomial(n-p,q-l)) &
                & * (-1)**l
        enddo
        overlap = overlap                            &
             & * cmplx(0.0_dp,1.0_dp,dp)**(-n+2*p+q) &
             & / sqrt(2.0_dp)**n
        output = output * overlap
      endif
      
      that_mode_accounted(j) = .true.
    endif
  enddo
  
  ! Check remaining modes in 'that'. If any have non-zero power then the
  !    overlap is zero.
  do i=1,size(that)
    if (.not. that_mode_accounted(j)) then
      if (that%power(j)/=0 .or. that%paired_power(j)/=0) then
        output = 0
        return
      endif
    endif
  enddo
  
  if (include_coefficients) then
    output = output * this%coefficient / that%coefficient
  endif
end procedure

module procedure element_RealMonomial_ComplexMonomial
  type(RealUnivariate)    :: this_mode
  type(ComplexUnivariate) :: that_mode
  
  logical, allocatable :: that_mode_accounted(:)
  
  complex(dp) :: overlap
  
  integer :: i,j,p,q,l,n
  
  if (this%total_power()/=that%total_power()) then
    output = 0
    return
  endif
  
  output = 1.0_dp
  
  that_mode_accounted = [(.false., i=1, size(that))]
  
  ! Loop over modes, multiplying the output by the overlap of each mode pair
  !    in 'this' with the eqivalent mode pair in 'that'.
  do i=1,size(this)
    this_mode = this%mode(i)
    
    ! Find the location of this%modes(i) in 'that'.
    j = first_equivalent(that%ids(), this%id(i), default=0, sorted=.true.)
    
    if (j==0) then
      ! If the modes doesn't exist in that, and the mode in 'this' has
      !    non-zero power, then the overlap between 'this' and 'that' is zero.
      ! If the mode in 'this' has zero power, then the mode has zero power
      !    in both, and can be neglected.
      if (this%power(i)/=0 .or. this%paired_power(i)/=0) then
        output = 0
        return
      endif
    else
      that_mode = that%mode(j)
      
      ! If the mode exists in both, account for it in the overlap.
      if (this_mode%total_power()/=that_mode%total_power()) then
        ! There is no overlap between this and that.
        output = 0
        return
      elseif (this_mode%id==this_mode%paired_id) then
        ! The mode is real. The overlap is 1.
        continue
      else
        ! this_mode = (c )^p * (s )^{n-p}.
        ! that_mode = (u+)^q * (u-)^{n-q}.
        !
        ! c^p * s^{n-p} = sum_{q=0}^n M(p,q,n) (u+)^q * (u-)^{n-q}.
        !
        ! c^p * s^{n-p} = (u+ + u-)^p * ((u+ - u-)/i)^{n-p} / sqrt(2)^n
        !
        !    = sum_{l=0}^p sum_{m=0}^{n-p}
        !      [ bin(p,l) bin(n-p,m) (u+)^{l+m} (u-)^{n-l-m} (-1)^{n-p-m} ]
        !    / (sqrt(2)^n * i^{n-p})
        !
        ! q = l+m
        !
        ! => sum_{l=0}^n sum_{m=0}^{n-p} = sum_{q=0}^n sum_{l=max(0,p+q-n)}^{min(p,q)}
        !
        ! => that = sum_{q,l}[ bin(p,l) bin(n-p,q-l) u+^q u-^{n-q} (-1)^{n-p-q+l} ]
        !      / (sqrt(2)^n * i^{n-p})
        !
        ! => M(p,q,n) = i^{n-p-2k} / sqrt(2)^n
        !             * sum_{l=max(0,p+q-n)}^{min(p,q)} bin(p,l)*bin(n-p,q-l)*(-1)^l
        p = this_mode%power
        q = that_mode%power
        n = this_mode%total_power()
        overlap = 0
        do l=max(0,p+q-n),min(p,q)
          overlap = overlap                                      &
                & + exp(log_binomial(p,l)+log_binomial(n-p,q-l)) &
                & * (-1)**l
        enddo
        overlap = overlap                             &
             & * cmplx(0.0_dp,1.0_dp,dp)**(n-p-2*q) &
             & / sqrt(2.0_dp)**n
        output = output * overlap
      endif
      
      that_mode_accounted(j) = .true.
    endif
  enddo
  
  ! Check remaining modes in 'that'. If any have non-zero power then the
  !    overlap is zero.
  do i=1,size(that)
    if (.not. that_mode_accounted(j)) then
      if (that%power(j)/=0 .or. that%paired_power(j)/=0) then
        output = 0
        return
      endif
    endif
  enddo
  
  if (include_coefficients) then
    output = output * this%coefficient / that%coefficient
  endif
end procedure

module procedure new_PairedMonomials_ComplexMonomials
  type(ComplexMonomial) :: conjugate
  
  integer, allocatable :: conjugates(:)
  
  integer :: i,j,ialloc
  
  conjugates = [(0, i=1, size(input))]
  do i=1,size(input)
    if (conjugates(i)/=0) then
      cycle
    endif
    
    conjugate = conjg(input(i))
    
    j = first_equivalent( input,                    &
                        & conjugate,                &
                        & compare_complex_monomials )
    
    conjugates(i) = j
    conjugates(j) = i
  enddo
  
  allocate(output(size(input)), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    if (conjugates(i)==i) then
      output(i) = PairedMonomial( real(input(i)%coefficient), &
                                & input(i)%modes(),           &
                                & pair=0                      )
    elseif (conjugates(i)>i) then
      output(i) = PairedMonomial( real(input(i)%coefficient), &
                                & input(i)%modes(),           &
                                & pair=1                      )
      output(conjugates(i)) = PairedMonomial( aimag(input(i)%coefficient), &
                                            & input(i)%modes(),            &
                                            & pair=-1                      )
    else
      continue
    endif
  enddo
end procedure

module procedure new_ComplexMonomials_PairedMonomials
  integer, allocatable :: matches(:)
  
  complex(dp) :: coefficient
  
  integer :: i,j,ialloc
  
  matches = [(0, i=1, size(input))]
  do i=1,size(input)
    if (matches(i)/=0) then
      cycle
    endif
    
    j = first_equivalent( input,        &
                        & input(i),     &
                        & matching_pair )
    
    matches(i) = j
    matches(j) = i
  enddo
  
  allocate(output(size(input)), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    if (matches(i)==i) then
      coefficient = cmplx(input(i)%coefficient, 0.0_dp, dp)
      output(i) = ComplexMonomial( coefficient,     &
                                 & input(i)%modes() )
    elseif (matches(i)>i) then
      coefficient = cmplx( input(i)%coefficient,          &
                         & input(matches(i))%coefficient, &
                         & dp                             )
      output(i) = ComplexMonomial( coefficient,     &
                                 & input(i)%modes() )
      output(matches(i)) = ComplexMonomial( conjg(coefficient),    &
                                          & conjg(input(i)%modes()) )
    else
      continue
    endif
  enddo
end procedure

module procedure new_PairedPolynomial_ComplexPolynomial
  output = PairedPolynomial(PairedMonomial(input%terms))
end procedure

module procedure new_ComplexPolynomial_PairedPolynomial
  output = ComplexPolynomial(ComplexMonomial(input%terms))
end procedure
end submodule
