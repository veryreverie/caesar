submodule (caesar_complex_polynomial_module) caesar_complex_polynomial_submodule
  use caesar_normal_mode_module
contains

module procedure new_ComplexUnivariate
  if (id<=paired_id) then
    if (id==paired_id .and. power/=paired_power) then
      call print_line(CODE_ERROR//': modes the same, but powers differ.')
      call err()
    endif
    this%id           = id
    this%paired_id    = paired_id
    this%power        = power
    this%paired_power = paired_power
  else
    this%id           = paired_id
    this%paired_id    = id
    this%power        = paired_power
    this%paired_power = power
  endif
end procedure

module procedure new_ComplexUnivariate_ComplexMode
  if (present(paired_power)) then
    if (mode%id==mode%paired_id .and. power/=paired_power) then
      call print_line(ERROR//': Mode is its own pair, but power does not &
         & match paired_power.')
      call err()
    endif
    this = ComplexUnivariate( id           = mode%id,        &
                            & paired_id    = mode%paired_id, &
                            & power        = power,          &
                            & paired_power = paired_power    )
  else
    if (mode%id==mode%paired_id) then
      this = ComplexUnivariate( id           = mode%id, &
                              & paired_id    = mode%id, &
                              & power        = power,   &
                              & paired_power = power    )
    else
      this = ComplexUnivariate( id           = mode%id,        &
                              & paired_id    = mode%paired_id, &
                              & power        = power,          &
                              & paired_power = 0               )
    endif
  endif
end procedure

module procedure new_ComplexMonomial
  this%coefficient = coefficient
  this%modes_      = modes(sort(modes%id))
end procedure

module procedure new_ComplexMonomial_ComplexMonomialable
  this = input%to_ComplexMonomial()
end procedure

module procedure new_ComplexPolynomial
  this%terms = terms
end procedure

module procedure new_ComplexPolynomial_ComplexPolynomialable
  this = input%to_ComplexPolynomial()
end procedure

module procedure to_ComplexMonomial_ComplexUnivariate
  select type(this); type is(ComplexUnivariate)
    output = ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp), &
                            & modes       = [this])
  class default
    call err()
  end select
end procedure

module procedure to_ComplexPolynomial_ComplexUnivariate
  select type(this); type is(ComplexUnivariate)
    output = ComplexPolynomial([this%to_ComplexMonomial()])
  class default
    call err()
  end select
end procedure

module procedure to_ComplexMonomial_ComplexMonomial
  select type(this); type is(ComplexMonomial)
    output = this
  class default
    call err()
  end select
end procedure

module procedure to_ComplexPolynomial_ComplexMonomial
  select type(this); type is(ComplexMonomial)
    output = ComplexPolynomial([this])
  class default
    call err()
  end select
end procedure

module procedure to_ComplexPolynomial_ComplexPolynomial
  select type(this); type is(ComplexPolynomial)
    output = this
  class default
    call err()
  end select
end procedure

module procedure size_ComplexMonomial
  output = size(this%modes_)
end procedure

module procedure size_ComplexPolynomial
  output = size(this%terms)
end procedure

module procedure mode_ComplexMonomial
  output = this%modes_(index)
end procedure

module procedure modes_ComplexMonomial
  integer, allocatable :: sort_key(:)
  
  integer :: i,j,k,ialloc
  
  if (count([present(indices),present(ids),present(exclude_ids)])>1) then
    call print_line(CODE_ERROR//': too many optional arguments present.')
    call err()
  elseif (present(ids) .neqv. present(paired_ids)) then
    call print_line(CODE_ERROR// &
       & ': Only one of ids and paired_ids is present.')
    call err()
  endif
  
  if (present(indices)) then
    ! Return the modes with the specified indices.
    output = this%modes_(indices)
  elseif (present(ids)) then
    ! Return only modes with given ids, in the order of the ids specified.
    ! If an id is not present in this monomial, instead include u^0.
    sort_key = sort(ids)
    allocate(output(size(ids)), stat=ialloc); call err(ialloc)
    j = 1
    ! Loop over the ids in ascending order (i.e. ids(sort_key)).
    do i=1,size(ids)
      ! Cycle through this%modes_ until the mode with id=ids(sort_key(i))
      !    is found.
      do
        if (j>size(this%modes_)) then
          exit
        elseif (this%modes_(j)%id>=ids(sort_key(i))) then
          exit
        else
          j = j+1
        endif
      enddo
      
      ! If such a mode exists, add it to the output.
      if (j<=size(this%modes_)) then
        if (this%modes_(j)%id == ids(sort_key(i))) then
          output(sort_key(i)) = this%modes_(j)
          cycle
        endif
      endif
      
      ! If no such mode exists, add u^0 to the output.
      output(sort_key(i)) = ComplexUnivariate(     &
         & id           = ids(sort_key(i)),        &
         & paired_id    = paired_ids(sort_key(i)), &
         & power        = 0,                       &
         & paired_power = 0                        )
    enddo
  elseif (present(exclude_ids)) then
    ! Return all modes apart from those in exclude_ids.
    sort_key = sort(exclude_ids)
    
    allocate(output(size(this%modes_)), stat=ialloc); call err(ialloc)
    j = 1
    k = 0
    ! Loop over this%modes_.
    do i=1,size(this%modes_)
      ! Cycle through exclude_ids in ascending order
      !    (i.e. exclude_ids(sort_key)) until
      !    exculde_ids(j)=this%modes_(i)%id.
      do
        if (j>size(exclude_ids)) then
          exit
        elseif (exclude_ids(sort_key(j))>=this%modes_(i)%id) then
          exit
        else
          j = j+1
        endif
      enddo
      
      ! If the mode is in exclude_ids, exclude it.
      if (j<=size(exclude_ids)) then
        if (exclude_ids(sort_key(j))==this%modes_(i)%id) then
          cycle
        endif
      endif
      
      ! If the mode is not in exclude_ids, add it to the output.
      k = k+1
      output(k) = this%modes_(i)
    enddo
    output = output(:k)
  else
    ! Return all modes in the monomial.
    output = this%modes_
  endif
end procedure

module procedure id_ComplexMonomial
  output = this%modes_(index)%id
end procedure

module procedure ids_ComplexMonomial
  if (present(indices)) then
    output = this%modes_(indices)%id
  else
    output = this%modes_%id
  endif
end procedure

module procedure paired_id_ComplexMonomial
  output = this%modes_(index)%paired_id
end procedure

module procedure paired_ids_ComplexMonomial
  if (present(indices)) then
    output = this%modes_(indices)%paired_id
  else
    output = this%modes_%paired_id
  endif
end procedure

module procedure power_ComplexMonomial
  output = this%modes_(index)%power
end procedure

module procedure powers_ComplexMonomial
  if (present(indices)) then
    output = this%modes_(indices)%power
  else
    output = this%modes_%power
  endif
end procedure

module procedure paired_power_ComplexMonomial
  output = this%modes_(index)%paired_power
end procedure

module procedure paired_powers_ComplexMonomial
  if (present(indices)) then
    output = this%modes_(indices)%paired_power
  else
    output = this%modes_%paired_power
  endif
end procedure

module procedure set_modes_ComplexMonomial
  this%modes_ = modes
end procedure

module procedure simplify_ComplexMonomial
  integer :: i
  
  ! Combine modes with the same ID, remove modes with power=paired_power=0.
  i = 1
  do while(i<=size(this))
    if (this%modes_(i)%power<0 .or. this%modes_(i)%paired_power<0) then
      call err()
    elseif (this%modes_(i)%power==0 .and. this%modes_(i)%paired_power==0) then
      this%modes_ = [this%modes_(:i-1), this%modes_(i+1:)]
      cycle
    endif
    
    if (i>1) then
      if (this%modes_(i)%id==this%modes_(i-1)%id) then
        this%modes_(i-1)%power = this%modes_(i-1)%power + this%modes_(i)%power
        this%modes_(i-1)%paired_power = this%modes_(i-1)%paired_power &
                                    & + this%modes_(i)%paired_power
        this%modes_ = [this%modes_(:i-1), this%modes_(i+1:)]
        cycle
      endif
    endif
    
    i = i+1
  enddo
end procedure

module procedure simplify_ComplexPolynomial
  integer,     allocatable :: unique_terms(:)
  complex(dp), allocatable :: coefficients(:)
  
  integer :: i,j
  
  call this%terms%simplify()
  
  ! Add together monomials with the same powers.
  unique_terms = set(this%terms, compare_complex_monomials)
  coefficients = [(cmplx(0.0_dp,0.0_dp,dp), i=1, size(unique_terms))]
  do i=1,size(this%terms)
    do j=1,size(unique_terms)
      if (compare_complex_monomials( this%terms(i),              &
                                   & this%terms(unique_terms(j)) )) then
        coefficients(j) = coefficients(j) + this%terms(i)%coefficient
      endif
    enddo
  enddo
  
  this%terms = this%terms(unique_terms)
  this%terms%coefficient = coefficients
end procedure

module procedure conjg_ComplexUnivariate
  output = ComplexUnivariate( id           = this%paired_id,   &
                            & paired_id    = this%id,          &
                            & power        = this%power,       &
                            & paired_power = this%paired_power )
end procedure

module procedure conjg_ComplexMonomial
  output = ComplexMonomial( coefficient = conjg(this%coefficient), &
                          & modes       = conjg(this%modes_)       )
end procedure

module procedure conjg_ComplexPolynomial
  output = ComplexPolynomial(conjg(this%terms))
end procedure

module procedure total_power_ComplexUnivariate
  if (this%id==this%paired_id) then
    output = this%power
  else
    output = this%power + this%paired_power
  endif
end procedure

module procedure total_power_ComplexMonomial
  output = sum(this%modes_%total_power())
end procedure

module procedure wavevector_ComplexUnivariate
  type(ComplexMode) :: mode
  type(QpointData)  :: qpoint
  
  mode = modes(first(modes%id==this%id))
  qpoint = qpoints(first(qpoints%id==mode%qpoint_id))
  if (this%id==this%paired_id) then
    output = qpoint%qpoint * this%power
  else
    output = qpoint%qpoint * (this%power-this%paired_power)
  endif
  
  output = vec(modulo(frac(output),1))
end procedure

module procedure wavevector_ComplexMonomial
  integer :: i
  
  output = fracvec(zeroes(3))
  do i=1,size(this%modes_)
    output = output + this%modes_(i)%wavevector(modes,qpoints)
  enddo
  
  output = vec(modulo(frac(output),1))
end procedure

module procedure energy_RealSingleDisplacement_ComplexUnivariate
  complex(dp) :: magnitude
  complex(dp) :: paired_magnitude
  
  if (.not. (present(displacement) .or. present(paired_displacement))) then
    call print_line(CODE_ERROR//': Neither displacement passed.')
    call err()
  elseif (present(paired_displacement) .and. this%id==this%paired_id) then
    call print_line(CODE_ERROR//': paired displacement passed to &
       &univariate with id=paired_id.')
    call err()
  endif
  
  if (this%id==this%paired_id) then
    magnitude = displacement%magnitude
  else
    ! Convert from real to complex co-ordinates.
    ! x_+ = (x_c + ix_s) / sqrt(2)
    ! x_- = (x_c - ix_s) / sqrt(2)
    magnitude = 0
    paired_magnitude = 0
    if (present(displacement)) then
      ! x_+ : x_c / sqrt(2)
      magnitude = magnitude &
              & + displacement%magnitude &
              & / sqrt(2.0_dp)
      ! x_- : x_c / sqrt(2)
      paired_magnitude = paired_magnitude &
                     & + displacement%magnitude &
                     & / sqrt(2.0_dp)
    endif
    if (present(paired_displacement)) then
      ! x_+ : ix_s / sqrt(2)
      magnitude = magnitude &
              & + paired_displacement%magnitude &
              & * cmplx(0,1,dp) / sqrt(2.0_dp)
      ! x_- : - ix_s / sqrt(2)
      paired_magnitude = paired_magnitude &
                     & - paired_displacement%magnitude &
                     & * cmplx(0,1,dp) / sqrt(2.0_dp)
    endif
  endif
  
  output = magnitude**this%power
  if (this%id/=this%paired_id) then
    output = output * paired_magnitude**this%paired_power
  endif
end procedure

module procedure energy_ComplexSingleDisplacement_ComplexUnivariate
  if (.not. (present(displacement) .or. present(paired_displacement))) then
    call print_line(CODE_ERROR//': Neither displacement passed.')
    call err()
  elseif (present(paired_displacement) .and. this%id==this%paired_id) then
    call print_line(CODE_ERROR//': paired displacement passed to &
       &univariate with id=paired_id.')
    call err()
  endif
  
  output = 1.0_dp
  
  if (present(displacement)) then
    output = output * displacement%magnitude**this%power
  endif
  
  if (present(paired_displacement)) then
    output = output * paired_displacement%magnitude**this%paired_power
  endif
end procedure

module procedure energy_RealModeDisplacement_ComplexMonomial
  integer              :: max_mode_id
  integer, allocatable :: indices(:)
  
  integer :: i,j,k,ialloc
  
  ! Generate `indices`, such that
  !    indices(i) = first(displacement%vectors%id==i, default=0)
  !    for all i in this%modes_%id and this%modes_%paired_id.
  ! N.B. this is only for optimisation purposes.
  ! As such, indices(i) is only initialised if `i` is in
  !    displacement%vectors%id, this%modes_%id or this%modes_%paired_id.
  ! This reduces a double loop across this%modes_ and displacement%vectors
  !    to two single loops.
  max_mode_id = maxval([ maxval(displacement%vectors%id), &
                       & maxval(this%modes_%id),          &
                       & maxval(this%modes_%paired_id)    ])
  allocate(indices(max_mode_id), stat=ialloc); call err(ialloc)
  do i=1,size(this%modes_)
    indices(this%modes_(i)%id) = 0
    indices(this%modes_(i)%paired_id) = 0
  enddo
  do i=1,size(displacement%vectors)
    indices(displacement%vectors(i)%id) = i
  enddo
  
  output = this%coefficient
  
  ! Loop over the modes in this monomial,
  !    find the mode (and paired mode if relevant) in the displacement
  !    corresponding to each mode,
  !    and calculate the energy term by term.
  do i=1,size(this)
    associate (mode => this%modes_(i))
      if (mode%id==mode%paired_id) then
        if (mode%power/=0) then
          j = indices(mode%id)
          
          if (j==0) then
            ! If the mode is not present in the displacement,
            !    then the displacement along that mode is zero.
            ! As such, the monomial is zero. (0^n=0 if n>0).
            output = 0
            return
          else
            output = output &
                 & * mode%energy_real(displacement=displacement%vectors(j))
          endif
        endif
      else
        if (mode%power/=0 .or. mode%paired_power/=0) then
          j = indices(mode%id)
          k = indices(mode%paired_id)
          
          if (j==0 .and. k==0) then
            ! If both modes are not present in the displacement,
            !    then the displacement along that mode is zero.
            ! As such, the monomial is zero. (0^n=0 if n>0).
            output = 0
            return
          elseif (j==0) then
            output = output * mode%energy_real(              &
               & paired_displacement=displacement%vectors(k) )
          elseif (k==0) then
            output = output * mode%energy_real(       &
               & displacement=displacement%vectors(j) )
          else
            output = output * mode%energy_real(              &
               & displacement=displacement%vectors(j),       &
               & paired_displacement=displacement%vectors(k) )
          endif
        endif
      endif
    end associate
  enddo
end procedure

module procedure energy_ComplexModeDisplacement_ComplexMonomial
  integer              :: max_mode_id
  integer, allocatable :: indices(:)
  
  integer :: i,j,k,ialloc
  
  ! Generate `indices`, such that
  !    indices(i) = first(displacement%vectors%id==i, default=0)
  !    for all i in this%modes_%id and this%modes_%paired_id.
  ! N.B. this is only for optimisation purposes.
  ! As such, indices(i) is only initialised if `i` is in
  !    displacement%vectors%id, this%modes_%id or this%modes_%paired_id.
  ! This reduces a double loop across this%modes_ and displacement%vectors
  !    to two single loops.
  max_mode_id = maxval([ maxval(displacement%vectors%id), &
                       & maxval(this%modes_%id),          &
                       & maxval(this%modes_%paired_id)    ])
  allocate(indices(max_mode_id), stat=ialloc); call err(ialloc)
  do i=1,size(this%modes_)
    indices(this%modes_(i)%id) = 0
    indices(this%modes_(i)%paired_id) = 0
  enddo
  do i=1,size(displacement%vectors)
    indices(displacement%vectors(i)%id) = i
  enddo
  
  output = this%coefficient
  
  ! Loop over the modes in this monomial,
  !    find the mode (and paired mode if relevant) in the displacement
  !    corresponding to each mode,
  !    and calculate the energy term by term.
  do i=1,size(this)
    associate (mode => this%modes_(i))
      if (mode%id==mode%paired_id) then
        if (mode%power/=0) then
          j = indices(mode%id)
          
          if (j==0) then
            ! If the mode is not present in the displacement,
            !    then the displacement along that mode is zero.
            ! As such, the monomial is zero. (0^n=0 if n>0).
            output = 0
            return
          else
            output = output * mode%energy_complex(    &
               & displacement=displacement%vectors(j) )
          endif
        endif
      else
        if (mode%power/=0 .or. mode%paired_power/=0) then
          j = indices(mode%id)
          k = indices(mode%paired_id)
          
          if (      (j==0 .and. mode%power/=0)        &
             & .or. (k==0 .and. mode%paired_power/=0) ) then
            ! If either mode is not present in the displacement,
            !    then the displacement along that mode is zero.
            ! As such, the monomial is zero. (0^n=0 if n>0).
            output = 0
            return
          elseif (j==0) then
            output = output * mode%energy_complex(           &
               & paired_displacement=displacement%vectors(k) )
          elseif (k==0) then
            output = output * mode%energy_complex(    &
               & displacement=displacement%vectors(j) )
          else
            output = output * mode%energy_complex(           &
               & displacement=displacement%vectors(j),       &
               & paired_displacement=displacement%vectors(k) )
          endif
        endif
      endif
    end associate
  enddo
end procedure

module procedure energy_RealModeDisplacement_ComplexPolynomial
  output = sum(this%terms%energy(displacement))
end procedure

module procedure energy_ComplexModeDisplacement_ComplexPolynomial
  output = sum(this%terms%energy(displacement))
end procedure

module procedure force_RealSingleDisplacement_ComplexUnivariate
  complex(dp) :: magnitude
  complex(dp) :: paired_magnitude
  
  complex(dp) :: values(2)
  complex(dp) :: derivatives(2)
  
  if (.not. (present(displacement) .or. present(paired_displacement))) then
    call print_line(CODE_ERROR//': Neither displacement passed.')
    call err()
  elseif (present(paired_displacement) .and. this%id==this%paired_id) then
    call print_line(CODE_ERROR//': paired displacement passed to &
       &univariate with id=paired_id.')
    call err()
  endif
  
  if (this%id==this%paired_id) then
    magnitude = displacement%magnitude
  else
    ! Convert from real to complex co-ordinates.
    ! x_+ = (x_c + ix_s) / sqrt(2)
    ! x_- = (x_c - ix_s) / sqrt(2)
    magnitude = 0
    paired_magnitude = 0
    if (present(displacement)) then
      ! x_+ : x_c / sqrt(2)
      magnitude = magnitude &
              & + displacement%magnitude &
              & / sqrt(2.0_dp)
      ! x_- : x_c / sqrt(2)
      paired_magnitude = paired_magnitude &
                     & + displacement%magnitude &
                     & / sqrt(2.0_dp)
    endif
    if (present(paired_displacement)) then
      ! x_+ : ix_s / sqrt(2)
      magnitude = magnitude &
              & + paired_displacement%magnitude &
              & * cmplx(0,1,dp) / sqrt(2.0_dp)
      ! x_- : - ix_s / sqrt(2)
      paired_magnitude = paired_magnitude &
                     & - paired_displacement%magnitude &
                     & * cmplx(0,1,dp) / sqrt(2.0_dp)
    endif
  endif
  
  if (this%power>1) then
    values(1) = magnitude**this%power
    derivatives(1) = this%power * magnitude**(this%power-1)
  elseif (this%power==1) then
    values(1) = magnitude
    derivatives(1) = 1
  else
    values(1) = 1
    derivatives(1) = 0
  endif
  
  if (this%id==this%paired_id) then
    output = [-derivatives(1)]
    return
  endif
  
  if (this%paired_power>1) then
    values(2) = paired_magnitude**this%paired_power
    derivatives(2) = this%paired_power &
                 & * paired_magnitude**(this%paired_power-1)
  elseif (this%paired_power==1) then
    values(2) = paired_magnitude
    derivatives(2) = 1
  else
    values(2) = 1
    derivatives(2) = 0
  endif
  
  ! Construct the output in complex co-ordinates...
  output = [-derivatives(1)*values(2), -derivatives(2)*values(1)]
  ! ... then convert back to real co-ordinates.
  ! f_c =   (f_+ + f_-) / sqrt(2)
  ! f_s = - (f_+ - f_-) / sqrt(2)i
  output = [ ( output(1)+output(2))/sqrt(2.0_dp),            &
           & (-output(1)+output(2))/cmplx(0,sqrt(2.0_dp),dp) ]
end procedure

module procedure force_ComplexSingleDisplacement_ComplexUnivariate
  complex(dp) :: magnitude
  complex(dp) :: paired_magnitude
  
  complex(dp) :: values(2)
  complex(dp) :: derivatives(2)
  
  if (.not. (present(displacement) .or. present(paired_displacement))) then
    call print_line(CODE_ERROR//': Neither displacement passed.')
    call err()
  elseif (present(paired_displacement) .and. this%id==this%paired_id) then
    call print_line(CODE_ERROR//': paired displacement passed to &
       &univariate with id=paired_id.')
    call err()
  endif
  
  if (present(displacement)) then
    magnitude = displacement%magnitude
  else
    magnitude = 0
  endif
  
  if (this%power>1) then
    values(1) = magnitude**this%power
    derivatives(1) = this%power * magnitude**(this%power-1)
  elseif (this%power==1) then
    values(1) = magnitude
    derivatives(1) = 1
  else
    values(1) = 1
    derivatives(1) = 0
  endif
  
  if (this%id==this%paired_id) then
    output = [-derivatives(1)]
    return
  endif
  
  if (present(paired_displacement)) then
    paired_magnitude = paired_displacement%magnitude
  else
    paired_magnitude = 0
  endif
  
  if (this%paired_power>1) then
    values(2) = paired_magnitude**this%paired_power
    derivatives(2) = this%paired_power &
                 & * paired_magnitude**(this%paired_power-1)
  elseif (this%paired_power==1) then
    values(2) = paired_magnitude
    derivatives(2) = 1
  else
    values(2) = 1
    derivatives(2) = 0
  endif
  
  output = [-derivatives(1)*values(2), -derivatives(2)*values(1)]
end procedure

module procedure force_RealModeDisplacement_ComplexMonomial
  integer              :: max_mode_id
  integer, allocatable :: indices(:)
  
  integer :: value_is_zero
  logical, allocatable :: force_is_zero(:)
  
  type(RealSingleDisplacement), allocatable :: displacements(:)
  type(RealSingleDisplacement), allocatable :: paired_displacements(:)
  
  complex(dp)              :: energy
  integer,     allocatable :: ids(:)
  complex(dp), allocatable :: forces(:)
  
  integer :: i,j,k,ialloc
  
  integer :: no_modes
  
  ! Generate `indices`, such that
  !    indices(i) = first(displacement%vectors%id==i, default=0)
  !    for all i in this%modes_%id and this%modes_%paired_id.
  ! N.B. this is only for optimisation purposes.
  ! As such, indices(i) is only initialised if `i` is in
  !    displacement%vectors%id, this%modes_%id or this%modes_%paired_id.
  ! This reduces a double loop across this%modes_ and displacement%vectors
  !    to two single loops.
  max_mode_id = maxval([ maxval(displacement%vectors%id), &
                       & maxval(this%modes_%id),          &
                       & maxval(this%modes_%paired_id)    ])
  allocate(indices(max_mode_id), stat=ialloc); call err(ialloc)
  do i=1,size(this%modes_)
    indices(this%modes_(i)%id) = 0
    indices(this%modes_(i)%paired_id) = 0
  enddo
  do i=1,size(displacement%vectors)
    indices(displacement%vectors(i)%id) = i
  enddo
  
  ! Match each mode in the monomial with the matching mode(s)
  !    in the displacement.
  allocate( displacements(size(this)),        &
          & paired_displacements(size(this)), &
          & force_is_zero(size(this)),        &
          & stat=ialloc); call err(ialloc)
  value_is_zero = 0
  no_modes = 0
  do i=1,size(this)
    associate (mode => this%modes_(i))
      j = indices(mode%id)
      if (j/=0) then
        displacements(i) = displacement%vectors(j)
      else
        displacements(i) = RealSingleDisplacement(mode%id, 0.0_dp)
      endif
      if (mode%id==mode%paired_id) then
        if (j==0 .and. mode%power/=0) then
          if (value_is_zero/=0) then
            ! If U_i is zero for more than one i, then
            !    prod_{j/=i}[ {U_j}^{n_j} ] is always zero,
            !    so all derivatives are zero.
            output = RealModeForce([RealSingleForce::])
            return
          endif
          value_is_zero = i
        endif
        force_is_zero(i) = j==0 .and. mode%power/=1
        if (.not. force_is_zero(i)) then
          no_modes = no_modes+1
        endif
      else
        k = indices(mode%paired_id)
        if (k/=0) then
          paired_displacements(i) = displacement%vectors(k)
        else
          paired_displacements(i) = RealSingleDisplacement( mode%paired_id, &
                                                          & 0.0_dp          )
        endif
        if (       (j==0 .and. k==0)                         &
           & .and. (mode%power/=0 .or. mode%paired_power/=0) ) then
          if (value_is_zero/=0) then
            ! If U_i is zero for more than one i, then
            !    prod_{j/=i}[ {U_j}^{n_j} ] is always zero,
            !    so all derivatives are zero.
            output = RealModeForce([RealSingleForce::])
            return
          endif
          value_is_zero = i
        endif
        force_is_zero(i) = (j==0 .and. mode%power/=1) &
                   & .and. (k==0 .and. mode%paired_power/=1)
        if (.not. force_is_zero(i)) then
          no_modes = no_modes+2
        endif
      endif
      if (force_is_zero(i) .and. value_is_zero==i) then
        ! U_i and the derivative are zero along mode i,
        !    so all derivatives are zero.
        output = RealModeForce([RealSingleForce::])
        return
      endif
    end associate
  enddo
  
  i = value_is_zero
  if (i/=0) then
    ! Only one U_i is zero.
    associate (mode => this%modes_(i))
    if (mode%id==mode%paired_id) then
      ids = [mode%id]
      forces = mode%force_real(displacements(i))
    else
      ids = [mode%id, mode%paired_id]
      forces = mode%force_real(displacements(i), paired_displacements(i))
    endif
    end associate
    do j=1,size(this)
      if (j/=i) then
        associate (mode => this%modes_(j))
          if (mode%id==mode%paired_id) then
            forces = forces &
                 & * mode%energy_real(displacements(j))
          else
            forces = forces                                    &
                 & * mode%energy_real( displacements(j),       &
                 &                     paired_displacements(j) )
          endif
        end associate
      endif
    enddo
  else
    ! The force is non-zero along multiple modes.
    allocate( ids(no_modes),    &
            & forces(no_modes), &
            & stat=ialloc); call err(ialloc)
    forces = 1.0_dp
    j = 0
    do i=1,size(this)
      associate (mode => this%modes_(i))
        if (mode%id==mode%paired_id) then
          energy = mode%energy_real(displacements(i))
        else
          energy = mode%energy_real( displacements(i),       &
                                   & paired_displacements(i) )
        endif
        
        if (force_is_zero(i)) then
          forces = forces * energy
        else
          if (mode%id==mode%paired_id) then
            j = j+1
            ids(j) = mode%id
            forces(j:j) = forces(j:j) * mode%force_real(displacements(i))
            forces(:j-1) = forces(:j-1) * energy
            forces(j+1:) = forces(j+1:) * energy
          else
            j = j+2
            ids(j-1) = mode%id
            ids(j) = mode%paired_id
            forces(j-1:j) = mode%force_real( displacements(i),       &
                                           & paired_displacements(i) )
            forces(:j-2) = forces(:j-2) * energy
            forces(j+1:) = forces(j+1:) * energy
          endif
        endif
      end associate
    enddo
  endif
  
  ! Construct output from components along each mode.
  output = RealModeForce(RealSingleForce(ids, real(forces*this%coefficient)))
end procedure

module procedure force_ComplexModeDisplacement_ComplexMonomial
  integer              :: max_mode_id
  integer, allocatable :: indices(:)
  
  integer :: value_is_zero
  logical, allocatable :: force_is_zero(:)
  
  type(ComplexSingleDisplacement), allocatable :: displacements(:)
  type(ComplexSingleDisplacement), allocatable :: paired_displacements(:)
  
  complex(dp)              :: energy
  integer,     allocatable :: ids(:)
  complex(dp), allocatable :: forces(:)
  
  integer :: i,j,k,ialloc
  
  integer :: no_modes
  
  ! Generate `indices`, such that
  !    indices(i) = first(displacement%vectors%id==i, default=0)
  !    for all i in this%modes_%id and this%modes_%paired_id.
  ! N.B. this is only for optimisation purposes.
  ! As such, indices(i) is only initialised if `i` is in
  !    displacement%vectors%id, this%modes_%id or this%modes_%paired_id.
  ! This reduces a double loop across this%modes_ and displacement%vectors
  !    to two single loops.
  max_mode_id = maxval([ maxval(displacement%vectors%id), &
                       & maxval(this%modes_%id),          &
                       & maxval(this%modes_%paired_id)    ])
  allocate(indices(max_mode_id), stat=ialloc); call err(ialloc)
  do i=1,size(this%modes_)
    indices(this%modes_(i)%id) = 0
    indices(this%modes_(i)%paired_id) = 0
  enddo
  do i=1,size(displacement%vectors)
    indices(displacement%vectors(i)%id) = i
  enddo
  
  ! Match each mode in the monomial with the matching mode(s)
  !    in the displacement.
  allocate( displacements(size(this)),        &
          & paired_displacements(size(this)), &
          & force_is_zero(size(this)),        &
          & stat=ialloc); call err(ialloc)
  value_is_zero = 0
  no_modes = 0
  do i=1,size(this)
    associate (mode => this%modes_(i))
      j = indices(mode%id)
      if (j/=0) then
        displacements(i) = displacement%vectors(j)
      else
        displacements(i) = ComplexSingleDisplacement( &
                            & mode%id,                &
                            & cmplx(0.0_dp,0.0_dp,dp) )
      endif
      if (mode%id==mode%paired_id) then
        if (j==0 .and. mode%power/=0) then
          if (value_is_zero/=0) then
            ! If U_i is zero for more than one i, then
            !    prod_{j/=i}[ {U_j}^{n_j} ] is always zero,
            !    so all derivatives are zero.
            output = ComplexModeForce([ComplexSingleForce::])
            return
          endif
          value_is_zero = i
        endif
        force_is_zero(i) = j==0 .and. mode%power/=1
        if (.not. force_is_zero(i)) then
          no_modes = no_modes+1
        endif
      else
        k = indices(mode%paired_id)
        if (k/=0) then
          paired_displacements(i) = displacement%vectors(k)
        else
          paired_displacements(i) = ComplexSingleDisplacement( &
                                     & mode%paired_id,         &
                                     & cmplx(0.0_dp,0.0_dp,dp) )
        endif
        if (      (j==0 .and. mode%power/=0)        &
           & .or. (k==0 .and. mode%paired_power/=0) ) then
          if (value_is_zero/=0) then
            ! If U_i is zero for more than one i, then
            !    prod_{j/=i}[ {U_j}^{n_j} ] is always zero,
            !    so all derivatives are zero.
            output = ComplexModeForce([ComplexSingleForce::])
            return
          endif
          value_is_zero = i
        endif
        force_is_zero(i) = (j==0 .and. mode%power/=1) &
                   & .and. (k==0 .and. mode%paired_power/=1)
        if (.not. force_is_zero(i)) then
          no_modes = no_modes+2
        endif
      endif
      if (force_is_zero(i) .and. value_is_zero==i) then
        ! U_i and the derivative are zero along mode i,
        !    so all derivatives are zero.
        output = ComplexModeForce([ComplexSingleForce::])
        return
      endif
    end associate
  enddo
  
  i = value_is_zero
  if (i/=0) then
    ! Only one U_i is zero.
    associate (mode => this%modes_(i))
    if (mode%id==mode%paired_id) then
      ids = [mode%id]
      forces = mode%force_complex(displacements(i))
    else
      ids = [mode%id, mode%paired_id]
      forces = mode%force_complex( displacements(i),       &
                                 & paired_displacements(i) )
    endif
    end associate
    do j=1,size(this)
      if (j/=i) then
        associate (mode => this%modes_(j))
          if (mode%id==mode%paired_id) then
            forces = forces &
                 & * mode%energy_complex(displacements(j))
          else
            forces = forces &
                 & * mode%energy_complex( displacements(j),       &
                 &                        paired_displacements(j) )
          endif
        end associate
      endif
    enddo
  else
    ! The force is non-zero along multiple modes.
    allocate( ids(no_modes),    &
            & forces(no_modes), &
            & stat=ialloc); call err(ialloc)
    forces = 1.0_dp
    j = 0
    do i=1,size(this)
      associate (mode => this%modes_(i))
        if (mode%id==mode%paired_id) then
          energy = mode%energy_complex(displacements(i))
        else
          energy = mode%energy_complex( displacements(i),       &
                                      & paired_displacements(i) )
        endif
        
        if (force_is_zero(i)) then
          forces = forces * energy
        else
          if (mode%id==mode%paired_id) then
            j = j+1
            ids(j) = mode%id
            forces(j:j) = forces(j:j) * mode%force_complex(displacements(i))
            forces(:j-1) = forces(:j-1) * energy
            forces(j+1:) = forces(j+1:) * energy
          else
            j = j+2
            ids(j-1) = mode%id
            ids(j) = mode%paired_id
            forces(j-1:j) = mode%force_complex( displacements(i),       &
                                              & paired_displacements(i) )
            forces(:j-2) = forces(:j-2) * energy
            forces(j+1:) = forces(j+1:) * energy
          endif
        endif
      end associate
    enddo
  endif
  
  ! Construct output from components along each mode.
  output = ComplexModeForce(ComplexSingleForce(ids, forces*this%coefficient))
end procedure

module procedure force_RealModeDisplacement_ComplexPolynomial
  type(RealSingleForce) :: zero_force(0)
  
  if (size(this)==0) then
    output = RealModeForce(zero_force)
  else
    output = sum(this%terms%force(displacement))
  endif
end procedure

module procedure force_ComplexModeDisplacement_ComplexPolynomial
  type(ComplexSingleForce) :: zero_force(0)
  
  if (size(this)==0) then
    output = ComplexModeForce(zero_force)
  else
    output = sum(this%terms%force(displacement))
  endif
end procedure

module procedure harmonic_expectation_ComplexPolynomial
  output = sum(this%terms%harmonic_expectation( frequency,      &
                                              & thermal_energy, &
                                              & supercell_size  ))
end procedure

module procedure harmonic_expectation_ComplexMonomial
  output = real( this%coefficient                                             &
       &       * product(this%modes_%harmonic_expectation( frequency,         &
       &                                                   thermal_energy,    &
       &                                                   supercell_size  )) )
end procedure

module procedure harmonic_expectation_ComplexUnivariate
  integer :: n
  
  real(dp) :: factor
  
  ! If the mode is its own conjugate, then
  !            f(n) (1+e^(-w/T))^n
  ! <(u)^2n> = -------------------
  !            (2Nw(1-e^(-w/T)))^n
  ! where f(n) is the odd factorial.
  
  ! If the mode is not its own conjugate, then
  !               n! (1+e^(-w/T))^n
  ! <(u+u-)^n> = -------------------
  !              (2Nw(1-e^(-w/T)))^n
  
  ! Calculate n, and return if the answer is 0 or 1.
  if (this%id==this%paired_id) then
    if (modulo(this%power,2)==1) then
      output = 0
      return
    endif
    n = this%power/2
  else
    if (this%power/=this%paired_power) then
      output = 0
      return
    endif
    n = this%power
  endif
  
  ! <X^0> = <1> = 1.
  if (n==0) then
    output = 1
    return
  endif
  
  ! Calculate (1+e^(-w/T))/(1-e^(-w/T)).
  if (frequency < 1e-10_dp*thermal_energy) then
    ! High temperature regime (N.B. this is unbounded for w=0).
    ! If w/T<1e-10, (w/T)^2 < 1e-20, smaller than floating point error.
    !
    ! e^(-w/T) = 1 - w/T + 0.5(w/T)^2 + ...
    !
    ! (1+e^(-w/T))/(1-e^(-w/T)) = (2 - w/T + ...) / (w/T - 0.5(w/T)^2 + ...)
    !                           = (2T/w - 1 + ...) / (1 - 0.5w/T + ...)
    !                           = (2T/w - 1 + ...) * (1 + 0.5w/T + ...)
    !                           = 2T/w + 0 + ...
    factor = 2*thermal_energy/frequency
  elseif (frequency < 23*thermal_energy) then
    ! Normal temperature regime.
    factor = (1+exp(-frequency/thermal_energy)) &
         & / (1-exp(-frequency/thermal_energy))
  elseif (frequency < 690*thermal_energy) then
    ! Low temperature regime.
    ! If w/T>23, exp(-w/T)<1e-9, so exp(-2w/T) < 1e-20.
    ! 
    ! (1+e^(-w/T))/(1-e^(-w/T)) = (1 + e^(-w/T)) * (1 + e^(-w/T) + ...)
    !                           = 1 + 2e^(-w/T) + ...
    factor = 1+2*exp(-frequency/thermal_energy)
  else
    ! Very low temperature regime.
    ! If w/T>690, exp(-w/T)<1e-300, smaller than floating point can store.
    factor = 1
  endif
  
  ! Calculate the output.
  if (this%id==this%paired_id) then
    output = odd_factorial(n) &
         & * (factor/(2.0_dp*supercell_size*frequency))**n
  else
    output = factorial(n) &
         & * (factor/(2.0_dp*supercell_size*frequency))**n
  endif
end procedure

module procedure multiply_ComplexMonomial_real
  output = ComplexMonomial(this%coefficient*that, this%modes_)
end procedure

module procedure multiply_real_ComplexMonomial
  output = ComplexMonomial(this*that%coefficient, that%modes_)
end procedure

module procedure multiply_ComplexMonomial_complex
  output = ComplexMonomial(this%coefficient*that, this%modes_)
end procedure

module procedure multiply_complex_ComplexMonomial
  output = ComplexMonomial(this*that%coefficient, that%modes_)
end procedure

module procedure multiply_ComplexPolynomial_real
  output = ComplexPolynomial(this%terms * that)
end procedure

module procedure multiply_real_ComplexPolynomial
  output = ComplexPolynomial(this * that%terms)
end procedure

module procedure multiply_ComplexPolynomial_complex
  output = ComplexPolynomial(this%terms * that)
end procedure

module procedure multiply_complex_ComplexPolynomial
  output = ComplexPolynomial(this * that%terms)
end procedure

module procedure divide_ComplexMonomial_real
  output = ComplexMonomial(this%coefficient/that, this%modes_)
end procedure

module procedure divide_ComplexMonomial_complex
  output = ComplexMonomial(this%coefficient/that, this%modes_)
end procedure

module procedure divide_ComplexPolynomial_real
  output = ComplexPolynomial(this%terms / that)
end procedure

module procedure divide_ComplexPolynomial_complex
  output = ComplexPolynomial(this%terms / that)
end procedure

module procedure multiply_ComplexMonomialable_ComplexMonomialable
  complex(dp)                          :: coefficient
  type(ComplexUnivariate), allocatable :: modes(:)
  
  type(ComplexMonomial) :: this_monomial
  type(ComplexMonomial) :: that_monomial
  
  integer :: i_this,i_that,i_out,ialloc
  
  this_monomial = this%to_ComplexMonomial()
  that_monomial = that%to_ComplexMonomial()
  
  coefficient = this_monomial%coefficient * that_monomial%coefficient
  
  if (size(this_monomial)==0) then
    modes = that_monomial%modes_
  elseif (size(that_monomial)==0) then
    modes = this_monomial%modes_
  else
    i_this = 1
    i_that = 1
    i_out = 0
    allocate( modes(size(this_monomial)+size(that_monomial)), &
            & stat=ialloc); call err(ialloc)
    do while(i_this<=size(this_monomial) .or. i_that<=size(that_monomial))
      i_out = i_out + 1
      if (i_this>size(this_monomial)) then
        modes(i_out) = that_monomial%modes_(i_that)
        i_that = i_that + 1
      elseif (i_that>size(that_monomial)) then
        modes(i_out) = this_monomial%modes_(i_this)
        i_this = i_this + 1
      elseif ( this_monomial%modes_(i_this)%id == &
             & that_monomial%modes_(i_that)%id    ) then
        modes(i_out) = ComplexUnivariate(                             &
           & id           = this_monomial%modes_(i_this)%id,          &
           & paired_id    = this_monomial%modes_(i_this)%paired_id,   &
           & power        = this_monomial%modes_(i_this)%power        &
           &              + that_monomial%modes_(i_that)%power,       &
           & paired_power = this_monomial%modes_(i_this)%paired_power &
           &              + that_monomial%modes_(i_that)%paired_power )
        i_this = i_this + 1
        i_that = i_that + 1
      elseif ( this_monomial%modes_(i_this)%id < &
             & that_monomial%modes_(i_that)%id   ) then
        modes(i_out) = this_monomial%modes_(i_this)
        i_this = i_this + 1
      elseif ( this_monomial%modes_(i_this)%id > &
             & that_monomial%modes_(i_that)%id   ) then
        modes(i_out) = that_monomial%modes_(i_that)
        i_that = i_that + 1
      else
        call err()
      endif
    enddo
    modes = modes(:i_out)
  endif
  
  output = ComplexMonomial(coefficient, modes)
end procedure

module procedure multiply_ComplexPolynomial_ComplexMonomialable
  output = ComplexPolynomial(this%terms * that)
end procedure

module procedure multiply_ComplexMonomialable_ComplexPolynomial
  output = ComplexPolynomial(this * that%terms)
end procedure

module procedure add_ComplexPolynomialable_ComplexPolynomialable
  type(ComplexPolynomial) :: this_polynomial
  type(ComplexPolynomial) :: that_polynomial
  
  integer :: no_terms
  
  integer :: i,j,ialloc
  
  this_polynomial = this%to_ComplexPolynomial()
  that_polynomial = that%to_ComplexPolynomial()
  
  allocate( output%terms(size(this_polynomial)+size(that_polynomial)), &
          & stat=ialloc); call err(ialloc)
  output%terms(:size(this_polynomial)) = this_polynomial%terms
  no_terms = size(this_polynomial)
  do i=1,size(that_polynomial)
    j = first_equivalent( this_polynomial%terms,     &
                        & that_polynomial%terms(i),  &
                        & compare_complex_monomials, &
                        & default=0                  )
    if (j==0) then
      no_terms = no_terms + 1
      output%terms(no_terms) = that_polynomial%terms(i)
    else
      output%terms(j)%coefficient = output%terms(j)%coefficient &
                                & + that_polynomial%terms(i)%coefficient
    endif
  enddo
  output%terms = output%terms(:no_terms)
end procedure

module procedure negative_ComplexMonomial
  output = ComplexMonomial( modes       = this%modes_,      &
                          & coefficient = -this%coefficient )
end procedure

module procedure negative_ComplexPolynomial
  output = ComplexPolynomial(-this%terms)
end procedure

module procedure subtract_ComplexPolynomialable_ComplexPolynomialable
  output = this + (-that%to_ComplexPolynomial())
end procedure

module procedure sum_ComplexPolynomialables
  type(ComplexMonomial) :: zero_monomial(0)
  
  integer :: i
  
  if (size(input)==0) then
    output = ComplexPolynomial(zero_monomial)
  else
    output = input(1)%to_ComplexPolynomial()
    do i=2,size(input)
      output = output + input(i)
    enddo
  endif
end procedure

module procedure select_modes_ComplexUnivariate
  if (input%id==input%paired_id) then
    output = [modes(first(modes%id==input%id))]
  else
    output = [ modes(first(modes%id==input%id)),       &
             & modes(first(modes%id==input%paired_id)) ]
  endif
end procedure

module procedure select_modes_ComplexUnivariates
  integer :: i,ialloc
  
  allocate(output(0), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    output = [output, select_modes(input(i), modes)]
  enddo
end procedure

module procedure select_displacements_ComplexUnivariate
  if (input%id==input%paired_id) then
    output = [displacements(first(displacements%id==input%id))]
  else
    output = [ displacements(first(displacements%id==input%id)),       &
             & displacements(first(displacements%id==input%paired_id)) ]
  endif
end procedure

module procedure select_displacements_ComplexUnivariates
  integer :: i,ialloc
  
  allocate(output(0), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    output = [output, select_displacements(input(i), displacements)]
  enddo
end procedure

module procedure select_forces_ComplexUnivariate
  if (input%id==input%paired_id) then
    output = [forces(first(forces%id==input%id))]
  else
    output = [ forces(first(forces%id==input%id)),       &
             & forces(first(forces%id==input%paired_id)) ]
  endif
end procedure

module procedure select_forces_ComplexUnivariates
  integer :: i,ialloc
  
  allocate(output(0), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    output = [output, select_forces(input(i), forces)]
  enddo
end procedure

module procedure compare_complex_monomials
  select type(this); type is(ComplexMonomial)
    select type(that); type is(ComplexMonomial)
      if (size(this)/=size(that)) then
        output = .false.
      else
        output = all( this%modes_%id==that%modes_%id .and.               &
                    & this%modes_%power==that%modes_%power .and.         &
                    & this%modes_%paired_power==that%modes_%paired_power )
      endif
    end select
  end select
end procedure

module procedure read_ComplexUnivariate
  type(String), allocatable :: line(:)
  integer                   :: id
  integer                   :: paired_id
  integer                   :: power
  integer                   :: paired_power
  
  select type(this); type is(ComplexUnivariate)
    ! If id=paired_id=5 and power=paired_power=3 then:
    !    input = '(u5^3)'.
    ! If id=5, paired_id=7, power=3 and paired_power=4 then:
    !    input = '(u5^3*u7^4)'.
    
    ! Strip off brackets, and split into mode and paired mode.
    line = split_line( slice(input,2,len(input)-1), &
                     & delimiter='*')
    
    if (size(line)==1) then
      ! line = [ 'u5^3' ]
      ! ID = paired ID.
      
      ! Split into ID and power.
      line = split_line(line(1), delimiter='^')
      id = int(slice(line(1),2,len(line(1))))
      power = int(line(2))
      paired_id = id
      paired_power = power
    elseif (size(line)==2) then
      ! line = [ 'u5^3', 'u7^4' ]
      ! ID /= paired ID.
      
      ! Split into ID, power, paired ID, paired power.
      line = [ split_line(line(1), delimiter='^'), &
             & split_line(line(2), delimiter='^')  ]
      id = int(slice(line(1),2,len(line(1))))
      power = int(line(2))
      paired_id = int(slice(line(3),2,len(line(3))))
      paired_power = int(line(4))
    else
      call print_line(ERROR//': unable to parse ComplexUnivariate.')
      call err()
    endif
    
    this = ComplexUnivariate( id           = id,          &
                            & paired_id    = paired_id,   &
                            & power        = power,       &
                            & paired_power = paired_power )
  end select
end procedure

module procedure write_ComplexUnivariate
  select type(this); type is(ComplexUnivariate)
    if (this%id==this%paired_id) then
      output = '(u'//this%id//'^'//this%power//')'
    else
      output = '(u'//this%id//'^'//this%power// &
             & '*u'//this%paired_id//'^'//this%paired_power//')'
    endif
  end select
end procedure

module procedure new_ComplexUnivariate_String
  call this%read(input)
end procedure

module procedure read_ComplexMonomial
  type(String),            allocatable :: line(:)
  complex(dp)                          :: coefficient
  type(ComplexUnivariate), allocatable :: modes(:)
  
  integer :: i,ialloc
  
  select type(this); type is(ComplexMonomial)
    ! Splitting the input by '*' separates the coefficient and the modes,
    !    but also splits some modes in two.
    line = split_line(input,delimiter='*')
    
    coefficient = cmplx(line(1))
    
    allocate(modes(0), stat=ialloc); call err(ialloc)
    i = 2
    do while (i<=size(line))
      ! Check if line(i) ends in a bracket.
      if (slice(line(i),len(line(i)),len(line(i)))==')') then
        ! line(i) is a mode on its own.
        modes = [modes, ComplexUnivariate(line(i))]
        i = i+1
      else
        ! line(i) and line(i+1) together make a mode.
        modes = [modes, ComplexUnivariate(line(i)//'*'//line(i+1))]
        i = i+2
      endif
    enddo
    
    this = ComplexMonomial(coefficient, modes)
  end select
end procedure

module procedure write_ComplexMonomial
  select type(this); type is(ComplexMonomial)
    if (size(this%modes_)>0) then
      output = this%coefficient//'*'//join(this%modes_, delimiter='*')
    else
      output = str(this%coefficient)
    endif
  end select
end procedure

module procedure new_ComplexMonomial_String
  call this%read(input)
end procedure

module procedure read_ComplexPolynomial
  type(String),          allocatable :: terms(:)
  type(String)                       :: plus
  type(ComplexMonomial), allocatable :: monomials(:)
  
  plus = '+'
  select type(this); type is(ComplexPolynomial)
    terms = split_line(input)
    terms = terms(filter(terms/=plus))
    
    monomials = ComplexMonomial(terms)
    
    this = ComplexPolynomial(monomials)
  end select
end procedure

module procedure write_ComplexPolynomial
  select type(this); type is(ComplexPolynomial)
    output = join(this%terms, delimiter=' + ')
  end select
end procedure

module procedure new_ComplexPolynomial_String
  call this%read(input)
end procedure
end submodule
