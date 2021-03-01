submodule (caesar_real_polynomial_module) caesar_real_polynomial_submodule
  use caesar_normal_mode_module
contains

module procedure new_RealUnivariate
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

module procedure new_RealUnivariate_RealMode
  if (present(paired_power)) then
    if (mode%id==mode%paired_id .and. power/=paired_power) then
      call print_line(ERROR//': Mode is its own pair, but power does not &
         & match paired_power.')
      call err()
    endif
    this = RealUnivariate( id           = mode%id,        &
                         & paired_id    = mode%paired_id, &
                         & power        = power,          &
                         & paired_power = paired_power    )
  else
    if (mode%id==mode%paired_id) then
      this = RealUnivariate( id           = mode%id, &
                           & paired_id    = mode%id, &
                           & power        = power,   &
                           & paired_power = power    )
    else
      this = RealUnivariate( id           = mode%id,        &
                           & paired_id    = mode%paired_id, &
                           & power        = power,          &
                           & paired_power = 0               )
    endif
  endif
end procedure

module procedure new_RealMonomial
  this%coefficient = coefficient
  this%modes_      = modes(sort(modes%id))
end procedure

module procedure new_RealMonomial_RealMonomialable
  this = input%to_RealMonomial()
end procedure

module procedure new_RealPolynomial
  this%terms = terms
end procedure

module procedure new_RealPolynomial_RealPolynomialable
  this = input%to_RealPolynomial()
end procedure

module procedure to_RealMonomial_RealUnivariate
  select type(this); type is(RealUnivariate)
    output = RealMonomial( coefficient = 1.0_dp, &
                         & modes       = [this])
  class default
    call err()
  end select
end procedure

module procedure to_RealPolynomial_RealUnivariate
  select type(this); type is(RealUnivariate)
    output = RealPolynomial([this%to_RealMonomial()])
  class default
    call err()
  end select
end procedure

module procedure to_RealMonomial_RealMonomial
  select type(this); type is(RealMonomial)
    output = this
  class default
    call err()
  end select
end procedure

module procedure to_RealPolynomial_RealMonomial
  select type(this); type is(RealMonomial)
    output = RealPolynomial([this])
  class default
    call err()
  end select
end procedure

module procedure to_RealPolynomial_RealPolynomial
  select type(this); type is(RealPolynomial)
    output = this
  class default
    call err()
  end select
end procedure

module procedure size_RealMonomial
  output = size(this%modes_)
end procedure

module procedure size_RealPolynomial
  output = size(this%terms)
end procedure

module procedure mode_RealMonomial
  output = this%modes_(index)
end procedure

module procedure modes_RealMonomial
  if (present(indices)) then
    output = this%modes_(indices)
  else
    output = this%modes_
  endif
end procedure

module procedure id_RealMonomial
  output = this%modes_(index)%id
end procedure

module procedure ids_RealMonomial
  if (present(indices)) then
    output = this%modes_(indices)%id
  else
    output = this%modes_%id
  endif
end procedure

module procedure paired_id_RealMonomial
  output = this%modes_(index)%paired_id
end procedure

module procedure paired_ids_RealMonomial
  if (present(indices)) then
    output = this%modes_(indices)%paired_id
  else
    output = this%modes_%paired_id
  endif
end procedure

module procedure power_RealMonomial
  output = this%modes_(index)%power
end procedure

module procedure powers_RealMonomial
  if (present(indices)) then
    output = this%modes_(indices)%power
  else
    output = this%modes_%power
  endif
end procedure

module procedure paired_power_RealMonomial
  output = this%modes_(index)%paired_power
end procedure

module procedure paired_powers_RealMonomial
  if (present(indices)) then
    output = this%modes_(indices)%paired_power
  else
    output = this%modes_%paired_power
  endif
end procedure

module procedure simplify_RealMonomial
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

module procedure simplify_RealPolynomial
  type(RealMonomial), allocatable :: equivalent_monomials(:)
  type(RealMonomial), allocatable :: monomials(:)
  
  integer :: i
  
  call this%terms%simplify()
  
  ! Add together monomials with the same powers.
  monomials = this%terms(set(this%terms, compare_real_monomials))
  do i=1,size(monomials)
    equivalent_monomials = this%terms(                          &
       & filter(this%terms,compare_real_monomials,monomials(i)) )
    monomials(i)%coefficient = sum(equivalent_monomials%coefficient)
  enddo
end procedure

module procedure total_power_RealUnivariate
  if (this%id==this%paired_id) then
    output = this%power
  else
    output = this%power + this%paired_power
  endif
end procedure

module procedure total_power_RealMonomial
  output = sum(this%modes_%total_power())
end procedure

module procedure energy_RealSingleDisplacement_RealUnivariate
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

module procedure energy_ComplexSingleDisplacement_RealUnivariate
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
    ! Convert from complex to real co-ordinates.
    ! x_c = (x_+ + x_-) / sqrt(2)
    ! x_s = (x_+ - x_-) / sqrt(2)i
    magnitude = 0
    paired_magnitude = 0
    if (present(displacement)) then
      ! x_c : x_+ / sqrt(2)
      magnitude = magnitude &
              & + displacement%magnitude &
              & / sqrt(2.0_dp)
      ! x_s : x_+ / sqrt(2)i
      paired_magnitude = paired_magnitude       &
                     & + displacement%magnitude &
                     & / cmplx(0.0_dp,sqrt(2.0_dp),dp)
    endif
    if (present(paired_displacement)) then
      ! x_c : x_- / sqrt(2)
      magnitude = magnitude &
              & + paired_displacement%magnitude &
              & / sqrt(2.0_dp)
      ! x_s : - x_- / sqrt(2)i
      paired_magnitude = paired_magnitude              &
                     & - paired_displacement%magnitude &
                     & / cmplx(0.0_dp,sqrt(2.0_dp),dp)
    endif
  endif
  
  output = magnitude**this%power
  if (this%id/=this%paired_id) then
    output = output * paired_magnitude**this%paired_power
  endif
end procedure

module procedure energy_RealModeDisplacement_RealMonomial
  integer :: i,j,k
  
  output = this%coefficient
  
  ! Loop over the modes in this monomial,
  !    find the mode (and paired mode if relevant) in the displacement
  !    corresponding to each mode,
  !    and calculate the energy term by term.
  do i=1,size(this)
    associate (mode => this%modes_(i))
      if (mode%id==mode%paired_id) then
        if (mode%power/=0) then
          j = first(displacement%vectors%id==mode%id, default=0)
          
          if (j==0) then
            ! If the mode is not present in the displacement,
            !    then the displacement along that mode is zero.
            ! As such, the monomial is zero. (0^n=0 if n>0).
            output = 0
            return
          else
            output = output * mode%energy_real(       &
               & displacement=displacement%vectors(j) )
          endif
        endif
      else
        if (mode%power/=0 .or. mode%paired_power/=0) then
          j = first(displacement%vectors%id==mode%id, default=0)
          k = first(displacement%vectors%id==mode%paired_id, default=0)
          
          if (      (j==0 .and. mode%power/=0)        &
             & .or. (k==0 .and. mode%paired_power/=0) ) then
            ! If either mode is not present in the displacement,
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

module procedure energy_ComplexModeDisplacement_RealMonomial
  integer :: i,j,k
  
  output = this%coefficient
  
  ! Loop over the modes in this monomial,
  !    find the mode (and paired mode if relevant) in the displacement
  !    corresponding to each mode,
  !    and calculate the energy term by term.
  do i=1,size(this)
    associate (mode => this%modes_(i))
      if (mode%id==mode%paired_id) then
        if (mode%power/=0) then
          j = first(displacement%vectors%id==mode%id, default=0)
          
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
          j = first(displacement%vectors%id==mode%id, default=0)
          k = first(displacement%vectors%id==mode%paired_id, default=0)
          
          if (j==0 .and. k==0) then
            ! If both modes are not present in the displacement,
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

module procedure energy_RealModeDisplacement_RealPolynomial
  output = sum(this%terms%energy(displacement))
end procedure

module procedure energy_ComplexModeDisplacement_RealPolynomial
  output = sum(this%terms%energy(displacement))
end procedure

module procedure force_RealSingleDisplacement_RealUnivariate
  real(dp) :: magnitude
  real(dp) :: paired_magnitude
  
  real(dp) :: values(2)
  real(dp) :: derivatives(2)
  
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

module procedure force_ComplexSingleDisplacement_RealUnivariate
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
    ! Convert from complex to real co-ordinates.
    ! x_c = (x_+ + x_-) / sqrt(2)
    ! x_s = (x_+ - x_-) / sqrt(2)i
    magnitude = 0
    paired_magnitude = 0
    if (present(displacement)) then
      ! x_c : x_+ / sqrt(2)
      magnitude = magnitude &
              & + displacement%magnitude &
              & / sqrt(2.0_dp)
      ! x_s : x_+ / sqrt(2)i
      paired_magnitude = paired_magnitude       &
                     & + displacement%magnitude &
                     & / cmplx(0.0_dp,sqrt(2.0_dp),dp)
    endif
    if (present(paired_displacement)) then
      ! x_c : x_- / sqrt(2)
      magnitude = magnitude &
              & + paired_displacement%magnitude &
              & / sqrt(2.0_dp)
      ! x_s : - x_- / sqrt(2)i
      paired_magnitude = paired_magnitude              &
                     & - paired_displacement%magnitude &
                     & / cmplx(0.0_dp,sqrt(2.0_dp),dp)
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
  
  ! Construct the output in real co-ordinates...
  output = [-derivatives(1)*values(2), -derivatives(2)*values(1)]
  ! ... then convert back to complex co-ordinates.
  ! f_+ = (f_c - if_s) / sqrt(2)
  ! f_- = (f_c + if_s) / sqrt(2)
  output = [ (output(1)-output(2)*cmplx(0,1,dp))/sqrt(2.0_dp), &
           & (output(1)+output(2)*cmplx(0,1,dp))/sqrt(2.0_dp)  ]
end procedure

module procedure force_RealModeDisplacement_RealMonomial
  logical, allocatable :: value_is_zero(:)
  logical, allocatable :: force_is_zero(:)
  
  type(RealSingleDisplacement), allocatable :: displacements(:)
  type(RealSingleDisplacement), allocatable :: paired_displacements(:)
  
  real(dp)              :: energy
  integer,  allocatable :: ids(:)
  real(dp), allocatable :: forces(:)
  
  integer :: i,j,k,ialloc
  
  integer :: no_modes
  
  ! Match each mode in the monomial with the matching mode(s)
  !    in the displacement.
  allocate( displacements(size(this)),        &
          & paired_displacements(size(this)), &
          & value_is_zero(size(this)),        &
          & force_is_zero(size(this)),        &
          & stat=ialloc); call err(ialloc)
  no_modes = 0
  do i=1,size(this)
    associate (mode => this%modes_(i))
      j = first(displacement%vectors%id==mode%id, default=0)
      if (j/=0) then
        displacements(i) = displacement%vectors(j)
      else
        displacements(i) = RealSingleDisplacement(mode%id, 0.0_dp)
      endif
      if (mode%id==mode%paired_id) then
        value_is_zero(i) = j==0 .and. mode%power/=0
        force_is_zero(i) = j==0 .and. mode%power/=1
        if (.not. force_is_zero(i)) then
          no_modes = no_modes+1
        endif
      else
        k = first(displacement%vectors%id==mode%paired_id, default=0)
        if (k/=0) then
          paired_displacements(i) = displacement%vectors(k)
        else
          paired_displacements(i) = RealSingleDisplacement( mode%paired_id, &
                                                          & 0.0_dp          )
        endif
        value_is_zero(i) = (j==0 .and. mode%power/=0) &
                    & .or. (k==0 .and. mode%paired_power/=0)
        force_is_zero(i) = (j==0 .and. mode%power/=1) &
                   & .and. (k==0 .and. mode%paired_power/=1)
        if (.not. force_is_zero(i)) then
          no_modes = no_modes+2
        endif
      endif
    end associate
  enddo
  
  if (count(value_is_zero)>1) then
    ! If U_i is zero for more than one i, then
    !    prod_{j/=i}[ {U_j}^{n_j} ] is always zero,
    !    so all derivatives are zero.
    allocate( ids(0),    &
            & forces(0), &
            & stat=ialloc); call err(ialloc)
  else
    i = first(value_is_zero, default=0)
    if (i/=0) then
      ! Only one U_i is zero.
      if (force_is_zero(i)) then
        ! The derivative is also zero along mode i.
        allocate( ids(0),    &
                & forces(0), &
                & stat=ialloc); call err(ialloc)
      else
        ! The derivative is not zero along mode i.
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
      endif
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
  endif
  
  ! Construct output from components along each mode.
  output = RealModeForce(RealSingleForce(ids, forces*this%coefficient))
end procedure

module procedure force_ComplexModeDisplacement_RealMonomial
  logical, allocatable :: value_is_zero(:)
  logical, allocatable :: force_is_zero(:)
  
  type(ComplexSingleDisplacement), allocatable :: displacements(:)
  type(ComplexSingleDisplacement), allocatable :: paired_displacements(:)
  
  complex(dp)              :: energy
  integer,     allocatable :: ids(:)
  complex(dp), allocatable :: forces(:)
  
  integer :: i,j,k,ialloc
  
  integer :: no_modes
  
  ! Match each mode in the monomial with the matching mode(s)
  !    in the displacement.
  allocate( displacements(size(this)),        &
          & paired_displacements(size(this)), &
          & value_is_zero(size(this)),        &
          & force_is_zero(size(this)),        &
          & stat=ialloc); call err(ialloc)
  no_modes = 0
  do i=1,size(this)
    associate (mode => this%modes_(i))
      j = first(displacement%vectors%id==mode%id, default=0)
      if (j/=0) then
        displacements(i) = displacement%vectors(j)
      else
        displacements(i) = ComplexSingleDisplacement( &
                            & mode%id,                &
                            & cmplx(0.0_dp,0.0_dp,dp) )
      endif
      if (mode%id==mode%paired_id) then
        value_is_zero(i) = j==0 .and. mode%power/=0
        force_is_zero(i) = j==0 .and. mode%power/=1
        if (.not. force_is_zero(i)) then
          no_modes = no_modes+1
        endif
      else
        k = first(displacement%vectors%id==mode%paired_id, default=0)
        if (k/=0) then
          paired_displacements(i) = displacement%vectors(k)
        else
          paired_displacements(i) = ComplexSingleDisplacement( &
                                     & mode%paired_id,         &
                                     & cmplx(0.0_dp,0.0_dp,dp) )
        endif
        value_is_zero(i) = (j==0 .and. k==0) &
                   & .and. (mode%power/=0 .or. mode%paired_power/=0)
        force_is_zero(i) = (j==0 .and. mode%power/=1) &
                   & .and. (k==0 .and. mode%paired_power/=1)
        if (.not. force_is_zero(i)) then
          no_modes = no_modes+2
        endif
      endif
    end associate
  enddo
  
  if (count(value_is_zero)>1) then
    ! If U_i is zero for more than one i, then
    !    prod_{j/=i}[ {U_j}^{n_j} ] is always zero,
    !    so all derivatives are zero.
    allocate( ids(0),    &
            & forces(0), &
            & stat=ialloc); call err(ialloc)
  else
    i = first(value_is_zero, default=0)
    if (i/=0) then
      ! Only one U_i is zero.
      if (force_is_zero(i)) then
        ! The derivative is also zero along mode i.
        allocate( ids(0),    &
                & forces(0), &
                & stat=ialloc); call err(ialloc)
      else
        ! The derivative is not zero along mode i.
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
                forces = forces                                       &
                     & * mode%energy_complex( displacements(j),       &
                     &                        paired_displacements(j) )
              endif
            end associate
          endif
        enddo
      endif
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
  endif
  
  ! Construct output from components along each mode.
  output = ComplexModeForce(ComplexSingleForce(ids, forces*this%coefficient))
end procedure

module procedure force_RealModeDisplacement_RealPolynomial
  type(RealSingleForce) :: zero_force(0)
  
  if (size(this)==0) then
    output = RealModeForce(zero_force)
  else
    output = sum(this%terms%force(displacement))
  endif
end procedure

module procedure force_ComplexModeDisplacement_RealPolynomial
  type(ComplexSingleForce) :: zero_force(0)
  
  if (size(this)==0) then
    output = ComplexModeForce(zero_force)
  else
    output = sum(this%terms%force(displacement))
  endif
end procedure

module procedure multiply_RealMonomial_real
  output = RealMonomial(this%coefficient*that, this%modes_)
end procedure

module procedure multiply_real_RealMonomial
  output = RealMonomial(this*that%coefficient, that%modes_)
end procedure

module procedure multiply_RealPolynomial_real
  output = RealPolynomial(this%terms * that)
end procedure

module procedure multiply_real_RealPolynomial
  output = RealPolynomial(this * that%terms)
end procedure

module procedure divide_RealMonomial_real
  output = RealMonomial(this%coefficient/that, this%modes_)
end procedure

module procedure divide_RealPolynomial_real
  output = RealPolynomial(this%terms / that)
end procedure

module procedure multiply_RealMonomialable_RealMonomialable
  real(dp)                          :: coefficient
  type(RealUnivariate), allocatable :: modes(:)
  
  type(RealMonomial) :: this_monomial
  type(RealMonomial) :: that_monomial
  
  integer :: i_this,i_that,i_out,ialloc
  
  this_monomial = this%to_RealMonomial()
  that_monomial = that%to_RealMonomial()
  
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
        modes(i_out) = RealUnivariate(                                 &
           & id           = this_monomial%modes_(i_this)%id,           &
           & paired_id    = this_monomial%modes_(i_this)%paired_id,    &
           & power        = this_monomial%modes_(i_this)%power         &
           &              + that_monomial%modes_(i_that)%power,        &
           & paired_power = this_monomial%modes_(i_this)%paired_power  &
           &              + that_monomial%modes_(i_that)%paired_power  )
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
  
  output = RealMonomial(coefficient, modes)
end procedure

module procedure multiply_RealPolynomial_RealMonomialable
  output = RealPolynomial(this%terms * that)
end procedure

module procedure multiply_RealMonomialable_RealPolynomial
  output = RealPolynomial(this * that%terms)
end procedure

module procedure add_RealPolynomialable_RealPolynomialable
  type(RealPolynomial) :: this_polynomial
  type(RealPolynomial) :: that_polynomial
  
  integer :: no_terms
  
  integer :: i,j,ialloc
  
  this_polynomial = this%to_RealPolynomial()
  that_polynomial = that%to_RealPolynomial()
  
  allocate( output%terms(size(this_polynomial)+size(that_polynomial)), &
          & stat=ialloc); call err(ialloc)
  output%terms(:size(this_polynomial)) = this_polynomial%terms
  no_terms = size(this_polynomial)
  do i=1,size(that_polynomial)
    j = first_equivalent( this_polynomial%terms,    &
                        & that_polynomial%terms(i), &
                        & compare_real_monomials,   &
                        & default=0)
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

module procedure negative_RealPolynomialable
  output = this%to_RealPolynomial()
  output%terms%coefficient = - output%terms%coefficient
end procedure

module procedure subtract_RealPolynomialable_RealPolynomialable
  output = this + (-that)
end procedure

module procedure sum_RealPolynomialables
  type(RealMonomial) :: zero_monomial(0)
  
  integer :: i
  
  if (size(input)==0) then
    output = RealPolynomial(zero_monomial)
  else
    output = input(1)%to_RealPolynomial()
    do i=2,size(input)
      output = output + input(i)
    enddo
  endif
end procedure

module procedure select_mode_RealUnivariate
  output = modes(first(modes%id==univariate%id))
end procedure

module procedure select_modes_RealUnivariates
  integer :: i,ialloc
  
  allocate(output(size(univariates)), stat=ialloc); call err(ialloc)
  do i=1,size(univariates)
    output(i) = select_mode(univariates(i), modes)
  enddo
end procedure

module procedure select_displacement_RealUnivariate
  output = displacements(first(displacements%id==univariate%id))
end procedure

module procedure select_displacements_RealUnivariates
  integer :: i,ialloc
  
  allocate(output(size(univariates)), stat=ialloc); call err(ialloc)
  do i=1,size(univariates)
    output(i) = select_displacement(univariates(i), displacements)
  enddo
end procedure

module procedure select_force_RealUnivariate
  output = forces(first(forces%id==univariate%id))
end procedure

module procedure select_forces_RealUnivariates
  integer :: i,ialloc
  
  allocate(output(size(univariates)), stat=ialloc); call err(ialloc)
  do i=1,size(univariates)
    output(i) = select_force(univariates(i), forces)
  enddo
end procedure

module procedure compare_real_monomials
  select type(this); type is(RealMonomial)
    select type(that); type is(RealMonomial)
      if (size(this)/=size(that)) then
        output = .false.
      else
        output = all( this%modes_%id           == that%modes_%id    .and.  &
                    & this%modes_%power        == that%modes_%power .and.  &
                    & this%modes_%paired_power == that%modes_%paired_power )
      endif
    end select
  end select
end procedure

module procedure read_RealUnivariate
  type(String), allocatable :: line(:)
  integer                   :: id
  integer                   :: paired_id
  integer                   :: power
  integer                   :: paired_power
  
  select type(this); type is(RealUnivariate)
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
      call print_line(ERROR//': unable to parse RealUnivariate.')
      call err()
    endif
    
    this = RealUnivariate( id           = id,          &
                         & paired_id    = paired_id,   &
                         & power        = power,       &
                         & paired_power = paired_power )
  end select
end procedure

module procedure write_RealUnivariate
  select type(this); type is(RealUnivariate)
    if (this%id==this%paired_id) then
      output = '(u'//this%id//'^'//this%power//')'
    else
      output = '(u'//this%id//'^'//this%power// &
             & '*u'//this%paired_id//'^'//this%paired_power//')'
    endif
  end select
end procedure

module procedure new_RealUnivariate_String
  call this%read(input)
end procedure

module procedure read_RealMonomial
  type(String),         allocatable :: line(:)
  real(dp)                          :: coefficient
  type(RealUnivariate), allocatable :: modes(:)
  
  integer :: i,ialloc
  
  select type(this); type is(RealMonomial)
    ! Splitting the input by '*' separates the coefficient and the modes,
    !    but also splits some modes in two.
    line = split_line(input,delimiter='*')
    
    coefficient = dble(line(1))
    
    allocate(modes(0), stat=ialloc); call err(ialloc)
    i = 2
    do while (i<=size(line))
      ! Check if line(i) ends in a bracket.
      if (slice(line(i),len(line(i)),len(line(i)))==')') then
        ! line(i) is a mode on its own.
        modes = [modes, RealUnivariate(line(i))]
        i = i+1
      else
        ! line(i) and line(i+1) together make a mode.
        modes = [modes, RealUnivariate(line(i)//'*'//line(i+1))]
        i = i+2
      endif
    enddo
    
    this = RealMonomial(coefficient, modes)
  end select
end procedure

module procedure write_RealMonomial
  select type(this); type is(RealMonomial)
    if (size(this%modes_)>0) then
      output = this%coefficient//'*'//join(this%modes_, delimiter='*')
    else
      output = str(this%coefficient)
    endif
  end select
end procedure

module procedure new_RealMonomial_String
  call this%read(input)
end procedure

module procedure read_RealPolynomial
  type(String),       allocatable :: terms(:)
  type(String)                    :: plus
  type(RealMonomial), allocatable :: monomials(:)
  
  plus = '+'
  select type(this); type is(RealPolynomial)
    terms = split_line(input)
    terms = terms(filter(terms/=plus))
    
    monomials = RealMonomial(terms)
    
    this = RealPolynomial(monomials)
  end select
end procedure

module procedure write_RealPolynomial
  select type(this); type is(RealPolynomial)
    output = join(this%terms, delimiter=' + ')
  end select
end procedure

module procedure new_RealPolynomial_String
  call this%read(input)
end procedure
end submodule
