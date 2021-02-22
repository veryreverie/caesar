submodule (caesar_random_module) caesar_random_submodule
contains

module procedure new_RandomReal
  integer              :: random_size
  integer, allocatable :: random_seeds(:)
  character(10)        :: time
  
  integer :: i,ialloc
  
  ! Set the seed, either to the given seed if present, or to the
  !    current time in milliseconds.
  if (present(seed)) then
    this%random_seed_ = seed
  else
    ! Get the time in seconds (to three decimal places).
    call date_and_time(time=time)
    ! Convert to an integer in milliseconds.
    this%random_seed_ = nint(1000 * dble(time))
  endif
  
  ! Initialise random number generator from the seed.
  call random_seed(size=random_size)
  allocate(random_seeds(random_size), stat=ialloc); call err(ialloc)
  random_seeds = this%random_seed_
  do i=1,size(random_seeds)
    random_seeds(i) = random_seeds(i) + i
  enddo
  call random_seed(put=random_seeds)
  
  ! Set this%initialised_ to .true.
  this%initialised_ = .true.
end procedure

module procedure get_seed
  output = this%random_seed_
end procedure

module procedure random_number_RandomReal
  logical :: log_distributed_
  
  if (.not. this%initialised_) then
    call print_line(CODE_ERROR//': Calling random number generator before it &
       &has been initialised.')
    call err()
  endif
  
  if (present(minimum).neqv.present(maximum)) then
    call print_line(CODE_ERROR//': If one bound is specified, &
       &the other must be too.')
    call err()
  endif
  
  if (present(minimum)) then
    if (minimum>maximum) then
      call print_line(CODE_ERROR//': minimum>maximum.')
      call err()
    elseif (minimum>=maximum) then
      ! minimum=maximum.
      output = minimum
    endif
  endif
  
  log_distributed_ = set_default(log_distributed, .false.)
  
  if (log_distributed_) then
    if (.not. present(minimum)) then
      call print_line(CODE_ERROR//': If a log distribution is &
         &requested then minimum and maximum must be specified.')
      call err()
    elseif (minimum<=0) then
      call print_line(CODE_ERROR//': If a log distribution is &
         &requested then the bounds must be >0.')
      call err()
    endif
  endif
  
  call random_number(output)
  
  if (log_distributed_) then
    output = minimum*(maximum/minimum)**output
  elseif (present(minimum)) then
    output = minimum+(maximum-minimum)*output
  endif
end procedure

module procedure random_numbers_RandomReal
  integer :: i,ialloc
  
  allocate(output(no_numbers), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    output(i) = this%random_number()
  enddo
end procedure
end submodule
