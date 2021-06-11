!> Provides the implementation of the [[QpointCombinations(type)]] methods.
submodule (caesar_qpoint_combination_module) &
   & caesar_qpoint_combinations_submodule
  use caesar_stars_module
contains

module procedure new_QpointCombinations
  this%power = power
  this%combinations = combinations
end procedure

module procedure read_QpointCombinations
  integer                              :: power
  type(QpointCombination), allocatable :: combinations(:)
  
  select type(this); type is(QpointCombinations)
    power = int(token(input(1), 6))
    combinations = QpointCombination(input(2:))
    this = QpointCombinations(power, combinations)
  class default
    call err()
  end select
end procedure

module procedure write_QpointCombinations
  select type(this); type is(QpointCombinations)
    output = [ 'q-point combinations with power = '//this%power//' :', &
             & str(this%combinations)                                  ]
  class default
    call err()
  end select
end procedure

module procedure new_QpointCombinations_Strings
  call this%read(input)
end procedure

module procedure new_QpointCombinations_StringArray
  call this%read(str(input))
end procedure

module procedure generate_qpoint_combinations
  type :: PowerArray
    type(QpointPower),    allocatable :: powers(:)
    type(FractionVector), allocatable :: wavevectors(:)
  end type
  
  type :: CombinationData
    type(QpointCombination) :: combination
    integer                 :: total_power
    integer                 :: no_qpoints
    type(FractionVector)    :: wavevector
  end type
  
  integer :: max_qpoints
  
  logical              :: qpoints_real
  integer, allocatable :: key(:)
  
  type(PowerArray), allocatable :: powers(:)
  
  type(CombinationData), allocatable :: combinations(:)
  type(CombinationData), allocatable :: new_combinations(:)
  
  integer              :: power
  integer, allocatable :: no_combinations(:)
  
  integer :: i,j,k,l,ialloc
  
  if (max_power<0) then
    call print_line(ERROR//': max_power must be > 0.')
    call err()
  endif
  
  max_qpoints = set_default(max_qpoint_coupling, max_power)
  
  allocate(output(max_power+1), stat=ialloc); call err(ialloc)
  output%power = [(i,i=0,max_power)]
  
  if (size(qpoints)==0) then
    output(1)%combinations = [QpointCombination([QpointPower::])]
    do i=1,max_power
      output(i+1)%combinations = [QpointCombination::]
    enddo
    return
  endif
  
  ! Determine if the q-points are real (id=paired_id) or not.
  ! If the q-points are not real, filter for the q-points with id<paired_id.
  if (all(qpoints%is_paired_qpoint())) then
    qpoints_real = .true.
    key = [(i,i=1,size(qpoints))]
  elseif (.not.any(qpoints%is_paired_qpoint())) then
    qpoints_real = .false.
    key = filter(qpoints%id<qpoints%paired_qpoint_id)
  else
    call print_line(CODE_ERROR//': Q-point stars may only contain q-points &
       &which are all real or all complex.')
    call err()
  endif
  
  ! For each q-point, list the possible QpointPowers, and calculate the
  !    wavevector of each power.
  allocate(powers(size(key)), stat=ialloc); call err(ialloc)
  do i=1,size(key)
    powers(i)%powers = generate_qpoint_powers(qpoints(key(i)), max_power)
    if (qpoints_real) then
      powers(i)%wavevectors = qpoints(key(i))%qpoint*powers(i)%powers%power()
    else
      powers(i)%wavevectors = qpoints(key(i))%qpoint            &
                          & * ( powers(i)%powers%power()        &
                          &   - powers(i)%powers%paired_power() )
    endif
  enddo
  
  ! Initialise `combinations` to the empty combination.
  combinations = [CombinationData(                       &
     & combination = QpointCombination([QpointPower::]), &
     & total_power = 0,                                  &
     & no_qpoints  = 0,                                  &
     & wavevector  = fracvec(zeroes(3))                  )]
  ! Construct all q-point combinations up to max_power.
  ! After iteration `i` of the loop, `combinations` contains the q-point
  !    combinations corresponding to the first `i` q-points, with
  !    total_power <= max_power and no_qpoints <= max_qpoints.
  ! If conserve_momentum is true then the final iteration prunes combinations
  !    which do not conserve momentum.
  do i=1,size(powers)
    ! Allocate `new_combinations` to the maximum size possibly needed.
    if (allocated(new_combinations)) then
      deallocate(new_combinations, stat=ialloc); call err(ialloc)
    endif
    allocate( new_combinations(                                   &
            &      size(combinations)                             &
            &    * (1+size(powers(i)%powers))                     &
            &    - count(combinations%no_qpoints==max_qpoints)    &
            &    * size(powers(i)%powers)                      ), &
            & stat=ialloc); call err(ialloc)
    
    ! Loop over the combinations from the previous `i` loop iteration,
    !    and the powers at q-point `i`.
    l = 0
    do j=1,size(combinations)
      l = l+1
      new_combinations(l) = combinations(j)
      
      if (i==size(powers)) then
        if ( set_default(conserve_momentum,.false.)             &
           & .and. .not. is_int(new_combinations(l)%wavevector) ) then
          l = l-1
        endif
      endif
      
      if (combinations(j)%no_qpoints==max_qpoints) then
        cycle
      endif
      
      do k=1,size(powers(i)%powers)
        if ( combinations(j)%total_power       &
         & + powers(i)%powers(k)%total_power() &
         & > max_power                         ) then
          exit
        endif
        
        l = l+1
        new_combinations(l)%combination%qpoints_ = [ &
             & combinations(j)%combination%qpoints_, &
             & powers(i)%powers(k)                   ]
        new_combinations(l)%total_power =  &
           &   combinations(j)%total_power &
           & + powers(i)%powers(k)%total_power()
        new_combinations(l)%no_qpoints = new_combinations(l)%no_qpoints + 1
        new_combinations(l)%wavevector = combinations(j)%wavevector &
                                     & + powers(i)%wavevectors(k)
        if (i==size(powers)) then
          if ( set_default(conserve_momentum,.false.)             &
             & .and. .not. is_int(new_combinations(l)%wavevector) ) then
            l = l-1
          endif
        endif
      enddo
    enddo
    
    combinations = new_combinations(:l)
  enddo
  
  ! Separate the combinations by total power.
  allocate(no_combinations(max_power+1), stat=ialloc); call err(ialloc)
  do power=0,max_power
    no_combinations(power+1) = count(combinations%total_power==power)
    allocate( output(power+1)%combinations(no_combinations(power+1)), &
            & stat=ialloc); call err(ialloc)
  enddo
  do i=1,size(combinations)
    power = combinations(i)%total_power
    output(power+1)%combinations(no_combinations(power+1)) = &
       & combinations(i)%combination
    no_combinations(power+1) = no_combinations(power+1)-1
  enddo
  
  if (any(no_combinations/=0)) then
    call print_line(CODE_ERROR//': Inconsistency counting combinations.')
    call err()
  endif
end procedure
end submodule
