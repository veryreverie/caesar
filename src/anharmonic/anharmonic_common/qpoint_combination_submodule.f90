submodule (caesar_qpoint_combination_module) &
   & caesar_qpoint_combination_submodule
  use caesar_anharmonic_common_module
contains

module procedure new_QpointCombination
  this%qpoints_ = qpoints(sort(qpoints%id()))
end procedure

module procedure qpoints_QpointCombination
  output = this%qpoints_
end procedure

module procedure total_power_QpointCombination
  output = sum(this%qpoints_%total_power())
end procedure

module procedure equality_QpointCombination_QpointCombination
  if (size(this%qpoints_)==size(that%qpoints_)) then
    output = all(this%qpoints_==that%qpoints_)
  else
    output = .false.
  endif
end procedure

module procedure non_equality_QpointCombination_QpointCombination
  output = .not. this==that
end procedure

module procedure lt_QpointCombination_QpointCombination
  integer :: i
  
  if (this%total_power()<that%total_power()) then
    output = .true.
  elseif (this%total_power()>that%total_power()) then
    output = .false.
  else
    do i=1,min(size(this%qpoints_),size(that%qpoints_))
      if (this%qpoints_(i)<that%qpoints_(i)) then
        output = .true.
        return
      elseif (this%qpoints_(i)>that%qpoints_(i)) then
        output = .false.
        return
      endif
      
      output = .false.
    enddo
  endif
end procedure

module procedure le_QpointCombination_QpointCombination
  output = this==that .or. this<that
end procedure

module procedure gt_QpointCombination_QpointCombination
  output = .not. this<=that
end procedure

module procedure ge_QpointCombination_QpointCombination
  output = .not. this<that
end procedure

module procedure read_QpointCombination
  type(String),      allocatable :: line(:)
  type(QpointPower), allocatable :: qpoints(:)
  
  integer :: i,ialloc
  
  select type(this); type is(QpointCombination)
    ! Splitting the input by '*' separates the q-points,
    !    but also splits q-point pairs in two.
    line = split_line(input,delimiter='*')
    
    allocate(qpoints(0), stat=ialloc); call err(ialloc)
    i = 1
    do while (i<=size(line))
      ! Check if line(i) ends in a bracket.
      if (slice(line(i),len(line(i)),len(line(i)))==')') then
        ! line(i) is a q-point on its own.
        qpoints = [qpoints, QpointPower(line(i))]
        i = i+1
      else
        ! line(i) and line(i+1) together make a q-point pair.
        qpoints = [qpoints, QpointPower(line(i)//'*'//line(i+1))]
        i = i+2
      endif
    enddo
    
    this = QpointCombination(qpoints)
  class default
    call err()
  end select
end procedure

module procedure write_QpointCombination
  select type(this); type is(QpointCombination)
    output = join(this%qpoints_, delimiter='*')
  class default
    call err()
  end select
end procedure

module procedure new_QpointCombination_String
  call this%read(input)
end procedure

module procedure conjg_QpointCombination
  output%qpoints_ = conjg(this%qpoints_)
end procedure

module procedure operate_Group_QpointCombination
  output = QpointCombination(qpoint_group*this%qpoints_)
end procedure

module procedure generate_qpoint_combinations
  type :: PowerArray
    type(QpointPower),    allocatable :: powers(:)
    type(FractionVector), allocatable :: wavevectors(:)
  end type
  
  type :: CombinationAndWavevector
    type(QpointCombination) :: combination
    integer                 :: combination_power
    type(FractionVector)    :: wavevector
  end type
  
  logical              :: qpoints_real
  integer, allocatable :: key(:)
  
  type(PowerArray), allocatable :: powers(:)
  
  type(CombinationAndWavevector), allocatable :: combinations(:)
  type(CombinationAndWavevector), allocatable :: new_combinations(:)
  
  integer :: i,j,k,l,ialloc
  
  if (size(qpoints)==0) then
    if (power==0) then
      output = [QpointCombination([QpointPower::])]
    else
      output = [QpointCombination::]
    endif
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
    powers(i)%powers = generate_qpoint_powers(qpoints(key(i)), power)
    if (qpoints_real) then
      powers(i)%wavevectors = qpoints(key(i))%qpoint*powers(i)%powers%power()
    else
      powers(i)%wavevectors = qpoints(key(i))%qpoint            &
                          & * ( powers(i)%powers%power()        &
                          &   - powers(i)%powers%paired_power() )
    endif
  enddo
  
  ! Initialise `combinations` to the empty combination.
  combinations = [CombinationAndWavevector(                    &
     & combination       = QpointCombination([QpointPower::]), &
     & combination_power = 0,                                  &
     & wavevector        = fracvec(zeroes(3))                  )]
  ! Construct all q-point combinations at the given power.
  ! After iteration `i` of the loop, `combinations` contains the q-point
  !    combinations corresponding to the first `i` q-points, with
  !    total_power <= power.
  ! The final iteration only includes combinations with total_power=power.
  do i=1,size(powers)
    ! Allocate `new_combinations` to the maximum size possibly needed.
    if (allocated(new_combinations)) then
      deallocate(new_combinations, stat=ialloc); call err(ialloc)
    endif
    allocate( new_combinations( size(combinations)            &
            &                 * (1+size(powers(i)%powers)) ), &
            & stat=ialloc); call err(ialloc)
    
    ! Loop over the combinations from the previous `i` loop iteration,
    !    and the powers at q-point `i`.
    l = 0
    do j=1,size(combinations)
      l = l+1
      new_combinations(l) = combinations(j)
      
      if (i==size(powers)) then
        if (      new_combinations(l)%combination_power<power            &
           & .or. (       set_default(conserve_momentum,.false.)         &
           &        .and. .not. is_int(new_combinations(l)%wavevector) ) ) then
          l = l-1
        endif
      endif
      
      do k=1,size(powers(i)%powers)
        if ( combinations(j)%combination_power &
         & + powers(i)%powers(k)%total_power() &
         & > power                             ) then
          exit
        endif
        
        l = l+1
        new_combinations(l)%combination%qpoints_ = [ &
             & combinations(j)%combination%qpoints_, &
             & powers(i)%powers(k)                   ]
        new_combinations(l)%combination_power =  &
           &   combinations(j)%combination_power &
           & + powers(i)%powers(k)%total_power()
        new_combinations(l)%wavevector = combinations(j)%wavevector &
                                     & + powers(i)%wavevectors(k)
        if (i==size(powers)) then
          if (      new_combinations(l)%combination_power<power              &
             & .or. (       set_default(conserve_momentum,.false.)           &
             &        .and. .not. is_int(new_combinations(l)%wavevector) ) ) &
             & then
            l = l-1
          endif
        endif
      enddo
    enddo
    
    combinations = new_combinations(:l)
  enddo
  
  output = combinations(l:1:-1)%combination
end procedure
end submodule
