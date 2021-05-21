!> Provides the implementation of the [[QpointCombination(type)]] methods.
submodule (caesar_qpoint_combination_module) &
   & caesar_qpoint_combination_submodule
  use caesar_stars_module
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

module procedure wavevector_QpointCombination
  integer :: i
  
  if (size(this%qpoints_)==0) then
    output = fracvec(zeroes(3))
  else
    output = sum([( this%qpoints_(i)%wavevector(qpoints), &
                  & i=1,                                  &
                  & size(this%qpoints_)                   )])
  endif
end procedure

module procedure complex_monomials_QpointCombination
  type(ComplexMonomial), allocatable :: old(:)
  type(ComplexMonomial), allocatable :: new(:)
  
  type(ComplexMonomial), allocatable :: qpoint_monomials(:)
  
  integer :: i,j,k,l,ialloc
  
  new = [ComplexMonomial((1.0_dp,0.0_dp), [ComplexUnivariate::])]
  do i=1,size(this%qpoints_)
    ! Calculate the monomials at qpoint i.
    qpoint_monomials = this%qpoints_(i)%complex_monomials(modes)
    
    ! Copy the array "new" to "old", and re-allocate "new" to be large enough
    !    to store the monomials for this iteration.
    old = new
    deallocate(new, stat=ialloc); call err(ialloc)
    allocate( new(size(old)*(size(qpoint_monomials))), &
            & stat=ialloc); call err(ialloc)
    
    ! Loop over "old" and "qpoint_monomials",
    !    and construct all monomials which are products of the two.
    l = 0
    do j=1,size(old)
      do k=1,size(qpoint_monomials)
        l = l+1
        new(l) = old(j)*qpoint_monomials(k)
      enddo
    enddo
  enddo
  
  output = new
  
  do i=1,size(output)
    output(i)%coefficient = sqrt(no_permutations(output(i)))
  enddo
end procedure

function no_permutations(input) result(output) 
  type(ComplexMonomial), intent(in) :: input
  real(dp)                          :: output
  
  integer, allocatable :: powers(:)
  
  integer :: i,ialloc
  
  allocate(powers(0), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    powers = [powers, input%power(i)]
    if (input%id(i)/=input%paired_id(i)) then
      powers = [powers, input%paired_power(i)]
    endif
  enddo
  
  output = real_multinomial(sum(powers), powers)
end function

module procedure read_QpointCombination
  type(String),      allocatable :: line(:)
  type(QpointPower), allocatable :: qpoints(:)
  
  integer :: i,ialloc
  
  select type(this); type is(QpointCombination)
    ! Splitting the input by '*' separates the q-points,
    !    but also splits q-point pairs in two.
    line = split_line(input,delimiter='*')
    
    if (size(line)==1) then
      if (line(1)=='()') then
        this = QpointCombination([QpointPower::])
        return
      endif
    endif
    
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
    if (size(this%qpoints_)==0) then
      output = '()'
    else
      output = join(this%qpoints_, delimiter='*')
    endif
  class default
    call err()
  end select
end procedure

module procedure new_QpointCombination_String
  call this%read(input)
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
    enddo
    
    output = .false.
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

module procedure conjg_QpointCombination
  output%qpoints_ = conjg(this%qpoints_)
end procedure

module procedure operate_Group_QpointCombination
  output = QpointCombination(qpoint_group*this%qpoints_)
end procedure
end submodule
