submodule (caesar_qpoint_star_product_module) &
   & caesar_qpoint_star_product_submodule
  use caesar_stars_module
contains

module procedure new_QpointStarProduct
  this%subspace_combination = subspace_combination
  this%qpoint_stars         = qpoint_stars
end procedure

module procedure read_QpointStarProduct
  type(SubspaceCombination)     :: subspace_combination
  type(QpointStar), allocatable :: qpoint_stars(:)
  
  select type(this); type is(QpointStarProduct)
    subspace_combination = SubspaceCombination(token(input(1), 7))
    qpoint_stars = QpointStar(split_into_sections( input(2:),          &
                                                 & separating_line='*' ))
    this = QpointStarProduct(subspace_combination, qpoint_stars)
  class default
    call err()
  end select
end procedure

module procedure write_QpointStarProduct
  integer :: i
  select type(this); type is(QpointStarProduct)
    output = [ 'q-point star product in subspace combination '// &
             & this%subspace_combination//' :'                   ]
    do i=1,size(this%qpoint_stars)
      output = [ output,                                   &
               & str(this%qpoint_stars(i))                 ]
      if (i/=size(this%qpoint_stars)) then
        output = [ output,  &
                 & str('*') ]
      endif
    enddo
  class default
    call err()
  end select
end procedure

module procedure new_QpointStarProduct_Strings
  call this%read(input)
end procedure

module procedure new_QpointStarProduct_StringArray
  call this%read(str(input))
end procedure

module procedure generate_qpoint_star_products
  type(QpointStars), allocatable :: power_stars(:)
  
  integer, allocatable :: products(:)
  
  integer :: i,j,k,ialloc
  
  if (size(subspace_combination)==0) then
    output = [QpointStarProduct(subspace_combination, [QpointStar::])]
    return
  endif
  
  allocate( power_stars(size(subspace_combination)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(subspace_combination)
    j = first(subspace_qpoint_stars%subspace_id==subspace_combination%ids(i))
    k = first(    subspace_qpoint_stars(j)%powers%power &
             & == subspace_combination%powers(i),       &
             & default=0                                )
    if (k==0) then
      output = [QpointStarProduct::]
      return
    elseif (size(subspace_qpoint_stars(j)%powers(k)%stars)==0) then
      output = [QpointStarProduct::]
      return
    endif
    power_stars(i) = subspace_qpoint_stars(j)%powers(k)
  enddo
  
  allocate( products(size(subspace_combination)+1), &
          & stat=ialloc); call err(ialloc)
  products(size(subspace_combination)+1) = 1
  do i=size(subspace_combination),1,-1
    products(i) = products(i+1)*size(power_stars(i)%stars)
  enddo
  
  allocate(output(products(1)), stat=ialloc); call err(ialloc)
  output%subspace_combination = subspace_combination
  do i=1,size(output)
    allocate( output(i)%qpoint_stars(size(subspace_combination)), &
            & stat=ialloc); call err(ialloc)
  enddo
  do i=1,size(subspace_combination)
    do j=1,size(output)
      k = modulo((j-1)/products(i+1), size(power_stars(i)%stars))+1
      output(j)%qpoint_stars(i) = power_stars(i)%stars(k)
    enddo
  enddo
end procedure
end submodule
