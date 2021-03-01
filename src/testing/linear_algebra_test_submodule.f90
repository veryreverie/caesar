submodule (caesar_linear_algebra_test_module) caesar_linear_algebra_test_submodule
  use caesar_testing_module
contains

module procedure linear_algebra_test
  ! Working variables.
  type(IntVector)  :: iv1,iv2
  type(RealVector) :: rv1,rv2
  type(IntMatrix)  :: im1,im2
  type(RealMatrix) :: rm1,rm2
  type(IntMatrix) :: matrix
  type(IntVector) :: vector
  type(IntVector) :: vec_of_vecs(2)
  
  integer :: i
  
  call print_line('')
  call print_line('Testing linear algebra.')
  
  iv1 = vec([ 2,3,-1 ])
  iv2 = vec([ 1,-2,3 ])
  
  if (iv1/=iv1) then
    call err()
  elseif (iv1==iv2) then
    call err()
  elseif (iv1*iv2/=-7) then
    call err()
  elseif (-iv1*iv2/=7) then
    call err()
  endif
  
  rv1 = vec([ 1.0_dp, 0.5_dp, -1.0_dp ])
  rv2 = vec([ 3.0_dp, -0.3_dp, 1.0_dp ])
  
  if (abs(rv1*rv2-1.85_dp)>1.0e-10_dp) then
    call err()
  elseif (abs(l2_norm(rv1)-sqrt(2.25_dp))>1.0e-10_dp) then
    call err()
  endif
  
  im1 = mat( [ 1,  1,  2, &
             & 2,  1, -3, &
             & 1, -2, -1  ], 3,3)
  
  if (im1*iv1 /= vec([ 3, 10, -3 ])) then
    call err()
  elseif (iv1*im1 /= vec([ 7, 7, -4])) then
    call err()
  elseif (determinant(im1) /= -18) then
    call err()
  endif
  
  im2 = mat( [  2,  2, 4, &
             &  4, -2, 4, &
             & -2,  4, 2  ], 3,3)
  
  if (im2/2 /= mat([  1, 1, 2, &
                   &  2,-1, 2, &
                   & -1, 2, 1  ],3,3)) then
    call err()
  elseif (im1*im2 /= mat([   2,   8,  12, &
                         &  14, -10,   6, &
                         & - 4,   2, - 6  ], 3,3)) then
    im1 = im1*im2
    call err()
  elseif (iv1*im1*iv2 /= -19) then
    call err()
  endif
  
  rm1 = mat( [ 1.0_dp, 0.0_dp, 1.0_dp, &
             & 0.0_dp, 2.0_dp, 0.0_dp, &
             & 0.0_dp, 1.0_dp, 1.0_dp  ], 3,3)
  
  rm2 = rm1 * mat( [ 0.0_dp, 1.0_dp, 0.0_dp, &
                   & 0.5_dp, 0.0_dp, 1.0_dp, &
                   & 1.0_dp, 0.5_dp, 0.0_dp  ], 3,3)
  
  if ( l2_norm( vec([1,0,0])*rm2 &
            & - vec([1.0_dp,1.5_dp,0.0_dp])) > 1.0e-10_dp ) then
    call err()
  elseif ( l2_norm( vec([0,1,0])*rm2 &
                & - vec([1.0_dp,0.0_dp,2.0_dp])) > 1.0e-10_dp ) then
    call err()
  elseif ( l2_norm( vec([0,0,1])*rm2 &
                & - vec([1.5_dp,0.5_dp,1.0_dp])) > 1.0e-10_dp ) then
    call err()
  endif
  
  rv1 = vec([ 1.0_dp, 0.5_dp, -1.0_dp ])
  rv2 = vec([ 3.0_dp, -0.3_dp, 1.0_dp ])
  rm2 = outer_product(rv1,rv2) - mat([  3.0_dp,  1.5_dp , -3.0_dp, &
                                     & -0.3_dp, -0.15_dp,  0.3_dp, &
                                     &  1.0_dp,  0.5_dp , -1.0_dp  ], 3,3)
  rm2 = outer_product(rv1,rv2)
  if ( l2_norm( vec([1,0,0])*rm2 &
            & - vec([3.0_dp,-0.3_dp,1.0_dp])) > 1.0e-10_dp ) then
    call err()
  elseif ( l2_norm( vec([0,1,0])*rm2 &
                & - vec([1.5_dp,-0.15_dp,0.5_dp])) > 1.0e-10_dp ) then
    call err()
  elseif ( l2_norm( vec([0,0,1])*rm2 &
                & - vec([-3.0_dp,0.3_dp,-1.0_dp])) > 1.0e-10_dp ) then
    call err()
  endif
  
  rm2 = rm1*invert(rm1)
  rv2 = vec([1.0_dp,0.0_dp,0.0_dp])
  if (abs(rv2*rm2*rv2-1.0_dp) > 1.0e-10_dp) then
    call err()
  endif
  rv2 = vec([0.0_dp,1.0_dp,0.0_dp])
  if (abs(rv2*rm2*rv2-1.0_dp) > 1.0e-10_dp) then
    call err()
  endif
  rv2 = vec([0.0_dp,0.0_dp,1.0_dp])
  if (abs(rv2*rm2*rv2-1.0_dp) > 1.0e-10_dp) then
    call err()
  endif
  
  matrix = mat([ 1,1,0, &
               & 0,1,0, &
               & 0,0,1],3,3)
  vector = vec([1,1,1])
  vec_of_vecs = [ vec([1,0,0]), vec([0,1,0]) ]
  
  call print_line('Matrix:')
  call print_lines(matrix)
  
  matrix = transpose(matrix)
  call print_line('Transposed(Matrix):')
  call print_lines(matrix)
  
  call print_line('Matrix*Vector:')
  vector = matrix*vector
  
  call print_line('transpose(Matrix)*Vector:')
  vector = transpose(matrix)*vector
  
  call print_line('Matrix*[Vectors]:')
  do i=1,2
    vec_of_vecs(i) = matrix*vec_of_vecs(i)
  enddo
  
  call print_line('transpose(Matrix)*[Vectors]:')
  do i=1,2
    vec_of_vecs(i) = transpose(matrix)*vec_of_vecs(i)
  enddo
  
  call print_line('')
  call print_line('All tests succesful.')
end procedure
end submodule
