
 MODULE symmetry
!----------!
! SYMMETRY !
!----------!
 USE constants
 USE utils
 USE linear_algebra

 IMPLICIT NONE

 CONTAINS

 LOGICAL FUNCTION lattice_point(cart,latt_vecs)
!------------------------------------------------------------------------------!
! Given a vector with Cartesian coordinates cart, this function is true if the !
! vector is a point on the lattice specified by latt_vecs.                     !
!------------------------------------------------------------------------------!
 IMPLICIT NONE
 REAL(dp),PARAMETER :: tol=1.d-4
 REAL(dp),INTENT(in) :: cart(3),latt_vecs(3,3)
 INTEGER :: i_cart
 REAL(dp) :: frac(3),rec_vecs(3,3)

 rec_vecs = transpose(invert(latt_vecs))

 do i_cart=1,3
  frac(i_cart)=dot_product(rec_vecs(i_cart,1:3),cart(1:3))
 enddo ! i_cart

 frac(1:3)=modulo(frac(1:3)+tol,1.d0)-tol

 lattice_point=all(abs(frac(1:3))<tol)

 END FUNCTION lattice_point


 SUBROUTINE kpoint_symmetry_maps(rec_vecs,no_points,points_cart,&
  &no_symms,point_symms,forwards,backwards)
!-----------------------------------------------------------------------------!
! Determine the mapping under each symmetry operation for each k-point on the !
! grid.                                                                       !
!-----------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: no_points,no_symms
 REAL(dp),INTENT(in) :: rec_vecs(3,3),points_cart(3,no_points),&
  &point_symms(3,3,no_symms)
 INTEGER,INTENT(out) :: forwards(no_points,no_symms),&
  &backwards(no_points,no_symms)
 INTEGER :: i_grid,i_cart,i_symm,j_grid
 REAL(dp) :: rot_point_cart(3)

 forwards=0
 backwards=0

 do i_symm=1,no_symms
  do i_grid=1,no_points
   do i_cart=1,3
    rot_point_cart(i_cart)=dot_product(point_symms(i_cart,1:3,i_symm),&
     &points_cart(1:3,i_grid))
   enddo ! i_cart
   do j_grid=1,no_points
    if(lattice_point(rot_point_cart-points_cart(1:3,j_grid),rec_vecs))then
     if(forwards(i_grid,i_symm)/=0)then
      call errstop('KPOINT_SYMMETRY_GROUPS','Grid point '//trim(i2s(i_grid))//&
       &' is transformed to more than one grid point by symmetry operation '&
       &//trim(i2s(i_symm))//'.')
     else
      forwards(i_grid,i_symm)=j_grid
     endif ! forwards
     if(backwards(j_grid,i_symm)/=0)then
      call errstop('KPOINT_SYMMETRY_GROUPS','More than one grid point is &
       &transformed to grid point '//trim(i2s(i_grid))//' by symmetry &
       &operation '//trim(i2s(i_symm))//'.')
     else
      backwards(j_grid,i_symm)=i_grid
     endif ! backwards
    endif ! lattice_point
   enddo ! j_grid
  enddo ! i_grid
 enddo ! i_symm

 END SUBROUTINE kpoint_symmetry_maps


 SUBROUTINE atom_symmetry_maps(latt_vecs,basis,atom_cart,no_symms,point_symms,&
  &trans_symms,forwards,backwards)
!---------------------------------------------------------------------!
! Construct the mapping of each atom in the primitive cell under each !
! symmetry operation.                                                 !
!---------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: basis,no_symms
 REAL(dp),INTENT(in) :: latt_vecs(3,3),atom_cart(3,basis),&
  &point_symms(3,3,no_symms),trans_symms(3,no_symms)
 INTEGER,INTENT(out) :: forwards(basis,no_symms),backwards(basis,no_symms)
 INTEGER :: i_symm,i_atom,i_cart,j_atom
 REAL(dp) :: symm_pos(3)
 LOGICAL :: found_atom(basis)

 forwards=0
 backwards=0

 do i_symm=1,no_symms
  found_atom=.false.
  do i_atom=1,basis
   do i_cart=1,3
    symm_pos(i_cart)=dot_product(point_symms(i_cart,1:3,i_symm),&
     &atom_cart(1:3,i_atom))+trans_symms(i_cart,i_symm)
   enddo ! i_cart
   do j_atom=1,basis
    if(lattice_point(symm_pos-atom_cart(1:3,j_atom),latt_vecs))then
     found_atom(i_atom)=.true.
     if(forwards(i_atom,i_symm)/=0)then
      call errstop('ATOM_SYMMETRY_MAPS','Atom '//trim(i2s(i_atom))//' is &
       &transformed to more than one atom by symmetry operation '&
       &//trim(i2s(i_symm))//'.')
     else
      forwards(i_atom,i_symm)=j_atom
     endif ! forwards
     if(backwards(j_atom,i_symm)/=0)then
      call errstop('ATOM_SYMMETRY_MAPS','More than one atom is mapped to atom '&
       &//trim(i2s(j_atom))//' by symmetry operation '//trim(i2s(i_symm))//'.')
     else
      backwards(j_atom,i_symm)=i_atom
     endif ! backwards
    endif ! lattice_point
   enddo ! j_atom
  enddo ! i_atom
  if(any(.not.found_atom))call errstop('ATOM_SYMMETRY_MAPS','Unable to &
   &map all atoms under symmetry operation '//trim(i2s(i_symm))//'.')
 enddo ! i_symm

 END SUBROUTINE atom_symmetry_maps


 SUBROUTINE match_kpoints(rec_latt_vecs,no_grid_points,grid_points_cart,&
  &no_ibz_points,ibz_points_cart,no_symms,point_symms,map_ibz,map_symm,&
  &time_reversal)
!-------------------------------------------------!
! Map each k-point on the grid to one in the IBZ. !
!-------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: no_grid_points,no_ibz_points,no_symms
 REAL(dp),INTENT(in) :: rec_latt_vecs(3,3),grid_points_cart(3,no_grid_points),&
  &ibz_points_cart(3,no_ibz_points),point_symms(3,3,no_symms)
 INTEGER,INTENT(out) :: map_ibz(no_grid_points),map_symm(no_grid_points)
 LOGICAL,INTENT(out) :: time_reversal(no_grid_points)
 INTEGER :: i_grid,i_cart,i_point,i_symm
 REAL(dp) :: rot_ibz_point(3)
 LOGICAL :: found_grid_point(no_grid_points)

 found_grid_point=.false.

 do i_point=1,no_ibz_points
  do i_symm=1,no_symms
   do i_cart=1,3
    rot_ibz_point(i_cart)=dot_product(point_symms(i_cart,1:3,i_symm),&
     &ibz_points_cart(1:3,i_point))
   enddo ! i_cart
   do i_grid=1,no_grid_points
    if(.not.found_grid_point(i_grid))then
     if(lattice_point(rot_ibz_point-grid_points_cart(1:3,i_grid),&
      &rec_latt_vecs))then
      found_grid_point(i_grid)=.true.
      map_ibz(i_grid)=i_point
      map_symm(i_grid)=i_symm
      time_reversal(i_grid)=lattice_point(2.d0*grid_points_cart(1:3,i_grid),&
       &rec_latt_vecs)
     endif ! lattice_point
    endif ! found_grid_point
   enddo ! i_grid
  enddo ! i_symm
 enddo ! i_point

 if(any(.not.found_grid_point))call errstop('MATCH_KPOINTS','Unable to map &
  &all k-points on grid to the IBZ.')

 END SUBROUTINE match_kpoints


 SUBROUTINE g_matrix_phases(kpoint_cart,point_symm,trans_symm,basis,atoms_cart,&
  &phase)
!------------------------------------------!
! Calculate phases used in gamma matrices. !
!------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: basis
 REAL(dp),INTENT(in) :: kpoint_cart(3),point_symm(3,3),trans_symm(3),&
  &atoms_cart(3,basis)
 COMPLEX(dp),INTENT(out) :: phase(basis,basis)
 INTEGER :: i_atom,i_cart,j_atom,j_cart
 REAL(dp) :: rot_kpoint(3),symm_pos(3),arg

 do i_cart=1,3
  rot_kpoint(i_cart)=dot_product(point_symm(i_cart,1:3),kpoint_cart(1:3))
 enddo ! i_cart

 do i_atom=1,basis
  do j_atom=1,basis
   do j_cart=1,3
    symm_pos(j_cart)=dot_product(point_symm(j_cart,1:3),&
     &atoms_cart(1:3,j_atom))+trans_symm(j_cart)
   enddo ! j_cart
   arg=dot_product(rot_kpoint,atoms_cart(1:3,i_atom)-symm_pos)
   phase(i_atom,j_atom)=cmplx(cos(arg),sin(arg),dp)
  enddo ! j_atom
 enddo ! i_atom

 END SUBROUTINE g_matrix_phases


 SUBROUTINE apply_symmetry_to_dyn_mat(basis,forwards,phase,dyn_mat_in,&
  &point_symm,dyn_mat_out)
!-----------------------------------------------------!
! Apply a symmetry operation to the dynamical matrix. !
!-----------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: basis,forwards(basis)
 REAL(dp),INTENT(in) :: point_symm(3,3)
 COMPLEX(dp),INTENT(in) :: phase(basis,basis),dyn_mat_in(basis,3,basis,3)
 COMPLEX(dp),INTENT(out) :: dyn_mat_out(basis,3,basis,3)
 INTEGER :: i_atom,j_atom,symm_i_atom,symm_j_atom
 REAL(dp) :: trans_symm(3,3)
 COMPLEX(dp) :: temp_mat(3,3)

 trans_symm=transpose(point_symm)

 dyn_mat_out=cmplx(0.d0,0.d0,dp)

 do i_atom=1,basis
  symm_i_atom=forwards(i_atom)
  do j_atom=1,basis
  symm_j_atom=forwards(j_atom)
  temp_mat=dyn_mat_in(i_atom,1:3,j_atom,1:3)
  temp_mat=matmul(matmul(point_symm,temp_mat),trans_symm)
  dyn_mat_out(symm_i_atom,1:3,symm_j_atom,1:3)=phase(symm_j_atom,j_atom)*&
   &temp_mat(1:3,1:3)*conjg(phase(symm_i_atom,i_atom))
  enddo ! j_atom
 enddo ! i_atom

! do i_atom=1,basis
!  inv_i_atom=inv(i_atom)
!  do j_atom=1,basis
!  inv_j_atom=inv(j_atom)
!  temp_mat=dyn_mat_in(inv_i_atom,1:3,inv_j_atom,1:3)
!  temp_mat=matmul(matmul(point_symm,temp_mat),trans_symm)
!  dyn_mat_out(i_atom,1:3,j_atom,1:3)=&
!   &phase(i_atom,inv_i_atom)*temp_mat*conjg(phase(j_atom,inv_j_atom))
!  enddo ! j_atom
! enddo ! i_atom

 END SUBROUTINE apply_symmetry_to_dyn_mat

 END MODULE symmetry
