! MODULE linear_algebra
!!----------------!
!! LINEAR_ALGEBRA !
!!----------------!
! USE constants
! USE utils
!
! IMPLICIT NONE
!
! INTERFACE
! SUBROUTINE zheev(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,INFO)
! CHARACTER(1),INTENT(in) :: JOBZ,UPLO
! INTEGER,INTENT(out) :: INFO
! INTEGER,INTENT(in) :: LDA,LWORK,N
! REAL(KIND(1.d0)),INTENT(out) :: W(*)
! REAL(KIND(1.d0)),INTENT(inout) :: RWORK(*)
! COMPLEX(KIND(1.d0)),INTENT(inout) :: A(LDA,*),WORK(*)
! END SUBROUTINE zheev
! END INTERFACE
!
! CONTAINS
!
!
! END MODULE linear_algebra


 MODULE minimum_image
!---------------!
! MINIMUM_IMAGE !
!---------------!
 USE constants
 USE utils
 USE linear_algebra

 IMPLICIT NONE

 CONTAINS 

 SUBROUTINE min_images_brute_force(a,lat_vec,b,nim)
! Compute the minimum image vector(s) b of vector a with respect to the 
! lattice specified by the columns of lat_vec. rec_vec are the reciprocal 
! lattice vectors (w/o 2pi). -b is the vector from a to its closest lattice 
! point.  nim is the number of image vectors.
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: a(3),lat_vec(3,3)
 REAL(dp),INTENT(out) :: b(3,8)
 INTEGER,INTENT(out) :: nim
 REAL(dp) :: rec_vec(3,3),delta1(3),delta2(3),delta3(3),mag_b_sq,dist2,tol_L2
 INTEGER :: n(3),i,j,k
 INTEGER,PARAMETER :: check_shell=3
 REAL(dp),PARAMETER :: tol=1.d-8

 call inv_33(lat_vec,rec_vec)
 rec_vec=transpose(rec_vec)

 tol_L2=tol*dot_product(lat_vec(1,1:3),lat_vec(1,1:3))
 n(1)=floor(dot_product(a(1:3),rec_vec(1,1:3)))
 n(2)=floor(dot_product(a(1:3),rec_vec(2,1:3)))
 n(3)=floor(dot_product(a(1:3),rec_vec(3,1:3)))

 mag_b_sq=-1.d0
 nim=-1

 do i=n(1)-check_shell,n(1)+check_shell+1
  delta1=a-dble(i)*lat_vec(1,1:3)
  do j=n(2)-check_shell,n(2)+check_shell+1
   delta2=delta1-dble(j)*lat_vec(2,1:3)
   do k=n(3)-check_shell,n(3)+check_shell+1
    delta3=delta2-dble(k)*lat_vec(3,1:3)
    dist2=dot_product(delta3,delta3)
    if(abs(dist2-mag_b_sq)<=tol_L2)then
     nim=nim+1
     if(nim>8)call errstop('MIN_IMAGES_BRUTE_FORCE','Need to increase &
      &maxim parameter.')
     b(1:3,nim)=delta3(1:3)
    elseif(dist2<mag_b_sq.or.nim==-1)then
     mag_b_sq=dist2
     nim=1
     b(1:3,1)=delta3(1:3)
    endif
   enddo ! k
  enddo ! j
 enddo ! i
 
 if(nim<=0)call errstop('MIN_IMAGES_BRUTE_FORCE','Bug.')
 
 END SUBROUTINE min_images_brute_force

 END MODULE minimum_image


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

 call inv_33(latt_vecs,rec_vecs)
 rec_vecs=transpose(rec_vecs)

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
 REAL(dp),PARAMETER :: tol=1.d-8
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


 MODULE phonon
!--------!
! PHONON !
!--------!
 USE constants
 USE utils
 USE linear_algebra
 USE symmetry

 IMPLICIT NONE

 CONTAINS

 SUBROUTINE construct_dyn_mat(kpoint,basis,mass,no_cells,no_ims,cell_vecs,&
  &force_consts,dyn_mat)
!--------------------------------------------------------------------------!
! Construct the dynamical matrix at an arbitrary wave vector. We take into !
! account all minimum image primitive cells in the supercell.              !
!--------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: basis,no_cells,no_ims(no_cells)
 REAL(dp),INTENT(in) :: kpoint(3),mass(basis),cell_vecs(3,8,no_cells),&
  &force_consts(basis,3,basis,3,no_cells)
 COMPLEX(dp),INTENT(out) :: dyn_mat(basis,3,basis,3)
 INTEGER :: i_atom,j_atom,i_cart,j_cart,i_cell,i_im
 REAL(dp) :: k_dot_r,prefactor
 COMPLEX(dp) :: phase
 
 dyn_mat=cmplx(0.d0,0.d0,dp)

 do i_atom=1,basis
  do i_cart=1,3
   do j_atom=1,basis
    do j_cart=1,3
     prefactor=1.d0/sqrt(mass(i_atom)*mass(j_atom))
     do i_cell=1,no_cells
      phase=cmplx(0.d0,0.d0,dp)
      do i_im=1,no_ims(i_cell)
       k_dot_r=dot_product(kpoint(1:3),cell_vecs(1:3,i_im,i_cell))
       phase=phase+cmplx(cos(k_dot_r),-sin(k_dot_r),dp)
      enddo ! i_im
      phase=phase/dble(no_ims(i_cell))
      dyn_mat(i_atom,i_cart,j_atom,j_cart)=&
       &dyn_mat(i_atom,i_cart,j_atom,j_cart)+&
       &force_consts(i_atom,i_cart,j_atom,j_cart,i_cell)*phase
     enddo ! i_cell
     dyn_mat(i_atom,i_cart,j_atom,j_cart)=prefactor*&
      &dyn_mat(i_atom,i_cart,j_atom,j_cart)
    enddo ! j_cart
   enddo ! j_atom
  enddo ! i_cart
 enddo ! i_atom

 END SUBROUTINE construct_dyn_mat


 SUBROUTINE calculate_frequencies(basis,dyn_mat,freqs)
!-----------------------------------------------------------------!
! Diagonalise the dynamical matrix and calculate its eigenvalues. !
!-----------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: basis
 COMPLEX(dp),INTENT(in) :: dyn_mat(basis,3,basis,3)
 REAL(dp),INTENT(out) :: freqs(3*basis)
 INTEGER :: i_atom,j_atom,i_cart,j_cart,i_index,j_index,info
 REAL(dp) :: rwork(9*basis-2),minus_freqs_sq(3*basis)
 COMPLEX(dp) :: temp,temp_mat(3*basis,3*basis),work(6*basis-1)

 i_index=0
 do i_atom=1,basis
  do i_cart=1,3
   i_index=i_index+1
   j_index=0
   do j_atom=1,basis
    do j_cart=1,3
     j_index=j_index+1
     temp_mat(i_index,j_index)=dyn_mat(i_atom,i_cart,j_atom,j_cart)
    enddo ! j_cart
   enddo ! j_atom
  enddo ! i_cart
 enddo ! i_atom

 do i_index=1,3*basis
  temp_mat(i_index,i_index)=cmplx(real(temp_mat(i_index,i_index)),0.d0,dp)
  do j_index=i_index+1,3*basis
   temp=0.5d0*(temp_mat(i_index,j_index)+conjg(temp_mat(j_index,i_index)))
   temp_mat(i_index,j_index)=temp
   temp_mat(j_index,i_index)=conjg(temp)
  enddo ! j_index
 enddo ! i_index

 call zheev('N','U',3*basis,temp_mat(1,1),3*basis,minus_freqs_sq(1),work(1),&
  &6*basis-1,rwork(1),info)
 if(info/=0)call errstop('CALCULATE_FREQUENCIES','ZHEEV failed.')

 i_index=3*basis
 do j_index=1,3*basis
  if(minus_freqs_sq(i_index)>=0.d0)then
   freqs(j_index)=-sqrt(minus_freqs_sq(i_index))
  else
   freqs(j_index)=sqrt(-minus_freqs_sq(i_index))
  endif ! minus_freqs_sq
  i_index=i_index-1
 enddo ! j_index

 ENDSUBROUTINE calculate_frequencies


 SUBROUTINE generate_dispersion(rec_vecs,basis,mass,no_cells,no_ims,cell_vecs,&
  &force_consts,no_points,path)
!---------------------!
! GENERATE_DISPERSION !
!---------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: basis,no_cells,no_ims(no_cells),no_points
 REAL(dp),INTENT(in) :: rec_vecs(3,3),mass(basis),cell_vecs(3,8,no_cells),&
  &force_consts(basis,3,basis,3,no_cells),path(3,no_points)
 INTEGER :: ialloc,ierr,i_path,i_cart,path_length,i_point,j_point,i_dof
 REAL(dp) :: k_start(3),k_stop(3),k_diff(3),k_dist,total_k_dist,delta_k,&
  &kpoint(3),omega(3*basis)
 COMPLEX(dp) :: dyn_mat(basis,3,basis,3)

 total_k_dist=0.d0
 do i_point=1,no_points-1
  do i_cart=1,3
   k_start(i_cart)=dot_product(path(1:3,i_point),rec_vecs(1:3,i_cart))
   k_stop(i_cart)=dot_product(path(1:3,i_point+1),rec_vecs(1:3,i_cart))
  enddo ! i_cart
  k_diff=k_stop-k_start
  total_k_dist=total_k_dist+sqrt(dot_product(k_diff,k_diff))
 enddo ! i_point

 delta_k=total_k_dist/1000.d0

 open(unit=14,file='phonon_dispersion_curve.dat',status='replace',iostat=ierr)
 if(ierr/=0)call errstop('GENERATE_DISPERSION','Problem opening &
  &phonon_dispersion_curve.dat file')

 open(unit=15,file='high_symmetry_points.dat',status='replace',iostat=ierr)
 if(ierr/=0)call errstop('GENERATE_DISPERSION','Problem opening &
  &high_symmetry_points.dat file.')

 total_k_dist=0.d0
 do i_point=1,no_points-1
  write(15,*)i_point,total_k_dist
  do i_cart=1,3
   k_start(i_cart)=dot_product(path(1:3,i_point),rec_vecs(1:3,i_cart))
   k_stop(i_cart)=dot_product(path(1:3,i_point+1),rec_vecs(1:3,i_cart))
  enddo ! i_cart
  k_diff=k_stop-k_start
  k_dist=sqrt(dot_product(k_diff,k_diff))
  path_length=int(k_dist/delta_k)
  do i_path=0,path_length-1
   kpoint=k_start+dble(i_path)*k_diff/dble(path_length)
   call construct_dyn_mat(kpoint,basis,mass,no_cells,no_ims,cell_vecs,&
    &force_consts,dyn_mat)
   call calculate_frequencies(basis,dyn_mat,omega)
   write(14,*)total_k_dist,omega
   total_k_dist=total_k_dist+k_dist/dble(path_length)
  enddo ! i_path
 enddo ! i_point
 call construct_dyn_mat(k_stop,basis,mass,no_cells,no_ims,cell_vecs,&
  &force_consts,dyn_mat)
 call calculate_frequencies(basis,dyn_mat,omega)
 write(14,*)total_k_dist,omega
 write(15,*)no_points,total_k_dist

 close(14)
 close(15)

 ENDSUBROUTINE generate_dispersion


 SUBROUTINE generate_dos(rec_vecs,basis,mass,no_cells,no_ims,cell_vecs,&
  &force_consts)
!--------------!
! GENERATE_DOS !
!--------------!
 IMPLICIT NONE
 INTEGER,PARAMETER :: no_bins=1000,no_prelims=10000,no_samples=1000000
 REAL(dp),PARAMETER :: freq_tol=1.d-8,safety_factor=1.1d0
 REAL(dp) :: T
 INTEGER,INTENT(in) :: basis,no_cells,no_ims(no_cells)
 REAL(dp),INTENT(in) :: rec_vecs(3,3),mass(basis),cell_vecs(3,8,no_cells),&
  &force_consts(basis,3,basis,3,no_cells)
 INTEGER :: ialloc,ierr,i_sample,i_cart,i_freq,i_bin
 REAL(dp) :: max_freq,min_freq,frac(3),kpoint(3),freqs(3*basis),bin_width,&
  &rec_bin_width,freq_dos(no_bins),free_energy,omega
 COMPLEX(dp) :: dyn_mat(basis,3,basis,3)
 LOGICAL :: soft_modes

! Initialise the random number generator
 call random_seed()

 ! Read in temperature
 open(1,file='temperature.dat')
 read(1,*)T
 close(1)

 max_freq=-1.d0
 min_freq=huge(1.d0)

 do i_sample=1,no_prelims
  call random_number(frac)
  do i_cart=1,3
   kpoint(i_cart)=dot_product(frac(1:3),rec_vecs(1:3,i_cart))
  enddo ! i_cart
  call construct_dyn_mat(kpoint,basis,mass,no_cells,no_ims,cell_vecs,&
    &force_consts,dyn_mat)
  call calculate_frequencies(basis,dyn_mat,freqs)
  if(freqs(1)<min_freq)min_freq=freqs(1)
  if(freqs(3*basis)>max_freq)max_freq=freqs(3*basis)
 enddo ! i_sample
 soft_modes=(min_freq<-freq_tol)
 if(max_freq<=0.d0)call errstop('GENERATE_DOS','The system is pathologically &
  &unstable.')

 bin_width=safety_factor*max_freq/dble(no_bins)
 rec_bin_width=1.d0/bin_width
 freq_dos=0.d0

 do i_sample=1,no_samples
  call random_number(frac)
  do i_cart=1,3
   kpoint(i_cart)=dot_product(frac(1:3),rec_vecs(1:3,i_cart))
  enddo ! i_cart
  call construct_dyn_mat(kpoint,basis,mass,no_cells,no_ims,cell_vecs,&
    &force_consts,dyn_mat)
  call calculate_frequencies(basis,dyn_mat,freqs)
  soft_modes=(freqs(1)<-freq_tol)
  do i_freq=1,3*basis
   if(freqs(i_freq)>-freq_tol)then
    i_bin=max(1,ceiling(rec_bin_width*freqs(i_freq)))
    if(i_bin>no_bins)call errstop('GENERATE_DOS','Frequency too high to be &
     &binned.')
    freq_dos(i_bin)=freq_dos(i_bin)+1.d0
   endif ! freqs
  enddo ! i_freq
 enddo ! i_sample

 free_energy=0.d0
 do i_bin=1,no_bins
  omega=bin_width*(dble(i_bin)-0.5d0)
  free_energy=freq_dos(i_bin)*harmonic_free_energy(T,omega)/dble(no_samples)+&
   &free_energy
 enddo ! i_bin
 open(unit=16,file='free_energy.dat',status='replace',iostat=ierr)
 if(ierr/=0)call errstop('GENERATE_DOS','Problem opening free_energy.dat file.')
 write(16,*)free_energy
 close(16)

 freq_dos=freq_dos*rec_bin_width/dble(no_samples)

 open(unit=16,file='freq_dos.dat',status='replace',iostat=ierr)
 if(ierr/=0)call errstop('GENERATE_DOS','Problem opening freq_dos.dat file.')
 do i_bin=1,no_bins
  write(16,*)bin_width*(dble(i_bin)-0.5d0),freq_dos(i_bin)
 enddo ! i_bin
 close(16)

 if(soft_modes)write(*,*)'Soft modes present.'

 END SUBROUTINE generate_dos


 REAL(dp) FUNCTION harmonic_free_energy(T,omega)
!----------------------!
! HARMONIC_FREE_ENERGY !
!----------------------!
 IMPLICIT NONE
 REAL(dp),PARAMETER :: tol=1.d-8
 REAL(dp),PARAMETER :: kB_au_per_K=3.16679002948702D-006
 REAL(dp),INTENT(in) :: T,omega
 REAL(dp) :: difference,kT

 if(T<tol)then
  harmonic_free_energy=0.5d0*omega
 else
  kT=kB_au_per_K*T
  difference=1.d0-exp(-omega/kT)
  if(difference>0.d0)then
   harmonic_free_energy=0.5d0*omega+kT*log(difference)
  else
   harmonic_free_energy=-huge(0.d0)
  endif
 endif

 END FUNCTION harmonic_free_energy

 END MODULE phonon

module fourier_interpolation_module
  use constants
  use utils
  implicit none
contains

 SUBROUTINE read_input_files(no_symm_ops,basis,grid,prim_latt_vecs,symm_ops)
!-------------------------!
! Read basic input files. !
!-------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: no_symm_ops
 INTEGER,INTENT(out) :: basis,grid(3)
 REAL(dp),INTENT(out) :: prim_latt_vecs(3,3),symm_ops(4,3,no_symm_ops)
 INTEGER :: ierr,i_symm,i_row

! Read number of atoms in primitive cell (basis)
 open(unit=11,file='equilibrium.dat',status='old',iostat=ierr)
 if(ierr/=0)call errstop('READ_INPUT_FILES','Problem opening basis.dat file.')
 read(11,*,iostat=ierr)basis
 if(ierr/=0)call errstop('READ_INPUT_FILES','Problem reading basis.dat file.')
 close(11)

! Read grid.dat file
 open(unit=11,file='grid.dat',status='old',iostat=ierr)
 if(ierr/=0)call errstop('READ_INPUT_FILES','Problem opening grid.dat file.')
 read(11,*,iostat=ierr)grid(1:3)
 if(ierr/=0)call errstop('READ_INPUT_FILES','Problem reading grid.dat file.')
 close(11)

! Read primitive lattice
 open(unit=11,file='lattice.dat',status='old',iostat=ierr)
 if(ierr/=0)call errstop('READ_INPUT_FILES','Problem opening prim.dat file.')
 read(11,*,iostat=ierr)prim_latt_vecs(1,1:3)
 if(ierr==0)read(11,*,iostat=ierr)prim_latt_vecs(2,1:3)
 if(ierr==0)read(11,*,iostat=ierr)prim_latt_vecs(3,1:3)
 if(ierr/=0)call errstop('READ_INPUT_FILES','Problem reading prim.dat file.')
 close(11)

! Read symmetry.dat file
 open(unit=11,file='symmetry.dat',status='old',iostat=ierr)
 if(ierr/=0)call errstop('READ_INPUT_FILES','Problem opening symmetry.dat &
  &file.')
 read(11,*,iostat=ierr)
 do i_symm=1,no_symm_ops
  do i_row=1,4
   read(11,*,iostat=ierr)symm_ops(i_row,1:3,i_symm)
   if(ierr/=0)call errstop('READ_INPUT_FILES','Problem reading symmetry.dat &
    &file.')
  enddo ! i_row
 enddo ! i_symm
 close(11)

 END SUBROUTINE read_input_files


 SUBROUTINE read_kpoints(no_kpoints,kpoints,multiplicity,kpoint_to_supercell)
!--------------------------------------------------!
! Read input files related to k-points in the IBZ. !
!--------------------------------------------------!
 IMPLICIT NONE
 REAL(dp),PARAMETER :: tol=1.d-8
 INTEGER,INTENT(in) :: no_kpoints
 INTEGER,INTENT(out) :: multiplicity(no_kpoints),&
  &kpoint_to_supercell(no_kpoints)
 REAL(dp),INTENT(out) :: kpoints(3,no_kpoints)
 INTEGER :: ierr,i_point
 REAL(dp) :: kpoint(3)

! Read ibz.dat file
 open(unit=11,file='ibz.dat',status='old',iostat=ierr)
 if(ierr/=0)call errstop('READ_KPOINTS','Problem opening ibz.dat file.')
 do i_point=1,no_kpoints
  read(11,*,iostat=ierr)kpoints(1:3,i_point),multiplicity(i_point)
  if(ierr/=0)call errstop('READ_KPOINTS','Problem reading ibz.dat file.')
 enddo ! i_point
 close(11)

! Read kpoint_to_supercell.dat file
 open(unit=11,file='kpoint_to_supercell.dat',status='old',iostat=ierr)
 if(ierr/=0)call errstop('READ_KPOINTS','Problem opening &
  &kpoint_to_supercell.dat file.')
 do i_point=1,no_kpoints
  read(11,*,iostat=ierr)kpoint(1:3),kpoint_to_supercell(i_point)
  if(ierr/=0)call errstop('READ_KPOINTS','Problem reading &
   &kpoint_to_supercell.dat file.')
  if(any(abs(kpoint(1:3)-kpoints(1:3,i_point))>tol))then
   call errstop('READ_KPOINTS','k-points in ibz.dat file and &
    &kpoints_to_supercell.dat file disagree.')
  endif ! tol
 enddo ! i_point
 close(11)
 
! Check k-points are in expected order
! do i_point=2,no_kpoints
!  if(.not.(kpoint_to_supercell(i_point)==kpoint_to_supercell(i_point-1)).and.&
!   &.not.(kpoint_to_supercell(i_point)==kpoint_to_supercell(i_point-1)+1))then
!   call errstop('READ_KPOINTS','k-points are not in expected order.')
!  endif ! kpoint_to_supercell
! enddo ! i_point

 do i_point=1,no_kpoints
  kpoints(1:3,i_point)=modulo(kpoints(1:3,i_point)+0.5d0+tol,1.d0)-0.5d0-tol
 enddo ! i_points

 END SUBROUTINE read_kpoints


 SUBROUTINE read_dyn_mats(basis,mass,atom_prim_frac,no_kpoints,dyn_mats,&
  &kpoint_to_supercell)
!--------------------------------------------------------!
! Read in dynamical matrices at each k-point in the IBZ. !
!--------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: basis,no_kpoints,kpoint_to_supercell(no_kpoints)
 REAL(dp),INTENT(out) :: mass(basis),atom_prim_frac(3,basis)
 COMPLEX(dp),INTENT(out) :: dyn_mats(basis,3,basis,3,no_kpoints)
 REAL(dp),PARAMETER :: mass_tol=1.d-4,frac_tol=1.d-8
 INTEGER :: ierr,i_atom,i_cart,j_atom,j_cart,atom1,cart1,atom2,cart2,ibz_point,&
  &supercell,atom_map(basis)
 REAL(dp) :: real_part,imag_part,temp_mass,temp_frac(3)
 LOGICAL :: found_atom(basis)

 open(unit=11,file='atoms_in_primitive_cell.1.dat',status='old',iostat=ierr)
 if(ierr/=0)call errstop('READ_DYN_MATS','Problem opening &
  &atoms_in_primitive_cell.1.dat file.')
 do i_atom=1,basis
  read(11,*,iostat=ierr)mass(i_atom),atom_prim_frac(1:3,i_atom)
  if(ierr/=0)call errstop('READ_DYN_MATS','Problem reading &
  &atoms_in_primitive_cell.1.dat file.')
 enddo ! i_atom
 if(any(atom_prim_frac(1:3,1:basis)<-1.d-4.or.&
  &atom_prim_frac(1:3,1:basis)>=1.d0))call errstop('READ_DYN_MATS',&
   &'Fractional atomic coordinates are not in range [0.0,1.0)')
 close(11)

 open(unit=11,file='dyn_mat.1.dat',status='old',iostat=ierr)
 if(ierr/=0)call errstop('READ_DYN_MATS','Problem opening dyn_mat.1.dat file.')
 do i_atom=1,basis
  do i_cart=1,3
   do j_atom=1,basis
    do j_cart=1,3
     read(11,*,iostat=ierr)atom1,cart1,atom2,cart2,real_part,imag_part
     if(ierr/=0)call errstop('READ_DYN_MATS','Problem reading dyn_mat.1.dat &
      &file.')
     if(atom1/=i_atom.or.cart1/=i_cart.or.atom2/=j_atom.or.cart2/=j_cart)call &
      errstop('READ_DYN_MATS','dyn_mat.1.dat file does not seem to be in the &
       &expected order.')
      dyn_mats(atom1,cart1,atom2,cart2,1)=cmplx(real_part,imag_part,dp)
    enddo ! j_cart
   enddo ! j_atom
  enddo ! i_cart
 enddo ! i_atom
 close(11)

 do ibz_point=2,no_kpoints
  supercell=kpoint_to_supercell(ibz_point)
  open(unit=11,file='atoms_in_primitive_cell.'//trim(i2s(supercell))//'.dat',&
   &status='old',iostat=ierr)
  if(ierr/=0)call errstop('READ_DYN_MATS','Problem opening &
   &atoms_in_primitive_cell.'//trim(i2s(supercell))//'.dat file.')
  found_atom(1:basis)=.false.
  atom_map(1:basis)=0
  do i_atom=1,basis
   read(11,*,iostat=ierr)temp_mass,temp_frac(1:3)
   if(ierr/=0)call errstop('READ_DYN_MATS','Problem reading &
    &atoms_in_primitive_cell.'//trim(i2s(supercell))//'.dat file.')
   do j_atom=1,basis
    if(abs(temp_mass-mass(j_atom))<mass_tol)then
     if(all(abs(temp_frac(1:3)-atom_prim_frac(1:3,j_atom))<frac_tol))then
      found_atom(j_atom)=.true.
      atom_map(i_atom)=j_atom
     endif ! frac_tol
    endif ! mass_tol
   enddo ! j_atom
  enddo ! i_atom
  if(.not.any(found_atom(1:basis)))call errstop('READ_DYN_MATS','Unable to &
   &find all atoms in supercell '//trim(i2s(supercell))//'.')
  close(11)
  open(unit=11,file='dyn_mat.'//trim(i2s(ibz_point))//'.dat',status='old',&
   &iostat=ierr)
  if(ierr/=0)call errstop('READ_DYN_MATS','Problem opening dyn_mat.'//&
   &trim(i2s(ibz_point))//'.dat file.')
  do i_atom=1,basis
   do i_cart=1,3
    do j_atom=1,basis
     do j_cart=1,3
      read(11,*,iostat=ierr)atom1,cart1,atom2,cart2,real_part,imag_part
      if(ierr/=0)call errstop('READ_DYN_MATS','Problem reading dyn_mat.'//&
       &trim(i2s(ibz_point))//'.dat file.')
      if(atom1/=i_atom.or.cart1/=i_cart.or.atom2/=j_atom.or.cart2/=j_cart)call &
       errstop('READ_DYN_MATS','dyn_mat.'//trim(i2s(ibz_point))//'.dat file &
        &does not seem to be in the expected order.')
       dyn_mats(atom_map(atom1),cart1,atom_map(atom2),cart2,ibz_point)=&
        &cmplx(real_part,imag_part,dp)
     enddo ! j_cart
    enddo ! j_atom
   enddo ! i_cart
  enddo ! i_atom
  close(11)
 enddo ! ibz_point

 END SUBROUTINE read_dyn_mats


 SUBROUTINE read_path(no_points,path)
!-----------------------------------------------------------------!
! Read in the high symmetry points on the phonon dispersion path. !
!-----------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: no_points
 REAL(dp),INTENT(out) :: path(3,no_points)
 INTEGER :: ierr,i_point

 open(unit=11,file='path.dat',status='old',iostat=ierr)
 if(ierr/=0)call errstop('READ_PATH','Problem opening path.dat file.')
 do i_point=1,no_points
  read(11,*,iostat=ierr)path(1:3,i_point)
  if(ierr/=0)call errstop('READ_PATH','Problem reading path.dat file.')
 enddo ! i_point
 close(11)

 END SUBROUTINE read_path


subroutine fourier_interpolation()
!-----------------------!
! FOURIER_INTERPOLATION !
!-----------------------!
 USE constants
 USE utils
 USE linear_algebra
 USE minimum_image
 USE symmetry
 USE phonon

 IMPLICIT NONE

 REAL(dp),PARAMETER :: tol=1.d-8 

 INTEGER,ALLOCATABLE :: atom_map_symm_backwards(:,:)
 INTEGER,ALLOCATABLE :: atom_map_symm_forwards(:,:)
 INTEGER,ALLOCATABLE :: grid_map_symm_backwards(:,:)
 INTEGER,ALLOCATABLE :: grid_map_symm_forwards(:,:)
 INTEGER,ALLOCATABLE :: grid_to_ibz_map(:)
 INTEGER,ALLOCATABLE :: ibz_to_grid_symm(:)
 INTEGER,ALLOCATABLE :: ibz_to_supercell_map(:)
 INTEGER,ALLOCATABLE :: identity_map(:)
 INTEGER,ALLOCATABLE :: multiplicity(:)
 INTEGER,ALLOCATABLE :: no_im_cells(:)

 REAL(dp),ALLOCATABLE :: atom_pos_cart(:,:)
 REAL(dp),ALLOCATABLE :: atom_pos_frac(:,:)
 REAL(dp),ALLOCATABLE :: cell_pos_cart(:,:)
 REAL(dp),ALLOCATABLE :: min_im_cell_pos(:,:,:)
 REAL(dp),ALLOCATABLE :: freqs(:,:)
 REAL(dp),ALLOCATABLE :: force_consts(:,:,:,:,:)
 REAL(dp),ALLOCATABLE :: grid_points_cart(:,:)
 REAL(dp),ALLOCATABLE :: grid_points_frac(:,:)
 REAL(dp),ALLOCATABLE :: ibz_points_cart(:,:)
 REAL(dp),ALLOCATABLE :: ibz_points_frac(:,:)
 REAL(dp),ALLOCATABLE :: mass(:)
 REAL(dp),ALLOCATABLE :: path(:,:)
 REAL(dp),ALLOCATABLE :: point_symms(:,:,:)
 REAL(dp),ALLOCATABLE :: symm_ops(:,:,:)
 REAL(dp),ALLOCATABLE :: trans_symms(:,:)

 COMPLEX(dp),ALLOCATABLE :: dyn_mats_grid(:,:,:,:,:)
 COMPLEX(dp),ALLOCATABLE :: dyn_mats_ibz(:,:,:,:,:)
 COMPLEX(dp),ALLOCATABLE :: phase(:,:)
 COMPLEX(dp),ALLOCATABLE :: temp_dyn_mat(:,:,:,:)
 COMPLEX(dp),ALLOCATABLE :: dyn_mats_symm(:,:,:,:,:)

 LOGICAL,ALLOCATABLE :: time_reversal(:)

 INTEGER :: basis
 INTEGER :: grid(3)
 INTEGER :: dof_prim
 INTEGER :: no_grid_points
 INTEGER :: no_ibz_points
 INTEGER :: no_kpoints_path
 INTEGER :: no_symm_ops

 INTEGER :: ialloc
 INTEGER :: ierr
 INTEGER :: istat

 INTEGER :: counter
 INTEGER :: i_atom,j_atom
 INTEGER :: i_back
 INTEGER :: i_cart,j_cart,k_cart
 INTEGER :: i_cell
 INTEGER :: i_grid
 INTEGER :: i_point
 INTEGER :: i_symm
 INTEGER :: m1,m2,m3

 REAL(dp) :: identity(3,3)
 REAL(dp) :: k_dot_r
 REAL(dp) :: kpoint(3)
 REAL(dp) :: prefactor
 REAL(dp) :: prim_latt_vecs(3,3)
 REAL(dp) :: prim_rec_vecs(3,3)
 REAL(dp) :: super_latt_vecs(3,3)
 REAL(dp) :: super_rec_vecs(3,3)

! Get total number of symmetry operations and allocate corresponding arrays
 open(unit=10,file='symmetry.dat',status='old',iostat=ierr)
 if(ierr/=0)then
  write(*,*)'Problem opening symmetry.dat file in main program.'
  stop
 endif ! ierr
 read(10,*,iostat=ierr)no_symm_ops
 if(ierr/=0)then
  write(*,*)'Problem reading symmetry.dat file in main program.'
  stop
 endif ! ierr
 close(10)
 
 allocate(symm_ops(4,3,no_symm_ops),stat=ialloc)
 if(ialloc/=0)call erralloc('SYMM_OPS')
 allocate(point_symms(3,3,no_symm_ops),stat=ialloc)
 if(ialloc/=0)call erralloc('POINT_SYMMS')
 allocate(trans_symms(3,no_symm_ops),stat=ialloc)
 if(ialloc/=0)call erralloc('TRANS_SYMMS')

! Read basic input files and allocate corresponding arrays
 call read_input_files(no_symm_ops,basis,grid,prim_latt_vecs,symm_ops)

 dof_prim=3*basis
 no_grid_points=product(grid(1:3))

 do i_symm=1,no_symm_ops
  point_symms(1:3,1:3,i_symm)=symm_ops(1:3,1:3,i_symm)
  do i_cart=1,3
   trans_symms(i_cart,i_symm)=dot_product(symm_ops(4,1:3,i_symm),&
    &prim_latt_vecs(1:3,i_cart))
  enddo ! i_cart
 enddo ! i_symm

 allocate(identity_map(basis),stat=ialloc)
 if(ialloc/=0)call erralloc('IDENTITY_MAP')

 allocate(grid_to_ibz_map(no_grid_points),stat=ialloc)
 if(ialloc/=0)call erralloc('GRID_TO_IBZ_MAP')

 allocate(ibz_to_grid_symm(no_grid_points),stat=ialloc)
 if(ialloc/=0)call erralloc('IBZ_TO_GRID_SYMM')

 allocate(time_reversal(no_grid_points),stat=ialloc)
 if(ialloc/=0)call erralloc('TIME_REVERSAL')

 allocate(grid_map_symm_backwards(no_grid_points,no_symm_ops),stat=ialloc)
 if(ialloc/=0)call erralloc('GRID_MAP_SYMM_BACKWARDS')

 allocate(grid_map_symm_forwards(no_grid_points,no_symm_ops),stat=ialloc)
 if(ialloc/=0)call erralloc('GRID_MAP_SYMM_FORWARDS')

 allocate(grid_points_cart(3,no_grid_points),stat=ialloc)
 if(ialloc/=0)call erralloc('GRID_POINTS_CART')

 allocate(grid_points_frac(3,no_grid_points),stat=ialloc)
 if(ialloc/=0)call erralloc('GRID_POINTS_FRAC')

 allocate(dyn_mats_grid(basis,3,basis,3,no_grid_points),stat=ialloc)
 if(ialloc/=0)call erralloc('DYN_MATS_GRID')

 allocate(temp_dyn_mat(basis,3,basis,3),stat=ialloc)
 if(ialloc/=0)call erralloc('TEMP_DYN_MAT')

 allocate(dyn_mats_symm(basis,3,basis,3,no_grid_points),stat=ialloc)
 if(ialloc/=0)call erralloc('DYN_MATS_SYMM')

 allocate(force_consts(basis,3,basis,3,no_grid_points),stat=ialloc)
 if(ialloc/=0)call erralloc('FORCE_CONSTS')

 identity=0.d0
 identity(1,1)=1.d0
 identity(2,2)=1.d0
 identity(3,3)=1.d0

 do i_atom=1,basis
  identity_map(i_atom)=i_atom
 enddo ! i_atom

 call inv_33(prim_latt_vecs,prim_rec_vecs)
 prim_rec_vecs=transpose(prim_rec_vecs)
 prim_rec_vecs=twopi*prim_rec_vecs

 super_latt_vecs(1,1:3)=dble(grid(1))*prim_latt_vecs(1,1:3)
 super_latt_vecs(2,1:3)=dble(grid(2))*prim_latt_vecs(2,1:3)
 super_latt_vecs(3,1:3)=dble(grid(3))*prim_latt_vecs(3,1:3)

 call inv_33(super_latt_vecs,super_rec_vecs)
 super_rec_vecs=transpose(super_rec_vecs)
 super_rec_vecs=twopi*super_rec_vecs

 i_grid=0
 do m1=0,grid(1)-1
  do m2=0,grid(2)-1
   do m3=0,grid(3)-1
    i_grid=i_grid+1
    if(i_grid>no_grid_points)then
     write(*,*)'Found more k-points than on grid.'
     stop
    endif ! i_grid
    grid_points_frac(1,i_grid)=dble(m1)/dble(grid(1))
    grid_points_frac(2,i_grid)=dble(m2)/dble(grid(2))
    grid_points_frac(3,i_grid)=dble(m3)/dble(grid(3))
    grid_points_frac(1:3,i_grid)=&
     &modulo(0.5d0+grid_points_frac(1:3,i_grid)+tol,1.d0)-0.5d0-tol
   enddo ! m3
  enddo ! m2
 enddo ! m1
 if(i_grid<no_grid_points)then
  write(*,*)'Not found all k-points on grid.'
  stop
 endif ! i_grid

! Convert grid points from fractional to Cartesian coodinates
 do i_grid=1,no_grid_points
  do i_cart=1,3
   grid_points_cart(i_cart,i_grid)=dot_product(grid_points_frac(1:3,i_grid),&
    &prim_rec_vecs(1:3,i_cart))
  enddo ! i_cart
 enddo ! i_grid

! Determine which symmetry operations are in the point group and inverse group
! for each wave vector
 call kpoint_symmetry_maps(prim_rec_vecs,no_grid_points,grid_points_cart,&
  &no_symm_ops,point_symms,grid_map_symm_forwards,grid_map_symm_backwards)

 allocate(atom_pos_frac(3,basis),stat=ialloc)
 if(ialloc/=0)call erralloc('ATOM_POS_FRAC')

 allocate(atom_pos_cart(3,basis),stat=ialloc)
 if(ialloc/=0)call erralloc('ATOM_POS_CART')

 allocate(mass(basis),stat=ialloc)
 if(ialloc/=0)call erralloc('MASS')

 allocate(atom_map_symm_forwards(basis,no_symm_ops),stat=ialloc)
 if(ialloc/=0)call erralloc('ATOM_MAP_SYMM_FORWARDS')

 allocate(atom_map_symm_backwards(basis,no_symm_ops),stat=ialloc)
 if(ialloc/=0)call erralloc('ATOM_MAP_SYMM_BACKWARDS')

 allocate(phase(basis,basis),stat=ialloc)
 if(ialloc/=0)call erralloc('PHASE')

 allocate(cell_pos_cart(3,no_grid_points),stat=ialloc)
 if(ialloc/=0)call erralloc('CELL_POS_CART')

 allocate(no_im_cells(no_grid_points),stat=ialloc)
 if(ialloc/=0)call erralloc('NO_IM_CELLS')

 allocate(min_im_cell_pos(3,8,no_grid_points),stat=ialloc)
 if(ialloc/=0)call erralloc('MIN_IM_CELL_POS')

 i_cell=0
 do m1=0,grid(1)-1
  do m2=0,grid(2)-1
   do m3=0,grid(3)-1
    i_cell=i_cell+1
    if(i_cell>no_grid_points)then
     write(*,*)'Found more primitive cells than in supercell.'
     stop
    endif ! i_cell
    cell_pos_cart(:,i_cell)=dble(m1)*prim_latt_vecs(1,:)+&
     &dble(m2)*prim_latt_vecs(2,:)+dble(m3)*prim_latt_vecs(3,:)
    call min_images_brute_force(cell_pos_cart(:,i_cell),super_latt_vecs,&
     &min_im_cell_pos(:,:,i_cell),no_im_cells(i_cell))
   enddo ! m3
  enddo ! m2
 enddo ! m1
 if(i_cell<no_grid_points)then
  write(*,*)'Not found all primitive cells in supercell.'
  stop
 endif ! i_cell

! Get the number of k-points in the IBZ and allocate corresponding arrays
 call system("echo $(wc -l ibz.dat | awk '{print $1}') > tempfile.dat",istat)
 if(istat/=0)then
  write(*,*)'Problem counting the number of lines in ibz.dat.'
  stop
 endif ! istat
 open(unit=10,file='tempfile.dat',status='old',iostat=ierr)
 if(ierr/=0)then
  write(*,*)'Problem opening tempfile.dat file.'
  stop
 endif ! ierr
 read(10,*,iostat=ierr)no_ibz_points
 if(ierr/=0)then
  write(*,*)'Problem reading tempfile.dat file.'
  stop
 endif ! ierr
 close(10,status='delete')

 allocate(ibz_points_cart(3,no_ibz_points),stat=ialloc)
 if(ialloc/=0)call erralloc('IBZ_POINTS_CART')

 allocate(ibz_points_frac(3,no_ibz_points),stat=ialloc)
 if(ialloc/=0)call erralloc('IBZ_POINTS_FRAC')

 allocate(ibz_to_supercell_map(no_ibz_points),stat=ialloc)
 if(ialloc/=0)call erralloc('IBZ_TO_SUPERCELL_MAP')

 allocate(multiplicity(no_ibz_points),stat=ialloc)
 if(ialloc/=0)call erralloc('MULTIPLICITY')

 allocate(dyn_mats_ibz(basis,3,basis,3,no_ibz_points),stat=ialloc)
 if(ialloc/=0)call erralloc('DYN_MATS_IBZ')

! Read input files related to k-points in the IBZ
 call read_kpoints(no_ibz_points,ibz_points_frac,multiplicity,&
  &ibz_to_supercell_map)

! Convert IBZ points from fractional to Cartesian coodinates
 do i_point=1,no_ibz_points
  do i_cart=1,3
   ibz_points_cart(i_cart,i_point)=dot_product(ibz_points_frac(1:3,i_point),&
    &prim_rec_vecs(1:3,i_cart))
  enddo ! i_cart
 enddo ! i_point

! Read in the dynamical matrix at each k-point in the IBZ
 call read_dyn_mats(basis,mass,atom_pos_frac,no_ibz_points,dyn_mats_ibz,&
  &ibz_to_supercell_map)

 do i_atom=1,basis
  do i_cart=1,3
   atom_pos_cart(i_cart,i_atom)=dot_product(atom_pos_frac(1:3,i_atom),&
    &prim_latt_vecs(1:3,i_cart))
  enddo ! i_cart
 enddo ! i_atom

! Determine the mapping of the atoms in the primitive cell under each symmetry
! operation
  call atom_symmetry_maps(prim_latt_vecs,basis,atom_pos_cart,no_symm_ops,&
  &point_symms,trans_symms,atom_map_symm_forwards,atom_map_symm_backwards)

! Map each k-point on the grid to one in the IBZ
 call match_kpoints(prim_rec_vecs,no_grid_points,grid_points_cart,&
  &no_ibz_points,ibz_points_cart,no_symm_ops,point_symms,grid_to_ibz_map,&
  &ibz_to_grid_symm,time_reversal)

! Determine the dynamical matrix at each point on the k-point grid by applying 
! the appropriate unitary trasformation to the dynamical matrix at a point in 
! the IBZ
 do i_grid=1,no_grid_points
  i_point=grid_to_ibz_map(i_grid)
  kpoint=ibz_points_cart(1:3,i_point)
  i_symm=ibz_to_grid_symm(i_grid)
  call g_matrix_phases(kpoint,point_symms(:,:,i_symm),trans_symms(:,i_symm),&
   &basis,atom_pos_cart,phase)
  call apply_symmetry_to_dyn_mat(basis,atom_map_symm_forwards(:,i_symm),phase,&
   &dyn_mats_ibz(:,:,:,:,i_point),point_symms(:,:,i_symm),&
   &dyn_mats_grid(:,:,:,:,i_grid))
  if(time_reversal(i_grid))then
   dyn_mats_grid(:,:,:,:,i_grid)=&
    &cmplx(real(dyn_mats_grid(:,:,:,:,i_grid)),0.d0,dp)
  endif ! time_reversal
 enddo ! i_grid

! Symmetrize the dynamical matrix at each wave vector with respect to the
! translation operations of the supercell
 dyn_mats_symm=cmplx(0.d0,0.d0,dp)
 do i_grid=1,no_grid_points
  kpoint=grid_points_cart(1:3,i_grid)
  do i_cell=1,no_grid_points
   call g_matrix_phases(kpoint,identity,cell_pos_cart(:,i_cell),basis,&
    &atom_pos_cart,phase)
   call apply_symmetry_to_dyn_mat(basis,identity_map,phase,&
    &dyn_mats_grid(:,:,:,:,i_grid),identity,temp_dyn_mat)
   dyn_mats_symm(:,:,:,:,i_grid)=temp_dyn_mat+dyn_mats_symm(:,:,:,:,i_grid)
  enddo ! i_cell
  dyn_mats_symm(:,:,:,:,i_grid)=1.d0/dble(no_grid_points)*&
   &dyn_mats_symm(:,:,:,:,i_grid)
  if(time_reversal(i_grid))then
   dyn_mats_symm(:,:,:,:,i_grid)=&
    &cmplx(real(dyn_mats_symm(:,:,:,:,i_grid)),0.d0,dp)
  endif ! time_reversal
 enddo ! i_grid
 dyn_mats_grid=dyn_mats_symm

! Symmetrize the dynamical matrix at each wave vector with respect to its
! point and inverse groups
 dyn_mats_symm=cmplx(0.d0,0.d0,dp)
 do i_grid=1,no_grid_points
  counter=0
  do i_symm=1,no_symm_ops
   if(grid_map_symm_backwards(i_grid,i_symm)>0)then
    counter=counter+1
    i_back=grid_map_symm_backwards(i_grid,i_symm)
    kpoint=grid_points_cart(1:3,i_back)
    call g_matrix_phases(kpoint,point_symms(:,:,i_symm),trans_symms(:,i_symm),&
     &basis,atom_pos_cart,phase)
    call apply_symmetry_to_dyn_mat(basis,atom_map_symm_forwards(:,i_symm),&
     &phase,dyn_mats_grid(:,:,:,:,i_back),point_symms(:,:,i_symm),temp_dyn_mat)
    dyn_mats_symm(:,:,:,:,i_grid)=temp_dyn_mat+dyn_mats_symm(:,:,:,:,i_grid)
   endif ! grid_map_symm_backwards
  enddo ! i_symm
  dyn_mats_symm(:,:,:,:,i_grid)=1.d0/dble(counter)*dyn_mats_symm(:,:,:,:,i_grid)
  if(time_reversal(i_grid))then
   dyn_mats_symm(:,:,:,:,i_grid)=&
    &cmplx(real(dyn_mats_symm(:,:,:,:,i_grid)),0.d0,dp)
  endif ! time_reversal
 enddo ! i_grid
 dyn_mats_grid=dyn_mats_symm

! do i_grid=1,no_grid_points
!  open(unit=10,file='dyn_mat_symm.'//trim(i2s(i_grid))//'.dat',status='replace')
!  do i_atom=1,basis
!   do i_cart=1,3
!    do j_atom=1,basis
!     do j_cart=1,3
!      write(10,*)i_atom,i_cart,j_atom,j_cart,&
!       &real(dyn_mats_grid(i_atom,i_cart,j_atom,j_cart,i_grid)),&
!       &aimag(dyn_mats_grid(i_atom,i_cart,j_atom,j_cart,i_grid))
!     enddo ! j_cart
!    enddo ! j_atom
!   enddo ! i_cart
!  enddo ! i_atom
!  close(10)
! enddo ! i_grid

! Construct the matrix of force constants
 force_consts=0.d0
 do i_cell=1,no_grid_points
  do i_atom=1,basis
   do i_cart=1,3
    do j_atom=1,basis
     do j_cart=1,3
      prefactor=sqrt(mass(i_atom)*mass(j_atom))/&
       &dble(no_grid_points)
      do i_grid=1,no_grid_points
       kpoint=grid_points_cart(:,i_grid)
       k_dot_r=dot_product(kpoint(:),cell_pos_cart(:,i_cell))
       force_consts(i_atom,i_cart,j_atom,j_cart,i_cell)=&
        &force_consts(i_atom,i_cart,j_atom,j_cart,i_cell)+&
        &real(dyn_mats_grid(i_atom,i_cart,j_atom,j_cart,i_grid)*&
        &cmplx(cos(k_dot_r),sin(k_dot_r),dp))
      enddo ! i_grid
      force_consts(i_atom,i_cart,j_atom,j_cart,i_cell)=prefactor*&
       &force_consts(i_atom,i_cart,j_atom,j_cart,i_cell)
     enddo ! j_cart
    enddo ! j_atom
   enddo ! i_cart
  enddo ! i_atom
!  open(unit=10,file='force_consts.'//trim(i2s(i_cell))//'.dat',status='replace')
!  do i_atom=1,basis
!   do i_cart=1,3
!    do j_atom=1,basis
!     do j_cart=1,3
!      write(10,*)i_atom,i_cart,j_atom,j_cart,&
!       &force_consts(i_atom,i_cart,j_atom,j_cart,i_cell)
!     enddo ! j_cart
!    enddo ! j_atom
!   enddo ! i_cart
!  enddo ! i_atom
!  close(10)
 enddo ! i_cell

! allocate(freqs(dof_prim,no_grid_points),stat=ialloc)
! if(ialloc/=0)call erralloc('FREQS')

! Get the number of high symmetry points on the dispersion path and allocate 
! corresponding arrays
 call system("echo $(wc -l path.dat | awk '{print $1}') > tempfile.dat",istat)
 if(istat/=0)then
  write(*,*)'Problem counting the number of lines in path.dat file.'
  stop
 endif ! istat
 open(unit=10,file='tempfile.dat',status='old',iostat=ierr)
 if(ierr/=0)then
  write(*,*)'Problem opening tempfile.dat file.'
  stop
 endif ! ierr
 read(10,*,iostat=ierr)no_kpoints_path
 if(ierr/=0)then
  write(*,*)'Problem reading tempfile.dat file.'
  stop
 endif ! ierr
 close(10,status='delete')

 allocate(path(3,no_kpoints_path),stat=ialloc)
 if(ialloc/=0)call erralloc('PATH')

 call read_path(no_kpoints_path,path)

 call generate_dispersion(prim_rec_vecs,basis,mass,no_grid_points,no_im_cells,&
  &min_im_cell_pos,force_consts,no_kpoints_path,path)

 call generate_dos(prim_rec_vecs,basis,mass,no_grid_points,no_im_cells,&
  &min_im_cell_pos,force_consts)
end subroutine
end module
