! Lattice Thermal Energy, Neil Drummond, 12/2004-1/2005.

! All modifications are labelled as: Modified by B. Monserrat

! This program enables the user to calculate the lattice thermal energy and
! free energy of a crystal lattice within the quasiharmonic approximation.
! It also enables the calculation of phonon dispersion curves.  Finally it
! can be used to calculate the longitudinal and transverse sound speeds for a
! monatomic crystal by averaging the gradients of the acoustic branches
! over all directions.

! The program requires an lte.dat file of the form:

!    Primitive lattice vectors (rows, in a.u.)
!    0.d0 3.25 3.25
!    3.25 0.d0 3.25
!    3.25 3.25 0.d0
!    Supercell lattice vectors (rows, in a.u.)
!    0.d0 6.5 6.5
!    6.5 0.d0 6.5
!    6.5 6.5 0.d0
!    Number of atoms in supercell
!    8
!    Species ; mass (a.u.) ; position of atom in supercell (in terms of SC LVs)
!    Ne  36784.8   0.0  0.0  0.0
!    Ne  36784.8   0.0  0.0  0.5
!      (etc)
!    Ne  36784.8   0.5  0.5  0.5
!    Number of point-symmetry operations
!    48
!    Rotation matrices (3 rows) and translations (1 row, in terms of SC LVs)
!    1.0   0.0   0.0
!    0.0   1.0   0.0
!    0.0   0.0   1.0
!    0.0 0.0 0.0
!      (etc)
!    0.0   0.0   1.0
!    0.0  -1.0   0.0
!   -1.0   0.0   0.0
!    0.0 0.0 0.0
!   Number of force constants supplied
!   24
!   Atom 1 ; Cartes'n direction ; Atom 2 ; C. direction ; force constant (a.u.)
!   1 1 1 1 -0.05
!   1 1 1 2 0.0
!     (etc)
!   1 1 8 3 0.0
!   Program fn: (1) calc thermal en.; (2) calc disp; (3) calc speed of snd
!   1
!   Temperature (K)
!   4.2
!   Number of lines in k space to plot
!   4
!   Points in journey through k space (in terms of 2pi * rec. prim. LVs)
!   0.0 0.0 0.0
!   0.5 0.5 0.0
!   1.0 0.5 0.5
!   0.0 0.0 0.0
!   0.5 0.5 0.5

! Only a limited number of force constants need to be supplied: the program
! will use translational symmetry and the point-symmetry matrices supplied
! in order to complete as many elements of the matrix of force constants as
! possible.

! The point symmetries to be supplied are for the unperturbed supercell.

! The symmetry of the matrix of force constants will be enforced, as will
! Newton's third law.  Inversion symmetry will be enforced for monatomic
! crystals.

! If a dispersion curve is to be plotted then a dummy temperature should be
! supplied, and if the thermal energies are to be calculated then dummy
! dispersion-curve data should be supplied.  Likewise, dummy data should
! be provided if the speed of sound is to be calculated.

! Can also calculate the frequencies and patterns of atomic displacement at the
! supercell G vectors by setting program function to 4.  This is useful if you
! are interested in following a mode.

! CHANGES TO CODE
! ===============
! 20/07/05 Added calculation of speed of sound for monatomic crystals.
!          (Needs eigenvectors of dynamical matrix.)
! 21/07/05 Bugfix: take min image of all prim-cell vectors w.r.t supercell.
! 12/11/05 Fix for nonzero translations: r->b+Rr instead of r->b+R(r-b)
! 30/01/08 Divided freq_dos into separate sets, to allow error estimates.
! 02/02/08 Bugfix: include multiple images where necessary.
! 29/02/08 Applied fix r->b=Rr instead of r->b+R(r-b) everywhere.
! 17/04/09 Eliminated binning of imaginary frequencies, to avoid large free
!          energy due to lowest-freq bin at T>0.  Tidied, introduced i2s, etc.
!          Inserted better min-image routine.
! 08/05/09 Sped up initialisation by eliminating min_image calls.
! 13/05/09 Added evaluation of frequencies on supercell G vectors.
! 18/05/09 Fixed min-image bug introduced on 17/04/09.  Removed more min-images.
! 20/04/10 Enabled calculation of pattern of atomic displacement at each G.
! 17/06/11 Introduced BLAS & LAPACK.  Fixed bug in randomisation of theta for
!          speed-of-sound calculation.

module lte_module
  use min_images,     only : is_lat_point, min_images_brute_force, maxim
  use constants,      only : dp, third, twopi, kB_au_per_K
  use utils,          only : i2s, errstop
  use file_io,        only : open_read_file, open_write_file
  use linear_algebra, only : determinant33, inv_33
  implicit none
  
  private
  public :: lte ! The main function.
contains

! ----------------------------------------------------------------------
! The reciprocal lattice vectors of the primitive cell and the supercell
! are calculated here.
! NOTE THAT RECIPROCAL LATTICE VECTORS do NOT CONTAIN THE FACTOR OF 2*pi.
! ----------------------------------------------------------------------
subroutine setup_geometry(tol,prim_lat_vec,sc_lat_vec,                 &
   & length_scale,no_prim_cells,prim_rec_vec,sc_rec_vec,small_k_scale, &
   & vol_scale)
  use linear_algebra
  implicit none
  
  real(dp), intent(in) :: tol
  real(dp), intent(in) :: prim_lat_vec(3,3)
  real(dp), intent(in) :: sc_lat_vec(3,3)
  
  real(dp), intent(out) :: length_scale
  integer,  intent(out) :: no_prim_cells
  real(dp), intent(out) :: prim_rec_vec(3,3)
  real(dp), intent(out) :: sc_rec_vec(3,3)
  real(dp), intent(out) :: small_k_scale
  real(dp), intent(out) :: vol_scale
  
  real(dp) :: prim_cell_volume
  real(dp) :: sc_cell_volume

  ! Define the "length scale" to be the average of the lengths of the
  ! primitive lattice vectors and the reciprocal length scale to be the
  ! inverse of this.
  length_scale=(SQRT(DOT_PRODUCT(prim_lat_vec(1:3,1), &
    &prim_lat_vec(1:3,1)))+SQRT(DOT_PRODUCT(prim_lat_vec(1:3,2), &
    &prim_lat_vec(1:3,2)))+SQRT(DOT_PRODUCT(prim_lat_vec(1:3,3), &
    &prim_lat_vec(1:3,3))))*third
  vol_scale=length_scale**3

  ! Primitive-cell reciprocal lattice vectors and volume.
  prim_cell_volume=ABS(determinant33(prim_lat_vec))
  if(prim_cell_volume<tol*vol_scale)call errstop('SETUP_GEOMETRY', &
    &'Primitive lattice vectors should be linearly independent.  Please &
    &check your lattice vectors.')
  prim_rec_vec = TRANSPOSE(inv_33(prim_lat_vec))

  ! Supercell reciprocal lattice vectors and volume.
  sc_cell_volume=ABS(determinant33(sc_lat_vec))
  if(sc_cell_volume<tol*vol_scale)call errstop('SETUP_GEOMETRY', &
    &'Supercell lattice vectors should be linearly independent.  Please &
    &check your lattice vectors.')
  sc_rec_vec = transpose(inv_33(sc_lat_vec))

  ! "Small" distance in k space, determined by size of supercell.
  ! Factor of 2pi included.
  small_k_scale=twopi*sc_cell_volume**(-third)

  ! Number of unit cells.
  no_prim_cells=NINT(sc_cell_volume/prim_cell_volume)
  if(ABS(DBLE(no_prim_cells)*prim_cell_volume-sc_cell_volume) &
    &>tol*vol_scale)call errstop('SETUP_GEOMETRY','Supercell volume should &
    &be an integer multiple of primitive-cell volume.  Please check your &
    &lattice vectors.')
end subroutine

! ----------------------------------------------------------------------
! This function returns the number of the atom at rvec (up to translations
! through a supercell lattice vector.
! ----------------------------------------------------------------------
integer function atom_at_pos(rvec,no_atoms_in_sc,atom_pos,sc_rec_vec)
  implicit none
  
  real(dp), intent(in) :: rvec(3)
  integer,  intent(in) :: no_atoms_in_sc
  real(dp), intent(in) :: atom_pos(:,:)
  real(dp), intent(in) :: sc_rec_vec(3,3)
  
  integer :: atom1
  
  do atom1=1,no_atoms_in_sc
    !write(*,*)is_lat_point(atom_pos(1:3,atom1)-rvec(1:3),sc_rec_vec),atom_pos(1:3,atom1),rvec(1:3)
    if(is_lat_point(atom_pos(1:3,atom1)-rvec(1:3),sc_rec_vec))then
      atom_at_pos=atom1
      return
    endif ! Separation of rvec and atom posn. Is a sc lattice vector.
  enddo ! atom1
  atom_at_pos=0
end function

! ----------------------------------------------------------------------
! Read in the data in lte.dat and allocate arrays, etc.
! ----------------------------------------------------------------------
subroutine read_lte(tol,lte_filename,prim_rec_vec,sc_rec_vec,fc_scale,      &
   & no_atoms_in_prim,no_DoF_prim,no_prim_cells,length_scale,small_k_scale, &
   & vol_scale,                                                             &
   & prim_lat_vec,sc_lat_vec,no_atoms_in_sc,species,mass,atom_pos,          &
   & force_const,defined,atom,rotation,offset,disp_kpoints,no_kspace_lines, &
   & no_point_symms,prog_function,temperature)
  implicit none
  
  real(dp),     intent(in) :: tol
  
  character(*), intent(in) :: lte_filename
  
  real(dp),                  intent(out) :: prim_rec_vec(3,3)
  real(dp),                  intent(out) :: sc_rec_vec(3,3)
  real(dp),                  intent(out) :: fc_scale
  integer,                   intent(out) :: no_atoms_in_prim
  integer,                   intent(out) :: no_DoF_prim
  integer,                   intent(out) :: no_prim_cells
  real(dp),                  intent(out) :: length_scale
  real(dp),                  intent(out) :: small_k_scale
  real(dp),                  intent(out) :: vol_scale
  real(dp),                  intent(out) :: prim_lat_vec(3,3)
  real(dp),                  intent(out) :: sc_lat_vec(3,3)
  integer,                   intent(out) :: no_atoms_in_sc
  character(2), allocatable, intent(out) :: species(:)
  real(dp),     allocatable, intent(out) :: mass(:)
  real(dp),     allocatable, intent(out) :: atom_pos(:,:)
  real(dp),     allocatable, intent(out) :: force_const(:,:,:,:)
  logical,      allocatable, intent(out) :: defined(:,:,:,:)
  integer,      allocatable, intent(out) :: atom(:,:)
  real(dp),     allocatable, intent(out) :: rotation(:,:,:)
  real(dp),     allocatable, intent(out) :: offset(:,:)
  real(dp),     allocatable, intent(out) :: disp_kpoints(:,:)
  integer,                   intent(out) :: no_kspace_lines
  integer,                   intent(out) :: no_point_symms
  integer,                   intent(out) :: prog_function
  real(dp),                  intent(out) :: temperature
  
  integer :: ierr,n,i,j,no_force_c_supplied,atom1,atom2,dir1,dir2,ialloc
  real(dp) :: fc,check_matrix

  open(unit=8,file=lte_filename,status='old',iostat=ierr)
  if(ierr/=0)call errstop('READ_LTE','Problem opening lte.dat.')

  ! Primitive lattice vectors
  READ(8,*,err=20,end=20)
  READ(8,*,err=20,end=20)prim_lat_vec(1:3,1)
  READ(8,*,err=20,end=20)prim_lat_vec(1:3,2)
  READ(8,*,err=20,end=20)prim_lat_vec(1:3,3)
  write(*,*)'Primitive lattice vectors (Cartesian components in rows, a.u.):'
  write(*,*)prim_lat_vec(1:3,1)
  write(*,*)prim_lat_vec(1:3,2)
  write(*,*)prim_lat_vec(1:3,3)
  write(*,*)

  ! Supercell lattice vectors
  READ(8,*,err=20,end=20)
  READ(8,*,err=20,end=20)sc_lat_vec(1:3,1)
  READ(8,*,err=20,end=20)sc_lat_vec(1:3,2)
  READ(8,*,err=20,end=20)sc_lat_vec(1:3,3)
  write(*,*)'Supercell lattice vectors (Cartesian components in rows, a.u.):'
  write(*,*)sc_lat_vec(1:3,1)
  write(*,*)sc_lat_vec(1:3,2)
  write(*,*)sc_lat_vec(1:3,3)
  write(*,*)

  ! Construct reciprocal lattice vectors, etc.
  call setup_geometry(tol,prim_lat_vec,sc_lat_vec,                       &
     & length_scale,no_prim_cells,prim_rec_vec,sc_rec_vec,small_k_scale, &
     & vol_scale)
  write(*,*)'Number of primitive unit cells     : '//TRIM(i2s(no_prim_cells))
  write(*,*)
  write(*,*)'Prim. rec. latt. vectors (Cart. cmpnts &
    &in rows, factor of 2pi inc., a.u.):'
  write(*,*)prim_rec_vec(1:3,1)*twopi
  write(*,*)prim_rec_vec(1:3,2)*twopi
  write(*,*)prim_rec_vec(1:3,3)*twopi
  write(*,*)

  write(*,*)'Supercell rec. latt. vectors(Cart. cmpnts &
   &in rows, factor of 2pi inc., a.u.):'
  write(*,*)sc_rec_vec(1:3,1)*twopi
  write(*,*)sc_rec_vec(1:3,2)*twopi
  write(*,*)sc_rec_vec(1:3,3)*twopi
  write(*,*)

  ! Number of atoms.
  READ(8,*,err=20,end=20)
  READ(8,*,err=20,end=20)no_atoms_in_sc
  write(*,*)'Number of atoms in supercell       : '//TRIM(i2s(no_atoms_in_sc))
  if(no_atoms_in_sc<no_prim_cells)call errstop('READ_LTE', &
    &'Need more atoms in the supercell!')
  if(MOD(no_atoms_in_sc,no_prim_cells)/=0)call errstop('READ_LTE', &
    &'Number of atoms in supercell should be a multiple of the number of &
    &primitive cells.')
  no_atoms_in_prim=no_atoms_in_sc/no_prim_cells
  allocate(species(no_atoms_in_sc),                   &
    & mass(no_atoms_in_sc),                           &
    & atom_pos(3,no_atoms_in_sc),                     &
    & force_const(no_atoms_in_sc,3,no_atoms_in_sc,3), &
    & defined(no_atoms_in_sc,3,no_atoms_in_sc,3),     &
    & atom(no_prim_cells,no_atoms_in_prim),           &
    & stat=ialloc)
  if(ialloc/=0)call errstop('READ_LTE','Allocation error: species, etc.')
  defined(:,:,:,:)=.FALSE.
  force_const(:,:,:,:)=0.d0
  no_DoF_prim=3*no_atoms_in_prim
  write(*,*)

  ! Read in atom positions, species and masses.
  ! Convert atom position from fractional coordinates (in file) to
  ! Cartesian coordinates (used in program).  Translate the atom
  ! coordinates into the supercell at the origin.
  READ(8,*,err=20,end=20)
  write(*,*)'Species ; Mass (a.u.) ; Position (Cartesian coordinates, a.u.)'
  do i=1,no_atoms_in_sc
    READ(8,*,err=20,end=20)species(i),mass(i),atom_pos(1:3,i)
    atom_pos(1:3,i)=atom_pos(1,i)*sc_lat_vec(1:3,1)+atom_pos(2,i) &
      &*sc_lat_vec(1:3,2)+atom_pos(3,i)*sc_lat_vec(1:3,3)
    write(*,'(" ",a,"  ",f14.6," ",3("  ",f14.6))')species(i),mass(i), &
      &atom_pos(1:3,i)
    if(mass(i)<=0.d0)call errstop('READ_LTE','Mass should be positive.')
  enddo ! i
  ! Check atoms aren't on top of each other.
  do i=1,no_atoms_in_sc-1
    do j=i+1,no_atoms_in_sc
      if(is_lat_point(atom_pos(1:3,j)-atom_pos(1:3,i),sc_rec_vec))call &
        &errstop('READ_LTE','Atoms '//TRIM(i2s(i))//' and '//TRIM(i2s(j)) &
        &//' appear to be on top of one another.')
    enddo ! j
  enddo ! i
  write(*,*)

  ! Read in point-symmetry operations.
  READ(8,*,err=20,end=20)
  READ(8,*,err=20,end=20)no_point_symms
  write(*,*)'Number of point symmetries         : '//TRIM(i2s(no_point_symms))
  if(no_point_symms<1)call errstop('READ_LTE','At least one point-symmetry &
    &rotation matrix (identity) must be supplied.')
  allocate(rotation(3,3,no_point_symms),offset(3,no_point_symms),stat=ialloc)
  if(ialloc/=0)call errstop('READ_LTE','Allocation error: rotation, etc.')
  READ(8,*,err=20,end=20)
  do n=1,no_point_symms
    READ(8,*,err=20,end=20)rotation(1:3,1,n)
    READ(8,*,err=20,end=20)rotation(1:3,2,n)
    READ(8,*,err=20,end=20)rotation(1:3,3,n)
    do i=1,3
      do j=i,3
        check_matrix=DOT_PRODUCT(rotation(1:3,i,n),rotation(1:3,j,n))
        if((i==j.AND.ABS(check_matrix-1.d0)>tol) &
          &.OR.(i/=j.AND.ABS(check_matrix)>tol))call &
          &errstop('READ_LTE','Rotation matrix '//TRIM(i2s(n)) &
          &//' is not orthogonal!')
      enddo ! j
    enddo ! i
    READ(8,*,err=20,end=20)offset(1:3,n)
    ! Convert translation to Cartesians.
    offset(1:3,n)=offset(1,n)*sc_lat_vec(1:3,1) &
      &+offset(2,n)*sc_lat_vec(1:3,2)+offset(3,n)*sc_lat_vec(1:3,3)
  enddo ! n
  write(*,*)'Have read in rotation matrices and translation vectors.'
  write(*,*)

  ! Read in force constants supplied.
  READ(8,*,err=20,end=20)
  READ(8,*,err=20,end=20)no_force_c_supplied
  write(*,*)'Number of force constants supplied : ' &
    &//TRIM(i2s(no_force_c_supplied))
  if(no_force_c_supplied<=0)call errstop('READ_LTE', &
    &'Need to supply more force data!')
  READ(8,*,err=20,end=20)
  fc_scale=0.d0
  do i=1,no_force_c_supplied
    READ(8,*,err=20,end=20)atom1,dir1,atom2,dir2,fc
    fc_scale=fc_scale+ABS(fc)
    call trans_symm(atom1,dir1,atom2,dir2,fc,tol,atom_pos,no_atoms_in_sc, &
       & no_prim_cells,prim_rec_vec,sc_rec_vec,mass,defined,force_const)
  enddo ! i
  fc_scale=fc_scale/DBLE(no_force_c_supplied)
  write(*,*)fc_scale
  write(*,*)'Have read in the force-constant data and applied &
    &translational symmetry.'
  write(*,*)

  READ(8,*,err=20,end=20)
  READ(8,*,err=20,end=20)prog_function
  if(prog_function==1)then
    write(*,*)'The mean thermal energy and the free energy will &
      &be calculated.'
  elseif(prog_function==2)then
    write(*,*)'A dispersion curve will be calculated.'
  elseif(prog_function==3)then
    write(*,*)'The speed of sound will be calculated.'
  elseif(prog_function==4)then
    write(*,*)'The phonon frequencies at the supercell G vectors will be &
      &calculated.'
  else
    call errstop('READ_LTE','Program function must be either 1, 2, 3 or 4.')
  endif ! prog_function

  READ(8,*,err=20,end=20)
  READ(8,*,err=20,end=20)temperature
  if(prog_function==1)then
    write(*,*)'Temperature (K)                    :',temperature
    if(temperature<0.d0)call errstop('READ_LTE', &
      &'Temperature should be non-negative.')
    if(temperature<=0.d0)write(*,*)'(i.e. the zero-point energy is to be &
      &calculated.)'
    write(*,*)
  endif ! LTE to be calculated.

  READ(8,*,err=20,end=20)
  READ(8,*,err=20,end=20)no_kspace_lines
  if(prog_function==2)then
    write(*,*)'Number of lines in k-space to plot     : ' &
      &//TRIM(i2s(no_kspace_lines))
    if(no_kspace_lines<1)call errstop('READ_LTE', &
      &'Need to supply more lines in k-space!')
    allocate(disp_kpoints(3,0:no_kspace_lines),stat=ialloc)
    if(ialloc/=0)call errstop('READ_LTE','Allocation error: disp_kpoints.')
  endif
  READ(8,*,err=20,end=20)
  if(prog_function==2)write(*,*)'Points along walk in reciprocal space &
    &(Cartesian components in a.u.):'
  do i=0,no_kspace_lines
    if(prog_function==2)then
      READ(8,*,err=20,end=20)disp_kpoints(1:3,i)
      disp_kpoints(1:3,i)=twopi*(disp_kpoints(1,i)*prim_rec_vec(1:3,1) &
        &+disp_kpoints(2,i)*prim_rec_vec(1:3,2) &
        &+disp_kpoints(3,i)*prim_rec_vec(1:3,3))
      write(*,'(3(" ",f16.8))')disp_kpoints(1:3,i)
    else
      READ(8,*,err=20,end=20)
    endif ! prog_function=2
  enddo ! i
  if(prog_function==2)write(*,*)'Have read in points for dispersion curve.'

  close(8)

  RETURN

  ! Error reading file...
20  call errstop('READ_LTE','Problem reading lte.dat.  Please check the format &
    &of the file.')
end subroutine

! ----------------------------------------------------------------------
! Take the force constant fc supplied for (atom1,dir1,atom2,dir2),
! place it in all translationally equivalent elements of the matrix of
! force constants.
! ----------------------------------------------------------------------
subroutine trans_symm(atom1,dir1,atom2,dir2,fc,tol,atom_pos,no_atoms_in_sc, &
   & no_prim_cells,prim_rec_vec,sc_rec_vec,mass,defined,force_const)
  implicit none
  
  real(dp), intent(in) :: tol
  integer , intent(in) :: atom1,dir1,atom2,dir2
  real(dp), intent(in) :: fc
  real(dp), intent(in) :: atom_pos(:,:)
  integer,  intent(in) :: no_atoms_in_sc
  integer,  intent(in) :: no_prim_cells
  real(dp), intent(in) :: prim_rec_vec(3,3)
  real(dp), intent(in) :: sc_rec_vec(3,3)
  real(dp), intent(in) :: mass(:)
  
  logical,  intent(inout) :: defined(:,:,:,:)
  real(dp), intent(inout) :: force_const(:,:,:,:)
  
  integer :: atom1p,atom2p,no_translations
  real(dp) :: pos_atom2p(3),relpos_atom2_atom1(3)
  
  relpos_atom2_atom1(1:3)=atom_pos(1:3,atom2)-atom_pos(1:3,atom1)
  no_translations=0
  do atom1p=1,no_atoms_in_sc
    if(is_lat_point(atom_pos(1:3,atom1p)-atom_pos(1:3,atom1), &
      &prim_rec_vec))then
      ! atom1p and atom1 are equivalent under trans. symm.
      if(ABS(mass(atom1p)-mass(atom1))>tol*mass(atom1))call &
        &errstop('TRANS_SYMM','Atoms '//TRIM(i2s(atom1))//' and ' &
        &//TRIM(i2s(atom1p))//' are equivalent by translational symmetry, &
        &but they have different masses.')
      pos_atom2p(1:3)=atom_pos(1:3,atom1p)+relpos_atom2_atom1(1:3)
      atom2p=atom_at_pos(pos_atom2p,no_atoms_in_sc,atom_pos,sc_rec_vec)
      if(atom2p<=0)call errstop('TRANS_SYMM','Please check that your atom &
        &coordinates satisfy the translational symmetry they should have.')
      ! atom2p and atom2 are related to each other by the same translation
      ! as atom1p and atom1.
      if(ABS(mass(atom2p)-mass(atom2))>tol*mass(atom2))call &
        &errstop('TRANS_SYMM','Atoms '//TRIM(i2s(atom2))//' and ' &
        &//TRIM(i2s(atom2p))//' are equivalent by translational &
        &symmetry, but they have different masses.')
      if(defined(atom1p,dir1,atom2p,dir2))then
        write(*,*)atom1p,dir1,atom2p,dir2
        call errstop('TRANS_SYMM', &
        &'Please check your atom coordinates -- you may have a repeated &
        &atom.  Alternatively, you may have conflicting force data for &
        &atoms that are equivalent under translational symmetry.')
      endif
      force_const(atom1p,dir1,atom2p,dir2)=fc
      defined(atom1p,dir1,atom2p,dir2)=.TRUE.
      no_translations=no_translations+1
    endif ! atom1p equiv to atom1
  enddo ! atom1p
  if(no_translations/=no_prim_cells)call errstop('TRANS_SYMM', &
    &'Number of translationally equivalent atoms found differs from the &
    &number of primitive cells.  Please check the atom coordinates!')
end subroutine

! ----------------------------------------------------------------------
! Apply all of the point symmetry operations in turn in order to complete
! the matrix of force constants.  Give a pair of atoms atom1 and atom2, we
! find the corresponding pair of atoms atom1p and atom2p after the 
! symmetry translation and rotation have been applied.  The matrix of
! force constants Phi then transforms as Phi' = R Phi R^T, where R is the
! rotation matrix and Phi and Phi' are the matrices of force constants
! between atom1 & atom2 and atom1p & atom2p, respectively.  Some elements
! of Phi may be unknown, but so long as they are multiplied by zero we can
! still work out the corresponding elements of Phi'.
! ----------------------------------------------------------------------
subroutine point_symm(tol,no_atoms_in_sc,no_point_symms,rotation, &
   & atom_pos,offset,sc_rec_vec,mass,defined,force_const)
  use linear_algebra, only : ddot
  implicit none
  
  real(dp), intent(in) :: tol
  integer,  intent(in) :: no_atoms_in_sc
  integer,  intent(in) :: no_point_symms
  real(dp), intent(in) :: rotation(:,:,:)
  real(dp), intent(in) :: atom_pos(:,:)
  real(dp), intent(in) :: offset(:,:)
  real(dp), intent(in) :: sc_rec_vec(3,3)
  real(dp), intent(in) :: mass(:)
  
  logical,  intent(inout) :: defined(:,:,:,:)
  real(dp), intent(inout) :: force_const(:,:,:,:)
  
  integer :: atom1,atom1p,atom2,atom2p,i,j,n,ip,jp,&
    &weight(no_atoms_in_sc,3,no_atoms_in_sc,3)
  real(dp) :: fc,product,pos_atom1p(3),pos_atom2p(3)
  logical :: well_defined

  ! We average over all the force constants that ought, by the full
  ! symmetry of the supercell, to be identical.  The weight array is
  ! used to perform this average.
  do j=1,3
    do atom2=1,no_atoms_in_sc
      do i=1,3
        do atom1=1,no_atoms_in_sc
          if(defined(atom1,i,atom2,j))then
            weight(atom1,i,atom2,j)=1
          else
            weight(atom1,i,atom2,j)=0
          endif ! defined
  !        write(*,*)atom1,i,atom2,j,weight(atom1,i,atom2,j)
        enddo ! atom1
      enddo ! i
    enddo ! atom2
  enddo ! j

  do n=1,no_point_symms

    do atom1=1,no_atoms_in_sc

      ! Rotate atom coordinates and identify equivalent atom.
      do i=1,3
        ! LAPACK commented out because it isn't working. 9/1/2017
        ! pos_atom1p(i) = offset(i,n) &
        !             & + ddot(3,rotation(i,1,n),3,atom_pos(1,atom1),1)
        pos_atom1p(i) = offset(i,n) &
                    & + dot_product(rotation(i,:,n),atom_pos(:,atom1))
      enddo ! i
      atom1p=atom_at_pos(pos_atom1p,no_atoms_in_sc,atom_pos,sc_rec_vec)
      !write(*,*)n,atom1,atom1p,atom_at_pos(pos_atom1p,atom_pos)
      if(atom1p<=0)call errstop('POINT_SYMM','Please check that &
        &your atom coordinates satisfy the rotational symmetries that you &
        &have supplied.  NB, I have assumed that r''=b+Rr, where R is the &
        &rotation matrix and b is the translation.  This may be wrong.')
      if(ABS(mass(atom1)-mass(atom1p))>tol*mass(atom1))call &
        &errstop('POINT_SYMM','Atoms '//TRIM(i2s(atom1))//' and ' &
        &//TRIM(i2s(atom1p))//' are equivalent by rotational symmetry, &
        &but they have different masses.')

      do atom2=1,no_atoms_in_sc

        ! Rotate atom coordinates and identify equivalent atom.
        do i=1,3
          ! LAPACK commented out because it isn't working. 9/1/2017
          ! pos_atom2p(i)=offset(i,n)+ddot(3,rotation(i,1,n),3, &
          !   &atom_pos(1,atom2),1)
          pos_atom2p(i) = offset(i,n) &
                      & + dot_product(rotation(i,:,n),atom_pos(:,atom2))
        enddo ! i
        atom2p=atom_at_pos(pos_atom2p,no_atoms_in_sc,atom_pos,sc_rec_vec)
        !write(*,*)n,atom2,atom2p,atom_at_pos(pos_atom2p,atom_pos),pos_atom2p(1),pos_atom2p(2),pos_atom2p(3)
        if(atom2p<=0)call errstop('POINT_SYMM','Please check that &
          &your atom coordinates satisfy the rotational symmetries that &
          &you have supplied.  NB, I have assumed that r''=b+Rr, where R &
          &is the rotation matrix and b is the translation.  This may be &
          &wrong.')

        ! Apply rotation to maxtrix of force constants.  Record whether or
        ! not each force constant is well-defined.
        do i=1,3
          do j=1,3
            fc=0.d0
            well_defined=.TRUE.
            ip_loop : do ip=1,3
              if(ABS(rotation(i,ip,n))>tol)then
                product=0.d0
                do jp=1,3
                  if(ABS(rotation(j,jp,n))>tol)then
                    if(defined(atom1,ip,atom2,jp))then
                      product=product+rotation(j,jp,n) &
                        &*force_const(atom1,ip,atom2,jp)
                    else
                      well_defined=.FALSE.
                      exit ip_loop
                    endif ! Force constant defined
                  endif ! Rotation non-zero.
                enddo ! jp
                fc=fc+product*rotation(i,ip,n)
              endif ! Rotation non-zero.
            enddo ip_loop ! ip
            if(well_defined)then
              if(.NOT.defined(atom1p,i,atom2p,j))then
                ! A previously undefined force constant...
                force_const(atom1p,i,atom2p,j)=fc
                defined(atom1p,i,atom2p,j)=.TRUE.
                weight(atom1p,i,atom2p,j)=1
              else
                ! A previously defined force constant.  Average.
                force_const(atom1p,i,atom2p,j)=(weight(atom1p,i, &
                  &atom2p,j)*force_const(atom1p,i,atom2p,j)+fc) &
                  &/DBLE(weight(atom1p,i,atom2p,j)+1)
                weight(atom1p,i,atom2p,j)=weight(atom1p,i,atom2p,j)+1
              endif ! Element already defined.
            endif ! product well-defined
          enddo ! j
        enddo ! i

      enddo ! atom2

    enddo ! atom1

  enddo ! n

!    do atom1=1,no_atoms_in_sc
!      do atom2=1,no_atoms_in_sc
!        do i=1,3
!          do j=1,3
!            if(.NOT.defined(atom1,i,atom2,j))then
!              open(unit=7,file='next.dat',status='new',iostat=ierr)
!              if(ierr/=0)call errstop('POINT_SYMM','Unable to open next.dat.')
!              write(7,*)atom1,i
!              close(7)
!              STOP
!            endif ! defined
!          enddo ! j
!        enddo ! i
!      enddo ! atom2
!    enddo ! atom1
end subroutine

! ----------------------------------------------------------------------
! Apply all of the point symmetry operations in turn in order to complete
! the matrix of force constants.  This is done iteratively until
! the changes to the matrix of force constants have converged to within
! some tolerance.  This approach is valid for any geometry.
! ----------------------------------------------------------------------
subroutine point_symm_brute_force(tol,fc_scale,no_atoms_in_sc, &
   & no_DoF_prim,no_point_symms,sc_rec_vec,rotation,atom_pos,offset,   &
   & force_const,defined)
  use linear_algebra,only : ddot
  implicit none
  
  real(dp), intent(in) :: tol
  real(dp), intent(in) :: fc_scale
  integer,  intent(in) :: no_atoms_in_sc
  integer,  intent(in) :: no_DoF_prim
  integer,  intent(in) :: no_point_symms
  real(dp), intent(in) :: sc_rec_vec(3,3)
  real(dp), intent(in) :: rotation(:,:,:)
  real(dp), intent(in) :: atom_pos(:,:)
  real(dp), intent(in) :: offset(:,:)
  
  real(dp), intent(inout) :: force_const(:,:,:,:)
  logical,  intent(inout) :: defined(:,:,:,:)
  
  integer,parameter :: max_t=10000,min_t=10 ! Max & min number of impositions.
  
  integer :: atom1,atom1p,atom2,atom2p,i,j,n,ip,t,three_noDoFprim_sq
  real(dp) :: fc,pos_atom1p(3),pos_atom2p(3),max_diff
  logical :: last_iteration

  three_noDoFprim_sq=3*no_DoF_prim**2

  last_iteration=.FALSE.

  do t=1,max_t

    max_diff=0.d0

    do n=1,no_point_symms

      do atom1=1,no_atoms_in_sc

        ! Rotate atom coordinates and identify equivalent atom.
        do i=1,3
          ! LAPACK commented out because it isn't working. 9/1/2017
          ! pos_atom1p(i)=offset(i,n)+ddot(3,rotation(i,1,n),3, &
          !   &atom_pos(1,atom1),1)
          pos_atom1p(i) = offset(i,n) &
                      & + dot_product(rotation(i,:,n),atom_pos(:,atom1))
        enddo ! i
        atom1p=atom_at_pos(pos_atom1p,no_atoms_in_sc,atom_pos,sc_rec_vec)
        if(atom1p<=0)call errstop('POINT_SYMM_BRUTE_FORCE','Please check &
          &that your atom coordinates satisfy the rotational symmetries that &
          &you have supplied.  NB, I have assumed that r''=b+Rr, where R is &
          &the rotation matrix and b is the translation.  This may be &
          &wrong.')

        do atom2=1,no_atoms_in_sc

          ! Rotate atom coordinates and identify equivalent atom.
          do i=1,3
            ! LAPACK commented out because it isn't working. 9/1/2017
            ! pos_atom2p(i)=offset(i,n)+ddot(3,rotation(i,1,n),3, &
            !   &atom_pos(1,atom2),1)
            pos_atom2p(i) = offset(i,n) &
                        & + dot_product(rotation(i,:,n),atom_pos(:,atom2))
          enddo ! i
          atom2p=atom_at_pos(pos_atom2p,no_atoms_in_sc,atom_pos,sc_rec_vec)
          if(atom2p<=0)call errstop('POINT_SYMM_BRUTE_FORCE','Please check &
            &that your atom coordinates satisfy the rotational symmetries &
            &that you have supplied.  NB, I have assumed that r''=b+Rr, &
            &where R is the rotation matrix and b is the translation.  &
            &This may be wrong.')

          do i=1,3
            do j=1,3
              if(.NOT.defined(atom1p,i,atom2p,j))then
                fc=0.d0
                do ip=1,3
                  ! LAPACK commented out because it isn't working. 9/1/2017
                  ! fc=fc+ddot(3,rotation(j,1,n),3, &
                  !   &force_const(atom1,ip,atom2,1),three_noDoFprim_sq) &
                  !   &*rotation(i,ip,n)
                  fc = fc                                        &
                   & + dot_product(rotation(j,:,n),              &
                   &               force_const(atom1,ip,atom2,:)) &
                   & * rotation(i,ip,n)
                enddo ! ip
                if(ABS(force_const(atom1p,i,atom2p,j)-fc)>max_diff)then 
                  max_diff=ABS(force_const(atom1p,i,atom2p,j)-fc)
                endif
                force_const(atom1p,i,atom2p,j)=fc
                if(last_iteration)defined(atom1p,i,atom2p,j)=.TRUE.
              endif ! undefined
            enddo ! j
          enddo ! i

        enddo ! atom2

      enddo ! atom1

    enddo ! n

    if(last_iteration)exit
    if(max_diff<tol*fc_scale.AND.t>=min_t)last_iteration=.TRUE.
    if(t==max_t)call errstop('POINT_SYMM_BRUTE_FORCE','Unable to impose &
      &point symmetry on the matrix of force constants.')

  enddo ! Iterative impositions of symmetry.
end subroutine

! ----------------------------------------------------------------------
! Impose Newton's third law on the matrix of force constants.
! ----------------------------------------------------------------------
subroutine newtons_law(tol,fc_scale,no_atoms_in_prim,no_atoms_in_sc, &
   & force_const)
  implicit none
  
  real(dp), intent(in) :: tol
  real(dp), intent(in) :: fc_scale
  integer,  intent(in) :: no_atoms_in_prim
  integer,  intent(in) :: no_atoms_in_sc
  
  real(dp), intent(inout) :: force_const(:,:,:,:)
  
  integer, parameter :: min_t = 10     ! min number of impositions.
  integer, parameter :: max_t = 100000 ! max number of impositions.
  
  integer :: atom1,atom2,i,j,t
  real(dp) :: fc,sum1,rescale,max_diff

  do t=1,max_t

    max_diff=0.d0

    ! Impose Newton's third law on matrix of force consts.
    do atom1=1,no_atoms_in_sc
      do i=1,3
        do j=1,3
          sum1=0.d0
          do atom2=1,atom1-1
            sum1=sum1+force_const(atom1,i,atom2,j)
          enddo
          do atom2=atom1+1,no_atoms_in_sc
            sum1=sum1+force_const(atom1,i,atom2,j)
          enddo
          rescale=(force_const(atom1,i,atom1,j)+sum1) &
            &/DBLE(no_atoms_in_sc-1)
          if(ABS(rescale)>max_diff)max_diff=ABS(rescale)
          do atom2=1,atom1-1
            force_const(atom1,i,atom2,j)=force_const(atom1,i,atom2,j) &
              &-rescale
          enddo ! atom2
          do atom2=atom1+1,no_atoms_in_sc
            force_const(atom1,i,atom2,j)=force_const(atom1,i,atom2,j) &
              &-rescale
          enddo ! atom2
        enddo ! j
      enddo ! i
    enddo ! atom1

    ! Impose symmetry on the matrix of force constants.
    do atom1=1,no_atoms_in_sc
      do i=1,3
        do atom2=1,no_atoms_in_sc
          do j=1,3
            fc=0.5d0*(force_const(atom1,i,atom2,j) &
              &+force_const(atom2,j,atom1,i))
            if(ABS(fc-force_const(atom1,i,atom2,j))>max_diff) &
              &max_diff=ABS(fc-force_const(atom1,i,atom2,j))
            if(ABS(fc-force_const(atom2,j,atom1,i))>max_diff) &
              &max_diff=ABS(fc-force_const(atom2,j,atom1,i))
            force_const(atom1,i,atom2,j)=fc
            force_const(atom2,j,atom1,i)=fc
          enddo ! j
        enddo ! atom2
      enddo ! i
    enddo ! atom1

    ! For monatomic crystals we have inversion symmetry too.
    ! See Ashcroft & Mermin, p438.
    if(no_atoms_in_prim==1)then
      do atom1=1,no_atoms_in_sc
        do i=1,3
          do atom2=1,no_atoms_in_sc
            do j=1,3
              fc=0.5d0*(force_const(atom1,i,atom2,j) &
                &+force_const(atom1,j,atom2,i))
              if(ABS(fc-force_const(atom1,i,atom2,j))>max_diff) &
                &max_diff=ABS(fc-force_const(atom1,i,atom2,j))
              if(ABS(fc-force_const(atom1,j,atom2,i))>max_diff) &
                &max_diff=ABS(fc-force_const(atom1,j,atom2,i))
              force_const(atom1,i,atom2,j)=fc
              force_const(atom1,j,atom2,i)=fc
            enddo ! j
          enddo ! atom2
        enddo ! i
      enddo ! atom1
    endif ! monatomic lattice

!      write(*,*)max_diff,tol,fc_scale,tol*fc_scale

    if(max_diff<(tol*fc_scale+tol*1.d-4).AND.t>=min_t)exit

    if(t==max_t)call errstop('NEWTONS_LAW', &
      &'Unable to impose Newton''s 3rd law. on matrix of force constants')

  enddo ! Iterative impositions of Newton's 3d law, etc.
end subroutine

! ----------------------------------------------------------------------
! Mass-reduce the matrix of force constants.
! ----------------------------------------------------------------------
subroutine mass_reduce(no_atoms_in_sc,mass,force_const)
  use linear_algebra, only : dscal
  implicit none
  
  integer,  intent(in) :: no_atoms_in_sc
  real(dp), intent(in) :: mass(:)
  
  real(dp), intent(inout) :: force_const(:,:,:,:)
  
  integer :: atom1,atom2,j
  real(dp) :: rec_root_mass(no_atoms_in_sc),rec_root_m1m2
  
  do atom1=1,no_atoms_in_sc
    rec_root_mass(atom1)=1.d0/sqrt(mass(atom1))
  enddo ! atom1
  do atom2=1,no_atoms_in_sc
    do atom1=1,no_atoms_in_sc
      rec_root_m1m2=rec_root_mass(atom1)*rec_root_mass(atom2)
      do j=1,3
        ! LAPACK commented out because it isn't working. 9/1/2017
        ! call dscal(3,rec_root_m1m2,force_const(atom1,1,atom2,j),no_atoms_in_sc)
        force_const(atom1,:,atom2,j) = force_const(atom1,:,atom2,j) &
                                   & * rec_root_m1m2
      enddo ! j
    enddo ! atom1
  enddo ! atom2
end subroutine

! ----------------------------------------------------------------------
! This subroutine evaluates and stores the primitive lattice vector
! associated with each atom.  It also evaluates an array holding the
! label of each atom as a function of the label of the primitive cell
! in which an atom is found and the label of the atom within the
! primitive cell.
! ----------------------------------------------------------------------
subroutine find_prim_cell(delta,no_atoms_in_sc,atom_pos,prim_lat_vec,   &
   & prim_rec_vec,sc_lat_vec,sc_rec_vec,no_prim_cells,no_atoms_in_prim, &
   & atom,atom_in_prim,prim_cell_for_atom,no_equiv_ims,delta_prim)
  implicit none
  
  real(dp), intent(in) :: delta
  integer,  intent(in) :: no_atoms_in_sc
  real(dp), intent(in) :: atom_pos(:,:)
  real(dp), intent(in) :: prim_lat_vec(3,3)
  real(dp), intent(in) :: prim_rec_vec(3,3)
  real(dp), intent(in) :: sc_lat_vec(3,3)
  real(dp), intent(in) :: sc_rec_vec(3,3)
  integer,  intent(in) :: no_prim_cells
  integer,  intent(in) :: no_atoms_in_prim
  
  integer,  intent(inout) :: atom(:,:)
  
  integer,  allocatable, intent(out) :: atom_in_prim(:)
  integer,  allocatable, intent(out) :: prim_cell_for_atom(:)
  integer,  allocatable, intent(out) :: no_equiv_ims(:,:,:)
  real(dp), allocatable, intent(out) :: delta_prim(:,:,:,:,:)
  
  integer :: n,n1,n1p,n1pp,n2,n2p,n2pp,n3,n3p,n3pp,m,label,im, &
    &prim_n(3,no_atoms_in_sc),p,atom1,ialloc
  real(dp) :: delta_vect(3,3),r_temp(3),delta_r_ims(3,maxim),delta_r_corr(3)

  allocate(atom_in_prim(no_atoms_in_sc),prim_cell_for_atom(no_atoms_in_sc), &
    &stat=ialloc)
  if(ialloc/=0)call errstop('FIND_PRIM_CELL','Allocation error: &
    &atom_in_prim, etc.')

  ! Calculate the primitive-lattice point corresponding to the 
  ! primitive cell in which each atom lies.
  ! We try to make sure that there is no uncertainty in the primitive
  ! cell due to atoms sitting on the boundary between two primitive cells.
  delta_vect(1:3,1:3)=delta*prim_lat_vec(1:3,1:3)
  do n=1,no_atoms_in_sc
    r_temp(1:3)=atom_pos(1:3,n)+2.d0*delta_vect(1:3,1) &
      &+2.d0*delta_vect(1:3,2)+2.d0*delta_vect(1:3,3)
    n1=FLOOR(DOT_PRODUCT(r_temp(1:3),prim_rec_vec(1:3,1)))
    n1p=FLOOR(DOT_PRODUCT(r_temp(1:3)+delta_vect(1:3,1), &
      &prim_rec_vec(1:3,1)))
    n1pp=FLOOR(DOT_PRODUCT(r_temp(1:3)-delta_vect(1:3,1), &
      &prim_rec_vec(1:3,1)))
    if(n1/=n1p.OR.n1/=n1pp)call errstop('FIND_PRIM_CELL','Problem &
      &identifying unit cell in which atom lies [1].  Please try increasing &
      &the "delta" parameter in subroutine FIND_PRIM_CELL.')
    n2=FLOOR(DOT_PRODUCT(r_temp(1:3),prim_rec_vec(1:3,2)))
    n2p=FLOOR(DOT_PRODUCT(r_temp(1:3)+delta_vect(1:3,2), &
      &prim_rec_vec(1:3,2)))
    n2pp=FLOOR(DOT_PRODUCT(r_temp(1:3)-delta_vect(1:3,2), &
      &prim_rec_vec(1:3,2)))
    if(n2/=n2p.OR.n2/=n2pp)call errstop('FIND_PRIM_CELL','Problem &
      &identifying unit cell in which atom lies [2].  Please try increasing &
      &the "delta" parameter in subroutine FIND_PRIM_CELL.')
    n3=FLOOR(DOT_PRODUCT(r_temp(1:3),prim_rec_vec(1:3,3)))
    n3p=FLOOR(DOT_PRODUCT(r_temp(1:3)+delta_vect(1:3,3), &
      &prim_rec_vec(1:3,3)))
    n3pp=FLOOR(DOT_PRODUCT(r_temp(1:3)-delta_vect(1:3,3), &
      &prim_rec_vec(1:3,3)))
    if(n3/=n3p.OR.n3/=n3pp)call errstop('FIND_PRIM_CELL','Problem &
      &identifying unit cell in which atom lies [3].  Please try increasing &
      &the "delta" parameter in subroutine FIND_PRIM_CELL.')
    prim_n(1,n)=n1 ; prim_n(2,n)=n2 ; prim_n(3,n)=n3
  enddo ! n

  ! Establish a label for each different atom in the primitive cell,
  ! and evaluate this label for each atom.
  atom_in_prim(1:no_atoms_in_sc)=-1
  label=0
  do n=1,no_atoms_in_sc
    if(atom_in_prim(n)==-1)then
      label=label+1
      atom_in_prim(n)=label
      do m=n+1,no_atoms_in_sc
        ! Is difference of atom positions an integer multiple of the
        ! primitive reciprocal lattice vectors?  If so, same label.
        if(is_lat_point(atom_pos(1:3,m)-atom_pos(1:3,n),prim_rec_vec)) &
          &atom_in_prim(m)=label
      enddo ! m
    endif ! Atom not yet labelled by number within prim cell.
  enddo ! n
  if(label/=no_atoms_in_prim)call errstop('FIND_PRIM_CELL','Problem &
    &labelling the atoms in the primitive cell.')

  ! Establish a label for each different primitive cell, and evaluate
  ! this label for each atom.
  prim_cell_for_atom(1:no_atoms_in_sc)=-1
  label=0
  do n=1,no_atoms_in_sc
    if(prim_cell_for_atom(n)==-1)then
      label=label+1
      prim_cell_for_atom(n)=label
      do m=n+1,no_atoms_in_sc
        r_temp=DBLE(prim_n(1,m)-prim_n(1,n))*prim_lat_vec(1:3,1) &
          &+DBLE(prim_n(2,m)-prim_n(2,n))*prim_lat_vec(1:3,2) &
          &+DBLE(prim_n(3,m)-prim_n(3,n))*prim_lat_vec(1:3,3)
        if(is_lat_point(r_temp,sc_rec_vec))prim_cell_for_atom(m)=label
      enddo ! m
    endif ! Atom not yet labelled by number of primitive cell.
  enddo ! n
  if(label/=no_prim_cells)call errstop('FIND_PRIM_CELL','Problem labelling &
    &the primitive cells.')

  ! Construct array holding atom number for a given primitive cell and
  ! atom within the primitive cell.
  atom(:,:)=-1
  do n=1,no_atoms_in_sc
    atom(prim_cell_for_atom(n),atom_in_prim(n))=n
  enddo ! n
  if(ANY(atom(:,:)==-1))call errstop('FIND_PRIM_CELL','Problem defining atom &
    &labels.')

  ! Work out number of equivalent images and Delta Prim. Lattice Vectors
  ! for pairs of atoms (used in evaluation of force-constant matrix).
  allocate(no_equiv_ims(no_prim_cells,no_atoms_in_prim,no_atoms_in_prim), &
    &delta_prim(3,maxim,no_prim_cells,no_atoms_in_prim, &
    &no_atoms_in_prim),stat=ialloc)
  if(ialloc/=0)call errstop('FIND_PRIM_CELL','Allocation error: &
    &no_equiv_ims, etc.')
  delta_prim=0.d0
  do n=1,no_atoms_in_prim
    atom1=atom(1,n)
    do m=1,no_atoms_in_prim
      delta_r_corr=atom_pos(1:3,atom1)-atom_pos(1:3,atom(1,m))
      do p=1,no_prim_cells
        ! Work out min. image distance(s) between atoms (1,n) and (p,m).
        call min_images_brute_force(atom_pos(1:3,atom(p,m)) &
          &-atom_pos(1:3,atom1),sc_lat_vec,sc_rec_vec,delta_r_ims, &
          &no_equiv_ims(p,m,n))
        ! Turn this into the corresponding difference(s) of latt. vects.
        do im=1,no_equiv_ims(p,m,n)
          delta_prim(1:3,im,p,m,n)=delta_r_ims(1:3,im)+delta_r_corr
        enddo ! im
      enddo ! p
    enddo ! m
  enddo ! n
end subroutine

! ----------------------------------------------------------------------
! Construct the dynamical matrix for a given k vector.
! ----------------------------------------------------------------------
subroutine construct_dyn_matrix(kvec,delta_prim,no_DoF_prim,no_prim_cells, &
   & no_atoms_in_prim,atom,no_equiv_ims,force_const,dyn_mat)
  implicit none
  
  real(dp), intent(in) :: kvec(3)
  real(dp), intent(in) :: delta_prim(:,:,:,:,:)
  integer,  intent(in) :: no_DOF_prim
  integer,  intent(in) :: no_prim_cells
  integer,  intent(in) :: no_atoms_in_prim
  integer,  intent(in) :: atom(:,:)
  integer,  intent(in) :: no_equiv_ims(:,:,:)
  real(dp), intent(in) :: force_const(:,:,:,:)
  
  complex(dp), intent(inout) :: dyn_mat(no_DoF_prim,no_DoF_prim)
  
  integer :: p,n,m,i,j,index1,index2,atom1,im
  complex(dp) :: dm,tempc,expikdotD(no_prim_cells,no_atoms_in_prim,&
    &no_atoms_in_prim)
  real(dp) :: k_dot_D

  ! Precompute exp(-ik.(R-R')) to go in the dynamical matrix.
  do n=1,no_atoms_in_prim
    atom1=atom(1,n)
    do m=1,no_atoms_in_prim
      do p=1,no_prim_cells
        tempc=CMPLX(0.d0,0.d0,dp)
        do im=1,no_equiv_ims(p,m,n)
          k_dot_D=-DOT_PRODUCT(kvec,delta_prim(1:3,im,p,m,n))
          tempc=tempc+CMPLX(COS(k_dot_D),SIN(k_dot_D),dp)
        enddo ! im
        if(no_equiv_ims(p,m,n)>1)then
          expikdotD(p,m,n)=tempc/DBLE(no_equiv_ims(p,m,n))
        else
          expikdotD(p,m,n)=tempc
        endif ! number of images > 1.
      enddo ! p
    enddo ! m
  enddo ! n

  ! Evaluate the dynamical matrix.
  index1=0
  do n=1,no_atoms_in_prim
    atom1=atom(1,n)
    do i=1,3
      index1=index1+1
      index2=0
      do m=1,no_atoms_in_prim
        do j=1,3
          index2=index2+1
          dm=CMPLX(0.d0,0.d0,dp)
          do p=1,no_prim_cells
            dm=dm+force_const(atom(p,m),j,atom1,i)*expikdotD(p,m,n)
          enddo ! p
          dyn_mat(index1,index2)=dm
        enddo ! j
      enddo ! i
    enddo ! m
  enddo ! n

  ! Enforce Hermiticity on the dynamical matrix.
  do index1=1,no_DoF_prim
    dyn_mat(index1,index1)=CMPLX(real(dyn_mat(index1,index1),dp),0.d0,dp)
    do index2=index1+1,no_DoF_prim
      dm=0.5d0*(dyn_mat(index1,index2)+CONJG(dyn_mat(index2,index1)))
      dyn_mat(index1,index2)=dm
      dyn_mat(index2,index1)=CONJG(dm)
    enddo ! index2
  enddo ! index 1
end subroutine

! ----------------------------------------------------------------------
! For a given k, construct and diagonalise the dynamical matrix.
! ----------------------------------------------------------------------
subroutine calculate_eigenfreqs(kvec,delta_prim,no_DoF_prim, &
   & no_atoms_in_prim,no_prim_cells,atom,no_equiv_ims,force_const,omega)
  use linear_algebra, only : zheev
  implicit none
  
  real(dp), intent(in) :: kvec(3)
  real(dp), intent(in) :: delta_prim(:,:,:,:,:)
  integer,  intent(in) :: no_DoF_prim
  integer,  intent(in) :: no_atoms_in_prim
  integer,  intent(in) :: no_prim_cells
  integer,  intent(in) :: atom(:,:)
  integer,  intent(in) :: no_equiv_ims(:,:,:)
  real(dp), intent(in) :: force_const(:,:,:,:)
  
  real(dp),intent(out) :: omega(no_DoF_prim)
  
  complex(dp) :: dyn_mat(no_DoF_prim,no_DoF_prim),work(2*no_DoF_prim-1)
  real(dp) :: rwork(3*no_DoF_prim-2),minusomegasq(no_DoF_prim)
  integer :: n,m,info
  
  ! Construct dynamical matrix.
  call construct_dyn_matrix(kvec,delta_prim,no_DoF_prim,no_prim_cells, &
     & no_atoms_in_prim,atom,no_equiv_ims,force_const,dyn_mat)

  call zheev('N','U',no_DoF_prim,dyn_mat(1,1),no_DoF_prim,minusomegasq(1), &
    &work(1),2*no_DoF_prim-1,rwork(1),info)
  if(info/=0)call errstop('CALCULATE_EIGENFREQS','ZHEEV failed (1).  Error &
    &code: '//TRIM(i2s(info))//'.')

  ! Eigenvalues of dynamical matrix are minus the frequencies squared.
  ! The eigenvalues are in ascending order, so the +ve frequencies
  ! will be in ascending order.
  m=no_DoF_prim
  do n=1,no_DoF_prim
    if(minusomegasq(m)>=0.d0)then
      omega(n)=-SQRT(minusomegasq(m)) ! Unstable mode.
    else
      omega(n)=SQRT(-minusomegasq(m)) ! Stable mode.
    endif
    m=m-1
  enddo ! n
end subroutine

! ----------------------------------------------------------------------
! For a given k, construct and diagonalise the dynamical matrix.
! This subroutine returns the polarisation vectors as well.
! It is not optimised for speed.
! ----------------------------------------------------------------------
subroutine calculate_eigenfreqs_and_vecs(kvec,no_DoF_prim,delta_prim, &
   & no_atoms_in_prim,no_prim_cells,atom,no_equiv_ims,force_const,    &
   & omega,pol_vec)
  use linear_algebra,only : zscal,dznrm2,zheev,zcopy
  implicit none
  
  real(dp), intent(in) :: kvec(3)
  integer,  intent(in) :: no_DoF_prim
  real(dp), intent(in) :: delta_prim(:,:,:,:,:)
  integer,  intent(in) :: no_atoms_in_prim
  integer,  intent(in) :: no_prim_cells
  integer,  intent(in) :: atom(:,:)
  integer,  intent(in) :: no_equiv_ims(:,:,:)
  real(dp), intent(in) :: force_const(:,:,:,:)
  
  real(dp),    intent(out) :: omega(no_DoF_prim)
  complex(dp), intent(out) :: pol_vec(no_DoF_prim,no_DoF_prim)
  
  real(dp) :: rwork(3*no_DoF_prim-2), &
    &minusomegasq(no_DoF_prim)
  complex(dp) :: dyn_mat(no_DoF_prim,no_DoF_prim),work(2*no_DoF_prim-1)
  integer :: n,m,info

  ! Construct dynamical matrix.
  call construct_dyn_matrix(kvec,delta_prim,no_DoF_prim,no_prim_cells, &
     & no_atoms_in_prim,atom,no_equiv_ims,force_const,dyn_mat)

  call zheev('V','U',no_DoF_prim,dyn_mat(1,1),no_DoF_prim,minusomegasq(1), &
    &work(1),2*no_DoF_prim-1,rwork(1),info)
  if(info/=0)call errstop('CALCULATE_EIGENFREQS_AND_VECS', &
    &'ZHEEV failed (1).  Error code: '//TRIM(i2s(info))//'.')

  m=no_DoF_prim
  do n=1,no_DoF_prim 
! Modified by B. Monserrat to output the correct 'omega' and 'pol_vec'
    if(minusomegasq(m)>=0.d0)then
      omega(n)=-SQRT(minusomegasq(m)) ! Unstable mode.
    else
      omega(n)=SQRT(-minusomegasq(m)) ! Stable mode.
    endif
    call zcopy(no_DoF_prim,dyn_mat(1,m),1,pol_vec(1,n),1)
    m=m-1
  enddo ! n

    ! Attempt to make the eigenvectors real by multiplying
    ! each eigenvector by the conjugate of the largest element.
    ! Then normalise the eigenvectors to the number of atoms in the
    ! primitive cell.  So the polarisation vector for a monatomic
    ! lattice should be normalised to unity.
   !do n=1,no_DoF_prim
   ! largest_mag2=0.d0
   ! largest_k=-1
   ! do k=1,no_DoF_prim
   !   mag2_element=DBLE(pol_vec(k,n))**2+AIMAG(pol_vec(k,n))**2
   !   if(mag2_element>largest_mag2)then
   !     largest_mag2=mag2_element
   !     largest_k=k
   !   endif
   ! enddo ! k
   ! if(largest_k==-1)call errstop('CALCULATE_EIGENFREQS_AND_VECS', &
   !   &'Eigenvector appears to be zero!')
   ! scalefactor=(DBLE(no_atoms_in_prim)/(ABS(pol_vec(largest_k,n)) &
   !   &*dznrm2(no_DoF_prim,pol_vec(1,n),1)))*CONJG(pol_vec(largest_k,n))
   ! call zscal(no_DoF_prim,scalefactor,pol_vec(1,n),1)
   !enddo ! n
end subroutine

! ----------------------------------------------------------------------
! Calculate the frequency density-of-states by Monte Carlo sampling of
! the Brillouin zone.
! ----------------------------------------------------------------------
subroutine calculate_freq_dos(tol,prim_rec_vec,no_DoF_prim,        &
   & delta_prim,no_atoms_in_prim,no_prim_cells,atom,no_equiv_ims,force_const,&
   & freq_dos_filename,bin_width,freq_dos)
  use constants,      only : max_bin, no_samples, no_fdos_sets
  use linear_algebra, only : dscal
  use rand_no_gen,    only : ranx
  implicit none
  
  real(dp), intent(in) :: tol
  real(dp), intent(in) :: prim_rec_vec(3,3)
  integer,  intent(in) :: no_DoF_prim
  real(dp), intent(in) :: delta_prim(:,:,:,:,:)
  integer,  intent(in) :: no_atoms_in_prim
  integer,  intent(in) :: no_prim_cells
  integer,  intent(in) :: atom(:,:)
  integer,  intent(in) :: no_equiv_ims(:,:,:)
  real(dp), intent(in) :: force_const(:,:,:,:)
  
  character(*), intent(in) :: freq_dos_filename
  
  real(dp),              intent(out) :: bin_width
  real(dp), allocatable, intent(out) :: freq_dos(:,:)
  
  ! Number of preliminary samples of Brillouin zone to make in order to
  ! establish maximum and minimum frequencies.
  integer,parameter :: no_samples_trial=10000
  ! Our preliminary sampling of the Brillouin zone is imperfect.
  ! Multiply the highest frequency found by this factor when choosing the
  ! highest frequency bin.
  real(dp),parameter :: safety_factor=1.15d0
  
  real(dp) :: omega(no_DoF_prim),kvec(3),rec_bin_width,max_freq,min_freq, &
    &rec_no_fdos_sets
  integer :: j,i,n,bin,ialloc,ierr
  logical :: soft_modes,soft_modes_prelim

  ! Establish (approximate) maximum and minimum frequencies and hence
  ! choose the bin width.
  max_freq=-1.d0
  min_freq=HUGE(1.d0)
  do i=1,no_samples_trial
    kvec(1:3)=twopi*(ranx()*prim_rec_vec(1:3,1) &
      &+ranx()*prim_rec_vec(1:3,2)+ranx()*prim_rec_vec(1:3,3))
    call calculate_eigenfreqs(kvec,delta_prim,no_DoF_prim,no_atoms_in_prim, &
      & no_prim_cells,atom,no_equiv_ims,force_const,omega)
    if(omega(1)<min_freq)min_freq=omega(1)
    if(omega(no_DoF_prim)>max_freq)max_freq=omega(no_DoF_prim)
  enddo ! i
  soft_modes_prelim=(min_freq<-tol)
  if(soft_modes_prelim)write(*,*)'WARNING: soft modes present.'
  write(*,*)'In preliminary sampling, largest frequency is : ',max_freq
  write(*,*)'and lowest frequency is                       : ',min_freq
  if(max_freq<=0)call errstop('CALCULATE_FREQ_DOS','The crystal lattice is &
    &pathologically unstable.')
  bin_width=safety_factor*max_freq/DBLE(max_bin)
  rec_bin_width=1.d0/bin_width

  write(*,*)'Number of random k vectors                    : ' &
    &//TRIM(i2s(no_samples))
  write(*,*)'Number of frequency bins                      : ' &
    &//TRIM(i2s(max_bin+1))
  write(*,*)'Frequency bin width                           : ',bin_width
  write(*,*)'Number of DoS sets (for computing error bars) : ' &
    &//TRIM(i2s(no_fdos_sets))

  allocate(freq_dos(0:max_bin,1:no_fdos_sets),stat=ialloc)
  if(ialloc/=0)call errstop('CALCULATE_FREQ_DOS','Allocation error: &
    &freq_dos.')
  freq_dos(0:max_bin,1:no_fdos_sets)=0.d0
  soft_modes=.FALSE.

  do j=1,no_fdos_sets
    do i=1,no_samples
      kvec(1:3)=twopi*(ranx()*prim_rec_vec(1:3,1) &
        &+ranx()*prim_rec_vec(1:3,2)+ranx()*prim_rec_vec(1:3,3))
      call calculate_eigenfreqs(kvec,delta_prim,no_DoF_prim, &
        & no_atoms_in_prim,no_prim_cells,atom,no_equiv_ims,force_const,omega)
      if(omega(1)<-tol)soft_modes=.TRUE.
      do n=1,no_DoF_prim
        if(omega(n)>0.d0)then ! Only bin positive frequencies.
          bin=MAX(0,FLOOR(rec_bin_width*omega(n)))
          if(bin>max_bin)call errstop('CALCULATE_FREQ_DOS', &
            &'Have encountered a frequency too high to be binned.  Try &
            &increasing the no_samples_trial or safety_factor parameters &
            &in the code.')
          freq_dos(bin,j)=freq_dos(bin,j)+1.d0
        endif ! positive frequency.
      enddo ! n
    enddo ! i
  enddo ! j
  if(soft_modes.AND..NOT.soft_modes_prelim)write(*,*)'WARNING: soft modes &
    &present.'

  ! Normalise frequency DoS so that its integral is the number of
  ! degrees of freedom in the primitive cell.  Note that the total
  ! number of frequencies sampled is 3*no_samples*no_atoms_in_prim.
  ! (Imaginary frequencies are ignored, however.)
  call dscal((max_bin+1)*no_fdos_sets,1.d0/(DBLE(no_samples)*bin_width), &
    &freq_dos(0,1),1)

  ! Write out the frequency DoS.
  rec_no_fdos_sets=1.d0/DBLE(no_fdos_sets)
  open(unit=8,file=freq_dos_filename,status='replace',iostat=ierr)
  if(ierr/=0)call errstop('CALCULATE_FREQ_DOS','Error opening freq_dos.dat.')
  do bin=0,max_bin
    write(8,*)(DBLE(bin)+0.5d0)*bin_width, &
      &SUM(freq_dos(bin,1:no_fdos_sets))*rec_no_fdos_sets
  enddo ! bin
  close(8)
end subroutine

! ----------------------------------------------------------------------
! Use the frequency density-of-states to evaluate the lattice thermal
! energy of the crystal as a function of the temperature in Kelvin.
! Repeat this for each set of frequency DoS data, to estimate the error
! in the LTFE.
! ----------------------------------------------------------------------
subroutine calc_lte(bin_width,temperature,freq_dos,tdependence1_filename)
  use constants,      only : max_bin, no_fdos_sets
  use linear_algebra, only : ddot
  implicit none
  
  real(dp), intent(in) :: bin_width
  real(dp), intent(in) :: temperature
  real(dp), intent(in) :: freq_dos(:,:)
  
  character(*), intent(in) :: tdependence1_filename
  
  integer :: bin,j
  real(dp) :: omega,lte_val,lte_sq,E_H(0:max_bin),lte,lte_err
  
  do bin=0,max_bin
    ! omega is the frequency in the middle of the corresponding bin.
    omega=(DBLE(bin)+0.5d0)*bin_width
    ! Array of harmonic energies at each frequency.
    E_H(bin)=harmonic_energy(temperature,omega) 
  enddo ! bin
  lte=0.d0 ; lte_sq=0.d0
  do j=1,no_fdos_sets
    ! LAPACK commented out because it isn't working. 9/1/2017
    ! lte_val=ddot(max_bin+1,freq_dos(0,j),1,E_H(0),1)
    lte_val = dot_product(freq_dos(:,j),E_H(:))
    lte=lte+lte_val ; lte_sq=lte_sq+lte_val**2
  enddo ! j
  lte=bin_width*lte/DBLE(no_fdos_sets)
  lte_sq=bin_width**2*lte_sq/DBLE(no_fdos_sets)
  lte_err=SQRT((lte_sq-lte**2)/DBLE(no_fdos_sets-1))
  write(*,'(1x,a,es18.10,a,es10.2)')'Done.  LTE per primitive cell : ', &
    &lte,' +/- ',lte_err
  write(*,'(1x,a,es18.10,a,es10.2)')'Done.  LTE per primitive cell (eV) : ', &
    &lte*27.211396132d0,' +/- ',lte_err*27.211396132d0
  open(1,FILE=tdependence1_filename)
  write(1,*)lte*27.211396132d0
  close(1)
end subroutine

! ----------------------------------------------------------------------
! This function returns the mean thermal energy of an isolated harmonic
! oscillator of frequency omega (in a.u.).  T is the temperature in Kelvin.
! ----------------------------------------------------------------------
real(dp) function harmonic_energy(t,omega)
  implicit none
  
  real(dp), intent(in) :: T
  real(dp), intent(in) :: omega
  
  real(dp) :: denominator
  
  if(T<=0.d0)then
    ! Zero-point energy.
    harmonic_energy=0.5d0*omega
  else
    denominator=EXP(omega/(kB_au_per_K*T))-1.d0
    if(denominator>0.d0)then
      ! General case.
      harmonic_energy=(1.d0/denominator+0.5d0)*omega
    else
      ! High-temperature limit.
      harmonic_energy=kB_au_per_K*T
    endif ! denominator>0
  endif ! T=0
end function

! ----------------------------------------------------------------------
! Use the frequency density-of-states to evaluate the lattice thermal
! free energy of the crystal as a function of the temperature in Kelvin.
! Repeat this for each set of frequency DoS data, to estimate the error
! in the LTFE.
! ----------------------------------------------------------------------
subroutine calc_ltfe(bin_width,temperature,freq_dos,tdependence2_filename)
  use constants, only : max_bin, no_fdos_sets
  implicit none
  
  real(dp), intent(in) :: bin_width
  real(dp), intent(in) :: temperature
  real(dp), intent(in) :: freq_dos(:,:)
  
  character(*), intent(in) :: tdependence2_filename
  
  integer :: bin,j
  real(dp) :: omega,ltfe_sq,ltfe_val,FE_H(0:max_bin),ltfe,ltfe_err
  
  do bin=0,max_bin
    ! omega is the frequency in the middle of the corresponding bin.
    omega=(DBLE(bin)+0.5d0)*bin_width
    ! Array of harmonic energies at each frequency.
    FE_H(bin)=harmonic_free_energy(temperature,omega)
  enddo ! bin
  ltfe=0.d0 ; ltfe_sq=0.d0
  do j=1,no_fdos_sets
    ltfe_val=DOT_PRODUCT(freq_dos(:,j),FE_H(:))
    ltfe=ltfe+ltfe_val ; ltfe_sq=ltfe_sq+ltfe_val**2
  enddo ! j
  ltfe=bin_width*ltfe/DBLE(no_fdos_sets)
  ltfe_sq=bin_width**2*ltfe_sq/DBLE(no_fdos_sets)
  ltfe_err=SQRT((ltfe_sq-ltfe**2)/DBLE(no_fdos_sets-1))
  write(*,'(1x,a,es18.10,a,es10.2)')'and LTFE per primitive cell   : ', &
    &ltfe,' +/- ',ltfe_err
  write(*,'(1x,a,es18.10,a,es10.2)')'and LTFE per primitive cell (eV)  : ', &
    &ltfe*27.211396132d0,' +/- ',ltfe_err*27.211396132d0
  open(1,FILE=tdependence2_filename)
  write(1,*)ltfe*27.211396132d0
  close(1)
end subroutine

! ----------------------------------------------------------------------
! This function returns the mean free energy of an isolated harmonic
! oscillator of frequency omega (in a.u.).  T is the temperature in Kelvin.
! ----------------------------------------------------------------------
real(dp) function harmonic_free_energy(T,omega)
  IMPLICIT NONE
  real(dp),INTENT(in) :: T,omega
  real(dp) :: difference,kT
  if(T<=0.d0)then
    ! Zero-point energy.
    harmonic_free_energy=0.5d0*omega
  else
    kT=kB_au_per_K*T
    difference=1.d0-EXP(-omega/kT)
    if(difference>0.d0)then
      harmonic_free_energy=0.5d0*omega+kT*LOG(difference)
    else
      ! High-temperature limit.
      harmonic_free_energy=-HUGE(0.d0)
    endif
  endif
end function

! ----------------------------------------------------------------------
! This subroutine generates a dispersion_curve.dat file, which contains
! all the branches of the dispersion curve in a format that xmgrace 
! can read.  The branches of the dispersion curve are plotted against
! the total distance travelled along the specified lines.
! ----------------------------------------------------------------------
subroutine generate_disp_curve(no_DoF_prim,delta_prim,no_kspace_lines, &
   & disp_kpoints,no_atoms_in_prim,no_prim_cells,atom,no_equiv_ims,    &
   & force_const,                                                      &
   & dispersion_curve_filename)
  implicit none
  
  integer,  intent(in) :: no_DoF_prim
  real(dp), intent(in) :: delta_prim(:,:,:,:,:)
  integer,  intent(in) :: no_kspace_lines
  real(dp), intent(in) :: disp_kpoints(:,:)
  integer,  intent(in) :: no_atoms_in_prim
  integer,  intent(in) :: no_prim_cells
  integer,  intent(in) :: atom(:,:)
  integer,  intent(in) :: no_equiv_ims(:,:,:)
  real(dp), intent(in) :: force_const(:,:,:,:)
  
  character(*), intent(in) :: dispersion_curve_filename
  
  real(dp) :: k_dist,kvec(3),delta_k(3),k_step,omega(no_DoF_prim)
  integer :: i,j,k,total_no_kpoints,ialloc,ierr
  real(dp),allocatable :: disp_k_dist(:),branch(:,:)
  integer,parameter :: no_kpoints_per_line=1000

  write(*,*)
  write(*,*)'Number of k points per line in k space : ' &
    &//TRIM(i2s(no_kpoints_per_line))
  write(*,*)

  ! Total number of k points at which the dispersion curve is to be calc'd.
  total_no_kpoints=no_kspace_lines*no_kpoints_per_line
  allocate(disp_k_dist(total_no_kpoints), &
    &branch(no_DoF_prim,total_no_kpoints),stat=ialloc)
  if(ialloc/=0)call errstop('GENERATE_DISP_CURVE','Allocation error: &
    &disp_k_dist, etc.')

  ! The step-size in all but the last line is |k_stop-k_start|/no_steps,
  ! so that k_stop is not included in that line (because k_stop is the
  ! first point on the next line).  For the last line the step-size
  ! is |k_stop-k_start|/(no_steps-1), because k_stop has to be the final
  ! point plotted.
  k=0
  k_dist=0.d0
  do i=1,no_kspace_lines
    write(*,*)'Start of new line at k-space distance : ',k_dist
    kvec(1:3)=disp_kpoints(1:3,i-1)
    if(i<no_kspace_lines)then
      delta_k(1:3)=(disp_kpoints(1:3,i)-disp_kpoints(1:3,i-1)) &
        &/DBLE(no_kpoints_per_line)
    else
      delta_k(1:3)=(disp_kpoints(1:3,i)-disp_kpoints(1:3,i-1)) &
        &/DBLE(no_kpoints_per_line-1)
    endif ! Last line or not
    k_step=SQRT(DOT_PRODUCT(delta_k,delta_k))
    do j=1,no_kpoints_per_line
      k=k+1
      call calculate_eigenfreqs(kvec,delta_prim,no_DoF_prim, &
        & no_atoms_in_prim,no_prim_cells,atom,no_equiv_ims,force_const,omega)
      branch(1:no_DoF_prim,k)=omega(1:no_DoF_prim)
      disp_k_dist(k)=k_dist
      kvec(1:3)=kvec(1:3)+delta_k(1:3)
      k_dist=k_dist+k_step
    enddo ! j
  enddo ! i
  k_dist=k_dist-k_step
  write(*,*)'Final line ends at k-space distance   : ',k_dist

  open(unit=8,file=dispersion_curve_filename,status='replace',iostat=ierr)
  if(ierr/=0)call errstop('GENERATE_DISP_CURVE', &
    &'Error opening dispersion_curve.dat.')
  do j=1,no_DoF_prim
    do k=1,total_no_kpoints
      write(8,*)disp_k_dist(k),branch(j,k)
    enddo ! k
    write(8,*)'&'
  enddo ! j
  close(8)

  deallocate(disp_k_dist,branch)
end subroutine

! ----------------------------------------------------------------------
! This subroutine calculates the mean speed of sound by evaluating
! domega/dk at Gamma and averaging over all directions.  The
! directions are uniformly distributed.  This subroutine only works
! for monatomic crystals at present.  The polarisation vectors are
! calculated for each k, and the dot product of k with the
! polarisation vectors are calculated.  The branch with the largest dot
! product is the longitudinal branch, whilst the other two branches
! are the transverse modes.
! ----------------------------------------------------------------------
subroutine calculate_speed_sound(no_atoms_in_prim,no_prim_cells,no_DoF_prim, &
   & small_k_scale,prim_lat_vec,delta_prim,atom,no_equiv_ims,force_const)
  use rand_no_gen,only : ranx
  implicit none
  
  integer,  intent(in) :: no_atoms_in_prim
  integer,  intent(in) :: no_prim_cells
  integer,  intent(in) :: no_DoF_prim
  real(dp), intent(in) :: small_k_scale
  real(dp), intent(in) :: prim_lat_vec(3,3)
  real(dp), intent(in) :: delta_prim(:,:,:,:,:)
  integer,  intent(in) :: atom(:,:)
  integer,  intent(in) :: no_equiv_ims(:,:,:)
  real(dp), intent(in) :: force_const(:,:,:,:)
  
  real(dp) :: kvec(3),kmag,omega(3),cos_theta,sin_theta,phi,c_tr_tot,c_tr, &
    &c2_tr_tot,c2_tr,c_ln_tot,c_ln,c2_ln_tot,c2_ln,err_tr,err_ln,c(3), &
    &kunit(3),pol_vec_real(3,3),dot_prod(3),temp,c_tr_old,c_ln_old
  complex(dp) :: pol_vec(3,3)
  integer :: i,no_samples,k,k2
  real(dp),parameter :: err_tol=1.d-3
  integer,parameter :: max_samples=1000000
  logical,parameter :: verbose=.FALSE.

  if(no_atoms_in_prim/=1)call errstop('CALCULATE_SPEED_SOUND', &
    &'At the moment this program can only work out the speed of sound in &
    &materials with a single atom per primitive cell.  Sorry about that.')

  ! First guess at a suitable radius of sphere in k-space for computing
  ! derivatives of omega w.r.t. k.
  kmag=0.75d0*small_k_scale

  c_tr_old=HUGE(1.d0) ; c_ln_old=HUGE(1.d0)

  ! Reduce kmag until the calculated sound speeds have converged.

  do

    c_tr_tot=0.d0  ; c_ln_tot=0.d0
    c2_tr_tot=0.d0 ; c2_ln_tot=0.d0
    no_samples=max_samples

    do i=1,max_samples

      ! Choose random k vector on sphere of radius kmag.
      cos_theta=1.d0-2.d0*ranx() ; phi=ranx()*twopi
      sin_theta=SQRT(1.d0-cos_theta**2)
      kunit(1:3)=(/sin_theta*COS(phi),sin_theta*SIN(phi),cos_theta/)
      kvec(1:3)=kmag*kunit

      ! Calculate corresponding eigenfrequencies.
      call calculate_eigenfreqs_and_vecs(kvec,no_DoF_prim,delta_prim,   &
        & no_atoms_in_prim,no_prim_cells,atom,no_equiv_ims,force_const, &
        & omega,pol_vec)
      if(ANY(omega<0.d0))then
        write(*,*)'Imaginary frequencies found.'
        write(*,*)'In terms of the primitive reciprocal lattice vectors, &
          &the k-point is:'
        write(*,*)DOT_PRODUCT(kvec,prim_lat_vec(1:3,1)), &
          &DOT_PRODUCT(kvec,prim_lat_vec(1:3,2)), &
          &DOT_PRODUCT(kvec,prim_lat_vec(1:3,3))
        write(*,*)'The frequencies are:'
        write(*,*)omega
        call errstop('CALCULATE_SPEED_SOUND','Cannot calculate speed of &
          &sound for unstable lattices.')
      endif ! Soft modes.

      ! Speed of sound corresponding to first three (acoustic) branches.
      c(1:3)=omega(1:3)/kmag

      ! Work out dot products of corresponding polarisation vectors
      ! with unit vector in direction of wave vector.
      pol_vec_real(:,:)=real(pol_vec(:,:),dp)
      dot_prod(1)=ABS(DOT_PRODUCT(kunit,pol_vec_real(1:3,1)))
      dot_prod(2)=ABS(DOT_PRODUCT(kunit,pol_vec_real(1:3,2)))
      dot_prod(3)=ABS(DOT_PRODUCT(kunit,pol_vec_real(1:3,3)))

      ! Arrange the three sound speeds in ascending order of dot
      ! product of k with polarisation vector.  The third component
      ! should be the longitudinal one, the first two should be
      ! the transverse components.
      do k=1,2
        do k2=k+1,3
          if(dot_prod(k)>dot_prod(k2))then
            temp=dot_prod(k) ; dot_prod(k)=dot_prod(k2) ; dot_prod(k2)=temp
            temp=c(k) ; c(k)=c(k2) ; c(k2)=temp
          endif ! Swap needed
        enddo ! k2
      enddo ! k

      ! Accumulate sound-speed statistics.
      c_tr_tot=c_tr_tot+c(1)+c(2)
      c_ln_tot=c_ln_tot+c(3)
      c2_tr_tot=c2_tr_tot+c(1)**2+c(2)**2
      c2_ln_tot=c2_ln_tot+c(3)**2

      ! Check whether we have desired accuracy level.
      if(i>=20)then
        c2_ln=c2_ln_tot/DBLE(i)
        c_ln=c_ln_tot/DBLE(i)
        err_ln=SQRT((c2_ln-c_ln**2)/DBLE(i-1))
        if(err_ln<err_tol*c_ln)then
          no_samples=i
          exit
        endif
      endif ! i>20

    enddo ! i

    ! Mean & standard error in mean for transverse speed.
    c_tr=c_tr_tot/DBLE(2*no_samples)
    c2_tr=c2_tr_tot/DBLE(2*no_samples)
    err_tr=SQRT((c2_tr-c_tr**2)/DBLE(2*no_samples-1))

    ! Mean & standard error in mean for longitudinal speed.
    c_ln=c_ln_tot/DBLE(no_samples)
    c2_ln=c2_ln_tot/DBLE(no_samples)
    err_ln=SQRT((c2_ln-c_ln**2)/DBLE(no_samples-1))

    if(verbose)then
      if(no_samples==max_samples)write(*,*)'Warning: have not reached &
        &desired error bar.'
      write(*,*)'Radius of k-space sphere : ',kmag
      write(*,*)'Speed of sound (a.u.)'
      write(*,*)'Transverse   : ',c_tr,' +/- ',err_tr
      write(*,*)'Longitudinal : ',c_ln,' +/- ',err_ln
      write(*,*)
    endif ! verbose

    if(ABS(c_tr-c_tr_old)<2.d0*err_tr.AND.ABS(c_ln-c_ln_old)<2.d0*err_ln)exit
    c_tr_old=c_tr
    c_ln_old=c_ln

    kmag=kmag*0.75d0

  enddo ! reduce kmag

  if(verbose)write(*,*)'Final results:'
  write(*,*)'Radius of k-space sphere : ',kmag
  write(*,*)'Please check this is sensible by examining a dispersion curve.'
  write(*,*)

  write(*,*)'Speed of sound (a.u.)'
  write(*,*)'Transverse   : ',c_tr,' +/- ',err_tr
  write(*,*)'Longitudinal : ',c_ln,' +/- ',err_ln
  write(*,*)
end subroutine

! ----------------------------------------------------------------------
! Evaluate the set of phonon frequencies on the supercell G vectors.
! Average the corresponding energies (for testing purposes).
! Write out the real part of the non-mass-reduced polarisation vector, which
! is the pattern of displacement corresponding to the normal mode.
! ----------------------------------------------------------------------
subroutine evaluate_freqs_on_grid(no_prim_cells,no_atoms_in_prim,    &
   & no_atoms_in_sc,no_DoF_prim,temperature,sc_rec_vec,prim_lat_vec, &
   & atom_pos,delta_prim,atom,mass,no_equiv_ims,force_const,         &
   & prim_cell_for_atom,atom_in_prim,                                &
   & kpairs_filename,freq_grids_filename,disp_patterns_filename,     &
   & kdisp_patterns_filename,pol_vec_filename,gvectors_filename,     &
   & gvectors_frac_filename,error_filename)
  implicit none
  
  integer,  intent(in) :: no_prim_cells
  integer,  intent(in) :: no_atoms_in_prim
  integer,  intent(in) :: no_atoms_in_sc
  integer,  intent(in) :: no_DoF_prim
  real(dp), intent(in) :: temperature
  real(dp), intent(in) :: sc_rec_vec(3,3)
  real(dp), intent(in) :: prim_lat_vec(3,3)
  real(dp), intent(in) :: atom_pos(:,:)
  real(dp), intent(in) :: delta_prim(:,:,:,:,:)
  integer,  intent(in) :: atom(:,:)
  real(dp), intent(in) :: mass(:)
  integer,  intent(in) :: no_equiv_ims(:,:,:)
  real(dp), intent(in) :: force_const(:,:,:,:)
  integer,  intent(in) :: prim_cell_for_atom(:)
  integer,  intent(in) :: atom_in_prim(:)
  
  ! filenames
  character(*), intent(in) :: kpairs_filename
  character(*), intent(in) :: freq_grids_filename
  character(*), intent(in) :: disp_patterns_filename
  character(*), intent(in) :: kdisp_patterns_filename
  character(*), intent(in) :: pol_vec_filename
  character(*), intent(in) :: gvectors_filename
  character(*), intent(in) :: gvectors_frac_filename
  character(*), intent(in) :: error_filename
  
  integer :: ng,i,j,k,ig,ierr,index1,index2,p,n,atom1
  logical :: found,soft_modes
  real(dp) :: gnew(3),gnew1(3),gnew2(3),gvec(3,no_prim_cells),R0(3), &
    &omega(no_DoF_prim),E,F,rec_root_mass(no_atoms_in_prim),GdotR, &
    &disp_pattern(3),kdisp_pattern(3),tot_disp_patt
  complex(dp) :: pol_vec(no_DoF_prim,no_DoF_prim), &
    &non_mr_pol_vec(3,no_atoms_in_prim),expiGdotR(no_prim_cells), &
    &kpol_vec(3,no_atoms_in_prim)
  real(dp),parameter :: tol_omega=1.d-6 ! For judging whether modes are soft.
  real(dp),parameter :: tol_g=1.d-8 ! For equivalent +/- G-points.
  integer :: reference(no_prim_cells) ! For equivalent +/- G-points.
  real(dp) :: gfrac(3,no_prim_cells) ! Fractional G vectors
  real(dp) :: prefactor

  write(*,*)'Frequencies at each supercell G vector:'
  write(*,*)

  ! Evaluate the set of supercell G vectors in the Brillouin zone of the
  ! primitive cell.
  ng=0
  do k=0,no_prim_cells-1
    gnew2=DBLE(k)*sc_rec_vec(1:3,3)
    do j=0,no_prim_cells-1
      gnew1=gnew2+DBLE(j)*sc_rec_vec(1:3,2)
      do i=0,no_prim_cells-1
        gnew=gnew1+DBLE(i)*sc_rec_vec(1:3,1)
        found=.TRUE.
        do ig=1,ng
          if(is_lat_point(gnew(1:3)-gvec(1:3,ig),prim_lat_vec))then
            found=.FALSE.
            exit
          endif ! ig
        enddo ! ig
        if(found)then
          ng=ng+1
          if(ng>no_prim_cells)call errstop('EVALUATE_FREQS_ON_GRID', &
            &'Bug: too many G vectors.')
          gvec(1:3,ng)=gnew(1:3)
        endif ! found
      enddo ! i
    enddo ! j
  enddo ! k
  if(ng/=no_prim_cells)call errstop('EVALUATE_FREQS_ON_GRID', &
    &'Bug: too few G vectors.')
  

! Calculate +/- G-vector pairs
! First, write G-vectors as fractions of rec. latt. vecs.
  do ig=1,no_prim_cells
    gfrac(1,ig)=DOT_PRODUCT(gvec(1:3,ig),prim_lat_vec(1:3,1))
    gfrac(2,ig)=DOT_PRODUCT(gvec(1:3,ig),prim_lat_vec(1:3,2))
    gfrac(3,ig)=DOT_PRODUCT(gvec(1:3,ig),prim_lat_vec(1:3,3))
    do while(gfrac(1,ig)+0.01*gfrac(1,ig)>0.5d0.AND.ABS(gfrac(1,ig)-0.5d0)>tol_g)
      gfrac(1,ig)=gfrac(1,ig)-1.d0
    enddo  
    do while(gfrac(1,ig)-0.01*gfrac(1,ig)<-0.5d0.AND.ABS(gfrac(1,ig)+0.5d0)>tol_g)
      gfrac(1,ig)=gfrac(1,ig)+1.d0
    enddo  
    do while(gfrac(2,ig)+0.01*gfrac(2,ig)>0.5d0.AND.ABS(gfrac(2,ig)-0.5d0)>tol_g)
      gfrac(2,ig)=gfrac(2,ig)-1.d0
    enddo  
    do while(gfrac(2,ig)-0.01*gfrac(2,ig)<-0.5d0.AND.ABS(gfrac(2,ig)+0.5d0)>tol_g)
      gfrac(2,ig)=gfrac(2,ig)+1.d0
    enddo  
    do while(gfrac(3,ig)+0.01*gfrac(3,ig)>0.5d0.AND.ABS(gfrac(3,ig)-0.5d0)>tol_g)
      gfrac(3,ig)=gfrac(3,ig)-1.d0
    enddo  
    do while(gfrac(3,ig)-0.01*gfrac(3,ig)<-0.5d0.AND.ABS(gfrac(3,ig)+0.5d0)>tol_g)
      gfrac(3,ig)=gfrac(3,ig)+1.d0
    enddo  
  enddo ! ig
  ! Second, pair them up
  reference=0
  do k=1,no_prim_cells
    do j=1,k-1
    if((ABS(gfrac(1,k)+gfrac(1,j))<tol_g.OR.ABS(gfrac(1,k)+gfrac(1,j)-1)<tol_g.OR.ABS(gfrac(1,k)+gfrac(1,j)+1)<tol_g).AND.&
      &(ABS(gfrac(2,k)+gfrac(2,j))<tol_g.OR.ABS(gfrac(2,k)+gfrac(2,j)-1)<tol_g.OR.ABS(gfrac(2,k)+gfrac(2,j)+1)<tol_g).AND.&
      &(ABS(gfrac(3,k)+gfrac(3,j))<tol_g.OR.ABS(gfrac(3,k)+gfrac(3,j)-1)<tol_g.OR.ABS(gfrac(3,k)+gfrac(3,j)+1)<tol_g))then
      reference(k)=j 
      reference(j)=k
      exit
    endif
    enddo ! j
  enddo ! k
  open(1,FILE=kpairs_filename)
  do i=1,no_prim_cells
    write(1,*)i,reference(i)
  enddo ! i
  close(1)

  open(unit=8,file=freq_grids_filename,status='replace',iostat=ierr)
  if(ierr/=0)call errstop('EVALUATE_FREQS_ON_GRID', &
    &'Problem opening freqs_grid.dat.')
  open(unit=9,file=disp_patterns_filename,status='replace',iostat=ierr)
  open(unit=10,file=kdisp_patterns_filename,status='replace',iostat=ierr)
  open(unit=11,file=pol_vec_filename,status='replace',iostat=ierr)
  if(ierr/=0)call errstop('EVALUATE_FREQS_ON_GRID', &
    &'Problem opening disp_patterns.dat.')

  do n=1,no_atoms_in_prim
    rec_root_mass(n)=1.d0/SQRT(mass(atom(1,n))) ! 1/sqrt(m) in prim. cell.
  enddo ! n

! Modified by B. Monserrat to output G vectors to file
  open(19,FILE=gvectors_filename)
  open(20,FILE=gvectors_frac_filename)
  write(19,*) no_prim_cells
  write(20,*) no_prim_cells

! Evaluate the frequencies at each supercell G vector.
  E=0.d0  ;  F=0.d0
  soft_modes=.FALSE.
  R0=atom_pos(1:3,atom(1,1))
  do ig=1,no_prim_cells
    write(*,'(" G = (",es20.12,",",es20.12,",",es20.12,")")')twopi*gvec(1:3,ig)
    write(*,'(" G = (",es20.12,",",es20.12,",",es20.12,")")')gfrac(1:3,ig)
    write(19,*)twopi*gvec(1:3,ig)
    ! G-vectors as a fraction of the primitive reciprocal lattice vectors
    write(20,*)ig,gfrac(1,ig),gfrac(2,ig),gfrac(3,ig)
    if(reference(ig)==0)then
      call calculate_eigenfreqs_and_vecs(twopi*gvec(1:3,ig),no_DoF_prim, &
        & delta_prim,no_atoms_in_prim,no_prim_cells,atom,no_equiv_ims,   &
        & force_const,omega,pol_vec)
    else
      if(reference(ig)>ig)then
        call calculate_eigenfreqs_and_vecs(twopi*gvec(1:3,ig),no_DoF_prim, &
          & delta_prim,no_atoms_in_prim,no_prim_cells,atom,no_equiv_ims,   &
          & force_const,omega,pol_vec)
      else
        call calculate_eigenfreqs_and_vecs(twopi*gvec(1:3,reference(ig)), &
          & no_DoF_prim,                                                  &
          & delta_prim,no_atoms_in_prim,no_prim_cells,atom,no_equiv_ims,  &
          & force_const,omega,pol_vec)
      endif
    endif

! The negative is used because the matrix of force constants is the transpose of
! the usual expression in derivations that lead to a positive exponential
    do p=1,no_prim_cells
      if(reference(ig)==0)then
        GdotR=-DOT_PRODUCT(twopi*gvec(1:3,ig),atom_pos(1:3,atom(p,1))-R0)
      else
        if(reference(ig)>ig)then
          GdotR=-DOT_PRODUCT(twopi*gvec(1:3,ig),atom_pos(1:3,atom(p,1))-R0)
        else
          GdotR=-DOT_PRODUCT(twopi*gvec(1:3,reference(ig)),atom_pos(1:3,atom(p,1))-R0)
        endif
      endif
      expiGdotR(p)=CMPLX(COS(GdotR),SIN(GdotR),dp) ! Store exp(iG.R_p).
    enddo ! p

    do index2=1,no_DoF_prim
      write(*,*)'  omega = ',omega(index2)
      if(omega(index2)>tol_omega)then
! Ignore contributions from imaginary or zero frequencies.
        E=E+harmonic_energy(temperature,omega(index2))
        F=F+harmonic_free_energy(temperature,omega(index2))
      elseif(omega(index2)<-tol_omega)then
        soft_modes=.TRUE.
      endif ! omega>0
      write(8,*)omega(index2),1.d0
! Compute the non-mass-reduced polarisation vector.
      index1=0
      do n=1,no_atoms_in_prim
        do i=1,3
          index1=index1+1
          non_mr_pol_vec(i,n)=pol_vec(index1,index2)*rec_root_mass(n)
          kpol_vec(i,n)=pol_vec(index1,index2)
        enddo ! i
      enddo ! n
      if(omega(index2)<-tol_omega)then
        write(9,*)'Frequency           : ',omega(index2),' (SOFT)'
        write(10,*)'Frequency           : ',omega(index2),' (SOFT)'
        write(11,*)'Mode number     :',index2, '   Frequency           : ',omega(index2),' (SOFT)'
      else
        write(9,*)'Frequency           : ',omega(index2)
        write(10,*)'Frequency           : ',omega(index2)
        write(11,*)'Mode number     :',index2, '   Frequency           : ',omega(index2)
      endif ! soft freq.
      write(9,*)gvec(1:3,ig)*twopi
      write(10,*)gvec(1:3,ig)*twopi
      write(11,*)gvec(1:3,ig)*twopi
      write(9,*)'Displacement pattern for each atom:'
      write(10,*)'Displacement pattern for each atom:'
      write(11,*)'Polarisation vector:'
      disp_pattern=0.d0
      kdisp_pattern=0.d0
      tot_disp_patt=0.d0
      do atom1=1,no_atoms_in_sc
! Displacement pattern: polarisation vector times exp(iG.R).
        if(reference(ig)==0)then
          disp_pattern=real(non_mr_pol_vec(1:3,atom_in_prim(atom1)) & ! Note only the real part is taken
            &*expiGdotR(prim_cell_for_atom(atom1)),dp)
          tot_disp_patt=tot_disp_patt+sqrt(disp_pattern(1)**2+disp_pattern(2)**2+disp_pattern(3)**2)
          kdisp_pattern=real(kpol_vec(1:3,atom_in_prim(atom1))) ! & 
            !&*expiGdotR(prim_cell_for_atom(atom1))
          prefactor=1.d0
        else
          if(reference(ig)>ig)then
            disp_pattern=real(non_mr_pol_vec(1:3,atom_in_prim(atom1)) &
             &*expiGdotR(prim_cell_for_atom(atom1)),dp)
            tot_disp_patt=tot_disp_patt+sqrt(disp_pattern(1)**2+disp_pattern(2)**2+disp_pattern(3)**2)
            kdisp_pattern=real(kpol_vec(1:3,atom_in_prim(atom1)) &
             &*expiGdotR(prim_cell_for_atom(atom1)),dp)
            prefactor=SQRT(2.d0)
          else
            disp_pattern=AIMAG(non_mr_pol_vec(1:3,atom_in_prim(atom1)) &
             &*expiGdotR(prim_cell_for_atom(atom1)))
            tot_disp_patt=tot_disp_patt+sqrt(disp_pattern(1)**2+disp_pattern(2)**2+disp_pattern(3)**2)
            kdisp_pattern=AIMAG(kpol_vec(1:3,atom_in_prim(atom1)) &
             &*expiGdotR(prim_cell_for_atom(atom1)))
            prefactor=SQRT(2.d0)
          endif
        endif
        write(9,*)disp_pattern,prefactor
        write(10,*)kdisp_pattern,prefactor
        write(11,*)real(non_mr_pol_vec(1:3,atom_in_prim(atom1)))
        write(11,*)AIMAG(non_mr_pol_vec(1:3,atom_in_prim(atom1)))
      enddo ! atom1
      write(9,*)
      write(10,*)
      write(11,*)
    enddo ! index2
    write(*,*)
    if(tot_disp_patt<1.d-8)then
      open(111,file=error_filename)
      write(111,*)'The total displacement is:',tot_disp_patt
      close(111)
    endif
  enddo ! ig
  E=E/DBLE(no_prim_cells)  ;  F=F/DBLE(no_prim_cells)

  close(8)
  close(9)
  close(10)
  close(11)
  close(19)
  close(20)

  write(*,*)'Mean LTE per primitive cell  : ',E
  write(*,*)'Mean LTE per primitive cell (eV)  : ',E*27.211396132d0
  write(*,*)'Mean LTFE per primitive cell : ',F
  write(*,*)'Mean LTFE per primitive cell (eV): ',F*27.211396132d0
  if(soft_modes)write(*,*)'WARNING: soft modes are present.'
  write(*,*)'Frequencies written to freqs_grid.dat.'
  write(*,*)'Displacement patterns written to disp_patterns.dat.'
  write(*,*)
end subroutine

! ----------------------------------------------------------------------
! Write out the dynamical matrix at each supercell G-vector to a file
! ----------------------------------------------------------------------
subroutine write_dynamical_matrix(sc_rec_vec,prim_lat_vec,no_DoF_prim, &
   & no_prim_cells,delta_prim,no_atoms_in_prim,atom,no_equiv_ims,      &
   & force_const,                                                      &
   & dyn_mat_fileroot)
  implicit none
  
  real(dp), intent(in) :: sc_rec_vec(3,3)
  real(dp), intent(in) :: prim_lat_vec(3,3)
  integer,  intent(in) :: no_DoF_prim
  integer,  intent(in) :: no_prim_cells
  real(dp), intent(in) :: delta_prim(:,:,:,:,:)
  integer,  intent(in) :: no_atoms_in_prim
  integer,  intent(in) :: atom(:,:)
  integer,  intent(in) :: no_equiv_ims(:,:,:)
  real(dp), intent(in) :: force_const(:,:,:,:)
  
  character(*), intent(in) :: dyn_mat_fileroot
  
  integer :: ierr,ng,k,j,i,ig,atom1,cart1,index1,atom2,cart2,index2
  real(dp) :: gnew2(3),gnew1(3),gnew(3),gvec(3,no_prim_cells)
  complex(dp) :: dyn_mat(no_DoF_prim,no_DoF_prim)
  logical :: found

  ng=0
  do k=0,no_prim_cells-1
    gnew2=DBLE(k)*sc_rec_vec(1:3,3)
    do j=0,no_prim_cells-1
      gnew1=gnew2+DBLE(j)*sc_rec_vec(1:3,2)
      do i=0,no_prim_cells-1
        gnew=gnew1+DBLE(i)*sc_rec_vec(1:3,1)
        found=.TRUE.
        do ig=1,ng
          if(is_lat_point(gnew(1:3)-gvec(1:3,ig),prim_lat_vec))then
            found=.FALSE.
            exit
          endif ! ig
        enddo ! ig
        if(found)then
          ng=ng+1
          if(ng>no_prim_cells)call errstop('WRITE_DYNAMICAL_MATRIX','Too &
           &many G-vectors.')
          gvec(1:3,ng)=gnew(1:3)
        endif ! found
      enddo ! i
    enddo ! j
  enddo ! k
  if(ng/=no_prim_cells)call errstop('WRITE_DYNAMICAL_MATRIX','Too few &
   &G-vectors.')

  dyn_mat(1:no_DoF_prim,1:no_DoF_prim)=CMPLX(0.d0,0.d0,dp)

  do ig=1,ng
    open(unit=101,file=dyn_mat_fileroot//TRIM(i2s(ig))//'.dat',status='replace',iostat=ierr)
    if(ierr/=0)call errstop('WRITE_DYNAMICAL_MATRIX','Problem opening &
     &dyn_mat.'//TRIM(i2s(ig))//'.dat file.')
    call construct_dyn_matrix(twopi*gvec(:,ig),delta_prim,no_DoF_prim, &
      & no_prim_cells,no_atoms_in_prim,atom,no_equiv_ims,force_const,dyn_mat)
    atom1=0
    do index1=1,no_DoF_prim
      atom2=0
      if(MOD(index1,3)==1)then
        atom1=atom1+1
        cart1=1
      endif ! MOD(index1,3)==1
      do index2=1,no_DoF_prim
        if(MOD(index2,3)==1)then
          atom2=atom2+1
          cart2=1
        endif ! MOD(index2,3)==1
        write(101,*)atom1,cart1,atom2,cart2,real(dyn_mat(index1,index2)),AIMAG(dyn_mat(index1,index2))
        cart2=cart2+1
      enddo ! index2
      cart1=cart1+1
    enddo ! index1
    close(101)
  end do ! ig
end subroutine

! ----------------------------------------------------------------------
! Write out atoms in primitive cell in order.
! ----------------------------------------------------------------------
subroutine write_atoms_in_primitive_cell(atom_pos,prim_rec_vec, &
   & no_atoms_in_prim,atom,mass,atoms_in_primitive_cell_filename)
  implicit none
  
  real(dp), intent(in) :: atom_pos(:,:)
  real(dp), intent(in) :: prim_rec_vec(3,3)
  integer,  intent(in) :: no_atoms_in_prim
  integer,  intent(in) :: atom(:,:)
  real(dp), intent(in) :: mass(:)
  
  character(*), intent(in) :: atoms_in_primitive_cell_filename
  
  real(dp),parameter :: tol=1.d-10
  
  integer :: ierr,n,atom1,i
  real(dp) :: pos(3),frac(3)

  open(unit=102,file=atoms_in_primitive_cell_filename,status='replace',&
   &iostat=ierr)
  if(ierr/=0)call errstop('write_atoms_in_primitive_cell','problem opening &
   &atoms_in_primitive_cell.dat file.')

  do n=1,no_atoms_in_prim
    atom1=atom(1,n)
    pos=atom_pos(1:3,atom1)
    do i=1,3
      frac(i)=dot_product(pos(1:3),prim_rec_vec(1:3,i))
    enddo
    frac(1:3)=modulo(frac(1:3)+tol,1.d0)-tol
    write(102,*)mass(atom1),frac(1:3)
  enddo

  close(102)
end subroutine

! ----------------------------------------------------------------------
! Deallocate arrays.
! ----------------------------------------------------------------------
subroutine finalise(species,mass,atom_pos,force_const,defined,atom, &
   & no_equiv_ims,delta_prim,rotation,offset,disp_kpoints,freq_dos, &
   & atom_in_prim,prim_cell_for_atom)
  implicit none
  
  character(2), allocatable, intent(inout) :: species(:)
  real(dp),     allocatable, intent(inout) :: mass(:)
  real(dp),     allocatable, intent(inout) :: atom_pos(:,:)
  real(dp),     allocatable, intent(inout) :: force_const(:,:,:,:)
  logical,      allocatable, intent(inout) :: defined(:,:,:,:)
  integer,      allocatable, intent(inout) :: atom(:,:)
  integer,      allocatable, intent(inout) :: no_equiv_ims(:,:,:)
  real(dp),     allocatable, intent(inout) :: delta_prim(:,:,:,:,:)
  real(dp),     allocatable, intent(inout) :: rotation(:,:,:)
  real(dp),     allocatable, intent(inout) :: offset(:,:)
  real(dp),     allocatable, intent(inout) :: disp_kpoints(:,:)
  real(dp),     allocatable, intent(inout) :: freq_dos(:,:)
  integer,      allocatable, intent(inout) :: atom_in_prim(:)
  integer,      allocatable, intent(inout) :: prim_cell_for_atom(:)
  
  if(allocated(species))deallocate(species,mass,atom_pos,force_const, &
    &defined,atom,no_equiv_ims,delta_prim)
  if(allocated(rotation))deallocate(rotation,offset)
  if(allocated(disp_kpoints))deallocate(disp_kpoints)
  if(allocated(freq_dos))deallocate(freq_dos)
  if(allocated(atom_in_prim))deallocate(atom_in_prim,prim_cell_for_atom)
end subroutine

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine lte(tol,tol2,delta,lte_filename,freq_dos_filename,                &
    & tdependence1_filename,tdependence2_filename,dispersion_curve_filename, &
    & kpairs_filename,freq_grids_filename,disp_patterns_filename,            &
    & kdisp_patterns_filename,pol_vec_filename,gvectors_filename,            &
    & gvectors_frac_filename,error_filename,dyn_mat_fileroot,                &
    & atoms_in_primitive_cell_filename)
  use constants, only : dp
  use utils,     only : errstop, wordwrap
  implicit none
  
  ! ----------------------------------------
  ! Inputs
  ! ----------------------------------------
  real(dp), intent(in) :: tol
  real(dp), intent(in) :: tol2
  real(dp), intent(in) :: delta
  
  ! ----------------------------------------
  ! filenames
  ! ----------------------------------------
  character(*), intent(in) :: lte_filename
  character(*), intent(in) :: freq_dos_filename
  character(*), intent(in) :: tdependence1_filename
  character(*), intent(in) :: tdependence2_filename
  character(*), intent(in) :: dispersion_curve_filename
  character(*), intent(in) :: kpairs_filename
  character(*), intent(in) :: freq_grids_filename
  character(*), intent(in) :: disp_patterns_filename
  character(*), intent(in) :: kdisp_patterns_filename
  character(*), intent(in) :: pol_vec_filename
  character(*), intent(in) :: gvectors_filename
  character(*), intent(in) :: gvectors_frac_filename
  character(*), intent(in) :: error_filename
  character(*), intent(in) :: dyn_mat_fileroot ! will have *.dat appended
  character(*), intent(in) :: atoms_in_primitive_cell_filename
  
  ! ----------------------------------------
  ! times
  ! ----------------------------------------
  real :: t1,t2
  
  ! ----------------------------------------
  ! previously global variables
  ! ----------------------------------------
  real(dp) :: prim_lat_vec(3,3)
  real(dp) :: sc_lat_vec(3,3)
  real(dp) :: prim_rec_vec(3,3)
  real(dp) :: sc_rec_vec(3,3)
  real(dp) :: length_scale
  real(dp) :: vol_scale
  real(dp) :: small_k_scale
  real(dp) :: fc_scale
  real(dp) :: bin_width
  real(dp) :: temperature
  integer  :: prog_function
  integer  :: no_atoms_in_sc
  integer  :: no_prim_cells
  integer  :: no_atoms_in_prim
  integer  :: no_point_symms
  integer  :: no_kspace_lines
  integer  :: no_DoF_prim
  
  character(2), allocatable :: species(:)
  real(dp),     allocatable :: mass(:)
  real(dp),     allocatable :: atom_pos(:,:)
  real(dp),     allocatable :: rotation(:,:,:) 
  real(dp),     allocatable ::  force_const(:,:,:,:)
  real(dp),     allocatable :: offset(:,:)
  real(dp),     allocatable :: freq_dos(:,:)
  real(dp),     allocatable :: disp_kpoints(:,:)
  real(dp),     allocatable :: delta_prim(:,:,:,:,:)
  integer,      allocatable :: atom(:,:)
  integer,      allocatable :: no_equiv_ims(:,:,:)
  integer,      allocatable :: atom_in_prim(:)
  integer,      allocatable :: prim_cell_for_atom(:)
  logical,      allocatable :: defined(:,:,:,:)

  call CPU_TIME(t1)

  write(*,*)
  write(*,*)'LATTICE THERMAL ENERGY'
  write(*,*)'======================'
  write(*,*)

  write(*,*)'Reading data from lte.dat...'
  write(*,*)
  call read_lte(tol,lte_filename,prim_rec_vec,sc_rec_vec,fc_scale,           &
     & no_atoms_in_prim,no_DoF_prim,no_prim_cells,length_scale,small_k_scale,&
     & vol_scale,                                                            &
     & prim_lat_vec,sc_lat_vec,no_atoms_in_sc,species,mass,atom_pos,         &
     & force_const,defined,atom,rotation,offset,disp_kpoints,no_kspace_lines,&
     & no_point_symms,prog_function,temperature)
  write(*,*)'Finished reading input data.'
  write(*,*)

  write(*,*)'Applying point symmetries to the matrix of force constants...'
  call point_symm(tol,no_atoms_in_sc,no_point_symms,rotation, &
     & atom_pos,offset,sc_rec_vec,mass,defined,force_const)
  write(*,*)'Done.'
  write(*,*)

  if(ANY(.NOT.defined))then
    call wordwrap('WARNING: will impose symmetries on the matrix of force &
      &constants iteratively...')
    call point_symm_brute_force(tol,fc_scale,no_atoms_in_sc,     &
       & no_DoF_prim,no_point_symms,sc_rec_vec,rotation,atom_pos,offset, &
       & force_const,defined)
    write(*,*)'Done.'
    write(*,*)
  endif
  if(ANY(.NOT.defined))call errstop('LTE','Some elements of the matrix of &
    &force constants are still undefined.')

  call wordwrap('Imposing Newton''s third law and symmetry on the matrix of &
    &force constants...')
  call newtons_law(tol2,fc_scale,no_atoms_in_prim,no_atoms_in_sc, &
    & force_const)
  write(*,*)'Done.'
  write(*,*)

  write(*,*)'Performing mass reduction on the matrix of force constants...'
  call mass_reduce(no_atoms_in_sc,mass,force_const)
  write(*,*)'Done.'
  write(*,*)

  write(*,*)'Establishing the primitive lattice vector associated with &
    &each atom...'
  call find_prim_cell(delta,no_atoms_in_sc,atom_pos,prim_lat_vec,         &
     & prim_rec_vec,sc_lat_vec,sc_rec_vec,no_prim_cells,no_atoms_in_prim, &
     & atom,atom_in_prim,prim_cell_for_atom,no_equiv_ims,                 &
     & delta_prim)
  write(*,*)'Done.'
  write(*,*)

  if(prog_function==1)then

    write(*,*)'Calculating the frequency density-of-states function...'
    call calculate_freq_dos(tol,prim_rec_vec,no_DoF_prim,   &
       & delta_prim,no_atoms_in_prim,no_prim_cells,atom,no_equiv_ims, &
       & force_const,                                                 &
       & freq_dos_filename,bin_width,freq_dos)
    call wordwrap('Done.  Frequency density-of-states function written to &
      &freq_dos.dat.  (Please view this file using XMGrace.)')
    write(*,*)

    write(*,*)'Calculating the lattice thermal energy (LTE) and free energy &
      &(LTFE)...'
    call calc_lte(bin_width,temperature,freq_dos,tdependence1_filename)
    call calc_ltfe(bin_width,temperature,freq_dos,tdependence2_filename)
    write(*,*)

  elseif(prog_function==2)then

    write(*,*)'Calculating the requested dispersion curve.'
    call generate_disp_curve(no_DoF_prim,delta_prim,no_kspace_lines,    &
       & disp_kpoints,no_atoms_in_prim,no_prim_cells,atom,no_equiv_ims, &
       & force_const,                                                   &
       & dispersion_curve_filename)
    call wordwrap('Done.  dispersion_curve.dat has been generated.  (Please &
      &view this file using XMGrace.)')
    write(*,*)

  elseif(prog_function==3)then

    write(*,*)'Calculating the speed of sound.'
    call calculate_speed_sound(no_atoms_in_prim,no_prim_cells,no_DoF_prim, &
       & small_k_scale,prim_lat_vec,delta_prim,atom,no_equiv_ims,force_const)
    write(*,*)'Done.  Speed of sound calculated.'
    write(*,*)

  elseif(prog_function==4)then
    
    write(*,*)'Calculating the frequencies and displacement patterns on the &
      &G-vector grid.'
    call evaluate_freqs_on_grid(no_prim_cells,no_atoms_in_prim,          &
       & no_atoms_in_sc,no_DoF_prim,temperature,sc_rec_vec,prim_lat_vec, &
       & atom_pos,delta_prim,atom,mass,no_equiv_ims,force_const,         &
       & prim_cell_for_atom,atom_in_prim,                                &
       & kpairs_filename,freq_grids_filename,disp_patterns_filename,     &
       & kdisp_patterns_filename,pol_vec_filename,gvectors_filename,     &
       & gvectors_frac_filename,error_filename)
    write(*,*)'Done.  Frequencies and displacement patterns calculated.'
    write(*,*)
    call write_dynamical_matrix(sc_rec_vec,prim_lat_vec,no_DoF_prim,  &
       & no_prim_cells,delta_prim,no_atoms_in_prim,atom,no_equiv_ims, &
       &force_const,                                                  &
       & dyn_mat_fileroot)
    call write_atoms_in_primitive_cell(atom_pos,prim_rec_vec, &
      & no_atoms_in_prim,atom,mass,atoms_in_primitive_cell_filename)
    
  else
    
    call errstop('LTE','Program function should be 1, 2, 3 or 4.')
    
  endif ! prog_function.
  
  call finalise(species,mass,atom_pos,force_const,defined,atom, &
     & no_equiv_ims,delta_prim,rotation,offset,disp_kpoints,freq_dos, &
     & atom_in_prim,prim_cell_for_atom)
  
  call CPU_TIME(t2)
  
  write(*,*)'Program finished.  Time taken: ',t2-t1
  write(*,*)
end subroutine
end module
