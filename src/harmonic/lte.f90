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


! Miscellaneous utilities etc.
module phonons
  use min_images,     only : is_lat_point, min_images_brute_force, maxim
  use constants,      only : dp, third, twopi, kB_au_per_K
  use utils,          only : i2s, errstop
  use file_io,        only : open_read_file, open_write_file
  use linear_algebra, only : determinant33, inv_33
  implicit none
  
  PRIVATE
  
  ! subroutines and functions
  PUBLIC defined,read_lte,point_symm,point_symm_brute_force,newtons_law, &
    &mass_reduce,find_prim_cell,calculate_freq_dos,calc_lte,calc_ltfe, &
    &prog_function,generate_disp_curve,calculate_speed_sound,finalise, &
    &evaluate_freqs_on_grid,write_dynamical_matrix,write_atoms_in_primitive_cell
  
  ! Global variables
  REAL(dp) :: prim_lat_vec(3,3),sc_lat_vec(3,3),prim_rec_vec(3,3), &
    &sc_rec_vec(3,3),length_scale,vol_scale,small_k_scale,fc_scale,bin_width, &
    &temperature
  INTEGER :: prog_function,no_atoms_in_sc,no_prim_cells,no_atoms_in_prim, &
    &no_point_symms,no_kspace_lines,no_DoF_prim
  CHARACTER(2),ALLOCATABLE :: species(:)
  REAL(dp),ALLOCATABLE :: mass(:),atom_pos(:,:),rotation(:,:,:), &
    &force_const(:,:,:,:),offset(:,:),freq_dos(:,:),disp_kpoints(:,:), &
    &delta_prim(:,:,:,:,:)
  INTEGER,ALLOCATABLE :: atom(:,:),no_equiv_ims(:,:,:),atom_in_prim(:), &
    &prim_cell_for_atom(:)
  LOGICAL,ALLOCATABLE :: defined(:,:,:,:)

contains


  SUBROUTINE setup_geometry(tol)
    ! The reciprocal lattice vectors of the primitive cell and the supercell
    ! are calculated here.
    ! NOTE THAT RECIPROCAL LATTICE VECTORS DO NOT CONTAIN THE FACTOR OF 2*pi.
    USE linear_algebra
    IMPLICIT NONE
    
    real(dp), intent(in) :: tol
    
    REAL(dp) :: prim_cell_volume,sc_cell_volume

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
    IF(prim_cell_volume<tol*vol_scale)CALL errstop('SETUP_GEOMETRY', &
      &'Primitive lattice vectors should be linearly independent.  Please &
      &check your lattice vectors.')
    CALL inv_33(prim_lat_vec,prim_rec_vec)
    prim_rec_vec=TRANSPOSE(prim_rec_vec)

    ! Supercell reciprocal lattice vectors and volume.
    sc_cell_volume=ABS(determinant33(sc_lat_vec))
    IF(sc_cell_volume<tol*vol_scale)CALL errstop('SETUP_GEOMETRY', &
      &'Supercell lattice vectors should be linearly independent.  Please &
      &check your lattice vectors.')
    CALL inv_33(sc_lat_vec,sc_rec_vec)
    sc_rec_vec=TRANSPOSE(sc_rec_vec)

    ! "Small" distance in k space, determined by size of supercell.
    ! Factor of 2pi included.
    small_k_scale=twopi*sc_cell_volume**(-third)

    ! Number of unit cells.
    no_prim_cells=NINT(sc_cell_volume/prim_cell_volume)
    IF(ABS(DBLE(no_prim_cells)*prim_cell_volume-sc_cell_volume) &
      &>tol*vol_scale)CALL errstop('SETUP_GEOMETRY','Supercell volume should &
      &be an integer multiple of primitive-cell volume.  Please check your &
      &lattice vectors.')

  END SUBROUTINE setup_geometry


  INTEGER FUNCTION atom_at_pos(rvec)
    ! This function returns the number of the atom at rvec (up to translations
    ! through a supercell lattice vector.
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: rvec(3)
    INTEGER :: atom1
    DO atom1=1,no_atoms_in_sc
      !WRITE(*,*)is_lat_point(atom_pos(1:3,atom1)-rvec(1:3),sc_rec_vec),atom_pos(1:3,atom1),rvec(1:3)
      IF(is_lat_point(atom_pos(1:3,atom1)-rvec(1:3),sc_rec_vec))THEN
        atom_at_pos=atom1
        RETURN
      ENDIF ! Separation of rvec and atom posn. is a SC lattice vector.
    ENDDO ! atom1
    atom_at_pos=0
  END FUNCTION atom_at_pos


  SUBROUTINE read_lte(tol,lte_filename)
    ! Read in the data in lte.dat and allocate arrays, etc.
    IMPLICIT NONE
    
    ! inputs
    real(dp)    , intent(in) :: tol
    character(*), intent(in) :: lte_filename
    
    INTEGER :: ierr,n,i,j,no_force_c_supplied,atom1,atom2,dir1,dir2,ialloc
    REAL(dp) :: fc,check_matrix

    OPEN(unit=8,file=lte_filename,status='old',iostat=ierr)
    IF(ierr/=0)CALL errstop('READ_LTE','Problem opening lte.dat.')

    ! Primitive lattice vectors
    READ(8,*,err=20,END=20)
    READ(8,*,err=20,END=20)prim_lat_vec(1:3,1)
    READ(8,*,err=20,END=20)prim_lat_vec(1:3,2)
    READ(8,*,err=20,END=20)prim_lat_vec(1:3,3)
    WRITE(*,*)'Primitive lattice vectors (Cartesian components in rows, a.u.):'
    WRITE(*,*)prim_lat_vec(1:3,1)
    WRITE(*,*)prim_lat_vec(1:3,2)
    WRITE(*,*)prim_lat_vec(1:3,3)
    WRITE(*,*)

    ! Supercell lattice vectors
    READ(8,*,err=20,END=20)
    READ(8,*,err=20,END=20)sc_lat_vec(1:3,1)
    READ(8,*,err=20,END=20)sc_lat_vec(1:3,2)
    READ(8,*,err=20,END=20)sc_lat_vec(1:3,3)
    WRITE(*,*)'Supercell lattice vectors (Cartesian components in rows, a.u.):'
    WRITE(*,*)sc_lat_vec(1:3,1)
    WRITE(*,*)sc_lat_vec(1:3,2)
    WRITE(*,*)sc_lat_vec(1:3,3)
    WRITE(*,*)

    ! Construct reciprocal lattice vectors, etc.
    CALL setup_geometry(tol)
    WRITE(*,*)'Number of primitive unit cells     : '//TRIM(i2s(no_prim_cells))
    WRITE(*,*)
    WRITE(*,*)'Prim. rec. latt. vectors (Cart. cmpnts &
      &in rows, factor of 2pi inc., a.u.):'
    WRITE(*,*)prim_rec_vec(1:3,1)*twopi
    WRITE(*,*)prim_rec_vec(1:3,2)*twopi
    WRITE(*,*)prim_rec_vec(1:3,3)*twopi
    WRITE(*,*)

    WRITE(*,*)'Supercell rec. latt. vectors(Cart. cmpnts &
     &in rows, factor of 2pi inc., a.u.):'
    WRITE(*,*)sc_rec_vec(1:3,1)*twopi
    WRITE(*,*)sc_rec_vec(1:3,2)*twopi
    WRITE(*,*)sc_rec_vec(1:3,3)*twopi
    WRITE(*,*)

    ! Number of atoms.
    READ(8,*,err=20,END=20)
    READ(8,*,err=20,END=20)no_atoms_in_sc
    WRITE(*,*)'Number of atoms in supercell       : '//TRIM(i2s(no_atoms_in_sc))
    IF(no_atoms_in_sc<no_prim_cells)CALL errstop('READ_LTE', &
      &'Need more atoms in the supercell!')
    IF(MOD(no_atoms_in_sc,no_prim_cells)/=0)CALL errstop('READ_LTE', &
      &'Number of atoms in supercell should be a multiple of the number of &
      &primitive cells.')
    no_atoms_in_prim=no_atoms_in_sc/no_prim_cells
    ALLOCATE(species(no_atoms_in_sc),mass(no_atoms_in_sc), &
      &atom_pos(3,no_atoms_in_sc), &
      &force_const(no_atoms_in_sc,3,no_atoms_in_sc,3), &
      &defined(no_atoms_in_sc,3,no_atoms_in_sc,3), &
      &atom(no_prim_cells,no_atoms_in_prim),stat=ialloc)
    IF(ialloc/=0)CALL errstop('READ_LTE','Allocation error: species, etc.')
    defined(:,:,:,:)=.FALSE.
    force_const(:,:,:,:)=0.d0
    no_DoF_prim=3*no_atoms_in_prim
    WRITE(*,*)

    ! Read in atom positions, species and masses.
    ! Convert atom position from fractional coordinates (in file) to
    ! Cartesian coordinates (used in program).  Translate the atom
    ! coordinates into the supercell at the origin.
    READ(8,*,err=20,END=20)
    WRITE(*,*)'Species ; Mass (a.u.) ; Position (Cartesian coordinates, a.u.)'
    DO i=1,no_atoms_in_sc
      READ(8,*,err=20,END=20)species(i),mass(i),atom_pos(1:3,i)
      atom_pos(1:3,i)=atom_pos(1,i)*sc_lat_vec(1:3,1)+atom_pos(2,i) &
        &*sc_lat_vec(1:3,2)+atom_pos(3,i)*sc_lat_vec(1:3,3)
      WRITE(*,'(" ",a,"  ",f14.6," ",3("  ",f14.6))')species(i),mass(i), &
        &atom_pos(1:3,i)
      IF(mass(i)<=0.d0)CALL errstop('READ_LTE','Mass should be positive.')
    ENDDO ! i
    ! Check atoms aren't on top of each other.
    DO i=1,no_atoms_in_sc-1
      DO j=i+1,no_atoms_in_sc
        IF(is_lat_point(atom_pos(1:3,j)-atom_pos(1:3,i),sc_rec_vec))CALL &
          &errstop('READ_LTE','Atoms '//TRIM(i2s(i))//' and '//TRIM(i2s(j)) &
          &//' appear to be on top of one another.')
      ENDDO ! j
    ENDDO ! i
    WRITE(*,*)

    ! Read in point-symmetry operations.
    READ(8,*,err=20,END=20)
    READ(8,*,err=20,END=20)no_point_symms
    WRITE(*,*)'Number of point symmetries         : '//TRIM(i2s(no_point_symms))
    IF(no_point_symms<1)CALL errstop('READ_LTE','At least one point-symmetry &
      &rotation matrix (identity) must be supplied.')
    ALLOCATE(rotation(3,3,no_point_symms),offset(3,no_point_symms),stat=ialloc)
    IF(ialloc/=0)CALL errstop('READ_LTE','Allocation error: rotation, etc.')
    READ(8,*,err=20,END=20)
    DO n=1,no_point_symms
      READ(8,*,err=20,END=20)rotation(1:3,1,n)
      READ(8,*,err=20,END=20)rotation(1:3,2,n)
      READ(8,*,err=20,END=20)rotation(1:3,3,n)
      DO i=1,3
        DO j=i,3
          check_matrix=DOT_PRODUCT(rotation(1:3,i,n),rotation(1:3,j,n))
          IF((i==j.AND.ABS(check_matrix-1.d0)>tol) &
            &.OR.(i/=j.AND.ABS(check_matrix)>tol))CALL &
            &errstop('READ_LTE','Rotation matrix '//TRIM(i2s(n)) &
            &//' is not orthogonal!')
        ENDDO ! j
      ENDDO ! i
      READ(8,*,err=20,END=20)offset(1:3,n)
      ! Convert translation to Cartesians.
      offset(1:3,n)=offset(1,n)*sc_lat_vec(1:3,1) &
        &+offset(2,n)*sc_lat_vec(1:3,2)+offset(3,n)*sc_lat_vec(1:3,3)
    ENDDO ! n
    WRITE(*,*)'Have read in rotation matrices and translation vectors.'
    WRITE(*,*)

    ! Read in force constants supplied.
    READ(8,*,err=20,END=20)
    READ(8,*,err=20,END=20)no_force_c_supplied
    WRITE(*,*)'Number of force constants supplied : ' &
      &//TRIM(i2s(no_force_c_supplied))
    IF(no_force_c_supplied<=0)CALL errstop('READ_LTE', &
      &'Need to supply more force data!')
    READ(8,*,err=20,END=20)
    fc_scale=0.d0
    DO i=1,no_force_c_supplied
      READ(8,*,err=20,END=20)atom1,dir1,atom2,dir2,fc
      fc_scale=fc_scale+ABS(fc)
      CALL trans_symm(atom1,dir1,atom2,dir2,fc,tol)
    ENDDO ! i
    fc_scale=fc_scale/DBLE(no_force_c_supplied)
    WRITE(*,*)fc_scale
    WRITE(*,*)'Have read in the force-constant data and applied &
      &translational symmetry.'
    WRITE(*,*)

    READ(8,*,err=20,END=20)
    READ(8,*,err=20,END=20)prog_function
    IF(prog_function==1)THEN
      WRITE(*,*)'The mean thermal energy and the free energy will &
        &be calculated.'
    ELSEIF(prog_function==2)THEN
      WRITE(*,*)'A dispersion curve will be calculated.'
    ELSEIF(prog_function==3)THEN
      WRITE(*,*)'The speed of sound will be calculated.'
    ELSEIF(prog_function==4)THEN
      WRITE(*,*)'The phonon frequencies at the supercell G vectors will be &
        &calculated.'
    ELSE
      CALL errstop('READ_LTE','Program function must be either 1, 2, 3 or 4.')
    ENDIF ! prog_function

    READ(8,*,err=20,END=20)
    READ(8,*,err=20,END=20)temperature
    IF(prog_function==1)THEN
      WRITE(*,*)'Temperature (K)                    :',temperature
      IF(temperature==0.d0)WRITE(*,*)'(i.e. the zero-point energy is to be &
        &calculated.)'
      IF(temperature<0.d0)CALL errstop('READ_LTE', &
        &'Temperature should be non-negative.')
      WRITE(*,*)
    ENDIF ! LTE to be calculated.

    READ(8,*,err=20,END=20)
    READ(8,*,err=20,END=20)no_kspace_lines
    IF(prog_function==2)THEN
      WRITE(*,*)'Number of lines in k-space to plot     : ' &
        &//TRIM(i2s(no_kspace_lines))
      IF(no_kspace_lines<1)CALL errstop('READ_LTE', &
        &'Need to supply more lines in k-space!')
      ALLOCATE(disp_kpoints(3,0:no_kspace_lines),stat=ialloc)
      IF(ialloc/=0)CALL errstop('READ_LTE','Allocation error: disp_kpoints.')
    ENDIF
    READ(8,*,err=20,END=20)
    IF(prog_function==2)WRITE(*,*)'Points along walk in reciprocal space &
      &(Cartesian components in a.u.):'
    DO i=0,no_kspace_lines
      IF(prog_function==2)THEN
        READ(8,*,err=20,END=20)disp_kpoints(1:3,i)
        disp_kpoints(1:3,i)=twopi*(disp_kpoints(1,i)*prim_rec_vec(1:3,1) &
          &+disp_kpoints(2,i)*prim_rec_vec(1:3,2) &
          &+disp_kpoints(3,i)*prim_rec_vec(1:3,3))
        WRITE(*,'(3(" ",f16.8))')disp_kpoints(1:3,i)
      ELSE
        READ(8,*,err=20,END=20)
      ENDIF ! prog_function=2
    ENDDO ! i
    IF(prog_function==2)WRITE(*,*)'Have read in points for dispersion curve.'

    CLOSE(8)

    RETURN

    ! Error reading file...
20  CALL errstop('READ_LTE','Problem reading lte.dat.  Please check the format &
      &of the file.')

  END SUBROUTINE read_lte


  SUBROUTINE trans_symm(atom1,dir1,atom2,dir2,fc,tol)
    ! Take the force constant fc supplied for (atom1,dir1,atom2,dir2),
    ! place it in all translationally equivalent elements of the matrix of
    ! force constants.
    IMPLICIT NONE
    
    real(dp), intent(in) :: tol
    
    INTEGER,INTENT(in) :: atom1,dir1,atom2,dir2
    REAL(dp),INTENT(in) :: fc
    INTEGER :: atom1p,atom2p,no_translations
    REAL(dp) :: pos_atom2p(3),relpos_atom2_atom1(3)
    relpos_atom2_atom1(1:3)=atom_pos(1:3,atom2)-atom_pos(1:3,atom1)
    no_translations=0
    DO atom1p=1,no_atoms_in_sc
      IF(is_lat_point(atom_pos(1:3,atom1p)-atom_pos(1:3,atom1), &
        &prim_rec_vec))THEN
        ! atom1p and atom1 are equivalent under trans. symm.
        IF(ABS(mass(atom1p)-mass(atom1))>tol*mass(atom1))CALL &
          &errstop('TRANS_SYMM','Atoms '//TRIM(i2s(atom1))//' and ' &
          &//TRIM(i2s(atom1p))//' are equivalent by translational symmetry, &
          &but they have different masses.')
        pos_atom2p(1:3)=atom_pos(1:3,atom1p)+relpos_atom2_atom1(1:3)
        atom2p=atom_at_pos(pos_atom2p)
        IF(atom2p<=0)CALL errstop('TRANS_SYMM','Please check that your atom &
          &coordinates satisfy the translational symmetry they should have.')
        ! atom2p and atom2 are related to each other by the same translation
        ! as atom1p and atom1.
        IF(ABS(mass(atom2p)-mass(atom2))>tol*mass(atom2))CALL &
          &errstop('TRANS_SYMM','Atoms '//TRIM(i2s(atom2))//' and ' &
          &//TRIM(i2s(atom2p))//' are equivalent by translational &
          &symmetry, but they have different masses.')
        IF(defined(atom1p,dir1,atom2p,dir2))THEN
          WRITE(*,*)atom1p,dir1,atom2p,dir2
          CALL errstop('TRANS_SYMM', &
          &'Please check your atom coordinates -- you may have a repeated &
          &atom.  Alternatively, you may have conflicting force data for &
          &atoms that are equivalent under translational symmetry.')
        ENDIF
        force_const(atom1p,dir1,atom2p,dir2)=fc
        defined(atom1p,dir1,atom2p,dir2)=.TRUE.
        no_translations=no_translations+1
      ENDIF ! atom1p equiv to atom1
    ENDDO ! atom1p
    IF(no_translations/=no_prim_cells)CALL errstop('TRANS_SYMM', &
      &'Number of translationally equivalent atoms found differs from the &
      &number of primitive cells.  Please check the atom coordinates!')
  END SUBROUTINE trans_symm


  SUBROUTINE point_symm(tol)
    ! Apply all of the point symmetry operations in turn in order to complete
    ! the matrix of force constants.  Give a pair of atoms atom1 and atom2, we
    ! find the corresponding pair of atoms atom1p and atom2p after the 
    ! symmetry translation and rotation have been applied.  The matrix of
    ! force constants Phi then transforms as Phi' = R Phi R^T, where R is the
    ! rotation matrix and Phi and Phi' are the matrices of force constants
    ! between atom1 & atom2 and atom1p & atom2p, respectively.  Some elements
    ! of Phi may be unknown, but so long as they are multiplied by zero we can
    ! still work out the corresponding elements of Phi'.
    USE linear_algebra,ONLY : ddot
    IMPLICIT NONE
    
    real(dp), intent(in) :: tol
    
    INTEGER :: atom1,atom1p,atom2,atom2p,i,j,n,ip,jp,ierr,&
      &weight(no_atoms_in_sc,3,no_atoms_in_sc,3)
    REAL(dp) :: fc,product,pos_atom1p(3),pos_atom2p(3)
    LOGICAL :: well_defined

    ! We average over all the force constants that ought, by the full
    ! symmetry of the supercell, to be identical.  The weight array is
    ! used to perform this average.
    DO j=1,3
      DO atom2=1,no_atoms_in_sc
        DO i=1,3
          DO atom1=1,no_atoms_in_sc
            IF(defined(atom1,i,atom2,j))THEN
              weight(atom1,i,atom2,j)=1
            ELSE
              weight(atom1,i,atom2,j)=0
            ENDIF ! defined
    !        WRITE(*,*)atom1,i,atom2,j,weight(atom1,i,atom2,j)
          ENDDO ! atom1
        ENDDO ! i
      ENDDO ! atom2
    ENDDO ! j

    DO n=1,no_point_symms

      DO atom1=1,no_atoms_in_sc

        ! Rotate atom coordinates and identify equivalent atom.
        DO i=1,3
          pos_atom1p(i)=offset(i,n)+ddot(3,rotation(i,1,n),3,atom_pos(1,atom1),1)
        ENDDO ! i
        atom1p=atom_at_pos(pos_atom1p)
        !WRITE(*,*)n,atom1,atom1p,atom_at_pos(pos_atom1p)
        IF(atom1p<=0)CALL errstop('POINT_SYMM','Please check that &
          &your atom coordinates satisfy the rotational symmetries that you &
          &have supplied.  NB, I have assumed that r''=b+Rr, where R is the &
          &rotation matrix and b is the translation.  This may be wrong.')
        IF(ABS(mass(atom1)-mass(atom1p))>tol*mass(atom1))CALL &
          &errstop('POINT_SYMM','Atoms '//TRIM(i2s(atom1))//' and ' &
          &//TRIM(i2s(atom1p))//' are equivalent by rotational symmetry, &
          &but they have different masses.')

        DO atom2=1,no_atoms_in_sc

          ! Rotate atom coordinates and identify equivalent atom.
          DO i=1,3
            pos_atom2p(i)=offset(i,n)+ddot(3,rotation(i,1,n),3, &
              &atom_pos(1,atom2),1)
          ENDDO ! i
          atom2p=atom_at_pos(pos_atom2p)
          !WRITE(*,*)n,atom2,atom2p,atom_at_pos(pos_atom2p),pos_atom2p(1),pos_atom2p(2),pos_atom2p(3)
          IF(atom2p<=0)CALL errstop('POINT_SYMM','Please check that &
            &your atom coordinates satisfy the rotational symmetries that &
            &you have supplied.  NB, I have assumed that r''=b+Rr, where R &
            &is the rotation matrix and b is the translation.  This may be &
            &wrong.')

          ! Apply rotation to maxtrix of force constants.  Record whether or
          ! not each force constant is well-defined.
          DO i=1,3
            DO j=1,3
              fc=0.d0
              well_defined=.TRUE.
              ip_loop : DO ip=1,3
                IF(ABS(rotation(i,ip,n))>tol)THEN
                  product=0.d0
                  DO jp=1,3
                    IF(ABS(rotation(j,jp,n))>tol)THEN
                      IF(defined(atom1,ip,atom2,jp))THEN
                        product=product+rotation(j,jp,n) &
                          &*force_const(atom1,ip,atom2,jp)
                      ELSE
                        well_defined=.FALSE.
                        EXIT ip_loop
                      ENDIF ! Force constant defined
                    ENDIF ! Rotation non-zero.
                  ENDDO ! jp
                  fc=fc+product*rotation(i,ip,n)
                ENDIF ! Rotation non-zero.
              ENDDO ip_loop ! ip
              IF(well_defined)THEN
                IF(.NOT.defined(atom1p,i,atom2p,j))THEN
                  ! A previously undefined force constant...
                  force_const(atom1p,i,atom2p,j)=fc
                  defined(atom1p,i,atom2p,j)=.TRUE.
                  weight(atom1p,i,atom2p,j)=1
                ELSE
                  ! A previously defined force constant.  Average.
                  force_const(atom1p,i,atom2p,j)=(weight(atom1p,i, &
                    &atom2p,j)*force_const(atom1p,i,atom2p,j)+fc) &
                    &/DBLE(weight(atom1p,i,atom2p,j)+1)
                  weight(atom1p,i,atom2p,j)=weight(atom1p,i,atom2p,j)+1
                ENDIF ! Element already defined.
              ENDIF ! product well-defined
            ENDDO ! j
          ENDDO ! i

        ENDDO ! atom2

      ENDDO ! atom1

    ENDDO ! n

!    DO atom1=1,no_atoms_in_sc
!      DO atom2=1,no_atoms_in_sc
!        DO i=1,3
!          DO j=1,3
!            IF(.NOT.defined(atom1,i,atom2,j))THEN
!              OPEN(unit=7,file='next.dat',status='new',iostat=ierr)
!              IF(ierr/=0)CALL errstop('POINT_SYMM','Unable to open next.dat.')
!              WRITE(7,*)atom1,i
!              CLOSE(7)
!              STOP
!            ENDIF ! defined
!          ENDDO ! j
!        ENDDO ! i
!      ENDDO ! atom2
!    ENDDO ! atom1

  END SUBROUTINE point_symm


  SUBROUTINE point_symm_brute_force(tol)
    ! Apply all of the point symmetry operations in turn in order to complete
    ! the matrix of force constants.  This is done iteratively until
    ! the changes to the matrix of force constants have converged to within
    ! some tolerance.  This approach is valid for any geometry.
    USE linear_algebra,ONLY : ddot
    IMPLICIT NONE
    
    real(dp), intent(in) :: tol
    
    INTEGER :: atom1,atom1p,atom2,atom2p,i,j,n,ip,t,three_noDoFprim_sq
    REAL(dp) :: fc,pos_atom1p(3),pos_atom2p(3),max_diff
    LOGICAL :: last_iteration
    INTEGER,PARAMETER :: max_t=10000,min_t=10 ! Max & min number of impositions.

    three_noDoFprim_sq=3*no_DoF_prim**2

    last_iteration=.FALSE.

    DO t=1,max_t

      max_diff=0.d0

      DO n=1,no_point_symms

        DO atom1=1,no_atoms_in_sc

          ! Rotate atom coordinates and identify equivalent atom.
          DO i=1,3
            pos_atom1p(i)=offset(i,n)+ddot(3,rotation(i,1,n),3, &
              &atom_pos(1,atom1),1)
          ENDDO ! i
          atom1p=atom_at_pos(pos_atom1p)
          IF(atom1p<=0)CALL errstop('POINT_SYMM_BRUTE_FORCE','Please check &
            &that your atom coordinates satisfy the rotational symmetries that &
            &you have supplied.  NB, I have assumed that r''=b+Rr, where R is &
            &the rotation matrix and b is the translation.  This may be &
            &wrong.')

          DO atom2=1,no_atoms_in_sc

            ! Rotate atom coordinates and identify equivalent atom.
            DO i=1,3
              pos_atom2p(i)=offset(i,n)+ddot(3,rotation(i,1,n),3, &
                &atom_pos(1,atom2),1)
            ENDDO ! i
            atom2p=atom_at_pos(pos_atom2p)
            IF(atom2p<=0)CALL errstop('POINT_SYMM_BRUTE_FORCE','Please check &
              &that your atom coordinates satisfy the rotational symmetries &
              &that you have supplied.  NB, I have assumed that r''=b+Rr, &
              &where R is the rotation matrix and b is the translation.  &
              &This may be wrong.')

            DO i=1,3
              DO j=1,3
                IF(.NOT.defined(atom1p,i,atom2p,j))THEN
                  fc=0.d0
                  DO ip=1,3
                    fc=fc+ddot(3,rotation(j,1,n),3, &
                      &force_const(atom1,ip,atom2,1),three_noDoFprim_sq) &
                      &*rotation(i,ip,n)
                  ENDDO ! ip
                  IF(ABS(force_const(atom1p,i,atom2p,j)-fc)>max_diff)THEN 
                    max_diff=ABS(force_const(atom1p,i,atom2p,j)-fc)
                  ENDIF
                  force_const(atom1p,i,atom2p,j)=fc
                  IF(last_iteration)defined(atom1p,i,atom2p,j)=.TRUE.
                ENDIF ! undefined
              ENDDO ! j
            ENDDO ! i

          ENDDO ! atom2

        ENDDO ! atom1

      ENDDO ! n

      IF(last_iteration)EXIT
      IF(max_diff<tol*fc_scale.AND.t>=min_t)last_iteration=.TRUE.
      IF(t==max_t)CALL errstop('POINT_SYMM_BRUTE_FORCE','Unable to impose &
        &point symmetry on the matrix of force constants.')

    ENDDO ! Iterative impositions of symmetry.

  END SUBROUTINE point_symm_brute_force


  SUBROUTINE newtons_law(tol)
    ! Impose Newton's third law on the matrix of force constants.
    IMPLICIT NONE
    
    real(dp), intent(in) :: tol
    
    INTEGER :: atom1,atom2,i,j,t
    REAL(dp) :: fc,sum1,rescale,max_diff
    INTEGER,PARAMETER :: min_t=10,max_t=100000 ! Max & min number of impositions.

    DO t=1,max_t

      max_diff=0.d0

      ! Impose Newton's third law on matrix of force consts.
      DO atom1=1,no_atoms_in_sc
        DO i=1,3
          DO j=1,3
            sum1=0.d0
            DO atom2=1,atom1-1
              sum1=sum1+force_const(atom1,i,atom2,j)
            ENDDO
            DO atom2=atom1+1,no_atoms_in_sc
              sum1=sum1+force_const(atom1,i,atom2,j)
            ENDDO
            rescale=(force_const(atom1,i,atom1,j)+sum1) &
              &/DBLE(no_atoms_in_sc-1)
            IF(ABS(rescale)>max_diff)max_diff=ABS(rescale)
            DO atom2=1,atom1-1
              force_const(atom1,i,atom2,j)=force_const(atom1,i,atom2,j) &
                &-rescale
            ENDDO ! atom2
            DO atom2=atom1+1,no_atoms_in_sc
              force_const(atom1,i,atom2,j)=force_const(atom1,i,atom2,j) &
                &-rescale
            ENDDO ! atom2
          ENDDO ! j
        ENDDO ! i
      ENDDO ! atom1
 
      ! Impose symmetry on the matrix of force constants.
      DO atom1=1,no_atoms_in_sc
        DO i=1,3
          DO atom2=1,no_atoms_in_sc
            DO j=1,3
              fc=0.5d0*(force_const(atom1,i,atom2,j) &
                &+force_const(atom2,j,atom1,i))
              IF(ABS(fc-force_const(atom1,i,atom2,j))>max_diff) &
                &max_diff=ABS(fc-force_const(atom1,i,atom2,j))
              IF(ABS(fc-force_const(atom2,j,atom1,i))>max_diff) &
                &max_diff=ABS(fc-force_const(atom2,j,atom1,i))
              force_const(atom1,i,atom2,j)=fc
              force_const(atom2,j,atom1,i)=fc
            ENDDO ! j
          ENDDO ! atom2
        ENDDO ! i
      ENDDO ! atom1

      ! For monatomic crystals we have inversion symmetry too.
      ! See Ashcroft & Mermin, p438.
      IF(no_atoms_in_prim==1)THEN
        DO atom1=1,no_atoms_in_sc
          DO i=1,3
            DO atom2=1,no_atoms_in_sc
              DO j=1,3
                fc=0.5d0*(force_const(atom1,i,atom2,j) &
                  &+force_const(atom1,j,atom2,i))
                IF(ABS(fc-force_const(atom1,i,atom2,j))>max_diff) &
                  &max_diff=ABS(fc-force_const(atom1,i,atom2,j))
                IF(ABS(fc-force_const(atom1,j,atom2,i))>max_diff) &
                  &max_diff=ABS(fc-force_const(atom1,j,atom2,i))
                force_const(atom1,i,atom2,j)=fc
                force_const(atom1,j,atom2,i)=fc
              ENDDO ! j
            ENDDO ! atom2
          ENDDO ! i
        ENDDO ! atom1
      ENDIF ! monatomic lattice

!      WRITE(*,*)max_diff,tol,fc_scale,tol*fc_scale

      IF(max_diff<(tol*fc_scale+tol*1.d-4).AND.t>=min_t)EXIT

      IF(t==max_t)CALL errstop('NEWTONS_LAW', &
        &'Unable to impose Newton''s 3rd law. on matrix of force constants')

    ENDDO ! Iterative impositions of Newton's 3d law, etc.

  END SUBROUTINE newtons_law


  SUBROUTINE mass_reduce
    ! Mass-reduce the matrix of force constants.
    USE linear_algebra,ONLY : dscal
    IMPLICIT NONE
    INTEGER :: atom1,atom2,j
    REAL(dp) :: rec_root_mass(no_atoms_in_sc),rec_root_m1m2
    DO atom1=1,no_atoms_in_sc
      rec_root_mass(atom1)=1.d0/SQRT(mass(atom1))
    ENDDO ! atom1
    DO atom2=1,no_atoms_in_sc
      DO atom1=1,no_atoms_in_sc
        rec_root_m1m2=rec_root_mass(atom1)*rec_root_mass(atom2)
        DO j=1,3
          CALL dscal(3,rec_root_m1m2,force_const(atom1,1,atom2,j),no_atoms_in_sc)
        ENDDO ! j
      ENDDO ! atom1
    ENDDO ! atom2
  END SUBROUTINE mass_reduce


  SUBROUTINE find_prim_cell(delta)
    ! This subroutine evaluates and stores the primitive lattice vector
    ! associated with each atom.  It also evaluates an array holding the
    ! label of each atom as a function of the label of the primitive cell
    ! in which an atom is found and the label of the atom within the
    ! primitive cell.
    IMPLICIT NONE
    
    real(dp), intent(in) :: delta
    
    INTEGER :: n,n1,n1p,n1pp,n2,n2p,n2pp,n3,n3p,n3pp,m,label,im, &
      &prim_n(3,no_atoms_in_sc),p,atom1,ialloc
    REAL(dp) :: delta_vect(3,3),r_temp(3),delta_r_ims(3,maxim),delta_r_corr(3)

    ALLOCATE(atom_in_prim(no_atoms_in_sc),prim_cell_for_atom(no_atoms_in_sc), &
      &stat=ialloc)
    IF(ialloc/=0)CALL errstop('FIND_PRIM_CELL','Allocation error: &
      &atom_in_prim, etc.')

    ! Calculate the primitive-lattice point corresponding to the 
    ! primitive cell in which each atom lies.
    ! We try to make sure that there is no uncertainty in the primitive
    ! cell due to atoms sitting on the boundary between two primitive cells.
    delta_vect(1:3,1:3)=delta*prim_lat_vec(1:3,1:3)
    DO n=1,no_atoms_in_sc
      r_temp(1:3)=atom_pos(1:3,n)+2.d0*delta_vect(1:3,1) &
        &+2.d0*delta_vect(1:3,2)+2.d0*delta_vect(1:3,3)
      n1=FLOOR(DOT_PRODUCT(r_temp(1:3),prim_rec_vec(1:3,1)))
      n1p=FLOOR(DOT_PRODUCT(r_temp(1:3)+delta_vect(1:3,1), &
        &prim_rec_vec(1:3,1)))
      n1pp=FLOOR(DOT_PRODUCT(r_temp(1:3)-delta_vect(1:3,1), &
        &prim_rec_vec(1:3,1)))
      IF(n1/=n1p.OR.n1/=n1pp)CALL errstop('FIND_PRIM_CELL','Problem &
        &identifying unit cell in which atom lies [1].  Please try increasing &
        &the "delta" parameter in subroutine FIND_PRIM_CELL.')
      n2=FLOOR(DOT_PRODUCT(r_temp(1:3),prim_rec_vec(1:3,2)))
      n2p=FLOOR(DOT_PRODUCT(r_temp(1:3)+delta_vect(1:3,2), &
        &prim_rec_vec(1:3,2)))
      n2pp=FLOOR(DOT_PRODUCT(r_temp(1:3)-delta_vect(1:3,2), &
        &prim_rec_vec(1:3,2)))
      IF(n2/=n2p.OR.n2/=n2pp)CALL errstop('FIND_PRIM_CELL','Problem &
        &identifying unit cell in which atom lies [2].  Please try increasing &
        &the "delta" parameter in subroutine FIND_PRIM_CELL.')
      n3=FLOOR(DOT_PRODUCT(r_temp(1:3),prim_rec_vec(1:3,3)))
      n3p=FLOOR(DOT_PRODUCT(r_temp(1:3)+delta_vect(1:3,3), &
        &prim_rec_vec(1:3,3)))
      n3pp=FLOOR(DOT_PRODUCT(r_temp(1:3)-delta_vect(1:3,3), &
        &prim_rec_vec(1:3,3)))
      IF(n3/=n3p.OR.n3/=n3pp)CALL errstop('FIND_PRIM_CELL','Problem &
        &identifying unit cell in which atom lies [3].  Please try increasing &
        &the "delta" parameter in subroutine FIND_PRIM_CELL.')
      prim_n(1,n)=n1 ; prim_n(2,n)=n2 ; prim_n(3,n)=n3
    ENDDO ! n

    ! Establish a label for each different atom in the primitive cell,
    ! and evaluate this label for each atom.
    atom_in_prim(1:no_atoms_in_sc)=-1
    label=0
    DO n=1,no_atoms_in_sc
      IF(atom_in_prim(n)==-1)THEN
        label=label+1
        atom_in_prim(n)=label
        DO m=n+1,no_atoms_in_sc
          ! Is difference of atom positions an integer multiple of the
          ! primitive reciprocal lattice vectors?  If so, same label.
          IF(is_lat_point(atom_pos(1:3,m)-atom_pos(1:3,n),prim_rec_vec)) &
            &atom_in_prim(m)=label
        ENDDO ! m
      ENDIF ! Atom not yet labelled by number within prim cell.
    ENDDO ! n
    IF(label/=no_atoms_in_prim)CALL errstop('FIND_PRIM_CELL','Problem &
      &labelling the atoms in the primitive cell.')

    ! Establish a label for each different primitive cell, and evaluate
    ! this label for each atom.
    prim_cell_for_atom(1:no_atoms_in_sc)=-1
    label=0
    DO n=1,no_atoms_in_sc
      IF(prim_cell_for_atom(n)==-1)THEN
        label=label+1
        prim_cell_for_atom(n)=label
        DO m=n+1,no_atoms_in_sc
          r_temp=DBLE(prim_n(1,m)-prim_n(1,n))*prim_lat_vec(1:3,1) &
            &+DBLE(prim_n(2,m)-prim_n(2,n))*prim_lat_vec(1:3,2) &
            &+DBLE(prim_n(3,m)-prim_n(3,n))*prim_lat_vec(1:3,3)
          IF(is_lat_point(r_temp,sc_rec_vec))prim_cell_for_atom(m)=label
        ENDDO ! m
      ENDIF ! Atom not yet labelled by number of primitive cell.
    ENDDO ! n
    IF(label/=no_prim_cells)CALL errstop('FIND_PRIM_CELL','Problem labelling &
      &the primitive cells.')

    ! Construct array holding atom number for a given primitive cell and
    ! atom within the primitive cell.
    atom(:,:)=-1
    DO n=1,no_atoms_in_sc
      atom(prim_cell_for_atom(n),atom_in_prim(n))=n
    ENDDO ! n
    IF(ANY(atom(:,:)==-1))CALL errstop('FIND_PRIM_CELL','Problem defining atom &
      &labels.')

    ! Work out number of equivalent images and Delta Prim. Lattice Vectors
    ! for pairs of atoms (used in evaluation of force-constant matrix).
    ALLOCATE(no_equiv_ims(no_prim_cells,no_atoms_in_prim,no_atoms_in_prim), &
      &delta_prim(3,maxim,no_prim_cells,no_atoms_in_prim, &
      &no_atoms_in_prim),stat=ialloc)
    IF(ialloc/=0)CALL errstop('FIND_PRIM_CELL','Allocation error: &
      &no_equiv_ims, etc.')
    delta_prim=0.d0
    DO n=1,no_atoms_in_prim
      atom1=atom(1,n)
      DO m=1,no_atoms_in_prim
        delta_r_corr=atom_pos(1:3,atom1)-atom_pos(1:3,atom(1,m))
        DO p=1,no_prim_cells
          ! Work out min. image distance(s) between atoms (1,n) and (p,m).
          CALL min_images_brute_force(atom_pos(1:3,atom(p,m)) &
            &-atom_pos(1:3,atom1),sc_lat_vec,sc_rec_vec,delta_r_ims, &
            &no_equiv_ims(p,m,n))
          ! Turn this into the corresponding difference(s) of latt. vects.
          DO im=1,no_equiv_ims(p,m,n)
            delta_prim(1:3,im,p,m,n)=delta_r_ims(1:3,im)+delta_r_corr
          ENDDO ! im
        ENDDO ! p
      ENDDO ! m
    ENDDO ! n

  END SUBROUTINE find_prim_cell


  SUBROUTINE construct_dyn_matrix(kvec,dyn_mat)
! Construct the dynamical matrix for a given k vector.
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: kvec(3)
    COMPLEX(dp),INTENT(inout) :: dyn_mat(no_DoF_prim,no_DoF_prim)
    INTEGER :: p,n,m,i,j,index1,index2,atom1,im
    COMPLEX(dp) :: dm,tempc,expikdotD(no_prim_cells,no_atoms_in_prim,&
      &no_atoms_in_prim)
    REAL(dp) :: k_dot_D

    ! Precompute exp(-ik.(R-R')) to go in the dynamical matrix.
    DO n=1,no_atoms_in_prim
      atom1=atom(1,n)
      DO m=1,no_atoms_in_prim
        DO p=1,no_prim_cells
          tempc=CMPLX(0.d0,0.d0,dp)
          DO im=1,no_equiv_ims(p,m,n)
            k_dot_D=-DOT_PRODUCT(kvec,delta_prim(1:3,im,p,m,n))
            tempc=tempc+CMPLX(COS(k_dot_D),SIN(k_dot_D),dp)
          ENDDO ! im
          IF(no_equiv_ims(p,m,n)>1)THEN
            expikdotD(p,m,n)=tempc/DBLE(no_equiv_ims(p,m,n))
          ELSE
            expikdotD(p,m,n)=tempc
          ENDIF ! number of images > 1.
        ENDDO ! p
      ENDDO ! m
    ENDDO ! n

    ! Evaluate the dynamical matrix.
    index1=0
    DO n=1,no_atoms_in_prim
      atom1=atom(1,n)
      DO i=1,3
        index1=index1+1
        index2=0
        DO m=1,no_atoms_in_prim
          DO j=1,3
            index2=index2+1
            dm=CMPLX(0.d0,0.d0,dp)
            DO p=1,no_prim_cells
              dm=dm+force_const(atom(p,m),j,atom1,i)*expikdotD(p,m,n)
            ENDDO ! p
            dyn_mat(index1,index2)=dm
          ENDDO ! j
        ENDDO ! i
      ENDDO ! m
    ENDDO ! n

    ! Enforce Hermiticity on the dynamical matrix.
    DO index1=1,no_DoF_prim
      dyn_mat(index1,index1)=CMPLX(REAL(dyn_mat(index1,index1),dp),0.d0,dp)
      DO index2=index1+1,no_DoF_prim
        dm=0.5d0*(dyn_mat(index1,index2)+CONJG(dyn_mat(index2,index1)))
        dyn_mat(index1,index2)=dm
        dyn_mat(index2,index1)=CONJG(dm)
      ENDDO ! index2
    ENDDO ! index 1

  END SUBROUTINE construct_dyn_matrix


  SUBROUTINE calculate_eigenfreqs(kvec,omega)
    ! For a given k, construct and diagonalise the dynamical matrix.
    USE linear_algebra,ONLY : zheev
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: kvec(3)
    REAL(dp),INTENT(out) :: omega(no_DoF_prim)
    COMPLEX(dp) :: dyn_mat(no_DoF_prim,no_DoF_prim),work(2*no_DoF_prim-1)
    REAL(dp) :: rwork(3*no_DoF_prim-2),minusomegasq(no_DoF_prim)
    INTEGER :: n,m,info
    
    ! Construct dynamical matrix.
    CALL construct_dyn_matrix(kvec,dyn_mat)

    CALL zheev('N','U',no_DoF_prim,dyn_mat(1,1),no_DoF_prim,minusomegasq(1), &
      &work(1),2*no_DoF_prim-1,rwork(1),info)
    IF(info/=0)CALL errstop('CALCULATE_EIGENFREQS','ZHEEV failed (1).  Error &
      &code: '//TRIM(i2s(info))//'.')

    ! Eigenvalues of dynamical matrix are minus the frequencies squared.
    ! The eigenvalues are in ascending order, so the +ve frequencies
    ! will be in ascending order.
    m=no_DoF_prim
    DO n=1,no_DoF_prim
      IF(minusomegasq(m)>=0.d0)THEN
        omega(n)=-SQRT(minusomegasq(m)) ! Unstable mode.
      ELSE
        omega(n)=SQRT(-minusomegasq(m)) ! Stable mode.
      ENDIF
      m=m-1
    ENDDO ! n

  END SUBROUTINE calculate_eigenfreqs


  SUBROUTINE calculate_eigenfreqs_and_vecs(kvec,omega,pol_vec)
    ! For a given k, construct and diagonalise the dynamical matrix.
    ! This subroutine returns the polarisation vectors as well.
    ! It is not optimised for speed.
    USE linear_algebra,ONLY : zscal,dznrm2,zheev,zcopy
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: kvec(3)
    REAL(dp),INTENT(out) :: omega(no_DoF_prim)
    COMPLEX(dp),INTENT(out) :: pol_vec(no_DoF_prim,no_DoF_prim)
    REAL(dp) :: largest_mag2,mag2_element,rwork(3*no_DoF_prim-2), &
      &minusomegasq(no_DoF_prim)
    COMPLEX(dp) :: dyn_mat(no_DoF_prim,no_DoF_prim),work(2*no_DoF_prim-1), &
      &scalefactor
    INTEGER :: n,m,k,largest_k,info

    ! Construct dynamical matrix.
    CALL construct_dyn_matrix(kvec,dyn_mat)

    CALL zheev('V','U',no_DoF_prim,dyn_mat(1,1),no_DoF_prim,minusomegasq(1), &
      &work(1),2*no_DoF_prim-1,rwork(1),info)
    IF(info/=0)CALL errstop('CALCULATE_EIGENFREQS_AND_VECS', &
      &'ZHEEV failed (1).  Error code: '//TRIM(i2s(info))//'.')

    m=no_DoF_prim
    DO n=1,no_DoF_prim 
! Modified by B. Monserrat to output the correct 'omega' and 'pol_vec'
      IF(minusomegasq(m)>=0.d0)THEN
        omega(n)=-SQRT(minusomegasq(m)) ! Unstable mode.
      ELSE
        omega(n)=SQRT(-minusomegasq(m)) ! Stable mode.
      ENDIF
      CALL zcopy(no_DoF_prim,dyn_mat(1,m),1,pol_vec(1,n),1)
      m=m-1
    ENDDO ! n

      ! Attempt to make the eigenvectors real by multiplying
      ! each eigenvector by the conjugate of the largest element.
      ! Then normalise the eigenvectors to the number of atoms in the
      ! primitive cell.  So the polarisation vector for a monatomic
      ! lattice should be normalised to unity.
     !DO n=1,no_DoF_prim
     ! largest_mag2=0.d0
     ! largest_k=-1
     ! DO k=1,no_DoF_prim
     !   mag2_element=DBLE(pol_vec(k,n))**2+AIMAG(pol_vec(k,n))**2
     !   IF(mag2_element>largest_mag2)THEN
     !     largest_mag2=mag2_element
     !     largest_k=k
     !   ENDIF
     ! ENDDO ! k
     ! IF(largest_k==-1)CALL errstop('CALCULATE_EIGENFREQS_AND_VECS', &
     !   &'Eigenvector appears to be zero!')
     ! scalefactor=(DBLE(no_atoms_in_prim)/(ABS(pol_vec(largest_k,n)) &
     !   &*dznrm2(no_DoF_prim,pol_vec(1,n),1)))*CONJG(pol_vec(largest_k,n))
     ! CALL zscal(no_DoF_prim,scalefactor,pol_vec(1,n),1)
     !ENDDO ! n

  END SUBROUTINE calculate_eigenfreqs_and_vecs


  ! Calculate the frequency density-of-states by Monte Carlo sampling of
  ! the Brillouin zone.
  subroutine calculate_freq_dos(tol,freq_dos_filename)
    use constants,      only : max_bin, no_samples, no_fdos_sets
    use linear_algebra, only : dscal
    use rand_no_gen,    only : ranx
    implicit none
    
    ! inputs
    real(dp),     intent(in) :: tol
    character(*), intent(in) :: freq_dos_filename
    
    REAL(dp) :: omega(no_DoF_prim),kvec(3),rec_bin_width,max_freq,min_freq, &
      &rec_no_fdos_sets
    INTEGER :: j,i,n,bin,ialloc,ierr
    LOGICAL :: soft_modes,soft_modes_prelim
    ! Number of preliminary samples of Brillouin zone to make in order to
    ! establish maximum and minimum frequencies.
    INTEGER,PARAMETER :: no_samples_trial=10000
    ! Our preliminary sampling of the Brillouin zone is imperfect.
    ! Multiply the highest frequency found by this factor when choosing the
    ! highest frequency bin.
    REAL(dp),PARAMETER :: safety_factor=1.15d0

    ! Establish (approximate) maximum and minimum frequencies and hence
    ! choose the bin width.
    max_freq=-1.d0
    min_freq=HUGE(1.d0)
    DO i=1,no_samples_trial
      kvec(1:3)=twopi*(ranx()*prim_rec_vec(1:3,1) &
        &+ranx()*prim_rec_vec(1:3,2)+ranx()*prim_rec_vec(1:3,3))
      CALL calculate_eigenfreqs(kvec,omega)
      IF(omega(1)<min_freq)min_freq=omega(1)
      IF(omega(no_DoF_prim)>max_freq)max_freq=omega(no_DoF_prim)
    ENDDO ! i
    soft_modes_prelim=(min_freq<-tol)
    IF(soft_modes_prelim)WRITE(*,*)'WARNING: soft modes present.'
    WRITE(*,*)'In preliminary sampling, largest frequency is : ',max_freq
    WRITE(*,*)'and lowest frequency is                       : ',min_freq
    IF(max_freq<=0)CALL errstop('CALCULATE_FREQ_DOS','The crystal lattice is &
      &pathologically unstable.')
    bin_width=safety_factor*max_freq/DBLE(max_bin)
    rec_bin_width=1.d0/bin_width

    WRITE(*,*)'Number of random k vectors                    : ' &
      &//TRIM(i2s(no_samples))
    WRITE(*,*)'Number of frequency bins                      : ' &
      &//TRIM(i2s(max_bin+1))
    WRITE(*,*)'Frequency bin width                           : ',bin_width
    WRITE(*,*)'Number of DoS sets (for computing error bars) : ' &
      &//TRIM(i2s(no_fdos_sets))

    ALLOCATE(freq_dos(0:max_bin,1:no_fdos_sets),stat=ialloc)
    IF(ialloc/=0)CALL errstop('CALCULATE_FREQ_DOS','Allocation error: &
      &freq_dos.')
    freq_dos(0:max_bin,1:no_fdos_sets)=0.d0
    soft_modes=.FALSE.

    DO j=1,no_fdos_sets
      DO i=1,no_samples
        kvec(1:3)=twopi*(ranx()*prim_rec_vec(1:3,1) &
          &+ranx()*prim_rec_vec(1:3,2)+ranx()*prim_rec_vec(1:3,3))
        CALL calculate_eigenfreqs(kvec,omega)
        IF(omega(1)<-tol)soft_modes=.TRUE.
        DO n=1,no_DoF_prim
          IF(omega(n)>0.d0)THEN ! Only bin positive frequencies.
            bin=MAX(0,FLOOR(rec_bin_width*omega(n)))
            IF(bin>max_bin)CALL errstop('CALCULATE_FREQ_DOS', &
              &'Have encountered a frequency too high to be binned.  Try &
              &increasing the no_samples_trial or safety_factor parameters &
              &in the code.')
            freq_dos(bin,j)=freq_dos(bin,j)+1.d0
          ENDIF ! positive frequency.
        ENDDO ! n
      ENDDO ! i
    ENDDO ! j
    IF(soft_modes.AND..NOT.soft_modes_prelim)WRITE(*,*)'WARNING: soft modes &
      &present.'

    ! Normalise frequency DoS so that its integral is the number of
    ! degrees of freedom in the primitive cell.  Note that the total
    ! number of frequencies sampled is 3*no_samples*no_atoms_in_prim.
    ! (Imaginary frequencies are ignored, however.)
    CALL dscal((max_bin+1)*no_fdos_sets,1.d0/(DBLE(no_samples)*bin_width), &
      &freq_dos(0,1),1)

    ! Write out the frequency DoS.
    rec_no_fdos_sets=1.d0/DBLE(no_fdos_sets)
    OPEN(unit=8,file=freq_dos_filename,status='replace',iostat=ierr)
    IF(ierr/=0)CALL errstop('CALCULATE_FREQ_DOS','Error opening freq_dos.dat.')
    DO bin=0,max_bin
      WRITE(8,*)(DBLE(bin)+0.5d0)*bin_width, &
        &SUM(freq_dos(bin,1:no_fdos_sets))*rec_no_fdos_sets
    ENDDO ! bin
    CLOSE(8)

  END SUBROUTINE calculate_freq_dos


  ! Use the frequency density-of-states to evaluate the lattice thermal
  ! energy of the crystal as a function of the temperature in Kelvin.
  ! Repeat this for each set of frequency DoS data, to estimate the error
  ! in the LTFE.
  SUBROUTINE calc_lte(tdependence1_filename)
    use constants,      only : max_bin, no_fdos_sets
    use linear_algebra, only : ddot
    implicit none
    
    ! inputs
    character(*), intent(in) :: tdependence1_filename
    
    INTEGER :: bin,j
    REAL(dp) :: omega,lte_val,lte_sq,E_H(0:max_bin),lte,lte_err
    
    DO bin=0,max_bin
      ! omega is the frequency in the middle of the corresponding bin.
      omega=(DBLE(bin)+0.5d0)*bin_width
      ! Array of harmonic energies at each frequency.
      E_H(bin)=harmonic_energy(temperature,omega) 
    ENDDO ! bin
    lte=0.d0 ; lte_sq=0.d0
    DO j=1,no_fdos_sets
      lte_val=ddot(max_bin+1,freq_dos(0,j),1,E_H(0),1)
      lte=lte+lte_val ; lte_sq=lte_sq+lte_val**2
    ENDDO ! j
    lte=bin_width*lte/DBLE(no_fdos_sets)
    lte_sq=bin_width**2*lte_sq/DBLE(no_fdos_sets)
    lte_err=SQRT((lte_sq-lte**2)/DBLE(no_fdos_sets-1))
    WRITE(*,'(1x,a,es18.10,a,es10.2)')'Done.  LTE per primitive cell : ', &
      &lte,' +/- ',lte_err
    WRITE(*,'(1x,a,es18.10,a,es10.2)')'Done.  LTE per primitive cell (eV) : ', &
      &lte*27.211396132d0,' +/- ',lte_err*27.211396132d0
    OPEN(1,FILE=tdependence1_filename)
    WRITE(1,*)lte*27.211396132d0
    CLOSE(1)
  END SUBROUTINE calc_lte


  REAL(dp) FUNCTION harmonic_energy(T,omega)
    ! This function returns the mean thermal energy of an isolated harmonic
    ! oscillator of frequency omega (in a.u.).  T is the temperature in Kelvin.
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: T,omega
    REAL(dp) :: denominator
    IF(T<=0.d0)THEN
      ! Zero-point energy.
      harmonic_energy=0.5d0*omega
    ELSE
      denominator=EXP(omega/(kB_au_per_K*T))-1.d0
      IF(denominator>0.d0)THEN
        ! General case.
        harmonic_energy=(1.d0/denominator+0.5d0)*omega
      ELSE
        ! High-temperature limit.
        harmonic_energy=kB_au_per_K*T
      ENDIF ! denominator>0
    ENDIF ! T=0
  END FUNCTION harmonic_energy


  ! Use the frequency density-of-states to evaluate the lattice thermal
  ! free energy of the crystal as a function of the temperature in Kelvin.
  ! Repeat this for each set of frequency DoS data, to estimate the error
  ! in the LTFE.
  SUBROUTINE calc_ltfe(tdependence2_filename)
    use constants, only : max_bin, no_fdos_sets
    IMPLICIT NONE
    
    ! inputs
    character(*), intent(in) :: tdependence2_filename
    
    INTEGER :: bin,j
    REAL(dp) :: omega,ltfe_sq,ltfe_val,FE_H(0:max_bin),ltfe,ltfe_err
    
    DO bin=0,max_bin
      ! omega is the frequency in the middle of the corresponding bin.
      omega=(DBLE(bin)+0.5d0)*bin_width
      ! Array of harmonic energies at each frequency.
      FE_H(bin)=harmonic_free_energy(temperature,omega)
    ENDDO ! bin
    ltfe=0.d0 ; ltfe_sq=0.d0
    DO j=1,no_fdos_sets
      ltfe_val=DOT_PRODUCT(freq_dos(0:max_bin,j),FE_H(0:max_bin))
      ltfe=ltfe+ltfe_val ; ltfe_sq=ltfe_sq+ltfe_val**2
    ENDDO ! j
    ltfe=bin_width*ltfe/DBLE(no_fdos_sets)
    ltfe_sq=bin_width**2*ltfe_sq/DBLE(no_fdos_sets)
    ltfe_err=SQRT((ltfe_sq-ltfe**2)/DBLE(no_fdos_sets-1))
    WRITE(*,'(1x,a,es18.10,a,es10.2)')'and LTFE per primitive cell   : ', &
      &ltfe,' +/- ',ltfe_err
    WRITE(*,'(1x,a,es18.10,a,es10.2)')'and LTFE per primitive cell (eV)  : ', &
      &ltfe*27.211396132d0,' +/- ',ltfe_err*27.211396132d0
    OPEN(1,FILE=tdependence2_filename)
    WRITE(1,*)ltfe*27.211396132d0
    CLOSE(1)
  END SUBROUTINE calc_ltfe


  REAL(dp) FUNCTION harmonic_free_energy(T,omega)
    ! This function returns the mean free energy of an isolated harmonic
    ! oscillator of frequency omega (in a.u.).  T is the temperature in Kelvin.
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: T,omega
    REAL(dp) :: difference,kT
    IF(T<=0.d0)THEN
      ! Zero-point energy.
      harmonic_free_energy=0.5d0*omega
    ELSE
      kT=kB_au_per_K*T
      difference=1.d0-EXP(-omega/kT)
      IF(difference>0.d0)THEN
        harmonic_free_energy=0.5d0*omega+kT*LOG(difference)
      ELSE
        ! High-temperature limit.
        harmonic_free_energy=-HUGE(0.d0)
      ENDIF
    ENDIF
  END FUNCTION harmonic_free_energy


  ! This subroutine generates a dispersion_curve.dat file, which contains
  ! all the branches of the dispersion curve in a format that xmgrace 
  ! can read.  The branches of the dispersion curve are plotted against
  ! the total distance travelled along the specified lines.
  subroutine generate_disp_curve(dispersion_curve_filename)
    implicit none
    
    ! inputs
    character(*), intent(in) :: dispersion_curve_filename
    
    REAL(dp) :: k_dist,kvec(3),delta_k(3),k_step,omega(no_DoF_prim)
    INTEGER :: i,j,k,total_no_kpoints,ialloc,ierr
    REAL(dp),ALLOCATABLE :: disp_k_dist(:),branch(:,:)
    INTEGER,PARAMETER :: no_kpoints_per_line=1000

    WRITE(*,*)
    WRITE(*,*)'Number of k points per line in k space : ' &
      &//TRIM(i2s(no_kpoints_per_line))
    WRITE(*,*)

    ! Total number of k points at which the dispersion curve is to be calc'd.
    total_no_kpoints=no_kspace_lines*no_kpoints_per_line
    ALLOCATE(disp_k_dist(total_no_kpoints), &
      &branch(no_DoF_prim,total_no_kpoints),stat=ialloc)
    IF(ialloc/=0)CALL errstop('GENERATE_DISP_CURVE','Allocation error: &
      &disp_k_dist, etc.')

    ! The step-size in all but the last line is |k_stop-k_start|/no_steps,
    ! so that k_stop is not included in that line (because k_stop is the
    ! first point on the next line).  For the last line the step-size
    ! is |k_stop-k_start|/(no_steps-1), because k_stop has to be the final
    ! point plotted.
    k=0
    k_dist=0.d0
    DO i=1,no_kspace_lines
      WRITE(*,*)'Start of new line at k-space distance : ',k_dist
      kvec(1:3)=disp_kpoints(1:3,i-1)
      IF(i<no_kspace_lines)THEN
        delta_k(1:3)=(disp_kpoints(1:3,i)-disp_kpoints(1:3,i-1)) &
          &/DBLE(no_kpoints_per_line)
      ELSE
        delta_k(1:3)=(disp_kpoints(1:3,i)-disp_kpoints(1:3,i-1)) &
          &/DBLE(no_kpoints_per_line-1)
      ENDIF ! Last line or not
      k_step=SQRT(DOT_PRODUCT(delta_k,delta_k))
      DO j=1,no_kpoints_per_line
        k=k+1
        CALL calculate_eigenfreqs(kvec,omega)
        branch(1:no_DoF_prim,k)=omega(1:no_DoF_prim)
        disp_k_dist(k)=k_dist
        kvec(1:3)=kvec(1:3)+delta_k(1:3)
        k_dist=k_dist+k_step
      ENDDO ! j
    ENDDO ! i
    k_dist=k_dist-k_step
    WRITE(*,*)'Final line ends at k-space distance   : ',k_dist

    OPEN(unit=8,file=dispersion_curve_filename,status='replace',iostat=ierr)
    IF(ierr/=0)CALL errstop('GENERATE_DISP_CURVE', &
      &'Error opening dispersion_curve.dat.')
    DO j=1,no_DoF_prim
      DO k=1,total_no_kpoints
        WRITE(8,*)disp_k_dist(k),branch(j,k)
      ENDDO ! k
      WRITE(8,*)'&'
    ENDDO ! j
    CLOSE(8)

    DEALLOCATE(disp_k_dist,branch)

  END SUBROUTINE generate_disp_curve


  SUBROUTINE calculate_speed_sound
    ! This subroutine calculates the mean speed of sound by evaluating
    ! domega/dk at Gamma and averaging over all directions.  The
    ! directions are uniformly distributed.  This subroutine only works
    ! for monatomic crystals at present.  The polarisation vectors are
    ! calculated for each k, and the dot product of k with the
    ! polarisation vectors are calculated.  The branch with the largest dot
    ! product is the longitudinal branch, whilst the other two branches
    ! are the transverse modes.
    USE rand_no_gen,ONLY : ranx
    IMPLICIT NONE
    REAL(dp) :: kvec(3),kmag,omega(3),cos_theta,sin_theta,phi,c_tr_tot,c_tr, &
      &c2_tr_tot,c2_tr,c_ln_tot,c_ln,c2_ln_tot,c2_ln,err_tr,err_ln,c(3), &
      &kunit(3),pol_vec_real(3,3),dot_prod(3),temp,c_tr_old,c_ln_old
    COMPLEX(dp) :: pol_vec(3,3)
    INTEGER :: i,no_samples,k,k2
    REAL(dp),PARAMETER :: err_tol=1.d-3
    INTEGER,PARAMETER :: max_samples=1000000
    LOGICAL,PARAMETER :: verbose=.FALSE.

    IF(no_atoms_in_prim/=1)CALL errstop('CALCULATE_SPEED_SOUND', &
      &'At the moment this program can only work out the speed of sound in &
      &materials with a single atom per primitive cell.  Sorry about that.')

    ! First guess at a suitable radius of sphere in k-space for computing
    ! derivatives of omega w.r.t. k.
    kmag=0.75d0*small_k_scale

    c_tr_old=HUGE(1.d0) ; c_ln_old=HUGE(1.d0)

    ! Reduce kmag until the calculated sound speeds have converged.

    DO

      c_tr_tot=0.d0  ; c_ln_tot=0.d0
      c2_tr_tot=0.d0 ; c2_ln_tot=0.d0
      no_samples=max_samples

      DO i=1,max_samples

        ! Choose random k vector on sphere of radius kmag.
        cos_theta=1.d0-2.d0*ranx() ; phi=ranx()*twopi
        sin_theta=SQRT(1.d0-cos_theta**2)
        kunit(1:3)=(/sin_theta*COS(phi),sin_theta*SIN(phi),cos_theta/)
        kvec(1:3)=kmag*kunit

        ! Calculate corresponding eigenfrequencies.
        CALL calculate_eigenfreqs_and_vecs(kvec,omega,pol_vec)
        IF(ANY(omega<0.d0))THEN
          WRITE(*,*)'Imaginary frequencies found.'
          WRITE(*,*)'In terms of the primitive reciprocal lattice vectors, &
            &the k-point is:'
          WRITE(*,*)DOT_PRODUCT(kvec,prim_lat_vec(1:3,1)), &
            &DOT_PRODUCT(kvec,prim_lat_vec(1:3,2)), &
            &DOT_PRODUCT(kvec,prim_lat_vec(1:3,3))
          WRITE(*,*)'The frequencies are:'
          WRITE(*,*)omega
          CALL errstop('CALCULATE_SPEED_SOUND','Cannot calculate speed of &
            &sound for unstable lattices.')
        ENDIF ! Soft modes.

        ! Speed of sound corresponding to first three (acoustic) branches.
        c(1:3)=omega(1:3)/kmag

        ! Work out dot products of corresponding polarisation vectors
        ! with unit vector in direction of wave vector.
        pol_vec_real(:,:)=REAL(pol_vec(:,:),dp)
        dot_prod(1)=ABS(DOT_PRODUCT(kunit,pol_vec_real(1:3,1)))
        dot_prod(2)=ABS(DOT_PRODUCT(kunit,pol_vec_real(1:3,2)))
        dot_prod(3)=ABS(DOT_PRODUCT(kunit,pol_vec_real(1:3,3)))

        ! Arrange the three sound speeds in ascending order of dot
        ! product of k with polarisation vector.  The third component
        ! should be the longitudinal one, the first two should be
        ! the transverse components.
        DO k=1,2
          DO k2=k+1,3
            IF(dot_prod(k)>dot_prod(k2))THEN
              temp=dot_prod(k) ; dot_prod(k)=dot_prod(k2) ; dot_prod(k2)=temp
              temp=c(k) ; c(k)=c(k2) ; c(k2)=temp
            ENDIF ! Swap needed
          ENDDO ! k2
        ENDDO ! k

        ! Accumulate sound-speed statistics.
        c_tr_tot=c_tr_tot+c(1)+c(2)
        c_ln_tot=c_ln_tot+c(3)
        c2_tr_tot=c2_tr_tot+c(1)**2+c(2)**2
        c2_ln_tot=c2_ln_tot+c(3)**2

        ! Check whether we have desired accuracy level.
        IF(i>=20)THEN
          c2_ln=c2_ln_tot/DBLE(i)
          c_ln=c_ln_tot/DBLE(i)
          err_ln=SQRT((c2_ln-c_ln**2)/DBLE(i-1))
          IF(err_ln<err_tol*c_ln)THEN
            no_samples=i
            EXIT
          ENDIF
        ENDIF ! i>20

      ENDDO ! i

      ! Mean & standard error in mean for transverse speed.
      c_tr=c_tr_tot/DBLE(2*no_samples)
      c2_tr=c2_tr_tot/DBLE(2*no_samples)
      err_tr=SQRT((c2_tr-c_tr**2)/DBLE(2*no_samples-1))

      ! Mean & standard error in mean for longitudinal speed.
      c_ln=c_ln_tot/DBLE(no_samples)
      c2_ln=c2_ln_tot/DBLE(no_samples)
      err_ln=SQRT((c2_ln-c_ln**2)/DBLE(no_samples-1))

      IF(verbose)THEN
        IF(no_samples==max_samples)WRITE(*,*)'Warning: have not reached &
          &desired error bar.'
        WRITE(*,*)'Radius of k-space sphere : ',kmag
        WRITE(*,*)'Speed of sound (a.u.)'
        WRITE(*,*)'Transverse   : ',c_tr,' +/- ',err_tr
        WRITE(*,*)'Longitudinal : ',c_ln,' +/- ',err_ln
        WRITE(*,*)
      ENDIF ! verbose

      IF(ABS(c_tr-c_tr_old)<2.d0*err_tr.AND.ABS(c_ln-c_ln_old)<2.d0*err_ln)EXIT
      c_tr_old=c_tr
      c_ln_old=c_ln

      kmag=kmag*0.75d0

    ENDDO ! reduce kmag

    IF(verbose)WRITE(*,*)'Final results:'
    WRITE(*,*)'Radius of k-space sphere : ',kmag
    WRITE(*,*)'Please check this is sensible by examining a dispersion curve.'
    WRITE(*,*)

    WRITE(*,*)'Speed of sound (a.u.)'
    WRITE(*,*)'Transverse   : ',c_tr,' +/- ',err_tr
    WRITE(*,*)'Longitudinal : ',c_ln,' +/- ',err_ln
    WRITE(*,*)

  END SUBROUTINE calculate_speed_sound


  ! Evaluate the set of phonon frequencies on the supercell G vectors.
  ! Average the corresponding energies (for testing purposes).
  ! Write out the real part of the non-mass-reduced polarisation vector, which
  ! is the pattern of displacement corresponding to the normal mode.
  subroutine evaluate_freqs_on_grid(kpairs_filename,freq_grids_filename, &
      & disp_patterns_filename,kdisp_patterns_filename,pol_vec_filename, &
      & gvectors_filename,gvectors_frac_filename,error_filename)
    implicit none
    
    ! inputs
    character(*), intent(in) :: kpairs_filename
    character(*), intent(in) :: freq_grids_filename
    character(*), intent(in) :: disp_patterns_filename
    character(*), intent(in) :: kdisp_patterns_filename
    character(*), intent(in) :: pol_vec_filename
    character(*), intent(in) :: gvectors_filename
    character(*), intent(in) :: gvectors_frac_filename
    character(*), intent(in) :: error_filename
    
    
    INTEGER :: ng,i,j,k,ig,ierr,index1,index2,p,n,atom1
    LOGICAL :: found,soft_modes
    REAL(dp) :: gnew(3),gnew1(3),gnew2(3),gvec(3,no_prim_cells),R0(3), &
      &omega(no_DoF_prim),E,F,rec_root_mass(no_atoms_in_prim),GdotR, &
      &disp_pattern(3),kdisp_pattern(3),tot_disp_patt
    COMPLEX(dp) :: pol_vec(no_DoF_prim,no_DoF_prim), &
      &non_mr_pol_vec(3,no_atoms_in_prim),expiGdotR(no_prim_cells), &
      &kpol_vec(3,no_atoms_in_prim)
    REAL(dp),PARAMETER :: tol_omega=1.d-6 ! For judging whether modes are soft.
    REAL(dp),PARAMETER :: tol_g=1.d-8 ! For equivalent +/- G-points.
    INTEGER :: reference(no_prim_cells) ! For equivalent +/- G-points.
    REAL(dp) :: gfrac(3,no_prim_cells),gdir(3) ! Fractional G vectors
    REAL(dp) :: prefactor

    WRITE(*,*)'Frequencies at each supercell G vector:'
    WRITE(*,*)

    ! Evaluate the set of supercell G vectors in the Brillouin zone of the
    ! primitive cell.
    ng=0
    DO k=0,no_prim_cells-1
      gnew2=DBLE(k)*sc_rec_vec(1:3,3)
      DO j=0,no_prim_cells-1
        gnew1=gnew2+DBLE(j)*sc_rec_vec(1:3,2)
        DO i=0,no_prim_cells-1
          gnew=gnew1+DBLE(i)*sc_rec_vec(1:3,1)
          found=.TRUE.
          DO ig=1,ng
            IF(is_lat_point(gnew(1:3)-gvec(1:3,ig),prim_lat_vec))THEN
              found=.FALSE.
              EXIT
            ENDIF ! ig
          ENDDO ! ig
          IF(found)THEN
            ng=ng+1
            IF(ng>no_prim_cells)CALL errstop('EVALUATE_FREQS_ON_GRID', &
              &'Bug: too many G vectors.')
            gvec(1:3,ng)=gnew(1:3)
          ENDIF ! found
        ENDDO ! i
      ENDDO ! j
    ENDDO ! k
    IF(ng/=no_prim_cells)CALL errstop('EVALUATE_FREQS_ON_GRID', &
      &'Bug: too few G vectors.')
    

! Calculate +/- G-vector pairs
! First, write G-vectors as fractions of rec. latt. vecs.
    DO ig=1,no_prim_cells
      gfrac(1,ig)=DOT_PRODUCT(gvec(1:3,ig),prim_lat_vec(1:3,1))
      gfrac(2,ig)=DOT_PRODUCT(gvec(1:3,ig),prim_lat_vec(1:3,2))
      gfrac(3,ig)=DOT_PRODUCT(gvec(1:3,ig),prim_lat_vec(1:3,3))
      DO WHILE(gfrac(1,ig)+0.01*gfrac(1,ig)>0.5d0.AND.ABS(gfrac(1,ig)-0.5d0)>tol_g)
        gfrac(1,ig)=gfrac(1,ig)-1.d0
      ENDDO  
      DO WHILE(gfrac(1,ig)-0.01*gfrac(1,ig)<-0.5d0.AND.ABS(gfrac(1,ig)+0.5d0)>tol_g)
        gfrac(1,ig)=gfrac(1,ig)+1.d0
      ENDDO  
      DO WHILE(gfrac(2,ig)+0.01*gfrac(2,ig)>0.5d0.AND.ABS(gfrac(2,ig)-0.5d0)>tol_g)
        gfrac(2,ig)=gfrac(2,ig)-1.d0
      ENDDO  
      DO WHILE(gfrac(2,ig)-0.01*gfrac(2,ig)<-0.5d0.AND.ABS(gfrac(2,ig)+0.5d0)>tol_g)
        gfrac(2,ig)=gfrac(2,ig)+1.d0
      ENDDO  
      DO WHILE(gfrac(3,ig)+0.01*gfrac(3,ig)>0.5d0.AND.ABS(gfrac(3,ig)-0.5d0)>tol_g)
        gfrac(3,ig)=gfrac(3,ig)-1.d0
      ENDDO  
      DO WHILE(gfrac(3,ig)-0.01*gfrac(3,ig)<-0.5d0.AND.ABS(gfrac(3,ig)+0.5d0)>tol_g)
        gfrac(3,ig)=gfrac(3,ig)+1.d0
      ENDDO  
    ENDDO ! ig
    ! Second, pair them up
    reference=0
    DO k=1,no_prim_cells
      DO j=1,k-1
      IF((ABS(gfrac(1,k)+gfrac(1,j))<tol_g.OR.ABS(gfrac(1,k)+gfrac(1,j)-1)<tol_g.OR.ABS(gfrac(1,k)+gfrac(1,j)+1)<tol_g).AND.&
        &(ABS(gfrac(2,k)+gfrac(2,j))<tol_g.OR.ABS(gfrac(2,k)+gfrac(2,j)-1)<tol_g.OR.ABS(gfrac(2,k)+gfrac(2,j)+1)<tol_g).AND.&
        &(ABS(gfrac(3,k)+gfrac(3,j))<tol_g.OR.ABS(gfrac(3,k)+gfrac(3,j)-1)<tol_g.OR.ABS(gfrac(3,k)+gfrac(3,j)+1)<tol_g))THEN
        reference(k)=j 
        reference(j)=k
        EXIT
      ENDIF
      ENDDO ! j
    ENDDO ! k
    OPEN(1,FILE=kpairs_filename)
    DO i=1,no_prim_cells
      WRITE(1,*)i,reference(i)
    ENDDO ! i
    CLOSE(1)

    OPEN(unit=8,file=freq_grids_filename,status='replace',iostat=ierr)
    IF(ierr/=0)CALL errstop('EVALUATE_FREQS_ON_GRID', &
      &'Problem opening freqs_grid.dat.')
    OPEN(unit=9,file=disp_patterns_filename,status='replace',iostat=ierr)
    OPEN(unit=10,file=kdisp_patterns_filename,status='replace',iostat=ierr)
    OPEN(unit=11,file=pol_vec_filename,status='replace',iostat=ierr)
    IF(ierr/=0)CALL errstop('EVALUATE_FREQS_ON_GRID', &
      &'Problem opening disp_patterns.dat.')

    DO n=1,no_atoms_in_prim
      rec_root_mass(n)=1.d0/SQRT(mass(atom(1,n))) ! 1/sqrt(m) in prim. cell.
    ENDDO ! n

! Modified by B. Monserrat to output G vectors to file
    OPEN(19,FILE=gvectors_filename)
    OPEN(20,FILE=gvectors_frac_filename)
    WRITE(19,*) no_prim_cells
    WRITE(20,*) no_prim_cells

! Evaluate the frequencies at each supercell G vector.
    E=0.d0  ;  F=0.d0
    soft_modes=.FALSE.
    R0=atom_pos(1:3,atom(1,1))
    DO ig=1,no_prim_cells
      WRITE(*,'(" G = (",es20.12,",",es20.12,",",es20.12,")")')twopi*gvec(1:3,ig)
      WRITE(*,'(" G = (",es20.12,",",es20.12,",",es20.12,")")')gfrac(1:3,ig)
      WRITE(19,*)twopi*gvec(1:3,ig)
      ! G-vectors as a fraction of the primitive reciprocal lattice vectors
      WRITE(20,*)ig,gfrac(1,ig),gfrac(2,ig),gfrac(3,ig)
      IF(reference(ig)==0)THEN
        CALL calculate_eigenfreqs_and_vecs(twopi*gvec(1:3,ig),omega,pol_vec)
      ELSE
        IF(reference(ig)>ig)THEN
          CALL calculate_eigenfreqs_and_vecs(twopi*gvec(1:3,ig),omega,pol_vec)
        ELSE
          CALL calculate_eigenfreqs_and_vecs(twopi*gvec(1:3,reference(ig)),omega,pol_vec)
        ENDIF
      ENDIF

! The negative is used because the matrix of force constants is the transpose of
! the usual expression in derivations that lead to a positive exponential
      DO p=1,no_prim_cells
        IF(reference(ig)==0)THEN
          GdotR=-DOT_PRODUCT(twopi*gvec(1:3,ig),atom_pos(1:3,atom(p,1))-R0)
        ELSE
          IF(reference(ig)>ig)THEN
            GdotR=-DOT_PRODUCT(twopi*gvec(1:3,ig),atom_pos(1:3,atom(p,1))-R0)
          ELSE
            GdotR=-DOT_PRODUCT(twopi*gvec(1:3,reference(ig)),atom_pos(1:3,atom(p,1))-R0)
          ENDIF
        ENDIF
        expiGdotR(p)=CMPLX(COS(GdotR),SIN(GdotR),dp) ! Store exp(iG.R_p).
      ENDDO ! p

      DO index2=1,no_DoF_prim
        WRITE(*,*)'  omega = ',omega(index2)
        IF(omega(index2)>tol_omega)THEN
! Ignore contributions from imaginary or zero frequencies.
          E=E+harmonic_energy(temperature,omega(index2))
          F=F+harmonic_free_energy(temperature,omega(index2))
        ELSEIF(omega(index2)<-tol_omega)THEN
          soft_modes=.TRUE.
        ENDIF ! omega>0
        WRITE(8,*)omega(index2),1.d0
! Compute the non-mass-reduced polarisation vector.
        index1=0
        DO n=1,no_atoms_in_prim
          DO i=1,3
            index1=index1+1
            non_mr_pol_vec(i,n)=pol_vec(index1,index2)*rec_root_mass(n)
            kpol_vec(i,n)=pol_vec(index1,index2)
          ENDDO ! i
        ENDDO ! n
        IF(omega(index2)<-tol_omega)THEN
          WRITE(9,*)'Frequency           : ',omega(index2),' (SOFT)'
          WRITE(10,*)'Frequency           : ',omega(index2),' (SOFT)'
          WRITE(11,*)'Mode number     :',index2, '   Frequency           : ',omega(index2),' (SOFT)'
        ELSE
          WRITE(9,*)'Frequency           : ',omega(index2)
          WRITE(10,*)'Frequency           : ',omega(index2)
          WRITE(11,*)'Mode number     :',index2, '   Frequency           : ',omega(index2)
        ENDIF ! soft freq.
        WRITE(9,*)gvec(1:3,ig)*twopi
        WRITE(10,*)gvec(1:3,ig)*twopi
        WRITE(11,*)gvec(1:3,ig)*twopi
        WRITE(9,*)'Displacement pattern for each atom:'
        WRITE(10,*)'Displacement pattern for each atom:'
        WRITE(11,*)'Polarisation vector:'
        disp_pattern=0.d0
        kdisp_pattern=0.d0
        tot_disp_patt=0.d0
        DO atom1=1,no_atoms_in_sc
! Displacement pattern: polarisation vector times exp(iG.R).
          IF(reference(ig)==0)THEN
            disp_pattern=REAL(non_mr_pol_vec(1:3,atom_in_prim(atom1)) & ! Note only the real part is taken
              &*expiGdotR(prim_cell_for_atom(atom1)),dp)
            tot_disp_patt=tot_disp_patt+sqrt(disp_pattern(1)**2+disp_pattern(2)**2+disp_pattern(3)**2)
            kdisp_pattern=kpol_vec(1:3,atom_in_prim(atom1))! & ! Note only the real part is taken
              !&*expiGdotR(prim_cell_for_atom(atom1))
            prefactor=1.d0
          ELSE
            IF(reference(ig)>ig)THEN
              disp_pattern=REAL(non_mr_pol_vec(1:3,atom_in_prim(atom1)) &
               &*expiGdotR(prim_cell_for_atom(atom1)),dp)
              tot_disp_patt=tot_disp_patt+sqrt(disp_pattern(1)**2+disp_pattern(2)**2+disp_pattern(3)**2)
              kdisp_pattern=REAL(kpol_vec(1:3,atom_in_prim(atom1)) &
               &*expiGdotR(prim_cell_for_atom(atom1)),dp)
              prefactor=SQRT(2.d0)
            ELSE
              disp_pattern=AIMAG(non_mr_pol_vec(1:3,atom_in_prim(atom1)) &
               &*expiGdotR(prim_cell_for_atom(atom1)))
              tot_disp_patt=tot_disp_patt+sqrt(disp_pattern(1)**2+disp_pattern(2)**2+disp_pattern(3)**2)
              kdisp_pattern=AIMAG(kpol_vec(1:3,atom_in_prim(atom1)) &
               &*expiGdotR(prim_cell_for_atom(atom1)))
              prefactor=SQRT(2.d0)
            ENDIF
          ENDIF
          WRITE(9,*)disp_pattern,prefactor
          WRITE(10,*)kdisp_pattern,prefactor
          WRITE(11,*)REAL(non_mr_pol_vec(1:3,atom_in_prim(atom1)))
          WRITE(11,*)AIMAG(non_mr_pol_vec(1:3,atom_in_prim(atom1)))
        ENDDO ! atom1
        WRITE(9,*)
        WRITE(10,*)
        WRITE(11,*)
      ENDDO ! index2
      WRITE(*,*)
      IF(tot_disp_patt<1.d-8)THEN
        OPEN(111,file=error_filename)
        WRITE(111,*)'The total displacement is:',tot_disp_patt
        CLOSE(111)
      ENDIF
    ENDDO ! ig
    E=E/DBLE(no_prim_cells)  ;  F=F/DBLE(no_prim_cells)

    CLOSE(8)
    CLOSE(9)
    CLOSE(10)
    CLOSE(11)
    CLOSE(19)
    CLOSE(20)

    WRITE(*,*)'Mean LTE per primitive cell  : ',E
    WRITE(*,*)'Mean LTE per primitive cell (eV)  : ',E*27.211396132d0
    WRITE(*,*)'Mean LTFE per primitive cell : ',F
    WRITE(*,*)'Mean LTFE per primitive cell (eV): ',F*27.211396132d0
    IF(soft_modes)WRITE(*,*)'WARNING: soft modes are present.'
    WRITE(*,*)'Frequencies written to freqs_grid.dat.'
    WRITE(*,*)'Displacement patterns written to disp_patterns.dat.'
    WRITE(*,*)

  END SUBROUTINE evaluate_freqs_on_grid


  ! Write out the dynamical matrix at each supercell G-vector to a file
  subroutine write_dynamical_matrix(dyn_mat_fileroot)
    implicit none
    
    ! inputs
    character(*), intent(in) :: dyn_mat_fileroot
    
    INTEGER :: ierr,ng,k,j,i,ig,atom1,cart1,index1,atom2,cart2,index2
    REAL(dp) :: gnew2(3),gnew1(3),gnew(3),gvec(3,no_prim_cells)
    COMPLEX(dp) :: dyn_mat(no_DoF_prim,no_DoF_prim)
    LOGICAL :: found

    ng=0
    DO k=0,no_prim_cells-1
      gnew2=DBLE(k)*sc_rec_vec(1:3,3)
      DO j=0,no_prim_cells-1
        gnew1=gnew2+DBLE(j)*sc_rec_vec(1:3,2)
        DO i=0,no_prim_cells-1
          gnew=gnew1+DBLE(i)*sc_rec_vec(1:3,1)
          found=.TRUE.
          DO ig=1,ng
            IF(is_lat_point(gnew(1:3)-gvec(1:3,ig),prim_lat_vec))THEN
              found=.FALSE.
              EXIT
            ENDIF ! ig
          ENDDO ! ig
          IF(found)THEN
            ng=ng+1
            IF(ng>no_prim_cells)CALL errstop('WRITE_DYNAMICAL_MATRIX','Too &
             &many G-vectors.')
            gvec(1:3,ng)=gnew(1:3)
          ENDIF ! found
        ENDDO ! i
      ENDDO ! j
    ENDDO ! k
    IF(ng/=no_prim_cells)CALL errstop('WRITE_DYNAMICAL_MATRIX','Too few &
     &G-vectors.')

    dyn_mat(1:no_DoF_prim,1:no_DoF_prim)=CMPLX(0.d0,0.d0,dp)

    DO ig=1,ng
      OPEN(unit=101,file=dyn_mat_fileroot//TRIM(i2s(ig))//'.dat',status='replace',iostat=ierr)
      IF(ierr/=0)CALL errstop('WRITE_DYNAMICAL_MATRIX','Problem opening &
       &dyn_mat.'//TRIM(i2s(ig))//'.dat file.')
      CALL construct_dyn_matrix(twopi*gvec(1:3,ig),dyn_mat)
      atom1=0
      DO index1=1,no_DoF_prim
        atom2=0
        IF(MOD(index1,3)==1)THEN
          atom1=atom1+1
          cart1=1
        ENDIF ! MOD(index1,3)==1
        DO index2=1,no_DoF_prim
          IF(MOD(index2,3)==1)THEN
            atom2=atom2+1
            cart2=1
          ENDIF ! MOD(index2,3)==1
          WRITE(101,*)atom1,cart1,atom2,cart2,REAL(dyn_mat(index1,index2)),AIMAG(dyn_mat(index1,index2))
          cart2=cart2+1
        ENDDO ! index2
        cart1=cart1+1
      ENDDO ! index1
      CLOSE(101)
    END DO ! ig

  END SUBROUTINE write_dynamical_matrix


  ! Write out atoms in primitive cell in order.
  subroutine write_atoms_in_primitive_cell(atoms_in_primitive_cell_filename)
    implicit none
    
    ! inputs
    character(*), intent(in) :: atoms_in_primitive_cell_filename
    
    REAL(dp),PARAMETER :: tol=1.d-10
    INTEGER :: ierr,n,atom1,i
    REAL(dp) :: pos(3),frac(3)

    OPEN(unit=102,file=atoms_in_primitive_cell_filename,status='replace',&
     &iostat=ierr)
    IF(ierr/=0)CALL errstop('WRITE_ATOMS_IN_PRIMITIVE_CELL','Problem opening &
     &atoms_in_primitive_cell.dat file.')

    DO n=1,no_atoms_in_prim
      atom1=atom(1,n)
      pos=atom_pos(1:3,atom1)
      DO i=1,3
        frac(i)=DOT_PRODUCT(pos(1:3),prim_rec_vec(1:3,i))
      ENDDO ! i
      frac(1:3)=MODULO(frac(1:3)+tol,1.d0)-tol
      WRITE(102,*)mass(atom1),frac(1:3)
    ENDDO ! n  

    CLOSE(102)

  END SUBROUTINE write_atoms_in_primitive_cell



  SUBROUTINE finalise
    ! Deallocate arrays.
    IMPLICIT NONE
    IF(ALLOCATED(species))DEALLOCATE(species,mass,atom_pos,force_const, &
      &defined,atom,no_equiv_ims,delta_prim)
    IF(ALLOCATED(rotation))DEALLOCATE(rotation,offset)
    IF(ALLOCATED(disp_kpoints))DEALLOCATE(disp_kpoints)
    IF(ALLOCATED(freq_dos))DEALLOCATE(freq_dos)
    IF(ALLOCATED(atom_in_prim))DEALLOCATE(atom_in_prim,prim_cell_for_atom)
  END SUBROUTINE finalise


END MODULE phonons


module lte_module
  implicit none
contains

! ------------------------------------------------------------
! Main program starts here.
! ------------------------------------------------------------
subroutine lte(tol,tol2,delta,lte_filename,freq_dos_filename,                &
    & tdependence1_filename,tdependence2_filename,dispersion_curve_filename, &
    & kpairs_filename,freq_grids_filename,disp_patterns_filename,            &
    & kdisp_patterns_filename,pol_vec_filename,gvectors_filename,            &
    & gvectors_frac_filename,error_filename,dyn_mat_fileroot,                &
    & atoms_in_primitive_cell_filename)
  use constants, only : dp
  use utils,     only : errstop, wordwrap
  use phonons
  IMPLICIT NONE
  
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
  REAL :: t1,t2

  CALL CPU_TIME(t1)

  WRITE(*,*)
  WRITE(*,*)'LATTICE THERMAL ENERGY'
  WRITE(*,*)'======================'
  WRITE(*,*)

  WRITE(*,*)'Reading data from lte.dat...'
  WRITE(*,*)
  CALL read_lte(tol, lte_filename)
  WRITE(*,*)'Finished reading input data.'
  WRITE(*,*)

  WRITE(*,*)'Applying point symmetries to the matrix of force constants...'
  CALL point_symm(tol)
  WRITE(*,*)'Done.'
  WRITE(*,*)

  IF(ANY(.NOT.defined))THEN
    CALL wordwrap('WARNING: will impose symmetries on the matrix of force &
      &constants iteratively...')
    CALL point_symm_brute_force(tol)
    WRITE(*,*)'Done.'
    WRITE(*,*)
  ENDIF
  IF(ANY(.NOT.defined))CALL errstop('LTE','Some elements of the matrix of &
    &force constants are still undefined.')

  CALL wordwrap('Imposing Newton''s third law and symmetry on the matrix of &
    &force constants...')
  CALL newtons_law(tol2)
  WRITE(*,*)'Done.'
  WRITE(*,*)

  WRITE(*,*)'Performing mass reduction on the matrix of force constants...'
  CALL mass_reduce
  WRITE(*,*)'Done.'
  WRITE(*,*)

  WRITE(*,*)'Establishing the primitive lattice vector associated with &
    &each atom...'
  CALL find_prim_cell(delta)
  WRITE(*,*)'Done.'
  WRITE(*,*)

  IF(prog_function==1)THEN

    WRITE(*,*)'Calculating the frequency density-of-states function...'
    CALL calculate_freq_dos(tol, freq_dos_filename)
    CALL wordwrap('Done.  Frequency density-of-states function written to &
      &freq_dos.dat.  (Please view this file using XMGrace.)')
    WRITE(*,*)

    WRITE(*,*)'Calculating the lattice thermal energy (LTE) and free energy &
      &(LTFE)...'
    CALL calc_lte(tdependence1_filename)
    CALL calc_ltfe(tdependence2_filename)
    WRITE(*,*)

  ELSEIF(prog_function==2)THEN

    WRITE(*,*)'Calculating the requested dispersion curve.'
    CALL generate_disp_curve(dispersion_curve_filename)
    CALL wordwrap('Done.  dispersion_curve.dat has been generated.  (Please &
      &view this file using XMGrace.)')
    WRITE(*,*)

  ELSEIF(prog_function==3)THEN

    WRITE(*,*)'Calculating the speed of sound.'
    CALL calculate_speed_sound
    WRITE(*,*)'Done.  Speed of sound calculated.'
    WRITE(*,*)

  ELSEIF(prog_function==4)THEN

    WRITE(*,*)'Calculating the frequencies and displacement patterns on the &
      &G-vector grid.'
    CALL evaluate_freqs_on_grid(kpairs_filename,freq_grids_filename,     &
      & disp_patterns_filename,kdisp_patterns_filename,pol_vec_filename, &
      & gvectors_filename,gvectors_frac_filename,error_filename)
    WRITE(*,*)'Done.  Frequencies and displacement patterns calculated.'
    WRITE(*,*)
    CALL write_dynamical_matrix(dyn_mat_fileroot)
    CALL write_atoms_in_primitive_cell(atoms_in_primitive_cell_filename)

  ELSE

    CALL errstop('LTE','Program function should be 1, 2, 3 or 4.')

  ENDIF ! prog_function.

  CALL finalise

  CALL CPU_TIME(t2)

  WRITE(*,*)'Program finished.  Time taken: ',t2-t1
  WRITE(*,*)
end subroutine
end module
