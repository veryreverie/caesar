! Lattice Thermal Energy, Neil Drummond, 12/2004-1/2005.

! All modifications are labelled as: Modified by B. Monserrat

! This program enables the user to calculate the lattice thermal energy and
! free energy of a crystal lattice within the quasiharmonic approximation.
! It also enables the calculation of phonon dispersion curves.  Finally it
! can be used to calculate the longitudinal and transverse sound speeds for a
! monatomic crystal by averaging the gradients of the acoustic branches
! over all directions.

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
  use constants,      only : dp, kB_au_per_K
  use utils,          only : i2s, errstop
  use linear_algebra, only : determinant33, inv_33
  use min_images,     only : min_images_brute_force, maxim
  use file_module
  implicit none
  
  private
  public :: lte_1,lte_2,lte_3,lte_4
contains

! ----------------------------------------------------------------------
! This function returns the number of the atom at rvec (up to translations
! through a supercell lattice vector.
! ----------------------------------------------------------------------
integer function atom_at_pos(rvec,structure_sc)
  use structure_module
  implicit none
  
  real(dp),            intent(in) :: rvec(3)
  type(StructureData), intent(in) :: structure_sc
  
  integer :: atom1
  
  atom_at_pos=0
  do atom1=1,structure_sc%no_atoms
    ! Separation of rvec and atom posn. Is a sc lattice vector.
    if (is_lat_point( structure_sc%cart_atoms(:,atom1)-rvec(:), &
                    & structure_sc%recip_lattice)) then
      atom_at_pos = atom1
      exit
    endif
  enddo
end function

! ----------------------------------------------------------------------
! This function returns T if and only if rvec is a lattice vector.
! rec_vec holds the reciprocal lattice vectors.
! ----------------------------------------------------------------------
function is_lat_point(rvec,rec_vec) result(output)
  implicit none
  
  real(dp), intent(in) :: rvec(3)
  real(dp), intent(in) :: rec_vec(3,3)
  logical              :: output
  
  real(dp), parameter :: tol=1.d-3
  
  real(dp) :: t(3)
  
  t = matmul(rvec,rec_vec)
  output = all(nint(t)-t<tol)
end function is_lat_point

! ----------------------------------------------------------------------
! Read in the data in lte.dat and allocate arrays, etc.
! ----------------------------------------------------------------------
subroutine read_lte(prog_function,tol,structure,structure_sc,atoms, &
   & displacements,forces,temperature,no_kspace_lines,disp_kpoints, &
   & fc_scale,      &
   & no_prim_cells,length_scale, &
   & vol_scale,                                                     &
   & force_const,defined,atom)
  use constants, only : pi
  use string_module
  use structure_module
  implicit none
  
  integer,             intent(in) :: prog_function
  real(dp),            intent(in) :: tol
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: structure_sc
  integer,             intent(in) :: atoms(:)
  integer,             intent(in) :: displacements(:)
  real(dp),            intent(in) :: forces(:,:,:)
  real(dp),            intent(in) :: temperature
  integer,             intent(in) :: no_kspace_lines
  real(dp),            intent(in) :: disp_kpoints(:,:)
  
  real(dp),              intent(out) :: fc_scale
  integer,               intent(out) :: no_prim_cells
  real(dp),              intent(out) :: length_scale
  real(dp),              intent(out) :: vol_scale
  real(dp), allocatable, intent(out) :: force_const(:,:,:,:)
  logical,  allocatable, intent(out) :: defined(:,:,:,:)
  integer,  allocatable, intent(out) :: atom(:,:)
  
  integer  :: n,i,j,k
  integer  :: atom1,atom2,dir1,dir2
  real(dp) :: check_matrix
  
  ! Define the "length scale" to be the average of the lengths of the
  ! primitive lattice vectors and the reciprocal length scale to be the
  ! inverse of this.
  length_scale = 0.d0
  do i=1,3
    length_scale = length_scale + norm2(structure%lattice(:,i))/3.0_dp
  enddo
  vol_scale=length_scale**3

  ! Number of unit cells.
  no_prim_cells = structure_sc%no_atoms/structure%no_atoms
  
  allocate( force_const(structure_sc%no_atoms,3,structure_sc%no_atoms,3), &
          & defined(structure_sc%no_atoms,3,structure_sc%no_atoms,3),     &
          & atom(no_prim_cells,structure%no_atoms))
  
  ! Read in atom positions, species and masses.
  
  fc_scale = 0.d0
  defined(:,:,:,:)=.FALSE.
  force_const(:,:,:,:)=0.d0
  do i=1,size(atoms)
    atom1 = atoms(i)
    dir1 = displacements(i)
    do j=1,structure_sc%no_atoms
      atom2 = j
      do k=1,3
        dir2 = k
        call trans_symm(atom1,dir1,atom2,dir2,forces(k,j,i), &
           & tol,structure_sc, &
           & no_prim_cells,structure,defined,force_const)
      enddo
    enddo
  enddo
  fc_scale=fc_scale/(size(atoms)*3*structure_sc%no_atoms)
  
  ! ----------------------------------------------------------------------
  ! Copy data to terminal
  ! ----------------------------------------------------------------------
  
  ! Primitive lattice vectors
  write(*,*) 'Primitive lattice vectors (Cartesian components in rows, a.u.):'
  write(*,*) structure%lattice(1:3,1)
  write(*,*) structure%lattice(1:3,2)
  write(*,*) structure%lattice(1:3,3)
  write(*,*)
  
  ! Supercell lattice vectors
  write(*,*) 'Supercell lattice vectors (Cartesian components in rows, a.u.):'
  write(*,*) structure_sc%lattice(1:3,1)
  write(*,*) structure_sc%lattice(1:3,2)
  write(*,*) structure_sc%lattice(1:3,3)
  write(*,*)
  
  ! Reciprocal lattice vectors, etc.
  write(*,*) 'Number of primitive unit cells     : '//TRIM(i2s(no_prim_cells))
  write(*,*)
  write(*,*) 'Prim. rec. latt. vectors (Cart. cmpnts &
    &in rows, factor of 2pi inc., a.u.):'
  write(*,*) structure%recip_lattice(1:3,1)*2*pi
  write(*,*) structure%recip_lattice(1:3,2)*2*pi
  write(*,*) structure%recip_lattice(1:3,3)*2*pi
  write(*,*)

  write(*,*) 'Supercell rec. latt. vectors(Cart. cmpnts &
   &in rows, factor of 2pi inc., a.u.):'
  write(*,*) structure_sc%recip_lattice(1:3,1)*2*pi
  write(*,*) structure_sc%recip_lattice(1:3,2)*2*pi
  write(*,*) structure_sc%recip_lattice(1:3,3)*2*pi
  write(*,*)
  
  if(structure%volume<tol*vol_scale)call errstop('READ_LTE', &
    &'Primitive lattice vectors should be linearly independent.  Please &
    &check your lattice vectors.')

  if(structure_sc%volume<tol*vol_scale)call errstop('READ_LTE', &
    &'Supercell lattice vectors should be linearly independent.  Please &
    &check your lattice vectors.')
  
  if(ABS(DBLE(no_prim_cells)*structure%volume-structure_sc%volume) &
    &>tol*vol_scale)call errstop('READ_LTE','Supercell volume should &
    &be an integer multiple of primitive-cell volume.  Please check your &
    &lattice vectors.')
  
  ! Number of atoms.
  write(*,*)'Number of atoms in supercell       : '//TRIM(i2s(structure_sc%no_atoms))
  if(structure_sc%no_atoms<no_prim_cells)call errstop('READ_LTE', &
    &'Need more atoms in the supercell!')
  if(MOD(structure_sc%no_atoms,no_prim_cells)/=0)call errstop('READ_LTE', &
    &'Number of atoms in supercell should be a multiple of the number of &
    &primitive cells.')
  write(*,*)
  
  ! Atom positions, species and masses.
  write(*,*)'Species ; Mass (a.u.) ; Position (Cartesian coordinates, a.u.)'
  do i=1,structure_sc%no_atoms
    write(*,'(" ",a,"  ",f14.6," ",3("  ",f14.6))')structure_sc%species(i),structure_sc%mass(i), &
      &structure_sc%cart_atoms(1:3,i)
    if(structure_sc%mass(i)<=0.d0)call errstop('READ_LTE','Mass should be positive.')
  enddo
  
  ! Check atoms aren't on top of each other.
  do i=1,structure_sc%no_atoms-1
    do j=i+1,structure_sc%no_atoms
      if(is_lat_point(structure_sc%cart_atoms(1:3,j)-structure_sc%cart_atoms(1:3,i),structure_sc%recip_lattice))call &
        &errstop('READ_LTE','Atoms '//TRIM(i2s(i))//' and '//TRIM(i2s(j)) &
        &//' appear to be on top of one another.')
    enddo ! j
  enddo ! i
  write(*,*)
  
  ! Point-symmetry operations.
  write(*,*)'Number of point symmetries         : '//TRIM(i2s(structure_sc%no_symmetries))
  if(structure_sc%no_symmetries<1)call errstop('READ_LTE','At least one point-symmetry &
    &structure_sc%rotation_matrices matrix (identity) must be supplied.')
  do n=1,structure_sc%no_symmetries
    do i=1,3
      do j=i,3
        check_matrix=DOT_PRODUCT(structure_sc%rotation_matrices(1:3,i,n),structure_sc%rotation_matrices(1:3,j,n))
        if((i==j.AND.ABS(check_matrix-1.d0)>tol) &
          &.OR.(i/=j.AND.ABS(check_matrix)>tol))call &
          &errstop('READ_LTE','Rotation matrix '//TRIM(i2s(n)) &
          &//' is not orthogonal!')
      enddo
    enddo
  enddo
  write(*,*)'Have read in structure_sc%rotation_matrices matrices and translation vectors.'
  write(*,*)

  ! Force constants supplied.
  write(*,*)'Number of force constants supplied : ' &
    &//TRIM(i2s(size(atoms)*3*structure_sc%no_atoms))
  if(size(atoms)*3*structure_sc%no_atoms<=0)call errstop('READ_LTE', &
    &'Need to supply more force data!')
  write(*,*)fc_scale
  write(*,*)'Have read in the force-constant data and applied &
    &translational symmetry.'
  write(*,*)

  if(prog_function==1)then
    write(*,*)'Temperature (K)                    :',temperature
    if(temperature<0.d0)call errstop('READ_LTE', &
      &'Temperature should be non-negative.')
    if(temperature<=0.d0)write(*,*)'(i.e. the zero-point energy is to be &
      &calculated.)'
    write(*,*)
  endif
  
  if(prog_function==2)then
    write(*,*)'Number of lines in k-space to plot     : ' &
      &//TRIM(i2s(no_kspace_lines))
    if(no_kspace_lines<1)call errstop('READ_LTE', &
      &'Need to supply more lines in k-space!')
    write(*,*)'Points along walk in reciprocal space &
      &(Cartesian components in a.u.):'
    do i=0,no_kspace_lines
      write(*,'(3(" ",f16.8))')disp_kpoints(1:3,i)
    enddo
    write(*,*)'Have read in points for dispersion curve.'
  endif
end subroutine

! ----------------------------------------------------------------------
! Take the force constant fc supplied for (atom1,dir1,atom2,dir2),
! place it in all translationally equivalent elements of the matrix of
! force constants.
! ----------------------------------------------------------------------
subroutine trans_symm(atom1,dir1,atom2,dir2,fc,tol,structure_sc, &
   & no_prim_cells,structure,defined,force_const)
  use structure_module
  implicit none
  
  real(dp), intent(in) :: tol
  integer , intent(in) :: atom1,dir1,atom2,dir2
  real(dp), intent(in) :: fc
  type(StructureData),  intent(in) :: structure_sc
  integer,  intent(in) :: no_prim_cells
  type(StructureData), intent(in) :: structure
  
  logical,  intent(inout) :: defined(:,:,:,:)
  real(dp), intent(inout) :: force_const(:,:,:,:)
  
  integer :: atom1p,atom2p,no_translations
  real(dp) :: pos_atom2p(3),relpos_atom2_atom1(3)
  
  relpos_atom2_atom1(1:3)=structure_sc%cart_atoms(1:3,atom2)-structure_sc%cart_atoms(1:3,atom1)
  no_translations=0
  do atom1p=1,structure_sc%no_atoms
    if(is_lat_point(structure_sc%cart_atoms(1:3,atom1p)-structure_sc%cart_atoms(1:3,atom1), &
      &structure%recip_lattice))then
      ! atom1p and atom1 are equivalent under trans. symm.
      if (abs(structure_sc%mass(atom1p)/structure_sc%mass(atom1)-1)>tol) then
        call errstop('TRANS_SYMM','Atoms '//TRIM(i2s(atom1))//' and ' &
        &//TRIM(i2s(atom1p))//' are equivalent by translational symmetry, &
        &but they have different masses.')
      endif
      pos_atom2p(1:3)=structure_sc%cart_atoms(1:3,atom1p)+relpos_atom2_atom1(1:3)
      atom2p=atom_at_pos(pos_atom2p,structure_sc)
      if(atom2p<=0)call errstop('TRANS_SYMM','Please check that your atom &
        &coordinates satisfy the translational symmetry they should have.')
      ! atom2p and atom2 are related to each other by the same translation
      ! as atom1p and atom1.
      if (abs(structure_sc%mass(atom2p)/structure_sc%mass(atom2)-1)>tol) then
        call errstop('TRANS_SYMM','Atoms '//TRIM(i2s(atom2))//' and ' &
        &//TRIM(i2s(atom2p))//' are equivalent by translational &
        &symmetry, but they have different masses.')
      endif
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
! symmetry translation and structure_sc%rotation_matrices have been applied.  The matrix of
! force constants Phi then transforms as Phi' = R Phi R^T, where R is the
! structure_sc%rotation_matrices matrix and Phi and Phi' are the matrices of force constants
! between atom1 & atom2 and atom1p & atom2p, respectively.  Some elements
! of Phi may be unknown, but so long as they are multiplied by zero we can
! still work out the corresponding elements of Phi'.
! ----------------------------------------------------------------------
subroutine point_symm(tol,structure_sc,defined,force_const)
  use linear_algebra, only : ddot
  use structure_module
  implicit none
  
  real(dp),            intent(in) :: tol
  type(StructureData), intent(in) :: structure_sc
  
  logical,  intent(inout) :: defined(:,:,:,:)
  real(dp), intent(inout) :: force_const(:,:,:,:)
  
  integer :: atom1,atom1p,atom2,atom2p,i,j,n,ip,jp,&
    &weight(structure_sc%no_atoms,3,structure_sc%no_atoms,3)
  real(dp) :: fc,product,pos_atom1p(3),pos_atom2p(3)
  logical :: well_defined

  ! We average over all the force constants that ought, by the full
  ! symmetry of the supercell, to be identical.  The weight array is
  ! used to perform this average.
  do j=1,3
    do atom2=1,structure_sc%no_atoms
      do i=1,3
        do atom1=1,structure_sc%no_atoms
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

  do n=1,structure_sc%no_symmetries

    do atom1=1,structure_sc%no_atoms

      ! Rotate atom coordinates and identify equivalent atom.
      do i=1,3
        ! LAPACK commented out because it isn't working. 9/1/2017
        ! pos_atom1p(i) = structure_sc%offsets_cart(i,n) &
        !             & + ddot(3,structure_sc%rotation_matrices(i,1,n),3,structure_sc%cart_atoms(1,atom1),1)
        pos_atom1p(i) = structure_sc%offsets_cart(i,n) &
                    & + dot_product(structure_sc%rotation_matrices(i,:,n),structure_sc%cart_atoms(:,atom1))
      enddo ! i
      atom1p=atom_at_pos(pos_atom1p,structure_sc)
      !write(*,*)n,atom1,atom1p,atom_at_pos(pos_atom1p,structure_sc%cart_atoms)
      if(atom1p<=0)call errstop('POINT_SYMM','Please check that &
        &your atom coordinates satisfy the rotational symmetries that you &
        &have supplied.  NB, I have assumed that r''=b+Rr, where R is the &
        &structure_sc%rotation_matrices matrix and b is the translation.  This may be wrong.')
      if(ABS(structure_sc%mass(atom1)-structure_sc%mass(atom1p))>tol*structure_sc%mass(atom1))call &
        &errstop('POINT_SYMM','Atoms '//TRIM(i2s(atom1))//' and ' &
        &//TRIM(i2s(atom1p))//' are equivalent by rotational symmetry, &
        &but they have different masses.')

      do atom2=1,structure_sc%no_atoms

        ! Rotate atom coordinates and identify equivalent atom.
        do i=1,3
          ! LAPACK commented out because it isn't working. 9/1/2017
          ! pos_atom2p(i)=structure_sc%offsets_cart(i,n)+ddot(3,structure_sc%rotation_matrices(i,1,n),3, &
          !   &structure_sc%cart_atoms(1,atom2),1)
          pos_atom2p(i) = structure_sc%offsets_cart(i,n) &
                      & + dot_product(structure_sc%rotation_matrices(i,:,n),structure_sc%cart_atoms(:,atom2))
        enddo ! i
        atom2p=atom_at_pos(pos_atom2p,structure_sc)
        !write(*,*)n,atom2,atom2p,atom_at_pos(pos_atom2p,structure_sc%cart_atoms),pos_atom2p(1),pos_atom2p(2),pos_atom2p(3)
        if(atom2p<=0)call errstop('POINT_SYMM','Please check that &
          &your atom coordinates satisfy the rotational symmetries that &
          &you have supplied.  NB, I have assumed that r''=b+Rr, where R &
          &is the structure_sc%rotation_matrices matrix and b is the translation.  This may be &
          &wrong.')

        ! Apply structure_sc%rotation_matrices to maxtrix of force constants.  Record whether or
        ! not each force constant is well-defined.
        do i=1,3
          do j=1,3
            fc=0.d0
            well_defined=.TRUE.
            ip_loop : do ip=1,3
              if(ABS(structure_sc%rotation_matrices(i,ip,n))>tol)then
                product=0.d0
                do jp=1,3
                  if(ABS(structure_sc%rotation_matrices(j,jp,n))>tol)then
                    if(defined(atom1,ip,atom2,jp))then
                      product=product+structure_sc%rotation_matrices(j,jp,n) &
                        &*force_const(atom1,ip,atom2,jp)
                    else
                      well_defined=.FALSE.
                      exit ip_loop
                    endif ! Force constant defined
                  endif ! Rotation non-zero.
                enddo ! jp
                fc=fc+product*structure_sc%rotation_matrices(i,ip,n)
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

!    do atom1=1,structure_sc%no_atoms
!      do atom2=1,structure_sc%no_atoms
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
subroutine point_symm_brute_force(tol,fc_scale,structure_sc, &
   & structure,   &
   & force_const,defined)
  use linear_algebra,only : ddot
  use structure_module
  implicit none
  
  real(dp), intent(in) :: tol
  real(dp), intent(in) :: fc_scale
  type(StructureData),  intent(in) :: structure_sc
  type(StructureData),  intent(in) :: structure
  
  real(dp), intent(inout) :: force_const(:,:,:,:)
  logical,  intent(inout) :: defined(:,:,:,:)
  
  integer,parameter :: max_t=10000,min_t=10 ! Max & min number of impositions.
  
  integer :: atom1,atom1p,atom2,atom2p,i,j,n,ip,t,three_noDoFprim_sq
  real(dp) :: fc,pos_atom1p(3),pos_atom2p(3),max_diff
  logical :: last_iteration

  three_noDoFprim_sq=3*structure%no_modes**2

  last_iteration=.FALSE.

  do t=1,max_t

    max_diff=0.d0

    do n=1,structure_sc%no_symmetries

      do atom1=1,structure_sc%no_atoms

        ! Rotate atom coordinates and identify equivalent atom.
        do i=1,3
          ! LAPACK commented out because it isn't working. 9/1/2017
          ! pos_atom1p(i)=structure_sc%offsets_cart(i,n)+ddot(3,structure_sc%rotation_matrices(i,1,n),3, &
          !   &structure_sc%cart_atoms(1,atom1),1)
          pos_atom1p(i) = structure_sc%offsets_cart(i,n) &
                      & + dot_product(structure_sc%rotation_matrices(i,:,n),structure_sc%cart_atoms(:,atom1))
        enddo ! i
        atom1p=atom_at_pos(pos_atom1p,structure_sc)
        if(atom1p<=0)call errstop('POINT_SYMM_BRUTE_FORCE','Please check &
          &that your atom coordinates satisfy the rotational symmetries that &
          &you have supplied.  NB, I have assumed that r''=b+Rr, where R is &
          &the structure_sc%rotation_matrices matrix and b is the translation.  This may be &
          &wrong.')

        do atom2=1,structure_sc%no_atoms

          ! Rotate atom coordinates and identify equivalent atom.
          do i=1,3
            ! LAPACK commented out because it isn't working. 9/1/2017
            ! pos_atom2p(i)=structure_sc%offsets_cart(i,n)+ddot(3,structure_sc%rotation_matrices(i,1,n),3, &
            !   &structure_sc%cart_atoms(1,atom2),1)
            pos_atom2p(i) = structure_sc%offsets_cart(i,n) &
                        & + dot_product(structure_sc%rotation_matrices(i,:,n),structure_sc%cart_atoms(:,atom2))
          enddo ! i
          atom2p=atom_at_pos(pos_atom2p,structure_sc)
          if(atom2p<=0)call errstop('POINT_SYMM_BRUTE_FORCE','Please check &
            &that your atom coordinates satisfy the rotational symmetries &
            &that you have supplied.  NB, I have assumed that r''=b+Rr, &
            &where R is the structure_sc%rotation_matrices matrix and b is the translation.  &
            &This may be wrong.')

          do i=1,3
            do j=1,3
              if(.NOT.defined(atom1p,i,atom2p,j))then
                fc=0.d0
                do ip=1,3
                  ! LAPACK commented out because it isn't working. 9/1/2017
                  ! fc=fc+ddot(3,structure_sc%rotation_matrices(j,1,n),3, &
                  !   &force_const(atom1,ip,atom2,1),three_noDoFprim_sq) &
                  !   &*structure_sc%rotation_matrices(i,ip,n)
                  fc = fc                                        &
                   & + dot_product(structure_sc%rotation_matrices(j,:,n),              &
                   &               force_const(atom1,ip,atom2,:)) &
                   & * structure_sc%rotation_matrices(i,ip,n)
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
subroutine newtons_law(tol,fc_scale,structure,structure_sc,force_const)
  use structure_module
  implicit none
  
  real(dp), intent(in) :: tol
  real(dp), intent(in) :: fc_scale
  type(StructureData),  intent(in) :: structure
  type(StructureData),  intent(in) :: structure_sc
  
  real(dp), intent(inout) :: force_const(:,:,:,:)
  
  integer, parameter :: min_t = 10     ! min number of impositions.
  integer, parameter :: max_t = 100000 ! max number of impositions.
  
  integer :: atom1,atom2,i,j,t
  real(dp) :: fc,sum1,rescale,max_diff

  do t=1,max_t

    max_diff=0.d0

    ! Impose Newton's third law on matrix of force consts.
    do atom1=1,structure_sc%no_atoms
      do i=1,3
        do j=1,3
          sum1=0.d0
          do atom2=1,atom1-1
            sum1=sum1+force_const(atom1,i,atom2,j)
          enddo
          do atom2=atom1+1,structure_sc%no_atoms
            sum1=sum1+force_const(atom1,i,atom2,j)
          enddo
          rescale=(force_const(atom1,i,atom1,j)+sum1) &
            &/DBLE(structure_sc%no_atoms-1)
          if(ABS(rescale)>max_diff)max_diff=ABS(rescale)
          do atom2=1,atom1-1
            force_const(atom1,i,atom2,j)=force_const(atom1,i,atom2,j) &
              &-rescale
          enddo ! atom2
          do atom2=atom1+1,structure_sc%no_atoms
            force_const(atom1,i,atom2,j)=force_const(atom1,i,atom2,j) &
              &-rescale
          enddo ! atom2
        enddo ! j
      enddo ! i
    enddo ! atom1

    ! Impose symmetry on the matrix of force constants.
    do atom1=1,structure_sc%no_atoms
      do i=1,3
        do atom2=1,structure_sc%no_atoms
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
    if(structure%no_atoms==1)then
      do atom1=1,structure_sc%no_atoms
        do i=1,3
          do atom2=1,structure_sc%no_atoms
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
subroutine mass_reduce(structure_sc,force_const)
  use linear_algebra, only : dscal
  use structure_module
  implicit none
  
  type(StructureData), intent(in) :: structure_sc
  
  real(dp), intent(inout) :: force_const(:,:,:,:)
  
  integer :: atom1,atom2,j
  real(dp) :: rec_root_mass(structure_sc%no_atoms),rec_root_m1m2
  
  do atom1=1,structure_sc%no_atoms
    rec_root_mass(atom1)=1.d0/sqrt(structure_sc%mass(atom1))
  enddo ! atom1
  do atom2=1,structure_sc%no_atoms
    do atom1=1,structure_sc%no_atoms
      rec_root_m1m2=rec_root_mass(atom1)*rec_root_mass(atom2)
      do j=1,3
        ! LAPACK commented out because it isn't working. 9/1/2017
        ! call dscal(3,rec_root_m1m2,force_const(atom1,1,atom2,j),structure_sc%no_atoms)
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
subroutine find_prim_cell(delta,structure_sc,structure,   &
   & no_prim_cells, &
   & atom,atom_in_prim,prim_cell_for_atom,no_equiv_ims,delta_prim)
  use structure_module
  implicit none
  
  real(dp), intent(in) :: delta
  type(StructureData), intent(in) :: structure_sc
  type(StructureData), intent(in) :: structure
  integer,  intent(in) :: no_prim_cells
  
  integer,  intent(inout) :: atom(:,:)
  
  integer,  allocatable, intent(out) :: atom_in_prim(:)
  integer,  allocatable, intent(out) :: prim_cell_for_atom(:)
  integer,  allocatable, intent(out) :: no_equiv_ims(:,:,:)
  real(dp), allocatable, intent(out) :: delta_prim(:,:,:,:,:)
  
  integer :: n,n1,n1p,n1pp,n2,n2p,n2pp,n3,n3p,n3pp,m,label,im, &
    &prim_n(3,structure_sc%no_atoms),p,atom1,ialloc
  real(dp) :: delta_vect(3,3),r_temp(3),delta_r_ims(3,maxim),delta_r_corr(3)

  allocate(atom_in_prim(structure_sc%no_atoms),prim_cell_for_atom(structure_sc%no_atoms), &
    &stat=ialloc)
  if(ialloc/=0)call errstop('FIND_PRIM_CELL','Allocation error: &
    &atom_in_prim, etc.')

  ! Calculate the primitive-lattice point corresponding to the 
  ! primitive cell in which each atom lies.
  ! We try to make sure that there is no uncertainty in the primitive
  ! cell due to atoms sitting on the boundary between two primitive cells.
  delta_vect(1:3,1:3)=delta*structure%lattice(1:3,1:3)
  do n=1,structure_sc%no_atoms
    r_temp(1:3)=structure_sc%cart_atoms(1:3,n)+2.d0*delta_vect(1:3,1) &
      &+2.d0*delta_vect(1:3,2)+2.d0*delta_vect(1:3,3)
    n1=FLOOR(DOT_PRODUCT(r_temp(1:3),structure%recip_lattice(1:3,1)))
    n1p=FLOOR(DOT_PRODUCT(r_temp(1:3)+delta_vect(1:3,1), &
      &structure%recip_lattice(1:3,1)))
    n1pp=FLOOR(DOT_PRODUCT(r_temp(1:3)-delta_vect(1:3,1), &
      &structure%recip_lattice(1:3,1)))
    if(n1/=n1p.OR.n1/=n1pp)call errstop('FIND_PRIM_CELL','Problem &
      &identifying unit cell in which atom lies [1].  Please try increasing &
      &the "delta" parameter in subroutine FIND_PRIM_CELL.')
    n2=FLOOR(DOT_PRODUCT(r_temp(1:3),structure%recip_lattice(1:3,2)))
    n2p=FLOOR(DOT_PRODUCT(r_temp(1:3)+delta_vect(1:3,2), &
      &structure%recip_lattice(1:3,2)))
    n2pp=FLOOR(DOT_PRODUCT(r_temp(1:3)-delta_vect(1:3,2), &
      &structure%recip_lattice(1:3,2)))
    if(n2/=n2p.OR.n2/=n2pp)call errstop('FIND_PRIM_CELL','Problem &
      &identifying unit cell in which atom lies [2].  Please try increasing &
      &the "delta" parameter in subroutine FIND_PRIM_CELL.')
    n3=FLOOR(DOT_PRODUCT(r_temp(1:3),structure%recip_lattice(1:3,3)))
    n3p=FLOOR(DOT_PRODUCT(r_temp(1:3)+delta_vect(1:3,3), &
      &structure%recip_lattice(1:3,3)))
    n3pp=FLOOR(DOT_PRODUCT(r_temp(1:3)-delta_vect(1:3,3), &
      &structure%recip_lattice(1:3,3)))
    if(n3/=n3p.OR.n3/=n3pp)call errstop('FIND_PRIM_CELL','Problem &
      &identifying unit cell in which atom lies [3].  Please try increasing &
      &the "delta" parameter in subroutine FIND_PRIM_CELL.')
    prim_n(1,n)=n1 ; prim_n(2,n)=n2 ; prim_n(3,n)=n3
  enddo ! n

  ! Establish a label for each different atom in the primitive cell,
  ! and evaluate this label for each atom.
  atom_in_prim(1:structure_sc%no_atoms)=-1
  label=0
  do n=1,structure_sc%no_atoms
    if(atom_in_prim(n)==-1)then
      label=label+1
      atom_in_prim(n)=label
      do m=n+1,structure_sc%no_atoms
        ! Is difference of atom positions an integer multiple of the
        ! primitive reciprocal lattice vectors?  If so, same label.
        if(is_lat_point(structure_sc%cart_atoms(1:3,m)-structure_sc%cart_atoms(1:3,n),structure%recip_lattice)) &
          &atom_in_prim(m)=label
      enddo ! m
    endif ! Atom not yet labelled by number within prim cell.
  enddo ! n
  if(label/=structure%no_atoms)call errstop('FIND_PRIM_CELL','Problem &
    &labelling the atoms in the primitive cell.')

  ! Establish a label for each different primitive cell, and evaluate
  ! this label for each atom.
  prim_cell_for_atom(1:structure_sc%no_atoms)=-1
  label=0
  do n=1,structure_sc%no_atoms
    if(prim_cell_for_atom(n)==-1)then
      label=label+1
      prim_cell_for_atom(n)=label
      do m=n+1,structure_sc%no_atoms
        r_temp=DBLE(prim_n(1,m)-prim_n(1,n))*structure%lattice(1:3,1) &
          &+DBLE(prim_n(2,m)-prim_n(2,n))*structure%lattice(1:3,2) &
          &+DBLE(prim_n(3,m)-prim_n(3,n))*structure%lattice(1:3,3)
        if(is_lat_point(r_temp,structure_sc%recip_lattice))prim_cell_for_atom(m)=label
      enddo ! m
    endif ! Atom not yet labelled by number of primitive cell.
  enddo ! n
  if(label/=no_prim_cells)call errstop('FIND_PRIM_CELL','Problem labelling &
    &the primitive cells.')

  ! Construct array holding atom number for a given primitive cell and
  ! atom within the primitive cell.
  atom(:,:)=-1
  do n=1,structure_sc%no_atoms
    atom(prim_cell_for_atom(n),atom_in_prim(n))=n
  enddo ! n
  if(ANY(atom(:,:)==-1))call errstop('FIND_PRIM_CELL','Problem defining atom &
    &labels.')

  ! Work out number of equivalent images and Delta Prim. Lattice Vectors
  ! for pairs of atoms (used in evaluation of force-constant matrix).
  allocate(no_equiv_ims(no_prim_cells,structure%no_atoms,structure%no_atoms), &
    &delta_prim(3,maxim,no_prim_cells,structure%no_atoms, &
    &structure%no_atoms),stat=ialloc)
  if(ialloc/=0)call errstop('FIND_PRIM_CELL','Allocation error: &
    &no_equiv_ims, etc.')
  delta_prim=0.d0
  do n=1,structure%no_atoms
    atom1=atom(1,n)
    do m=1,structure%no_atoms
      delta_r_corr=structure_sc%cart_atoms(1:3,atom1)-structure_sc%cart_atoms(1:3,atom(1,m))
      do p=1,no_prim_cells
        ! Work out min. image distance(s) between atoms (1,n) and (p,m).
        call min_images_brute_force(structure_sc%cart_atoms(1:3,atom(p,m)) &
          &-structure_sc%cart_atoms(1:3,atom1),structure_sc%lattice,delta_r_ims, &
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
subroutine construct_dyn_matrix(kvec,delta_prim,structure,no_prim_cells, &
   & atom,no_equiv_ims,force_const,dyn_mat)
  use structure_module
  implicit none
  
  real(dp), intent(in) :: kvec(3)
  real(dp), intent(in) :: delta_prim(:,:,:,:,:)
  type(StructureData), intent(in) :: structure
  integer,  intent(in) :: no_prim_cells
  integer,  intent(in) :: atom(:,:)
  integer,  intent(in) :: no_equiv_ims(:,:,:)
  real(dp), intent(in) :: force_const(:,:,:,:)
  
  complex(dp), intent(inout) :: dyn_mat(structure%no_modes,structure%no_modes)
  
  integer :: p,n,m,i,j,index1,index2,atom1,im
  complex(dp) :: dm,tempc,expikdotD(no_prim_cells,structure%no_atoms,&
    &structure%no_atoms)
  real(dp) :: k_dot_D

  ! Precompute exp(-ik.(R-R')) to go in the dynamical matrix.
  do n=1,structure%no_atoms
    atom1=atom(1,n)
    do m=1,structure%no_atoms
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
  do n=1,structure%no_atoms
    atom1=atom(1,n)
    do i=1,3
      index1=index1+1
      index2=0
      do m=1,structure%no_atoms
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
  do index1=1,structure%no_modes
    dyn_mat(index1,index1)=CMPLX(real(dyn_mat(index1,index1),dp),0.d0,dp)
    do index2=index1+1,structure%no_modes
      dm=0.5d0*(dyn_mat(index1,index2)+CONJG(dyn_mat(index2,index1)))
      dyn_mat(index1,index2)=dm
      dyn_mat(index2,index1)=CONJG(dm)
    enddo ! index2
  enddo ! index 1
end subroutine

! ----------------------------------------------------------------------
! For a given k, construct and diagonalise the dynamical matrix.
! ----------------------------------------------------------------------
subroutine calculate_eigenfreqs(kvec,delta_prim,structure, &
   & no_prim_cells,atom,no_equiv_ims,force_const,omega)
  use linear_algebra, only : zheev
  use structure_module
  implicit none
  
  real(dp), intent(in) :: kvec(3)
  real(dp), intent(in) :: delta_prim(:,:,:,:,:)
  type(StructureData), intent(in) :: structure
  integer,  intent(in) :: no_prim_cells
  integer,  intent(in) :: atom(:,:)
  integer,  intent(in) :: no_equiv_ims(:,:,:)
  real(dp), intent(in) :: force_const(:,:,:,:)
  
  real(dp),intent(out) :: omega(structure%no_modes)
  
  complex(dp) :: dyn_mat(structure%no_modes,structure%no_modes),work(2*structure%no_modes-1)
  real(dp) :: rwork(3*structure%no_modes-2),minusomegasq(structure%no_modes)
  integer :: n,m,info
  
  ! Construct dynamical matrix.
  call construct_dyn_matrix(kvec,delta_prim,structure,no_prim_cells, &
     & atom,no_equiv_ims,force_const,dyn_mat)

  call zheev('N','U',structure%no_modes,dyn_mat(1,1),structure%no_modes,minusomegasq(1), &
    &work(1),2*structure%no_modes-1,rwork(1),info)
  if(info/=0)call errstop('CALCULATE_EIGENFREQS','ZHEEV failed (1).  Error &
    &code: '//TRIM(i2s(info))//'.')

  ! Eigenvalues of dynamical matrix are minus the frequencies squared.
  ! The eigenvalues are in ascending order, so the +ve frequencies
  ! will be in ascending order.
  m=structure%no_modes
  do n=1,structure%no_modes
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
subroutine calculate_eigenfreqs_and_vecs(kvec,structure,delta_prim, &
   & no_prim_cells,atom,no_equiv_ims,force_const,    &
   & omega,pol_vec)
  use linear_algebra,only : zscal,dznrm2,zheev,zcopy
  use structure_module
  implicit none
  
  real(dp), intent(in) :: kvec(3)
  type(StructureData), intent(in) :: structure
  real(dp), intent(in) :: delta_prim(:,:,:,:,:)
  integer,  intent(in) :: no_prim_cells
  integer,  intent(in) :: atom(:,:)
  integer,  intent(in) :: no_equiv_ims(:,:,:)
  real(dp), intent(in) :: force_const(:,:,:,:)
  
  real(dp),    intent(out) :: omega(structure%no_modes)
  complex(dp), intent(out) :: pol_vec(structure%no_modes,structure%no_modes)
  
  real(dp) :: rwork(3*structure%no_modes-2), &
    &minusomegasq(structure%no_modes)
  complex(dp) :: dyn_mat(structure%no_modes,structure%no_modes),work(2*structure%no_modes-1)
  integer :: n,m,info

  ! Construct dynamical matrix.
  call construct_dyn_matrix(kvec,delta_prim,structure,no_prim_cells, &
     & atom,no_equiv_ims,force_const,dyn_mat)

  call zheev('V','U',structure%no_modes,dyn_mat(1,1),structure%no_modes,minusomegasq(1), &
    &work(1),2*structure%no_modes-1,rwork(1),info)
  if(info/=0)call errstop('CALCULATE_EIGENFREQS_AND_VECS', &
    &'ZHEEV failed (1).  Error code: '//TRIM(i2s(info))//'.')

  m=structure%no_modes
  do n=1,structure%no_modes 
! Modified by B. Monserrat to output the correct 'omega' and 'pol_vec'
    if(minusomegasq(m)>=0.d0)then
      omega(n)=-SQRT(minusomegasq(m)) ! Unstable mode.
    else
      omega(n)=SQRT(-minusomegasq(m)) ! Stable mode.
    endif
    call zcopy(structure%no_modes,dyn_mat(1,m),1,pol_vec(1,n),1)
    m=m-1
  enddo ! n

    ! Attempt to make the eigenvectors real by multiplying
    ! each eigenvector by the conjugate of the largest element.
    ! Then normalise the eigenvectors to the number of atoms in the
    ! primitive cell.  So the polarisation vector for a monatomic
    ! lattice should be normalised to unity.
   !do n=1,structure%no_modes
   ! largest_mag2=0.d0
   ! largest_k=-1
   ! do k=1,structure%no_modes
   !   mag2_element=DBLE(pol_vec(k,n))**2+AIMAG(pol_vec(k,n))**2
   !   if(mag2_element>largest_mag2)then
   !     largest_mag2=mag2_element
   !     largest_k=k
   !   endif
   ! enddo ! k
   ! if(largest_k==-1)call errstop('CALCULATE_EIGENFREQS_AND_VECS', &
   !   &'Eigenvector appears to be zero!')
   ! scalefactor=(DBLE(structure%no_atoms)/(ABS(pol_vec(largest_k,n)) &
   !   &*dznrm2(structure%no_modes,pol_vec(1,n),1)))*CONJG(pol_vec(largest_k,n))
   ! call zscal(structure%no_modes,scalefactor,pol_vec(1,n),1)
   !enddo ! n
end subroutine

! ----------------------------------------------------------------------
! Calculate the frequency density-of-states by Monte Carlo sampling of
! the Brillouin zone.
! ----------------------------------------------------------------------
subroutine calculate_freq_dos(tol,structure,        &
   & delta_prim,no_prim_cells,atom,no_equiv_ims,force_const,&
   & freq_dos_filename,bin_width,freq_dos)
  use constants,      only : max_bin, no_samples, no_fdos_sets, pi
  use linear_algebra, only : dscal
  use rand_no_gen,    only : ranx
  use string_module
  use structure_module
  implicit none
  
  real(dp), intent(in) :: tol
  type(StructureData), intent(in) :: structure
  real(dp), intent(in) :: delta_prim(:,:,:,:,:)
  integer,  intent(in) :: no_prim_cells
  integer,  intent(in) :: atom(:,:)
  integer,  intent(in) :: no_equiv_ims(:,:,:)
  real(dp), intent(in) :: force_const(:,:,:,:)
  
  ! Filename
  type(String), intent(in) :: freq_dos_filename
  
  ! File unit
  integer :: freq_dos_file
  
  real(dp),              intent(out) :: bin_width
  real(dp), allocatable, intent(out) :: freq_dos(:,:)
  
  ! Number of preliminary samples of Brillouin zone to make in order to
  ! establish maximum and minimum frequencies.
  integer,parameter :: no_samples_trial=10000
  ! Our preliminary sampling of the Brillouin zone is imperfect.
  ! Multiply the highest frequency found by this factor when choosing the
  ! highest frequency bin.
  real(dp),parameter :: safety_factor=1.15d0
  
  real(dp) :: omega(structure%no_modes),kvec(3),rec_bin_width,max_freq,min_freq, &
    &rec_no_fdos_sets
  integer :: j,i,n,bin,ialloc
  logical :: soft_modes,soft_modes_prelim

  ! Establish (approximate) maximum and minimum frequencies and hence
  ! choose the bin width.
  max_freq=-1.d0
  min_freq=HUGE(1.d0)
  do i=1,no_samples_trial
    kvec(1:3)=2*pi*(ranx()*structure%recip_lattice(1:3,1) &
      &+ranx()*structure%recip_lattice(1:3,2)+ranx()*structure%recip_lattice(1:3,3))
    call calculate_eigenfreqs(kvec,delta_prim,structure, &
      & no_prim_cells,atom,no_equiv_ims,force_const,omega)
    if(omega(1)<min_freq)min_freq=omega(1)
    if(omega(structure%no_modes)>max_freq)max_freq=omega(structure%no_modes)
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
      kvec(1:3)=2*pi*(ranx()*structure%recip_lattice(1:3,1) &
        &+ranx()*structure%recip_lattice(1:3,2)+ranx()*structure%recip_lattice(1:3,3))
      call calculate_eigenfreqs(kvec,delta_prim,structure, &
        & no_prim_cells,atom,no_equiv_ims,force_const,omega)
      if(omega(1)<-tol)soft_modes=.TRUE.
      do n=1,structure%no_modes
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
  ! number of frequencies sampled is 3*no_samples*structure%no_atoms.
  ! (Imaginary frequencies are ignored, however.)
  call dscal((max_bin+1)*no_fdos_sets,1.d0/(DBLE(no_samples)*bin_width), &
    &freq_dos(0,1),1)

  ! Write out the frequency DoS.
  rec_no_fdos_sets=1.d0/DBLE(no_fdos_sets)
  freq_dos_file = open_write_file(freq_dos_filename)
  do bin=0,max_bin
    write(freq_dos_file,*) (bin+0.5_dp)*bin_width, &
                         & sum(freq_dos(bin,1:no_fdos_sets))*rec_no_fdos_sets
  enddo
  close(freq_dos_file)
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
  use string_module
  implicit none
  
  real(dp), intent(in) :: bin_width
  real(dp), intent(in) :: temperature
  real(dp), intent(in) :: freq_dos(:,:)
  
  ! File name
  type(String), intent(in) :: tdependence1_filename
  
  ! File unit
  integer :: tdependence1_file
  
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
   
  tdependence1_file = open_write_file(tdependence1_filename)
  write(tdependence1_file,*) lte*27.211396132d0
  close(tdependence1_file)
end subroutine

! ----------------------------------------------------------------------
! This function returns the mean thermal energy of an isolated harmonic
! oscillator of frequency omega (in a.u.).  T is the temperature in Kelvin.
! ----------------------------------------------------------------------
real(dp) function harmonic_energy(T,omega)
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
  use string_module
  implicit none
  
  real(dp), intent(in) :: bin_width
  real(dp), intent(in) :: temperature
  real(dp), intent(in) :: freq_dos(:,:)
  
  ! File name
  type(String), intent(in) :: tdependence2_filename
  
  ! File unit
  integer :: tdependence2_file
  
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
  
  tdependence2_file = open_write_file(tdependence2_filename)
  write(tdependence2_file,*) ltfe*27.211396132d0
  close(tdependence2_file)
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
subroutine generate_disp_curve(structure,delta_prim,no_kspace_lines, &
   & disp_kpoints,no_prim_cells,atom,no_equiv_ims,    &
   & force_const,                                                      &
   & dispersion_curve_filename)
  use string_module
  use structure_module
  implicit none
  
  type(StructureData), intent(in) :: structure
  real(dp), intent(in) :: delta_prim(:,:,:,:,:)
  integer,  intent(in) :: no_kspace_lines
  real(dp), intent(in) :: disp_kpoints(:,:)
  integer,  intent(in) :: no_prim_cells
  integer,  intent(in) :: atom(:,:)
  integer,  intent(in) :: no_equiv_ims(:,:,:)
  real(dp), intent(in) :: force_const(:,:,:,:)
  
  ! File name
  type(String), intent(in) :: dispersion_curve_filename
  
  ! File unit
  integer :: dispersion_curve_file
  
  real(dp) :: k_dist,kvec(3),delta_k(3),k_step,omega(structure%no_modes)
  integer :: i,j,k,total_no_kpoints,ialloc
  real(dp),allocatable :: disp_k_dist(:),branch(:,:)
  integer,parameter :: no_kpoints_per_line=1000

  write(*,*)
  write(*,*)'Number of k points per line in k space : ' &
    &//TRIM(i2s(no_kpoints_per_line))
  write(*,*)

  ! Total number of k points at which the dispersion curve is to be calc'd.
  total_no_kpoints=no_kspace_lines*no_kpoints_per_line
  allocate(disp_k_dist(total_no_kpoints), &
    &branch(structure%no_modes,total_no_kpoints),stat=ialloc)
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
      call calculate_eigenfreqs(kvec,delta_prim,structure, &
        & no_prim_cells,atom,no_equiv_ims,force_const,omega)
      branch(1:structure%no_modes,k)=omega(1:structure%no_modes)
      disp_k_dist(k)=k_dist
      kvec(1:3)=kvec(1:3)+delta_k(1:3)
      k_dist=k_dist+k_step
    enddo ! j
  enddo ! i
  k_dist=k_dist-k_step
  write(*,*)'Final line ends at k-space distance   : ',k_dist
  
  dispersion_curve_file = open_write_file(dispersion_curve_filename)
  do j=1,structure%no_modes
    do k=1,total_no_kpoints
      write(dispersion_curve_file,*) disp_k_dist(k), branch(j,k)
    enddo
    write(dispersion_curve_file,*) '&'
  enddo
  close(dispersion_curve_file)

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
subroutine calculate_speed_sound(no_prim_cells, &
   & structure,structure_sc,delta_prim,atom,no_equiv_ims,force_const)
  use constants,   only : pi
  use rand_no_gen, only : ranx
  use structure_module
  implicit none
  
  integer,  intent(in) :: no_prim_cells
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: structure_sc
  real(dp), intent(in) :: delta_prim(:,:,:,:,:)
  integer,  intent(in) :: atom(:,:)
  integer,  intent(in) :: no_equiv_ims(:,:,:)
  real(dp), intent(in) :: force_const(:,:,:,:)
  
  real(dp) :: kvec(3),omega(3),cos_theta,sin_theta,phi,c_tr_tot,c_tr, &
    &c2_tr_tot,c2_tr,c_ln_tot,c_ln,c2_ln_tot,c2_ln,err_tr,err_ln,c(3), &
    &kunit(3),pol_vec_real(3,3),dot_prod(3),temp,c_tr_old,c_ln_old
  complex(dp) :: pol_vec(3,3)
  integer :: i,no_samples,k,k2
  real(dp),parameter :: err_tol=1.d-3
  integer,parameter :: max_samples=1000000
  logical,parameter :: verbose=.FALSE.
  
  real(dp) :: small_k_scale
  real(dp) :: kmag

  if(structure%no_atoms/=1)call errstop('CALCULATE_SPEED_SOUND', &
    &'At the moment this program can only work out the speed of sound in &
    &materials with a single atom per primitive cell.  Sorry about that.')
  
  ! "Small" distance in k space, determined by size of supercell.
  ! Factor of 2pi included.
  small_k_scale=2*pi*structure_sc%volume**(-1/3.0_dp)

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
      cos_theta=1.d0-2.d0*ranx()
      phi=ranx()*2*pi
      sin_theta=SQRT(1.d0-cos_theta**2)
      kunit(1:3)=(/sin_theta*COS(phi),sin_theta*SIN(phi),cos_theta/)
      kvec(1:3)=kmag*kunit

      ! Calculate corresponding eigenfrequencies.
      call calculate_eigenfreqs_and_vecs(kvec,structure,delta_prim,   &
        & no_prim_cells,atom,no_equiv_ims,force_const, &
        & omega,pol_vec)
      if(ANY(omega<0.d0))then
        write(*,*)'Imaginary frequencies found.'
        write(*,*)'In terms of the primitive reciprocal lattice vectors, &
          &the k-point is:'
        write(*,*)DOT_PRODUCT(kvec,structure%lattice(1:3,1)), &
          &DOT_PRODUCT(kvec,structure%lattice(1:3,2)), &
          &DOT_PRODUCT(kvec,structure%lattice(1:3,3))
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
subroutine evaluate_freqs_on_grid(no_prim_cells,    &
   & structure_sc,temperature,structure, &
   & delta_prim,atom,no_equiv_ims,force_const,         &
   & prim_cell_for_atom,atom_in_prim,                                &
   & kpairs_filename,freq_grids_filename,disp_patterns_filename,     &
   & kdisp_patterns_filename,pol_vec_filename,gvectors_filename,     &
   & gvectors_frac_filename,error_filename)
  use constants, only : pi
  use utils,     only : reduce_interval
  use string_module
  use structure_module
  implicit none
  
  integer,  intent(in) :: no_prim_cells
  type(StructureData),  intent(in) :: structure_sc
  real(dp), intent(in) :: temperature
  type(StructureData), intent(in) :: structure
  real(dp), intent(in) :: delta_prim(:,:,:,:,:)
  integer,  intent(in) :: atom(:,:)
  integer,  intent(in) :: no_equiv_ims(:,:,:)
  real(dp), intent(in) :: force_const(:,:,:,:)
  integer,  intent(in) :: prim_cell_for_atom(:)
  integer,  intent(in) :: atom_in_prim(:)
  
  ! File names
  type(String), intent(in) :: kpairs_filename
  type(String), intent(in) :: freq_grids_filename
  type(String), intent(in) :: disp_patterns_filename
  type(String), intent(in) :: kdisp_patterns_filename
  type(String), intent(in) :: pol_vec_filename
  type(String), intent(in) :: gvectors_filename
  type(String), intent(in) :: gvectors_frac_filename
  type(String), intent(in) :: error_filename
  
  ! File units
  integer :: kpairs_file
  integer :: freq_grids_file
  integer :: disp_patterns_file
  integer :: kdisp_patterns_file
  integer :: pol_vec_file
  integer :: gvectors_file
  integer :: gvectors_frac_file
  integer :: error_file
  
  integer :: ng,i,j,k,ig,index1,index2,p,n,atom1
  logical :: found,soft_modes
  real(dp) :: gnew(3),gvec(3,no_prim_cells),R0(3), &
    &omega(structure%no_modes),E,F,rec_root_mass(structure%no_atoms),GdotR, &
    &disp_pattern(3),kdisp_pattern(3),tot_disp_patt
  complex(dp) :: pol_vec(structure%no_modes,structure%no_modes), &
    &non_mr_pol_vec(3,structure%no_atoms),expiGdotR(no_prim_cells), &
    &kpol_vec(3,structure%no_atoms)
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
    do j=0,no_prim_cells-1
      do i=0,no_prim_cells-1
        gnew = matmul(structure_sc%recip_lattice,(/i,j,k/))
        found=.TRUE.
        do ig=1,ng
          if(is_lat_point(gnew-gvec(:,ig),structure%lattice))then
            found=.FALSE.
            exit
          endif ! ig
        enddo ! ig
        if(found)then
          ng=ng+1
          if(ng>no_prim_cells)call errstop('EVALUATE_FREQS_ON_GRID', &
            &'Bug: too many G vectors.')
          gvec(:,ng) = gnew
        endif ! found
      enddo ! i
    enddo ! j
  enddo ! k
  if(ng/=no_prim_cells)call errstop('EVALUATE_FREQS_ON_GRID', &
    &'Bug: too few G vectors.')
  

  ! Calculate +/- G-vector pairs
  ! First, write G-vectors as fractions of rec. latt. vecs.
  gfrac = reduce_interval(matmul(transpose(structure%lattice),gvec),tol_g)
  ! Second, pair them up
  reference=0
  do k=1,no_prim_cells
    do j=1,k-1
      if (all(abs(reduce_interval(gfrac(:,j)+gfrac(:,k),tol_g))<tol_g)) then
        reference(k)=j 
        reference(j)=k
        exit
      endif
    enddo
  enddo
  
  ! Scale gvec by 2 pi
  gvec = 2*pi*gvec
  
  kpairs_file = open_write_file(kpairs_filename)
  do i=1,no_prim_cells
    write(kpairs_file,*)i,reference(i)
  enddo
  close(kpairs_file)
  
  freq_grids_file = open_write_file(freq_grids_filename)
  disp_patterns_file = open_write_file(disp_patterns_filename)
  kdisp_patterns_file = open_write_file(kdisp_patterns_filename)
  pol_vec_file = open_write_file(pol_vec_filename)

  do n=1,structure%no_atoms
    rec_root_mass(n)=1.d0/SQRT(structure_sc%mass(atom(1,n))) ! 1/sqrt(m) in prim. cell.
  enddo ! n

  ! Modified by B. Monserrat to output G vectors to file
  gvectors_file = open_write_file(gvectors_filename)
  gvectors_frac_file = open_write_file(gvectors_frac_filename)
  write(gvectors_file,*) no_prim_cells
  write(gvectors_frac_file,*) no_prim_cells
  do ig=1,no_prim_cells
    write(*,'(" G = (",es20.12,",",es20.12,",",es20.12,")")') gvec(:,ig)
    write(*,'(" G = (",es20.12,",",es20.12,",",es20.12,")")') gfrac(:,ig)
    write(gvectors_file,*) gvec(:,ig)
    ! G-vectors as a fraction of the primitive reciprocal lattice vectors
    write(gvectors_frac_file,*) ig, gfrac(:,ig)
  enddo
  close(gvectors_file)
  close(gvectors_frac_file)

  ! Evaluate the frequencies at each supercell G vector.
  E=0.d0  ;  F=0.d0
  soft_modes=.FALSE.
  R0=structure_sc%cart_atoms(:,atom(1,1))
  do ig=1,no_prim_cells
    if(reference(ig)==0)then
      call calculate_eigenfreqs_and_vecs(gvec(:,ig),structure, &
        & delta_prim,no_prim_cells,atom,no_equiv_ims,   &
        & force_const,omega,pol_vec)
    else
      if(reference(ig)>ig)then
        call calculate_eigenfreqs_and_vecs(gvec(:,ig),structure, &
          & delta_prim,no_prim_cells,atom,no_equiv_ims,   &
          & force_const,omega,pol_vec)
      else
        call calculate_eigenfreqs_and_vecs(gvec(:,reference(ig)), &
          & structure,                                                  &
          & delta_prim,no_prim_cells,atom,no_equiv_ims,  &
          & force_const,omega,pol_vec)
      endif
    endif

    ! The negative is used because the matrix of force constants is the transpose of
    ! the usual expression in derivations that lead to a positive exponential
    do p=1,no_prim_cells
      if(reference(ig)==0)then
        GdotR=-DOT_PRODUCT(gvec(:,ig),structure_sc%cart_atoms(1:3,atom(p,1))-R0)
      else
        if(reference(ig)>ig)then
          GdotR=-DOT_PRODUCT(gvec(:,ig),structure_sc%cart_atoms(1:3,atom(p,1))-R0)
        else
          GdotR=-DOT_PRODUCT(gvec(:,reference(ig)),structure_sc%cart_atoms(1:3,atom(p,1))-R0)
        endif
      endif
      expiGdotR(p)=CMPLX(COS(GdotR),SIN(GdotR),dp) ! Store exp(iG.R_p).
    enddo ! p

    do index2=1,structure%no_modes
      write(*,*)'  omega = ',omega(index2)
      if(omega(index2)>tol_omega)then
    ! Ignore contributions from imaginary or zero frequencies.
        E=E+harmonic_energy(temperature,omega(index2))
        F=F+harmonic_free_energy(temperature,omega(index2))
      elseif(omega(index2)<-tol_omega)then
        soft_modes=.TRUE.
      endif ! omega>0
      write(freq_grids_file,*)omega(index2),1.d0
! Compute the non-mass-reduced polarisation vector.
      index1=0
      do n=1,structure%no_atoms
        do i=1,3
          index1=index1+1
          non_mr_pol_vec(i,n)=pol_vec(index1,index2)*rec_root_mass(n)
          kpol_vec(i,n)=pol_vec(index1,index2)
        enddo ! i
      enddo ! n
      if(omega(index2)<-tol_omega)then
        write(disp_patterns_file,*)'Frequency           : ',omega(index2),' (SOFT)'
        write(kdisp_patterns_file,*)'Frequency           : ',omega(index2),' (SOFT)'
        write(pol_vec_file,*)'Mode number     :',index2, '   Frequency           : ',omega(index2),' (SOFT)'
      else
        write(disp_patterns_file,*)'Frequency           : ',omega(index2)
        write(kdisp_patterns_file,*)'Frequency           : ',omega(index2)
        write(pol_vec_file,*)'Mode number     :',index2, '   Frequency           : ',omega(index2)
      endif ! soft freq.
      write(disp_patterns_file,*) gvec(:,ig)
      write(kdisp_patterns_file,*) gvec(:,ig)
      write(pol_vec_file,*) gvec(:,ig)
      write(disp_patterns_file,*)'Displacement pattern for each atom:'
      write(kdisp_patterns_file,*)'Displacement pattern for each atom:'
      write(pol_vec_file,*)'Polarisation vector:'
      disp_pattern=0.d0
      kdisp_pattern=0.d0
      tot_disp_patt=0.d0
      do atom1=1,structure_sc%no_atoms
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
        write(disp_patterns_file,*)disp_pattern,prefactor
        write(kdisp_patterns_file,*)kdisp_pattern,prefactor
        write(pol_vec_file,*)real(non_mr_pol_vec(1:3,atom_in_prim(atom1)))
        write(pol_vec_file,*)AIMAG(non_mr_pol_vec(1:3,atom_in_prim(atom1)))
      enddo ! atom1
      write(disp_patterns_file,*)
      write(kdisp_patterns_file,*)
      write(pol_vec_file,*)
    enddo ! index2
    write(*,*)
    if(tot_disp_patt<1.d-8)then
      error_file = open_write_file(error_filename)
      write(error_file,*)'The total displacement is:',tot_disp_patt
      close(error_file)
    endif
  enddo ! ig
  E=E/DBLE(no_prim_cells)  ;  F=F/DBLE(no_prim_cells)

  close(freq_grids_file)
  close(disp_patterns_file)
  close(kdisp_patterns_file)
  close(pol_vec_file)

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
subroutine write_dynamical_matrix(structure_sc,structure, &
   & no_prim_cells,delta_prim,atom,no_equiv_ims,      &
   & force_const,                                                      &
   & dyn_mat_fileroot)
  use constants, only : pi
  use string_module
  use structure_module
  implicit none
  
  type(StructureData), intent(in) :: structure_sc
  type(StructureData), intent(in) :: structure
  integer,  intent(in) :: no_prim_cells
  real(dp), intent(in) :: delta_prim(:,:,:,:,:)
  integer,  intent(in) :: atom(:,:)
  integer,  intent(in) :: no_equiv_ims(:,:,:)
  real(dp), intent(in) :: force_const(:,:,:,:)
  
  ! File root (file name minus ending)
  type(String), intent(in) :: dyn_mat_fileroot
  
  ! File unit
  integer :: dyn_mat_file
  
  integer :: ng,k,j,i,ig,atom1,cart1,index1,atom2,cart2,index2
  real(dp) :: gnew(3),gvec(3,no_prim_cells)
  complex(dp) :: dyn_mat(structure%no_modes,structure%no_modes)
  logical :: found

  ng=0
  do k=0,no_prim_cells-1
    do j=0,no_prim_cells-1
      do i=0,no_prim_cells-1
        gnew = matmul(structure_sc%recip_lattice,(/i,j,k/))
        found=.TRUE.
        do ig=1,ng
          if(is_lat_point(gnew-gvec(:,ig),structure%lattice))then
            found=.FALSE.
            exit
          endif ! ig
        enddo ! ig
        if(found)then
          ng=ng+1
          if(ng>no_prim_cells)call errstop('WRITE_DYNAMICAL_MATRIX','Too &
           &many G-vectors.')
          gvec(:,ng)=gnew(:)
        endif ! found
      enddo ! i
    enddo ! j
  enddo ! k
  if(ng/=no_prim_cells)call errstop('WRITE_DYNAMICAL_MATRIX','Too few &
   &G-vectors.')

  dyn_mat = cmplx(0.0_dp,0.0_dp,dp)
  gvec = 2*pi*gvec

  do ig=1,ng
    dyn_mat_file = open_write_file(dyn_mat_fileroot//'.'//ig//'.dat')
    call construct_dyn_matrix(gvec(:,ig),delta_prim,structure, &
      & no_prim_cells,atom,no_equiv_ims,force_const,dyn_mat)
    atom1=0
    do index1=1,structure%no_modes
      atom2=0
      if(MOD(index1,3)==1)then
        atom1=atom1+1
        cart1=1
      endif ! MOD(index1,3)==1
      do index2=1,structure%no_modes
        if(MOD(index2,3)==1)then
          atom2=atom2+1
          cart2=1
        endif ! MOD(index2,3)==1
        write(dyn_mat_file,*) atom1,                        &
                            & cart1,                        &
                            & atom2,                        &
                            & cart2,                        &
                            & real(dyn_mat(index1,index2)), &
                            & AIMAG(dyn_mat(index1,index2))
        cart2=cart2+1
      enddo ! index2
      cart1=cart1+1
    enddo ! index1
    close(dyn_mat_file)
  end do ! ig
end subroutine

! ----------------------------------------------------------------------
! Write out atoms in primitive cell in order.
! ----------------------------------------------------------------------
subroutine write_atoms_in_primitive_cell(structure,structure_sc, &
   & atom,atoms_in_primitive_cell_filename)
  use string_module
  use structure_module
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: structure_sc
  integer,  intent(in) :: atom(:,:)
  
  ! File name
  type(String), intent(in) :: atoms_in_primitive_cell_filename
  
  ! File unit
  integer :: atoms_in_primitive_cell_file
  
  real(dp),parameter :: tol=1.d-10
  
  integer :: n,atom1,i
  real(dp) :: pos(3),frac(3)
  
  atoms_in_primitive_cell_file = open_write_file(atoms_in_primitive_cell_filename)

  do n=1,structure%no_atoms
    atom1=atom(1,n)
    pos=structure_sc%cart_atoms(1:3,atom1)
    do i=1,3
      frac(i)=dot_product(pos(1:3),structure%recip_lattice(1:3,i))
    enddo
    frac(1:3)=modulo(frac(1:3)+tol,1.d0)-tol
    write(atoms_in_primitive_cell_file,*) structure_sc%mass(atom1), frac(1:3)
  enddo

  close(atoms_in_primitive_cell_file)
end subroutine

! ----------------------------------------------------------------------
! Set up lte
! ----------------------------------------------------------------------
subroutine initialise(prog_function,tol,tol2,delta,structure,structure_sc, &
   & atoms,displacements,forces,temperature,no_kspace_lines,disp_kpoints, &
   & length_scale,vol_scale,fc_scale, &
   & no_prim_cells,                &
   & force_const,delta_prim,     &
   & atom,no_equiv_ims,atom_in_prim,prim_cell_for_atom,defined)
  use constants, only : dp
  use utils,     only : errstop, wordwrap
  use string_module
  use structure_module
  implicit none
  
  ! ----------------------------------------
  ! Inputs
  ! ----------------------------------------
  integer,  intent(in) :: prog_function
  real(dp), intent(in) :: tol
  real(dp), intent(in) :: tol2
  real(dp), intent(in) :: delta
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: structure_sc
  integer,             intent(in) :: atoms(:)
  integer,             intent(in) :: displacements(:)
  real(dp),            intent(in) :: forces(:,:,:)
  real(dp),            intent(in) :: temperature
  integer,             intent(in) :: no_kspace_lines
  real(dp),            intent(in) :: disp_kpoints(:,:)
  
  ! ----------------------------------------
  ! previously global variables
  ! ----------------------------------------
  real(dp), intent(out) :: length_scale
  real(dp), intent(out) :: vol_scale
  real(dp), intent(out) :: fc_scale
  integer,  intent(out)  :: no_prim_cells
  
  real(dp),     allocatable, intent(out) :: force_const(:,:,:,:)
  real(dp),     allocatable, intent(out) :: delta_prim(:,:,:,:,:)
  integer,      allocatable, intent(out) :: atom(:,:)
  integer,      allocatable, intent(out) :: no_equiv_ims(:,:,:)
  integer,      allocatable, intent(out) :: atom_in_prim(:)
  integer,      allocatable, intent(out) :: prim_cell_for_atom(:)
  logical,      allocatable, intent(out) :: defined(:,:,:,:)
  
  write(*,*)
  write(*,*)'LATTICE THERMAL ENERGY'
  write(*,*)'======================'
  write(*,*)

  write(*,*)'Reading data from lte.dat...'
  write(*,*)
  call read_lte(prog_function,tol,structure,structure_sc,atoms,displacements, &
     & forces,temperature,no_kspace_lines,disp_kpoints, &
     & fc_scale,           &
     & no_prim_cells,length_scale,&
     & vol_scale,                                                            &
     & force_const,defined,atom)
  write(*,*)'Finished reading input data.'
  write(*,*)

  write(*,*)'Applying point symmetries to the matrix of force constants...'
  call point_symm(tol,structure_sc,defined,force_const)
  write(*,*)'Done.'
  write(*,*)

  if(ANY(.NOT.defined))then
    call wordwrap('WARNING: will impose symmetries on the matrix of force &
      &constants iteratively...')
    call point_symm_brute_force(tol,fc_scale,structure_sc,     &
       & structure,force_const,defined)
    write(*,*)'Done.'
    write(*,*)
  endif
  if(ANY(.NOT.defined))call errstop('LTE','Some elements of the matrix of &
    &force constants are still undefined.')

  call wordwrap('Imposing Newton''s third law and symmetry on the matrix of &
    &force constants...')
  call newtons_law(tol2,fc_scale,structure,structure_sc,force_const)
  write(*,*)'Done.'
  write(*,*)

  write(*,*)'Performing mass reduction on the matrix of force constants...'
  call mass_reduce(structure_sc,force_const)
  write(*,*)'Done.'
  write(*,*)

  write(*,*)'Establishing the primitive lattice vector associated with &
    &each atom...'
  call find_prim_cell(delta,structure_sc,structure,         &
     & no_prim_cells, &
     & atom,atom_in_prim,prim_cell_for_atom,no_equiv_ims,                 &
     & delta_prim)
  write(*,*)'Done.'
  write(*,*)
end subroutine

! ----------------------------------------------------------------------
! Main program. (Split into four modes)
! ----------------------------------------------------------------------
subroutine lte_1(tol,tol2,delta,structure,structure_sc,atoms, &
   & displacements,forces,temperature,no_kspace_lines,disp_kpoints,freq_dos_filename, &
   & tdependence1_filename,tdependence2_filename)
  use constants, only : dp
  use utils,     only : errstop, wordwrap
  use string_module
  use structure_module
  implicit none
  
  ! ----------------------------------------
  ! Inputs
  ! ----------------------------------------
  real(dp),            intent(in) :: tol
  real(dp),            intent(in) :: tol2
  real(dp),            intent(in) :: delta
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: structure_sc
  integer,             intent(in) :: atoms(:)
  integer,             intent(in) :: displacements(:)
  real(dp),            intent(in) :: forces(:,:,:)
  real(dp),            intent(in) :: temperature
  integer,             intent(in) :: no_kspace_lines
  real(dp),            intent(in) :: disp_kpoints(:,:)
  
  ! ----------------------------------------
  ! filenames
  ! ----------------------------------------
  type(String), intent(in) :: freq_dos_filename
  type(String), intent(in) :: tdependence1_filename
  type(String), intent(in) :: tdependence2_filename
  
  ! ----------------------------------------
  ! Program mode
  ! ----------------------------------------
  integer :: prog_function = 1
  
  ! ----------------------------------------
  ! times
  ! ----------------------------------------
  real :: t1,t2
  
  ! ----------------------------------------
  ! previously global variables
  ! ----------------------------------------
  real(dp) :: length_scale
  real(dp) :: vol_scale
  real(dp) :: fc_scale
  real(dp) :: bin_width
  integer  :: no_prim_cells
  
  real(dp),     allocatable :: force_const(:,:,:,:)
  real(dp),     allocatable :: freq_dos(:,:)
  real(dp),     allocatable :: delta_prim(:,:,:,:,:)
  integer,      allocatable :: atom(:,:)
  integer,      allocatable :: no_equiv_ims(:,:,:)
  integer,      allocatable :: atom_in_prim(:)
  integer,      allocatable :: prim_cell_for_atom(:)
  logical,      allocatable :: defined(:,:,:,:)
  
  call cpu_time(t1)
  
  call initialise(prog_function,tol,tol2,delta,structure,structure_sc,atoms, &
     & displacements,forces,temperature,no_kspace_lines,disp_kpoints, &
     & length_scale,vol_scale,fc_scale,&
     & no_prim_cells,               &
     & force_const,delta_prim,    &
     & atom,no_equiv_ims,atom_in_prim,prim_cell_for_atom,defined)

  write(*,*)'The mean thermal energy and the free energy will &
    &be calculated.'
  write(*,*)'Calculating the frequency density-of-states function...'
  call calculate_freq_dos(tol,structure,   &
     & delta_prim,no_prim_cells,atom,no_equiv_ims, &
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
  
  call cpu_time(t2)
  
  write(*,*)'Program finished.  Time taken: ',t2-t1
  write(*,*)
end subroutine

subroutine lte_2(tol,tol2,delta,structure,structure_sc,atoms, &
   & displacements,forces,temperature,no_kspace_lines,disp_kpoints, &
   & dispersion_curve_filename)
  use constants, only : dp
  use utils,     only : errstop, wordwrap
  use string_module
  use structure_module
  implicit none
  
  ! ----------------------------------------
  ! Inputs
  ! ----------------------------------------
  real(dp),            intent(in) :: tol
  real(dp),            intent(in) :: tol2
  real(dp),            intent(in) :: delta
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: structure_sc
  integer,             intent(in) :: atoms(:)
  integer,             intent(in) :: displacements(:)
  real(dp),            intent(in) :: forces(:,:,:)
  real(dp),            intent(in) :: temperature
  integer,             intent(in) :: no_kspace_lines
  real(dp),            intent(in) :: disp_kpoints(:,:)
  
  ! ----------------------------------------
  ! filenames
  ! ----------------------------------------
  type(String), intent(in) :: dispersion_curve_filename
  
  ! ----------------------------------------
  ! Program mode
  ! ----------------------------------------
  integer :: prog_function = 2
  
  ! ----------------------------------------
  ! times
  ! ----------------------------------------
  real :: t1,t2
  
  ! ----------------------------------------
  ! previously global variables
  ! ----------------------------------------
  real(dp) :: length_scale
  real(dp) :: vol_scale
  real(dp) :: fc_scale
  integer  :: no_prim_cells
  
  real(dp),     allocatable :: force_const(:,:,:,:)
  real(dp),     allocatable :: delta_prim(:,:,:,:,:)
  integer,      allocatable :: atom(:,:)
  integer,      allocatable :: no_equiv_ims(:,:,:)
  integer,      allocatable :: atom_in_prim(:)
  integer,      allocatable :: prim_cell_for_atom(:)
  logical,      allocatable :: defined(:,:,:,:)
  
  call cpu_time(t1)
  
  call initialise(prog_function,tol,tol2,delta,structure,structure_sc,atoms, &
     & displacements,forces,temperature,no_kspace_lines,disp_kpoints, &
     & length_scale,vol_scale,fc_scale,&
     & no_prim_cells,               &
     & force_const,delta_prim,    &
     & atom,no_equiv_ims,atom_in_prim,prim_cell_for_atom,defined)

  write(*,*)'A dispersion curve will be calculated.'
  write(*,*)'Calculating the requested dispersion curve.'
  call generate_disp_curve(structure,delta_prim,no_kspace_lines,    &
     & disp_kpoints,no_prim_cells,atom,no_equiv_ims, &
     & force_const,                                                   &
     & dispersion_curve_filename)
  call wordwrap('Done.  dispersion_curve.dat has been generated.  (Please &
    &view this file using XMGrace.)')
  write(*,*)
  
  call cpu_time(t2)
  
  write(*,*)'Program finished.  Time taken: ',t2-t1
  write(*,*)
end subroutine

subroutine lte_3(tol,tol2,delta,structure,structure_sc,atoms, &
   & displacements,forces,temperature,no_kspace_lines,disp_kpoints)
  use constants, only : dp
  use utils,     only : errstop, wordwrap
  use string_module
  use structure_module
  implicit none
  
  ! ----------------------------------------
  ! Inputs
  ! ----------------------------------------
  real(dp),            intent(in) :: tol
  real(dp),            intent(in) :: tol2
  real(dp),            intent(in) :: delta
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: structure_sc
  integer,             intent(in) :: atoms(:)
  integer,             intent(in) :: displacements(:)
  real(dp),            intent(in) :: forces(:,:,:)
  real(dp),            intent(in) :: temperature
  integer,             intent(in) :: no_kspace_lines
  real(dp),            intent(in) :: disp_kpoints(:,:)
  
  ! ----------------------------------------
  ! Program mode
  ! ----------------------------------------
  integer :: prog_function = 3
  
  ! ----------------------------------------
  ! times
  ! ----------------------------------------
  real :: t1,t2
  
  ! ----------------------------------------
  ! previously global variables
  ! ----------------------------------------
  real(dp) :: length_scale
  real(dp) :: vol_scale
  real(dp) :: fc_scale
  integer  :: no_prim_cells
  
  real(dp),     allocatable :: force_const(:,:,:,:)
  real(dp),     allocatable :: delta_prim(:,:,:,:,:)
  integer,      allocatable :: atom(:,:)
  integer,      allocatable :: no_equiv_ims(:,:,:)
  integer,      allocatable :: atom_in_prim(:)
  integer,      allocatable :: prim_cell_for_atom(:)
  logical,      allocatable :: defined(:,:,:,:)
  
  call cpu_time(t1)
  
  call initialise(prog_function,tol,tol2,delta,structure,structure_sc,atoms, &
     & displacements,forces,temperature,no_kspace_lines,disp_kpoints, &
     & length_scale,vol_scale,fc_scale,&
     & no_prim_cells,               &
     & force_const,delta_prim,    &
     & atom,no_equiv_ims,atom_in_prim,prim_cell_for_atom,defined)

  write(*,*)'The speed of sound will be calculated.'
  write(*,*)'Calculating the speed of sound.'
  call calculate_speed_sound(no_prim_cells, &
     & structure,structure_sc,delta_prim,atom,no_equiv_ims,force_const)
  write(*,*)'Done.  Speed of sound calculated.'
  write(*,*)

  
  call cpu_time(t2)
  
  write(*,*)'Program finished.  Time taken: ',t2-t1
  write(*,*)
end subroutine

subroutine lte_4(tol,tol2,delta,structure,structure_sc,atoms, &
   & displacements,forces,temperature,no_kspace_lines,disp_kpoints, &
   & kpairs_filename,freq_grids_filename,disp_patterns_filename,            &
   & kdisp_patterns_filename,pol_vec_filename,gvectors_filename,            &
   & gvectors_frac_filename,error_filename,dyn_mat_fileroot,                &
   & atoms_in_primitive_cell_filename)
  use constants, only : dp
  use utils,     only : errstop, wordwrap
  use string_module
  use structure_module
  implicit none
  
  ! ----------------------------------------
  ! Inputs
  ! ----------------------------------------
  real(dp),            intent(in) :: tol
  real(dp),            intent(in) :: tol2
  real(dp),            intent(in) :: delta
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: structure_sc
  integer,             intent(in) :: atoms(:)
  integer,             intent(in) :: displacements(:)
  real(dp),            intent(in) :: forces(:,:,:)
  real(dp),            intent(in) :: temperature
  integer,             intent(in) :: no_kspace_lines
  real(dp),            intent(in) :: disp_kpoints(:,:)
  
  ! ----------------------------------------
  ! filenames
  ! ----------------------------------------
  type(String), intent(in) :: kpairs_filename
  type(String), intent(in) :: freq_grids_filename
  type(String), intent(in) :: disp_patterns_filename
  type(String), intent(in) :: kdisp_patterns_filename
  type(String), intent(in) :: pol_vec_filename
  type(String), intent(in) :: gvectors_filename
  type(String), intent(in) :: gvectors_frac_filename
  type(String), intent(in) :: error_filename
  type(String), intent(in) :: dyn_mat_fileroot ! will have *.dat appended
  type(String), intent(in) :: atoms_in_primitive_cell_filename
  
  ! ----------------------------------------
  ! Program mode
  ! ----------------------------------------
  integer :: prog_function = 4
  
  ! ----------------------------------------
  ! times
  ! ----------------------------------------
  real :: t1,t2
  
  ! ----------------------------------------
  ! previously global variables
  ! ----------------------------------------
  real(dp) :: length_scale
  real(dp) :: vol_scale
  real(dp) :: fc_scale
  integer  :: no_prim_cells
  
  real(dp),     allocatable :: force_const(:,:,:,:)
  real(dp),     allocatable :: delta_prim(:,:,:,:,:)
  integer,      allocatable :: atom(:,:)
  integer,      allocatable :: no_equiv_ims(:,:,:)
  integer,      allocatable :: atom_in_prim(:)
  integer,      allocatable :: prim_cell_for_atom(:)
  logical,      allocatable :: defined(:,:,:,:)
  
  call cpu_time(t1)
  
  call initialise(prog_function,tol,tol2,delta,structure,structure_sc,atoms, &
     & displacements,forces,temperature,no_kspace_lines,disp_kpoints, &
     & length_scale,vol_scale,fc_scale,&
     & no_prim_cells,               &
     & force_const,delta_prim,    &
     & atom,no_equiv_ims,atom_in_prim,prim_cell_for_atom,defined)
  
    
  write(*,*)'The phonon frequencies at the supercell G vectors will be &
    &calculated.'
  write(*,*)'Calculating the frequencies and displacement patterns on the &
    &G-vector grid.'
  call evaluate_freqs_on_grid(no_prim_cells,         &
     & structure_sc,temperature,structure, &
     & delta_prim,atom,no_equiv_ims,force_const,         &
     & prim_cell_for_atom,atom_in_prim,                                &
     & kpairs_filename,freq_grids_filename,disp_patterns_filename,     &
     & kdisp_patterns_filename,pol_vec_filename,gvectors_filename,     &
     & gvectors_frac_filename,error_filename)
  write(*,*)'Done.  Frequencies and displacement patterns calculated.'
  write(*,*)
  call write_dynamical_matrix(structure_sc,structure,  &
     & no_prim_cells,delta_prim,atom,no_equiv_ims, &
     &force_const,                                                  &
     & dyn_mat_fileroot)
  call write_atoms_in_primitive_cell(structure_sc,structure, &
    & atom,atoms_in_primitive_cell_filename)
  
  call cpu_time(t2)
  
  write(*,*)'Program finished.  Time taken: ',t2-t1
  write(*,*)
end subroutine
end module
