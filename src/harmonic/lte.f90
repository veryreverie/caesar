! Lattice Thermal Energy, Neil Drummond, 12/2004-1/2005.
! Modified by B. Monserrat

! Calculates the frequencies and patterns of atomic displacement at the 
!    supercell G vectors. This is useful if you are interested in following
!    a mode.

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
! 2016     Integrated into Caesar. See Caesar git history.

module lte_module
  use constants,      only : dp, kB_au_per_K
  use utils,          only : i2s, errstop
  use min_images,     only : min_images_brute_force, maxim
  use file_module
  implicit none
  
  private
  public :: lte
contains

! ----------------------------------------------------------------------
! Construct the dynamical matrix for a given k vector.
! ----------------------------------------------------------------------
function construct_dyn_matrix(kvec,structure,structure_sc, &
   & atom,force_constants) result(dynamical_matrix)
  use structure_module
  implicit none
  
  real(dp),            intent(in) :: kvec(3)
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: structure_sc
  integer,             intent(in) :: atom(:,:)
  real(dp),            intent(in) :: force_constants(:,:,:,:)
  
  complex(dp), allocatable :: dynamical_matrix(:,:)
  
  real(dp) :: delta_r_ims(3,maxim),delta_r_corr(3)
  integer,  allocatable :: no_equiv_ims(:,:,:)
  real(dp), allocatable :: delta_prim(:,:,:,:,:)
  
  integer :: p,n,m,i,j,index1,index2,atom1,im
  complex(dp) :: dm,tempc,expikdotD(structure_sc%supercell%sc_size,structure%no_atoms,&
    &structure%no_atoms)
  real(dp) :: k_dot_D
  
  allocate(dynamical_matrix(structure%no_modes,structure%no_modes))

  ! Work out number of equivalent images and Delta Prim. Lattice Vectors
  ! for pairs of atoms (used in evaluation of force-constant matrix).
  allocate(no_equiv_ims(structure_sc%supercell%sc_size,structure%no_atoms,structure%no_atoms),&
    &delta_prim(3,maxim,structure_sc%supercell%sc_size,structure%no_atoms,structure%no_atoms))
  delta_prim=0.d0
  do n=1,structure%no_atoms
    do m=1,structure%no_atoms
      delta_r_corr = structure_sc%atoms(1:3,atom(1,n)) &
                 & - structure_sc%atoms(1:3,atom(1,m))
      do p=1,structure_sc%supercell%sc_size
        ! Work out min. image distance(s) between atoms (1,n) and (p,m).
        call min_images_brute_force( structure_sc%atoms(:,atom(p,m))  &
                                   & -structure_sc%atoms(:,atom(1,n)),&
                                   & transpose(structure_sc%lattice), &
                                   & delta_r_ims,                     &
                                   & no_equiv_ims(p,m,n))
        ! Turn this into the corresponding difference(s) of latt. vects.
        do im=1,no_equiv_ims(p,m,n)
          delta_prim(1:3,im,p,m,n)=delta_r_ims(1:3,im)+delta_r_corr
        enddo ! im
      enddo ! p
    enddo ! m
  enddo ! n

  ! Precompute exp(-ik.(R-R')) to go in the dynamical matrix.
  do n=1,structure%no_atoms
    atom1=atom(1,n)
    do m=1,structure%no_atoms
      do p=1,structure_sc%supercell%sc_size
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
          do p=1,structure_sc%supercell%sc_size
            dm=dm+force_constants(j,i,atom1,atom(p,m))*expikdotD(p,m,n)
          enddo ! p
          dynamical_matrix(index1,index2)=dm
        enddo ! j
      enddo ! i
    enddo ! m
  enddo ! n

  ! Enforce Hermiticity on the dynamical matrix.
  dynamical_matrix = (dynamical_matrix + conjg(transpose(dynamical_matrix)))/2
end function

! ----------------------------------------------------------------------
! Diagonalise the dynamical matrix.
! This subroutine returns the polarisation vectors as well.
! It is not optimised for speed.
! ----------------------------------------------------------------------
subroutine calculate_eigenfreqs_and_vecs(dynamical_matrix,omega,pol_vec)
  use linear_algebra
  use structure_module
  implicit none
  
  complex(dp), intent(in) :: dynamical_matrix(:,:)
  
  real(dp),    allocatable, intent(out) :: omega(:)
  complex(dp), allocatable, intent(out) :: pol_vec(:,:)
  
  
  complex(dp), allocatable :: dyn_mat(:,:) ! Diagonalised dynamical matrix.
  real(dp),    allocatable :: minusomegasq(:)
  integer :: no_modes
  integer :: n,m
  
  type(ComplexEigenstuff) :: estuff
  
  no_modes = size(dynamical_matrix,1)
  
  allocate(minusomegasq(no_modes))
  allocate(dyn_mat(no_modes,no_modes))
  allocate(omega(no_modes))
  allocate(pol_vec(no_modes,no_modes))
  
  ! Diagonalise dynamical matrix.
  estuff       = calculate_eigenstuff(dynamical_matrix)
  minusomegasq = estuff%evals
  dyn_mat      = estuff%evecs
  
  ! Calculate frequencies and polarisation vectors.
  m=no_modes ! no_modes - n + 1
  do n=1,no_modes
  ! Modified by B. Monserrat to output the correct 'omega' and 'pol_vec'
    if(minusomegasq(m)>=0.d0)then
      omega(n)=-SQRT(minusomegasq(m)) ! Unstable mode.
    else
      omega(n)=SQRT(-minusomegasq(m)) ! Stable mode.
    endif
    pol_vec(:,n) = dyn_mat(:,m)
    m=m-1
  enddo
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
! Evaluate the set of phonon frequencies on the supercell G vectors.
! Average the corresponding energies (for testing purposes).
! Write out the real part of the non-mass-reduced polarisation vector, which
! is the pattern of displacement corresponding to the normal mode.
! ----------------------------------------------------------------------
subroutine evaluate_freqs_on_grid(    &
   & structure_sc,temperature,structure, &
   & atom,         &
   & prim_cell_for_atom,atom_in_prim,frequencies,polarisation_vectors, &
   & freq_grids_filename,disp_patterns_filename,     &
   & kdisp_patterns_filename,pol_vec_filename,     &
   & error_filename)
  use constants, only : pi
  use utils,     only : reduce_interval
  use string_module
  use structure_module
  implicit none
  
  type(StructureData),  intent(in) :: structure_sc
  real(dp), intent(in) :: temperature
  type(StructureData), intent(in) :: structure
  integer,  intent(in) :: atom(:,:)
  integer,  intent(in) :: prim_cell_for_atom(:)
  integer,  intent(in) :: atom_in_prim(:)
  real(dp), intent(in) :: frequencies(:,:)
  complex(dp), intent(in) :: polarisation_vectors(:,:,:)
  
  ! File names
  type(String), intent(in) :: freq_grids_filename
  type(String), intent(in) :: disp_patterns_filename
  type(String), intent(in) :: kdisp_patterns_filename
  type(String), intent(in) :: pol_vec_filename
  type(String), intent(in) :: error_filename
  
  ! File units
  integer :: freq_grids_file
  integer :: disp_patterns_file
  integer :: kdisp_patterns_file
  integer :: pol_vec_file
  integer :: error_file
  
  integer :: i,j,ig,index1,index2,p,n,atom1
  logical :: soft_modes
  real(dp) :: gvec(3,structure_sc%supercell%sc_size),R0(3), &
    &omega(structure%no_modes),E,F,GdotR, &
    &disp_pattern(3),kdisp_pattern(3),tot_disp_patt
  complex(dp) :: pol_vec(structure%no_modes,structure%no_modes), &
    &non_mr_pol_vec(3,structure%no_atoms),expiGdotR(structure_sc%supercell%sc_size), &
    &kpol_vec(3,structure%no_atoms)
  real(dp),parameter :: tol_omega=1.d-6 ! For judging whether modes are soft.
  integer :: reference(structure_sc%supercell%sc_size) ! For equivalent +/- G-points.
  real(dp) :: prefactor
  
  ! Transform G-vectors into reciprocal space.
  do i=1,structure_sc%supercell%sc_size
    gvec(:,i) = 2*pi*matmul( transpose(structure_sc%recip_lattice), &
                           & modulo(structure_sc%supercell%gvectors(:,i),structure_sc%supercell%sc_size))
  enddo
  
  ! Calculate +/- G-vector pairs
  reference=0
  do i=1,structure_sc%supercell%sc_size
    do j=1,i-1
      if (all( structure_sc%supercell%gvectors(:,i) &
           & + structure_sc%supercell%gvectors(:,j)==0)) then
        reference(i)=j 
        reference(j)=i
        exit
      endif
    enddo
  enddo
  
  freq_grids_file = open_write_file(freq_grids_filename)
  disp_patterns_file = open_write_file(disp_patterns_filename)
  kdisp_patterns_file = open_write_file(kdisp_patterns_filename)
  pol_vec_file = open_write_file(pol_vec_filename)

  ! Evaluate the frequencies at each supercell G vector.
  E=0.d0  ;  F=0.d0
  soft_modes=.FALSE.
  R0=structure_sc%atoms(:,atom(1,1))
  do ig=1,structure_sc%supercell%sc_size
    if (reference(ig)==0 .or. reference(ig)>ig) then
      omega = frequencies(:,ig)
      pol_vec = polarisation_vectors(:,:,ig)
    else
      omega = frequencies(:,reference(ig))
      pol_vec = polarisation_vectors(:,:,ig)
    endif

    ! The negative is used because the matrix of force constants is the transpose of
    ! the usual expression in derivations that lead to a positive exponential
    do p=1,structure_sc%supercell%sc_size
      if(reference(ig)==0)then
        GdotR=-DOT_PRODUCT(gvec(:,ig),structure_sc%atoms(1:3,atom(p,1))-R0)
      else
        if(reference(ig)>ig)then
          GdotR=-DOT_PRODUCT(gvec(:,ig),structure_sc%atoms(1:3,atom(p,1))-R0)
        else
          GdotR=-DOT_PRODUCT(gvec(:,reference(ig)),structure_sc%atoms(1:3,atom(p,1))-R0)
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
          non_mr_pol_vec(i,n) = pol_vec(index1,index2) &
                            & / sqrt(structure_sc%mass(atom(1,n)))
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
      call print_line(disp_patterns_file,'G-vector '//ig)
      call print_line(kdisp_patterns_file,'G-vector '//ig)
      call print_line(pol_vec_file,'G-vector '//ig)
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
  E=E/DBLE(structure_sc%supercell%sc_size)  ;  F=F/DBLE(structure_sc%supercell%sc_size)

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
subroutine write_dynamical_matrix(dynamical_matrices,dyn_mat_fileroot)
  use string_module
  implicit none
  
  complex(dp), intent(in) :: dynamical_matrices(:,:,:)
  
  ! File root (file name minus ending)
  type(String), intent(in) :: dyn_mat_fileroot
  
  ! File unit
  integer :: dyn_mat_file
  
  integer :: ig,atom1,cart1,index1,atom2,cart2,index2
  complex(dp), allocatable :: dyn_mat(:,:)
  
  integer :: no_gvectors
  integer :: no_modes
  
  no_gvectors = size(dynamical_matrices,3)
  no_modes = size(dynamical_matrices,2)
  
  allocate(dyn_mat(no_modes,no_modes))

  do ig=1,no_gvectors
    dyn_mat = dynamical_matrices(:,:,ig)
    dyn_mat_file = open_write_file(dyn_mat_fileroot//ig//'.dat')
    atom1=0
    do index1=1,no_modes
      atom2=0
      if(MOD(index1,3)==1)then
        atom1=atom1+1
        cart1=1
      endif ! MOD(index1,3)==1
      do index2=1,no_modes
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
! Set up lte
! ----------------------------------------------------------------------
subroutine initialise(structure,structure_sc,force_constants,atom,atom_in_prim, &
   & prim_cell_for_atom,dynamical_matrices,frequencies,polarisation_vectors)
  use constants, only : pi
  use string_module
  use structure_module
  implicit none
  
  ! Inputs.
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: structure_sc
  real(dp),            intent(in) :: force_constants(:,:,:,:)
  
  ! Atom placement information.
  integer,  allocatable, intent(out) :: atom(:,:)
  integer,  allocatable, intent(out) :: atom_in_prim(:)
  integer,  allocatable, intent(out) :: prim_cell_for_atom(:)
  
  ! Dynamical matrices.
  complex(dp), allocatable, intent(out) :: dynamical_matrices(:,:,:)
  real(dp),    allocatable, intent(out) :: frequencies(:,:)
  complex(dp), allocatable, intent(out) :: polarisation_vectors(:,:,:)
  
  ! Temporary variables.
  integer  :: atom_id
  integer  :: i,j
  real(dp) :: gvec(3)
  
  real(dp),    allocatable :: omega(:)
  complex(dp), allocatable :: pol_vec(:,:)
  
  ! Calculate atom placements.
  allocate(atom(structure_sc%supercell%sc_size,structure%no_atoms))
  allocate(atom_in_prim(structure_sc%no_atoms))
  allocate(prim_cell_for_atom(structure_sc%no_atoms))
  
  do i=1,structure_sc%supercell%sc_size
    do j=1,structure%no_atoms
      atom_id = (j-1)*structure_sc%supercell%sc_size + i
      atom(i,j) = atom_id
      atom_in_prim(atom_id) = j
      prim_cell_for_atom = i
    enddo
  enddo
  
  ! Calculate dynamical matrices, and their eigenstuff.
  allocate(dynamical_matrices( structure%no_modes, &
                             & structure%no_modes, &
                             & structure_sc%supercell%sc_size))
  allocate(frequencies(structure%no_modes,structure_sc%supercell%sc_size))
  allocate(polarisation_vectors( structure%no_modes, &
                               & structure%no_modes, &
                               & structure_sc%supercell%sc_size))
  do i=1,structure_sc%supercell%sc_size
    gvec = 2*pi*matmul( structure_sc%recip_lattice, &
                      & modulo(structure_sc%supercell%gvectors(:,i),structure_sc%supercell%sc_size))
    dynamical_matrices(:,:,i) = construct_dyn_matrix(gvec,structure, &
      & structure_sc,atom,force_constants)
    
    call calculate_eigenfreqs_and_vecs(dynamical_matrices(:,:,i),omega,pol_vec)
    frequencies(:,i) = omega
    polarisation_vectors(:,:,i) = pol_vec
  enddo
end subroutine

subroutine lte(structure,structure_sc,force_constants, &
   & temperature, &
   & freq_grids_filename,disp_patterns_filename,            &
   & kdisp_patterns_filename,pol_vec_filename,            &
   & error_filename,dyn_mat_fileroot)
  use constants, only : dp
  use string_module
  use structure_module
  implicit none
  
  ! ----------------------------------------
  ! Inputs
  ! ----------------------------------------
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: structure_sc
  real(dp),            intent(in) :: force_constants(:,:,:,:)
  real(dp),            intent(in) :: temperature
  
  ! ----------------------------------------
  ! filenames
  ! ----------------------------------------
  type(String), intent(in) :: freq_grids_filename
  type(String), intent(in) :: disp_patterns_filename
  type(String), intent(in) :: kdisp_patterns_filename
  type(String), intent(in) :: pol_vec_filename
  type(String), intent(in) :: error_filename
  type(String), intent(in) :: dyn_mat_fileroot ! will have *.dat appended
  
  ! ----------------------------------------
  ! previously global variables
  ! ----------------------------------------
  integer,     allocatable :: atom(:,:)
  integer,     allocatable :: atom_in_prim(:)
  integer,     allocatable :: prim_cell_for_atom(:)
  complex(dp), allocatable :: dynamical_matrices(:,:,:)
  real(dp),    allocatable :: frequencies(:,:)
  complex(dp), allocatable :: polarisation_vectors(:,:,:)
  
  call initialise(structure,structure_sc,force_constants,atom,atom_in_prim, &
     & prim_cell_for_atom,dynamical_matrices,frequencies,polarisation_vectors)
    
  write(*,*)'The phonon frequencies at the supercell G vectors will be &
    &calculated.'
  write(*,*)'Calculating the frequencies and displacement patterns on the &
    &G-vector grid.'
  call evaluate_freqs_on_grid(        &
     & structure_sc,temperature,structure, &
     & atom,         &
     & prim_cell_for_atom,atom_in_prim,frequencies,polarisation_vectors, &
     & freq_grids_filename,disp_patterns_filename,     &
     & kdisp_patterns_filename,pol_vec_filename,     &
     & error_filename)
  write(*,*)'Done.  Frequencies and displacement patterns calculated.'
  write(*,*)
  call write_dynamical_matrix(dynamical_matrices,dyn_mat_fileroot)
end subroutine
end module
