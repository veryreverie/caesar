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
! 2016     Integrated into Caesar. See Caesar git history.

module lte_module
  use constants,      only : dp, kB_au_per_K
  use utils,          only : i2s, errstop
  use min_images,     only : min_images_brute_force, maxim
  use file_module
  implicit none
  
  private
  public :: lte_1,lte_2,lte_3,lte_4
contains

! ----------------------------------------------------------------------
! Construct the dynamical matrix for a given k vector.
! ----------------------------------------------------------------------
subroutine construct_dyn_matrix(kvec,structure,structure_sc, &
   & atom,force_constants,dyn_mat)
  use structure_module
  implicit none
  
  real(dp), intent(in) :: kvec(3)
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: structure_sc
  integer,  intent(in) :: atom(:,:)
  real(dp), intent(in) :: force_constants(:,:,:,:)
  
  complex(dp), intent(inout) :: dyn_mat(structure%no_modes,structure%no_modes)
  
  real(dp) :: delta_r_ims(3,maxim),delta_r_corr(3)
  integer,  allocatable :: no_equiv_ims(:,:,:)
  real(dp), allocatable :: delta_prim(:,:,:,:,:)
  
  integer :: p,n,m,i,j,index1,index2,atom1,im
  complex(dp) :: dm,tempc,expikdotD(structure_sc%supercell%sc_size,structure%no_atoms,&
    &structure%no_atoms)
  real(dp) :: k_dot_D

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
subroutine calculate_eigenfreqs(kvec,structure, &
   & structure_sc,atom,force_constants,omega)
  use linear_algebra
  use structure_module
  implicit none
  
  real(dp), intent(in) :: kvec(3)
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: structure_sc
  integer,  intent(in) :: atom(:,:)
  real(dp), intent(in) :: force_constants(:,:,:,:)
  
  real(dp),intent(out) :: omega(structure%no_modes)
  
  complex(dp) :: dyn_mat(structure%no_modes,structure%no_modes)
  real(dp) :: minusomegasq(structure%no_modes)
  integer :: n,m
  
  type(ComplexEigenstuff) :: estuff
  
  ! Construct dynamical matrix.
  call construct_dyn_matrix(kvec,structure,structure_sc, &
     & atom,force_constants,dyn_mat)
  
  estuff = calculate_eigenstuff(dyn_mat)
  minusomegasq = estuff%evals

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
subroutine calculate_eigenfreqs_and_vecs(kvec,structure, &
   & structure_sc,atom,force_constants,    &
   & omega,pol_vec)
  use linear_algebra
  use structure_module
  implicit none
  
  real(dp), intent(in) :: kvec(3)
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: structure_sc
  integer,  intent(in) :: atom(:,:)
  real(dp), intent(in) :: force_constants(:,:,:,:)
  
  real(dp),    intent(out) :: omega(structure%no_modes)
  complex(dp), intent(out) :: pol_vec(structure%no_modes,structure%no_modes)
  
  real(dp) :: minusomegasq(structure%no_modes)
  complex(dp) :: dyn_mat(structure%no_modes,structure%no_modes)
  integer :: n,m
  
  type(ComplexEigenstuff) :: estuff

  ! Construct dynamical matrix.
  call construct_dyn_matrix(kvec,structure,structure_sc, &
     & atom,force_constants,dyn_mat)
  
  estuff = calculate_eigenstuff(dyn_mat)
  minusomegasq = estuff%evals
  dyn_mat = estuff%evecs

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
subroutine calculate_freq_dos(structure,        &
   & structure_sc,atom,force_constants,&
   & freq_dos_filename,bin_width,freq_dos)
  use constants,      only : max_bin, no_samples, no_fdos_sets, pi
  use linear_algebra, only : dscal
  use rand_no_gen,    only : ranx
  use string_module
  use structure_module
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: structure_sc
  integer,  intent(in) :: atom(:,:)
  real(dp), intent(in) :: force_constants(:,:,:,:)
  
  real(dp), parameter :: tol = 1.0e-5_dp
  
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
    kvec = 2*pi*( ranx()*structure%recip_lattice(1,:) &
              & + ranx()*structure%recip_lattice(2,:) &
              & + ranx()*structure%recip_lattice(3,:))
    call calculate_eigenfreqs(kvec,structure, &
      & structure_sc,atom,force_constants,omega)
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
      kvec = 2*pi*( ranx()*structure%recip_lattice(1,:) &
                & + ranx()*structure%recip_lattice(2,:) &
                & + ranx()*structure%recip_lattice(3,:))
      call calculate_eigenfreqs(kvec,structure, &
        & structure_sc,atom,force_constants,omega)
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
subroutine generate_disp_curve(structure,no_kspace_lines, &
   & disp_kpoints,structure_sc,atom,    &
   & force_constants,                                                      &
   & dispersion_curve_filename)
  use string_module
  use structure_module
  implicit none
  
  type(StructureData), intent(in) :: structure
  integer,  intent(in) :: no_kspace_lines
  real(dp), intent(in) :: disp_kpoints(:,:)
  type(StructureData), intent(in) :: structure_sc
  integer,  intent(in) :: atom(:,:)
  real(dp), intent(in) :: force_constants(:,:,:,:)
  
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
      call calculate_eigenfreqs(kvec,structure, &
        & structure_sc,atom,force_constants,omega)
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
subroutine calculate_speed_sound( &
   & structure,structure_sc,atom,force_constants)
  use constants,   only : pi
  use rand_no_gen, only : ranx
  use structure_module
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: structure_sc
  integer,  intent(in) :: atom(:,:)
  real(dp), intent(in) :: force_constants(:,:,:,:)
  
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
      call calculate_eigenfreqs_and_vecs(kvec,structure,   &
        & structure_sc,atom,force_constants, &
        & omega,pol_vec)
      if(ANY(omega<0.d0))then
        write(*,*)'Imaginary frequencies found.'
        write(*,*)'In terms of the primitive reciprocal lattice vectors, &
          &the k-point is:'
        write(*,*) dot_product(kvec,structure%lattice(1,:)), &
                 & dot_product(kvec,structure%lattice(2,:)), &
                 & dot_product(kvec,structure%lattice(3,:))
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
subroutine evaluate_freqs_on_grid(    &
   & structure_sc,temperature,structure, &
   & atom,force_constants,         &
   & prim_cell_for_atom,atom_in_prim,                                &
   & kpairs_filename,freq_grids_filename,disp_patterns_filename,     &
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
  real(dp), intent(in) :: force_constants(:,:,:,:)
  integer,  intent(in) :: prim_cell_for_atom(:)
  integer,  intent(in) :: atom_in_prim(:)
  
  ! File names
  type(String), intent(in) :: kpairs_filename
  type(String), intent(in) :: freq_grids_filename
  type(String), intent(in) :: disp_patterns_filename
  type(String), intent(in) :: kdisp_patterns_filename
  type(String), intent(in) :: pol_vec_filename
  type(String), intent(in) :: error_filename
  
  ! File units
  integer :: kpairs_file
  integer :: freq_grids_file
  integer :: disp_patterns_file
  integer :: kdisp_patterns_file
  integer :: pol_vec_file
  integer :: error_file
  
  integer :: i,j,ig,index1,index2,p,n,atom1
  logical :: soft_modes
  real(dp) :: gvec(3,structure_sc%supercell%sc_size),R0(3), &
    &omega(structure%no_modes),E,F,rec_root_mass(structure%no_atoms),GdotR, &
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
                           & structure_sc%supercell%gvectors(:,i))
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
  
  kpairs_file = open_write_file(kpairs_filename)
  do i=1,structure_sc%supercell%sc_size
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

  ! Evaluate the frequencies at each supercell G vector.
  E=0.d0  ;  F=0.d0
  soft_modes=.FALSE.
  R0=structure_sc%atoms(:,atom(1,1))
  do ig=1,structure_sc%supercell%sc_size
    if(reference(ig)==0 .or. reference(ig)>ig)then
      call calculate_eigenfreqs_and_vecs(gvec(:,ig),structure, &
        & structure_sc,atom,   &
        & force_constants,omega,pol_vec)
    else
      call calculate_eigenfreqs_and_vecs(gvec(:,reference(ig)), &
        & structure,                                                  &
        & structure_sc,atom,  &
        & force_constants,omega,pol_vec)
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
subroutine write_dynamical_matrix(structure_sc,structure, &
   & atom,      &
   & force_constants,                                                      &
   & dyn_mat_fileroot)
  use constants, only : pi
  use string_module
  use structure_module
  implicit none
  
  type(StructureData), intent(in) :: structure_sc
  type(StructureData), intent(in) :: structure
  integer,  intent(in) :: atom(:,:)
  real(dp), intent(in) :: force_constants(:,:,:,:)
  
  ! File root (file name minus ending)
  type(String), intent(in) :: dyn_mat_fileroot
  
  ! File unit
  integer :: dyn_mat_file
  
  integer :: ig,atom1,cart1,index1,atom2,cart2,index2
  real(dp) :: gvec(3)
  complex(dp) :: dyn_mat(structure%no_modes,structure%no_modes)

  do ig=1,structure_sc%supercell%sc_size
    gvec = 2*pi*matmul( structure_sc%recip_lattice, &
                      & structure_sc%supercell%gvectors(:,ig))
    dyn_mat_file = open_write_file(dyn_mat_fileroot//'.'//ig//'.dat')
    call construct_dyn_matrix(gvec,structure, &
      & structure_sc,atom,force_constants,dyn_mat)
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
! Set up lte
! ----------------------------------------------------------------------
subroutine initialise(structure,structure_sc, &
   & atom,atom_in_prim,prim_cell_for_atom)
  use utils,     only : errstop, wordwrap
  use string_module
  use structure_module
  implicit none
  
  ! ----------------------------------------
  ! Inputs
  ! ----------------------------------------
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: structure_sc
  
  ! ----------------------------------------
  ! previously global variables
  ! ----------------------------------------
  integer,  allocatable, intent(out) :: atom(:,:)
  integer,  allocatable, intent(out) :: atom_in_prim(:)
  integer,  allocatable, intent(out) :: prim_cell_for_atom(:)
  
  integer :: atom_id
  integer :: i,j
  
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
end subroutine

! ----------------------------------------------------------------------
! Main program. (Split into four modes)
! ----------------------------------------------------------------------
subroutine lte_1(structure,structure_sc,force_constants, &
   & temperature,freq_dos_filename, &
   & tdependence1_filename,tdependence2_filename)
  use constants, only : dp
  use utils,     only : errstop, wordwrap
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
  type(String), intent(in) :: freq_dos_filename
  type(String), intent(in) :: tdependence1_filename
  type(String), intent(in) :: tdependence2_filename
  
  ! ----------------------------------------
  ! previously global variables
  ! ----------------------------------------
  real(dp) :: bin_width
  
  real(dp),     allocatable :: freq_dos(:,:)
  integer,      allocatable :: atom(:,:)
  integer,      allocatable :: atom_in_prim(:)
  integer,      allocatable :: prim_cell_for_atom(:)
  
  call initialise(structure,structure_sc, &
     & atom,atom_in_prim,prim_cell_for_atom)

  write(*,*)'Temperature (K)                    :',temperature
  if(temperature<0.d0)call errstop('READ_LTE', &
    &'Temperature should be non-negative.')
  if(temperature<=0.d0)write(*,*)'(i.e. the zero-point energy is to be &
    &calculated.)'
  write(*,*)

  write(*,*)'The mean thermal energy and the free energy will &
    &be calculated.'
  write(*,*)'Calculating the frequency density-of-states function...'
  call calculate_freq_dos(structure,   &
     & structure_sc,atom, &
     & force_constants,                                                 &
     & freq_dos_filename,bin_width,freq_dos)
  call wordwrap('Done.  Frequency density-of-states function written to &
    &freq_dos.dat.  (Please view this file using XMGrace.)')
  write(*,*)

  write(*,*)'Calculating the lattice thermal energy (LTE) and free energy &
    &(LTFE)...'
  call calc_lte(bin_width,temperature,freq_dos,tdependence1_filename)
  call calc_ltfe(bin_width,temperature,freq_dos,tdependence2_filename)
  write(*,*)
end subroutine

subroutine lte_2(structure,structure_sc,force_constants, &
   & no_kspace_lines,disp_kpoints, &
   & dispersion_curve_filename)
  use constants, only : dp
  use utils,     only : errstop, wordwrap
  use string_module
  use structure_module
  implicit none
  
  ! ----------------------------------------
  ! Inputs
  ! ----------------------------------------
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: structure_sc
  real(dp),            intent(in) :: force_constants(:,:,:,:)
  integer,             intent(in) :: no_kspace_lines
  real(dp),            intent(in) :: disp_kpoints(:,:)
  
  ! ----------------------------------------
  ! filenames
  ! ----------------------------------------
  type(String), intent(in) :: dispersion_curve_filename
  
  ! ----------------------------------------
  ! previously global variables
  ! ----------------------------------------
  
  integer,      allocatable :: atom(:,:)
  integer,      allocatable :: atom_in_prim(:)
  integer,      allocatable :: prim_cell_for_atom(:)
  
  integer :: i
  
  call initialise(structure,structure_sc, &
     & atom,atom_in_prim,prim_cell_for_atom)

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
  
  write(*,*)'A dispersion curve will be calculated.'
  write(*,*)'Calculating the requested dispersion curve.'
  call generate_disp_curve(structure,no_kspace_lines,    &
     & disp_kpoints,structure_sc,atom, &
     & force_constants,                                                   &
     & dispersion_curve_filename)
  call wordwrap('Done.  dispersion_curve.dat has been generated.  (Please &
    &view this file using XMGrace.)')
  write(*,*)
end subroutine

subroutine lte_3(structure,structure_sc,force_constants)
  use constants, only : dp
  use utils,     only : errstop, wordwrap
  use string_module
  use structure_module
  implicit none
  
  ! ----------------------------------------
  ! Inputs
  ! ----------------------------------------
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: structure_sc
  real(dp),            intent(in) :: force_constants(:,:,:,:)
  
  ! ----------------------------------------
  ! previously global variables
  ! ----------------------------------------
  
  integer,      allocatable :: atom(:,:)
  integer,      allocatable :: atom_in_prim(:)
  integer,      allocatable :: prim_cell_for_atom(:)
  
  call initialise(structure,structure_sc, &
     & atom,atom_in_prim,prim_cell_for_atom)

  write(*,*)'The speed of sound will be calculated.'
  write(*,*)'Calculating the speed of sound.'
  call calculate_speed_sound( &
     & structure,structure_sc,atom,force_constants)
  write(*,*)'Done.  Speed of sound calculated.'
  write(*,*)
end subroutine

subroutine lte_4(structure,structure_sc,force_constants, &
   & temperature, &
   & kpairs_filename,freq_grids_filename,disp_patterns_filename,            &
   & kdisp_patterns_filename,pol_vec_filename,            &
   & error_filename,dyn_mat_fileroot)
  use constants, only : dp
  use utils,     only : errstop, wordwrap
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
  type(String), intent(in) :: kpairs_filename
  type(String), intent(in) :: freq_grids_filename
  type(String), intent(in) :: disp_patterns_filename
  type(String), intent(in) :: kdisp_patterns_filename
  type(String), intent(in) :: pol_vec_filename
  type(String), intent(in) :: error_filename
  type(String), intent(in) :: dyn_mat_fileroot ! will have *.dat appended
  
  ! ----------------------------------------
  ! previously global variables
  ! ----------------------------------------
  
  integer,      allocatable :: atom(:,:)
  integer,      allocatable :: atom_in_prim(:)
  integer,      allocatable :: prim_cell_for_atom(:)
  
  call initialise(structure,structure_sc, &
     & atom,atom_in_prim,prim_cell_for_atom)
    
  write(*,*)'The phonon frequencies at the supercell G vectors will be &
    &calculated.'
  write(*,*)'Calculating the frequencies and displacement patterns on the &
    &G-vector grid.'
  call evaluate_freqs_on_grid(        &
     & structure_sc,temperature,structure, &
     & atom,force_constants,         &
     & prim_cell_for_atom,atom_in_prim,                                &
     & kpairs_filename,freq_grids_filename,disp_patterns_filename,     &
     & kdisp_patterns_filename,pol_vec_filename,     &
     & error_filename)
  write(*,*)'Done.  Frequencies and displacement patterns calculated.'
  write(*,*)
  call write_dynamical_matrix(structure_sc,structure,  &
     & atom, &
     &force_constants,                                                  &
     & dyn_mat_fileroot)
end subroutine
end module
