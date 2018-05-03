! ======================================================================
! Calculates the vibrational self-consistent field (VSCF) approximation
!    to the phonon Hamiltonian
! ======================================================================
module vscf_module
  use common_module
  
  use scf_old_module
  implicit none
contains

! ----------------------------------------------------------------------
! vibrational self-consistent field main loop
! ----------------------------------------------------------------------
subroutine vscf()
  implicit none
  
  integer :: first_mode, last_mode, no_modes, no_indep_params, &
     & no_coupled_params
  type(String) :: functional
  real(dp),allocatable :: indep_params(:,:), coupled_params(:,:,:)
  real(dp),allocatable :: indep_pot(:,:), indep_pot_bare(:,:), &
     &coupled_pot(:,:,:,:), scf_pot(:,:)
  real(dp),allocatable :: max_amplitude(:),harmonic_freq(:),qs(:,:),dq(:)
  real(dp),allocatable :: shift_harmonic_freq(:,:)
  real(dp),allocatable :: eigenvalues(:,:), eigenvalues_old(:,:), &
     &eigenvectors(:,:,:), basis(:,:,:)
  real(dp),allocatable :: shift_eigenvalues(:,:)
  real(dp),allocatable :: omega(:),x0(:),v0(:)
  real(dp) :: density
  real(dp) :: amp_ratio
  real(dp),allocatable :: root_two_over_i(:),root_i_over_two(:)
  integer :: nbasis, nstates
  real(dp) :: bfp
  integer :: integration_points
  real(dp) :: q
  integer,parameter :: min_t=5,max_t=1000
  real(dp) :: max_diff
  real(dp),allocatable :: scf_indep(:)
  real(dp),allocatable :: scf_indep2(:,:),scf_coupled2(:,:,:,:)
  integer,allocatable :: state(:)
  integer :: alpha,beta,delta,epsilon1
  real(dp) :: int_temp1,int_temp2
  real(dp) :: final_energy,final_energy0
  real(dp) :: tol=1.d-8,tolerance=1.d-18
  real(dp) :: mp2
  real(dp),allocatable :: single_integral(:,:,:),coupled_integral(:,:,:,:,:,:)
  real(dp),allocatable :: partial_coupled_integral(:,:,:,:,:,:), &
     & temp_coupled_integral(:,:,:,:,:,:), &
     & temp_partial_coupled_integral(:,:,:,:,:)
  integer,allocatable :: symmetry_mode(:),symmetry_ref(:)
  integer :: no_unit_cells,execution_mode
  real(dp) :: harmonic_energy,anharmonic_energy,mp2_energy
  real(dp),allocatable :: vec_11(:),vec_21(:),vec_22(:),vec_23(:)
  real(dp) :: temperature, thermal_energy, fenergy, hfenergy, int_energy, &
     & hint_energy, temp_int_energy, temp_hint_energy
  real(dp),allocatable :: part_fn(:),har_part_fn(:)
  
  type(ScfOutput) :: scf_output
  
  ! Input files.
  type(IFile) :: input_file
  type(IFile) :: fit_input_file
  type(IFile) :: symmetry_file
  type(IFile) :: amplitude_ratio_file
  type(IFile) :: amplitude_file
  type(IFile) :: basis_file
  type(IFile) :: vec_in_file
  type(IFile) :: coupled_file
  type(IFile) :: potential_file
  type(IFile) :: coupled_potential_file
  
  ! Output files.
  type(OFile) :: scf_pot_file
  type(OFile) :: convergence_file
  type(OFile) :: output_file
  type(OFile) :: result_file
  type(OFile) :: density_file
  type(OFile) :: eigenvalues_file
  type(OFile) :: independent_mode_file
  type(OFile) :: internal_file
  type(OFile) :: temperature_file
  type(OFile) :: vec_out_file
  
  ! Temporary variables
  integer :: i,j,k,l,ialloc,t,m,n
  integer :: i2,j2
  type(String), allocatable :: line(:)
  
  ! --------------------------------------------------
  ! Read in basis information.
  ! --------------------------------------------------
  input_file = 'vscf_input.dat'
  
  line = split(input_file%line(2))
  first_mode = int(line(1))
  last_mode = int(line(2))
  no_modes = last_mode-first_mode+1
  no_indep_params = int(line(3))
  no_coupled_params = int(line(4))
  
  line = split(input_file%line(4))
  nstates = int(line(1))
  nbasis = int(line(2))
  no_unit_cells = int(line(3))
  temperature = dble(line(4))
  
  line = split(input_file%line(6))
  execution_mode = int(line(1))
  
  ! Read in the functional.
  fit_input_file = 'fit_input.dat'
  line = split(fit_input_file%line(1))
  functional = line(1)
  integration_points = int(line(2)) + 1
  
  if(functional=='polynomial')then
    no_indep_params=no_indep_params-1
  else
    no_indep_params=no_indep_params+1
  endif
  
  ! Write out inputs.
  call print_line('the total number of modes is '//no_modes)
  call print_line('the total number of parameters for the independent &
     &potential is '//no_indep_params)
  call print_line('the total number of parameters for the coupled potential &
     &is '//no_coupled_params)
  if(execution_mode==0)then
    call print_line('calculation with independent modes only.')
  elseif(execution_mode==1)then
    call print_line('calculation with coupled modes.')
  else
    call print_line('execution mode error: invalid input.')
    call err()
  endif
  call print_line('the total number of integration points is '// &
     & integration_points)
  
  ! --------------------------------------------------
  ! Allocate arrays.
  ! --------------------------------------------------
  allocate( basis(no_modes,integration_points,nbasis),                 &
          & max_amplitude(no_modes),                                   &
          & qs(integration_points,no_modes),                           &
          & dq(no_modes),                                              &
          & x0(no_modes),                                              &
          & coupled_params(no_modes,no_modes,no_coupled_params),       &
          & vec_11(no_modes),                                          &
!          & vec_12(no_modes),                                          &
!          & vec_13(no_modes),                                          &
          & vec_21(no_modes),                                          &
          & vec_22(no_modes),                                          &
          & vec_23(no_modes),                                          &
          & part_fn(no_modes),                                         &
          & har_part_fn(no_modes),                                     &
          & indep_params(no_modes,no_indep_params),                    &
          & harmonic_freq(no_modes),                                   &
          & shift_harmonic_freq(no_modes,nbasis),                      &
          & indep_pot(no_modes,integration_points),                    &
          & indep_pot_bare(no_modes,integration_points),               &
          & eigenvalues(no_modes,nbasis),                              &
          & eigenvectors(no_modes,nbasis,nstates),                     &
          & shift_eigenvalues(no_modes,nbasis),                        &
          & eigenvalues_old(no_modes,nbasis),                          &
          & omega(no_modes),                                           &
          & v0(no_modes),                                              &
          & root_two_over_i(nbasis),                                   &
          & root_i_over_two(nbasis),                                   &
          & scf_indep(no_modes),                                       &
          & scf_indep2(no_modes,nstates),                              &
          & scf_coupled2(no_modes,integration_points,nstates,nstates), &
          & state(no_modes),                                           &
          & single_integral(no_modes,nbasis,nbasis),                   &
! &coupled_integral(no_modes,no_modes,nbasis,nbasis,nbasis,nbasis),     &
          & symmetry_mode(no_modes),                                   &
          & symmetry_ref(no_modes),                                    &
          & stat=ialloc); call err(ialloc)

  if(execution_mode==1)then
    allocate( coupled_pot( no_modes, no_modes,                           &
            &              integration_points, integration_points),      &
            & partial_coupled_integral( nbasis, nbasis, nbasis, nbasis,  &
            &                           no_modes, no_modes),             &
            & temp_partial_coupled_integral( nbasis, nbasis,             &
            &                                no_modes, no_modes,         &
            &                                integration_points),        &
            & temp_coupled_integral( nbasis, nbasis, nbasis, nbasis,     &
            &                        no_modes, no_modes),                &
            & coupled_integral( nbasis, nbasis, nbasis, nbasis,          &
            &                   no_modes, no_modes),                     &
            & stat=ialloc); call err(ialloc)
  endif

  ! --------------------------------------------------
  ! Read in data.
  ! --------------------------------------------------
  ! Read in symmetry-related modes.
  symmetry_file = 'symmetry.dat'
  do i=first_mode,last_mode
    i2 = i-first_mode+1
    line = split(symmetry_file%line(i))
    symmetry_mode(i2) = int(line(1))
    symmetry_ref(i2) = int(line(2))
  enddo

  ! Read in maximum amplitude. Calculate sampling points.
  amplitude_ratio_file = 'amplitude_ratio.dat'
  line = split(amplitude_ratio_file%line(1))
  amp_ratio = dble(line(1))/dble(line(2))

  do i=first_mode,last_mode
    i2 = i-first_mode+1
    
    if (symmetry_mode(i2)==symmetry_ref(i2)) then
      amplitude_file = 'max_amplitude.'//i//'.dat'
      line = split(amplitude_file%line(1))
      max_amplitude(i2) = dble(line(1)) * amp_ratio
      harmonic_freq(i2) = dble(line(2))
    else
      max_amplitude(i2) = max_amplitude(symmetry_ref(i2)-(first_mode-1))
      harmonic_freq(i2) = harmonic_freq(symmetry_ref(i2)-(first_mode-1))
    endif
    
    dq(i2) = max_amplitude(i2)*2.0_dp/integration_points
    do j=1,integration_points
      qs(j,i2) = max_amplitude(i2)*(-1+(j-1)*2.0_dp/integration_points)
    enddo
  enddo

  ! Read in eigenvectors
  do i=first_mode,last_mode
    i2 = i-first_mode+1
    if (symmetry_mode(i2)==symmetry_ref(i2)) then
      basis_file = 'functions_basis/basis_'//i//'.dat'
      vec_11(i2) = dble(basis_file%line(1))
      vec_21(i2) = dble(basis_file%line(2))
      vec_22(i2) = dble(basis_file%line(3))
      vec_23(i2) = dble(basis_file%line(4))
      if (vec_23(i2)<1.0e-6_dp) then
        vec_23(i2) = abs(vec_11(i2))
      endif
      
      if (execution_mode==0) then
        do j=1,nstates
          do k=1,nbasis
            if(j==k)then
              eigenvectors(i2,k,j) = 1
            else
              eigenvectors(i2,k,j) = 0
            endif
          enddo
        enddo
      elseif (execution_mode==1) then
        vec_in_file = 'anharmonic_eigenvectors_mode_'//i//'.dat'
        do j=1,nstates
          do k=1,nbasis
            eigenvectors(i2,k,j) = dble( &
               & vec_in_file%line((j-1)*(nbasis+1)+k+3))
          enddo
        enddo
      endif
    else
      do j=1,nstates
        do k=1,nbasis
          eigenvectors(i2,k,j) = eigenvectors(symmetry_ref(i2) &
                             & - (first_mode-1),k,j)
          vec_11(i2) = vec_11(symmetry_ref(i2)-(first_mode-1))
          vec_21(i2) = vec_21(symmetry_ref(i2)-(first_mode-1))
          vec_22(i2) = vec_22(symmetry_ref(i2)-(first_mode-1))
          vec_23(i2) = vec_23(symmetry_ref(i2)-(first_mode-1))
        enddo
      enddo
    endif
  enddo

  ! read in independent potential parameters
  !do i=first_mode,last_mode
  !  i2 = i-first_mode+1
  !  if(symmetry_mode(i2)==symmetry_ref(i2))then
  !    indep_file = 'indep_fit_parameters_'//i//'.dat'
  !    do j=1,no_indep_params
  !      indep_params(i2,j) = int(indep_file%line(j))
  !    enddo
  !  else
  !    do j=1,no_indep_params
  !      indep_params(i2,j)=indep_params(symmetry_ref(i2)-(first_mode-1),j)
  !    enddo
  !  endif
  !enddo

  ! Read in coupled potential parameters.
  if (execution_mode==1) then
    coupled_params = 0
    do i=first_mode,last_mode
      i2 = i-first_mode+1
      if (symmetry_mode(i2)==symmetry_ref(i2)) then
        do j=first_mode,last_mode
          j2 = j-first_mode+1
          if(j>i.or.(j<i.and.symmetry_mode(j2)/=symmetry_ref(j2)))then
            if (file_exists('coupled_fit_params_'//i//'.'//j//'.dat')) then
              coupled_file = 'coupled_fit_params_'//i//'.'//'.dat'
              do k=1,no_coupled_params
                if (j>i) then
                  coupled_params(i2,j2,k) = dble(coupled_file%line(k))
                else
                  coupled_params(j2,i2,k) = dble(coupled_file%line(k))
                endif
              enddo
            endif
          endif
        enddo
      !else
      !  do j=i+1,last_mode
      !    j2 = j-first_mode+1
      !    do k=1,no_coupled_params
      !      coupled_params(i2,j2,k)=&
      !       &coupled_params(symmetry_ref(i2)-(first_mode-1),j2,k)
      !    enddo ! k
      !  enddo ! j 
      endif
    enddo
  endif
  
  !call print_line('coupled parameters')
  !call print_line(coupled_params(1,2,1))
  !call print_line(coupled_params(10,18,1))
  !call print_line(coupled_params(12,21,1))
  !call print_line(coupled_params(21,29,1))
  
  ! Read in basis function defining potential.
  do i=first_mode,last_mode
    i2 = i-first_mode+1
    if(symmetry_mode(i2)==symmetry_ref(i2))then
      basis_file = 'functions_basis/basis_'//i//'.dat'
      v0(i2) = dble(basis_file%line(2))
      x0(i2) = dble(basis_file%line(3))
      omega(i2) = dble(basis_file%line(4))
      if (omega(i2)<1.0e-5_dp) then
        omega(i2) = dble(basis_file%line(1))
      endif
    else
      v0(i2)=v0(symmetry_ref(i2)-(first_mode-1))
      x0(i2)=x0(symmetry_ref(i2)-(first_mode-1))
      omega(i2)=omega(symmetry_ref(i2)-(first_mode-1))
    endif
  enddo

  ! Calculate basis functions.
  ! Harmonic eigenbasis wavefunctions.
  call print_line('calculating basis functions...')
  do i=first_mode,last_mode
    i2 = i-first_mode+1
    bfp=(omega(i2)/PI)**0.25_dp
    
    do j=1,nbasis
      root_i_over_two(j)=sqrt(0.5_dp*j)
      root_two_over_i(j)=1/root_i_over_two(j)
    enddo ! j
    
    do j=1,integration_points
      q = qs(j,i2)
      basis(i2,j,1)=bfp*exp(-0.5_dp*q*q*omega(i2))
      basis(i2,j,2)=sqrt(omega(i2))*q*root_two_over_i(1)*basis(i2,j,1)
      do k=3,nbasis
        basis(i2,j,k)=root_two_over_i(k-1)*(sqrt(omega(i2))*q*basis(i2,j,k-1) &
         &-root_i_over_two(k-2)*basis(i2,j,k-2))
      enddo
    enddo
  enddo
  call print_line('done.')
  call print_line('')

  ! Calculate potentials.
  !ortuzar_file = 'ortuzar.dat'
  call print_line('reading in potentials...')
!    do i=1,no_modes
!      do j=1,integration_points
!        ! anharmonic potential.  nb, instead of shifting basis functions and
!        ! integration grid, we shift the potential.
!        q = qs(j,i)
!        indep_pot(i,j) = v_indep(indep_params(i,:), q-x0(i)) &
!                     & - 0.5_dp*omega(i)*omega(i)*q*q-v0(i)
!        !if(i==9)call print_line(q,indep_pot(i,j))
!      enddo ! j
!    !  call ortuzar_file%print_line('')
!    !  call ortuzar_file%print_line('')
!    enddo ! i

  do i=1,no_modes
    if(symmetry_mode(i)==symmetry_ref(i))then
      l=i+(first_mode-1)
      potential_file = 'indep_potential/indep_pot_'//l//'.dat'
      do j=1,integration_points
        q = qs(j,i)
        line = split(potential_file%line(j))
        indep_pot_bare(i,j) = dble(line(2))
        indep_pot(i,j) = indep_pot_bare(i,j) &
                     & - 0.5_dp*omega(i)*omega(i)*q*q-v0(i)
      enddo
    else 
      do j=1,integration_points
        indep_pot(i,j)=indep_pot(symmetry_ref(i)-(first_mode-1),j)
      enddo
    endif
  enddo

  ! read in coupled modes potential
  if(execution_mode==1)then
    do i=1,no_modes
      if(symmetry_mode(i)==symmetry_ref(i))then
        l=i+3
        do j=1,no_modes
          if(j>i.or.(j<i.and.symmetry_mode(j)/=symmetry_ref(j)))then
            k=j+3
            if (file_exists('coupled_potential/coupled_pot_'//l//'.'//k//'.dat')) then
              coupled_potential_file = 'coupled_potential/coupled_pot_'//l//'.'//k//'.dat'
              do m=1,integration_points
                do n=1,integration_points
                  line = split(coupled_potential_file%line((m-1)*integration_points+n))
                  coupled_pot(i,j,m,n) = dble(line(1))
                enddo ! n
              enddo ! m
            endif ! file exists
          endif ! j symmetry_mode
        enddo ! j
      endif ! symmetry_mode
    enddo ! i
  endif ! execution_mode
  
  ! calculate coupled modes potential
  !if(execution_mode==1)then
  !  do i=1,no_modes
  !  !  if(symmetry_mode(i)==symmetry_ref(i))then
  !  !    do j=1,i-1
  !  !      if(symmetry_mode(j)/=symmetry_ref(j))then
  !  !        do k=1,integration_points
  !  !          q1 = qs(k,i)
  !  !          do l=1,integration_points
  !  !            q2 = qs(l,j)
  !  !            coupled_pot(j,i,k,l) = v_coupled( coupled_params(j,i,:), &
  !  !                                            & q2+x0(j),              &
  !  !                                            & q1+x0(i))
  !  !          enddo
  !  !        enddo
  !  !      endif
  !  !    enddo 
  !      do j=i+1,no_modes
  !        do k=1,integration_points
  !          q1 = qs(k,i)
  !          do l=1,integration_points
  !            q2 = qs(l,j)
  !            ! nb, no shift of potential because
  !            !    v_{ppa}(q) = sum_i v_i(q_i) + sum_ij v_ij(q_i,q_j),
  !            !    and the shifted potential is v_{ppa}(q-q0) + v0
  !            !    and we add v0 to the independent potential above. 
  !            coupled_pot(i,j,k,l) = v_coupled( coupled_params(i,j), &
  !                                            & q1-x0(i),            &
  !                                            & q2-x0(j)) 
  !          enddo
  !        enddo
  !      enddo
  !  !  endif
  !  enddo
  !endif
  call print_line('done.')
  call print_line('')

  ! state to consider
  !state_file = 'state.dat'
  !do i=1,no_modes
  !  state(i) = int(state_file%line(i))
  !enddo
  state=1
  
  ! --------------------------------------------------
  ! Calculate couplings.
  ! --------------------------------------------------

  ! calculate matrix elements
  call print_line('calculating matrix elements...')

  if(execution_mode==1) then
    temp_partial_coupled_integral=0
    do i=1,no_modes
      do j=i+1,no_modes
        do alpha=1,nbasis
          do beta=1,nbasis
            do l=1,integration_points
              temp_partial_coupled_integral(alpha,beta,i,j,l) = &
                 & sum(   basis(i,:,alpha)                      &
                 &      * basis(i,:,beta)                       &
                 &      * coupled_pot(i,j,:,l))                 &
                 & * dq(i)
            enddo
          enddo
        enddo
      enddo
    enddo

    partial_coupled_integral=0
    do i=1,no_modes
      do j=i+1,no_modes
        do alpha=1,nbasis
          do beta=1,nbasis
            do delta=1,nbasis
              do epsilon1=1,nbasis
                partial_coupled_integral(alpha,beta,delta,epsilon1,i,j) =    &
                   & sum(   basis(j,:,delta)                                 &
                   &      * basis(j,:,epsilon1)                              &
                   &      * temp_partial_coupled_integral(alpha,beta,i,j,:)) &
                   & * dq(j)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  endif

  final_energy = sum( (0.5_dp+(state-1)) * harmonic_freq )
  harmonic_energy=final_energy
  
  ! --------------------------------------------------
  ! Self-consistent field loop 
  ! --------------------------------------------------
  ! Open output files.
  scf_pot_file = 'scf_pot.dat'
  convergence_file = 'convergence.dat'
  output_file = 'caesar.output'
  
  do t=1,max_t
    ! Run scf loop.
    scf_output = scf(basis,coupled_pot,indep_pot,execution_mode,final_energy,&
       & no_modes,integration_points,nbasis,                                 &
       & symmetry_mode,symmetry_ref,dq,state,                                &
       & eigenvectors,omega,v0,coupled_integral)
    
    scf_pot = scf_output%scf_pot
    eigenvalues = scf_output%eigenvalues
    eigenvectors = scf_output%eigenvectors
    final_energy = scf_output%final_energy
    final_energy0 = scf_output%final_energy0
    max_diff = scf_output%max_diff
    
    ! Write out self-consistent potential.
    do k=1,integration_points
      call scf_pot_file%print_line( qs(k,1)        //' '// &
                                  & scf_pot(1,k)   //' '// &
                                  & indep_pot(1,k) //' '// &
                                  & scf_output%scf_lin_coupled(1,k))
    enddo
    call scf_pot_file%print_line('')
    call scf_pot_file%print_line('')
    
    ! Print convergence data.
    call convergence_file%print_line(t//' '//final_energy//' '//max_diff)
    call output_file%print_line(t//' '//final_energy//' '//max_diff)
    
    ! Break condition.
    if ((max_diff<tol.and.t>=min_t).or.execution_mode==0) then
      exit
    elseif (t==max_t) then
      call print_line('the scf calculation did not converge')
      call err()
    endif
  enddo
  
  anharmonic_energy = final_energy
  
  ! --------------------------------------------------
  ! Write out scf results.
  ! --------------------------------------------------
  ! Write out eigenvalues.
  eigenvalues_file = 'anharmonic_eigenvalues.dat'
  do i=1,no_modes
    call eigenvalues_file%print_line( 'eigenvalues for mod '//i)
    do j=1,nbasis
      call eigenvalues_file%print_line( eigenvalues(i,j))
    enddo
    call eigenvalues_file%print_line( '')
    call eigenvalues_file%print_line( '')
  enddo
  
  ! Write out independent modes.
  independent_mode_file = 'indep_mode.dat'
  call independent_mode_file%print_line('independent modes frequnecies')
  do i=1,no_modes
    call independent_mode_file%print_line( i+3              //' '// &
                                         & harmonic_freq(i) //' '// &
                                         & eigenvalues(i,1))
  enddo
  call independent_mode_file%print_line('')
  
  ! Write out eigenvectors.
  do i=1,no_modes
    if(symmetry_mode(i)==symmetry_ref(i))then
      l=i+3
      vec_out_file = 'anharmonic_eigenvectors_mode_'//l//'.dat'
      call vec_out_file%print_line( vec_11(i) //' '// &
                                  & nstates   //' '// &
                                  & nbasis)
      call vec_out_file%print_line( vec_21(i) //' '// &
                                  & vec_22(i) //' '// &
                                  & vec_23(i))
      call vec_out_file%print_line('')
      do j=1,nstates
        do k=1,nbasis
          call vec_out_file%print_line(eigenvectors(i,k,j))
        enddo
        call vec_out_file%print_line('')
      enddo
      call vec_out_file%print_line('')
    endif
  enddo
  
  ! Write out convergence data.
  call output_file%print_line( '')
  call output_file%print_line( 'converged in '//t//' iterations.')
  call output_file%print_line( 'final energy (a.u., ev) '// &
     & final_energy//' '//final_energy*EV_PER_HARTREE)
  call output_file%print_line( 'final energy puc (a.u., ev) '// &
     & final_energy/no_unit_cells                       //' '// &
     & final_energy*EV_PER_HARTREE/no_unit_cells)
  call output_file%print_line( '')
  
  result_file = 'vscf_results.dat'
  call result_file%print_line( 'converged in '//t//' iterations.')
  call result_file%print_line( 'final energy (a.u., ev, cm-1) ' // &
                             & final_energy                //' '// &
                             & final_energy*EV_PER_HARTREE //' '// &
                             & final_energy*ev_per_inverse_cm)
  call result_file%print_line( 'energy difference is '//max_diff)
  
  ! --------------------------------------------------
  ! Perturbation theory.
  ! --------------------------------------------------
  if(execution_mode==1)then
    ! calculate second order perturbation.
    call print_line('calculating second order perturbation...')

    ! Independent terms correction.
    call print_line(' - calculating independent terms...')
    scf_indep2=0
    do i=1,no_modes 
      do j=1,nstates
        if (j/=state(i)) then
          do alpha=1,nbasis
            do beta=1,nbasis
              scf_indep2(i,j) = scf_indep2(i,j)                      &
                            & + eigenvectors(i,alpha,state(i))       &
                            & * sum(   basis(i,:,alpha)              &
                            &        * (indep_pot(i,:)-scf_pot(i,:)) &
                            &        * basis(i,:,beta))              &
                            & * eigenvectors(i,beta,j)               &
                            & * dq(i)
            enddo
          enddo
        endif
      enddo
    enddo
    
    ! Coupled terms correction.
    call print_line(' - calculating coupled terms...')
    scf_coupled2=0
    do i=1,no_modes
      do j=i+1,no_modes
        do m=1,nstates
          do n=1,nstates
            if (m/=state(i) .or. n/=state(j)) then
              do alpha=1,nbasis
                do beta=1,nbasis
                  do delta=1,nbasis
                    do epsilon1=1,nbasis
                      scf_coupled2(i,j,m,n) = scf_coupled2(i,j,m,n) &
                         & + eigenvectors(i,alpha,state(i))         &
                         & * eigenvectors(i,beta,state(i))          &
                         & * eigenvectors(j,delta,m)                &
                         & * eigenvectors(j,epsilon1,n)             &
                         & * coupled_integral(alpha,beta,delta,epsilon1,i,j)
                    enddo
                  enddo
                enddo
              enddo
            endif
          enddo
        enddo
      enddo
    enddo
    
    ! Second order perturbation correction.
    call print_line(' - calculating perturbation correction...')
    mp2=0
    ! Single excited modes.
    do i=1,no_modes
      do j=1,nstates
        if (j/=state(i)) then
          int_temp1=0
          do k=1,i-1
            int_temp1=int_temp1+scf_coupled2(k,i,state(k),j)
          enddo
          do k=i+1,no_modes
            int_temp1=int_temp1+scf_coupled2(i,k,j,state(k))
          enddo 
          
          int_temp2 = eigenvalues(i,j)
          do k=1,no_modes
            if (k/=i) then
              int_temp2=int_temp2+eigenvalues(k,state(k))
            endif
          enddo
          
          mp2 = mp2                         &
            & + (scf_indep2(i,j)+int_temp1) &
            & * (scf_indep2(i,j)+int_temp1) &
            & / (final_energy0-int_temp2)
        endif
      enddo
    enddo
    
    ! Double excited modes.
    do i=1,no_modes
      do j=i+1,no_modes
        do m=1,nstates
          do n=1,nstates
            if(m/=state(i).and.n/=state(j))then
              int_temp1=0
              do l=1,no_modes
                if(l/=i.and.l/=j)then
                  int_temp1=int_temp1+eigenvalues(l,state(l))
                endif
                int_temp1=int_temp1+eigenvalues(i,m)+eigenvalues(j,n)
              enddo
              mp2 = mp2                   &
                & + scf_coupled2(i,j,m,n) &
                & * scf_coupled2(i,j,m,n) &
                & / (final_energy0-int_temp1)
            endif
          enddo
        enddo
      enddo
    enddo
    
    ! Write out perturbation theory results.
    call print_line( '')
    call print_line( '     mp2 correction (a.u., ev, cm-1) '// &
                   & mp2                               //' '// &
                   & mp2*EV_PER_HARTREE                //' '// &
                   & mp2*ev_per_inverse_cm)
    call print_line( '     mp2 energy (a.u., ev, cm-1) '    // &
                   & final_energy+mp2                  //' '// &
                   & (final_energy+mp2)*EV_PER_HARTREE //' '// &
                   & (final_energy+mp2)*ev_per_inverse_cm)
    call print_line( '     mp2 energy puc (a.u., ev, cm-1) '  // &
       & (final_energy+mp2)/no_unit_cells                //' '// &
       & (final_energy+mp2)*EV_PER_HARTREE/no_unit_cells //' '// &
       & (final_energy+mp2)*ev_per_inverse_cm/no_unit_cells)
    
    call result_file%print_line( '')
    call result_file%print_line( 'mp2 energy (a.u., ev, cm-1): '// &
       & final_energy+mp2                                  //' '// &
       & (final_energy+mp2)*EV_PER_HARTREE                 //' '// &
       & (final_energy+mp2)*ev_per_inverse_cm)
    
    mp2_energy=final_energy+mp2

  endif ! execution mode

  ! --------------------------------------------------
  ! Write out energies.
  ! --------------------------------------------------
  call print_line('')
  call print_line('-------summary of results-------')
  call print_line( 'harmonic energy (ev): '// &
                 & harmonic_energy*EV_PER_HARTREE/no_unit_cells)
  call print_line( 'anharmonic energy (ev): '// &
                 & anharmonic_energy*EV_PER_HARTREE/no_unit_cells)
  call output_file%print_line( '')
  call output_file%print_line( '-------summary of results-------')
  call output_file%print_line( 'harmonic energy (ev): '// &
     & harmonic_energy*EV_PER_HARTREE/no_unit_cells)
  call output_file%print_line( 'anharmonic energy (ev): '// &
     & anharmonic_energy*EV_PER_HARTREE/no_unit_cells)
  if (execution_mode==1) then
    call print_line( 'mp2 energy (ev): '// &
                   & mp2_energy*EV_PER_HARTREE/no_unit_cells)
  endif

  ! --------------------------------------------------
  ! Calculate thermodynamic quantities.
  ! --------------------------------------------------
  if(temperature>tol)then
    shift_eigenvalues=0
    shift_harmonic_freq=0
    do i=1,no_modes
      if(eigenvalues(i,1)<0)then
        shift_eigenvalues(i,:) = eigenvalues(i,:) - 2*eigenvalues(i,1)
      !elseif(eigenvalues(i,1)<eigen_tol)then
      !  shift_eigenvalues(i,:) = eigenvalues(i,:) - 2*eigenvalues(i,1)
      else 
        shift_eigenvalues(i,:) = eigenvalues(i,:)
      endif
      
      if(harmonic_freq(i)<0)then
        do j=1,nstates
          shift_harmonic_freq(i,j) = harmonic_freq(i)*((j-1)+0.5_dp) &
                                 & - harmonic_freq(i)
        enddo
      else
        do j=1,nstates
          shift_harmonic_freq(i,j) = harmonic_freq(i)*((j-1)+0.5_dp)
        enddo
      endif
    enddo

    ! Calculate partition function.
    thermal_energy = temperature * KB_IN_AU
    part_fn = sum(exp(-shift_eigenvalues/thermal_energy), 2)
    har_part_fn = sum(exp(-shift_harmonic_freq/thermal_energy), 2)

    ! Calculate internal energy.
    int_energy = 0
    hint_energy = 0
    if (temperature>tolerance) then
      internal_file = 'indep_mode_internal.dat'
      do i=1,no_modes
        temp_int_energy = sum( shift_eigenvalues(i,:)       &
                      &      * exp( -shift_eigenvalues(i,:) &
                      &           / thermal_energy))        &
                      & / part_fn(i)
        if (eigenvalues(i,1)<0) then
          temp_int_energy = temp_int_energy-(-2*eigenvalues(i,1))
        endif
        int_energy = int_energy + temp_int_energy
        
        temp_hint_energy = sum( shift_harmonic_freq(i,:)       &
                       &      * exp( -shift_harmonic_freq(i,:) &
                       &           / thermal_energy))          &
                       & / har_part_fn(i)
        if (harmonic_freq(i)<0) then
          temp_hint_energy = temp_hint_energy-(-harmonic_freq(i))
        endif
        hint_energy = hint_energy + temp_hint_energy
        
        call internal_file%print_line( i+3              //' '// &
                                     & temp_hint_energy //' '// &
                                     & temp_int_energy)
      enddo
    endif

    ! calculate free energy
    fenergy=0
    hfenergy=0
    if(temperature>tolerance)then
      temperature_file = 'indep_mode_temperature.dat'
      do i=1,no_modes
        fenergy=fenergy-thermal_energy*log(part_fn(i))
        hfenergy=hfenergy-thermal_energy*log(har_part_fn(i))
        if(eigenvalues(i,1)<0)then
          fenergy=fenergy-(-2*eigenvalues(i,1))
        endif
        if(harmonic_freq(i)<0)then
          hfenergy=hfenergy-(-harmonic_freq(i))
        endif
        call temperature_file%print_line(                    &
           & -(  thermal_energy*log(part_fn(i))              &
           &   - thermal_energy*log(har_part_fn(i))) //' '// &
           & -thermal_energy*log(part_fn(i))         //' '// &
           & -thermal_energy*log(har_part_fn(i)))
      enddo
    else
      fenergy = sum(eigenvalues(:,1))
      hfenergy = 0.5_dp*sum(harmonic_freq)
    endif

    ! Calculate vibrational density.
    density_file = 'density.dat'
    do i=1,integration_points
      density=0
      do j=1,nstates
        do alpha=1,nbasis
          do beta=1,nbasis
            density = density                               &
                  & + eigenvectors(1,alpha,j)               &
                  & * eigenvectors(1,beta,j)                &
                  & * basis(1,i,alpha)                      &
                  & * basis(1,i,beta)                       &
                  & * exp(-eigenvalues(1,j)/thermal_energy) &
                  & / part_fn(1)
          enddo
        enddo
      enddo
      call density_file%print_line(qs(i,1)//' '//density)
    enddo
    
    ! Write out thermodynamic quantities.
    call output_file%print_line( '')
    call output_file%print_line( &
       & '-------summary of free energy results-------')
    call output_file%print_line( 'temperature (k): '//temperature)
    call output_file%print_line( 'harmonic free energy (ev): '// &
       & hfenergy/no_unit_cells*EV_PER_HARTREE)
    call output_file%print_line( 'anharmonic free energy (ev): '// &
       & fenergy/no_unit_cells*EV_PER_HARTREE)
    call output_file%print_line( '')
    call output_file%print_line( &
       & '-------summary of int energy results-------')
    call output_file%print_line( 'temperature (k): '//temperature)
    call output_file%print_line( 'harmonic internal energy (ev): '// &
       & hint_energy/no_unit_cells*EV_PER_HARTREE)
    call output_file%print_line( 'anharmonic internal energy (ev): '// &
       & int_energy/no_unit_cells*EV_PER_HARTREE)
  endif
end subroutine

! ======================================================================
! Functional forms for the one- and two-body potentials
! ======================================================================

! ----------------------------------------------------------------------
! Single mode potential.
! ----------------------------------------------------------------------
function v_indep(indep_params,q) result(output)
  implicit none
  
  real(dp),intent(in) :: indep_params(:),q
  real(dp) :: output
  integer :: order
  type(String) :: functional
  
  type(IFile) :: fit_input_file
  
  ! Temporary variables.
  integer                   :: j
  type(String), allocatable :: line(:)
  
  fit_input_file = 'fit_input.dat'
  line = split(fit_input_file%line(1))
  functional = line(1)
  order = int(line(2))
    
  output = 0
  do j=1,order-1
    output = output + indep_params(j)*q**(j+1)
  enddo
  
  if(functional=='sine')then
    output = output + indep_params(order)*sin(indep_params(order+1)*q)**2
  elseif (functional/='polynomial') then
    call print_line('error: the functional form of the independent potential &
       &is incorrect!')
    call err()
  endif
end function

! ----------------------------------------------------------------------
! Two-mode coupling potential.
! ----------------------------------------------------------------------
function v_coupled(coupled_params,q1,q2) result(output)
  implicit none
  
  real(dp),intent(in) :: q1,q2,coupled_params(:)
  real(dp) :: output
  
  integer :: no_coupled_params
  
  no_coupled_params = size(coupled_params)

  if(no_coupled_params==1)then
    output = coupled_params(1)*q1**2*q2**2
  elseif(no_coupled_params==2)then
    output = coupled_params(1)*q1*q2 &
         & + coupled_params(2)*q1**2*q2**2
  elseif(no_coupled_params==3)then
    output = coupled_params(1)*q1**2*q2**2 &
         & + coupled_params(2)*q1**2*q2    &
         & + coupled_params(3)*q1*q2**2
  elseif(no_coupled_params==4)then
    output = coupled_params(1)*q1*q2       &
         & + coupled_params(2)*q1**2*q2**2 &
         & + coupled_params(3)*q1**2*q2    &
         & + coupled_params(4)*q1*q2**2
  elseif(no_coupled_params==5)then
    output = coupled_params(1)            &
         & * sin(coupled_params(2)*q1)**2 &
         & * sin(coupled_params(3)*q2)**2 &
         & + coupled_params(4)*q1**3*q2   &
         & + coupled_params(5)*q1*q2**3 
  elseif(no_coupled_params==6)then
    output = coupled_params(1)*q1*q2       &
         & + coupled_params(2)*q1**2*q2    &
         & + coupled_params(3)*q1*q2**2    &
         & + coupled_params(4)*q1**3*q2    &
         & + coupled_params(5)*q1**2*q2**2 &
         & + coupled_params(6)*q1*q2**3
  elseif(no_coupled_params==15)then
    output = coupled_params(1)*q1*q2        &
         & + coupled_params(2)*q1**2*q2     &
         & + coupled_params(3)*q1*q2**2     &
         & + coupled_params(4)*q1**3*q2     &
         & + coupled_params(5)*q1**2*q2**2  &
         & + coupled_params(6)*q1*q2**3     &
         & + coupled_params(7)*q1**4*q2     &
         & + coupled_params(8)*q1**3*q2**2  &
         & + coupled_params(9)*q1**2*q2**3  &
         & + coupled_params(10)*q1*q2**4    &
         & + coupled_params(11)*q1**5*q2    &
         & + coupled_params(12)*q1**4*q2**2 &
         & + coupled_params(13)*q1**3*q2**3 &
         & + coupled_params(14)*q1**2*q2**4 &
         & + coupled_params(15)*q1*q2**5
  elseif(no_coupled_params==45)then
    output = coupled_params(1)*q1*q2        &
         & +(coupled_params(2)*q1**2*q2     &
         & + coupled_params(3)*q1*q2**2)    &
         & +(coupled_params(4)*q1**3*q2     &
         & + coupled_params(5)*q1**2*q2**2  &
         & + coupled_params(6)*q1*q2**3)    &
         & +(coupled_params(7)*q1**4*q2     &
         & + coupled_params(8)*q1**3*q2**2  &
         & + coupled_params(9)*q1**2*q2**3  &
         & + coupled_params(10)*q1*q2**4)   &
         & +(coupled_params(11)*q1**5*q2    &
         & + coupled_params(12)*q1**4*q2**2 &
         & + coupled_params(13)*q1**3*q2**3 &
         & + coupled_params(14)*q1**2*q2**4 &
         & + coupled_params(15)*q1*q2**5)   &
         & +(coupled_params(16)*q1**6*q2    &
         & + coupled_params(17)*q1**5*q2**2 &
         & + coupled_params(18)*q1**4*q2**3 &
         & + coupled_params(19)*q1**3*q2**4 &
         & + coupled_params(20)*q1**2*q2**5 &
         & + coupled_params(21)*q1*q2**6)   &
         & +(coupled_params(22)*q1**7*q2    &
         & + coupled_params(23)*q1**6*q2**2 &
         & + coupled_params(24)*q1**5*q2**3 &
         & + coupled_params(25)*q1**4*q2**4 &
         & + coupled_params(26)*q1**3*q2**5 &
         & + coupled_params(27)*q1**2*q2**6 &
         & + coupled_params(28)*q1*q2**7)   &
         & +(coupled_params(29)*q1**8*q2    &
         & + coupled_params(30)*q1**7*q2**2 &
         & + coupled_params(31)*q1**6*q2**3 &
         & + coupled_params(32)*q1**5*q2**4 &
         & + coupled_params(33)*q1**4*q2**5 &
         & + coupled_params(34)*q1**3*q2**6 &
         & + coupled_params(35)*q1**2*q2**7 &
         & + coupled_params(36)*q1*q2**8)   &
         & +(coupled_params(37)*q1**9*q2    &
         & + coupled_params(38)*q1**8*q2**2 &
         & + coupled_params(39)*q1**7*q2**3 &
         & + coupled_params(40)*q1**6*q2**4 &
         & + coupled_params(41)*q1**5*q2**5 &
         & + coupled_params(42)*q1**4*q2**6 &
         & + coupled_params(43)*q1**3*q2**7 &
         & + coupled_params(44)*q1**2*q2**8 &
         & + coupled_params(45)*q1*q2**9)
  else
    call print_line('error: the functional form of the coupled potential is &
       &incorrect!')
    call err()
  endif
end function
end module
