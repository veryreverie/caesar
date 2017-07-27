! Calculates the vibrational self-consistent field (VSCF) approximation
!    to the phonon Hamiltonian
! Author: B. Monserrat
! Created: 07 February 2012
! Modifications:
!   - 14/02/2012 : added MP2
!   - 29/02/2012 : parallelised coupled matrix element calculation
module vscf_module
  use string_module
  use io_module
  use constants_module, only : dp
contains

! ----------------------------------------------------------------------
! vibrational self-consistent field main loop
! ----------------------------------------------------------------------
subroutine vscf()
  use constants_module, only : pi, ev_per_hartree, kb_in_au, ev_per_inverse_cm
  use scf_module
  !use parallel
  implicit none
  
  integer :: first_mode, last_mode, no_modes, no_indep_params, &
     & no_coupled_params
  type(String) :: functional
  real(dp),allocatable :: indep_params(:,:), coupled_params(:,:,:)
  real(dp) :: mireia,crispin
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
  real(dp) :: q,q1
  integer,parameter :: min_t=5,max_t=1000
  real(dp) :: max_diff
  real(dp),allocatable :: scf_indep(:)
  real(dp),allocatable :: scf_indep2(:,:),scf_coupled2(:,:,:,:)
  integer,allocatable :: state(:)
  integer :: alpha,beta,gamma1,delta,epsilon1
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
  type(String), allocatable :: input_file(:)
  type(String), allocatable :: fit_input_file(:)
  type(String), allocatable :: symmetry_file(:)
  type(String), allocatable :: amplitude_ratio_file(:)
  type(String), allocatable :: amplitude_file(:)
  type(String), allocatable :: basis_file(:)
  type(String), allocatable :: vec_in_file(:)
  type(String), allocatable :: coupled_file(:)
  type(String), allocatable :: potential_file(:)
  type(String), allocatable :: coupled_potential_file(:)
  
  ! Output files.
  integer :: convergence_file
  integer :: output_file
  integer :: result_file
  integer :: density_file
  integer :: eigenvalues_file
  integer :: independent_mode_file
  integer :: internal_file
  integer :: temperature_file
  integer :: vec_out_file
  
  ! Temporary variables
  integer :: i,j,k,l,ialloc,t,m,n
  integer :: i2,j2
  type(String), allocatable :: line(:)

  ! read in the total number of modes and order of p and coupled potentials 
  input_file = read_lines('vscf_input.dat')
  
  line = split(input_file(2))
  first_mode = int(line(1))
  last_mode = int(line(2))
  no_modes = last_mode-first_mode+1
  no_indep_params = int(line(3))
  no_coupled_params = int(line(4))
  
  line = split(input_file(4))
  nstates = int(line(1))
  nbasis = int(line(2))
  no_unit_cells = int(line(3))
  temperature = dble(line(4))
  
  line = split(input_file(6))
  execution_mode = int(line(1))
  
  ! Read in the functional.
  fit_input_file = read_lines('fit_input.dat')
  line = split(fit_input_file(1))
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

  ! allocate various arrays
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
          & scf_pot(no_modes,integration_points),                      &
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

  ! Read in symmetry-related modes.
  symmetry_file = read_lines('symmetry.dat')
  do i=first_mode,last_mode
    i2 = i-first_mode+1
    line = split(symmetry_file(i))
    symmetry_mode(i2) = int(line(1))
    symmetry_ref(i2) = int(line(2))
  enddo

  ! Read in maximum amplitude.
  amplitude_ratio_file = read_lines('amplitude_ratio.dat')
  line = split(amplitude_ratio_file(1))
  amp_ratio = dble(line(1))/dble(line(2))

  do i=first_mode,last_mode
    i2 = i-first_mode+1
    
    if (symmetry_mode(i2)==symmetry_ref(i2)) then
      amplitude_file = read_lines('max_amplitude.'//i//'.dat')
      line = split(amplitude_file(1))
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

  ! read in eigenvectors
  do i=first_mode,last_mode
    i2 = i-first_mode+1
    if (symmetry_mode(i2)==symmetry_ref(i2)) then
      basis_file = read_lines('functions_basis/basis_'//i//'.dat')
      vec_11(i2) = dble(basis_file(1))
      vec_21(i2) = dble(basis_file(2))
      vec_22(i2) = dble(basis_file(3))
      vec_23(i2) = dble(basis_file(4))
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
        vec_in_file = read_lines('anharmonic_eigenvectors_mode_'//i//'.dat')
        do j=1,nstates
          do k=1,nbasis
            eigenvectors(i2,k,j) = dble(vec_in_file((j-1)*(nbasis+1)+k+3))
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
  !    indep_file = read_lines('indep_fit_parameters_'//i//'.dat')
  !    do j=1,no_indep_params
  !      indep_params(i2,j) = int(indep_file(j))
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
              coupled_file = read_lines('coupled_fit_params_'//i//'.'//'.dat')
              do k=1,no_coupled_params
                if (j>i) then
                  coupled_params(i2,j2,k) = dble(coupled_file(k))
                else
                  coupled_params(j2,i2,k) = dble(coupled_file(k))
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
      basis_file = read_lines('functions_basis/basis_'//i//'.dat')
      v0(i2) = dble(basis_file(2))
      x0(i2) = dble(basis_file(3))
      omega(i2) = dble(basis_file(4))
      if (omega(i2)<1.0e-5_dp) then
        omega(i2) = dble(basis_file(1))
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
    bfp=(omega(i2)/pi)**0.25_dp
    
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
  !ortuzar_file = open_write_file('ortuzar.dat')
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
!    !  call print_line(ortuzar_file, '')
!    !  call print_line(ortuzar_file, '')
!    enddo ! i
!    !close(ortuzar_file)

  do i=1,no_modes
    if(symmetry_mode(i)==symmetry_ref(i))then
      l=i+(first_mode-1)
      potential_file = read_lines('indep_potential/indep_pot_'//l//'.dat')
      do j=1,integration_points
        q = qs(j,i)
        line = split(potential_file(j))
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
              coupled_potential_file = read_lines('coupled_potential/coupled_pot_'//l//'.'//k//'.dat')
              do m=1,integration_points
                do n=1,integration_points
                  line = split(coupled_potential_file((m-1)*integration_points+n))
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
  !state_file = read_lines('state.dat')
  !do i=1,no_modes
  !  state(i) = int(state_file(i))
  !enddo ! i
  !close(state_file)
  state=1
  
  ! --------------------------------------------------
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
  
  ! Open output files.
  convergence_file = open_write_file('convergence.dat')
  output_file = open_append_file('caesar.output')
  
  ! --------------------------------------------------
  ! self-consistent field loop 
  ! --------------------------------------------------
  do t=1,max_t
    ! Run scf loop.
    scf_output = scf(basis,coupled_pot,indep_pot,execution_mode,final_energy,&
       & no_modes,integration_points,nbasis,nstates,                         &
       & symmetry_mode,symmetry_ref,qs,dq,state,                             &
       & eigenvectors,omega,v0,coupled_integral)
    
    scf_pot = scf_output%scf_pot
    eigenvalues = scf_output%eigenvalues
    final_energy0 = scf_output%final_energy0
    max_diff = scf_output%max_diff
    
    ! Print convergence data.
    call print_line(convergence_file, t//' '//final_energy//' '//max_diff)
    call print_line(output_file, t//' '//final_energy//' '//max_diff)
    
    ! Break condition.
    if ((max_diff<tol.and.t>=min_t).or.execution_mode==0) then
      exit
    elseif (t==max_t) then
      call print_line('the scf calculation did not converge')
      call err()
    endif
  enddo
  
  close(convergence_file)
  
  ! --------------------------------------------------
  ! After SCF loop.
  ! --------------------------------------------------
  eigenvalues_file = open_write_file('anharmonic_eigenvalues.dat')
  independent_mode_file = open_write_file('indep_mode.dat')
  call print_line(independent_mode_file,'independent modes frequnecies')
  mireia = 0
  do i=1,no_modes
    mireia = mireia + harmonic_freq(i)/2
    crispin = crispin + eigenvalues(i,1)
    call print_line(independent_mode_file, i+3              //' '// &
                                         & harmonic_freq(i) //' '// &
                                         & eigenvalues(i,1))
    call print_line(eigenvalues_file, 'eigenvalues for mod '//i)
    do j=1,nbasis
      call print_line(eigenvalues_file, eigenvalues(i,j))
    enddo
    call print_line(eigenvalues_file, '')
    call print_line(eigenvalues_file, '')
  enddo
  call print_line(independent_mode_file, '')
  close(eigenvalues_file)
  close(independent_mode_file)

  do i=1,no_modes
    if(symmetry_mode(i)==symmetry_ref(i))then
      l=i+3
      vec_out_file = open_write_file('anharmonic_eigenvectors_mode_'//l//'.dat')
      call print_line(vec_out_file, vec_11(i) //' '// &
                                  & nstates   //' '// &
                                  & nbasis)
      call print_line(vec_out_file, vec_21(i) //' '// &
                                  & vec_22(i) //' '// &
                                  & vec_23(i))
      call print_line(vec_out_file, '')
      do j=1,nstates
        do k=1,nbasis
          call print_line(vec_out_file, eigenvectors(i,k,j))
        enddo ! k
        call print_line(vec_out_file, '')
      enddo ! j
      call print_line(vec_out_file, '')
      close(vec_out_file)
    endif ! symmetry
  enddo ! i

  call print_line('')
  call print_line('converged in '//t//' iterations.')
  call print_line('final energy (a.u., ev) '//final_energy//' '// &
     & final_energy*ev_per_hartree)
  call print_line('final energy puc (a.u., ev) '// &
     & final_energy/no_unit_cells          //' '// &
     & final_energy*ev_per_hartree/no_unit_cells)
  call print_line('')
  
  call print_line(output_file, '')
  call print_line(output_file, 'converged in '//t//' iterations.')
  call print_line(output_file, 'final energy (a.u., ev) '// &
     & final_energy//' '//final_energy*ev_per_hartree)
  call print_line(output_file, 'final energy puc (a.u., ev) '// &
     & final_energy/no_unit_cells                       //' '// &
     & final_energy*ev_per_hartree/no_unit_cells)
  call print_line(output_file, '')
  
  result_file = open_write_file('vscf_results.dat')
  call print_line(result_file, 'converged in '//t//' iterations.')
  call print_line(result_file, 'final energy (a.u., ev, cm-1) ' // &
                             & final_energy                //' '// &
                             & final_energy*ev_per_hartree //' '// &
                             & final_energy*ev_per_inverse_cm)
  call print_line(result_file, 'energy difference is '//max_diff)
  
  anharmonic_energy=final_energy
  
  ! --------------------------------------------------
  ! Perturbation theory.
  ! --------------------------------------------------
  if(execution_mode==1)then
    ! calculate second order perturbation
    call print_line('calculating second order perturbation...')

    ! independent terms correction
    call print_line(' - calculating independent terms...')
    scf_indep2=0
    do i=1,no_modes 
      do j=1,state(i)-1
        do alpha=1,nbasis
          do beta=1,nbasis
            int_temp1=0
            do gamma1=1,integration_points
              int_temp1=int_temp1+basis(i,gamma1,alpha)*(indep_pot(i,gamma1)-scf_pot(i,gamma1))*&
               &basis(i,gamma1,beta)*dq(i)
            enddo ! gamma
            scf_indep2(i,j)=scf_indep2(i,j)+eigenvectors(i,alpha,state(i))*eigenvectors(i,beta,j)*&
             &int_temp1
          enddo ! beta
        enddo ! alpha
      enddo ! j
      do j=state(i)+1,nstates
        do alpha=1,nbasis
          do beta=1,nbasis
            int_temp1=0
            do gamma1=1,integration_points
              int_temp1=int_temp1+basis(i,gamma1,alpha)*(indep_pot(i,gamma1)-scf_pot(i,gamma1))*&
               &basis(i,gamma1,beta)*dq(i)
            enddo ! gamma
            scf_indep2(i,j)=scf_indep2(i,j)+eigenvectors(i,alpha,state(i))*eigenvectors(i,beta,j)*int_temp1
          enddo ! beta
        enddo ! alpha
      enddo ! j
    enddo ! i
    
    ! coupled terms correction
    call print_line(' - calculating coupled terms...')
    scf_coupled2=0
    do i=1,no_modes
      do j=i+1,no_modes
        do m=1,nstates
          do n=1,nstates
            if(m/=state(i).or.n/=state(j))then
              do alpha=1,nbasis
                do beta=1,nbasis
                  do delta=1,nbasis
                    do epsilon1=1,nbasis
                      scf_coupled2(i,j,m,n)=scf_coupled2(i,j,m,n)+eigenvectors(i,alpha,state(i))*&
                       &eigenvectors(i,beta,state(i))*eigenvectors(j,delta,m)*&
                       &eigenvectors(j,epsilon1,n)*coupled_integral(alpha,beta,delta,epsilon1,i,j)
                    enddo ! epsilon1
                 enddo ! delta
                enddo ! beta
              enddo ! alpha
            endif ! m,n/=reference state
          enddo ! n 
        enddo ! m
      enddo ! j
    enddo ! i
    
    ! second order perturbation correction
    call print_line(' - calculating perturbation correction...')
    mp2=0
    int_temp1=0
    int_temp2=0
    ! single excited modes
    do i=1,no_modes
      do j=1,nstates
        int_temp1=0
        int_temp2=0
        if(j/=state(i))then
          do k=1,i-1
            int_temp1=int_temp1+scf_coupled2(k,i,state(k),j)
          enddo ! k
          do k=i+1,no_modes
            int_temp1=int_temp1+scf_coupled2(i,k,j,state(k))
          enddo 
          do k=1,i-1
            int_temp2=int_temp2+eigenvalues(k,state(k))
          enddo ! k
          do k=i+1,no_modes
            int_temp2=int_temp2+eigenvalues(k,state(k))
          enddo ! k 
          int_temp2=int_temp2+eigenvalues(i,j)
          mp2=mp2+(scf_indep2(i,j)+int_temp1)*(scf_indep2(i,j)+int_temp1)/(final_energy0-int_temp2)
        endif ! j/=reference state
      enddo ! j 
    enddo ! i
    ! double excited modes
    int_temp1=0
    do i=1,no_modes
      do j=i+1,no_modes
        do m=1,nstates
          do n=1,nstates
            int_temp1=0
            if(m/=state(i).and.n/=state(j))then
              do l=1,no_modes
                if(l/=i.and.l/=j)then
                  int_temp1=int_temp1+eigenvalues(l,state(l))
                endif ! l/= reference state
                int_temp1=int_temp1+eigenvalues(i,m)+eigenvalues(j,n)
              enddo ! l
              mp2=mp2+scf_coupled2(i,j,m,n)*scf_coupled2(i,j,m,n)/(final_energy0-int_temp1)
            endif ! m,n/= reference state
          enddo ! n
        enddo ! m
      enddo ! j
    enddo ! i

    call print_line('')
    call print_line( '     mp2 correction (a.u., ev, cm-1) '// &
                   & mp2                               //' '// &
                   & mp2*ev_per_hartree                //' '// &
                   & mp2*ev_per_inverse_cm)
    call print_line( '     mp2 energy (a.u., ev, cm-1) '    // &
                   & final_energy+mp2                  //' '// &
                   & (final_energy+mp2)*ev_per_hartree //' '// &
                   & (final_energy+mp2)*ev_per_inverse_cm)
    call print_line( '     mp2 energy puc (a.u., ev, cm-1) '  // &
       & (final_energy+mp2)/no_unit_cells                //' '// &
       & (final_energy+mp2)*ev_per_hartree/no_unit_cells //' '// &
       & (final_energy+mp2)*ev_per_inverse_cm/no_unit_cells)
    
    call print_line(result_file, '')
    call print_line(result_file, 'mp2 energy (a.u., ev, cm-1): '// &
       & final_energy+mp2                                  //' '// &
       & (final_energy+mp2)*ev_per_hartree                 //' '// &
       & (final_energy+mp2)*ev_per_inverse_cm)
    
    mp2_energy=final_energy+mp2

  endif ! execution mode
  close(result_file)

  call print_line('')
  call print_line('-------summary of results-------')
  call print_line( 'harmonic energy (ev): '// &
                 & harmonic_energy*ev_per_hartree/no_unit_cells)
  call print_line( 'anharmonic energy (ev): '// &
                 & anharmonic_energy*ev_per_hartree/no_unit_cells)
  call print_line(output_file, '')
  call print_line(output_file, '-------summary of results-------')
  call print_line(output_file, 'harmonic energy (ev): '// &
     & harmonic_energy*ev_per_hartree/no_unit_cells)
  call print_line(output_file, 'anharmonic energy (ev): '// &
     & anharmonic_energy*ev_per_hartree/no_unit_cells)
  if (execution_mode==1) then
    call print_line( 'mp2 energy (ev): '// &
                   & mp2_energy*ev_per_hartree/no_unit_cells)
  endif

  if(temperature>tol)then
    shift_eigenvalues=0
    shift_harmonic_freq=0
    do i=1,no_modes
      if(eigenvalues(i,1)<0)then
        do j=1,nstates
          shift_eigenvalues(i,j)=eigenvalues(i,j)-2*eigenvalues(i,1)
        enddo ! j
      !elseif(eigenvalues(i,1)<eigen_tol)then
      !  do j=1,nstates
      !    shift_eigenvalues(i,j)=eigenvalues(i,j)-2*eigenvalues(i,1)
      !  enddo ! j
      else 
        do j=1,nstates
          shift_eigenvalues(i,j)=eigenvalues(i,j)
        enddo ! j
      endif ! eigenval < 0
      if(i==9)then
        do j=1,nstates
        enddo 
      endif 
      if(harmonic_freq(i)<0)then
        do j=1,nstates
          shift_harmonic_freq(i,j)=harmonic_freq(i)*((j-1)+0.5_dp)-harmonic_freq(i)
        enddo !j
      else
        do j=1,nstates
          shift_harmonic_freq(i,j)=harmonic_freq(i)*((j-1)+0.5_dp)
        enddo !j
      endif ! harmonic_freq < 0
    enddo ! i

    ! calculate partition function
    thermal_energy=temperature * kb_in_au
    part_fn=0
    har_part_fn=0
    do i=1,no_modes
      do j=1,nstates
        part_fn(i)=part_fn(i)+exp(-shift_eigenvalues(i,j)/thermal_energy)
        har_part_fn(i)=har_part_fn(i)+exp(-shift_harmonic_freq(i,j)/thermal_energy)
      enddo ! j 
    enddo ! i

    ! calculate internal energy
    int_energy=0
    hint_energy=0
    internal_file = open_write_file('indep_mode_internal.dat')
    if(temperature>tolerance)then
      do i=1,no_modes
        temp_int_energy=0
        temp_hint_energy=0
        do j=1,nstates
          int_energy = int_energy                                  &
                   & + shift_eigenvalues(i,j)                      &
                   & * exp(-shift_eigenvalues(i,j)/thermal_energy) &
                   & / part_fn(i)
          temp_int_energy = temp_int_energy                            &
                        & + shift_eigenvalues(i,j)                     &
                        & * exp(-shift_eigenvalues(i,j)/thermal_energy)&
                        & / part_fn(i)
          hint_energy = hint_energy                                   &
                    & + shift_harmonic_freq(i,j)                      &
                    & * exp(-shift_harmonic_freq(i,j)/thermal_energy) &
                    & / har_part_fn(i)
          temp_hint_energy = temp_hint_energy               &
                         & + shift_harmonic_freq(i,j)       &
                         & * exp( -shift_harmonic_freq(i,j) &
                         &      / thermal_energy)           &
                         & / har_part_fn(i)
        enddo ! j
        if(eigenvalues(i,1)<0)then
          int_energy=int_energy-(-2*eigenvalues(i,1))
          temp_int_energy=temp_int_energy-(-2*eigenvalues(i,1))
        endif ! eigenval < 0
        if(harmonic_freq(i)<0)then
          hint_energy=hint_energy-(-2*harmonic_freq(i))
          temp_hint_energy=temp_hint_energy-(-harmonic_freq(i))
        endif ! eigenval < 0
        call print_line(internal_file, i+3              //' '// &
                                     & temp_hint_energy //' '// &
                                     & temp_int_energy)
      enddo ! i
    endif ! temp > tolerance
    close(internal_file)

    ! calculate free energy
    fenergy=0
    hfenergy=0
    temperature_file = open_write_file('indep_mode_temperature.dat')
    if(temperature>tolerance)then
      do i=1,no_modes
        fenergy=fenergy-thermal_energy*log(part_fn(i))
        hfenergy=hfenergy-thermal_energy*log(har_part_fn(i))
        if(eigenvalues(i,1)<0)then
          fenergy=fenergy-(-2*eigenvalues(i,1))
        endif ! eigenval < 0
        if(harmonic_freq(i)<0)then
          hfenergy=hfenergy-(-harmonic_freq(i))
        endif ! eigenval < 0
        call print_line(temperature_file,                    &
           & -(  thermal_energy*log(part_fn(i))              &
           &   - thermal_energy*log(har_part_fn(i))) //' '// &
           & -thermal_energy*log(part_fn(i))         //' '// &
           & -thermal_energy*log(har_part_fn(i)))
      enddo ! i
    else
      do i=1,no_modes
        fenergy=fenergy+eigenvalues(i,1)
        hfenergy=hfenergy+0.5_dp*harmonic_freq(i)
      enddo ! i
    endif
    close(temperature_file)

    ! calculate vibrational density
    density_file = open_write_file('density.dat')
    do i=1,integration_points
      density=0
      do j=1,nstates
        do alpha=1,nbasis
          do beta=1,nbasis
            density=density+eigenvectors(1,alpha,j)*eigenvectors(1,beta,j)*&
             &basis(1,i,alpha)*basis(1,i,beta)*exp(-eigenvalues(1,j)/thermal_energy)/part_fn(1)
          enddo ! beta
        enddo ! alpha
      enddo ! j
      q1 = qs(i,1)
      call print_line(density_file, q1//' '//density)
    enddo ! i
    close(density_file)
    
    call print_line(output_file, '')
    call print_line(output_file, &
       & '-------summary of free energy results-------')
    call print_line(output_file, 'temperature (k): '//temperature)
    call print_line(output_file, 'harmonic free energy (ev): '// &
       & hfenergy/no_unit_cells*ev_per_hartree)
    call print_line(output_file, 'anharmonic free energy (ev): '// &
       & fenergy/no_unit_cells*ev_per_hartree)
    call print_line(output_file, '')
    call print_line(output_file, &
       & '-------summary of int energy results-------')
    call print_line(output_file, 'temperature (k): '//temperature)
    call print_line(output_file, 'harmonic internal energy (ev): '// &
       & hint_energy/no_unit_cells*ev_per_hartree)
    call print_line(output_file, 'anharmonic internal energy (ev): '// &
       & int_energy/no_unit_cells*ev_per_hartree)
  endif
      
  close(output_file)
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
  
  type(String), allocatable :: fit_input_file(:)
  
  ! Temporary variables.
  integer                   :: j
  type(String), allocatable :: line(:)
  
  fit_input_file = read_lines('fit_input.dat')
  line = split(fit_input_file(1))
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
