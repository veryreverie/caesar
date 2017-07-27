! ======================================================================
! The scf loop for vscfCoupled
! ======================================================================
module scf_module
  use string_module
  use io_module
  use constants_module, only : dp
  
  type ScfOutput
    real(dp), allocatable :: scf_pot(:,:)
    real(dp), allocatable :: eigenvalues(:,:)
    real(dp)              :: final_energy0
    real(dp)              :: max_diff
  end type
contains

function scf(basis,coupled_pot,indep_pot,execution_mode,final_energy, &
   & no_modes,integration_points,nbasis,nstates,                      &
   & symmetry_mode,symmetry_ref,qs,dq,state,                          &
   & eigenvectors,omega,v0,coupled_integral) result(output)
  implicit none
  
  real(dp), intent(in) :: basis(:,:,:)
  real(dp), intent(in) :: coupled_pot(:,:,:,:)
  real(dp), intent(in) :: indep_pot(:,:)
  integer,  intent(in) :: execution_mode
  real(dp), intent(inout) :: final_energy
  integer,  intent(in) :: no_modes
  integer,  intent(in) :: integration_points
  integer,  intent(in) :: nbasis
  integer,  intent(in) :: nstates
  integer,  intent(in) :: symmetry_mode(:)
  integer,  intent(in) :: symmetry_ref(:)
  real(dp), intent(in) :: qs(:,:)
  real(dp), intent(in) :: dq(:)
  integer,  intent(in) :: state(:)
  real(dp), intent(inout) :: eigenvectors(:,:,:)
  real(dp), intent(in) :: omega(:)
  real(dp), intent(in) :: v0(:)
  real(dp), intent(in) :: coupled_integral(:,:,:,:,:,:)
  type(ScfOutput)      :: output
  
  ! Outputs.
  real(dp), allocatable :: scf_pot(:,:)
  real(dp), allocatable :: eigenvalues(:,:)
  real(dp)              :: final_energy0
  real(dp)              :: max_diff
  
  ! Working variables.
  real(dp), allocatable :: scf_lin_coupled(:,:)
  real(dp), allocatable :: hamiltonian(:,:)
  real(dp), allocatable :: e(:)
  real(dp), allocatable :: scf_coupled(:,:)
  real(dp)              :: final_energy_pre
  real(dp), allocatable :: eigenvectors_old(:,:,:)
  
  ! Files.
  integer :: scf_pot_file
  
  ! Temporary variables.
  real(dp),allocatable :: work(:)
  integer :: lwork,info
  
  integer :: alpha, beta, delta, epsilon1
  
  integer :: i,j,k,ialloc
  
  allocate( scf_lin_coupled(no_modes,integration_points), &
          & scf_pot(no_modes,integration_points),         &
          & hamiltonian(nbasis,nbasis),                   &
          & eigenvalues(no_modes,nbasis),                 &
          & e(nbasis),                                    &
          & scf_coupled(no_modes,no_modes),               &
          & eigenvectors_old(no_modes,nbasis,nbasis),     &
          & stat=ialloc); call err(ialloc)
  
  ! linear coupled correction
  if(execution_mode==1)then
    scf_lin_coupled=0
    do i=1,no_modes
      if(symmetry_mode(i)==symmetry_ref(i))then
        do j=1,no_modes
          if (j==i) then
            cycle
          endif
          
          do k=1,integration_points
            do alpha=1,nbasis
              do beta=1,nbasis
                ! note here is v_ji, also change integration points order
                scf_lin_coupled(i,k) = scf_lin_coupled(i,k)           &
                                   & + eigenvectors(j,alpha,state(j)) &
                                   & * eigenvectors(j,beta,state(j))  &
                                   & * sum(   basis(j,:,alpha)        &
                                   &        * coupled_pot(j,i,:,k)    &
                                   &        * basis(j,:,beta) )       &
                                   & * dq(j)
             enddo
            enddo
          enddo
        enddo
      else
        scf_lin_coupled(i,:) = scf_lin_coupled(symmetry_ref(i),:)
      endif
    enddo
  endif

  ! effective self-consistent potential
  scf_pot = indep_pot
  if (execution_mode==1) then
    scf_pot = scf_pot + scf_lin_coupled
  endif
  
  scf_pot_file = open_append_file('scf_pot.dat')
  do k=1,integration_points
    call print_line(scf_pot_file, qs(k,1)        //' '// &
                                & scf_pot(1,k)   //' '// &
                                & indep_pot(1,k) //' '// &
                                & scf_lin_coupled(1,k))
  enddo
  call print_line(scf_pot_file, '')
  call print_line(scf_pot_file, '')
  close(scf_pot_file)

  ! construct and diagonalise hamiltonian matrix
  !call print_line(' - constructing and diagonalising hamiltonian matrix...')
  do i=1,no_modes
    ! construct hamiltonian matrix      
    hamiltonian=0
    ! harmonic diagonal part: on diagonal we have sho eigenvalues.  
    ! add on matrix elements of anharmonic potential w.r.t. sho eigenfunctions.
    do j=1,nbasis
      hamiltonian(j,j) = hamiltonian(j,j) + (j-1+0.5_dp)*omega(i) + v0(i)
    enddo
    do j=1,nbasis
      do k=1,nbasis
        hamiltonian(j,k) = hamiltonian(j,k) &
                       & + sum(basis(i,:,j)*basis(i,:,k)*scf_pot(i,:)) &
                       !& * (scf_pot(i,:)-0.5_dp*omega(i)*qs(:,i)*qs(:,i)))&
                       & * dq(i)
      enddo
    enddo
    ! diagonalise hamiltonian matrix
    lwork=3*nbasis-1
    allocate(work(lwork), stat=ialloc); call err(ialloc)
    call dsyev('v','u',nbasis,hamiltonian(1,1),nbasis,e(1),work(1),lwork,info)
    if(info/=0) then
      call print_line('dsyev failed. info='//info//'.')
      call err()
    endif
    deallocate(work)
    ! Write out the eigenvalues and eigenvectors obtained by diagonalisation.
    !eigenvalues_old(i,:)=eigenvalues(i,:)
    !eigenvalues(i,:)=0.3_dp*e(:)+0.7_dp*eigenvalues_old(i,:)
    eigenvalues(i,:)=e(:)
    !if(hamiltonian(1,1)*eigenvectors_old(i,1,1)>0)then
      do j=1,nstates
        do k=1,nbasis
          eigenvectors_old(i,k,j)=eigenvectors(i,k,j)
    !      ! mixing scheme for stable convergence
    !      eigenvectors(i,k,j)=(0.5_dp*hamiltonian(k,j)+0.5_dp*eigenvectors_old(i,k,j)) 
          eigenvectors(i,k,j)=hamiltonian(k,j)
        enddo ! k
      enddo ! j
    !else
    !  do j=1,nstates
    !    do k=1,nbasis
    !      eigenvectors_old(i,k,j)=eigenvectors(i,k,j)
    !      ! mixing scheme for stable convergence
    !      eigenvectors(i,k,j)=(0.5_dp*hamiltonian(k,j)-0.5_dp*eigenvectors_old(i,k,j))
    !    enddo ! k
    !  enddo ! j
    !endif
  enddo ! i 

  ! total final energy
  final_energy_pre=final_energy
  final_energy=0
  final_energy0=0
  do i=1,no_modes
    final_energy=final_energy+eigenvalues(i,state(i))
    final_energy0=final_energy0+eigenvalues(i,state(i))
  enddo

  ! coupling correction
  if(execution_mode==1)then
    scf_coupled=0
    do i=1,no_modes
      do j=i+1,no_modes
        do alpha=1,nbasis
          do beta=1,nbasis
            do delta=1,nbasis
              do epsilon1=1,nbasis
               scf_coupled(i,j) = scf_coupled(i,j)          &
                  & + eigenvectors_old(i,alpha,state(i))    &
                  & * eigenvectors_old(i,beta,state(i))     &
                  & * eigenvectors_old(j,delta,state(j))    &
                  & * eigenvectors_old(j,epsilon1,state(j)) &
                  & * coupled_integral(alpha,beta,delta,epsilon1,i,j)
              enddo ! epsilon1
            enddo ! delta
          enddo ! beta
        enddo ! alpha
      enddo ! j
    enddo ! i
  
    do i=1,no_modes
      do j=i+1,no_modes
        final_energy=final_energy-scf_coupled(i,j)
      enddo
    enddo
    !call print_line(final_energy)
  endif

  max_diff=abs(final_energy-final_energy_pre)
  
  output%scf_pot = scf_pot
  output%eigenvalues = eigenvalues
  output%final_energy0 = final_energy0
  output%max_diff = max_diff
end function
end module
