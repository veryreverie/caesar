! ======================================================================
! The scf loop for vscfCoupled
! ======================================================================
module caesar_scf_old_module
  use caesar_common_module
  implicit none
  
  type ScfOutput
    real(dp), allocatable :: scf_pot(:,:)
    real(dp), allocatable :: scf_lin_coupled(:,:)
    real(dp), allocatable :: eigenvalues(:,:)
    real(dp), allocatable :: eigenvectors(:,:,:)
    real(dp)              :: final_energy
    real(dp)              :: final_energy0
    real(dp)              :: max_diff
  end type
contains

function scf(basis,coupled_pot,indep_pot,execution_mode,final_energy_pre, &
   & no_modes,integration_points,nbasis,                                  &
   & symmetry_mode,symmetry_ref,dq,state,                                 &
   & eigenvectors_old,omega,v0,coupled_integral) result(output)
  implicit none
  
  real(dp), intent(in) :: basis(:,:,:)
  real(dp), intent(in) :: coupled_pot(:,:,:,:)
  real(dp), intent(in) :: indep_pot(:,:)
  integer,  intent(in) :: execution_mode
  real(dp), intent(in) :: final_energy_pre
  integer,  intent(in) :: no_modes
  integer,  intent(in) :: integration_points
  integer,  intent(in) :: nbasis
  integer,  intent(in) :: symmetry_mode(:)
  integer,  intent(in) :: symmetry_ref(:)
  real(dp), intent(in) :: dq(:)
  integer,  intent(in) :: state(:)
  real(dp), intent(in) :: eigenvectors_old(:,:,:)
  real(dp), intent(in) :: omega(:)
  real(dp), intent(in) :: v0(:)
  real(dp), intent(in) :: coupled_integral(:,:,:,:,:,:)
  type(ScfOutput)      :: output
  
  ! Outputs.
  real(dp), allocatable :: scf_pot(:,:)
  real(dp), allocatable :: scf_lin_coupled(:,:)
  real(dp), allocatable :: eigenvalues(:,:)
  real(dp), allocatable :: eigenvectors(:,:,:)
  real(dp)              :: final_energy
  real(dp)              :: final_energy0
  real(dp)              :: max_diff
  
  ! Working variables.
  real(dp), allocatable :: hamiltonian(:,:)
  real(dp), allocatable :: scf_coupled(:,:)
  type(RealEigenstuff)  :: eigenstuff
  
  ! Temporary variables.
  integer :: alpha, beta, delta, epsilon1
  integer :: i,j,k,ialloc
  
  allocate( scf_lin_coupled(no_modes,integration_points), &
          & scf_pot(no_modes,integration_points),         &
          & hamiltonian(nbasis,nbasis),                   &
          & eigenvalues(no_modes,nbasis),                 &
          & eigenvectors(no_modes,nbasis,nbasis),         &
          & scf_coupled(no_modes,no_modes),               &
          & stat=ialloc); call err(ialloc)
  
  ! linear coupled correction
  if(execution_mode==1)then
    scf_lin_coupled=0
    do i=1,no_modes
      if(symmetry_mode(i)==symmetry_ref(i))then
        do j=1,no_modes
          if (j/=i) then
            do k=1,integration_points
              do alpha=1,nbasis
                do beta=1,nbasis
                  ! note here is v_ji, also change integration points order
                  scf_lin_coupled(i,k) = scf_lin_coupled(i,k)               &
                                     & + eigenvectors_old(j,alpha,state(j)) &
                                     & * eigenvectors_old(j,beta,state(j))  &
                                     & * sum(   basis(j,:,alpha)            &
                                     &        * coupled_pot(j,i,:,k)        &
                                     &        * basis(j,:,beta) )       &
                                     & * dq(j)
                enddo
              enddo
            enddo
          endif
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

  ! Construct and diagonalise Hamiltonian.
  do i=1,no_modes
    ! Construct Hamiltonian.
    ! harmonic diagonal part: on diagonal we have sho eigenvalues.  
    ! add on matrix elements of anharmonic potential w.r.t. sho eigenfunctions.
    do j=1,nbasis
      do k=1,nbasis
        hamiltonian(j,k) = sum(basis(i,:,j)*basis(i,:,k)*scf_pot(i,:)) &
                       & * dq(i)
      enddo
      hamiltonian(j,j) = hamiltonian(j,j) + (j-1+0.5_dp)*omega(i) + v0(i)
    enddo
    
    ! Diagonalise Hamiltonian.
    eigenstuff = calculate_eigenstuff(hamiltonian)
    eigenvalues(i,:) = eigenstuff%evals
    eigenvectors(i,:,:) = eigenstuff%evecs
  enddo

  ! total final energy
  final_energy0=0
  do i=1,no_modes
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
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    
    final_energy = final_energy0 - sum(scf_coupled)
  endif

  max_diff=abs(final_energy-final_energy_pre)
  
  output%scf_pot = scf_pot
  output%scf_lin_coupled = scf_lin_coupled
  output%eigenvalues = eigenvalues
  output%eigenvectors = eigenvectors
  output%final_energy = final_energy
  output%final_energy0 = final_energy0
  output%max_diff = max_diff
end function
end module
