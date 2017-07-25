! Calculates the vibrational self-consistent field (VSCF) approximation
!    to the phonon Hamiltonian
! Author: B. Monserrat
! Created: 07 February 2012
! Modifications:
!   - 14/02/2012 : added MP2
!   - 29/02/2012 : parallelised coupled matrix element calculation

! Compile and run with:
! mpifort -llapack -lpthread vscf.f90 -o vscf
! mpirun -np <nnodes> vscf

! ----------------------------------------------------------------------
! Module for handling parallelism.
! See MPI specification at e.g. http://www.mpi-forum.org/docs/docs.html.
! ----------------------------------------------------------------------
!module parallel
!  implicit none
!  include 'mpif.h'
!  integer :: nnodes ! number of nodes.
!  integer :: my_node ! node (between 0 and nnodes-1).
!  logical :: am_master  ! node zero.  only node that writes to standard out.
!  integer :: ierror,mpistatus(mpi_status_size)
!contains
!
!  ! set up mpi and decide whether we are the master node.
!  subroutine initialise_parallel
!    implicit none
!    call mpi_init(ierror)
!    call mpi_comm_size(mpi_comm_world,nnodes,ierror)
!    call mpi_comm_rank(mpi_comm_world,my_node,ierror)
!    am_master=(my_node==0)
!  end subroutine
!
!  ! call mpi_finalize.
!  subroutine finalise_parallel
!    implicit none
!    call mpi_finalize(ierror)
!  end subroutine
!
!  ! abort the run (ending the program on all nodes).
!  subroutine abort_parallel
!    implicit none
!    call mpi_abort(mpi_comm_world,-1,ierror)
!    stop
!  end subroutine
!end module parallel

module vscf_module
  use string_module
  use io_module
  use constants_module, only : dp
contains

! ----------------------------------------------------------------------
! vibrational self-consistent field main loop
! ----------------------------------------------------------------------
subroutine vscf()
  use constants_module, only : pi, ev_per_hartree, kb_in_au
  !use parallel
  implicit none
  
  real(dp) :: v_indep,v_coupled
  integer :: first_mode, last_mode, no_modes, no_indep_params, &
     & no_coupled_params
  type(String) :: functional
  real(dp),allocatable :: indep_params(:,:), coupled_params(:,:,:)
  real(dp) :: mireia,crispin
  real(dp),allocatable :: indep_pot(:,:), indep_pot_bare(:,:), &
     &coupled_pot(:,:,:,:), coupled_pot_bare(:,:,:,:), scf_pot(:,:)
  real(dp),allocatable :: max_amplitude(:),harmonic_freq(:)
  real(dp),allocatable :: shift_harmonic_freq(:,:)
  real(dp),allocatable :: eigenvalues(:,:), eigenvalues_old(:,:), &
     &eigenvectors(:,:,:), basis(:,:,:), eigenvectors_old(:,:,:)
  real(dp),allocatable :: shift_eigenvalues(:,:)
  real(dp),allocatable :: omega(:),x0(:),v0(:)
  real(dp) :: omega_temp,density
  real(dp) :: amp1,amp2,amp_ratio
  real(dp),allocatable :: root_two_over_i(:),root_i_over_two(:)
  integer :: nbasis, nstates
  real(dp) :: bfp,dump1,dump2
  type(String) :: vec_name, basis_name, amp_name, file_type, pot_name, &
     & coupled_pot_name, indep_name, coupled_name
  integer :: integration_points
  real(dp) :: q,dq,q1,dq1,q2,dq2
  integer,parameter :: min_t=5,max_t=1000
  integer,allocatable :: master(:),source(:),msource(:)
  real(dp) :: max_diff
  real(dp),allocatable :: scf_indep(:),scf_lin_coupled(:,:),scf_coupled(:,:)
  real(dp),allocatable :: scf_indep2(:,:),scf_coupled2(:,:,:,:)
  real(dp),allocatable :: hamiltonian(:,:),e(:)
  integer,allocatable :: state(:)
  integer :: alpha,beta,gamma1,delta,epsilon1
  real(dp) :: int_temp1,int_temp2
  real(dp),allocatable :: work(:)
  integer :: lwork,info
  real(dp) :: final_energy,final_energy_pre,final_energy0
  real(dp) :: tol=1.d-8,tolerance=1.d-18,eigen_tol=1.d-5
  integer :: counter,total_counter
  real(dp) :: mp2
  real(dp),allocatable :: single_integral(:,:,:),coupled_integral(:,:,:,:,:,:)
  real(dp),allocatable :: partial_coupled_integral(:,:,:,:,:,:), &
     & temp_coupled_integral(:,:,:,:,:,:), &
     & temp_partial_coupled_integral(:,:,:,:,:)
  integer,allocatable :: symmetry_mode(:),symmetry_ref(:)
  integer :: inode,avg_modes_per_node,message_id
  integer :: no_unit_cells,execution_mode
  real(dp) :: harmonic_energy,anharmonic_energy,mp2_energy
  real(dp),allocatable :: vec_11(:),vec_21(:),vec_22(:),vec_23(:)
  real(dp) :: temperature, thermal_energy, fenergy, hfenergy, int_energy, &
     & hint_energy, temp_int_energy, temp_hint_energy
  real(dp),allocatable :: part_fn(:),har_part_fn(:)
  
  ! Input files.
  type(String), allocatable :: input_file(:)
  type(String), allocatable :: fit_input_file(:)
  type(String), allocatable :: symmetry_file(:)
  type(String), allocatable :: amplitude_ratio_file(:)
  type(String), allocatable :: amplitude_file(:)
  type(String), allocatable :: basis_file(:)
  type(String), allocatable :: vec_file(:)
  type(String), allocatable :: coupled_file(:)
  type(String), allocatable :: potential_file(:)
  type(String), allocatable :: coupled_potential_file(:)
  
  ! Output files.
  integer :: mireia_file
  integer :: scf_pot_file
  integer :: convergence_file
  integer :: output_file
  integer :: result_file
  
  ! Temporary variables
  integer :: i,j,k,l,ialloc,t,m,n,imax,imin
  integer :: i2,j2
  type(String), allocatable :: line(:)
  
  !call initialise_parallel

  !if(am_master)then

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
  else if(execution_mode==1)then
    call print_line('calculation with coupled modes.')
  else
    call print_line('execution mode error: invalid input.')
    call err()
  endif
  call print_line('the total number of integration points is '// &
     & integration_points)
  !endif ! am_master

  ! broadcast parameters
  !call mpi_bcast(no_modes,1,mpi_integer,0,mpi_comm_world,ierror)
  !call mpi_bcast(nbasis,1,mpi_integer,0,mpi_comm_world,ierror)
  !call mpi_bcast(integration_points,1,mpi_integer,0,mpi_comm_world,ierror)
  !call mpi_bcast(no_coupled_params,1,mpi_integer,0,mpi_comm_world,ierror)
  !call mpi_bcast(first_mode,1,mpi_integer,0,mpi_comm_world,ierror)
  !call mpi_bcast(no_unit_cells,1,mpi_integer,0,mpi_comm_world,ierror)
  !call mpi_bcast(execution_mode,1,mpi_integer,0,mpi_comm_world,ierror)

  ! allocate various arrays
  allocate( basis(no_modes,integration_points,nbasis),           &
          & max_amplitude(no_modes),                             &
          & x0(no_modes),                                        &
          & coupled_params(no_modes,no_modes,no_coupled_params), &
          & master(no_modes),                                    &
          & source(no_modes),                                    &
          & msource(no_modes),                                   &
          & vec_11(no_modes),                                    &
!          & vec_12(no_modes),                                    &
!          & vec_13(no_modes),                                    &
          & vec_21(no_modes),                                    &
          & vec_22(no_modes),                                    &
          & vec_23(no_modes),                                    &
          & part_fn(no_modes),                                   &
          & har_part_fn(no_modes),                               &
          & stat=ialloc); call err(ialloc)
  
  if(execution_mode==1)then
    allocate( coupled_pot( no_modes, no_modes,                           &
            &              integration_points, integration_points),      &
            & coupled_pot_bare( no_modes, no_modes,                      &
            &                   integration_points, integration_points), &
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

  !if(am_master)then
    allocate( indep_params(no_modes,no_indep_params),       &
            & harmonic_freq(no_modes),                      &
            & shift_harmonic_freq(no_modes,nbasis),         &
            & indep_pot(no_modes,integration_points),       &
            & indep_pot_bare(no_modes,integration_points),  &
            & scf_pot(no_modes,integration_points),         &
            & eigenvalues(no_modes,nbasis),                 &
            & eigenvectors(no_modes,nbasis,nstates),        &
            & shift_eigenvalues(no_modes,nbasis),           &
            & eigenvalues_old(no_modes,nbasis),             &
            & eigenvectors_old(no_modes,nbasis,nstates),    &
            & omega(no_modes),                              &
            & v0(no_modes),                                 &
            & hamiltonian(nbasis,nbasis),                   &
            & e(nbasis),                                    &
            & root_two_over_i(nbasis),                      &
            & root_i_over_two(nbasis),                      &
            & scf_indep(no_modes),                          &
            & scf_lin_coupled(no_modes,integration_points), &
            & scf_coupled(no_modes,no_modes),               &
            & scf_indep2(no_modes,nstates),                 &
            & scf_coupled2(no_modes,integration_points,nstates,nstates), &
            & state(no_modes),                                           &
            & single_integral(no_modes,nbasis,nbasis),                   &
! &coupled_integral(no_modes,no_modes,nbasis,nbasis,nbasis,nbasis),       &
            & symmetry_mode(no_modes),                                   &
            & symmetry_ref(no_modes),                                    &
            & stat=ialloc); call err(ialloc)

    ! read in the various input files
    indep_name = "indep_fit_parameters_"
    coupled_name = "coupled_fit_params_"
    amp_name = "max_amplitude."
    vec_name = "anharmonic_eigenvectors_mode_"
    pot_name = "indep_potential/indep_pot_"
    coupled_pot_name = "coupled_potential/coupled_pot_"
    basis_name = "functions_basis/basis_"
    file_type = ".dat"

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
    amp1 = dble(line(1))
    amp2 = dble(line(2))
    amp_ratio = amp1/amp2

    do i=first_mode,last_mode
      i2 = i-first_mode+1
      if (symmetry_mode(i2)==symmetry_ref(i2)) then
        amplitude_file = read_lines(amp_name//i//file_type)
        line = split(amplitude_file(1))
        max_amplitude(i2) = dble(line(1)) * amp_ratio
        harmonic_freq(i2) = dble(line(2))
      else
        max_amplitude(i2) = max_amplitude(symmetry_ref(i2)-(first_mode-1))
        harmonic_freq(i2) = harmonic_freq(symmetry_ref(i2)-(first_mode-1))
      endif
    enddo

    ! read in eigenvectors
    do i=first_mode,last_mode
      i2 = i-first_mode+1
      if (symmetry_mode(i2)==symmetry_ref(i2)) then
        basis_file = read_lines(basis_name//i//file_type)
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
          vec_file = read_lines(vec_name//i//file_type)
          do j=1,nstates
            do k=1,nbasis
              eigenvectors(i2,k,j) = dble(vec_file((j-1)*(nbasis+1)+k+3))
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
    !    open(1,file=indep_name//i//file_type)
    !    do j=1,no_indep_params
    !      read(1,*)indep_params(i2,j)
    !    enddo ! j
    !    close(1)
    !  else
    !    do j=1,no_indep_params
    !      indep_params(i2,j)=indep_params(symmetry_ref(i2)-(first_mode-1),j)
    !    enddo ! j
    !  endif ! symmetry_mode==symmetry_ref
    !enddo ! i

    ! Read in coupled potential parameters.
    if (execution_mode==1) then
      coupled_params = 0
      do i=first_mode,last_mode
        i2 = i-first_mode+1
        if (symmetry_mode(i2)==symmetry_ref(i2)) then
          do j=first_mode,last_mode
            j2 = j-first_mode+1
            if(j>i.or.(j<i.and.symmetry_mode(j2)/=symmetry_ref(j2)))then
              if (file_exists(coupled_name//i//'.'//j//file_type)) then
                coupled_file = read_lines(coupled_name//i//'.'//file_type)
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
    
    !write(*,*)'coupled parameters'
    !write(*,*)coupled_params(1,2,1)
    !write(*,*)coupled_params(10,18,1)
    !write(*,*)coupled_params(12,21,1)
    !write(*,*)coupled_params(21,29,1)
    
    ! Read in basis function defining potential.
    do i=first_mode,last_mode
      i2 = i-first_mode+1
      if(symmetry_mode(i2)==symmetry_ref(i2))then
        basis_file = read_lines(basis_name//i//file_type)
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
    write(*,*)'calculating basis functions...'
    do i=first_mode,last_mode
      i2 = i-first_mode+1
      bfp=(omega(i2)/pi)**0.25_dp
      do j=1,nbasis
        root_i_over_two(j)=sqrt(0.5_dp*j)
        root_two_over_i(j)=1/root_i_over_two(j)
      enddo ! j
      q=-max_amplitude(i2)
      dq=2*max_amplitude(i2)/integration_points
      do j=1,integration_points
        basis(i2,j,1)=bfp*exp(-0.5_dp*q*q*omega(i2))
        basis(i2,j,2)=sqrt(omega(i2))*q*root_two_over_i(1)*basis(i2,j,1)
        do k=3,nbasis
          basis(i2,j,k)=root_two_over_i(k-1)*(sqrt(omega(i2))*q*basis(i2,j,k-1) &
           &-root_i_over_two(k-2)*basis(i2,j,k-2))
        enddo
        q=q+dq 
      enddo
    enddo
    write(*,*)'done.'
    write(*,*)

    ! Calculate potentials.
    !open(1,file='ortuzar.dat')
    write(*,*)'reading in potentials...'
!    do i=1,no_modes
!      q = -max_amplitude(i)
!      dq = 2*max_amplitude(i)/integration_points
!      do j=1,integration_points
!        ! anharmonic potential.  nb, instead of shifting basis functions and
!        ! integration grid, we shift the potential.
!        indep_pot(i,j) = v_indep(i,no_modes,no_indep_params,indep_params,q-x0(i))-0.5_dp*omega(i)*omega(i)*q*q-v0(i)
!        !if(i==9)write(*,*)q,indep_pot(i,j)
!        q = q+dq
!      enddo ! j
!    !  write(1,*)
!    !  write(1,*)
!    enddo ! i
!    !close(1)
    do i=1,no_modes
      if(symmetry_mode(i)==symmetry_ref(i))then
        l=i+(first_mode-1)
        q = -max_amplitude(i)
        dq = 2*max_amplitude(i)/integration_points
        potential_file = read_lines(pot_name//l//file_type)
        do j=1,integration_points
          line = split(potential_file(j))
          indep_pot_bare(i,j) = dble(line(2))
          indep_pot(i,j) = indep_pot_bare(i,j) &
                       & - 0.5_dp*omega(i)*omega(i)*q*q-v0(i)
          q=q+dq
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
              if (file_exists(coupled_pot_name//l//'.'//k//file_type)) then
                coupled_potential_file = read_lines(coupled_pot_name//l//'.'//k//file_type)
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
    !  !        q1 = -max_amplitude(i)
    !  !        dq1 = 2*max_amplitude(i)/integration_points
    !  !        dq2 = 2*max_amplitude(j)/integration_points
    !  !        do k=1,integration_points
    !  !          q2=-max_amplitude(j)
    !  !          do l=1,integration_points
    !  !            coupled_pot(j,i,k,l) = v_coupled(j,i,no_modes,no_coupled_params,coupled_params,q2+x0(j),q1+x0(i))
    !  !            q2=q2+dq2
    !  !          enddo ! l
    !  !          q1=q1+dq1
    !  !        enddo ! k
    !  !      endif
    !  !    enddo 
    !      do j=i+1,no_modes
    !        q1 = -max_amplitude(i)
    !        dq1 = 2*max_amplitude(i)/integration_points
    !        dq2 = 2*max_amplitude(j)/integration_points
    !        do k=1,integration_points
    !          q2=-max_amplitude(j)
    !          do l=1,integration_points
    !            ! nb, no shift of potential because v_{ppa}(q) = sum_i v_i(q_i) + sum_ij v_ij(q_i,q_j), and the 
    !            ! shifted potential is v_{ppa}(q-q0) + v0 and we add v0 to the independent potential above. 
    !            coupled_pot(i,j,k,l) = v_coupled(i,j,no_modes,no_coupled_params,coupled_params,q1-x0(i),q2-x0(j)) 
    !            q2=q2+dq2
    !          enddo ! l
    !          q1=q1+dq1
    !        enddo ! k
    !      enddo ! j
    !  !  endif ! symmetry
    !  enddo ! i
    !endif ! execution mode
    write(*,*)'done.'
    write(*,*)

    ! state to consider
    !open(2,file='state.dat')
    !do i=1,no_modes
    !  read(2,*)state(i)
    !enddo ! i
    !close(2)
    state=1

    ! calculate matrix elements
    write(*,*)'calculating matrix elements...'

  !endif ! am_master

  ! broadcast basis, maximum amplitude and coupled potential
  !call mpi_bcast(basis,no_modes*nbasis*integration_points,mpi_double_precision,0,&
  ! &mpi_comm_world,ierror)
  !if(execution_mode==1)call mpi_bcast(coupled_pot,no_modes*no_modes*integration_points*integration_points,&
  ! &mpi_double_precision,0,mpi_comm_world,ierror)
  !call mpi_bcast(max_amplitude,no_modes,mpi_double_precision,0,mpi_comm_world,ierror)
  !call mpi_bcast(x0,no_modes,mpi_double_precision,0,mpi_comm_world,ierror)
  !if(execution_mode==1)call mpi_bcast(coupled_params,no_modes*no_modes*no_coupled_params,mpi_double_precision,&
  ! &0,mpi_comm_world,ierror)

  !avg_modes_per_node=no_modes/nnodes
  !if(am_master)write(*,*)'     the average number of modes per node is',avg_modes_per_node

!  write(*,*)'     processor ',my_node,'starts at mode',my_node+1,'in steps of',nnodes

  if(execution_mode==1)temp_partial_coupled_integral=0
!  master=0
  source=0

  

  !do i=my_node+1,no_modes,nnodes
  do i=1,no_modes
    !if(my_node==0)then
    !  master(i)=i
    !endif ! my_node==0
    !source(i)=my_node
    if(execution_mode==1)then
      do j=i+1,no_modes
        do alpha=1,nbasis
          do beta=1,nbasis
            do l=1,integration_points
              int_temp1=0
              dq1=2*max_amplitude(i)/integration_points
              do k=1,integration_points
                int_temp1=int_temp1+basis(i,k,alpha)*basis(i,k,beta)*coupled_pot(i,j,k,l)*dq1
              enddo ! k
              temp_partial_coupled_integral(alpha,beta,i,j,l)=int_temp1
            enddo ! l
          enddo ! beta
        enddo ! alpha
      enddo ! j
    endif ! execution mode
  enddo ! i

  if(execution_mode==1)partial_coupled_integral=0
  !master=0
  !source=0

  !do i=my_node+1,no_modes,nnodes
  do i=1,no_modes
    !if(my_node==0)then
    !  master(i)=i
    !endif ! my_node==0
    !source(i)=my_node
    if(execution_mode==1)then
      do j=i+1,no_modes
        do alpha=1,nbasis
          do beta=1,nbasis
            do delta=1,nbasis
              do epsilon1=1,nbasis
                int_temp1=0
                dq2=2*max_amplitude(j)/integration_points
                do l=1,integration_points
                  int_temp1=int_temp1+basis(j,l,delta)*basis(j,l,epsilon1)*dq2*&
                   &temp_partial_coupled_integral(alpha,beta,i,j,l)
                enddo ! l
                partial_coupled_integral(alpha,beta,delta,epsilon1,i,j)=int_temp1
              enddo ! epsilon1
            enddo ! delta
          enddo ! beta
        enddo ! alpha
      enddo ! j
    endif ! execution mode
  enddo ! i

  !call mpi_reduce(source,msource,no_modes,mpi_integer,mpi_sum,0,mpi_comm_world,ierror)

  !call mpi_barrier(mpi_comm_world,ierror)


!  if(execution_mode==1)then
!    if(am_master)write(*,*)'     starting transmission of data...' 
!    do i=1,no_modes
!      do j=i+1,no_modes
!        call mpi_reduce(partial_coupled_integral(1,1,1,1,i,j),coupled_integral(1,1,1,1,i,j),&
!              &nbasis*nbasis*nbasis*nbasis,mpi_double_precision,&
!              &mpi_sum,0,mpi_comm_world,ierror)
!      enddo ! j
!    enddo ! i
!    if(am_master)write(*,*)'     finished transmission of data.' 
!  endif ! execution_mode

!  if(am_master)then
!    !write(*,*) 
!    write(*,*)'     starting transmission of data...' 
!    coupled_integral=partial_coupled_integral
!    do i=1,no_modes
!      if(master(i)/=i)then
!        !!call mpi_recv(temp_coupled_integral(1,1,1,1,1,i),&
!        !! &no_modes*nbasis*nbasis*nbasis*nbasis,mpi_double_precision,&
!        !! &msource(i),i,mpi_comm_world,mpistatus,ierror)
!        !!coupled_integral=coupled_integral+temp_coupled_integral
!        !!temp_coupled_integral=0
!        !! if large number of modes (>100?) uncomment following and comment mpi_recv above
!        do j=i+1,no_modes
!        !  call mpi_recv(temp_coupled_integral(1,1,1,1,i,j),&
!        !   &nbasis*nbasis*nbasis*nbasis,mpi_double_precision,&
!        !   &msource(i),i,mpi_comm_world,mpistatus,ierror)
!        !  coupled_integral=coupled_integral+temp_coupled_integral
!        !  temp_coupled_integral=0
!        enddo ! j
!      endif ! not master
!    enddo ! i
!    write(*,*)'     finished transmission of data.' 
!  else
!    do i=my_node+1,no_modes,nnodes
!      !!call mpi_ssend(partial_coupled_integral(1,1,1,1,1,i),&
!      !! &no_modes*nbasis*nbasis*nbasis*nbasis,mpi_double_precision,&
!      !! &0,i,mpi_comm_world,ierror)
!      !! if large number of modes (>100?) uncomment following and comment mpi_ssend above
!      do j=i+1,no_modes
!      !  call mpi_ssend(partial_coupled_integral(1,1,1,1,i,j),&
!      !   &nbasis*nbasis*nbasis*nbasis,mpi_double_precision,&
!      !   &0,i,mpi_comm_world,ierror)
!      enddo ! j
!    enddo ! i
!  endif ! am_master

  !call mpi_barrier(mpi_comm_world,ierror)

  !if(am_master)then

    final_energy=0
    do i=1,no_modes
      final_energy=final_energy+(0.5_dp+(state(i)-1))*harmonic_freq(i)
    enddo ! i
    harmonic_energy=final_energy
    
    ! Open output files.
    mireia_file = open_write_file('mireia.dat')
    call print_line(mireia_file, '0 '//final_energy)
    
    scf_pot_file = open_write_file('scf_pot.dat')
    
    convergence_file = open_write_file('convergence.dat')
    
    output_file = open_append_file('caesar.output')

    ! self-consistent field loop 
    do t=1,max_t
      ! linear coupled correction
      if(execution_mode==1)then
        scf_lin_coupled=0
        do i=1,no_modes
          if(symmetry_mode(i)==symmetry_ref(i))then
            do j=1,i-1
              if(i==1.and.t==1)q1=-max_amplitude(i);dq1=2*max_amplitude(i)/integration_points
              do k=1,integration_points
                do alpha=1,nbasis
                  do beta=1,nbasis
                    int_temp1=0
                    dq=2*max_amplitude(j)/integration_points
                    do gamma1=1,integration_points
                      ! note here is v_ji, also change integration points order
                      int_temp1=int_temp1+basis(j,gamma1,alpha)*coupled_pot(j,i,gamma1,k)*basis(j,gamma1,beta)*dq
                    enddo ! gamma
                    scf_lin_coupled(i,k)=scf_lin_coupled(i,k)+eigenvectors(j,alpha,state(j))*eigenvectors(j,beta,state(j))*int_temp1
                 enddo ! beta
                enddo ! alpha
                if (i==1.and.t==1) then
                  call print_line(mireia_file, q1//' '//scf_lin_coupled(i,k))
                  q1=q1+dq1
                endif
              enddo ! k
            enddo ! j
            do j=i+1,no_modes
              if(i==1.and.t==1)q1=-max_amplitude(i);dq1=2*max_amplitude(i)/integration_points
              do k=1,integration_points
                do alpha=1,nbasis
                  do beta=1,nbasis
                    int_temp1=0
                    dq=2*max_amplitude(j)/integration_points
                    do gamma1=1,integration_points
                      ! note here is v_ij, also change integration points order
                      int_temp1=int_temp1+basis(j,gamma1,alpha)*coupled_pot(i,j,k,gamma1)*basis(j,gamma1,beta)*dq
                    enddo ! gamma
                    scf_lin_coupled(i,k)=scf_lin_coupled(i,k)+eigenvectors(j,alpha,state(j))*eigenvectors(j,beta,state(j))*int_temp1
                  enddo ! beta
                enddo ! alpha
                if (i==1.and.t==1) then
                  call print_line(mireia_file, q1//' '//scf_lin_coupled(i,k))
                  q1=q1+dq1
                endif
              enddo ! k
              if (i==1.and.t==1) then
                call print_line(mireia_file, '')
                call print_line(mireia_file, '')
              endif
            enddo ! j
          else
            do k=1,integration_points
              scf_lin_coupled(i,k)=scf_lin_coupled(symmetry_ref(i),k) 
            enddo ! k
          endif ! symmetry
        enddo ! i
      endif ! execution mode

      ! effective self-consistent potential
      scf_pot=0
      do i=1,no_modes
        if(i==1)q=-max_amplitude(i);dq=2*max_amplitude(i)/integration_points
        do k=1,integration_points
          if(execution_mode==1)scf_pot(i,k)=scf_pot(i,k)+indep_pot(i,k)+scf_lin_coupled(i,k)
          if(execution_mode==0)scf_pot(i,k)=scf_pot(i,k)+indep_pot(i,k)
          if (i==1) then
            call print_line(scf_pot_file, q              //' '// &
                                        & scf_pot(i,k)   //' '// &
                                        & indep_pot(i,k) //' '// &
                                        & scf_lin_coupled(i,k))
          endif
          q=q+dq
        enddo ! k
        if (i==1) then
          call print_line(scf_pot_file, '')
          call print_line(scf_pot_file, '')
        endif
      enddo ! i 

      ! construct and diagonalise hamiltonian matrix
      !write(*,*)' - constructing and diagonalising hamiltonian matrix...'
      do i=1,no_modes
        ! construct hamiltonian matrix      
        hamiltonian=0
        ! harmonic diagonal part: on diagonal we have sho eigenvalues.  
        ! add on matrix elements of anharmonic potential w.r.t. sho eigenfunctions.
        do j=1,nbasis
          hamiltonian(j,j)=hamiltonian(j,j)+(real(j-1,dp)+0.5_dp)*omega(i)+v0(i)
        enddo ! j
        counter=0
        do j=1,nbasis
          do k=1,nbasis
            q=-max_amplitude(i);dq=2*max_amplitude(i)/integration_points
            int_temp1=0
            do l=1,integration_points
              int_temp1=int_temp1+basis(i,l,j)*basis(i,l,k)&
               &*scf_pot(i,l)*dq
               !&*(scf_pot(i,l)-0.5_dp*omega(i)*omega(i)*q*q)*dq
              q=q+dq
            enddo !l
            hamiltonian(j,k)=hamiltonian(j,k)+int_temp1
          enddo ! k
        enddo ! j
        ! diagonalise hamiltonian matrix
        lwork=3*nbasis-1
        allocate(work(lwork), stat=ialloc); call err(ialloc)
        call dsyev('v','u',nbasis,hamiltonian(1,1),nbasis,e(1),work(1),lwork,info)
        if(info/=0) then
          call print_line('dsyev failed. info='//info//'.')
          call err()
        endif
        deallocate(work)
        ! write out the eigenvalues and eigenvectors obtained by diagonalisation.
        do l=1,nbasis
          !eigenvalues_old(i,l)=eigenvalues(i,l)
          !eigenvalues(i,l)=0.3_dp*e(l)+0.7_dp*eigenvalues_old(i,l)
          eigenvalues(i,l)=e(l)
        enddo ! l
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
                   scf_coupled(i,j)=scf_coupled(i,j)+eigenvectors_old(i,alpha,state(i))*eigenvectors_old(i,beta,state(i))*&
                    &eigenvectors_old(j,delta,state(j))*eigenvectors_old(j,epsilon1,state(j))*&
                    &coupled_integral(alpha,beta,delta,epsilon1,i,j)
                  enddo ! epsilon1
                enddo ! delta
              enddo ! beta
            enddo ! alpha
          enddo ! j
        enddo ! i
      
        do i=1,no_modes
          do j=i+1,no_modes
            final_energy=final_energy-scf_coupled(i,j)
          enddo ! j
        enddo ! i
        !write(*,*)final_energy
      endif ! execution mode

      max_diff=abs(final_energy-final_energy_pre)
      call print_line(convergence_file, t//' '//final_energy//' '//max_diff)
      
      call print_line(output_file, t//' '//final_energy//' '//max_diff)


      if((max_diff<tol.and.t>=min_t).or.execution_mode==0)then

        open(3,file='anharmonic_eigenvalues.dat')
        open(21,file='indep_mode.dat')
        write(21,*)'independent modes frequencies'
        mireia=0
        do i=1,no_modes
          mireia=mireia+harmonic_freq(i)/2
          crispin=crispin+eigenvalues(i,1)
          write(21,*)i+3,harmonic_freq(i)/2,eigenvalues(i,1)!,eigenvalues(i,1)-harmonic_freq(i)/2,mireia,crispin
          write(3,*)'eigenvalues for mode',i
          do j=1,nbasis
            write(3,*)eigenvalues(i,j)
          enddo ! j
          write(3,*)
          write(3,*)
        enddo ! i
        close(3)
        write(21,*)
        close(21)

        do i=1,no_modes
          if(symmetry_mode(i)==symmetry_ref(i))then
            l=i+3
            vec_file = open_write_file(vec_name//l//file_type)
            call print_line(vec_file, vec_11(i) //' '// &
                                    & nstates   //' '// &
                                    & nbasis)
            call print_line(vec_file, vec_21(i) //' '// &
                                    & vec_22(i) //' '// &
                                    & vec_23(i))
            call print_line(vec_file, '')
            do j=1,nstates
              do k=1,nbasis
                call print_line(vec_file, eigenvectors(i,k,j))
              enddo ! k
              call print_line(vec_file, '')
            enddo ! j
            call print_line(vec_file, '')
            close(vec_file)
          endif ! symmetry
        enddo ! i

        open(2,file='vscf_results.dat')
        write(*,*)
        write(*,*)'converged in',t,'iterations.'
        write(*,*)'final energy (a.u., ev)',final_energy,final_energy*ev_per_hartree
        write(*,*)'final energy puc (a.u., ev)',final_energy/no_unit_cells,final_energy*ev_per_hartree/no_unit_cells
        
        call print_line(output_file, '')
        call print_line(output_file, 'converged in '//t//' iterations.')
        call print_line(output_file, 'final energy (a.u., ev) '// &
           & final_energy//' '//final_energy*ev_per_hartree)
        call print_line(output_file, 'final energy puc (a.u., ev) '// &
           & final_energy/no_unit_cells                       //' '// &
           & final_energy*ev_per_hartree/no_unit_cells)
        call print_line(output_file, '')
        
        write(2,*)'converged in',t,'iterations.'
        write(2,*)'final energy (a.u., ev, cm-1)',final_energy,final_energy*ev_per_hartree,&
         &final_energy*219472.62
        write(2,*)'energy difference is',max_diff
        close(2)
        write(*,*)
        anharmonic_energy=final_energy
 
        if(execution_mode==1)then

          ! calculate second order perturbation
          write(*,*)'calculating second order perturbation...'

          ! independent terms correction
          write(*,*)' - calculating independent terms...'
          scf_indep2=0
          do i=1,no_modes 
            do j=1,state(i)-1
              do alpha=1,nbasis
                do beta=1,nbasis
                  int_temp1=0
                  dq=2*max_amplitude(i)/integration_points
                  do gamma1=1,integration_points
                    int_temp1=int_temp1+basis(i,gamma1,alpha)*(indep_pot(i,gamma1)-scf_pot(i,gamma1))*&
                     &basis(i,gamma1,beta)*dq
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
                  dq=2*max_amplitude(i)/integration_points
                  do gamma1=1,integration_points
                    int_temp1=int_temp1+basis(i,gamma1,alpha)*(indep_pot(i,gamma1)-scf_pot(i,gamma1))*&
                     &basis(i,gamma1,beta)*dq
                  enddo ! gamma
                  scf_indep2(i,j)=scf_indep2(i,j)+eigenvectors(i,alpha,state(i))*eigenvectors(i,beta,j)*int_temp1
                enddo ! beta
              enddo ! alpha
            enddo ! j
          enddo ! i
          ! coupled terms correction
          write(*,*)' - calculating coupled terms...'
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
          write(*,*)' - calculating perturbation correction...'
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

          write(*,*)
          write(*,*)'     mp2 correction (a.u., ev, cm-1)', mp2,mp2*ev_per_hartree,mp2*219472.62
          write(*,*)'     mp2 energy (a.u., ev, cm-1)' ,final_energy+mp2,(final_energy+mp2)*ev_per_hartree,&
           &(final_energy+mp2)*219472.62
          write(*,*)'     mp2 energy puc (a.u., ev, cm-1)' ,(final_energy+mp2)/no_unit_cells,&
           &(final_energy+mp2)*ev_per_hartree/no_unit_cells,&
           &(final_energy+mp2)*219472.62/no_unit_cells
          
          result_file = open_append_file('vscf_results.dat')
          call print_line(result_file, '')
          call print_line(result_file, 'mp2 energy (a.u., ev, cm-1): '// &
             & final_energy+mp2                                  //' '// &
             & (final_energy+mp2)*ev_per_hartree                 //' '// &
             & (final_energy+mp2)*219472.62)
          
          mp2_energy=final_energy+mp2

        endif ! execution mode

        write(*,*)
        write(*,*)'-------summary of results-------'
        write(*,*)'harmonic energy (ev):',harmonic_energy*ev_per_hartree/no_unit_cells
        write(*,*)'anharmonic energy (ev):',anharmonic_energy*ev_per_hartree/no_unit_cells
        call print_line(output_file, '')
        call print_line(output_file, '-------summary of results-------')
        call print_line(output_file, 'harmonic energy (ev): '// &
           & harmonic_energy*ev_per_hartree/no_unit_cells)
        call print_line(output_file, 'anharmonic energy (ev): '// &
           & anharmonic_energy*ev_per_hartree/no_unit_cells)
        if(execution_mode==1)write(*,*)'mp2 energy (ev):',mp2_energy*ev_per_hartree/no_unit_cells


        if(temperature>tol)then
         
          shift_eigenvalues=0
          shift_harmonic_freq=0
          do i=1,no_modes
            if(eigenvalues(i,1)<0)then
              do j=1,nstates
                shift_eigenvalues(i,j)=eigenvalues(i,j)-2*eigenvalues(i,1)
              enddo ! j
            !else if(eigenvalues(i,1)<eigen_tol)then
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
          open(1,file='indep_mode_internal.dat')
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
              write(1,*)i+3,temp_hint_energy,temp_int_energy
            enddo ! i
          endif ! temp > tolerance
          close(1)

          ! calculate free energy
          fenergy=0
          hfenergy=0
          open(1,file='indep_mode_temperature.dat')
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
              write(1,*)-(thermal_energy*log(part_fn(i))-thermal_energy*log(har_part_fn(i))),-thermal_energy*log(part_fn(i)),&
               &-thermal_energy*log(har_part_fn(i))
            enddo ! i
          else
            do i=1,no_modes
              fenergy=fenergy+eigenvalues(i,1)
              hfenergy=hfenergy+0.5_dp*harmonic_freq(i)
            enddo ! i
          endif
          close(1)

          ! calculate vibrational density
          open(1,file='density.dat')
          q1=-max_amplitude(1)
          dq1=2*max_amplitude(1)/integration_points
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
            write(1,*)q1,density
            q1=q1+dq1
          enddo ! i
           close(1)

          write(*,*)
          write(*,*)'-------summary of free energy results-------'
          write(*,*)'temperature (k):',temperature
          write(*,*)'harmonic free energy (ev):',hfenergy/no_unit_cells*ev_per_hartree
          write(*,*)'anharmonic free energy (ev):',fenergy/no_unit_cells*ev_per_hartree
          call print_line(output_file, '')
          call print_line(output_file, '-------summary of free energy results-------')
          call print_line(output_file, 'temperature (k): '//temperature)
          call print_line(output_file, 'harmonic free energy (ev): '// &
             & hfenergy/no_unit_cells*ev_per_hartree)
          call print_line(output_file, 'anharmonic free energy (ev): '// &
             & fenergy/no_unit_cells*ev_per_hartree)
          write(*,*)
          write(*,*)'-------summary of int energy results-------'
          write(*,*)'temperature (k):',temperature
          write(*,*)'harmonic internal energy (ev):',hint_energy/no_unit_cells*ev_per_hartree
          write(*,*)'anharmonic internal energy (ev):',int_energy/no_unit_cells*ev_per_hartree
          call print_line(output_file, '')
          call print_line(output_file, &
             & '-------summary of int energy results-------')
          call print_line(output_file, 'temperature (k): '//temperature)
          call print_line(output_file, 'harmonic internal energy (ev): '// &
             & hint_energy/no_unit_cells*ev_per_hartree)
          call print_line(output_file, 'anharmonic internal energy (ev): '// &
             & int_energy/no_unit_cells*ev_per_hartree)
        endif ! temperature > 0
        exit
      endif 
      
      if (t==max_t) then
        !call finalise_parallel
        call print_line('the scf calculation did not converge')
        call err()
      endif
    enddo
    close(mireia_file)
    close(scf_pot_file)
    close(convergence_file)
    close(output_file)
  !endif ! am_master
  !call finalise_parallel
end subroutine

! ======================================================================
! functional forms for the one- and two-body potentials
! ======================================================================

! ----------------------------------------------------------------------
! single mode potential
! ----------------------------------------------------------------------
function v_indep(i,no_modes,no_indep_params,indep_params,q) result(output)
  implicit none
  
  integer,intent(in) :: i,no_modes,no_indep_params
  real(dp),intent(in) :: indep_params(no_modes,no_indep_params),q
  real(dp) :: output
  integer :: order
  character(100) :: functional

  open(1,file='fit_input.dat')
  read(1,*)functional,order
  close(1)

  if(functional=='polynomial')then
    if(order==6)then
      output = indep_params(i,1)*q**2 + indep_params(i,2)*q**3 + &
       &indep_params(i,3)*q**4 + indep_params(i,4)*q**5 + &
       &indep_params(i,5)*q**6
    else if(order==8)then
      output = indep_params(i,1)*q**2 + indep_params(i,2)*q**3 + &
       &indep_params(i,3)*q**4 + indep_params(i,4)*q**5 + &
       &indep_params(i,5)*q**6 + indep_params(i,6)*q**7 + &
       &indep_params(i,7)*q**8
    else if(order==10)then
      output = indep_params(i,1)*q**2 + indep_params(i,2)*q**3 + &
       &indep_params(i,3)*q**4 + indep_params(i,4)*q**5 + &
       &indep_params(i,5)*q**6 + indep_params(i,6)*q**7 + &
       &indep_params(i,7)*q**8 + indep_params(i,8)*q**9 + &
       &indep_params(i,9)*q**10
    else if(order==12)then
      output = indep_params(i,1)*q**2 + indep_params(i,2)*q**3 + &
       &indep_params(i,3)*q**4 + indep_params(i,4)*q**5 + &
       &indep_params(i,5)*q**6 + indep_params(i,6)*q**7 + &
       &indep_params(i,7)*q**8 + indep_params(i,8)*q**9 + &
       &indep_params(i,9)*q**10 + indep_params(i,10)*q**11 + &
       &indep_params(i,11)*q**12
    endif
  else if(functional=='sine')then
    if(order==6)then
      output = indep_params(i,1)*q**2 + indep_params(i,2)*q**3 + &
       &indep_params(i,3)*q**4 + indep_params(i,4)*q**5 + &
       &indep_params(i,5)*q**6 + &
       &indep_params(i,6)*sin(indep_params(i,7)*q)**2
    else if(order==8)then
      output = indep_params(i,1)*q**2 + indep_params(i,2)*q**3 + &
       &indep_params(i,3)*q**4 + indep_params(i,4)*q**5 + &
       &indep_params(i,5)*q**6 + indep_params(i,6)*q**7 + &
       &indep_params(i,7)*q**8 + &
       &indep_params(i,8)*sin(indep_params(i,9)*q)**2
    else if(order==10)then
      output = indep_params(i,1)*q**2 + indep_params(i,2)*q**3 + &
       &indep_params(i,3)*q**4 + indep_params(i,4)*q**5 + &
       &indep_params(i,5)*q**6 + indep_params(i,6)*q**7 + &
       &indep_params(i,7)*q**8 + indep_params(i,8)*q**9 + &
       &indep_params(i,9)*q**10 + &
       &indep_params(i,10)*sin(indep_params(i,11)*q)**2
    endif
  else
    call print_line('error: the functional form of the independent potential &
       &is incorrect!')
    call err()
  endif
end function

! ----------------------------------------------------------------------
! two-mode coupling potential
! ----------------------------------------------------------------------
function v_coupled(i,j,no_modes,no_coupled_params,coupled_params,q1,q2) &
   & result(output)
  implicit none
  
  integer,intent(in) :: i,j,no_modes,no_coupled_params
  real(dp),intent(in) :: q1,q2,coupled_params(no_modes,no_modes,no_coupled_params)   
  real(dp) :: output

  if(no_coupled_params==1)then
    output = coupled_params(i,j,1)*q1**2*q2**2
  else if(no_coupled_params==2)then
    output = coupled_params(i,j,1)*q1*q2 + coupled_params(i,j,2)*q1**2*q2**2
  else if(no_coupled_params==3)then
   output = coupled_params(i,j,1)*q1**2*q2**2 + &
    &(coupled_params(i,j,2)*q1**2*q2+coupled_params(i,j,3)*q1*q2**2)
  else if(no_coupled_params==4)then
    output = coupled_params(i,j,1)*q1*q2 + coupled_params(i,j,2)*q1**2*q2**2 +&
     &coupled_params(i,j,3)*q1**2*q2 + coupled_params(i,j,4)*q1*q2**2
  else if(no_coupled_params==5)then
    output = coupled_params(i,j,1)*sin(coupled_params(i,j,2)*q1)**2*&
               &sin(coupled_params(i,j,3)*q2)**2 +&
               &coupled_params(i,j,4)*q1**3*q2 +&
               &coupled_params(i,j,5)*q1*q2**3 
  else if(no_coupled_params==6)then
    output = coupled_params(i,j,1)*q1*q2 + coupled_params(i,j,2)*q1**2*q2 +&
     &coupled_params(i,j,3)*q1*q2**2 + coupled_params(i,j,4)*q1**3*q2 +&
     &coupled_params(i,j,5)*q1**2*q2**2 + coupled_params(i,j,6)*q1*q2**3
  else if(no_coupled_params==15)then
    output=coupled_params(i,j,1)*q1*q2 + coupled_params(i,j,2)*q1**2*q2 +&
     & coupled_params(i,j,3)*q1*q2**2 + coupled_params(i,j,4)*q1**3*q2 +&
     & coupled_params(i,j,5)*q1**2*q2**2 + coupled_params(i,j,6)*q1*q2**3 +&
     &coupled_params(i,j,7)*q1**4*q2 + coupled_params(i,j,8)*q1**3*q2**2 +&
     & coupled_params(i,j,9)*q1**2*q2**3 + coupled_params(i,j,10)*q1*q2**4 +&
     &coupled_params(i,j,11)*q1**5*q2+coupled_params(i,j,12)*q1**4*q2**2+&
     &coupled_params(i,j,13)*q1**3*q2**3+coupled_params(i,j,14)*q1**2*q2**4+&
     &coupled_params(i,j,15)*q1*q2**5
  else if(no_coupled_params==45)then
   output = coupled_params(i,j,1)*q1*q2 + &
    &(coupled_params(i,j,2)*q1**2*q2+coupled_params(i,j,3)*q1*q2**2) + &
    &(coupled_params(i,j,4)*q1**3*q2+coupled_params(i,j,5)*q1**2*q2**2+ &
    &coupled_params(i,j,6)*q1*q2**3) + &
    &(coupled_params(i,j,7)*q1**4*q2+coupled_params(i,j,8)*q1**3*q2**2+coupled_params(i,j,9)*q1**2*q2**3 + &
    &coupled_params(i,j,10)*q1*q2**4) + &
    &(coupled_params(i,j,11)*q1**5*q2+coupled_params(i,j,12)*q1**4*q2**2+coupled_params(i,j,13)*q1**3*q2**3 + &
    &coupled_params(i,j,14)*q1**2*q2**4+coupled_params(i,j,15)*q1*q2**5) +&
    &(coupled_params(i,j,16)*q1**6*q2+coupled_params(i,j,17)*q1**5*q2**2+coupled_params(i,j,18)*q1**4*q2**3 + &
    &coupled_params(i,j,19)*q1**3*q2**4+coupled_params(i,j,20)*q1**2*q2**5+&
    &coupled_params(i,j,21)*q1*q2**6) + &
    &(coupled_params(i,j,22)*q1**7*q2+coupled_params(i,j,23)*q1**6*q2**2+coupled_params(i,j,24)*q1**5*q2**3 + &
    &coupled_params(i,j,25)*q1**4*q2**4+coupled_params(i,j,26)*q1**3*q2**5+coupled_params(i,j,27)*q1**2*q2**6 + &
    &coupled_params(i,j,28)*q1*q2**7) + &
    &(coupled_params(i,j,29)*q1**8*q2+coupled_params(i,j,30)*q1**7*q2**2+coupled_params(i,j,31)*q1**6*q2**3 + &
    &coupled_params(i,j,32)*q1**5*q2**4+coupled_params(i,j,33)*q1**4*q2**5+coupled_params(i,j,34)*q1**3*q2**6 + &
    &coupled_params(i,j,35)*q1**2*q2**7+coupled_params(i,j,36)*q1*q2**8) +&
    &(coupled_params(i,j,37)*q1**9*q2+coupled_params(i,j,38)*q1**8*q2**2+coupled_params(i,j,39)*q1**7*q2**3 + &
    &coupled_params(i,j,40)*q1**6*q2**4+coupled_params(i,j,41)*q1**5*q2**5+coupled_params(i,j,42)*q1**4*q2**6 + &
    &coupled_params(i,j,43)*q1**3*q2**7+coupled_params(i,j,44)*q1**2*q2**8+& 
    &coupled_params(i,j,45)*q1*q2**9)
  else
    call print_line('error: the functional form of the coupled potential is &
       &incorrect!')
    call err()
  endif
end function
end module
