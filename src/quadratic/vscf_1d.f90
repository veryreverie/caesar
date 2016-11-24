module vscf_1d_module
  use constants, only : dp
  implicit none
  
  ! holds output data for vscf_1d
  type VscfReturn
    real(dp), allocatable :: anh_pot(:,:)
    real(dp), allocatable :: hamiltonian(:,:)
    real(dp), allocatable :: eigenvals(:,:)
    real(dp), allocatable :: eigenvecs(:,:)
  end type
  
contains

subroutine vscf_1d()
  use constants,      only : dp, pi, eV, thermal
  use linear_algebra, only : dsyev
  implicit none
 
  ! Hard coded numerical parameters
  ! If this value is changed, it also needs to be changed in 'calculate_anharmonic.f90'
  integer :: Nbasis=20 

  ! Input variables
  integer :: integration_points
  real(dp) :: frequency,max_amplitude
 
  ! Working variables
  integer :: i,j,k
  integer :: lwork,info
  real(dp) :: bfp,q,dq,dump_real,mat_element,basis_frequency
  real(dp),allocatable :: root_two_over_i(:),root_i_over_two(:)
  real(dp),allocatable :: basis(:,:),potential(:),potential_bare(:)
  real(dp),allocatable :: Hamiltonian(:,:),E(:)
  real(dp),allocatable :: work(:)

  ! Allocate various arrays
  allocate(root_two_over_i(Nbasis),root_i_over_two(Nbasis))
  allocate(Hamiltonian(Nbasis,Nbasis),E(Nbasis))

  ! Read in frequency (note that it is in eV)
  open(1,file='frequency.dat')
  read(1,*)frequency
  close(1)
  ! Convert frequency to a.u.
  frequency=frequency/eV 

  ! Read in maximum amplitude
  open(1,file='max_amplitude.dat')
  read(1,*)max_amplitude
  close(1)
  ! The maximum amplitude is read in as negative, turn positive
  max_amplitude=-max_amplitude

  ! Read in integration points
  open(1,file='integration_points.dat')
  read(1,*)integration_points
  close(1)
  allocate(basis(integration_points,Nbasis))
  allocate(potential(integration_points),potential_bare(integration_points))

  ! Calculate basis functions
  basis_frequency=abs(frequency)
  bfp=(basis_frequency/pi)**0.25
  do i=1,Nbasis
    root_i_over_two(i)=SQRT(0.5d0*REAL(i,dp))
    root_two_over_i(i)=1.d0/root_i_over_two(i)
  enddo ! i
  q=-max_amplitude
  dq=2.d0*max_amplitude/integration_points
  do i=1,integration_points
    basis(i,1)=bfp*exp(-0.5d0*q*q*basis_frequency)
    basis(i,2)=sqrt(basis_frequency)*q*root_two_over_i(1)*basis(i,1)
    do j=3,Nbasis
      basis(i,j)=root_two_over_i(j-1)*(sqrt(basis_frequency)*q*basis(i,j-1) &
       &-root_i_over_two(j-2)*basis(i,j-2))
    enddo ! j
    q=q+dq
  enddo ! i

  ! Read in potential
  open(1,file='interp_energy.dat')
  open(2,file='anh_pot.dat')
  q=-max_amplitude
  dq=2.d0*max_amplitude/integration_points
  do i=1,integration_points
    read(1,*)dump_real,potential_bare(i)
    potential_bare(i)=potential_bare(i)/eV
    potential(i)=potential_bare(i)-0.5d0*frequency*frequency*q*q
    write(2,*)q,potential(i),potential_bare(i),0.5d0*frequency*frequency*q*q
    q=q+dq
  enddo ! i
  close(1)
  close(2)

  

  ! Construct and diagonalise Hamiltonian matrix
  Hamiltonian=0.d0
  ! Harmonic diagonal part: On diagonal we have SHO eigenvalues.  
  do i=1,Nbasis
    Hamiltonian(i,i)=Hamiltonian(i,i)+(real(i-1,dp)+0.5d0)*basis_frequency
  enddo ! i
  ! Add on matrix elements of anharmonic potential w.r.t. SHO eigenfunctions.
  do i=1,Nbasis
    do j=1,Nbasis
      dq=2.d0*max_amplitude/integration_points
      mat_element=0.d0
      do k=1,integration_points
        mat_element=mat_element+basis(k,i)*basis(k,j)*potential(k)*dq
        !write(*,*)mat_element,basis(k,i),basis(k,j),potential(k),dq
      enddo ! k
      Hamiltonian(i,j)=Hamiltonian(i,j)+mat_element
    enddo ! j
  enddo ! i
 open(1,file='hamiltonian.dat')
 write(1,*)Hamiltonian
 close(1)
 ! Diagonalise Hamiltonian matrix
 lwork=3*Nbasis-1
 allocate(work(lwork))
 call dsyev('V','U',Nbasis,Hamiltonian(1,1),Nbasis,E(1),work(1),lwork,info)
 deallocate(work)

 ! Output anharmonic eigenvalues and eigenvectors
 open(1,file='eigenvals.dat')
 open(2,file='eigenvecs.dat')
 do i=1,Nbasis
   write(1,*)i,frequency*(real(i-1,dp)+0.5d0)*eV,E(i)*eV
   do j=1,Nbasis
     write(2,*)i,j,Hamiltonian(j,i)
   enddo ! j
   write(2,*)
 enddo ! i
 close(1)
 close(2)
end subroutine

! as above, but takes arguments rather than reading files
pure function vscf_1d_2(frequency_ev,max_amplitude_nve,spline_potential, Nbasis) &
  &result(output)
  use constants,      only : dp, pi, eV, thermal
  use linear_algebra, only : dsyev
  implicit none
  
  real(dp), intent(in) :: frequency_ev        ! frequency in eV
  real(dp), intent(in) :: max_amplitude_nve   ! max amplitude (negative)
  real(dp), intent(in) :: spline_potential(:) ! 
  integer,  intent(in) :: Nbasis
  type(VscfReturn)     :: output
  
  real(dp) :: frequency          ! frequency in a.u.
  real(dp) :: max_amplitude      ! -1*max_amplitude_nve
  integer  :: integration_points ! size(spline_potential)
  
  ! Working variables
  integer :: i,j,k
  integer :: lwork,info
  real(dp) :: bfp,q,dq,dump_real,mat_element,basis_frequency
  real(dp),allocatable :: root_two_over_i(:),root_i_over_two(:)
  real(dp),allocatable :: basis(:,:),potential(:),potential_bare(:)
  real(dp),allocatable :: Hamiltonian(:,:),E(:)
  real(dp),allocatable :: work(:)
  
  ! Convert frequency to a.u.
  frequency=frequency_ev/eV 

  ! The maximum amplitude is read in as negative, turn positive
  max_amplitude=-max_amplitude_nve
  
  ! get integration_points
  integration_points = size(spline_potential)
  
  ! Allocate output
  allocate(output%anh_pot(4,integration_points))
  allocate(output%hamiltonian(Nbasis,Nbasis))
  allocate(output%eigenvals(3,Nbasis))
  allocate(output%eigenvecs(3,Nbasis*Nbasis))

  ! Allocate various arrays
  allocate(root_two_over_i(Nbasis),root_i_over_two(Nbasis))
  allocate(Hamiltonian(Nbasis,Nbasis),E(Nbasis))
  allocate(basis(integration_points,Nbasis))
  allocate(potential(integration_points),potential_bare(integration_points))

  ! Calculate basis functions
  basis_frequency=abs(frequency)
  bfp=(basis_frequency/pi)**0.25
  do i=1,Nbasis
    root_i_over_two(i)=SQRT(0.5d0*REAL(i,dp))
    root_two_over_i(i)=1.d0/root_i_over_two(i)
  enddo ! i
  q=-max_amplitude
  dq=2.d0*max_amplitude/integration_points
  do i=1,integration_points
    basis(i,1)=bfp*exp(-0.5d0*q*q*basis_frequency)
    basis(i,2)=sqrt(basis_frequency)*q*root_two_over_i(1)*basis(i,1)
    do j=3,Nbasis
      basis(i,j)=root_two_over_i(j-1)*(sqrt(basis_frequency)*q*basis(i,j-1) &
       &-root_i_over_two(j-2)*basis(i,j-2))
    enddo ! j
    q=q+dq
  enddo ! i

  ! write anh_pot
  q=-max_amplitude
  dq=2.d0*max_amplitude/integration_points
  do i=1,integration_points
    potential_bare(i) = spline_potential(i)/eV
    potential(i)=potential_bare(i)-0.5d0*frequency*frequency*q*q
    output%anh_pot(1,i) = q
    output%anh_pot(2,i) = potential(i)
    output%anh_pot(3,i) = potential_bare(i)
    output%anh_pot(4,i) = 0.5d0*frequency*frequency*q*q
    q=q+dq
  enddo ! i

  

  ! Construct and diagonalise Hamiltonian matrix
  Hamiltonian=0.d0
  ! Harmonic diagonal part: On diagonal we have SHO eigenvalues.  
  do i=1,Nbasis
    Hamiltonian(i,i)=Hamiltonian(i,i)+(real(i-1,dp)+0.5d0)*basis_frequency
  enddo ! i
  ! Add on matrix elements of anharmonic potential w.r.t. SHO eigenfunctions.
  do i=1,Nbasis
    do j=1,Nbasis
      dq=2.d0*max_amplitude/integration_points
      mat_element=0.d0
      do k=1,integration_points
        mat_element=mat_element+basis(k,i)*basis(k,j)*potential(k)*dq
        !write(*,*)mat_element,basis(k,i),basis(k,j),potential(k),dq
      enddo ! k
      Hamiltonian(i,j)=Hamiltonian(i,j)+mat_element
    enddo ! j
  enddo ! i
  output%hamiltonian = Hamiltonian
  ! Diagonalise Hamiltonian matrix
  lwork=3*Nbasis-1
  allocate(work(lwork))
  call dsyev('V','U',Nbasis,Hamiltonian(1,1),Nbasis,E(1),work(1),lwork,info)
  deallocate(work)

  ! Output anharmonic eigenvalues and eigenvectors
  do i=1,Nbasis
    output%eigenvals(1,i) = i
    output%eigenvals(2,i) = frequency*(real(i-1,dp)+0.5d0)*eV
    output%eigenvals(3,i) = E(i)*eV
    do j=1,Nbasis
      output%eigenvecs(1,(i-1)*Nbasis+j) = i
      output%eigenvecs(2,(i-1)*Nbasis+j) = j
      output%eigenvecs(3,(i-1)*Nbasis+j) = Hamiltonian(j,i)
    enddo ! j
  enddo ! i
end function
end module
