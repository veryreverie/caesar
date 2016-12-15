module vscf_1d_module
  use constants, only : dp
  implicit none
  
  ! holds anharmonic potential data
  type Potential
    real(dp) :: q
    real(dp) :: harmonic   ! harmonic f(q) = 0.5*(frequency*q)^2
    real(dp) :: anharmonic ! anharmonic f(q)
    real(dp) :: difference ! anharmonic f(q) - harmonic f(q)
  end type
  
  ! holds output data for vscf_1d
  type VscfReturn
    type(Potential), allocatable :: anh_pot(:)
    real(dp),        allocatable :: hamiltonian(:,:)
    real(dp),        allocatable :: harmonic(:)
    real(dp),        allocatable :: eigenvals(:)
    real(dp),        allocatable :: eigenvecs(:,:)
  end type
  
contains

function vscf_1d(frequency_ev, potential, Nbasis) result(output)
  use constants,      only : dp, pi, eV, thermal
  use linear_algebra, only : Eigenstuff, calculate_eigenstuff, size
  implicit none
  
  real(dp), intent(in) :: frequency_ev   ! frequency in eV
  real(dp), intent(in) :: potential(:,:) ! anharmonic potential, {q,V(q)}
  integer,  intent(in) :: Nbasis
  type(VscfReturn)     :: output
  
  integer  :: Npoints   ! size(potential,2)
  real(dp) :: frequency ! frequency in a.u.
  
  ! Working variables
  integer :: i,j,k
  real(dp) :: bfp,q,dq,basis_frequency
  real(dp),allocatable :: basis(:,:)
  real(dp),allocatable :: Hamiltonian(:,:)
  type(Eigenstuff)     :: estuff ! evals and evecs of Hamiltonian
  
  Npoints = size(potential,2)
  dq = (potential(1,Npoints)-potential(1,1))/Npoints
  
  ! Convert frequency to a.u.
  frequency=frequency_ev/eV 

  ! Allocate output
  allocate(output%anh_pot(Npoints))
  allocate(output%hamiltonian(Nbasis,Nbasis))
  allocate(output%harmonic(Nbasis))
  allocate(output%eigenvals(Nbasis))
  allocate(output%eigenvecs(Nbasis,Nbasis))

  ! Allocate various arrays
  allocate(basis(Npoints,Nbasis))

  ! Calculate basis functions
  basis_frequency=dabs(frequency)
  bfp=(basis_frequency/pi)**0.25
  
  do i=1,Npoints
    q = potential(1,i)
    basis(i,1)=bfp*dexp(-q*q*basis_frequency/2)
    basis(i,2)=dsqrt(2*basis_frequency)*q*basis(i,1)
    do j=3,Nbasis
      basis(i,j) = dsqrt(2*basis_frequency/(j-1))*q * basis(i,j-1) &
               & - dsqrt(dble(j-2)/(j-1))           * basis(i,j-2)
    enddo ! j
  enddo ! i

  ! write anh_pot
  do i=1,Npoints
    q = potential(1,i)
    output%anh_pot(i)%q = q
    output%anh_pot(i)%harmonic = 0.5d0*frequency*frequency*q*q
    output%anh_pot(i)%anharmonic = potential(2,i)/eV
    output%anh_pot(i)%difference = output%anh_pot(i)%anharmonic &
                               & - output%anh_pot(i)%harmonic
  enddo ! i

  ! Construct and diagonalise Hamiltonian matrix
  output%hamiltonian=0.d0
  ! Harmonic diagonal part: On diagonal we have SHO eigenvalues.  
  do i=1,Nbasis
    output%hamiltonian(i,i) = (i-0.5d0)*basis_frequency
  enddo
  
  ! Add on matrix elements of anharmonic potential w.r.t. SHO eigenfunctions.
  do i=1,Nbasis
    do j=1,Nbasis
      do k=1,Npoints
        output%hamiltonian(i,j) = output%hamiltonian(i,j) &
                              & + basis(k,i)*basis(k,j)*potential(2,k)*dq/eV
      enddo
    enddo
  enddo
  
  estuff = calculate_eigenstuff(output%hamiltonian)

  ! Output anharmonic eigenvalues and eigenvectors
  do i=1,Nbasis
    output%harmonic(i) = frequency*(i-0.5d0)*eV
    output%eigenvals(i) = estuff%evals(i)*eV
  enddo
  
  output%eigenvecs = estuff%evecs
  
end function

! deallocate(VscfReturn)
pure subroutine drop(this)
  implicit none
  
  type(VscfReturn), intent(inout) :: this
  
  deallocate(this%anh_pot)
  deallocate(this%hamiltonian)
  deallocate(this%harmonic)
  deallocate(this%eigenvals)
  deallocate(this%eigenvecs)
end subroutine

end module
