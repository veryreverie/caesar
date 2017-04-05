module vscf_1d_module
  use constants, only : dp
  implicit none
  
  ! Holds anharmonic potential data
  type Potential
    real(dp) :: q
    real(dp) :: harmonic   ! harmonic f(q) = 0.5*(frequency*q)^2
    real(dp) :: anharmonic ! anharmonic f(q)
    real(dp) :: difference ! anharmonic f(q) - harmonic f(q)
  end type
  
  ! Holds output data for vscf_1d
  type VscfData
    type(Potential), allocatable :: anh_pot(:)
    real(dp),        allocatable :: hamiltonian(:,:)
    real(dp),        allocatable :: harmonic(:)
    real(dp),        allocatable :: eigenvals(:)
    real(dp),        allocatable :: eigenvecs(:,:)
  end type

  interface new
    module procedure new_VscfData
  end interface
  
  interface drop
    module procedure drop_VscfData
  end interface
  
contains

function vscf_1d(frequency, potential, Nbasis) result(output)
  use constants,      only : dp, pi
  use linear_algebra, only : RealEigenstuff, calculate_eigenstuff, size
  implicit none
  
  real(dp), intent(in) :: frequency
  real(dp), intent(in) :: potential(:,:) ! Anharmonic potential, {q,V(q)}
  integer,  intent(in) :: Nbasis
  type(VscfData)       :: output
  
  integer  :: Npoints   ! Size(potential,2)
  
  ! Working variables
  integer :: i,j,k
  real(dp) :: bfp,q,dq,basis_frequency
  real(dp),allocatable :: basis(:,:)
  type(RealEigenstuff) :: estuff ! evals and evecs of Hamiltonian
  
  Npoints = size(potential,2)
  dq = (potential(1,Npoints)-potential(1,1))/Npoints
  
  ! Allocate output
  call new(output,Npoints,Nbasis)

  ! Calculate basis functions
  basis_frequency=dabs(frequency)
  bfp=(basis_frequency/pi)**0.25
  
  allocate(basis(Npoints,Nbasis))
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
    output%anh_pot(i)%anharmonic = potential(2,i)
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
                              & + basis(k,i)*basis(k,j)*potential(2,k)*dq
      enddo
    enddo
  enddo
  
  estuff = calculate_eigenstuff(output%hamiltonian)

  ! Output anharmonic eigenvalues and eigenvectors
  do i=1,Nbasis
    output%harmonic(i) = frequency*(i-0.5d0)
    output%eigenvals(i) = estuff%evals(i)
  enddo
  
  output%eigenvecs = estuff%evecs
  
end function

! allocate(VscfData)
pure subroutine new_VscfData(this,Npoints,Nbasis)
  implicit none
  
  type(VscfData), intent(out) :: this
  integer,        intent(in)  :: Npoints
  integer,        intent(in)  :: Nbasis
  
  allocate(this%anh_pot(Npoints))
  allocate(this%hamiltonian(Nbasis,Nbasis))
  allocate(this%harmonic(Nbasis))
  allocate(this%eigenvals(Nbasis))
  allocate(this%eigenvecs(Nbasis,Nbasis))
end subroutine

! deallocate(VscfData)
pure subroutine drop_VscfData(this)
  implicit none
  
  type(VscfData), intent(inout) :: this
  
  deallocate(this%anh_pot)
  deallocate(this%hamiltonian)
  deallocate(this%harmonic)
  deallocate(this%eigenvals)
  deallocate(this%eigenvecs)
end subroutine

end module
