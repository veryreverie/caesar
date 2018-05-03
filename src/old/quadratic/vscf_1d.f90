module vscf_1d_module
  use common_module
  implicit none
  
  ! Holds anharmonic potential data
  type Potential
    real(dp) :: u
    real(dp) :: harmonic   ! harmonic f(u) = 0.5*(frequency*u)^2
    real(dp) :: anharmonic ! anharmonic f(u)
    real(dp) :: difference ! anharmonic f(u) - harmonic f(u)
  end type
  
  ! Holds output data for vscf_1d
  type VscfData
    type(Potential), allocatable :: anh_pot(:)
    real(dp),        allocatable :: hamiltonian(:,:)
    real(dp),        allocatable :: harmonic(:)
    real(dp),        allocatable :: eigenvals(:)
    real(dp),        allocatable :: eigenvecs(:,:)
  end type

  interface VscfData
    module procedure new_VscfData
  end interface
  
contains

function vscf_1d(frequency, potential, Nbasis) result(output)
  implicit none
  
  real(dp), intent(in) :: frequency
  real(dp), intent(in) :: potential(:,:) ! Anharmonic potential, {u,V(u)}
  integer,  intent(in) :: Nbasis
  type(VscfData)       :: output
  
  integer  :: Npoints   ! Size(potential,2)
  
  ! Working variables
  integer :: i,j,k
  real(dp) :: bfp,u,du,basis_frequency
  real(dp),allocatable :: basis(:,:)
  type(RealEigenstuff) :: estuff ! evals and evecs of Hamiltonian
  
  Npoints = size(potential,2)
  du = (potential(1,Npoints)-potential(1,1))/Npoints
  
  ! Allocate output
  output = VscfData(Npoints,Nbasis)

  ! Calculate basis functions (1-d Harmonic Oscillator wavefunctions).
  basis_frequency = abs(frequency)
  bfp=(basis_frequency/PI)**0.25
  
  allocate(basis(Npoints,Nbasis))
  do i=1,Npoints
    u = potential(1,i)
    basis(i,1) = bfp*exp(-u*u*basis_frequency/2)
    basis(i,2) = sqrt(2*basis_frequency)*u*basis(i,1)
    do j=3,Nbasis
      basis(i,j) = sqrt(2*basis_frequency/(j-1))*u * basis(i,j-1) &
               & - sqrt(dble(j-2)/(j-1))           * basis(i,j-2)
    enddo ! j
  enddo ! i

  ! write anh_pot
  do i=1,Npoints
    u = potential(1,i)
    output%anh_pot(i)%u = u
    output%anh_pot(i)%harmonic = 0.5d0*frequency*frequency*u*u
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
                              & + basis(k,i)*basis(k,j)*potential(2,k)*du
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
function new_VscfData(Npoints,Nbasis) result(this)
  implicit none
  
  integer, intent(in) :: Npoints
  integer, intent(in) :: Nbasis
  type(VscfData)      :: this
  
  integer :: ialloc
  
  allocate( this%anh_pot(Npoints),           &
          & this%hamiltonian(Nbasis,Nbasis), &
          & this%harmonic(Nbasis),           &
          & this%eigenvals(Nbasis),          &
          & this%eigenvecs(Nbasis,Nbasis),   &
          & stat=ialloc); call err(ialloc)
end function

end module
