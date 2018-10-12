! ======================================================================
! Calculate the LO/TO correction to the energy and force
!    at a given displacement,
!    using the electric field response calculated by DFPT.
! ======================================================================
module loto_splitting_module
  use utils_module
  
  use normal_mode_module
  use structure_module
  implicit none
  
  private
  
  public :: LotoCorrection
  
  type, extends(NoDefaultConstructor) :: LotoCorrection
    real(dp)                      :: energy
    type(RealVector), allocatable :: force(:)
  end type
  
  interface LotoCorrection
    module procedure new_LotoCorrection
    module procedure new_LotoCorrection_dfpt
  end interface
contains

! Constructor.
function new_LotoCorrection(energy,force) result(this)
  implicit none
  
  real(dp),         intent(in) :: energy
  type(RealVector), intent(in) :: force(:)
  type(LotoCorrection)         :: this
  
  this%energy = energy
  this%force = force
end function

! N.B. The force correction is only valid under the assumption that
!    the Born effective charges and the permitivity are constants.
! born_charges are the Born effective charges.
! permittivity is the Dielectric permittivity.
! loto_direction is the direction along which q approaches 0.
function new_LotoCorrection_dfpt(born_charges,permittivity,loto_direction, &
   & displacement,structure) result(this)
  implicit none
  
  type(RealMatrix),            intent(in) :: born_charges(:)
  type(RealMatrix),            intent(in) :: permittivity
  type(FractionVector),        intent(in) :: loto_direction
  type(CartesianDisplacement), intent(in) :: displacement
  type(StructureData),         intent(in) :: structure
  type(LotoCorrection)                    :: this
  
  ! Electric susceptibility and unit cell polarisation.
  type(RealMatrix) :: susceptibility
  type(RealVector) :: polarisation
  
  ! Normalised LO/TO direction.
  type(RealVector) :: direction
  
  ! Output variables.
  real(dp)                      :: energy
  type(RealVector), allocatable :: force(:)
  
  ! Temporary variables.
  type(RealMatrix) :: identity
  type(RealMatrix) :: a
  type(RealMatrix) :: b
  type(RealMatrix) :: c
  
  ! Normalise the LO/TO direction.
  direction = dblevec(loto_direction)
  direction = direction / l2_norm(direction)
  
  ! Calculate the susceptibility from the permittivity.
  identity = dblemat(make_identity_matrix(3))
  susceptibility = permittivity - identity
  
  ! If S is susceptibility,
  !    V is the unit cell volume and
  !    |q> is the normalised LO/TO vector, then:
  ! A = -4pi/V * |q><q|
  ! B = (I-A.S)^-1 . A
  ! C = B + 1/2 * (B^T).S.B
  a = (-4*PI/structure%volume)  * outer_product(direction, direction)
  b = invert(identity - a*susceptibility) * a
  c = b + 0.5_dp*transpose(b)*susceptibility*b
  
  ! Calculate polarisation, energy correction and force correction.
  polarisation = sum(born_charges*displacement%vectors)
  energy = - polarisation * c * polarisation
  force = polarisation * (c + transpose(c)) * born_charges
  
  ! Construct output.
  this = LotoCorrection(energy,force)
end function

! Calculate the dynamical matrix correction.
function dynamical_matrix_correction(born_charges,permittivity, &
   & loto_direction,structure) result(output)
  implicit none
  
  type(RealMatrix),     intent(in) :: born_charges(:)
  type(RealMatrix),     intent(in) :: permittivity
  type(FractionVector), intent(in) :: loto_direction
  type(StructureData),  intent(in) :: structure
  type(RealMatrix), allocatable    :: output(:,:)
  
  ! Normalised LO/TO direction.
  type(RealVector) :: direction
  
  ! Temporary variables.
  real(dp)                      :: scaling
  type(RealVector), allocatable :: a(:)
  integer                       :: i,j,ialloc
  
  ! If Z(i) is the Born effective charge on atom i,
  !    |q>  is the normalised LO/TO direction,
  !    V    is the volume of the primitive cell, and
  !    e    is the permittivity,
  ! The correction to the dynamical matrix D(i,j) is given by:
  !
  ! D(i,j) = (Z(i)^T)|q><q|Z(j) * 4*pi/(V<q|e|q>)

  ! Normalise the LO/TO direction.
  direction = dblevec(loto_direction)
  direction = direction / l2_norm(direction)
  
  ! Calculate the scaling factor, 4*pi/(V<q|e|q>).
  scaling = 4*PI / (structure%volume*direction*permittivity*direction)
  
  ! Calculate (Z(i)^T)|q>.
  a = direction * born_charges
  
  ! Calculate the dynamical matrix correction.
  allocate( output(size(born_charges),size(born_charges)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(born_charges)
    do j=1,size(born_charges)
      output(i,j) = outer_product(a(i), a(j)) * scaling
    enddo
  enddo
end function
end module
