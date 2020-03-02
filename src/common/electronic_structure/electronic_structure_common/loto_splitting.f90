! ======================================================================
! Calculate the LO/TO correction to the energy and forces
!    at a given displacement,
!    using the electric field response calculated by DFPT.
! ======================================================================
module loto_splitting_module
  use utils_module
  
  use normal_mode_module
  use structure_module
  use electronic_structure_data_module
  implicit none
  
  private
  
  public :: LotoCorrection
  public :: calculate_loto_correction
  public :: dynamical_matrix_correction
  
  type, extends(NoDefaultConstructor) :: LotoCorrection
    real(dp)             :: energy
    type(CartesianForce) :: forces
  end type
  
  interface LotoCorrection
    module procedure new_LotoCorrection
    module procedure new_LotoCorrection_LinearResponse
  end interface
contains

! Constructor.
function new_LotoCorrection(energy,forces) result(this)
  implicit none
  
  real(dp),             intent(in) :: energy
  type(CartesianForce), intent(in) :: forces
  type(LotoCorrection)             :: this
  
  this%energy = energy
  this%forces = forces
end function

! N.B. The force correction is only valid under the assumption that
!    the Born effective charges and the permitivity are constants.
function new_LotoCorrection_LinearResponse(linear_response,loto_direction, &
   & displacement,structure) result(this)
  implicit none
  
  type(LinearResponse),        intent(in) :: linear_response
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
  type(RealVector), allocatable :: forces(:)
  
  ! Temporary variables.
  type(RealMatrix) :: identity
  type(RealMatrix) :: a
  type(RealMatrix) :: b
  type(RealMatrix) :: c
  
  integer :: i
  
  ! Normalise the LO/TO direction.
  direction = dblevec(loto_direction)
  direction = direction / l2_norm(direction)
  
  ! Calculate the susceptibility from the permittivity.
  identity = dblemat(make_identity_matrix(3))
  susceptibility = linear_response%permittivity - identity
  
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
  polarisation = sum(linear_response%born_charges*displacement%vectors)
  energy = - polarisation * c * polarisation
  forces = [(                                                             &
     & polarisation * (c+transpose(c)) * linear_response%born_charges(i), &
     & i=1,                                                               &
     & size(linear_response%born_charges)                                 )]
  
  ! Construct output.
  this = LotoCorrection(energy, CartesianForce(forces))
end function

! Calculate the dynamical matrix correction.
function dynamical_matrix_correction(linear_response,loto_direction, &
   & structure) result(output)
  implicit none
  
  type(LinearResponse), intent(in) :: linear_response
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
  scaling = 4*PI &
        & / (structure%volume*direction*linear_response%permittivity*direction)
  
  ! Calculate (Z(i)^T)|q>.
  a = direction * linear_response%born_charges
  
  ! Calculate the dynamical matrix correction.
  allocate( output(structure%no_atoms,structure%no_atoms), &
          & stat=ialloc); call err(ialloc)
  do i=1,structure%no_atoms
    do j=1,structure%no_atoms
      output(i,j) = outer_product(a(i), a(j)) * scaling
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Uses calculated LO/TO correction to update electronic structure output.
! ----------------------------------------------------------------------
function calculate_loto_correction(electronic_structure,loto_correction) &
   & result(output)
  implicit none
  
  type(ElectronicStructure), intent(in) :: electronic_structure
  type(LotoCorrection),      intent(in) :: loto_correction
  type(ElectronicStructure)             :: output
  
  real(dp)             :: energy
  type(CartesianForce) :: forces
  
  energy = electronic_structure%energy() + loto_correction%energy
  forces = electronic_structure%forces() + loto_correction%forces
  
  output = ElectronicStructure(energy, forces)
end function
end module
