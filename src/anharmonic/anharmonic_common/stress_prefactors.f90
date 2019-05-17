! ======================================================================
! Prefactors for stress calculations
! ======================================================================
! The stress prefactor between mode ui at q-point q and uj at q-point -q is
! I_{q,i,j} = sum_k ui_k ^ uj_k*, where ui_k is the mode at atom k.
module stress_prefactors_module
  use common_module
  implicit none
  
  private
  
  public :: StressPrefactors
  
  type, extends(NoDefaultConstructor) :: StressPrefactors
    integer,                       private :: subspace_id_
    integer,          allocatable, private :: mode_ids_(:)
    type(RealMatrix), allocatable, private :: prefactors_(:,:)
  contains
    procedure, public :: prefactor => prefactor_StressPrefactors
  end type
  
  interface StressPrefactors
    module procedure new_StressPrefactors
    module procedure new_StressPrefactors_subspace
  end interface
contains

! Constructor.
function new_StressPrefactors(subspace_id,mode_ids,prefactors) result(this)
  implicit none
  
  integer,          intent(in) :: subspace_id
  integer,          intent(in) :: mode_ids(:)
  type(RealMatrix), intent(in) :: prefactors(:,:)
  type(StressPrefactors)       :: this
  
  this%subspace_id_ = subspace_id
  this%mode_ids_    = mode_ids
  this%prefactors_  = prefactors
end function

! ----------------------------------------------------------------------
! Calculates the stress prefactors for the modes in the subspace.
! ----------------------------------------------------------------------
function new_StressPrefactors_subspace(subspace,modes) result(this)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: subspace
  type(ComplexMode),        intent(in) :: modes(:)
  type(StressPrefactors)               :: this
  
  type(ComplexMode), allocatable :: subspace_modes(:)
  
  type(RealMatrix), allocatable :: prefactors(:,:)
  
  integer :: i,j,k,ialloc
  
  subspace_modes = subspace%modes(modes)
  
  allocate( prefactors(size(subspace_modes),size(subspace_modes)), &
          & stat=ialloc); call err(ialloc)
  prefactors = dblemat(zeroes(3,3))
  do i=1,size(subspace_modes)
    do j=1,size(subspace_modes)
      if (subspace_modes(i)%qpoint_id==subspace_modes(j)%paired_qpoint_id) then
        do k=1,size(subspace_modes(i)%unit_vector)
          prefactors(i,j) = prefactors(i,j)                            &
                        & + real(outer_product(                        &
                        &      subspace_modes(i)%unit_vector(k),       &
                        &      conjg(subspace_modes(j)%unit_vector(k)) ))
        enddo
      endif
    enddo
  enddo
  
  this = StressPrefactors(subspace%id, subspace_modes%id, prefactors)
end function

! --------------------------------------------------
! Return the prefactor for a given pair of modes.
! --------------------------------------------------
impure elemental function prefactor_StressPrefactors(this,mode1,mode2) &
   & result(output)
  implicit none
  
  class(StressPrefactors), intent(in) :: this
  type(ComplexMode),       intent(in) :: mode1
  type(ComplexMode),       intent(in) :: mode2
  type(RealMatrix)                    :: output
  
  output = this%prefactors_( first(this%mode_ids_==mode1%id), &
                           & first(this%mode_ids_==mode2%id)  )
end function
end module
