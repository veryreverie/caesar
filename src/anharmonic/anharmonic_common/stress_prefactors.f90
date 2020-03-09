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
    integer,          allocatable, private :: pair_locs_(:)
    type(RealMatrix), allocatable, private :: prefactors_(:,:)
  contains
    generic,   public  :: prefactor => &
                        & prefactor_StressPrefactors_modes, &
                        & prefactor_StressPrefactors_ids
    procedure, private :: prefactor_StressPrefactors_modes
    procedure, private :: prefactor_StressPrefactors_ids
    
    procedure, public :: average_prefactor
  end type
  
  interface StressPrefactors
    module procedure new_StressPrefactors
    module procedure new_StressPrefactors_subspace
  end interface
contains

! Constructor.
function new_StressPrefactors(subspace_id,mode_ids,pair_locs,prefactors) &
   & result(this) 
  implicit none
  
  integer,          intent(in) :: subspace_id
  integer,          intent(in) :: mode_ids(:)
  integer,          intent(in) :: pair_locs(:)
  type(RealMatrix), intent(in) :: prefactors(:,:)
  type(StressPrefactors)       :: this
  
  this%subspace_id_ = subspace_id
  this%mode_ids_    = mode_ids
  this%pair_locs_   = pair_locs
  this%prefactors_  = prefactors
end function

! Calculates the stress prefactors for the modes in the subspace.
function new_StressPrefactors_subspace(subspace,modes) result(this)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: subspace
  type(ComplexMode),        intent(in) :: modes(:)
  type(StressPrefactors)               :: this
  
  type(ComplexMode), allocatable :: subspace_modes(:)
  
  integer, allocatable :: pair_locs(:)
  
  type(RealMatrix), allocatable :: prefactors(:,:)
  
  integer :: i,j,ialloc
  
  subspace_modes = subspace%modes(modes)
  
  allocate( prefactors(size(subspace_modes),size(subspace_modes)), &
          & stat=ialloc); call err(ialloc)
  prefactors = dblemat(zeroes(3,3))
  do i=1,size(subspace_modes)
    do j=1,size(subspace_modes)
      if (subspace_modes(i)%qpoint_id==subspace_modes(j)%paired_qpoint_id) then
        prefactors(i,j) = subspace_modes(i)%stress_prefactor(subspace_modes(j))
      endif
    enddo
  enddo
  
  pair_locs = [( first(subspace_modes%id==subspace_modes(i)%paired_id), &
               & i=1,                                                   &
               & size(subspace_modes)                                   )]
  
  this = StressPrefactors( subspace%id,       &
                         & subspace_modes%id, &
                         & pair_locs,         &
                         & prefactors         )
end function

! Return the prefactor for a given pair of modes.
impure elemental function prefactor_StressPrefactors_modes(this,mode1,mode2) &
   & result(output)
  implicit none
  
  class(StressPrefactors), intent(in) :: this
  type(ComplexMode),       intent(in) :: mode1
  type(ComplexMode),       intent(in) :: mode2
  type(RealMatrix)                    :: output
  
  output = this%prefactors_( first(this%mode_ids_==mode1%id), &
                           & first(this%mode_ids_==mode2%id)  )
end function

impure elemental function prefactor_StressPrefactors_ids(this,mode1,mode2) &
   & result(output)
  implicit none
  
  class(StressPrefactors), intent(in) :: this
  integer,                 intent(in) :: mode1
  integer,                 intent(in) :: mode2
  type(RealMatrix)                    :: output
  
  output = this%prefactors_( first(this%mode_ids_==mode1), &
                           & first(this%mode_ids_==mode2)  )
end function

! Return the average prefactor.
! Useful for harmonic calculations, where all modes are treated equally.
impure elemental function average_prefactor(this) result(output)
  implicit none
  
  class(StressPrefactors), intent(in) :: this
  type(RealMatrix)                    :: output
  
  integer :: no_modes
  
  integer :: i
  
  no_modes = size(this%mode_ids_)
  
  output = sum([(this%prefactors_(i,this%pair_locs_(i)), i=1, no_modes)]) &
       & / no_modes
end function
end module
