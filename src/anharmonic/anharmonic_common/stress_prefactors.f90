! ======================================================================
! Prefactors for stress calculations
! ======================================================================
! The stress prefactor between mode ui at q-point q and uj at q-point -q is
! I_{q,i,j} = sum_k ui_k ^ uj_k*, where ui_k is the mode at atom k.
module caesar_stress_prefactors_module
  use caesar_common_module
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
    ! Constructor.
    module function new_StressPrefactors(subspace_id,mode_ids,pair_locs, &
       & prefactors) result(this) 
      integer,          intent(in) :: subspace_id
      integer,          intent(in) :: mode_ids(:)
      integer,          intent(in) :: pair_locs(:)
      type(RealMatrix), intent(in) :: prefactors(:,:)
      type(StressPrefactors)       :: this
    end function
  
    ! Calculates the stress prefactors for the modes in the subspace.
    module function new_StressPrefactors_subspace(subspace,modes) result(this) 
      type(DegenerateSubspace), intent(in) :: subspace
      type(ComplexMode),        intent(in) :: modes(:)
      type(StressPrefactors)               :: this
    end function
  end interface
  
  interface
    ! Return the prefactor for a given pair of modes.
    impure elemental module function prefactor_StressPrefactors_modes(this, &
       & mode1,mode2) result(output) 
      class(StressPrefactors), intent(in) :: this
      type(ComplexMode),       intent(in) :: mode1
      type(ComplexMode),       intent(in) :: mode2
      type(RealMatrix)                    :: output
    end function
  end interface
  
  interface
    impure elemental module function prefactor_StressPrefactors_ids(this, &
       & mode1,mode2) result(output) 
      class(StressPrefactors), intent(in) :: this
      integer,                 intent(in) :: mode1
      integer,                 intent(in) :: mode2
      type(RealMatrix)                    :: output
    end function
  end interface
  
  interface
    ! Return the average prefactor.
    ! Useful for harmonic calculations, where all modes are treated equally.
    impure elemental module function average_prefactor(this) result(output) 
      class(StressPrefactors), intent(in) :: this
      type(RealMatrix)                    :: output
    end function
  end interface
end module
