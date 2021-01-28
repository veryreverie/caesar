! ======================================================================
! A pair of SubspaceStates, in order to perform operations of the form
!    <bra|X|ket>, for some operator X.
! ======================================================================
! This module exists since Fortran does not support multiple passed-object
!    dummy arguments, so functions which take two states can only be made
!    fast if they are polymorphic over a type which points to both states.
module caesar_subspace_braket_module
  use caesar_common_module
  
  use caesar_stress_prefactors_module
  use caesar_anharmonic_data_module
  use caesar_sparse_monomial_module
  use caesar_subspace_state_module
  implicit none
  
  private
  
  public :: SubspaceBraKet
  
  type, abstract, extends(NoDefaultConstructor) :: SubspaceBraKet
    ! The id of the subspace.
    integer :: subspace_id
    ! The ids and paired ids of the modes across which the bra and ket
    !    operate.
    integer, allocatable :: mode_ids(:)
    integer, allocatable :: paired_mode_ids(:)
  contains
    procedure(set_bra_pointer_SubspaceBraKet), public, deferred :: &
       & set_bra_pointer
    procedure(set_ket_pointer_SubspaceBraKet), public, deferred :: &
       & set_ket_pointer
    
    ! Whether <bra|ket> is finite.
    procedure(finite_overlap_SubspaceBraKet), public, deferred :: &
       & finite_overlap
    
    ! <bra|ket>.
    procedure(inner_product_SubspaceBraKet),  public, deferred :: &
       & inner_product
    
    ! <bra|V|ket>.
    ! This operates on an intent(inout) ComplexMonomial.
    procedure(integrate_SubspaceBraKet),      public, deferred :: &
       & integrate
    
    ! <bra|T|ket>.
    procedure(kinetic_energy_SubspaceBraKet), public, deferred :: &
       & kinetic_energy
    
    ! <bra|V|ket>, where V is the harmonic potential energy.
    procedure(harmonic_potential_energy_SubspaceBraKet), public, deferred :: &
       & harmonic_potential_energy
    
    ! <bra|stress|ket>, where stress is the kinetic stress.
    procedure(kinetic_stress_SubspaceBraKet), public, deferred :: &
       & kinetic_stress
  end type
  
  abstract interface
    subroutine set_bra_pointer_SubspaceBraKet(this,bra)
      import SubspaceBraKet
      import SubspaceState
      implicit none
      
      class(SubspaceBraKet), intent(inout)      :: this
      class(SubspaceState),  intent(in), target :: bra
    end subroutine
    
    subroutine set_ket_pointer_SubspaceBraKet(this,ket)
      import SubspaceBraKet
      import SubspaceState
      implicit none
      
      class(SubspaceBraKet), intent(inout)      :: this
      class(SubspaceState),  intent(in), target :: ket
    end subroutine
    
    impure elemental function finite_overlap_SubspaceBraKet(this, &
       & anharmonic_data) result(output)
      import SubspaceBraKet
      import AnharmonicData
      implicit none
      
      class(SubspaceBraKet), intent(in) :: this
      type(AnharmonicData),  intent(in) :: anharmonic_data
      logical                           :: output
    end function
    
    impure elemental function inner_product_SubspaceBraKet(this, &
       & anharmonic_data) result(output)
      import SubspaceBraKet
      import AnharmonicData
      import dp
      implicit none
      
      class(SubspaceBraKet), intent(in) :: this
      type(AnharmonicData),  intent(in) :: anharmonic_data
      real(dp)                          :: output
    end function
    
    impure elemental function integrate_SubspaceBraKet(this,monomial, &
       & anharmonic_data) result(output)
      import SubspaceBraKet
      import SparseMonomial
      import AnharmonicData
      import dp
      implicit none
      
      class(SubspaceBraKet), intent(in) :: this
      type(SparseMonomial),  intent(in) :: monomial
      type(AnharmonicData),  intent(in) :: anharmonic_data
      complex(dp)                       :: output
    end function
    
    impure elemental function kinetic_energy_SubspaceBraKet(this, &
       & anharmonic_data) result(output)
      import SubspaceBraKet
      import AnharmonicData
      import dp
      implicit none
      
      class(SubspaceBraKet), intent(in) :: this
      type(AnharmonicData),  intent(in) :: anharmonic_data
      real(dp)                          :: output
    end function
    
    impure elemental function harmonic_potential_energy_SubspaceBraKet(this, &
       & anharmonic_data) result(output)
      import SubspaceBraKet
      import AnharmonicData
      import dp
      implicit none
      
      class(SubspaceBraKet), intent(in) :: this
      type(AnharmonicData),  intent(in) :: anharmonic_data
      real(dp)                          :: output
    end function
    
    impure elemental function kinetic_stress_SubspaceBraKet(this, &
       & stress_prefactors,anharmonic_data) result(output)
      import SubspaceBraKet
      import StressPrefactors
      import AnharmonicData
      import RealMatrix
      implicit none
      
      class(SubspaceBraKet),  intent(in) :: this
      type(StressPrefactors), intent(in) :: stress_prefactors
      type(AnharmonicData),   intent(in) :: anharmonic_data
      type(RealMatrix)                   :: output
    end function
  end interface
end module
