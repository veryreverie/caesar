! ======================================================================
! A pair of HarmonicStateComplex.
! ======================================================================
module caesar_harmonic_braket_complex_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  
  use caesar_harmonic_state_complex_module
  implicit none
  
  private
  
  public :: HarmonicBraKetComplex
  
  type, extends(SubspaceBraKet) :: HarmonicBraKetComplex
    type(HarmonicStateComplex), pointer :: bra_
    type(HarmonicStateComplex), pointer :: ket_
    
    real(dp), private :: frequency_
    real(dp), private :: log_2nw_
    integer,  private :: maximum_power_
    integer,  private :: expansion_order_
  contains
    procedure, public :: set_bra_pointer => &
                       & set_bra_pointer_HarmonicBraKetComplex
    procedure, public :: set_ket_pointer => &
                       & set_ket_pointer_HarmonicBraKetComplex
    
    procedure, public :: finite_overlap => &
                       & finite_overlap_HarmonicBraKetComplex
    procedure, public :: inner_product => &
                       & inner_product_HarmonicBraKetComplex
    procedure, public :: integrate => &
                       & integrate_HarmonicBraKetComplex
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_HarmonicBraKetComplex
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_HarmonicBraKetComplex
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_HarmonicBraKetComplex
  end type
  
  interface HarmonicBraKetComplex
    module function new_HarmonicBraKetComplex(subspace_id,mode_ids, &
       & paired_mode_ids,frequency,supercell_size,maximum_power,    &
       & expansion_order) result(this) 
      integer,  intent(in)        :: subspace_id
      integer,  intent(in)        :: mode_ids(:)
      integer,  intent(in)        :: paired_mode_ids(:)
      real(dp), intent(in)        :: frequency
      integer,  intent(in)        :: supercell_size
      integer,  intent(in)        :: maximum_power
      integer,  intent(in)        :: expansion_order
      type(HarmonicBraKetComplex) :: this
    end function
  end interface
  
  interface
    module subroutine set_bra_pointer_HarmonicBraKetComplex(this,bra) 
      class(HarmonicBraKetComplex), intent(inout)      :: this
      class(SubspaceState),         intent(in), target :: bra
    end subroutine
  end interface
  
  interface
    module subroutine set_ket_pointer_HarmonicBraKetComplex(this,ket) 
      class(HarmonicBraKetComplex), intent(inout)      :: this
      class(SubspaceState),         intent(in), target :: ket
    end subroutine
  end interface
  
  interface
    impure elemental module function finite_overlap_HarmonicBraKetComplex(this,anharmonic_data) result(output) 
      class(HarmonicBraKetComplex), intent(in) :: this
      type(AnharmonicData),         intent(in) :: anharmonic_data
      logical                                  :: output
    end function
  end interface
  
  interface
    impure elemental module function inner_product_HarmonicBraKetComplex(this,anharmonic_data) result(output) 
      class(HarmonicBraKetComplex), intent(in) :: this
      type(AnharmonicData),         intent(in) :: anharmonic_data
      real(dp)                                 :: output
    end function
  end interface
  
  interface
    impure elemental module function integrate_HarmonicBraKetComplex(this, &
       & monomial,anharmonic_data) result(output) 
      class(HarmonicBraKetComplex), intent(in) :: this
      type(SparseMonomial),         intent(in) :: monomial
      type(AnharmonicData),         intent(in) :: anharmonic_data
      complex(dp)                              :: output
    end function
  end interface
  
  interface
    impure elemental module function kinetic_energy_HarmonicBraKetComplex(this,anharmonic_data) result(output) 
      class(HarmonicBraKetComplex), intent(in) :: this
      type(AnharmonicData),         intent(in) :: anharmonic_data
      real(dp)                                 :: output
    end function
  end interface
  
  interface
    impure elemental module function harmonic_potential_energy_HarmonicBraKetComplex(   this,anharmonic_data) result(output) 
      class(HarmonicBraKetComplex), intent(in) :: this
      type(AnharmonicData),         intent(in) :: anharmonic_data
      real(dp)                                 :: output
    end function
  end interface
  
  interface
    impure elemental module function kinetic_stress_HarmonicBraKetComplex(this,stress_prefactors,anharmonic_data) result(output) 
      class(HarmonicBraKetComplex), intent(in) :: this
      type(StressPrefactors),       intent(in) :: stress_prefactors
      type(AnharmonicData),         intent(in) :: anharmonic_data
      type(RealMatrix)                         :: output
    end function
  end interface
end module
