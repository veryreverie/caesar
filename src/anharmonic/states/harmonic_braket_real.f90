! ======================================================================
! A pair of HarmonicStateReal.
! ======================================================================
module caesar_harmonic_braket_real_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  
  use caesar_harmonic_state_real_module
  implicit none
  
  private
  
  public :: HarmonicBraKetReal
  
  type, extends(SubspaceBraKet) :: HarmonicBraKetReal
    type(HarmonicStateReal), pointer :: bra_
    type(HarmonicStateReal), pointer :: ket_
    
    real(dp), private :: frequency_
    real(dp), private :: log_2nw_
    integer,  private :: maximum_power_
    integer,  private :: expansion_order_
  contains
    procedure, public :: set_bra_pointer => set_bra_pointer_HarmonicBraKetReal
    procedure, public :: set_ket_pointer => set_ket_pointer_HarmonicBraKetReal
    
    procedure, public :: finite_overlap => &
                       & finite_overlap_HarmonicBraKetReal
    procedure, public :: inner_product => &
                       & inner_product_HarmonicBraKetReal
    procedure, public :: integrate => &
                       & integrate_HarmonicBraKetReal
    procedure, public :: kinetic_energy => &
                       & kinetic_energy_HarmonicBraKetReal
    procedure, public :: harmonic_potential_energy => &
                       & harmonic_potential_energy_HarmonicBraKetReal
    procedure, public :: kinetic_stress => &
                       & kinetic_stress_HarmonicBraKetReal
  end type
  
  interface HarmonicBraKetReal
    module function new_HarmonicBraKetReal(subspace_id,mode_ids, &
       & paired_mode_ids,frequency,supercell_size,maximum_power, &
       & expansion_order) result(this) 
      integer,  intent(in)     :: subspace_id
      integer,  intent(in)     :: mode_ids(:)
      integer,  intent(in)     :: paired_mode_ids(:)
      real(dp), intent(in)     :: frequency
      integer,  intent(in)     :: supercell_size
      integer,  intent(in)     :: maximum_power
      integer,  intent(in)     :: expansion_order
      type(HarmonicBraKetReal) :: this
    end function
  end interface
  
  interface
    module subroutine set_bra_pointer_HarmonicBraKetReal(this,bra) 
      class(HarmonicBraKetReal), intent(inout)      :: this
      class(SubspaceState),      intent(in), target :: bra
    end subroutine
  end interface
  
  interface
    module subroutine set_ket_pointer_HarmonicBraKetReal(this,ket) 
      class(HarmonicBraKetReal), intent(inout)      :: this
      class(SubspaceState),      intent(in), target :: ket
    end subroutine
  end interface
  
  interface
    impure elemental module function finite_overlap_HarmonicBraKetReal(this, &
       & anharmonic_data) result(output) 
      class(HarmonicBraKetReal), intent(in) :: this
      type(AnharmonicData),      intent(in) :: anharmonic_data
      logical                               :: output
    end function
  end interface
  
  interface
    impure elemental module function inner_product_HarmonicBraKetReal(this, &
       & anharmonic_data) result(output) 
      class(HarmonicBraKetReal), intent(in) :: this
      type(AnharmonicData),      intent(in) :: anharmonic_data
      real(dp)                              :: output
    end function
  end interface
  
  interface
    impure elemental module function integrate_HarmonicBraKetReal(this, &
       & monomial,anharmonic_data) result(output) 
      class(HarmonicBraKetReal), intent(in) :: this
      type(SparseMonomial),      intent(in) :: monomial
      type(AnharmonicData),      intent(in) :: anharmonic_data
      complex(dp)                           :: output
    end function
  end interface
  
  interface
    impure elemental module function kinetic_energy_HarmonicBraKetReal(this, &
       & anharmonic_data) result(output) 
      class(HarmonicBraKetReal), intent(in) :: this
      type(AnharmonicData),      intent(in) :: anharmonic_data
      real(dp)                              :: output
    end function
  end interface
  
  interface
    impure elemental module function harmonic_potential_energy_HarmonicBraKetReal(this,anharmonic_data) result(output) 
      class(HarmonicBraKetReal), intent(in) :: this
      type(AnharmonicData),      intent(in) :: anharmonic_data
      real(dp)                              :: output
    end function
  end interface
  
  interface
    impure elemental module function kinetic_stress_HarmonicBraKetReal(this, &
       & stress_prefactors,anharmonic_data) result(output) 
      class(HarmonicBraKetReal), intent(in) :: this
      type(StressPrefactors),    intent(in) :: stress_prefactors
      type(AnharmonicData),      intent(in) :: anharmonic_data
      type(RealMatrix)                      :: output
    end function
  end interface
end module
