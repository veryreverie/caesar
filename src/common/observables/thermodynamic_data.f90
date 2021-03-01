! ======================================================================
! Calculates thermodynamic properties, either for energy spectra or for
!    uncoupled harmonic oscillators.
! ======================================================================
! N.B. Stress is stored with the same extensivity as U, F, S, H and G.
module caesar_thermodynamic_data_module
  use caesar_utils_module
  implicit none
  
  private
  
  public :: ThermodynamicData
  public :: calculate_bose_factor
  public :: calculate_state_weight
  
  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: operator(/)
  
  public :: sum
  public :: min
  
  type, extends(Stringable) :: ThermodynamicData
    ! T, U, F and S.
    real(dp) :: thermal_energy
    real(dp) :: energy
    real(dp) :: free_energy
    real(dp) :: entropy
    
    ! P, V, H and G.
    type(RealMatrix), allocatable :: stress
    real(dp),         allocatable :: primitive_volume
    real(dp),         allocatable :: enthalpy
    real(dp),         allocatable :: gibbs
  contains
    procedure, public :: set_stress => set_stress_ThermodynamicData
    
    procedure, public :: add_energy  => add_energy_ThermodynamicData
    procedure, public :: add_entropy => add_entropy_ThermodynamicData
    procedure, public :: add_stress  => add_stress_ThermodynamicData
    
    ! I/O.
    procedure, public :: read  => read_ThermodynamicData
    procedure, public :: write => write_ThermodynamicData
  end type
  
  interface ThermodynamicData
    impure elemental module function new_ThermodynamicData(thermal_energy,  &
       & energy,free_energy,entropy,stress,primitive_volume,enthalpy,gibbs) &
       & result(this) 
      real(dp),         intent(in)           :: thermal_energy
      real(dp),         intent(in)           :: energy
      real(dp),         intent(in)           :: free_energy
      real(dp),         intent(in)           :: entropy
      type(RealMatrix), intent(in), optional :: stress
      real(dp),         intent(in), optional :: primitive_volume
      real(dp),         intent(in), optional :: enthalpy
      real(dp),         intent(in), optional :: gibbs
      type(ThermodynamicData)                :: this
    end function
  
    module function new_ThermodynamicData_spectrum(thermal_energy,energies, &
       & stresses,primitive_volume) result(output) 
      real(dp),         intent(in)           :: thermal_energy
      real(dp),         intent(in)           :: energies(:)
      type(RealMatrix), intent(in), optional :: stresses(:)
      real(dp),         intent(in), optional :: primitive_volume
      type(ThermodynamicData)                :: output
    end function
  
    module function new_ThermodynamicData_harmonic(thermal_energy,frequency, &
       & stress_prefactor,potential_stress,primitive_volume) result(output) 
      real(dp),         intent(in)           :: thermal_energy ! T.
      real(dp),         intent(in)           :: frequency      ! w.
      type(RealMatrix), intent(in), optional :: stress_prefactor
      type(RealMatrix), intent(in), optional :: potential_stress
      real(dp),         intent(in), optional :: primitive_volume
      type(ThermodynamicData) :: output
    end function
  end interface
  
  interface
    ! Set the stress of a ThermodynamicData.
    impure elemental module subroutine set_stress_ThermodynamicData(this, &
       & stress,primitive_volume) 
      class(ThermodynamicData), intent(inout) :: this
      type(RealMatrix),         intent(in)    :: stress
      real(dp),                 intent(in)    :: primitive_volume
    end subroutine
  end interface
  
  interface
    ! Add energy, entropy or stress.
    impure elemental module subroutine add_energy_ThermodynamicData(this, &
       & energy) 
      class(ThermodynamicData), intent(inout) :: this
      real(dp),                 intent(in)    :: energy
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine add_entropy_ThermodynamicData(this, &
       & entropy) 
      class(ThermodynamicData), intent(inout) :: this
      real(dp),                 intent(in)    :: entropy
    end subroutine
  end interface
  
  interface
    impure elemental module subroutine add_stress_ThermodynamicData(this, &
       & stress) 
      class(ThermodynamicData), intent(inout) :: this
      type(RealMatrix),         intent(in)    :: stress
    end subroutine
  end interface
  
  interface operator(+)
    ! ----------------------------------------------------------------------
    ! Algebra with ThermodynamicData.
    ! N.B. + and - assume that the thermal_energy and volume are the same
    !    for both arguments.
    ! ----------------------------------------------------------------------
    impure elemental module function add_ThermodynamicData_ThermodynamicData(this,that) result(output) 
      type(ThermodynamicData), intent(in) :: this
      type(ThermodynamicData), intent(in) :: that
      type(ThermodynamicData)             :: output
    end function
  end interface
  
  interface operator(-)
    impure elemental module function negative_ThermodynamicData(input) &
       & result(output) 
      type(ThermodynamicData), intent(in) :: input
      type(ThermodynamicData)             :: output
    end function
  
    impure elemental module function subtract_ThermodynamicData_ThermodynamicData(this,that) result(output) 
      type(ThermodynamicData), intent(in) :: this
      type(ThermodynamicData), intent(in) :: that
      type(ThermodynamicData)             :: output
    end function
  end interface
  
  interface operator(*)
    impure elemental module function multiply_ThermodynamicData_integer(this,that) result(output) 
      type(ThermodynamicData), intent(in) :: this
      integer,                 intent(in) :: that
      type(ThermodynamicData)             :: output
    end function
  
    impure elemental module function multiply_integer_ThermodynamicData(this,that) result(output) 
      integer,                 intent(in) :: this
      type(ThermodynamicData), intent(in) :: that
      type(ThermodynamicData)             :: output
    end function
  
    impure elemental module function multiply_ThermodynamicData_real(this, &
       & that) result(output) 
      type(ThermodynamicData), intent(in) :: this
      real(dp),                intent(in) :: that
      type(ThermodynamicData)             :: output
    end function
  
    impure elemental module function multiply_real_ThermodynamicData(this, &
       & that) result(output) 
      real(dp),                intent(in) :: this
      type(ThermodynamicData), intent(in) :: that
      type(ThermodynamicData)             :: output
    end function
  end interface
  
  interface operator(/)
    impure elemental module function divide_ThermodynamicData_integer(this, &
       & that) result(output) 
      type(ThermodynamicData), intent(in) :: this
      integer,                 intent(in) :: that
      type(ThermodynamicData)             :: output
    end function
  
    impure elemental module function divide_ThermodynamicData_real(this,that) &
       & result(output) 
      type(ThermodynamicData), intent(in) :: this
      real(dp),                intent(in) :: that
      type(ThermodynamicData)             :: output
    end function
  end interface
  
  interface sum
    module function sum_ThermodynamicData(input) result(output) 
      type(ThermodynamicData), intent(in) :: input(:)
      type(ThermodynamicData)             :: output
    end function
  end interface
  
  interface min
    module function min_ThermodynamicData(input) result(output) 
      type(ThermodynamicData), intent(in) :: input(:)
      type(ThermodynamicData)             :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns n(T,w) = 1/(e^(w/T)-1), using numerically stable strategies.
    ! ----------------------------------------------------------------------
    impure elemental module function calculate_bose_factor(thermal_energy, &
       & frequency) result(output) 
      real(dp), intent(in) :: thermal_energy
      real(dp), intent(in) :: frequency
      real(dp)             :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Returns P(T,w,n) = (1-e^(-w/T))^d * e^(-nw/T),
    !    using numerically stable strategies.
    ! ----------------------------------------------------------------------
    impure elemental module function calculate_state_weight(thermal_energy, &
       & frequency,occupation,state_dimension) result(output) 
      real(dp), intent(in) :: thermal_energy
      real(dp), intent(in) :: frequency
      integer,  intent(in) :: occupation
      integer,  intent(in) :: state_dimension
      real(dp)             :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_ThermodynamicData(this,input) 
      class(ThermodynamicData), intent(out) :: this
      type(String),             intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_ThermodynamicData(this) result(output) 
      class(ThermodynamicData), intent(in) :: this
      type(String)                         :: output
    end function
  end interface
  
  interface ThermodynamicData
    impure elemental module function new_ThermodynamicData_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(ThermodynamicData)  :: this
    end function
  end interface
end module
