! ======================================================================
! Routines to calculate a number of harmonic properties of the crystal from
!    the Hessian matrix.
! ======================================================================
module caesar_harmonic_properties_module
  use caesar_utils_module
  
  use caesar_structure_module
  use caesar_normal_mode_module
  use caesar_dynamical_matrices_module
  use caesar_harmonic_stress_module
  
  use caesar_thermodynamic_data_module
  use caesar_qpoint_path_module
  implicit none
  
  private
  
  public :: PhononDispersion
  public :: PhononDos
  
  public :: PathFrequencies
  public :: SampledQpoint
  public :: PdosBin
  
  ! The minor points along the q-point path.
  type, extends(StringsWriteable) :: PathFrequencies
    real(dp)              :: path_fraction
    type(RealVector)      :: qpoint
    real(dp), allocatable :: frequencies(:)
    type(DynamicalMatrix) :: dynamical_matrix
  contains
    procedure, public :: write => write_PathFrequencies
  end type
  
  ! The number of frequencies below the threshold at a given q-point.
  type, extends(Stringable) :: SampledQpoint
    type(RealVector) :: qpoint
    integer          :: no_soft_frequencies
  contains
    procedure, public :: read  => read_SampledQpoint
    procedure, public :: write => write_SampledQpoint
  end type
  
  ! The phonon density of states (PDOS).
  type, extends(Stringable) :: PdosBin
    real(dp) :: min_frequency
    real(dp) :: max_frequency
    real(dp) :: occupation
  contains
    procedure, public :: read  => read_PdosBin
    procedure, public :: write => write_PdosBin
  end type
  
  ! Return types.
  type, extends(NoDefaultConstructor) :: PhononDispersion
    type(QpointPath),                   private :: path_
    type(PathFrequencies), allocatable, private :: frequencies_(:)
  contains
    procedure, public :: path => path_PhononDispersion
    procedure, public :: frequencies => frequencies_PhononDispersion
    procedure, public :: json => json_PhononDispersion
    
    procedure, public :: write_files => write_files_PhononDispersion
  end type
  
  type, extends(NoDefaultConstructor) :: PhononDos
    type(SampledQpoint),     allocatable :: qpoints(:)
    type(PdosBin),           allocatable :: pdos(:)
    type(ThermodynamicData), allocatable :: thermodynamic_data(:)
  contains
    procedure, public :: write_files => write_files_PhononDos
  end type
  
  interface PathFrequencies
    ! Constructors.
    module function new_PathFrequencies(path_fraction,qpoint,frequencies, &
       & dynamical_matrix) result(this) 
      real(dp),              intent(in) :: path_fraction
      type(RealVector),      intent(in) :: qpoint
      real(dp),              intent(in) :: frequencies(:)
      type(DynamicalMatrix), intent(in) :: dynamical_matrix
      type(PathFrequencies)             :: this
    end function
  end interface
  
  interface SampledQpoint
    impure elemental module function new_SampledQpoint(qpoint, &
       & no_soft_frequencies) result(this) 
      type(RealVector), intent(in) :: qpoint
      integer,          intent(in) :: no_soft_frequencies
      type(SampledQpoint)          :: this
    end function
  end interface
  
  interface PdosBin
    impure elemental module function new_PdosBin(min_frequency, &
       & max_frequency,occupation) result(this) 
      real(dp), intent(in) :: min_frequency
      real(dp), intent(in) :: max_frequency
      real(dp), intent(in) :: occupation
      type(PdosBin)        :: this
    end function
  end interface
  
  interface PhononDispersion
    module function new_PhononDispersion(path,frequencies) result(this) 
      type(QpointPath),      intent(in) :: path
      type(PathFrequencies), intent(in) :: frequencies(:)
      type(PhononDispersion)            :: this
    end function
  end interface
  
  interface PhononDos
    module function new_PhononDos(qpoints,pdos,thermodynamic_data) &
       & result(this) 
      type(SampledQpoint),     intent(in) :: qpoints(:)
      type(PdosBin),           intent(in) :: pdos(:)
      type(ThermodynamicData), intent(in) :: thermodynamic_data(:)
      type(PhononDos)                     :: this
    end function
  end interface
  
  interface PhononDispersion
    ! ----------------------------------------------------------------------
    ! Generates the phonon dispersion curve.
    ! ----------------------------------------------------------------------
    module function new_PhononDispersion_CartesianHessian(supercell, &
       & min_images,hessian,stress_hessian,stress_supercell,         &
       & stress_min_images,path_string,no_path_points,logfile) result(this) 
      type(StructureData),    intent(in)           :: supercell
      type(MinImages),        intent(in)           :: min_images(:,:)
      type(CartesianHessian), intent(in)           :: hessian
      type(StressHessian),    intent(in), optional :: stress_hessian
      type(StructureData),    intent(in), optional :: stress_supercell
      type(MinImages),        intent(in), optional :: stress_min_images(:,:)
      type(String),           intent(in)           :: path_string
      integer,                intent(in)           :: no_path_points
      type(OFile),            intent(inout)        :: logfile
      type(PhononDispersion)                       :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Various printouts of the phonon dispersion.
    ! ----------------------------------------------------------------------
    module function path_PhononDispersion(this) result(output) 
      class(PhononDispersion), intent(in) :: this
      type(String), allocatable           :: output(:)
    end function
  end interface
  
  interface
    module function frequencies_PhononDispersion(this) result(output) 
      class(PhononDispersion), intent(in) :: this
      type(String), allocatable           :: output(:)
    end function
  end interface
  
  interface
    module function json_PhononDispersion(this,seedname,structure) &
       & result(output) 
      class(PhononDispersion), intent(in) :: this
      type(String),            intent(in) :: seedname
      type(StructureData),     intent(in) :: structure
      type(String), allocatable           :: output(:)
    end function
  end interface
  
  interface PhononDos
    ! ----------------------------------------------------------------------
    ! Calculate the frequency density-of-states by random sampling of
    !    the Brillouin zone.
    ! ----------------------------------------------------------------------
    module function new_PhononDos_CartesianHessian(supercell,min_images, &
       & hessian,stress_hessian,stress_supercell,stress_min_images,      &
       & thermal_energies,min_frequency,no_dos_samples,logfile,          &
       & random_generator) result(this) 
      type(StructureData),    intent(in)           :: supercell
      type(MinImages),        intent(in)           :: min_images(:,:)
      type(CartesianHessian), intent(in)           :: hessian
      type(StressHessian),    intent(in), optional :: stress_hessian
      type(StructureData),    intent(in), optional :: stress_supercell
      type(MinImages),        intent(in), optional :: stress_min_images(:,:)
      real(dp),               intent(in)           :: thermal_energies(:)
      real(dp),               intent(in)           :: min_frequency
      integer,                intent(in)           :: no_dos_samples
      type(OFile),            intent(inout)        :: logfile
      type(RandomReal),       intent(in)           :: random_generator
      type(PhononDos)                              :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module function write_PathFrequencies(this) result(output) 
      class(PathFrequencies), intent(in) :: this
      type(String), allocatable          :: output(:)
    end function
  end interface
  
  interface
    module subroutine read_SampledQpoint(this,input) 
      class(SampledQpoint), intent(out) :: this
      type(String),         intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_SampledQpoint(this) result(output) 
      class(SampledQpoint), intent(in) :: this
      type(String)                     :: output
    end function
  end interface
  
  interface SampledQpoint
    module function new_SampledQpoint_String(input) result(this) 
      type(String), intent(in) :: input
      type(SampledQpoint)      :: this
    end function
  end interface
  
  interface
    module subroutine read_PdosBin(this,input) 
      class(PdosBin), intent(out) :: this
      type(String),   intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_PdosBin(this) result(output) 
      class(PdosBin), intent(in) :: this
      type(String)                         :: output
    end function
  end interface
  
  interface PdosBin
    module function new_PdosBin_String(input) result(this) 
      type(String), intent(in) :: input
      type(PdosBin)            :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Multiple-file I/O.
    ! ----------------------------------------------------------------------
    module subroutine write_files_PhononDispersion(this,directory,seedname, &
       & structure) 
      class(PhononDispersion), intent(in) :: this
      type(String),            intent(in) :: directory
      type(String),            intent(in) :: seedname
      type(StructureData),     intent(in) :: structure
    end subroutine
  end interface
  
  interface
    module subroutine write_files_PhononDos(this,directory) 
      class(PhononDos), intent(in) :: this
      type(String),     intent(in) :: directory
    end subroutine
  end interface
end module
