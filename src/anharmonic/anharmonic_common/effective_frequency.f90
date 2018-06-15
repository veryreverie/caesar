! ======================================================================
! In order to generate the harmonic basis states needed to run VSCF,
!    the anharmonic potential is first fit to an effective harmonic
!    potential.
! ======================================================================
module effective_frequency_module
  use common_module
  
  use potential_module
  implicit none
  
  private
  
  public :: EffectiveFrequency
  
  type, extends(Stringsable) :: EffectiveFrequency
    integer                  :: mode_id
    real(dp)                 :: harmonic_frequency
    real(dp)                 :: effective_frequency
    real(dp),    allocatable :: displacements(:)
    complex(dp), allocatable :: energies(:)
    real(dp),    allocatable :: harmonic_energies(:)
    real(dp),    allocatable :: effective_energies(:)
  contains
    procedure, public :: read  => read_EffectiveFrequency
    procedure, public :: write => write_EffectiveFrequency
  end type
  
  interface EffectiveFrequency
    module procedure new_EffectiveFrequency
    module procedure new_EffectiveFrequency_potential
    module procedure new_EffectiveFrequency_StringArray
  end interface
contains

function new_EffectiveFrequency(mode_id,harmonic_frequency,        &
   & effective_frequency,displacements,energies,harmonic_energies, &
   & effective_energies) result(this)
  implicit none
  
  integer,     intent(in)  :: mode_id
  real(dp),    intent(in)  :: harmonic_frequency
  real(dp),    intent(in)  :: effective_frequency
  real(dp),    intent(in)  :: displacements(:)
  complex(dp), intent(in)  :: energies(:)
  real(dp),    intent(in)  :: harmonic_energies(:)
  real(dp),    intent(in)  :: effective_energies(:)
  type(EffectiveFrequency) :: this
  
  this%mode_id             = mode_id
  this%harmonic_frequency  = harmonic_frequency
  this%effective_frequency = effective_frequency
  this%displacements       = displacements
  this%energies            = energies
  this%harmonic_energies   = harmonic_energies
  this%effective_energies  = effective_energies
end function

function new_EffectiveFrequency_potential(displacements,mode,potential) &
   & result(this)
  implicit none
  
  real(dp),             intent(in) :: displacements(:)
  type(ComplexMode),    intent(in) :: mode
  class(PotentialData), intent(in) :: potential
  type(EffectiveFrequency)         :: this
  
  ! Output variables.
  integer                  :: mode_id
  real(dp)                 :: harmonic_frequency
  real(dp)                 :: effective_frequency
  complex(dp), allocatable :: energies(:)
  real(dp),    allocatable :: harmonic_energies(:)
  real(dp),    allocatable :: effective_energies(:)
  
  ! Displacement in complex mode co-ordinates.
  type(ComplexSingleModeVector) :: single_mode_displacement
  type(ComplexModeDisplacement) :: mode_displacement
  
  ! Temporary variables.
  integer :: i,ialloc
  
  mode_id = mode%id
  harmonic_frequency = mode%frequency
  
  allocate( energies(size(displacements)),           &
          & harmonic_energies(size(displacements)),  &
          & stat=ialloc); call err(ialloc)
  do i=1,size(displacements)
    single_mode_displacement = ComplexSingleModeVector( &
                                 & id=mode_id,          &
                                 & magnitude=cmplx(displacements(i),0.0_dp,dp))
    mode_displacement = ComplexModeDisplacement([single_mode_displacement])
    energies(i) = potential%energy(mode_displacement)
    harmonic_energies(i) = 0.5_dp             &
                       & * harmonic_frequency &
                       & * displacements(i)   &
                       & * displacements(i)
  enddo
  
  ! Fit effective harmonic potential, using 1-D linear least squares.
  ! The effective harmonic potential is x times the harmonic potential.
  ! If the harmonic and anharmonic potentials at displacement i are given by
  !    h_i and a_i respectively, then x is found by minimizing
  !    L = sum_i[ (x*h_i - a_i)^2 ]
  !      = sum_i[ x^2*{h_i}^2 - 2*x*h_i*a_i + {a_i}^2 ]
  ! -> 0 = dL/dx
  !      = sum_i[ 2*x*{h_i}^2 - 2*h_i*a_i ]
  ! -> x = sum_i[ h_i*a_i ] / sum_i[ {h_i}^2 ]
  effective_frequency = harmonic_frequency &
                    & * sum(harmonic_energies*energies) &
                    & / sum(harmonic_energies*harmonic_energies)
  allocate( effective_energies(size(displacements)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(displacements)
    effective_energies(i) = 0.5_dp              &
                        & * effective_frequency &
                        & * displacements(i)    &
                        & * displacements(i)
  enddo
  
  this = EffectiveFrequency( mode_id,             &
                           & harmonic_frequency,  &
                           & effective_frequency, &
                           & displacements,       &
                           & energies,            &
                           & harmonic_energies,   &
                           & effective_energies)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_EffectiveFrequency(this,input)
  implicit none
  
  class(EffectiveFrequency), intent(out) :: this
  type(String),              intent(in)  :: input(:)
  
  select type(this); type is(EffectiveFrequency)
    ! TODO
  end select
end subroutine

function write_EffectiveFrequency(this) result(output)
  implicit none
  
  class(EffectiveFrequency), intent(in) :: this
  type(String), allocatable             :: output(:)
  
  integer :: i
  
  select type(this); type is(EffectiveFrequency)
    output = [ 'Mode ID: '//this%mode_id,                            &
             & 'Harmonic frequency: '//this%harmonic_frequency,      &
             & 'Effective frequency: '//this%effective_frequency,    &
             & str('Displacement | Anharmonic energy | Harmonic energy | &
             &Effective energy') ]
    do i=1,size(this%energies)
      output = [ output,                           &
               & this%displacements(i)     //' '// &
               & this%energies(i)          //' '// &
               & this%harmonic_energies(i) //' '// &
               & this%effective_energies(i)]
    enddo
  end select
end function

impure elemental function new_EffectiveFrequency_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(EffectiveFrequency)      :: this
  
  this = input
end function
end module
