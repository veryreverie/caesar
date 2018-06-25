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
    integer               :: mode_id
    real(dp)              :: harmonic_frequency
    real(dp)              :: effective_frequency
    real(dp), allocatable :: displacements(:)
    real(dp), allocatable :: energies(:)
    real(dp), allocatable :: harmonic_energies(:)
    real(dp), allocatable :: effective_energies(:)
    real(dp), allocatable :: sampled_energies(:)
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
   & effective_energies,sampled_energies) result(this)
  implicit none
  
  integer,  intent(in)           :: mode_id
  real(dp), intent(in)           :: harmonic_frequency
  real(dp), intent(in)           :: effective_frequency
  real(dp), intent(in)           :: displacements(:)
  real(dp), intent(in)           :: energies(:)
  real(dp), intent(in)           :: harmonic_energies(:)
  real(dp), intent(in)           :: effective_energies(:)
  real(dp), intent(in), optional :: sampled_energies(:)
  type(EffectiveFrequency)       :: this
  
  this%mode_id             = mode_id
  this%harmonic_frequency  = harmonic_frequency
  this%effective_frequency = effective_frequency
  this%displacements       = displacements
  this%energies            = energies
  this%harmonic_energies   = harmonic_energies
  this%effective_energies  = effective_energies
  if (present(sampled_energies)) then
    this%sampled_energies = sampled_energies
  endif
end function

function new_EffectiveFrequency_potential(displacements,mode,real_modes, &
   & potential) result(this)
  implicit none
  
  real(dp),             intent(in) :: displacements(:)
  type(ComplexMode),    intent(in) :: mode
  type(RealMode),       intent(in) :: real_modes(:)
  class(PotentialData), intent(in) :: potential
  type(EffectiveFrequency)         :: this
  
  ! Output variables.
  integer               :: mode_id
  real(dp)              :: harmonic_frequency
  real(dp)              :: harmonic_coefficient
  real(dp)              :: effective_coefficient
  real(dp)              :: effective_frequency
  real(dp), allocatable :: energies(:)
  real(dp), allocatable :: harmonic_energies(:)
  real(dp), allocatable :: effective_energies(:)
  
  real(dp) :: relative_energy
  
  ! Displacement in real mode co-ordinates.
  type(RealMode)             :: real_mode
  type(RealModeDisplacement) :: mode_displacement
  
  ! Temporary variables.
  integer :: i,ialloc
  
  mode_id = mode%id
  harmonic_frequency = mode%frequency
  
  if (harmonic_frequency>=0) then
    harmonic_coefficient =  0.5_dp * harmonic_frequency * harmonic_frequency
  else
    harmonic_coefficient = -0.5_dp * harmonic_frequency * harmonic_frequency
  endif
  
  relative_energy = potential%energy(RealModeDisplacement( &
     & [RealSingleModeVector::]))
  
  real_mode = real_modes(first(real_modes%id==min(mode%id,mode%paired_id)))
  
  allocate( energies(size(displacements)),           &
          & harmonic_energies(size(displacements)),  &
          & stat=ialloc); call err(ialloc)
  do i=1,size(displacements)
    mode_displacement = RealModeDisplacement([real_mode],[displacements(i)])
    energies(i) = potential%energy(mode_displacement) - relative_energy
    harmonic_energies(i) = harmonic_coefficient &
                       & * displacements(i)     &
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
  effective_coefficient = harmonic_coefficient            &
                      & * sum(harmonic_energies*energies) &
                      & / sum(harmonic_energies*harmonic_energies)
  if (effective_coefficient>=0) then
    effective_frequency =  sqrt( 2*effective_coefficient)
  else
    effective_frequency = -sqrt(-2*effective_coefficient)
  endif
  allocate( effective_energies(size(displacements)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(displacements)
    effective_energies(i) = 0.5_dp                &
                        & * effective_coefficient &
                        & * displacements(i)      &
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
  
  integer               :: mode_id
  real(dp)              :: harmonic_frequency
  real(dp)              :: effective_frequency
  real(dp), allocatable :: displacements(:)
  real(dp), allocatable :: energies(:)
  real(dp), allocatable :: harmonic_energies(:)
  real(dp), allocatable :: effective_energies(:)
  real(dp), allocatable :: sampled_energies(:)
  
  type(String), allocatable :: line(:)
  
  integer :: i,ialloc
  
  select type(this); type is(EffectiveFrequency)
    line = split_line(input(1))
    mode_id = int(line(3))
    
    line = split_line(input(3))
    harmonic_frequency = dble(line(3))
    
    line = split_line(input(4))
    effective_frequency = dble(line(3))
    
    allocate( displacements(size(input)-4),      &
            & energies(size(input)-4),           &
            & harmonic_energies(size(input)-4),  &
            & effective_energies(size(input)-4), &
            & stat=ialloc); call err(ialloc)
    
    if (size(line)==13) then
      allocate(sampled_energies(size(input)-4), stat=ialloc); call err(ialloc)
    endif
    
    do i=1,size(input)-4
      line = split_line(input(i+4))
      displacements(i) = dble(line(1))
      energies(i) = dble(line(2))
      harmonic_energies(i) = dble(line(3))
      effective_energies(i) = dble(line(4))
      if (allocated(sampled_energies)) then
        sampled_energies(i) = dble(line(5))
      endif
    enddo
    
    if (allocated(sampled_energies)) then
      this = EffectiveFrequency( mode_id,             &
                               & harmonic_frequency,  &
                               & effective_frequency, &
                               & displacements,       &
                               & energies,            &
                               & harmonic_energies,   &
                               & effective_energies,  &
                               & sampled_energies)
    else
      this = EffectiveFrequency( mode_id,             &
                               & harmonic_frequency,  &
                               & effective_frequency, &
                               & displacements,       &
                               & energies,            &
                               & harmonic_energies,   &
                               & effective_energies)
    endif
  end select
end subroutine

function write_EffectiveFrequency(this) result(output)
  implicit none
  
  class(EffectiveFrequency), intent(in) :: this
  type(String), allocatable             :: output(:)
  
  integer :: i
  
  select type(this); type is(EffectiveFrequency)
    output = [ 'Mode ID: '//this%mode_id,                        &
             & 'Harmonic frequency: '//this%harmonic_frequency,  &
             & 'Effective frequency: '//this%effective_frequency ]
    if (allocated(this%sampled_energies)) then
      output = [ output,                                                   &
               & str('Displacement | Anharmonic energy | Harmonic energy | &
               &Effective energy | Sampled energy') ]
      do i=1,size(this%displacements)
        output = [ output,                            &
                 & this%displacements(i)      //' '// &
                 & this%energies(i)           //' '// &
                 & this%harmonic_energies(i)  //' '// &
                 & this%effective_energies(i) //' '// &
                 & this%sampled_energies(i)           ]
      enddo
    else
      output = [ output,                                                   &
               & str('Displacement | Anharmonic energy | Harmonic energy | &
               &Effective energy') ]
      do i=1,size(this%displacements)
        output = [ output,                           &
                 & this%displacements(i)     //' '// &
                 & this%energies(i)          //' '// &
                 & this%harmonic_energies(i) //' '// &
                 & this%effective_energies(i)        ]
      enddo
    endif
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
