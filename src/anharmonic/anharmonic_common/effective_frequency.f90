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
    real(dp), allocatable :: anharmonic_energies(:)
    real(dp), allocatable :: anharmonic_forces(:)
    real(dp), allocatable :: harmonic_energies(:)
    real(dp), allocatable :: harmonic_forces(:)
    real(dp), allocatable :: effective_energies(:)
    real(dp), allocatable :: effective_forces(:)
    real(dp), allocatable :: sampled_energies(:)
    real(dp), allocatable :: sampled_forces(:)
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

function new_EffectiveFrequency(mode_id,harmonic_frequency,                   &
   & effective_frequency,displacements,anharmonic_energies,anharmonic_forces, &
   & harmonic_energies,harmonic_forces,effective_energies,effective_forces,   &
   & sampled_energies,sampled_forces) result(this)
  implicit none
  
  integer,  intent(in)           :: mode_id
  real(dp), intent(in)           :: harmonic_frequency
  real(dp), intent(in)           :: effective_frequency
  real(dp), intent(in)           :: displacements(:)
  real(dp), intent(in)           :: anharmonic_energies(:)
  real(dp), intent(in)           :: anharmonic_forces(:)
  real(dp), intent(in)           :: harmonic_energies(:)
  real(dp), intent(in)           :: harmonic_forces(:)
  real(dp), intent(in)           :: effective_energies(:)
  real(dp), intent(in)           :: effective_forces(:)
  real(dp), intent(in), optional :: sampled_energies(:)
  real(dp), intent(in), optional :: sampled_forces(:)
  type(EffectiveFrequency)       :: this
  
  this%mode_id             = mode_id
  this%harmonic_frequency  = harmonic_frequency
  this%effective_frequency = effective_frequency
  this%displacements       = displacements
  this%anharmonic_energies = anharmonic_energies
  this%anharmonic_forces   = anharmonic_forces
  this%harmonic_energies   = harmonic_energies
  this%harmonic_forces     = harmonic_forces
  this%effective_energies  = effective_energies
  this%effective_forces    = effective_forces
  
  if (present(sampled_energies)) then
    this%sampled_energies = sampled_energies
  endif
  if (present(sampled_forces)) then
    this%sampled_forces = sampled_forces
  endif
end function

function new_EffectiveFrequency_potential(displacements,mode,real_modes, &
   & potential,structure) result(this)
  implicit none
  
  real(dp),             intent(in) :: displacements(:)
  type(ComplexMode),    intent(in) :: mode
  type(RealMode),       intent(in) :: real_modes(:)
  class(PotentialData), intent(in) :: potential
  type(StructureData),  intent(in) :: structure
  type(EffectiveFrequency)         :: this
  
  ! Output variables.
  real(dp)              :: effective_spring_constant
  real(dp)              :: effective_scaling
  real(dp)              :: effective_frequency
  real(dp), allocatable :: anharmonic_energies(:)
  real(dp), allocatable :: anharmonic_forces(:)
  real(dp), allocatable :: harmonic_energies(:)
  real(dp), allocatable :: harmonic_forces(:)
  real(dp), allocatable :: effective_energies(:)
  real(dp), allocatable :: effective_forces(:)
  
  ! The energy of the structure in equilibrium (displacement=0).
  real(dp) :: relative_energy
  
  ! Displacement in real mode co-ordinates.
  type(RealMode)             :: real_mode
  type(RealModeDisplacement) :: mode_displacement
  
  ! Force in real mode co-ordinates.
  type(RealModeForce) :: anharmonic_force
  
  ! Temporary variables.
  integer :: i,ialloc
  
  real_mode = real_modes(first(real_modes%id==min(mode%id,mode%paired_id)))
  
  ! Calculate the energy at displacement=0.
  relative_energy = potential%energy(RealModeDisplacement( &
                             & [RealSingleDisplacement::]) )
  
  allocate( anharmonic_energies(size(displacements)), &
          & anharmonic_forces(size(displacements)),   &
          & harmonic_energies(size(displacements)),   &
          & harmonic_forces(size(displacements)),     &
          & stat=ialloc); call err(ialloc)
  do i=1,size(displacements)
    mode_displacement = RealModeDisplacement([real_mode],[displacements(i)])
    anharmonic_energies(i) = potential%energy(mode_displacement) &
                         & - relative_energy
    harmonic_energies(i) = 0.5_dp               &
                       & * mode%spring_constant &
                       & * displacements(i)     &
                       & * displacements(i)
    anharmonic_force = potential%force(mode_displacement)
    anharmonic_forces(i) = anharmonic_force%force(real_mode)
    harmonic_forces(i) = - mode%spring_constant &
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
  effective_spring_constant = mode%spring_constant                   &
                      & * sum(harmonic_energies*anharmonic_energies) &
                      & / sum(harmonic_energies*harmonic_energies)
  if (effective_spring_constant>=0) then
    effective_frequency =  sqrt( effective_spring_constant)
  else
    effective_frequency = -sqrt(-effective_spring_constant)
  endif
  allocate( effective_energies(size(displacements)), &
          & effective_forces(size(displacements)),   &
          & stat=ialloc); call err(ialloc)
  do i=1,size(displacements)
    effective_energies(i) = 0.5_dp                    &
                        & * effective_spring_constant &
                        & * displacements(i)          &
                        & * displacements(i)
    effective_forces(i) = - effective_spring_constant &
                      & * displacements(i)
  enddo
  
  this = EffectiveFrequency( mode%id,             &
                           & mode%frequency,      &
                           & effective_frequency, &
                           & displacements,       &
                           & anharmonic_energies, &
                           & anharmonic_forces,   &
                           & harmonic_energies,   &
                           & harmonic_forces,     &
                           & effective_energies,  &
                           & effective_forces     )
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
  real(dp), allocatable :: anharmonic_energies(:)
  real(dp), allocatable :: anharmonic_forces(:)
  real(dp), allocatable :: harmonic_energies(:)
  real(dp), allocatable :: harmonic_forces(:)
  real(dp), allocatable :: effective_energies(:)
  real(dp), allocatable :: effective_forces(:)
  real(dp), allocatable :: sampled_energies(:)
  real(dp), allocatable :: sampled_forces(:)
  
  type(String), allocatable :: line(:)
  
  integer :: i,ialloc
  
  select type(this); type is(EffectiveFrequency)
    line = split_line(input(1))
    mode_id = int(line(3))
    
    line = split_line(input(3))
    harmonic_frequency = dble(line(3))
    
    line = split_line(input(4))
    effective_frequency = dble(line(3))
    
    allocate( displacements(size(input)-4),       &
            & anharmonic_energies(size(input)-4), &
            & harmonic_energies(size(input)-4),   &
            & effective_energies(size(input)-4),  &
            & stat=ialloc); call err(ialloc)
    
    if (line(size(line)-1)=='Sampled') then
      allocate( sampled_energies(size(input)-4), &
              & sampled_forces(size(input)-4),   &
              & stat=ialloc); call err(ialloc)
    endif
    
    do i=1,size(input)-4
      line = split_line(input(i+4))
      displacements(i) = dble(line(1))
      anharmonic_energies(i) = dble(line(2))
      anharmonic_forces(i) = dble(line(3))
      harmonic_energies(i) = dble(line(4))
      harmonic_forces(i) = dble(line(5))
      effective_energies(i) = dble(line(6))
      effective_forces(i) = dble(line(7))
      if (allocated(sampled_energies)) then
        sampled_energies(i) = dble(line(8))
      endif
      if (allocated(sampled_forces)) then
        sampled_forces(i) = dble(line(9))
      endif
    enddo
    
    if (allocated(sampled_energies)) then
      this = EffectiveFrequency( mode_id,             &
                               & harmonic_frequency,  &
                               & effective_frequency, &
                               & displacements,       &
                               & anharmonic_energies, &
                               & anharmonic_forces,   &
                               & harmonic_energies,   &
                               & harmonic_forces,     &
                               & effective_energies,  &
                               & effective_forces,    &
                               & sampled_energies,    &
                               & sampled_forces       )
    else
      this = EffectiveFrequency( mode_id,             &
                               & harmonic_frequency,  &
                               & effective_frequency, &
                               & displacements,       &
                               & anharmonic_energies, &
                               & anharmonic_forces,   &
                               & harmonic_energies,   &
                               & harmonic_forces,     &
                               & effective_energies,  &
                               & effective_forces)
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
    if ( allocated(this%sampled_energies) .and. &
       & allocated(this%sampled_forces)         ) then
      output = [ output,                                                    &
               & str('Displacement | Anharmonic energy | Anharmonic force | &
               &Harmonic energy | Harmonic force | Effective energy | &
               &Effective force | Sampled energy | Sampled force')          ]
      do i=1,size(this%displacements)
        output = [ output,                             &
                 & this%displacements(i)       //' '// &
                 & this%anharmonic_energies(i) //' '// &
                 & this%anharmonic_forces(i)   //' '// &
                 & this%harmonic_energies(i)   //' '// &
                 & this%harmonic_forces(i)     //' '// &
                 & this%effective_energies(i)  //' '// &
                 & this%effective_forces(i)    //' '// &
                 & this%sampled_energies(i)    //' '// &
                 & this%sampled_forces(i)              ]
      enddo
    else
      output = [ output,                                                    &
               & str('Displacement | Anharmonic energy | Anharmonic force | &
               &Harmonic energy | Harmonic force | Effective energy | &
               &Effective force')                                           ]
      do i=1,size(this%displacements)
        output = [ output,                             &
                 & this%displacements(i)       //' '// &
                 & this%anharmonic_energies(i) //' '// &
                 & this%anharmonic_forces(i)   //' '// &
                 & this%harmonic_energies(i)   //' '// &
                 & this%harmonic_forces(i)     //' '// &
                 & this%effective_energies(i)  //' '// &
                 & this%effective_forces(i)            ]
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
