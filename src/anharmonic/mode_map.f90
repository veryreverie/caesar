! ======================================================================
! A map of the potential along a given mode.
! ======================================================================
module mode_map_module
  use common_module
  
  use anharmonic_common_module
  implicit none
  
  private
  
  public :: ModeMap
  
  type, extends(Stringsable) :: ModeMap
    integer               :: mode_id
    real(dp)              :: harmonic_frequency
    real(dp), allocatable :: displacements(:)
    real(dp), allocatable :: harmonic_energies(:)
    real(dp), allocatable :: harmonic_forces(:)
    real(dp), allocatable :: anharmonic_cos_energies(:)
    real(dp), allocatable :: anharmonic_cos_forces(:)
    real(dp), allocatable :: anharmonic_sin_energies(:)
    real(dp), allocatable :: anharmonic_sin_forces(:)
    real(dp), allocatable :: sampled_cos_energies(:)
    real(dp), allocatable :: sampled_cos_forces(:)
    real(dp), allocatable :: sampled_sin_energies(:)
    real(dp), allocatable :: sampled_sin_forces(:)
  contains
    procedure, public :: read  => read_ModeMap
    procedure, public :: write => write_ModeMap
  end type
  
  interface ModeMap
    module procedure new_ModeMap
    module procedure new_ModeMap_potential
    module procedure new_ModeMap_Strings
    module procedure new_ModeMap_StringArray
  end interface
contains

function new_ModeMap(mode_id,harmonic_frequency,displacements, &
   & harmonic_energies,harmonic_forces,anharmonic_cos_energies, &
   & anharmonic_cos_forces,anharmonic_sin_energies,anharmonic_sin_forces, &
   & sampled_cos_energies,sampled_cos_forces,sampled_sin_energies, &
   & sampled_sin_forces) result(this)
  implicit none
  
  integer,  intent(in)           :: mode_id
  real(dp), intent(in)           :: harmonic_frequency
  real(dp), intent(in)           :: displacements(:)
  real(dp), intent(in)           :: harmonic_energies(:)
  real(dp), intent(in)           :: harmonic_forces(:)
  real(dp), intent(in)           :: anharmonic_cos_energies(:)
  real(dp), intent(in)           :: anharmonic_cos_forces(:)
  real(dp), intent(in)           :: anharmonic_sin_energies(:)
  real(dp), intent(in)           :: anharmonic_sin_forces(:)
  real(dp), intent(in), optional :: sampled_cos_energies(:)
  real(dp), intent(in), optional :: sampled_cos_forces(:)
  real(dp), intent(in), optional :: sampled_sin_energies(:)
  real(dp), intent(in), optional :: sampled_sin_forces(:)
  type(ModeMap)                  :: this
  
  this%mode_id                 = mode_id
  this%harmonic_frequency      = harmonic_frequency
  this%displacements           = displacements
  this%harmonic_energies       = harmonic_energies
  this%harmonic_forces         = harmonic_forces
  this%anharmonic_cos_energies = anharmonic_cos_energies
  this%anharmonic_cos_forces   = anharmonic_cos_forces
  this%anharmonic_sin_energies = anharmonic_sin_energies
  this%anharmonic_sin_forces   = anharmonic_sin_forces
  
  if (present(sampled_cos_energies)) then
    this%sampled_cos_energies = sampled_cos_energies
  endif
  if (present(sampled_cos_forces)) then
    this%sampled_cos_forces = sampled_cos_forces
  endif
  if (present(sampled_sin_energies)) then
    this%sampled_sin_energies = sampled_sin_energies
  endif
  if (present(sampled_sin_forces)) then
    this%sampled_sin_forces = sampled_sin_forces
  endif
end function

function new_ModeMap_potential(displacements,mode,real_modes,potential) &
   & result(this)
  implicit none
  
  real(dp),             intent(in) :: displacements(:)
  type(ComplexMode),    intent(in) :: mode
  type(RealMode),       intent(in) :: real_modes(:)
  class(PotentialData), intent(in) :: potential
  type(ModeMap)                    :: this
  
  ! Output variables.
  real(dp), allocatable :: harmonic_energies(:)
  real(dp), allocatable :: harmonic_forces(:)
  real(dp), allocatable :: anharmonic_cos_energies(:)
  real(dp), allocatable :: anharmonic_cos_forces(:)
  real(dp), allocatable :: anharmonic_sin_energies(:)
  real(dp), allocatable :: anharmonic_sin_forces(:)
  
  ! Displacement in real mode co-ordinates.
  type(RealMode)             :: cos_mode
  type(RealMode)             :: sin_mode
  
  ! Force in real mode co-ordinates.
  type(RealModeDisplacement) :: displacement
  type(RealModeForce)        :: anharmonic_force
  
  ! Temporary variables.
  integer :: i,ialloc
  
  cos_mode = real_modes(first(real_modes%id==min(mode%id,mode%paired_id)))
  sin_mode = real_modes(first(real_modes%id==max(mode%id,mode%paired_id)))
  
  allocate( harmonic_energies(size(displacements)),       &
          & harmonic_forces(size(displacements)),         &
          & anharmonic_cos_energies(size(displacements)), &
          & anharmonic_cos_forces(size(displacements)),   &
          & anharmonic_sin_energies(size(displacements)), &
          & anharmonic_sin_forces(size(displacements)),   &
          & stat=ialloc); call err(ialloc)
  do i=1,size(displacements)
    displacement = RealModeDisplacement([cos_mode],[displacements(i)])
    
    harmonic_energies(i) = 0.5_dp               &
                       & * mode%spring_constant &
                       & * displacements(i)     &
                       & * displacements(i)
    harmonic_forces(i) = - mode%spring_constant &
                     & * displacements(i)
    
    displacement = RealModeDisplacement([cos_mode],[displacements(i)])
    anharmonic_cos_energies(i) = potential%energy(displacement) &
                               - potential%undisplaced_energy()
    anharmonic_force = potential%force(displacement)
    anharmonic_cos_forces(i) = anharmonic_force%force(cos_mode)
    
    displacement = RealModeDisplacement([sin_mode],[displacements(i)])
    anharmonic_sin_energies(i) = potential%energy(displacement) &
                               - potential%undisplaced_energy()
    anharmonic_force = potential%force(displacement)
    anharmonic_sin_forces(i) = anharmonic_force%force(sin_mode)
  enddo
  
  this = ModeMap( mode%id,                 &
                & mode%frequency,          &
                & displacements,           &
                & harmonic_energies,       &
                & harmonic_forces,         &
                & anharmonic_cos_energies, &
                & anharmonic_cos_forces,   &
                & anharmonic_sin_energies, &
                & anharmonic_sin_forces    )
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_ModeMap(this,input)
  implicit none
  
  class(ModeMap), intent(out) :: this
  type(String),   intent(in)  :: input(:)
  
  integer               :: mode_id
  real(dp)              :: harmonic_frequency
  real(dp), allocatable :: displacements(:)
  real(dp), allocatable :: harmonic_energies(:)
  real(dp), allocatable :: harmonic_forces(:)
  real(dp), allocatable :: anharmonic_cos_energies(:)
  real(dp), allocatable :: anharmonic_cos_forces(:)
  real(dp), allocatable :: anharmonic_sin_energies(:)
  real(dp), allocatable :: anharmonic_sin_forces(:)
  real(dp), allocatable :: sampled_cos_energies(:)
  real(dp), allocatable :: sampled_cos_forces(:)
  real(dp), allocatable :: sampled_sin_energies(:)
  real(dp), allocatable :: sampled_sin_forces(:)
  
  type(String), allocatable :: line(:)
  
  integer :: i,ialloc
  
  select type(this); type is(ModeMap)
    line = split_line(input(1))
    mode_id = int(line(3))
    
    line = split_line(input(2))
    harmonic_frequency = dble(line(3))
    
    allocate( displacements(size(input)-3),           &
            & harmonic_energies(size(input)-3),       &
            & harmonic_forces(size(input)-3),         &
            & anharmonic_cos_energies(size(input)-3), &
            & anharmonic_cos_forces(size(input)-3),   &
            & anharmonic_sin_energies(size(input)-3), &
            & anharmonic_sin_forces(size(input)-3),   &
            & stat=ialloc); call err(ialloc)
    
    line = split_line(input(3))
    if (line(size(line)-2)=='Sampled') then
      allocate( sampled_cos_energies(size(input)-3), &
              & sampled_cos_forces(size(input)-3),   &
              & sampled_sin_energies(size(input)-3), &
              & sampled_sin_forces(size(input)-3),   &
              & stat=ialloc); call err(ialloc)
    endif
    
    do i=1,size(input)-4
      line = split_line(input(i+4))
      displacements(i) = dble(line(1))
      harmonic_energies(i) = dble(line(2))
      harmonic_forces(i) = dble(line(3))
      anharmonic_cos_energies(i) = dble(line(4))
      anharmonic_cos_forces(i) = dble(line(5))
      anharmonic_sin_energies(i) = dble(line(6))
      anharmonic_sin_forces(i) = dble(line(7))
      if (allocated(sampled_cos_energies)) then
        sampled_cos_energies(i) = dble(line(8))
        sampled_cos_forces(i) = dble(line(9))
        sampled_sin_energies(i) = dble(line(10))
        sampled_sin_forces(i) = dble(line(11))
      endif
    enddo
    
    if (allocated(sampled_cos_energies)) then
      this = ModeMap( mode_id,                 &
                    & harmonic_frequency,      &
                    & displacements,           &
                    & harmonic_energies,       &
                    & harmonic_forces,         &
                    & anharmonic_cos_energies, &
                    & anharmonic_cos_forces,   &
                    & anharmonic_sin_energies, &
                    & anharmonic_sin_forces,   &
                    & sampled_cos_energies,    &
                    & sampled_cos_forces,      &
                    & sampled_sin_energies,    &
                    & sampled_sin_forces       )
    else
      this = ModeMap( mode_id,                 &
                    & harmonic_frequency,      &
                    & displacements,           &
                    & harmonic_energies,       &
                    & harmonic_forces,         &
                    & anharmonic_cos_energies, &
                    & anharmonic_cos_forces,   &
                    & anharmonic_sin_energies, &
                    & anharmonic_sin_forces    )
    endif
  class default
    call err()
  end select
end subroutine

function write_ModeMap(this) result(output)
  implicit none
  
  class(ModeMap), intent(in) :: this
  type(String), allocatable  :: output(:)
  
  integer :: i
  
  select type(this); type is(ModeMap)
    output = [ 'Mode ID: '//this%mode_id,                      &
             & 'Harmonic frequency: '//this%harmonic_frequency ]
    if (allocated(this%sampled_cos_energies)) then
      output = [ output,                                                    &
               & str('Displacement | Harmonic energy | Harmonic force | &
               &Anharmonic cos energy | Anharmonic cos force | Anharmonic sin &
               &energy | Anharmonic sin force | Sampled cos energy | Sampled &
               &cos force | Sampled sin energy | Sampled sin force')          ]
      do i=1,size(this%displacements)
        output = [ output,                                 &
                 & this%displacements(i)           //' '// &
                 & this%harmonic_energies(i)       //' '// &
                 & this%harmonic_forces(i)         //' '// &
                 & this%anharmonic_cos_energies(i) //' '// &
                 & this%anharmonic_cos_forces(i)   //' '// &
                 & this%anharmonic_sin_energies(i) //' '// &
                 & this%anharmonic_sin_forces(i)   //' '// &
                 & this%sampled_cos_energies(i)    //' '// &
                 & this%sampled_cos_forces(i)      //' '// &
                 & this%sampled_sin_energies(i)    //' '// &
                 & this%sampled_sin_forces(i)              ]
      enddo
    else
      output = [ output,                                                    &
               & str('Displacement | Harmonic energy | Harmonic force | &
               &Anharmonic cos energy | Anharmonic cos force | Anharmonic sin &
               &energy | Anharmonic sin force')                               ]
      do i=1,size(this%displacements)
        output = [ output,                                 &
                 & this%displacements(i)           //' '// &
                 & this%harmonic_energies(i)       //' '// &
                 & this%harmonic_forces(i)         //' '// &
                 & this%anharmonic_cos_energies(i) //' '// &
                 & this%anharmonic_cos_forces(i)   //' '// &
                 & this%anharmonic_sin_energies(i) //' '// &
                 & this%anharmonic_sin_forces(i)           ]
      enddo
    endif
  class default
    call err()
  end select
end function

function new_ModeMap_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(ModeMap)            :: this
  
  call this%read(input)
end function

impure elemental function new_ModeMap_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(ModeMap)                 :: this
  
  this = ModeMap(str(input))
end function
end module
