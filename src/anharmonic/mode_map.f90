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
    real(dp), allocatable :: anharmonic_energies(:)
    real(dp), allocatable :: anharmonic_forces(:)
    real(dp), allocatable :: sampled_energies(:)
    real(dp), allocatable :: sampled_forces(:)
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

function new_ModeMap(mode_id,harmonic_frequency,displacements,                &
   & harmonic_energies,harmonic_forces,anharmonic_energies,anharmonic_forces, &
   & sampled_energies,sampled_forces) result(this)
  implicit none
  
  integer,  intent(in)           :: mode_id
  real(dp), intent(in)           :: harmonic_frequency
  real(dp), intent(in)           :: displacements(:)
  real(dp), intent(in)           :: harmonic_energies(:)
  real(dp), intent(in)           :: harmonic_forces(:)
  real(dp), intent(in)           :: anharmonic_energies(:)
  real(dp), intent(in)           :: anharmonic_forces(:)
  real(dp), intent(in), optional :: sampled_energies(:)
  real(dp), intent(in), optional :: sampled_forces(:)
  type(ModeMap)                  :: this
  
  this%mode_id             = mode_id
  this%harmonic_frequency  = harmonic_frequency
  this%displacements       = displacements
  this%harmonic_energies   = harmonic_energies
  this%harmonic_forces     = harmonic_forces
  this%anharmonic_energies = anharmonic_energies
  this%anharmonic_forces   = anharmonic_forces
  
  if (present(sampled_energies)) then
    this%sampled_energies = sampled_energies
  endif
  if (present(sampled_forces)) then
    this%sampled_forces = sampled_forces
  endif
end function

function new_ModeMap_potential(displacements,mode,potential,anharmonic_data) &
   & result(this)
  implicit none
  
  real(dp),             intent(in) :: displacements(:)
  type(RealMode),       intent(in) :: mode
  class(PotentialData), intent(in) :: potential
  type(AnharmonicData), intent(in) :: anharmonic_data
  type(ModeMap)                    :: this
  
  ! Output variables.
  real(dp), allocatable :: harmonic_energies(:)
  real(dp), allocatable :: harmonic_forces(:)
  real(dp), allocatable :: anharmonic_energies(:)
  real(dp), allocatable :: anharmonic_forces(:)
  
  ! Force in real mode co-ordinates.
  type(RealModeDisplacement) :: displacement
  type(RealModeForce)        :: anharmonic_force
  
  ! Temporary variables.
  integer :: i,ialloc
  
  allocate( harmonic_energies(size(displacements)),   &
          & harmonic_forces(size(displacements)),     &
          & anharmonic_energies(size(displacements)), &
          & anharmonic_forces(size(displacements)),   &
          & stat=ialloc); call err(ialloc)
  do i=1,size(displacements)
    displacement = RealModeDisplacement([mode],[displacements(i)])
    
    harmonic_energies(i) = 0.5_dp               &
                       & * mode%spring_constant &
                       & * displacements(i)     &
                       & * displacements(i)
    harmonic_forces(i) = - mode%spring_constant &
                     & * displacements(i)
    
    displacement = RealModeDisplacement([mode],[displacements(i)])
    anharmonic_energies(i) = ( potential%energy(displacement)   &
                         &   - potential%undisplaced_energy() ) &
                         & / anharmonic_data%anharmonic_supercell%sc_size
    anharmonic_force = potential%force(displacement)
    anharmonic_forces(i) = anharmonic_force%force(mode) &
                       & / anharmonic_data%anharmonic_supercell%sc_size
  enddo
  
  this = ModeMap( mode%id,             &
                & mode%frequency,      &
                & displacements,       &
                & harmonic_energies,   &
                & harmonic_forces,     &
                & anharmonic_energies, &
                & anharmonic_forces    )
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
  real(dp), allocatable :: anharmonic_energies(:)
  real(dp), allocatable :: anharmonic_forces(:)
  real(dp), allocatable :: sampled_energies(:)
  real(dp), allocatable :: sampled_forces(:)
  
  type(String), allocatable :: line(:)
  
  integer :: i,ialloc
  
  select type(this); type is(ModeMap)
    line = split_line(input(1))
    mode_id = int(line(3))
    
    line = split_line(input(2))
    harmonic_frequency = dble(line(3))
    
    allocate( displacements(size(input)-3),       &
            & harmonic_energies(size(input)-3),   &
            & harmonic_forces(size(input)-3),     &
            & anharmonic_energies(size(input)-3), &
            & anharmonic_forces(size(input)-3),   &
            & stat=ialloc); call err(ialloc)
    
    line = split_line(input(3))
    if (line(size(line)-2)=='Sampled') then
      allocate( sampled_energies(size(input)-3), &
              & sampled_forces(size(input)-3),   &
              & stat=ialloc); call err(ialloc)
    endif
    
    do i=1,size(input)-4
      line = split_line(input(i+4))
      displacements(i) = dble(line(1))
      harmonic_energies(i) = dble(line(2))
      harmonic_forces(i) = dble(line(3))
      anharmonic_energies(i) = dble(line(4))
      anharmonic_forces(i) = dble(line(5))
      if (allocated(sampled_energies)) then
        sampled_energies(i) = dble(line(6))
        sampled_forces(i) = dble(line(7))
      endif
    enddo
    
    if (allocated(sampled_energies)) then
      this = ModeMap( mode_id,             &
                    & harmonic_frequency,  &
                    & displacements,       &
                    & harmonic_energies,   &
                    & harmonic_forces,     &
                    & anharmonic_energies, &
                    & anharmonic_forces,   &
                    & sampled_energies,    &
                    & sampled_forces       )
    else
      this = ModeMap( mode_id,             &
                    & harmonic_frequency,  &
                    & displacements,       &
                    & harmonic_energies,   &
                    & harmonic_forces,     &
                    & anharmonic_energies, &
                    & anharmonic_forces    )
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
    if (allocated(this%sampled_energies)) then
      output = [ output,                                                &
               & str('Displacement | Harmonic energy | Harmonic force | &
               &Anharmonic energy | Anharmonic force | Sampled energy | &
               &Sampled force')                                         ]
      do i=1,size(this%displacements)
        output = [ output,                             &
                 & this%displacements(i)       //' '// &
                 & this%harmonic_energies(i)   //' '// &
                 & this%harmonic_forces(i)     //' '// &
                 & this%anharmonic_energies(i) //' '// &
                 & this%anharmonic_forces(i)   //' '// &
                 & this%sampled_energies(i)    //' '// &
                 & this%sampled_forces(i)              ]
      enddo
    else
      output = [ output,                                                &
               & str('Displacement | Harmonic energy | Harmonic force | &
               &Anharmonic energy | Anharmonic force')                  ]
      do i=1,size(this%displacements)
        output = [ output,                             &
                 & this%displacements(i)       //' '// &
                 & this%harmonic_energies(i)   //' '// &
                 & this%harmonic_forces(i)     //' '// &
                 & this%anharmonic_energies(i) //' '// &
                 & this%anharmonic_forces(i)           ]
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
