! ======================================================================
! A storage class, for all the variables which are passed to the
!    various potential methods.
! ======================================================================
module anharmonic_data_module
  use common_module
  
  use degenerate_symmetry_module
  use subspace_coupling_module
  implicit none
  
  private
  
  public :: AnharmonicData
  
  type, extends(Stringsable) :: AnharmonicData
    type(StructureData)                   :: structure
    type(StructureData)                   :: anharmonic_supercell
    type(QpointData),         allocatable :: qpoints(:)
    type(ComplexMode),        allocatable :: complex_modes(:)
    type(RealMode),           allocatable :: real_modes(:)
    type(DegenerateSubspace), allocatable :: degenerate_subspaces(:)
    type(DegenerateSymmetry), allocatable :: degenerate_symmetries(:)
    type(SubspaceCoupling),   allocatable :: subspace_couplings(:)
    logical                               :: vscf_basis_functions_only
    real(dp)                              :: maximum_weighted_displacement
    real(dp)                              :: frequency_of_max_displacement
  contains
    procedure, public :: read  => read_AnharmonicData
    procedure, public :: write => write_AnharmonicData
  end type
  
  interface AnharmonicData
    module procedure new_AnharmonicData
    module procedure new_AnharmonicData_Strings
    module procedure new_AnharmonicData_StringArray
  end interface
contains

function new_AnharmonicData(structure,anharmonic_supercell,qpoints,       &
   & complex_modes,real_modes,degenerate_subspaces,degenerate_symmetries, &
   & subspace_couplings,vscf_basis_functions_only,                        &
   & maximum_weighted_displacement,frequency_of_max_displacement) result(this)
  implicit none
  
  type(StructureData),      intent(in) :: structure
  type(StructureData),      intent(in) :: anharmonic_supercell
  type(QpointData),         intent(in) :: qpoints(:)
  type(ComplexMode),        intent(in) :: complex_modes(:)
  type(RealMode),           intent(in) :: real_modes(:)
  type(DegenerateSubspace), intent(in) :: degenerate_subspaces(:)
  type(DegenerateSymmetry), intent(in) :: degenerate_symmetries(:)
  type(SubspaceCoupling),   intent(in) :: subspace_couplings(:)
  logical,                  intent(in) :: vscf_basis_functions_only
  real(dp),                 intent(in) :: maximum_weighted_displacement
  real(dp),                 intent(in) :: frequency_of_max_displacement
  type(AnharmonicData)                 :: this
  
  this%structure                     = structure
  this%anharmonic_supercell          = anharmonic_supercell
  this%qpoints                       = qpoints
  this%complex_modes                 = complex_modes
  this%real_modes                    = real_modes
  this%degenerate_subspaces          = degenerate_subspaces
  this%degenerate_symmetries         = degenerate_symmetries
  this%subspace_couplings            = subspace_couplings
  this%vscf_basis_functions_only     = vscf_basis_functions_only
  this%maximum_weighted_displacement = maximum_weighted_displacement
  this%frequency_of_max_displacement = frequency_of_max_displacement
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_AnharmonicData(this,input)
  implicit none
  
  class(AnharmonicData), intent(out) :: this
  type(String),          intent(in)  :: input(:)
  
  type(StructureData)                   :: structure
  type(StructureData)                   :: anharmonic_supercell
  type(QpointData),         allocatable :: qpoints(:)
  type(ComplexMode),        allocatable :: complex_modes(:)
  type(RealMode),           allocatable :: real_modes(:)
  type(DegenerateSubspace), allocatable :: degenerate_subspaces(:)
  type(DegenerateSymmetry), allocatable :: degenerate_symmetries(:)
  type(SubspaceCoupling),   allocatable :: subspace_couplings(:)
  logical                               :: vscf_basis_functions_only
  real(dp)                              :: maximum_weighted_displacement
  real(dp)                              :: frequency_of_max_displacement
  
  type(StringArray), allocatable :: sections(:)
  
  type(String) :: separator
  
  integer :: i
  
  separator = repeat('=',70)
  
  select type(this); type is(AnharmonicData)
    ! Split the file into sections.
    sections = split_into_sections(input, separating_line=separator)
    
    ! Trim the header line from each section.
    do i=1,size(sections)
      sections(i)%strings = sections(i)%strings(2:)
    enddo
    
    ! Parse each section into the relevant variable.
    structure = StructureData(sections(1))
    anharmonic_supercell = StructureData(sections(2))
    qpoints = QpointData(split_into_sections(sections(3)))
    complex_modes = ComplexMode(split_into_sections(sections(4)))
    real_modes = RealMode(split_into_sections(sections(5)))
    degenerate_subspaces = DegenerateSubspace(split_into_sections(sections(6)))
    degenerate_symmetries = &
       & DegenerateSymmetry(split_into_sections(sections(7)))
    subspace_couplings = SubspaceCoupling(sections(8)%strings)
    vscf_basis_functions_only = lgcl(sections(9)%strings(1))
    maximum_weighted_displacement = dble(sections(10)%strings(1))
    frequency_of_max_displacement = dble(sections(11)%strings(1))
    
    ! Construct the output.
    this = AnharmonicData( structure,                     &
                         & anharmonic_supercell,          &
                         & qpoints,                       &
                         & complex_modes,                 &
                         & real_modes,                    &
                         & degenerate_subspaces,          &
                         & degenerate_symmetries,         &
                         & subspace_couplings,            &
                         & vscf_basis_functions_only,     &
                         & maximum_weighted_displacement, &
                         & frequency_of_max_displacement  )
  class default
    call err()
  end select
end subroutine

function write_AnharmonicData(this) result(output)
  implicit none
  
  class(AnharmonicData), intent(in) :: this
  type(String), allocatable         :: output(:)
  
  type(String) :: separator
  
  separator = repeat('=',70)
  
  select type(this); type is(AnharmonicData)
    output = [ str('Structure'),                                    &
             & str(this%structure),                                 &
             & separator,                                           &
             & str('Anharmonic Supercell'),                         &
             & str(this%anharmonic_supercell),                      &
             & separator,                                           &
             & str('q-points'),                                     &
             & str(this%qpoints, separating_line=''),               &
             & separator,                                           &
             & str('Complex modes'),                                &
             & str(this%complex_modes, separating_line=''),         &
             & separator,                                           &
             & str('Real modes'),                                   &
             & str(this%real_modes, separating_line=''),            &
             & separator,                                           &
             & str('Degenerate subspaces'),                         &
             & str(this%degenerate_subspaces, separating_line=''),  &
             & separator,                                           &
             & str('Degenerate symmetries'),                        &
             & str(this%degenerate_symmetries, separating_line=''), &
             & separator,                                           &
             & str('Subspace couplings'),                           &
             & str(this%subspace_couplings),                        &
             & separator,                                           &
             & str('VSCF basis functions only'),                    &
             & str(this%vscf_basis_functions_only),                 &
             & separator,                                           &
             & str('Maximum weighted displacement'),                &
             & str(this%maximum_weighted_displacement),             &
             & separator,                                           &
             & str('Frequency of maximum displacement'),            &
             & str(this%frequency_of_max_displacement)              ]
       
  class default
    call err()
  end select
end function

function new_AnharmonicData_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(AnharmonicData)     :: this
  
  call this%read(input)
end function

impure elemental function new_AnharmonicData_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(AnharmonicData)          :: this
  
  this = AnharmonicData(str(input))
end function
end module
