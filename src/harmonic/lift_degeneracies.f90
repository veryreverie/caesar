module lift_degeneracies_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  use normal_mode_module
  implicit none
  
  type :: LiftDegeneraciesReturn
    type(NormalMode), allocatable :: normal_modes(:)
  end type
contains

function lift_degeneracies(normal_modes,supercell,atom_symmetry_group) &
   & result(output)
  use normal_mode_module
  use structure_module
  use linear_algebra_module
  use group_module
  implicit none
  
  type(NormalMode),    intent(in) :: normal_modes(:)
  type(StructureData), intent(in) :: supercell
  type(Group),         intent(in) :: atom_symmetry_group(:)
  type(LiftDegeneraciesReturn)    :: output
  
  type(RealMatrix), allocatable :: rotations(:)
  
  logical, allocatable :: translations(:)
  type(Group) :: symmetry_group
  
  integer :: i,j,ialloc
  
  ! Identify purely translational symmetries.
  allocate( translations(supercell%no_symmetries), &
          & stat=ialloc); call err(ialloc)
  translations = .false.
  do i=1,supercell%no_symmetries
    if (supercell%rotations(i)==identity(3)) then
      do j=1,size(atom_symmetry_group(i))
        if (atom_symmetry_group(i)*j/=j) then
          translations(i) = .true.
          exit
        endif
      enddo
    endif
  enddo
  
  output%normal_modes = normal_modes
  
  rotations = calculate_cartesian_rotations(supercell)
  
  call print_line('')
  call print_line('R-vector translations:')
  do i=1,size(rotations)
    if (translations(i)) then
      call print_line(rotations(i))
      call print_line('')
    endif
  enddo
  
  call print_line('')
  call print_line('Rotations:')
  do i=1,size(rotations)
    call print_line(rotations(i))
    call print_line('')
  enddo
end function
end module
