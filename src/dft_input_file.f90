module dft_input_file_module
contains

function castep_input_file_to_structure(filename) result(output)
  use constants, only : dp,pi,bohr,amu_to_me,kg_to_me,identity
  use string_module
  use structure_module
  use file_module
  use err_module
  implicit none
  
  type(String), intent(in) :: filename
  type(StructureData)      :: output
  
  type(String), allocatable :: cell_file(:)
  
  ! Line numbers.
  integer :: lattice_start_line
  integer :: lattice_end_line
  integer :: positions_start_line
  integer :: positions_end_line
  integer :: masses_start_line
  integer :: masses_end_line
  
  ! Lattice variables.
  logical  :: lattice_is_cart
  real(dp) :: conversion
  real(dp) :: lattice(3,3)
  real(dp) :: lengths(3)
  real(dp) :: angles(3)
  
  ! Atomic variables.
  logical                   :: positions_are_abs
  type(String), allocatable :: species(:)
  real(dp),     allocatable :: positions(:,:)
  integer                   :: no_atoms
  logical,      allocatable :: masses_found(:)
  
  ! Temporary variables.
  type(String), allocatable :: line(:)
  character(1)              :: first_char
  integer                   :: i,j
  
  cell_file = read_lines(filename)
  
  lattice_start_line = 0
  lattice_end_line = 0
  positions_start_line = 0
  positions_end_line = 0
  masses_start_line = 0
  masses_end_line = 0
  
  ! Work out line numbers.
  do i=1,size(cell_file)
    line = split(lower_case(cell_file(i)))
    if (size(line) >= 2) then
      if (line(1)=='%block' .and. line(2)=='lattice_cart') then
        lattice_start_line = i
        lattice_is_cart = .true.
      elseif (line(1)=='%block' .and. line(2)=='lattice_abc') then
        lattice_start_line = i
        lattice_is_cart = .false.
      elseif (line(1)=='%endblock' .and. ( line(2)=='lattice_cart' .or. &
                                         & line(2)=='lattice_abc')) then
        lattice_end_line = i
      elseif (line(1)=='%block' .and. line(2)=='positions_abs') then
        positions_start_line = i
        positions_are_abs = .true.
      elseif (line(1)=='%block' .and. line(2)=='positions_frac') then
        positions_start_line = i
        positions_are_abs = .false.
      elseif (line(1)=='%endblock' .and. ( line(2)=='positions_abs' .or. &
                                         & line(2)=='positions_frac')) then
        positions_end_line = i
      elseif (line(1)=='%block' .and. line(2)=='species_mass') then
        masses_start_line = i
      elseif (line(1)=='%endblock' .and. line(2)=='species_mass') then
        masses_end_line = i
      endif
    endif
  enddo
  
  ! Check all blocks exist and are of reasonable sizes.
  if (lattice_start_line==0) then
    call print_line('Error: lattice block not present in '//filename)
    call err()
  elseif (lattice_end_line-lattice_start_line <= 1) then
    call print_line('Error: lattice block of unexpected size in '//filename)
    call err()
  elseif (positions_start_line==0) then
    call print_line('Error: positions block not present in '//filename)
    call err()
  elseif (positions_end_line-positions_start_line <= 1) then
    call print_line('Error: positions block of unexpected size in '//filename)
    call err()
  elseif (masses_start_line==0) then
    call print_line('Error: mass block not present in '//filename)
    call err()
  elseif (masses_end_line-masses_start_line <= 1) then
    call print_line('Error: mass block of unexpected size in '//filename)
    call err()
  endif
  
  ! Parse lattice.
  conversion = bohr
  j=1
  
  do i=lattice_start_line+1,lattice_end_line-1
    line = split(lower_case(cell_file(i)))
    
    ! Ignore comments.
    first_char = char(line(1))
    if (first_char=='!') then
      cycle
    
    ! Read units if present.
    elseif (line(1) == 'bohr') then
      conversion = 1.0_dp
    elseif (line(1) == 'a0') then
      conversion = 1.0_dp
    elseif (line(1) == 'm') then
      conversion = 1e10_dp*bohr
    elseif (line(1) == 'cm') then
      conversion = 1e8_dp*bohr
    elseif (line(1) == 'nm') then
      conversion = 1e1_dp*bohr
    elseif (line(1) == 'ang') then
      conversion = bohr
      
    elseif (lattice_is_cart) then
      if (j>3) then
        call print_line('Error: too many lattice lines found in '//filename)
        call err()
      endif
      
      lattice(j,:) = dble(line(1:3))*conversion
      j = j+1
    
    else
      if (j==1) then
        lengths = dble(line(1:3)) * conversion
      elseif (j==2) then
        angles = dble(line(1:3)) * pi/180
      else
        call print_line('Error: too many lattice lines found in '//filename)
        call err()
      endif
      j = j+1
    endif
  enddo
  
  if (.not. lattice_is_cart) then
    lattice = 0.0_dp
    lattice(1,1) = lengths(1)
    lattice(2,1) = lengths(2)*dcos(angles(3))
    lattice(2,2) = lengths(2)*dsin(angles(3))
    lattice(3,1) = lengths(3)*dcos(angles(2))
    lattice(3,2) = lengths(3)* ( dcos(angles(1)) &
                           &   - dcos(angles(2))*dcos(angles(3))) &
                           & / dsin(angles(3))
    lattice(3,3) = dsqrt(lengths(3)**2 - lattice(3,1)**2 - lattice(3,2)**2)
  endif
  
  ! Parse atoms.
  conversion = bohr
  j=1
  allocate(species(positions_end_line-positions_start_line))
  allocate(positions(3,positions_end_line-positions_start_line-1))
  
  do i=positions_start_line+1,positions_end_line-1
    line = split(lower_case(cell_file(i)))
    
    ! Ignore comments.
    first_char = char(line(1))
    if (first_char=='!') then
      cycle
    
    ! Read units if present.
    elseif (line(1) == 'bohr' .or. line(1) == 'a0') then
      conversion = 1.0_dp
    elseif (line(1) == 'm') then
      conversion = 1e10_dp*bohr
    elseif (line(1) == 'cm') then
      conversion = 1e8_dp*bohr
    elseif (line(1) == 'nm') then
      conversion = 1e1_dp*bohr
    elseif (line(1) == 'ang') then
      conversion = bohr
    
    elseif (j>size(positions,2)) then
      call print_line('Error: too many atom lines found in '//filename)
    
    ! Read in atomic positions.
    else
      line = split(cell_file(i)) ! N.B. no lower_case
      species(j) = line(1)
      positions(:,j) = dble(line(2:4))
      j = j+1
    endif
  enddo
  
  no_atoms = j-1
  
  if (positions_are_abs) then
    positions(:,:no_atoms) = positions(:,:no_atoms) * conversion
  else
    positions(:,:no_atoms) = matmul(transpose(lattice),positions(:,:no_atoms))
  endif
  
  ! Make output.
  call new(output,j-1,0,1)
  
  output%lattice = lattice
  output%species = species(:no_atoms)
  output%atoms = positions(:,:no_atoms)
  
  conversion = amu_to_me
  allocate(masses_found(no_atoms))
  masses_found = .false.
  do i=masses_start_line+1,masses_end_line-1
    line = split(lower_case(cell_file(i)))
    
    ! Ignore comments.
    first_char = char(line(1))
    if (first_char=='!') then
      cycle
    
    ! Read units if present.
    elseif (line(1) == 'me') then
      conversion = 1.0_dp
    elseif (line(1) == 'amu') then
      conversion = amu_to_me
    elseif (line(1) == 'kg') then
      conversion = kg_to_me
    elseif (line(1) == 'g') then
      conversion = 1e-3_dp*kg_to_me
    
    ! Read in masses.
    else
      line = split(cell_file(i)) ! N.B. no lower_case
      do j=1,no_atoms
        if (line(1)==species(j)) then
          masses_found(j) = .true.
          output%mass(j) = dble(line(2))*conversion
        endif
      enddo
    endif
  enddo
  
  if (.not. all(masses_found)) then
    call print_line('Error: not all masses specified in '//filename)
    call err()
  endif
  
  output%supercell = identity
  output%gvectors(:,1) = 0
  
  call calculate_derived_atom_quantities(output)
  call calculate_derived_supercell_quantities(output)
end function

function dft_input_file_to_structure(dft_code,filename) result(output)
  use string_module
  use structure_module
  use err_module
  implicit none
  
  type(String), intent(in) :: dft_code
  type(String), intent(in) :: filename
  type(StructureData)      :: output
  
  if (dft_code == 'castep') then
    output = castep_input_file_to_structure(filename)
  else
    call print_line('Conversion from dft file to structure not yet supported.')
    call err()
  endif
end function
end module
