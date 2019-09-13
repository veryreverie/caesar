! ======================================================================
! Reads and writes extended .xyz files.
! ======================================================================
module xyz_module
  use utils_module
  use structure_module
  implicit none
  
  private
  
  public :: make_input_filename_xyz
  public :: read_input_file_xyz
  public :: write_input_file_xyz
contains

function make_input_filename_xyz(seedname) result(output)
  implicit none
  
  type(String), intent(in) :: seedname
  type(String)             :: output
  
  output = seedname//'.xyz'
end function

function read_input_file_xyz(filename) result(output)
  implicit none
  
  type(String), intent(in)  :: filename
  type(BasicStructure)      :: output
  
  ! The xyz file.
  type(IFile) :: xyz_file
  
  ! Variables for parsing the second line.
  character(:), allocatable :: line
  
  type(String) :: key
  type(String) :: val
  
  logical :: lattice_found
  logical :: properties_found
  logical :: species_found
  logical :: positions_found
  logical :: masses_found
  
  real(dp)                  :: lattice(3,3)
  type(String), allocatable :: properties(:)
  integer                   :: species_column
  integer                   :: positions_column
  integer                   :: masses_column
  
  ! Variables for storing atomic information.
  type(String)                 :: species
  type(RealVector)             :: position
  real(dp)                     :: mass
  type(BasicAtom), allocatable :: atoms(:)
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  xyz_file = IFile(filename)
  
  ! Check the no_atoms line.
  if (size(xyz_file)/=2+int(xyz_file%line(1))) then
    call print_line(ERROR//': '//filename//' is inconsistent with its first &
       &line, no. atoms.')
    call err()
  endif
  
  ! Parse the comment line.
  line = char(xyz_file%line(2))
  lattice_found = .false.
  properties_found = .false.
  do
    ! Find the first '=' in the line.
    i = first([(line(k:k)=='=',k=1,len(line))], default=0)
    
    ! Everything before the '=' is the key.
    if (i==0 .or. i==len(line)) then
      call print_line(ERROR//': Unable to parse line 2 of '//filename//'.')
      call err()
    else
      key = lower_case(trim(str(line(:i-1))))
    endif
    
    ! The text after the '=' is the value.
    ! Check if this is quoted or terminated by a space.
    if (line(i+1:i+1)=='"') then
      ! If the value is quoted with ", find the closing quote.
      i = i+2
      j = first( [(line(k:k)=='"',k=1,len(line))], &
               & mask=[(k>i,k=1,len(line))],       &
               & default=0                         )
    elseif (line(i+1:i+1)=="'") then
      ! If the value is quoted with ', find the closing quote.
      i = i+2
      j = first( [(line(k:k)=="'",k=1,len(line))], &
               & mask=[(k>i,k=1,len(line))],       &
               & default=0                         )
    else
      ! If the value is not quoted, find the separating space.
      i = i+1
      j = first( [(line(k:k)==' ',k=1,len(line))], &
               & mask=[(k>i,k=1,len(line))],       &
               & default=0                         )
    endif
    
    if (j==0) then
      if (.not. lattice_found) then
        call print_line(ERROR//': Unable to find "Lattice" tag on line 2 of &
           &'//filename//'.')
        call err()
      elseif (.not. properties_found) then
        call print_line(ERROR//': Unable to find "Properties" tag on line 2 &
           &of '//filename//'.')
        call err()
      endif
    endif
    
    ! Parse the value.
    val = trim(str(line(i:j-1)))
    
    if (key=='lattice') then
      lattice_found = .true.
      lattice(1,:) = dble(tokens(val,1,3)) / ANGSTROM_PER_BOHR
      lattice(2,:) = dble(tokens(val,4,6)) / ANGSTROM_PER_BOHR
      lattice(3,:) = dble(tokens(val,7,9)) / ANGSTROM_PER_BOHR
    elseif (key=='properties') then
      properties_found = .true.
      properties = tokens(val,delimiter=':')
      
      ! Check the properties tag is formated as expected;
      !    a series of "name", "type", "no. columns" trios,
      !    separated by colons.
      if (modulo(size(properties),3)/=0) then
        call print_line(ERROR//': Unable to parse "Properties" tag in '// &
           & filename//'.')
        call err()
      endif
      
      ! Search for the column names "species", "pos" and "masses".
      species_found = .false.
      positions_found = .false.
      masses_found = .false.
      k = 0
      do j=1,size(properties)/3
        if (properties(3*j-2)=='species') then
          species_found = .true.
          if (properties(3*j-1)/='S') then
            call print_line(ERROR//': "species" type in "Properties" tag of &
               &'//filename//' is not "S".')
            call err()
          elseif (properties(3*j)/='1') then
            call print_line(ERROR//': "species" in "Properties" tag of '// &
               & filename//' should have exactly one column.')
            call err()
          endif
          species_column = k+1
        elseif (properties(3*j-2)=='pos') then
          positions_found = .true.
          if (properties(3*j-1)/='R') then
            call print_line(ERROR//': "pos" type in "Properties" tag of &
               &'//filename//' is not "R".')
            call err()
          elseif (properties(3*j)/='3') then
            call print_line(ERROR//': "pos" in "Properties" tag of '// &
               & filename//' should have exactly three columns.')
            call err()
          endif
          positions_column = k+1
        elseif (properties(3*j-2)=='masses') then
          masses_found = .true.
          if (properties(3*j-1)/='R') then
            call print_line(ERROR//': "masses" type in "Properties" tag of &
               &'//filename//' is not "R".')
            call err()
          elseif (properties(3*j)/='1') then
            call print_line(ERROR//': "masses" in "Properties" tag of '// &
               & filename//' should have exactly one column.')
            call err()
          endif
          masses_column = k+1
        endif
        
        if (properties(3*j-2)=='C') then
          k = k+int(properties(3*j))*2
        else
          k = k+int(properties(3*j))
        endif
      enddo
      
      if (.not. all([species_found,positions_found,masses_found])) then
        call print_line(ERROR//': Unable to find all of "species", "pos" and &
           &"masses" in "Properties" tag of '//filename//'.')
        call err()
      endif
    endif
    
    ! Stop parsing the second line if both the lattice and properties have
    !    been found.
    if (lattice_found .and. properties_found) then
      exit
    endif
    
    ! Remove the parsed key and value from the line,
    !    and move on to parsing the next key and value.
    line = line(j+1:)
  enddo
  
  ! Parse atom lines.
  allocate(atoms(size(xyz_file)-2), stat=ialloc); call err(ialloc)
  do i=3,size(xyz_file)
    species = token(xyz_file%line(i),species_column)
    position = vec(dble(tokens( xyz_file%line(i),  &
                              & positions_column,  &
                              & positions_column+2 ))) / ANGSTROM_PER_BOHR
    mass = dble(token(xyz_file%line(i),masses_column)) / AMU_PER_ME
    atoms(i-2) = BasicAtom(species,mass,position)
  enddo
  
  ! Construct the output.
  output = BasicStructure(mat(lattice), atoms)
end function

subroutine write_input_file_xyz(structure,input_filename,output_filename)
  implicit none
  
  type(BasicStructure), intent(in)           :: structure
  type(String),         intent(in), optional :: input_filename
  type(String),         intent(in)           :: output_filename
  
  type(String) :: lattice_string
  
  type(OFile) :: xyz_file
  
  integer :: i,j
  
  xyz_file = OFile(output_filename)
  
  lattice_string = ''
  do i=1,3
    do j=1,3
      lattice_string = lattice_string// &
                   & structure%lattice_matrix%element(i,j) * ANGSTROM_PER_BOHR
      if (i/=3 .or. j/=3) then
        lattice_string = lattice_string//' '
      endif
    enddo
  enddo
  
  call xyz_file%print_line(size(structure%atoms))
  call xyz_file%print_line('Lattice="'//lattice_string//'" &
     &Properties="species:S:1:pos:R:3:masses:R:1"')
  do i=1,size(structure%atoms)
    associate (atom => structure%atoms(i))
      call xyz_file%print_line(                                &
         & atom%species                                //' '// &
         & atom%cartesian_position * ANGSTROM_PER_BOHR //' '// &
         & atom%mass * AMU_PER_ME                              )
    end associate
  enddo
end subroutine
end module
