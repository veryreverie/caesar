module test_copy_harmonic_module
contains
subroutine test_copy_harmonic(wd,cwd)
  use constants, only : directions
  use utils,     only : format_path, make_dft_output_filename
  use string_module
  use file_module
  use structure_module
  use unique_directions_module
  use dft_output_file_module
  use group_module
  use atom_mapping_module
  implicit none
  
  ! Current working directory.
  type(String), intent(in) :: wd
  type(String), intent(in) :: cwd
  
  ! Terminal inputs.
  type(String) :: copy_dir
  
  ! Previous settings.
  type(String), allocatable :: user_inputs_file(:)
  type(String)              :: dft_code
  type(String)              :: seedname
  
  type(String), allocatable :: no_supercells_file(:)
  integer                   :: no_supercells
  
  ! Directory and file names.
  type(String) :: sdir
  type(String) :: copy_dirs(2)
  type(String) :: new_dirs(2)
  type(String) :: dft_output_filename
  
  ! Structure information.
  type(StructureData)   :: structure_copy
  type(StructureData)   :: structure_new
  type(Group)           :: new_to_copy
  
  ! Unique direction information.
  type(UniqueDirections) :: unique_directions
  integer                :: atom
  character(1)           :: direction
  
  ! DFT output files.
  type(DftOutputFile) :: dft_output
  integer             :: copied_output
  
  ! Temporary variables.
  integer :: i,j,k,l,m
  
  ! ----------------------------------------------------------------------
  ! Read in directory to compare against and copy from.
  ! ----------------------------------------------------------------------
  call print_line('')
  call print_line('This test will check dft input files against a previous &
     &calculations,')
  call print_line('   and copy over dft output files.')
  call print_line('')
  call print_line('Where is the harmonic directory for comparison?')
  copy_dir = format_path(read_line_from_user(),cwd)
  
  ! ----------------------------------------------------------------------
  ! Read in previous settings.
  ! ----------------------------------------------------------------------
  user_inputs_file = read_lines(wd//'/user_input.txt')
  dft_code = user_inputs_file(1)
  seedname = user_inputs_file(2)
  
  no_supercells_file = read_lines(wd//'/no_sc.dat')
  no_supercells = int(no_supercells_file(1))
  
  ! ----------------------------------------------------------------------
  ! Loop over supercells.
  ! ----------------------------------------------------------------------
  do i=1,no_supercells
    call print_line('')
    call print_line('Supercell '//i//':')
    sdir = wd//'/Supercell_'//i
    
    ! Read in both structure files.
    structure_copy = read_structure_file(copy_dir//'/'//sdir//'/structure.dat')
    structure_new = read_structure_file(sdir//'/structure.dat')
    
    new_to_copy = atom_mapping(structure_new,structure_copy)
    
    call print_line('Atom positions correct.')
    call print_line('Atoms in new calculation map to those in &
       &old calculation as:')
    call print_line(new_to_copy%operation)
    
    unique_directions = &
       & read_unique_directions_file(sdir//'/unique_directions.dat')
    do j=1,size(unique_directions)
      atom = unique_directions%unique_atoms(j)
      
      do k=1,3
        if (k==2 .and. unique_directions%xy_symmetry(j)/=0) then
          cycle
        endif
        
        if (k==3 .and. ( unique_directions%yz_symmetry(j)/=0 .or. &
                       & unique_directions%xz_symmetry(j)/=0)) then
          cycle
        endif
        
        direction = directions(k)
        
        copy_dirs = &
           & (/copy_dir//'/'//sdir//'/atom.'//atom//'.disp.'//k//'/positive', &
           &   copy_dir//'/'//sdir//'/atom.'//atom//'.disp.'//k//'/negative' /)
        
        new_dirs = (/sdir//'/atom.'//atom//'.+d'//direction, &
                   & sdir//'/atom.'//atom//'.-d'//direction /)
        do l=1,2
          ! Read in the old castep file.
          dft_output_filename = make_dft_output_filename(dft_code,seedname)
          dft_output_filename = copy_dirs(l)//'/'//dft_output_filename
          dft_output = read_dft_output_file(dft_code,dft_output_filename)
        
          ! Make a fake castep output file.
          copied_output = open_write_file( &
             & new_dirs(l)//'/'//seedname//'.castep')
          call print_line( copied_output,&
                         & 'final energy, E = '//dft_output%energy)
          call print_line(copied_output,'*********************** forces')
          do m=1,5
            call print_line(copied_output,'')
          enddo
          do m=1,structure_copy%no_atoms
            call print_line(copied_output, &
               & ': '//structure_copy%species(operate(new_to_copy,m))//' : '// &
               & dft_output%forces(:,operate(new_to_copy,m)))
          enddo
          call print_line(copied_output,'')
          call print_line(copied_output,'*************************************&
                          &*****************')
          close(copied_output)
        enddo
      enddo
    enddo
  enddo
end subroutine
end module
