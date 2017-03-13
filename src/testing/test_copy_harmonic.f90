module test_copy_harmonic_module
contains
subroutine test_copy_harmonic()
  use constants, only : directions
  use string_module
  use file_module
  use structure_module
  use unique_directions_module
  use dft_output_file_module
  implicit none
  
  ! Parameters
  real(dp), parameter :: tol = 1.0e-10_dp
  
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
  
  ! Structure information.
  type(StructureData)   :: structure_copy
  type(StructureData)   :: structure_new
  real(dp), allocatable :: copy_frac_pos(:,:)
  real(dp), allocatable :: new_frac_pos(:,:)
  real(dp)              :: difference(3)
  integer, allocatable  :: new_to_copy(:)
  
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
  copy_dir = read_line_from_user()
  
  ! ----------------------------------------------------------------------
  ! Read in previous settings.
  ! ----------------------------------------------------------------------
  user_inputs_file = read_lines('user_input.txt')
  dft_code = user_inputs_file(1)
  seedname = user_inputs_file(2)
  
  no_supercells_file = read_lines('no_sc.dat')
  no_supercells = int(no_supercells_file(1))
  
  ! ----------------------------------------------------------------------
  ! Loop over supercells.
  ! ----------------------------------------------------------------------
  do i=1,no_supercells
    call print_line('')
    call print_line('Supercell '//i//':')
    sdir = 'Supercell_'//i
    
    ! Read in both structure files.
    structure_copy = read_structure_file(copy_dir//'/'//sdir//'/structure.dat')
    structure_new = read_structure_file(sdir//'/structure.dat')
    
    ! Check no_atoms is the same.
    if (structure_copy%no_atoms/=structure_new%no_atoms) then
      call print_line('Atom counts do not match')
      stop
    endif
    
    ! Check lattices are the same.
    if (.not. all(structure_copy%lattice-structure_new%lattice < tol)) then
      call print_line('Lattices do not match.')
    endif
    
    ! Calculate fractional atomic positions.
    allocate(copy_frac_pos(3,structure_copy%no_atoms))
    allocate(new_frac_pos(3,structure_copy%no_atoms))
    copy_frac_pos = matmul(structure_copy%recip_lattice,structure_copy%atoms)
    new_frac_pos = matmul(structure_new%recip_lattice,structure_new%atoms)
    
    ! Find mapping between copy and new atoms.
    allocate(new_to_copy(structure_copy%no_atoms))
    new_to_copy = 0
    do j=1,structure_new%no_atoms
      do k=1,structure_copy%no_atoms
        ! Check if new atom j and copy atom k are the same, down to translation
        !    by lattice vectors.
        difference = new_frac_pos(:,j) - copy_frac_pos(:,k)
        if (all(abs(difference-nint(difference)) < tol)) then
          if (new_to_copy(j)/=0) then
            call print_line('Duplicate atom: atom '//j//' in new matches &
               &atom '//new_to_copy(j)//' and atom '//k//' in copy.')
            stop
          endif
          new_to_copy(j) = k
        endif
      enddo
    enddo
    
    do j=1,structure_copy%no_atoms
      if (new_to_copy(j)==0) then
        call print_line('Atom '//j//' in new has no equivalent in copy.')
        stop
      endif
      
      do k=1,j-1
        if (new_to_copy(j)==new_to_copy(k)) then
          call print_line('Duplicate atom: atom '//new_to_copy(j)//' in copy &
             &matches atom '//j//' and atom '//k//' in new.')
        endif
      enddo
    enddo
    
    call print_line('Atom positions correct.')
    call print_line('Atoms in new (1 2...) match to atoms in copy:')
    call print_line(new_to_copy)
      
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
          dft_output = read_dft_output_file(dft_code,copy_dirs(l),seedname)
        
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
               & ': '//structure_copy%species(new_to_copy(m))//' : '// &
               & dft_output%forces(:,new_to_copy(m)))
          enddo
          call print_line(copied_output,'')
          call print_line(copied_output,'*************************************&
                          &*****************')
          close(copied_output)
        enddo
      enddo
    enddo
    deallocate(copy_frac_pos)
    deallocate(new_frac_pos)
    deallocate(new_to_copy)
  enddo
end subroutine
end module
