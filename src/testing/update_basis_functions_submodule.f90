submodule (caesar_update_basis_functions_module) caesar_update_basis_functions_submodule
  use caesar_testing_module
contains

module procedure update_basis_functions_mode
  integer :: ialloc
  
  output%mode_name = 'update_basis_functions'
  output%description = 'Updates the format of files to allow calculations &
     &performed using Caesar from December 2018 and earlier to be used with &
     &Caesar from January 2019 and later. If this update is run, neither &
     &setup_anharmonic nor run_anharmonic need be run again, but &
     &calculate_potential and later steps will need re-running.'
  allocate(output%keywords(0), stat=ialloc); call err(ialloc)
  output%main_subroutine => update_basis_functions_subroutine
end procedure

module procedure update_basis_functions_subroutine
  type(AnharmonicData) :: anharmonic_data
  type(IFile)          :: anharmonic_data_file
  
  type(String) :: coupling_dir
  
  type(IFile)                    :: basis_function_ifile
  type(StringArray), allocatable :: sections(:)
  type(OFile)                    :: basis_function_ofile
  type(String)                   :: new_line
  
  type(String), allocatable :: line(:)
  
  integer :: i
  
  ! Read in anharmonic data.
  anharmonic_data_file = IFile('anharmonic_data.dat')
  anharmonic_data = AnharmonicData(anharmonic_data_file%lines())
  
  do i=1,size(anharmonic_data%subspace_couplings)
    coupling_dir = 'sampling_points/coupling_'// &
       & left_pad(i, str(size(anharmonic_data%subspace_couplings)))
    
    new_line = 'Subspace Coupling: '//anharmonic_data%subspace_couplings(i)
    
    ! Read in old file.
    basis_function_ifile = IFile(coupling_dir//'/basis_functions.dat')
    sections = basis_function_ifile%sections()
    
    ! Check the file hasn't already been updated.
    line = split_line(lower_case(basis_function_ifile%line(1)))
    if (line(1)=='subspace') then
      call print_line(ERROR//': Basis function file already contains subspace &
         &coupling.')
      call quit()
    endif
    
    ! Remove the unique monomials section, and add a subspace coupling section.
    sections = [StringArray([new_line]), sections(:size(sections)-1)]
    
    ! Write out the new file.
    basis_function_ofile = OFile(coupling_dir//'/basis_functions.dat')
    call basis_function_ofile%print_lines(str(sections, separating_line=''))
  enddo
  
end procedure
end submodule
