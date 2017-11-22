! ======================================================================
! The third stage of Caesar.
! Uses the forces calculated previously to calculate harmonic normal modes.
! ======================================================================
module calculate_normal_modes_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function calculate_normal_modes_keywords() result(keywords)
  use keyword_module
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [KeywordData::]
end function

function calculate_normal_modes_mode() result(output)
  use caesar_modes_module
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'calculate_normal_modes'
  output%description = 'Finds harmonic normal modes. Should be called &
     &after run_harmonic.'
  output%keywords = calculate_normal_modes_keywords()
  output%main_subroutine => calculate_normal_modes
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine calculate_normal_modes(arguments)
  use utils_module, only : mkdir
  use ifile_module
  use linear_algebra_module
  use setup_harmonic_module
  use structure_module
  use dft_output_file_module
  use lte_module
  use unique_directions_module
  use group_module
  use qpoints_module
  use dictionary_module
  use normal_mode_module
  use dynamical_matrix_module
  use force_constants_module
  use lift_degeneracies_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Working directory.
  type(String) :: wd
  
  ! Input arguments.
  type(Dictionary) :: setup_harmonic_arguments
  
  ! No. supercells file.
  type(IFile) :: no_supercells_file
  
  ! Setup data.
  integer             :: no_supercells
  type(String)        :: dft_code
  type(String)        :: seedname
  type(StructureData) :: structure
  type(StructureData) :: supercell
  
  ! Force constant data.
  type(UniqueDirections)             :: unique_directions
  type(RealVector),      allocatable :: forces(:,:)
  type(RealMatrix),      allocatable :: force_constants(:,:,:)
  
  ! q-point data.
  type(StructureData)           :: large_supercell
  type(QpointData), allocatable :: qpoints(:)
  
  ! Lte output data.
  type(DynamicalMatrixAndMode), allocatable :: qpoint_modes(:)
  type(DynamicalMatrixAndMode), allocatable :: supercell_modes(:)
  logical,                      allocatable :: modes_calculated(:)
  integer                                   :: mode
  integer                                   :: atom
  integer                                   :: gvector
  
  ! Normal modes and their symmetries.
  type(LiftDegeneraciesReturn)   :: lifted_degeneracies
  
  ! Temporary variables.
  integer      :: i,j,k,ialloc
  type(String) :: sdir,qdir
  
  ! --------------------------------------------------
  ! Read in arguments from user.
  ! --------------------------------------------------
  wd = arguments%value('working_directory')
  
  ! --------------------------------------------------
  ! Read in previous arguments.
  ! --------------------------------------------------
  setup_harmonic_arguments = setup_harmonic_keywords()
  call setup_harmonic_arguments%read_file(wd//'/setup_harmonic.used_settings')
  dft_code = setup_harmonic_arguments%value('dft_code')
  seedname = setup_harmonic_arguments%value('seedname')
  
  no_supercells_file = wd//'/no_supercells.dat'
  no_supercells = int(no_supercells_file%line(1))
  
  structure = read_structure_file(wd//'/structure.dat')
  
  large_supercell = read_structure_file(wd//'/large_supercell.dat')
  
  qpoints = read_qpoints_file(wd//'/qpoints.dat')
  
  ! --------------------------------------------------
  ! Loop over supercells, calculating dynamical matrices and normal modes.
  ! Transfer information to the relevant q-points.
  ! --------------------------------------------------
  allocate( qpoint_modes(size(qpoints)),     &
          & modes_calculated(size(qpoints)), &
          & stat=ialloc); call err(ialloc)
  modes_calculated = .false.
  do i=1,no_supercells
    sdir = wd//'/Supercell_'//i
    
    ! Read in supercell structure data.
    supercell = read_structure_file(sdir//'/structure.dat')
    
    ! Read in symmetry group and unique atoms.
    unique_directions = read_unique_directions_file( &
       & sdir//'/unique_directions.dat')
    
    ! Calculate force constants.
    forces = read_forces(supercell,unique_directions,sdir,dft_code, &
       & seedname)
    
    ! Mass-reduce forces.
    do j=1,size(unique_directions)
      atom = unique_directions%atoms(j)
      do k=1,supercell%no_atoms
        forces(k,j) = forces(k,j)                          &
                    & / sqrt( supercell%atoms(atom)%mass() &
                    &       * supercell%atoms(k)%mass())
      enddo
    enddo
    
    force_constants = construct_force_constants(forces,supercell, &
       & unique_directions)
    
    ! Run normal mode analysis to find all normal modes of the supercell.
    supercell_modes = evaluate_normal_modes(supercell, force_constants)
    
    deallocate(force_constants)
    
    ! Pick out the normal modes corresponding to desired q-points.
    do_j : do j=1,size(qpoints)
      if (qpoints(j)%supercell_matrix==supercell%supercell) then
        do gvector=1,supercell%sc_size
          if (qpoints(j)%supercell_gvector==supercell%gvectors(gvector)) then
            qpoint_modes(j) = supercell_modes(gvector)
            
            lifted_degeneracies = lift_degeneracies( &
                    & qpoint_modes(j)%complex_modes, &
                    & supercell)
            qpoint_modes(j)%complex_modes = lifted_degeneracies%complex_modes
            modes_calculated(j) = .true.
            cycle do_j
          endif
        enddo
        
        call print_line(CODE_ERROR//': No matching G-vector found.')
        call err()
      endif
    enddo do_j
  enddo
  
  ! --------------------------------------------------
  ! Calculate information at remaining q-points from symmetry.
  ! --------------------------------------------------
  do_i : do i=1,size(qpoints)
    if (.not. modes_calculated(i)) then
      do j=1,i-1
        do k=1,size(structure%symmetries)
          if (structure%symmetries(j)%rotation * qpoints(i)%scaled_qpoint &
                                            & == qpoints(j)%scaled_qpoint) then
            qpoint_modes(j) = rotate_modes( qpoint_modes(k),         &
                                          & structure%symmetries(j), &
                                          & structure,               &
                                          & qpoints(i))
            modes_calculated(j) = .true.
            cycle do_i
          endif
        enddo
      enddo
      
      call print_line(CODE_ERROR//': No rotationally equivalent q-point &
         &found.')
      call err()
    endif
  enddo do_i
  
  ! --------------------------------------------------
  ! Write out output.
  ! --------------------------------------------------
  do i=1,size(qpoints)
    qdir = wd//'/qpoint_'//j
    call mkdir(qdir)
    
    ! Write out dynamical matrix.
    call write_dynamical_matrix_file(      &
       & qpoint_modes(i)%dynamical_matrix, &
       & qdir//'/dynamical_matrix.dat')
    
    ! Write out normal modes.
    do mode=1,structure%no_modes
      call qpoint_modes(i)%complex_modes(mode)%write_file( &
         & qdir//'/mode_'//mode//'.dat')
    enddo
  enddo
end subroutine
end module
