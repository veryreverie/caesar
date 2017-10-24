! ======================================================================
! The third stage of Caesar.
! Uses the forces calculated previously to calculate harmonic normal modes.
! ======================================================================
module calculate_harmonic_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function calculate_harmonic_keywords() result(keywords)
  use keyword_module
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [KeywordData::]
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine calculate_harmonic(arguments)
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
  type(QpointData), allocatable :: qpoints_ibz(:)
  
  ! Lte output data.
  type(LteReturn) :: lte_result
  integer         :: mode
  integer         :: atom
  integer         :: gvector
  
  ! Normal modes and their symmetries.
  type(LiftDegeneraciesReturn) :: normal_modes
  
  ! Temporary variables.
  integer                   :: i,j,k
  type(String)              :: sdir,qdir
  
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
  
  qpoints_ibz = read_qpoints_file(wd//'/qpoints_ibz.dat')
  
  ! --------------------------------------------------
  ! Loop over supercells
  ! --------------------------------------------------
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
        forces(k,j) = forces(k,j) &
                    & / sqrt(supercell%mass(atom)*supercell%mass(k))
      enddo
    enddo
    force_constants = construct_force_constants(forces,supercell, &
       & unique_directions)
    
    ! Run normal mode analysis to find all normal modes of the supercell.
    lte_result = evaluate_freqs_on_grid(supercell, force_constants)
    
    deallocate(force_constants)
    
    ! Pick out the normal modes corresponding to desired q-points.
    do j=1,size(qpoints_ibz)
      if (qpoints_ibz(j)%sc_id/=i) then
        cycle
      endif
      
      qdir = wd//'/qpoint_'//j
      
      gvector = qpoints_ibz(j)%gvector_id
      
      ! Write out dynamical matrix.
      call write_dynamical_matrix_file(            &
         & lte_result%dynamical_matrices(gvector), &
         & qdir//'/dynamical_matrix.dat')
      
      ! Lift degeneracies using symmetry operators.
      normal_modes = lift_degeneracies( lte_result%normal_modes(:,gvector), &
                                      & supercell)
      
      ! Write out normal modes.
      qdir = wd//'/qpoint_'//j
      call mkdir(qdir)
      do mode=1,structure%no_modes
        call write_normal_mode_file( normal_modes%normal_modes(mode), &
                                   & qdir//'/mode_'//mode//'.dat')
      enddo
    enddo
  enddo
end subroutine
end module
