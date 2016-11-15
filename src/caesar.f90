program caesar
  use utils, only : command_line_args
  implicit none
  
  integer                        :: i             ! loop index
  character(len=32), allocatable :: args(:)       ! command line arguments
  character(len=:),  allocatable :: arg           ! command line argument
  integer                        :: return_status ! system() status

  ! read in command line arguments
  args = command_line_args()
  
  if (size(args) == 0) then
    write(*,*) "No arguments given. For help, call caesar -h"
    deallocate(args)
    stop
  else
    allocate(character(len=len(trim(args(1)))) :: arg)
    arg = trim(args(1))
  endif
  
  if (arg == "-h" .or. arg == "--help") then
    write(*,*) "caesar [-h] [option]"
    write(*,*) ""
    write(*,*) "-h :"
    write(*,*) "  Displays this help text"
    write(*,*) ""
    write(*,*) "option : utilities :"
    write(*,*) "  hartree_to_eV :"
    write(*,*) "    Provides a Hartree to eV calculator"
    write(*,*) ""
    write(*,*) "option : harmonic calculations :"
    write(*,*) "  setup_harmonic :"
    write(*,*) "    Sets up calculation"
    write(*,*) "  convert_harmonic :"
    write(*,*) "    Converts calculation to specific code"
    write(*,*) "    Choices are castep, vasp and quantum espresso"
    write(*,*) "    Should be called after setup_harmonic"
    write(*,*) "  tcm_cluster_run_harmonic :"
    write(*,*) "    Runs calculation on the TCM cluster"
    write(*,*) "    Should be called after convert_harmonic"
    write(*,*) "  rutgers_run_harmonic :"
    write(*,*) "    Runs calculations"
    write(*,*) "    Should be called after convert_harmonic"
    write(*,*) "  lte_harmonic :"
    write(*,*) "    Runs harmonic calculations"
    write(*,*) "    Should be run after one of the run_harmonic options"
    write(*,*) "  clear_all :"
    write(*,*) "    Deletes all temporary files and folders"
    write(*,*) ""
    write(*,*) "option : quadratic calculations :"
    write(*,*) "  [quadratic help text yet to be written]"
  ! Wrappers for Fortran 
  elseif (arg == "band_folding") then
    return_status = system(arg)
  elseif (arg == "calculate_anharmonic") then
    return_status = system(arg)
  elseif (arg == "calculate_bs") then
    return_status = system(arg)
  elseif (arg == "calculate_gap") then
    return_status = system(arg)
  elseif (arg == "combine_forces") then
    return_status = system(arg)
  elseif (arg == "compare_kpoints") then
    return_status = system(arg)
  elseif (arg == "construct_finite_displacement") then
    return_status = system(arg)
  elseif (arg == "construct_matrix_force_cnsts") then
    return_status = system(arg)
  elseif (arg == "construct_supercell") then
    return_status = system(arg)
  elseif (arg == "convert_forces_from_Rybohr_to_eVang") then
    return_status = system(arg)
  elseif (arg == "equilibrium_frac") then
    return_status = system(arg)
  elseif (arg == "fourier_interpolation") then
    return_status = system(arg)
  elseif (arg == "generate_amplitudes") then
    return_status = system(arg)
  elseif (arg == "generate_kgrid") then
    return_status = system(arg)
  elseif (arg == "generate_quadratic_configurations") then
    return_status = system(arg)
  elseif (arg == "generate_sc_path") then
    return_status = system(arg)
  elseif (arg == "generate_supercell_kpoint_mesh_qe") then
    return_status = system(arg)
  elseif (arg == "generate_supercells") then
    return_status = system(arg)
  elseif (arg == "lte") then
    return_status = system(arg)
  elseif (arg == "lte_lower") then
    return_status = system(arg)
  elseif (arg == "quadratic_spline") then
    return_status = system(arg)
  elseif (arg == "vscf_1d") then
    return_status = system(arg)
  ! wrappers for shell scripts
  elseif (arg == "anharmonic") then
    return_status = system(arg//".sh")
  elseif (arg == "bs_quadratic") then
    return_status = system(arg//".sh")
  elseif (arg == "clear_all") then
    return_status = system(arg//".sh")
  elseif (arg == "convert_harmonic") then
    return_status = system(arg//".sh")
  elseif (arg == "convert_quadratic") then
    return_status = system(arg//".sh")
  elseif (arg == "dyn_mats") then
    return_status = system(arg//".sh")
  elseif (arg == "eigenval_castep_to_bands") then
    return_status = system(arg//".sh")
  elseif (arg == "eigenval_vasp_to_bands") then
    return_status = system(arg//".sh")
  elseif (arg == "fetch_forces_castep") then
    return_status = system(arg//".sh")
  elseif (arg == "fetch_forces_qe") then
    return_status = system(arg//".sh")
  elseif (arg == "hartree_to_eV") then
    return_status = system(arg//".sh")
  elseif (arg == "lte_harmonic") then
    return_status = system(arg//".sh")
  elseif (arg == "rutgers_run_harmonic") then
    return_status = system(arg//".sh")
  elseif (arg == "setup_harmonic") then
    return_status = system(arg//".sh")
  elseif (arg == "setup_quadratic") then
    return_status = system(arg//".sh")
  elseif (arg == "structure_to_castep") then
    return_status = system(arg//".sh")
  elseif (arg == "structure_to_qe") then
    return_status = system(arg//".sh")
  elseif (arg == "structure_to_vasp") then
    return_status = system(arg//".sh")
  elseif (arg == "tcm_cleanup_anharmonic") then
    return_status = system(arg//".sh")
  elseif (arg == "tcm_cleanup_bs") then
    return_status = system(arg//".sh")
  elseif (arg == "tcm_cluster_run_harmonic") then
    return_status = system(arg//".sh")
  elseif (arg == "tcm_cluster_run_quadratic") then
    return_status = system(arg//".sh")
  ! wrappers for python scripts
  elseif (arg == "get_kpoints") then
    return_status = system(arg//".py")
  ! unrecognised argument
  else
    write(*,*) "Unrecognised argument : "//arg
  endif
  
  deallocate(arg)
  deallocate(args)

end program
