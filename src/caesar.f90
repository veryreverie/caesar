module utils
contains

  ! returns an array containing the command line arguments.
  function command_line_args() result(o_args)
    implicit none

    integer                        :: i         ! loop index
    integer                        :: arg_count ! no. of command line args
    character(len=32), allocatable :: o_args(:) ! return value
    
    ! read in the number of arguments
    arg_count = iargc()
    ! allocate the return array
    allocate (o_args(arg_count))
    ! read the arguments into the array
    do i=1,arg_count
      call getarg(i, o_args(i))
    enddo
    
    return
  end function

end module

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
    write(*,*) "option selected: ", arg
    return_status = system(arg)
  elseif (arg == "calculate_anharmonic") then
    write(*,*) "option selected: ", arg
    return_status = system(arg)
  elseif (arg == "calculate_bs") then
    write(*,*) "option selected: ", arg
    return_status = system(arg)
  elseif (arg == "calculate_gap") then
    write(*,*) "option selected: ", arg
    return_status = system(arg)
  elseif (arg == "combine_forces") then
    write(*,*) "option selected: ", arg
    return_status = system(arg)
  elseif (arg == "compare_kpoints") then
    write(*,*) "option selected: ", arg
    return_status = system(arg)
  elseif (arg == "construct_finite_displacement") then
    write(*,*) "option selected: ", arg
    return_status = system(arg)
  elseif (arg == "construct_matrix_force_cnsts") then
    write(*,*) "option selected: ", arg
    return_status = system(arg)
  elseif (arg == "construct_supercell") then
    write(*,*) "option selected: ", arg
    return_status = system(arg)
  elseif (arg == "convert_forces_from_Rybohr_to_eVang") then
    write(*,*) "option selected: ", arg
    return_status = system(arg)
  elseif (arg == "equilibrium_frac") then
    write(*,*) "option selected: ", arg
    return_status = system(arg)
  elseif (arg == "fourier_interpolation") then
    write(*,*) "option selected: ", arg
    return_status = system(arg)
  elseif (arg == "generate_amplitudes") then
    write(*,*) "option selected: ", arg
    return_status = system(arg)
  elseif (arg == "generate_kgrid") then
    write(*,*) "option selected: ", arg
    return_status = system(arg)
  elseif (arg == "generate_quadratic_configurations") then
    write(*,*) "option selected: ", arg
    return_status = system(arg)
  elseif (arg == "generate_sc_path") then
    write(*,*) "option selected: ", arg
    return_status = system(arg)
  elseif (arg == "generate_supercells") then
    write(*,*) "option selected: ", arg
    return_status = system(arg)
  elseif (arg == "lte") then
    write(*,*) "option selected: ", arg
    return_status = system(arg)
  elseif (arg == "lte_lower") then
    write(*,*) "option selected: ", arg
    return_status = system(arg)
  elseif (arg == "quadratic_spline") then
    write(*,*) "option selected: ", arg
    return_status = system(arg)
  elseif (arg == "vscf_1d") then
    write(*,*) "option selected: ", arg
    return_status = system(arg)
  ! wrappers for shell scripts
  elseif (arg == "anharmonic") then
    write(*,*) "option selected: ", arg
    return_status = system(arg//".sh")
  elseif (arg == "bs_quadratic") then
    write(*,*) "option selected: ", arg
    return_status = system(arg//".sh")
  elseif (arg == "clear_all") then
    write(*,*) "option selected: ", arg
    return_status = system(arg//".sh")
  elseif (arg == "convert_harmonic") then
    write(*,*) "option selected: ", arg
    return_status = system(arg//".sh")
  elseif (arg == "convert_quadratic") then
    write(*,*) "option selected: ", arg
    return_status = system(arg//".sh")
  elseif (arg == "dyn_mats") then
    write(*,*) "option selected: ", arg
    return_status = system(arg//".sh")
  elseif (arg == "eigenval_castep_to_bands") then
    write(*,*) "option selected: ", arg
    return_status = system(arg//".sh")
  elseif (arg == "eigenval_vasp_to_bands") then
    write(*,*) "option selected: ", arg
    return_status = system(arg//".sh")
  elseif (arg == "fetch_forces_castep") then
    write(*,*) "option selected: ", arg
    return_status = system(arg//".sh")
  elseif (arg == "fetch_forces_qe") then
    write(*,*) "option selected: ", arg
    return_status = system(arg//".sh")
  elseif (arg == "hartree_to_eV") then
    write(*,*) "option selected: ", arg
    return_status = system(arg//".sh")
  elseif (arg == "lte_harmonic") then
    write(*,*) "option selected: ", arg
    return_status = system(arg//".sh")
  elseif (arg == "rutgers_run_harmonic") then
    write(*,*) "option selected: ", arg
    return_status = system(arg//".sh")
  elseif (arg == "setup_harmonic") then
    write(*,*) "option selected: ", arg
    return_status = system(arg//".sh")
  elseif (arg == "setup_quadratic") then
    write(*,*) "option selected: ", arg
    return_status = system(arg//".sh")
  elseif (arg == "structure_to_castep") then
    write(*,*) "option selected: ", arg
    return_status = system(arg//".sh")
  elseif (arg == "structure_to_qe") then
    write(*,*) "option selected: ", arg
    return_status = system(arg//".sh")
  elseif (arg == "structure_to_vasp") then
    write(*,*) "option selected: ", arg
    return_status = system(arg//".sh")
  elseif (arg == "tcm_cleanup_anharmonic") then
    write(*,*) "option selected: ", arg
    return_status = system(arg//".sh")
  elseif (arg == "tcm_cleanup_bs") then
    write(*,*) "option selected: ", arg
    return_status = system(arg//".sh")
  elseif (arg == "tcm_cluster_run_harmonic") then
    write(*,*) "option selected: ", arg
    return_status = system(arg//".sh")
  elseif (arg == "tcm_cluster_run_quadratic") then
    write(*,*) "option selected: ", arg
    return_status = system(arg//".sh")
  ! wrappers for python scripts
  elseif (arg == "get_kpoints") then
    write(*,*) "option selected: ", arg
    write(*,*) arg//".py"
    return_status = system(arg//".py")
  ! unrecognised argument
  else
    write(*,*) "Unrecognised argument : "//arg
  endif
  
  deallocate(arg)
  deallocate(args)

end program
