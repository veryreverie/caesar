module utils
contains

  ! returns an array containing the command line arguments.
  function command_line_args() result(o_args)
    implicit none

    integer                        :: i
    integer                        :: arg_count
    character(len=32), allocatable :: o_args(:)
    
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
  
  integer           :: i
  character(len=32), allocatable :: args(:)

  ! read in command line arguments
  args = command_line_args()
  
  if (size(args) == 0) then
    write(*,*) "No arguments given. For help, call caesar -h"
    stop
  elseif (args(1) == "-h" .or. args(1) == "--help") then
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
    stop
  elseif (args(1) == "combine_forces") then
    write(*,*) "option selected"
    stop
  elseif (args(1) == "compare_kpoints") then
    write(*,*) "option selected"
    stop
  elseif (args(1) == "construct_finite_displacement") then
    write(*,*) "option selected"
    stop
  elseif (args(1) == "construct_matrix_force_cnsts") then
    write(*,*) "option selected"
    stop
  elseif (args(1) == "construct_supercell") then
    write(*,*) "option selected"
    stop
  elseif (args(1) == "") then
    write(*,*) "option selected"
    stop
  elseif (args(1) == "") then
    write(*,*) "option selected"
    stop
  elseif (args(1) == "") then
    write(*,*) "option selected"
    stop
  elseif (args(1) == "") then
    write(*,*) "option selected"
    stop
  elseif (args(1) == "") then
    write(*,*) "option selected"
    stop
  elseif (args(1) == "") then
    write(*,*) "option selected"
    stop
  elseif (args(1) == "") then
    write(*,*) "option selected"
    stop
  elseif (args(1) == "") then
    write(*,*) "option selected"
    stop
  else
    write(*,*) "Unrecognised argument : "//args(1)
    stop
  endif
  
  deallocate (args)

end program
