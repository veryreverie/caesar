! ======================================================================
! Writes a Castep .phonon file.
! ======================================================================
module phonon_file_module
  use utils_module
  use structure_module
  use normal_mode_module
  implicit none
  
  private
  
  public :: make_castep_phonon_filename
  public :: write_castep_phonon_file
contains

function make_castep_phonon_filename(seedname) result(output)
  implicit none
  
  type(String), intent(in) :: seedname
  type(String)             :: output
  
  output = seedname//'.phonon'
end function

subroutine write_castep_phonon_file(phonon_file,complex_modes,qpoints, &
   & structure)
  implicit none
  
  type(OFile),           intent(inout) :: phonon_file
  type(ComplexMode),     intent(in)    :: complex_modes(:,:)
  type(QpointData),      intent(in)    :: qpoints(:)
  type(StructureData),   intent(in)    :: structure
  
  type(String), allocatable :: lines(:)
  
  real(dp)     :: real_part(3)
  real(dp)     :: imag_part(3)
  
  integer :: i,j,k
  
  call set_print_settings(PrintSettings( floating_point_format = str('f'), &
                                       & integer_digits        = 4,        &
                                       & decimal_places        = 6         ))
  
  lines = [ str('BEGIN header'),                           &
          & 'Number of ions         '//structure%no_atoms, &
          & 'Number of branches     '//structure%no_modes, &
          & 'Number of wavevectors  '//size(qpoints),      &
          & str('Frequencies in         cm-1'),            &
          & str('Unit cell vectors (A)'),                  &
          & str(structure%lattice * ANGSTROM_PER_BOHR),    &
          & str('Fractional Co-ordinates')                 ]
  do i=1,structure%no_atoms
    lines = [ lines,                                             &
            & left_pad(i,'     ',' ')//' '               //' '// &
            & structure%atoms(i)%fractional_position()   //' '// &
            & '  '//structure%atoms(i)%species()         //' '// &
            & structure%atoms(i)%mass() * KG_PER_ME / KG_PER_AMU ]
  enddo
  lines = [ lines,            &
          & str('END header') ]
  
  do i=1,size(lines)
    lines(i) = ' '//lines(i)
  enddo
  
  do i=1,size(qpoints)
    lines = [ lines,                             &
            & '    q-pt= '                    // &
            & left_pad(i,'    ',' ')     //' '// &
            & dblevec(qpoints(i)%qpoint) //' '// &
            & 1.0_dp/size(qpoints)               ]
    do j=1,size(complex_modes,1)
      lines = [ lines,                                                &
              & left_pad(j,'       ',' ')                  //'    '// &
              & complex_modes(j,i)%frequency * INVERSE_CM_PER_HARTREE ]
    enddo
    lines = [ lines,                      &
            & str('Phonon Eigenvectors'), &
            & str('Mode Ion X Y Z')       ]
    do j=1,size(complex_modes,1)
      do k=1,structure%no_atoms
        real_part = dble(real( complex_modes(j,i)%unit_vector(k)))
        imag_part = dble(aimag(complex_modes(j,i)%unit_vector(k)))
        lines = [ lines,               &
                & j            //' '// &
                & k            //' '// &
                & real_part(1) //' '// &
                & imag_part(1) //' '// &
                & real_part(2) //' '// &
                & imag_part(2) //' '// &
                & real_part(3) //' '// &
                & imag_part(3)         ]
      enddo
    enddo
  enddo
  
  call phonon_file%print_lines(lines)
  
  call unset_print_settings()
end subroutine
end module
