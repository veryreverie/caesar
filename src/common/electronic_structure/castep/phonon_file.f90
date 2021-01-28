! ======================================================================
! Writes a Castep .phonon file.
! ======================================================================
module caesar_phonon_file_module
  use caesar_utils_module
  use caesar_structure_module
  use caesar_normal_mode_module
  implicit none
  
  private
  
  public :: make_castep_phonon_filename
  public :: write_castep_phonon_file
  public :: read_castep_phonon_file
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
        ! The eigenvector is listed as [Re(x),Im(x),Re(y),Im(y),Re(z),Im(z)].
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

function read_castep_phonon_file(phonon_file,structure,qpoints) result(output)
  implicit none
  
  type(IFile),         intent(in) :: phonon_file
  type(StructureData), intent(in) :: structure
  type(QpointData),    intent(in) :: qpoints(:)
  type(ComplexMode), allocatable  :: output(:,:)
  
  type(String), allocatable :: lines(:)
  
  integer :: block_length
  integer :: header_length
  
  type(RealVector), allocatable :: input_qpoints(:)
  type(RealVector), allocatable :: file_qpoints(:)
  real(dp)                      :: difference(3)
  integer,          allocatable :: file_to_input(:)
  
  real(dp)                         :: frequency
  real(dp)                         :: spring_constant
  real(dp)                         :: atom_vector(6)
  type(ComplexVector), allocatable :: unit_vector(:)
  
  integer :: i,j,k,ialloc
  
  lines = lower_case(phonon_file%lines())
  
  ! Check the header is consistent with qpoints and structure.
  if (any(tokens(lines(1),1,2)/=tokens('begin  header'))) then
    call print_line(ERROR//': .phonon file line 1 is not "BEGIN header".')
    call err()
  elseif (any(tokens(lines(2),1,3)/=tokens('number of ions'))) then
    call print_line(ERROR//': .phonon file line 2 does not begin &
       &"Number of ions".')
    call err()
  elseif (int(token(lines(2),4))/=structure%no_atoms) then
    call print_line(ERROR//': .phonon file contains an unexpected number of &
       &atoms ("ions").')
    call err()
  elseif (any(tokens(lines(3),1,3)/=tokens('number of branches'))) then
    call print_line(ERROR//': .phonon file line 3 does not begin &
       &"Number of branches".')
    call err()
  elseif (int(token(lines(3),4))/=structure%no_modes) then
    call print_line(ERROR//': .phonon file contains an unexpected number of &
       &modes ("branches").')
    call err()
  elseif (any(tokens(lines(4),1,3)/=tokens('number of wavevectors'))) then
    call print_line(ERROR//': .phonon file line 4 does not begin &
       &"Number of wavevectors".')
    call err()
  elseif (int(token(lines(4),4))/=size(qpoints)) then
    call print_line(ERROR//': .phonon file contains an unexpected number of &
       &q-points ("wavevectors").')
    call err()
  elseif (any(tokens(lines(5),1,2)/=tokens('frequencies in'))) then
    call print_line(ERROR//': .phonon file line 5 does not begin &
       &"Frequencies in".')
    call err()
  elseif (token(lines(5),3)/='cm-1') then
    call print_line(ERROR//': .phonon file has frequencies in a unit other &
       &than inverse cm. This is not currently supported.')
    call err()
  endif
  
  header_length = 11+structure%no_atoms
  if (any(tokens(lines(header_length),1,2)/=tokens('end header'))) then
    call print_line(ERROR//': .phonon file line '//header_length//' is not &
       &"END header".')
    call err()
  endif
  
  ! Check the file length is as expected.
  block_length = (structure%no_atoms+1)*structure%no_modes+3
  if (size(lines)/=header_length+block_length*size(qpoints)) then
    call print_line(ERROR//': .phonon file does not contain the expected &
       &number of lines.')
  endif
  
  ! Match q-points between qpoints and the .phonon file.
  input_qpoints = dblevec(qpoints%qpoint)
  file_qpoints = [(                                                        &
     & vec(dble(tokens(lines(header_length+block_length*(i-1)+1), 3, 5))), &
     & i=1,                                                                &
     & size(qpoints)                                                       )]
  file_to_input = [(0,i=1,size(qpoints))]
  do i=1,size(file_qpoints)
    do j=1,size(input_qpoints)
      difference = dble(input_qpoints(j)-file_qpoints(i))
      if (sum(abs(difference-nint(difference)))<0.5_dp/size(qpoints)) then
        file_to_input(i) = j
        cycle
      endif
    enddo
  enddo
  
  if (size(set(file_to_input))/=size(file_to_input)) then
    call print_line(ERROR//': Unable to establish mapping between q-points in &
       &.phonon file and input q-points.')
    call err()
  endif
  
  ! Read in normal mode frequencies and eigenvectors.
  allocate( unit_vector(structure%no_atoms),          &
          & output(structure%no_modes,size(qpoints)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(qpoints)
    do j=1,structure%no_modes
      frequency = dble(token( lines(header_length+block_length*(i-1)+1+j),    &
              &               2                                            )) &
              & / INVERSE_CM_PER_HARTREE
      
      if (frequency>=0) then
        spring_constant = frequency**2
      else
        spring_constant = -frequency**2
      endif
      
      do k=1,structure%no_atoms
        ! The eigenvector is listed as [Re(x),Im(x),Re(y),Im(y),Re(z),Im(z)].
        atom_vector = dble(tokens( lines( header_length               &
                                 &      + block_length*(i-1)          &
                                 &      + structure%no_atoms*(j-1)    &
                                 &      + 3+structure%no_modes+k   ), &
                                 & 3,                                 &
                                 & 8                                  ))
        unit_vector(k) = vec([ cmplx(atom_vector(1),atom_vector(2),dp), &
                             & cmplx(atom_vector(3),atom_vector(4),dp), &
                             & cmplx(atom_vector(5),atom_vector(6),dp)  ])
      enddo
      
      output(j,file_to_input(i)) = ComplexMode(       &
         & id                 = 0,                    &
         & paired_id          = 0,                    &
         & frequency          = frequency,            &
         & spring_constant    = spring_constant,      &
         & soft_mode          = frequency<-1.0e-6_dp, &
         & translational_mode = .false.,              &
         & unit_vector        = unit_vector,          &
         & qpoint_id          = 0,                    &
         & paired_qpoint_id   = 0,                    &
         & subspace_id        = 0                     )
    enddo
  enddo
end function
end module
