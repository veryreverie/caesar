! ======================================================================
! Writes a Quantum Espresso .fc file.
! ======================================================================
module force_constants_file_module
  use utils_module
  use structure_module
  use normal_mode_module
  
  use cartesian_hessian_module
  implicit none
  
  private
  
  public :: make_qe_force_constants_filename
  public :: write_qe_force_constants_file
  public :: read_qe_force_constants_file
contains
  
function make_qe_force_constants_filename(seedname) result(output)
  implicit none
  
  type(String), intent(in) :: seedname
  type(String)             :: output
  
  output = seedname//'.fc'
end function

subroutine write_qe_force_constants_file(fc_file,hessian,structure,supercell)
  implicit none
  
  type(OFile),            intent(inout) :: fc_file
  type(CartesianHessian), intent(in)    :: hessian
  type(StructureData),    intent(in)    :: structure
  type(StructureData),    intent(in)    :: supercell
  
  real(dp) :: alat
  
  type(String)              :: atom_species
  type(String), allocatable :: species(:)
  real(dp),     allocatable :: masses(:)
  
  integer :: grid(3)
  
  integer :: rvector(3)
  
  real(dp) :: matrix(3,3)
  real(dp) :: none_mass_reduced_force
  
  integer, allocatable :: ids(:,:,:,:,:,:)
  
  ! Cartesian directions.
  integer :: i1,i2
  
  ! Atom labels and atoms.
  integer        :: a1,a2
  type(AtomData) :: atom1,atom2
  
  ! Primitive atom labels.
  integer :: p1,p2
  
  ! R-vector labels.
  integer :: r1,r2,r3
  
  integer :: i,ialloc
  
  ! Construct supercell grid.
  grid = [(supercell%supercell%element(i,i), i=1, 3)]
  
  ! Construct mapping from r1,r2,r3,p1,p2 to a1,a2, where:
  !    - r1,r2 and r3 are R-vector labels,
  !    - p1 and p2 are primitive cell atom labels,
  !    - a1 and a2 are supercell atom labels.
  ids = construct_ids(supercell)
  
  ! Calculate 'alat', the length of the 'a' lattice vector.
  alat = l2_norm(vec([1,0,0])*structure%lattice)
  
  ! Construct species mapping.
  species = [String::]
  masses  = [real::]
  do p1=1,size(structure%atoms)
    atom_species = structure%atoms(p1)%species()
    if (.not. any(atom_species==species)) then
      species = [species, structure%atoms(p1)%species()]
      masses  = [masses, structure%atoms(p1)%mass()]
    endif
  enddo
  
  ! Write output file.
  call fc_file%print_line( size(species)         //' '// &
                         & size(structure%atoms) //' '// &
                         & 0                     //' '// &
                         & alat                  //' '// &
                         & 0.0_dp                //' '// &
                         & 0.0_dp                //' '// &
                         & 0.0_dp                //' '// &
                         & 0.0_dp                //' '// &
                         & 0.0_dp                        )
  call fc_file%print_lines(structure%lattice/alat)
  do i=1,size(species)
    call fc_file%print_line( i                    //" '"// &
                           & species(i)           //"' "// &
                           & masses(i)*RYDBERG_MASS_PER_ME )
  enddo
  do p1=1,size(structure%atoms)
    atom_species = structure%atoms(p1)%species()
    call fc_file%print_line(                                   &
       & p1                                            //' '// &
       & first(atom_species==species)                  //' '// &
       & structure%atoms(p1)%cartesian_position()/alat         )
  enddo
  call fc_file%print_line('F')
  call fc_file%print_line(grid)
  
  do i1=1,3
    do i2=1,3
      do p1=1,supercell%no_atoms_prim
        do p2=1,supercell%no_atoms_prim
          do r3=1,grid(3)
            do r2=1,grid(2)
              call fc_file%print_line(i1//' '//i2//' '//p1//' '//p2)
              do r1=1,grid(1)
                a1 = ids(1,r1,r2,r3,p1,p2)
                a2 = ids(2,r1,r2,r3,p1,p2)
                matrix = dble(hessian%elements( supercell%atoms(a1), &
                                              & supercell%atoms(a2)  ))
                none_mass_reduced_force = -matrix(i1,i2)                   &
                                      & * sqrt( structure%atoms(p1)%mass() &
                                      &       * structure%atoms(p2)%mass() )
                call fc_file%print_line(                            &
                   & r1                                     //' '// &
                   & r2                                     //' '// &
                   & r3                                     //' '// &
                   &  none_mass_reduced_force * RYDBERG_PER_HARTREE )
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
end subroutine

function read_qe_force_constants_file(fc_file,structure,supercell) &
   & result(output)
  implicit none
  
  type(IFile),         intent(in) :: fc_file
  type(StructureData), intent(in) :: structure
  type(StructureData), intent(in) :: supercell
  type(CartesianHessian)          :: output
  
  integer :: grid(3)
  
  integer, allocatable :: ids(:,:,:,:,:,:)
  
  type(String), allocatable :: lines(:)
  
  integer :: no_species
  integer :: header_length
  integer :: line_no
  
  real(dp),         allocatable :: elements(:,:,:,:)
  type(RealMatrix), allocatable :: matrix_elements(:,:)
  
  ! Atom labels and atoms.
  integer        :: a1,a2
  type(AtomData) :: atom1,atom2
  
  ! Cartesian component labels.
  integer :: i1,i2
  
  ! Primitive atom labels.
  integer :: p1,p2
  
  ! R-vector labels.
  integer :: r1,r2,r3
  
  integer :: i,ialloc
  
  ! Construct supercell grid.
  grid = [(supercell%supercell%element(i,i), i=1, 3)]
  
  ! Construct mapping from r1,r2,r3,p1,p2 to a1,a2, where:
  !    - r1,r2 and r3 are R-vector labels,
  !    - p1 and p2 are primitive cell atom labels,
  !    - a1 and a2 are supercell atom labels.
  ids = construct_ids(supercell)
  
  lines = lower_case(fc_file%lines())
  
  ! Check no. atoms and grid.
  if (int(token(lines(1),2))/=structure%no_atoms) then
    call print_line(ERROR//': .fc file contains the wrong number of atoms.')
    call err()
  endif
  
  no_species = int(token(lines(1),1))
  header_length = 6 + no_species + structure%no_atoms
  
  if (any(int(tokens(lines(header_length),1,3))/=grid)) then
    call print_line(ERROR//': .fc file contains unexpected q-point grid.')
  endif
  
  ! Read force constants.
  allocate( elements(3,3,supercell%no_atoms,supercell%no_atoms), &
          & stat=ialloc); call err(ialloc)
  line_no = header_length
  do i1=1,3
    do i2=1,3
      do p1=1,supercell%no_atoms_prim
        do p2=1,supercell%no_atoms_prim
          do r3=1,grid(3)
            do r2=1,grid(2)
              line_no = line_no+1
              do r1=1,grid(1)
                line_no = line_no + 1
                elements(i1,i2,a1,a2) = -dble(token(lines(line_no),4))     &
                                    & / sqrt( structure%atoms(p1)%mass()   &
                                    &       * structure%atoms(p2)%mass() ) &
                                    & / RYDBERG_PER_HARTREE
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  
  allocate( matrix_elements(supercell%no_atoms,supercell%no_atoms), &
          & stat=ialloc); call err(ialloc)
  do a1=1,supercell%no_atoms
    do a2=1,supercell%no_atoms
      matrix_elements(a1,a2) = mat(elements(:,:,a1,a2))
    enddo
  enddo
  
  output = CartesianHessian(structure,matrix_elements)
end function

! Construct mapping from r1,r2,r3,p1,p2 to a1,a2, where:
!    - r1,r2 and r3 are R-vector labels,
!    - p1 and p2 are primitive cell atom labels,
!    - a1 and a2 are supercell atom labels.
function construct_ids(supercell) result(output)
  implicit none
  
  type(StructureData), intent(in) :: supercell
  integer, allocatable            :: output(:,:,:,:,:,:)
  
  integer :: grid(3)
  
  type(AtomData) :: atom1,atom2
  
  integer :: rvector(3)
  
  integer :: a1,a2,i,ialloc
  
  grid = [(supercell%supercell%element(i,i), i=1, 3)]
  
  allocate( output( 2,                          &
          &         grid(1),                    &
          &         grid(2),                    &
          &         grid(3),                    &
          &         supercell%no_atoms_prim,    &
          &         supercell%no_atoms_prim  ), &
          & stat=ialloc); call err(ialloc)
  do a1=1,size(supercell%atoms)
    atom1 = supercell%atoms(a1)
    if (atom1%prim_id()/=atom1%id()) then
      cycle
    endif
    
    do a2=1,size(supercell%atoms)
      atom2 = supercell%atoms(a2)
      rvector = modulo(int(supercell%rvectors(atom2%rvec_id()))-1, grid) + 1
      
      output( 1,               &
            & rvector(1),      &
            & rvector(2),      &
            & rvector(3),      &
            & atom1%prim_id(), &
            & atom2%prim_id()  ) = a1
      output( 2,               &
            & rvector(1),      &
            & rvector(2),      &
            & rvector(3),      &
            & atom1%prim_id(), &
            & atom2%prim_id()  ) = a2
    enddo
  enddo
end function
end module
