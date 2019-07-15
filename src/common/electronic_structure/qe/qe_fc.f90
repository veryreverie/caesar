! ======================================================================
! Reads and Writes Quantum Espresso .fc files.
! ======================================================================
module qe_fc_module
  use utils_module
  use structure_module
  use normal_mode_module
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
  
  real(dp) :: matrix(3,3)
  real(dp) :: none_mass_reduced_force
  
  integer, allocatable :: ids(:,:,:,:)
  
  ! Cartesian directions.
  integer :: i1,i2
  
  ! Atom labels and atoms.
  integer        :: id
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
  ids = construct_atom_ids(supercell)
  
  ! Calculate 'alat', the length of the 'a' lattice vector.
  alat = l2_norm(vec([1,0,0])*structure%lattice)
  
  ! Construct species mapping.
  allocate( species(0), &
          & masses(0),  &
          & stat=ialloc); call err(ialloc)
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
          call fc_file%print_line(i1//' '//i2//' '//p1//' '//p2)
          do r3=1,grid(3)
            do r2=1,grid(2)
              do r1=1,grid(1)
                id = ids(r1,r2,r3,p2)
                matrix = dble(hessian%elements( supercell%atoms(p1),    &
                     &                          supercell%atoms(id)  )) &
                     & / supercell%sc_size
                none_mass_reduced_force = -matrix(i1,i2)                   &
                                      & * sqrt( structure%atoms(p1)%mass() &
                                      &       * structure%atoms(p2)%mass() )
                call fc_file%print_line(                           &
                   & r1                                    //' '// &
                   & r2                                    //' '// &
                   & r3                                    //' '// &
                   & none_mass_reduced_force * RYDBERG_PER_HARTREE )
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
end subroutine

function read_qe_force_constants_file(directory,seedname,supercell) &
   & result(output)
  implicit none
  
  type(String),        intent(in) :: directory
  type(String),        intent(in) :: seedname
  type(StructureData), intent(in) :: supercell
  type(CartesianHessian)          :: output
  
  type(IFile) :: fc_file
  
  integer :: grid(3)
  
  integer, allocatable :: ids(:,:,:,:)
  
  type(String), allocatable :: lines(:)
  
  integer :: no_species
  integer :: header_length
  integer :: line_no
  
  real(dp),         allocatable :: elements(:,:,:,:)
  type(RealMatrix), allocatable :: matrix_elements(:,:)
  
  ! Atom labels and atoms.
  integer        :: id1,id2
  type(AtomData) :: atom1,atom2
  
  ! Cartesian component labels.
  integer :: i1,i2
  
  ! Primitive atom labels.
  integer :: p1,p2
  
  ! R-vector labels.
  integer :: r1,r2,r3
  
  integer :: i,ialloc
  
  ! Open file.
  fc_file = IFile(directory//'/'//make_qe_force_constants_filename(seedname))
  lines = lower_case(fc_file%lines())
  
  ! Construct supercell grid.
  grid = [(supercell%supercell%element(i,i), i=1, 3)]
  
  ! Construct mapping from r1,r2,r3,p1,p2 to a1,a2, where:
  !    - r1,r2 and r3 are R-vector labels,
  !    - p1 and p2 are primitive cell atom labels,
  !    - a1 and a2 are supercell atom labels.
  ids = construct_atom_ids(supercell)
  
  ! Check no. atoms and grid.
  if (int(token(lines(1),2))/=supercell%no_atoms_prim) then
    call print_line(ERROR//': .fc file contains the wrong number of atoms.')
    call err()
  endif
  
  no_species = int(token(lines(1),1))
  header_length = 6 + no_species + supercell%no_atoms_prim
  
  if (any(int(tokens(lines(header_length),1,3))/=grid)) then
    call print_line(ERROR//': .fc file contains unexpected q-point grid.')
  endif
  
  ! Read force constants.
  allocate( elements(3,3,supercell%no_atoms_prim,supercell%no_atoms), &
          & stat=ialloc); call err(ialloc)
  line_no = header_length
  do i1=1,3
    do i2=1,3
      do p1=1,supercell%no_atoms_prim
        do p2=1,supercell%no_atoms_prim
          line_no = line_no+1
          do r3=1,grid(3)
            do r2=1,grid(2)
              do r1=1,grid(1)
                line_no = line_no + 1
                id2 = ids(r1,r2,r3,p2)
                elements(i1,i2,p1,id2) = -dble(token(lines(line_no),4))     &
                                     & / sqrt( supercell%atoms(p1)%mass()   &
                                     &       * supercell%atoms(p2)%mass() ) &
                                     & / RYDBERG_PER_HARTREE                &
                                     & * supercell%sc_size
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  
  allocate( matrix_elements(supercell%no_atoms_prim,supercell%no_atoms), &
          & stat=ialloc); call err(ialloc)
  do id1=1,supercell%no_atoms_prim
    do id2=1,supercell%no_atoms
      matrix_elements(id1,id2) = mat(elements(:,:,id1,id2))
    enddo
  enddo
  
  output = CartesianHessian( supercell      = supercell,       &
                           & elements       = matrix_elements, &
                           & check_symmetry = .false.          )
end function

! Construct mapping from R-vector and primitive id to atom id.
function construct_atom_ids(supercell) result(output)
  implicit none
  
  type(StructureData), intent(in) :: supercell
  integer, allocatable            :: output(:,:,:,:)
  
  integer :: grid(3)
  
  type(AtomData) :: atom
  
  integer :: rvector(3)
  
  integer :: id,i,ialloc
  
  grid = [(supercell%supercell%element(i,i), i=1, 3)]
  
  allocate( output(grid(1), grid(2), grid(3), supercell%no_atoms_prim), &
          & stat=ialloc); call err(ialloc)
  do id=1,size(supercell%atoms)
    atom = supercell%atoms(id)
    
    ! Quantum Espresso 1-indexes R-vector components, so the primitive cell is
    !    at (1,1,1).
    ! Quantum Espresso also adopts the opposite convention to Caesar for
    !    the sign of the R-vector.
    ! If the supercell grid is (A,B,C), then the R-vector at (a,b,c)
    !    with a in [0,A), b in [0,B) and c in [0,C) is stored in the .fc file
    !    as (A-a+1,B-b+1,C-c+1).
    rvector = modulo(grid-int(supercell%rvectors(atom%rvec_id())), grid) + 1
    
    output(rvector(1), rvector(2), rvector(3), atom%prim_id()) = id
  enddo
end function
end module
