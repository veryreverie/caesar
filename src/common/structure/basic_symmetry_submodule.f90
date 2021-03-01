submodule (caesar_basic_symmetry_module) caesar_basic_symmetry_submodule
  use caesar_structure_module
contains

module procedure new_BasicSymmetry
  this%id          = id
  this%tensor      = tensor
  this%translation = translation
  this%atom_group  = atom_group
  this%rvectors    = rvectors
end procedure

module procedure new_BasicSymmetry_SpglibSymmetries
  ! Variables for calculating atom groups and R-vectors.
  type(Group),      allocatable :: atom_groups(:)
  type(IntVector),  allocatable :: rvectors(:,:)
  type(RealVector)              :: transformed_position
  type(RealVector), allocatable :: kj_displacements(:)
  type(IntVector),  allocatable :: kj_rvectors(:)
  real(dp),         allocatable :: kj_distances(:)
  integer,          allocatable :: atom_group(:)
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  ! --------------------------------------------------
  ! Output symmetry information.
  ! --------------------------------------------------
  call print_line('Space group: '//input%international_symbol)
  
  ! --------------------------------------------------
  ! Calculate atom group and R-vectors.
  ! --------------------------------------------------
  ! The symmetry S transforms atomic equilibrium positions, {r_i} in
  !    fractional co-ordinates as:
  ! r_i -> r_j + R
  !
  ! atom_group * i = j
  ! rvectors(i) = R
  allocate( atom_groups(size(input)),           &
          & rvectors(size(atoms), size(input)), &
          & kj_displacements(size(atoms)),      &
          & kj_rvectors(size(atoms)),           &
          & kj_distances(size(atoms)),          &
          & atom_group(size(atoms)),            &
          & stat=ialloc); call err(ialloc)
  do i=1,size(input)
    do j=1,size(atoms)
      ! Calculate the position of the transformed atom.
      transformed_position = input%tensors(i)               &
                         & * atoms(j)%fractional_position() &
                         & + input%translations(i)
      
      do k=1,size(atoms)
        ! Calculate the displacement from atom k to the transfomed atom j.
        kj_displacements(k) = transformed_position &
                          & - atoms(k)%fractional_position()
        ! Calculate the nearest R-vector to the displacement.
        kj_rvectors(k) = nint(kj_displacements(k))
        ! Calculate the distance between the displacement and the nearest
        !    R-vector.
        kj_distances(k) = l2_norm(kj_displacements(k)-kj_rvectors(k))
      enddo
      
      ! Check that exactly one k->j distance is small.
      if (count(kj_distances<symmetry_precision)==0) then
        call print_line(ERROR//': A symmetry maps an atom onto space containing no atoms. Try increasing &
                       &symmetry_precision.')
        call err()
      elseif (count(kj_distances<symmetry_precision)>1) then
        call print_line(ERROR//': A symmetry maps an atom onto more than one atom. Try reducing &
                       &symmetry_precision.')
        call err()
      endif
      
      ! Identify the atom which symmetry i maps atom j to.
      k = minloc(kj_distances,1)
      
      if (atoms(k)%species()/=atoms(j)%species()) then
        call print_line(ERROR//': symmetry operation maps between atoms of &
           &different species.')
        call err()
      endif
      
      atom_group(j) = atoms(k)%id()
      rvectors(j,i) = kj_rvectors(k)
    enddo
    atom_groups(i) = Group(atom_group)
  enddo
  
  ! --------------------------------------------------
  ! Construct output.
  ! --------------------------------------------------
  
  allocate(output(size(input)), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    output(i) = BasicSymmetry( id          = i,                     &
                             & tensor      = input%tensors(i),      &
                             & translation = input%translations(i), &
                             & atom_group  = atom_groups(i),        &
                             & rvectors    = rvectors(:,i)          )
  enddo
end procedure

module procedure read_BasicSymmetry
  integer                      :: id
  type(IntMatrix)              :: tensor
  type(RealVector)             :: translation
  integer,         allocatable :: atom_group(:)
  type(IntVector), allocatable :: rvectors(:)
  
  type(String), allocatable :: line(:)
  
  integer :: i,ialloc
  
  select type(this); type is(BasicSymmetry)
    line = split_line(input(1))
    id = int(line(2))
    tensor = IntMatrix(input(3:5))
    translation = RealVector(input(7))
    
    allocate( atom_group(size(input)-8), &
            & rvectors(size(input)-8),   &
            & stat=ialloc); call err(ialloc)
    do i=1,size(atom_group)
      line = split_line(input(8+i))
      atom_group(i) = int(line(5))
      rvectors(i) = vec(int(line(8:10)))
    enddo
    
    this = BasicSymmetry(id,tensor,translation,Group(atom_group),rvectors)
  class default
    call err()
  end select
end procedure

module procedure write_BasicSymmetry
  type(String), allocatable :: atom_strings(:)
  
  integer :: i,ialloc
  
  allocate( atom_strings(size(this%rvectors)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(this%rvectors)
    atom_strings(i) = 'Atom '//i//' -> Atom '//this%atom_group*i// &
                    & ' + R-vector '//this%rvectors(i)
  enddo
  
  select type(this); type is(BasicSymmetry)
    output = [ 'Operation '//this%id,                          &
             & str('Tensor in fractional co-ordinates:'),      &
             & str(this%tensor),                               &
             & str('Translation in fractional co-ordinates:'), &
             & str(this%translation),                          &
             & str('Effect on atoms:'),                        &
             & atom_strings                                    ]
  class default
    call err()
  end select
end procedure

module procedure new_BasicSymmetry_Strings
  call this%read(input)
end procedure

module procedure new_BasicSymmetry_StringArray
  this = BasicSymmetry(str(input))
end procedure
end submodule
