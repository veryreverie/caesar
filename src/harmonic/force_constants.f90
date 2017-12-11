! ======================================================================
! Reads harmonic forces, and generates the matrix of force constants.
! ======================================================================
module force_constants_module
  use constants_module, only : dp
  use string_module
  use io_module
contains
  
! ----------------------------------------------------------------------
! Read in |f>, average over +/- |x>, and divide by | |x> |.
! ----------------------------------------------------------------------
! |x> is the collective vector of displacements, {xi}.
! |f> is the collective vector of forces, {fi}.
! Both are supercell%no_modes long.
!
! Under the harmonic approximation, the force constants, F, are defined as:
!    U = - 1/2 <x|F|x>
!
!    |f> = -dU/d|x> = F|x>
!
! Under the symmetry s:
!    |x> -> |xs> = Rs|x>
!    |f> -> |fs> = Rs|f>
!
! F can be found by minimising L:
!    L = sum(s)[ (|fs>-F|xs>)^2 ]
!      = sum(s)[ <fs|fs> - 2<fs|F|xs> + <xs|FF|xs> ]
!
! => 0 = dL/dF = -2 * sum(s)[ |fs><xs| - F|xs><xs| ]
! => F = sum(s)[|fs><xs|] . (sum(s)[|xs><xs|])^-1
!
! sum(s)[|xs><xs|] is block diagonal, so can be inverted in 3x3 blocks.
function read_forces(supercell,unique_directions,sdir,file_type,seedname) &
   & result(output)
  use structure_module
  use unique_directions_module
  use output_file_module
  use linear_algebra_module
  implicit none
  
  type(StructureData),    intent(in) :: supercell
  type(UniqueDirections), intent(in) :: unique_directions
  type(String),           intent(in) :: sdir
  type(String),           intent(in) :: file_type
  type(String),           intent(in) :: seedname
  type(RealVector), allocatable      :: output(:,:)
  
  ! DFT output data.
  type(String)     :: output_filename
  type(OutputFile) :: positive
  type(OutputFile) :: negative
  
  ! Direction information.
  integer      :: atom
  character(1) :: direction
  type(String) :: atom_string
  
  ! Temporary variables.
  integer          :: i,j,ialloc
  type(RealVector) :: total
  
  output_filename = make_output_filename(file_type,seedname)
  
  allocate( output(supercell%no_atoms, size(unique_directions)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(unique_directions)
    atom = unique_directions%atoms(i)
    direction = unique_directions%directions_char(i)
    atom_string = left_pad(atom,str(maxval(unique_directions%atoms)))
    
    positive = read_output_file(file_type,               &
       & sdir//'/atom.'//atom_string//'.+d'//direction// &
       & '/'//output_filename)
    negative = read_output_file(file_type,            &
       & sdir//'/atom.'//atom_string//'.-d'//direction// &
       & '/'//output_filename)
    
    do j=1,supercell%no_atoms
      output(j,i) = (positive%forces(j)-negative%forces(j)) / 0.02_dp
    enddo
    
    ! Enforce continuous translational symmetry,
    !    i.e. ensure sum(output,1)=0.
    total = vec([0.0_dp,0.0_dp,0.0_dp])
    do j=1,supercell%no_atoms
      total = total+output(j,i)
    enddo
    do j=1,supercell%no_atoms
      output(j,i) = output(j,i)-total/supercell%no_atoms
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Uses symmetry operations to construct force constants.
! ----------------------------------------------------------------------
function construct_force_constants(forces,supercell,unique_directions) &
   & result(output)
  use linear_algebra_module
  use structure_module
  use unique_directions_module
  use group_module
  implicit none
  
  type(RealVector),       intent(in) :: forces(:,:)
  type(StructureData),    intent(in) :: supercell
  type(UniqueDirections), intent(in) :: unique_directions
  type(RealMatrix), allocatable      :: output(:,:,:)
  
  ! Atom ids.
  integer :: atom_1,atom_1p,atom_2,atom_2p
  
  ! Parts of |x> and |f>.
  type(RealVector) :: x
  type(RealVector) :: f
  
  ! sum(s)[ |fs><xs| ].
  type(RealMatrix), allocatable :: fx(:,:)
  
  ! sum(s)[ |xs><xs| ] (diagonal blocks only).
  type(RealMatrix), allocatable :: xx(:)
  type(RealMatrix), allocatable :: xx_inverse(:)
  
  ! Force constants, F.
  type(RealMatrix), allocatable :: force_constants(:,:)
  
  ! R-vector information.
  type(Group), allocatable :: rvector_group(:)
  integer                  :: rvector_1
  integer                  :: rvector_2
  integer                  :: rvector
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  ! --------------------------------------------------
  ! Construct xx and fx.
  ! --------------------------------------------------
  allocate( xx(supercell%no_atoms),                    &
            fx(supercell%no_atoms,supercell%no_atoms), &
          & stat=ialloc); call err(ialloc)
  xx = mat([ 0.0_dp,0.0_dp,0.0_dp, &
           & 0.0_dp,0.0_dp,0.0_dp, &
           & 0.0_dp,0.0_dp,0.0_dp], 3,3)
  fx = mat([ 0.0_dp,0.0_dp,0.0_dp, &
           & 0.0_dp,0.0_dp,0.0_dp, &
           & 0.0_dp,0.0_dp,0.0_dp], 3,3)
  do i=1,size(supercell%symmetries)
    do j=1,size(unique_directions)
      atom_1 = unique_directions%atoms(j)
      atom_1p = supercell%symmetries(i)%atom_group * atom_1
      
      if (unique_directions%directions_char(j)=='x') then
        x = [ 1.0_dp, 0.0_dp, 0.0_dp ]
      elseif (unique_directions%directions_char(j)=='y') then
        x = [ 0.0_dp, 1.0_dp, 0.0_dp ]
      else
        x = [ 0.0_dp, 0.0_dp, 1.0_dp ]
      endif
      x = supercell%symmetries(i)%cartesian_rotation * x
      
      xx(atom_1p) = xx(atom_1p) + outer_product(x,x)
      
      do atom_2=1,supercell%no_atoms
        atom_2p = supercell%symmetries(i)%atom_group * atom_2
        
        f = supercell%symmetries(i)%cartesian_rotation * forces(atom_2,j)
        
        fx(atom_1p,atom_2p) = fx(atom_1p,atom_2p) + outer_product(f,x)
      enddo
    enddo
  enddo
  
  ! --------------------------------------------------
  ! Construct xx_inverse.
  ! --------------------------------------------------
  allocate(xx_inverse(supercell%no_atoms),stat=ialloc); call err(ialloc)
  do i=1,supercell%no_atoms
    xx_inverse(i) = invert(xx(i))
  enddo
  
  ! --------------------------------------------------
  ! Construct F.
  ! --------------------------------------------------
  allocate(force_constants(supercell%no_atoms,supercell%no_atoms), &
     & stat=ialloc); call err(ialloc)
  do i=1,supercell%no_atoms
    do j=1,supercell%no_atoms
      force_constants(j,i) = fx(j,i) * xx_inverse(i)
    enddo
  enddo
  
  ! --------------------------------------------------
  ! Symmetrise F.
  ! --------------------------------------------------
  do i=1,supercell%no_atoms
    do j=1,supercell%no_atoms
      force_constants(j,i) = ( force_constants(j,i)            &
                         & + transpose(force_constants(i,j)) &
                         & ) / 2.0_dp
      force_constants(i,j) = transpose(force_constants(i,j))
    enddo
  enddo
  
  ! --------------------------------------------------
  ! Average F across primitive lattice R-vectors, and
  !    convert to a mode-mode-Rvector representation.
  ! --------------------------------------------------
  allocate( output( supercell%sc_size,        &
          &         supercell%no_atoms_prim,  &
          &         supercell%no_atoms_prim), &
          & stat=ialloc); call err(ialloc)
  output = mat([ 0.0_dp,0.0_dp,0.0_dp, &
               & 0.0_dp,0.0_dp,0.0_dp, &
               & 0.0_dp,0.0_dp,0.0_dp  ], 3,3)
  
  rvector_group = supercell%calculate_rvector_group()
  
  do atom_1=1,supercell%no_atoms
    atom_1p = supercell%atom_to_prim(atom_1)
    rvector_1 = supercell%atom_to_rvec(atom_1)
    do atom_2=1,supercell%no_atoms
      atom_2p = supercell%atom_to_prim(atom_2)
      rvector_2 = supercell%atom_to_rvec(atom_2)
      
      rvector = rvector_group(supercell%paired_rvec(rvector_1)) * rvector_2
      
      output(rvector,atom_1p,atom_2p) = output(rvector,atom_1p,atom_2p) &
                                    & + force_constants(atom_1,atom_2)  &
                                    & / supercell%sc_size
    enddo
  enddo
end function
end module
