! ---------------------------------------------------------------------------
! Calculates Regular and Moller-Plesset perturbations to energies and states.
! ---------------------------------------------------------------------------
!
! let x be some perturbational parameter
!
! The unperturbed Hamiltonian, h, has states {|p>} and energies {e}
! The total Hamiltonian, H, has states {|P>} and energies {E}
!
! For brevity, ei=(e)i, and |i>=(|p>)i, where {i} label states
!
! H|P> = E|P>
!
! H = h + x.v
! E = e + x.e1 + x^2.e2 + x^3.e3 + x^4.e4 + ... 0(x^5)
! |P> = |p> + x.|p1> + x^2.|p2> + x^3.|p3> + x^4.|p4> + ... O(x^5)
! 
! => 0 =  h|p>          - e|p>
!      + (h|p1> + v|p>  - e|p1> - e1|p>).x
!      + (h|p2> + v|p1> - e|p2> - e1|p1> - e2|p>).x^2
!      + (h|p3> + v|p2> - e|p3> - e1|p2> - e2|p1> - e3|p>).x^3
!      + ...
!
! For convenience, the powers of x are dropped from here on.
module perturbation_theory
  use constants, only : dp
  use linear_algebra, only : RealEigenstuff, calculate_eigenstuff, size
  implicit none
  
  ! Moller-Plesset corrections
  public :: mp1_energy ! first order energy correction
  public :: mp1_state  ! first order state correction
  public :: mp2_energy ! second order energy correction
  public :: mp2_state  ! second order state correction
  public :: mp3_energy ! third order energy correction
  
  ! Helper functions
  private :: calculate_d ! d(i,j) = 1/(ei-ej)
  private :: calculate_m ! m(i,j) = <i|v|j>

contains

! --------------------------------------------------
! Calculates d, d(i,j) = 1/(ei-ej)
! if i and j are degenerate, d(i,j) is set to 1
! (in such cases, m(i,j) is set to 0)
! --------------------------------------------------
function calculate_d(estuff) result(d)
  implicit none
  
  type(RealEigenstuff), intent(in) :: estuff ! unperturbed evals and evecs
  real(dp), allocatable            :: d(:,:) ! denominators, d(i,j) = 1/(ei-ej)
  
  integer :: i,j ! loop indices
  
  allocate(d(size(estuff),size(estuff)))
  
  do i=1,size(estuff)
    do j=i+1,size(estuff)
      if (allocated(estuff%degeneracy)                  .and. &
        & estuff%degeneracy(i) /= 0                     .and. &
        & estuff%degeneracy(i) == estuff%degeneracy(j)) then
        ! i and j are degenerate with one another. prevent NaNs.
        d(i,j) = 1
      else
        ! i and j are not degenerate with one another
        d(i,j) = 1/(estuff%evals(i)-estuff%evals(j))
      endif
      d(j,i) = - d(i,j)
    enddo
  enddo
end function

! --------------------------------------------------
! Calculates m, m(i,j) = <i|v|j>
! if i and j are degenerate with one another, m(i,j) = 0 to prevent mixing
! --------------------------------------------------
function calculate_m(estuff,v) result(m)
  implicit none
  
  type(RealEigenstuff), intent(in) :: estuff ! unperturbed evals and evecs
  real(dp),             intent(in) :: v(:,:) ! perturbation
  real(dp), allocatable            :: m(:,:) ! matrix elements, m(i,j)=<i|v|j>
  
  integer :: i,j ! loop indices
  
  allocate(m(size(estuff),size(estuff)))
  
  ! calculate all elements
  m = matmul(matmul(estuff%evecs,v),estuff%evecs)
  
  ! zero degenerate elements to prevent mixing
  do i=1,size(estuff)
    do j=i+1,size(estuff)
      if (allocated(estuff%degeneracy)                  .and. &
        & estuff%degeneracy(i) /= 0                     .and. &
        & estuff%degeneracy(i) == estuff%degeneracy(j)) then
        m(i,j) = 0
        m(j,i) = 0
      endif
    enddo
  enddo
  
end function

! --------------------------------------------------
! Calculates non-degenerate basis
! --------------------------------------------------
function lift_degeneracy(estuff,v,threshold) result(output)
  implicit none
  
  type(RealEigenstuff), intent(in) :: estuff    ! unperturbed evals and evecs
  real(dp),             intent(in) :: v(:,:)    ! perturbation
  real(dp),             intent(in) :: threshold ! degeneracy threshold
  type(RealEigenstuff)             :: output    ! output evals and evecs
  
  integer :: no_states ! number of states. = size(estuff%evals)
  
  ! degeneracy lists
  integer, allocatable :: lists(:,:)    ! lists of degenerate eigenstates
  integer, allocatable :: lengths(:)    ! size(lists(i,:))
  integer, allocatable :: which_list(:) ! the list in which ei is a member
  
  ! degenerate objects
  real(dp), allocatable     :: m(:,:) ! m(i,j) = <i|v|j>
  real(dp), allocatable     :: h(:,:) ! degenerate hamiltonian
  type(RealEigenstuff)      :: hstuff ! evals and evecs of h
  
  ! temporary variables
  integer :: i,j,k ! loop indices
  integer :: list_id
  integer :: list_j_id
  
  no_states = size(estuff)
  
  ! allocate and initialise output
  allocate(output%evals(no_states))
  output%evals = 0.d0
  allocate(output%evecs(no_states,no_states))
  output%evecs = 0.d0
  allocate(output%degeneracy(no_states))
  output%degeneracy = 0
  
  ! allocate and initialise lists, lengths and which_list
  allocate(lists(no_states,no_states))
  lists = 0
  allocate(lengths(no_states))
  lengths = 0
  allocate(which_list(no_states))
  which_list = 0
  
  ! identify degeneracies
  do i=1,no_states
    do j=i+1,no_states
      if (abs(estuff%evals(i)-estuff%evals(j)) < threshold) then
        ! check if either eigenstate is already in a list
        if (which_list(i) /= 0 .and. which_list(j) /= 0) then
          ! move list containing j to list containing i
          list_id = which_list(i)
          list_j_id = which_list(j)
          do k=1,lengths(list_j_id)
            lists(list_id,lengths(list_id)+k) = lists(list_j_id,k)
            which_list(lists(list_j_id,k)) = list_id
            lists(list_j_id,k) = 0
          enddo
          lengths(list_id) = lengths(list_id)+lengths(which_list(j))
          lengths(which_list(j)) = 0
        elseif(which_list(i) /= 0) then
          ! add j to list containing i
          list_id = which_list(i)
          lists(list_id,lengths(list_id)+1) = j
          lengths(list_id) = lengths(list_id)+1
          which_list(j) = list_id
        elseif(which_list(j) /= 0) then
          ! add i to list containing j
          list_id = which_list(j)
          lists(list_id,lengths(list_id)+1) = i
          lengths(list_id) = lengths(list_id)+1
          which_list(i) = list_id
        else
          lists(i,1) = i
          lists(i,2) = j
          lengths(i) = 2
          which_list(i) = i
          which_list(j) = i
        endif
      endif
    enddo
  enddo
  
  ! copy accross non-degenerate states
  do i=1,no_states
    if (which_list(i) == 0) then
      output%evals(i) = estuff%evals(i)
      output%evecs(i,:) = estuff%evecs(i,:)
      output%degeneracy(i) = 0
    endif
  enddo
  
  ! calculate m = <i|v|j>
  m = calculate_m(estuff,v)
  
  ! calculate non-degenerate states
  do i=1,no_states
    if (lengths(i) == 0) cycle
    
    ! assemble degenerate hamiltonian
    allocate(h(lengths(i),lengths(i)))
    h = 0.d0
    do j=1,lengths(i)
      h(j,j) = estuff%evals(lists(i,j)) ! add in unperturbed elements
      do k=1,lengths(i)
        h(j,k) = h(j,k) + m(lists(i,j),lists(i,k))
      enddo
    enddo
    
    ! diagonalise h
    hstuff = calculate_eigenstuff(h)
    
    ! copy across new states
    do j=1,lengths(i)
      output%evals(lists(i,j)) = hstuff%evals(j)
      do k=1,lengths(i)
        output%evecs(lists(i,j),:) = output%evecs(lists(i,j),:) &
                                 & + hstuff%evecs(j,k) &
                                 & * estuff%evecs(lists(i,k),:)
      enddo
      output%degeneracy(lists(i,j)) = i
    enddo
    
    deallocate(hstuff%evals)
    deallocate(hstuff%evecs)
    deallocate(h)
  enddo
  
  ! deallocate arrays
  deallocate(lists)
  deallocate(lengths)
  deallocate(which_list)
end function

! --------------------------------------------------
! Calculates the first-order energy correction
! e1 = <p|v|p>
! --------------------------------------------------
function mp1_energy(estuff,v) result(e1)
  implicit none
  
  type(RealEigenstuff), intent(in) :: estuff ! unperturbed eigenvals and eigenvecs
  real(dp),             intent(in) :: v(:,:) ! perturbation
  real(dp), allocatable            :: e1(:)
  
  integer :: i ! loop index
  
  ! allocate output
  allocate(e1(size(estuff)))
  
  do i=1,size(estuff)
    if (estuff%degeneracy(i) /= 0) then
      ! the state is degenerate; no correction needed
      e1(i) = 0
    else
      ! the state is not degenerate; calculate correction
      e1(i) = dot_product(matmul(estuff%evecs(i,:),v),estuff%evecs(i,:))
    endif
  enddo
end function


! --------------------------------------------------
! Calculates the first-order state correction
! |p1>i = sum(j/=i) <j|v|i>/(ei-ej) = sum(j/=i) mij*dij
! --------------------------------------------------
function mp1_state(estuff,v) result(p1)
  implicit none
  
  type(RealEigenstuff), intent(in) :: estuff ! unperturbed evals and evecs
  real(dp),             intent(in) :: v(:,:) ! perturbation, v
  real(dp), allocatable            :: p1(:,:)
  
  real(dp), allocatable :: m(:,:) ! m(i,j) = <i|v|j>
  real(dp), allocatable :: d(:,:) ! d(i,j) = 1/(ei-ej)
  
  integer :: i,j ! loop indices
  
  ! allocate p1
  allocate(p1(size(estuff),size(estuff)))
  
  d = calculate_d(estuff)
  m = calculate_m(estuff,v)
  
  do i=1,size(estuff)
    do j=i+1,size(estuff)
      p1(i,:) = p1(i,:) + m(i,j)*d(i,j)*estuff%evecs(j,:)
      p1(j,:) = p1(j,:) - m(i,j)*d(i,j)*estuff%evecs(i,:)
    enddo
  enddo
  
  deallocate(d)
  deallocate(m)
end function

! --------------------------------------------------
! calculates second-order energy correction
! (e2)i = sum(j/=i) |<j|v|i>|**2/(ei-ej)
! --------------------------------------------------
function mp2_energy(estuff,v) result(e2)
  implicit none
  
  type(RealEigenstuff), intent(in) :: estuff ! unperturbed evals and evecs
  real(dp),             intent(in) :: v(:,:) ! perturbation
  real(dp), allocatable            :: e2(:)
  
  real(dp), allocatable :: m(:,:) ! m(i,j) = <i|v|j>
  real(dp), allocatable :: d(:,:) ! d(i,j) = 1/(ei-ej)
  
  integer :: i,j ! loop indices
  
  allocate(e2(size(estuff)))
  e2 = 0.d0
  
  d = calculate_d(estuff)
  m = calculate_m(estuff,v)
  
  do i=1,size(estuff)
    do j=i+1,size(estuff)
      e2(i) = e2(i) + m(i,j)*m(i,j)*d(i,j)
      e2(j) = e2(j) - m(i,j)*m(i,j)*d(i,j)
    enddo
  enddo
  
  deallocate(d)
  deallocate(m)
end function

! --------------------------------------------------
! calculates second-order state correction
! --------------------------------------------------
function mp2_state(estuff,v) result(p2)
  implicit none
  
  type(RealEigenstuff), intent(in) :: estuff  ! unperturbed evals and evecs
  real(dp),             intent(in) :: v(:,:)  ! perturbation
  real(dp), allocatable            :: p2(:,:)
  
  real(dp), allocatable :: m(:,:) ! m(i,j) = <i|v|j>
  real(dp), allocatable :: d(:,:) ! d(i,j) = 1/(ei-ej)
  
  integer :: i,j,k ! loop indices
  
  allocate(p2(size(estuff),size(estuff)))
  p2 = 0.d0
  
  d = calculate_d(estuff)
  m = calculate_m(estuff,v)
  
  ! three-state term
  do i=1,size(estuff)
    do j=1,size(estuff)
      if (j==i) cycle
      
      ! three-state term
      do k=1,size(estuff)
        if (k==i) cycle
        p2(i,:) = p2(i,:) + m(k,j)*m(j,i)*d(i,k)*d(i,j)*estuff%evecs(k,:)
      enddo
      
      ! two-state term
      p2(i,:) = p2(i,:) - m(i,i)*m(j,i)*d(i,j)*d(i,j)*estuff%evecs(j,:)
      
      ! self-mixing term
      p2(i,:) = p2(i,:) - 0.5d0*m(i,j)*m(i,j)*d(i,j)*d(i,j)*estuff%evecs(i,:)
    enddo
  enddo
  
  deallocate(d)
  deallocate(m)
end function

! --------------------------------------------------
! calculates third-order energy correction
! --------------------------------------------------
function mp3_energy(estuff,v) result(e3)
  implicit none
  
  type(RealEigenstuff), intent(in) :: estuff ! unperturbed evals and evecs
  real(dp),             intent(in) :: v(:,:) ! perturbation
  real(dp), allocatable            :: e3(:)
  
  real(dp), allocatable :: m(:,:) ! m(i,j) = <i|v|j>
  real(dp), allocatable :: d(:,:) ! d(i,j) = 1/(ei-ej)
  
  integer :: i,j,k     ! loop indices
  
  allocate(e3(size(estuff)))
  e3 = 0.d0
  
  d = calculate_d(estuff)
  m = calculate_m(estuff,v)
  
  ! three-state term
  do i=1,size(estuff)
    do j=1,size(estuff)
      if (j==i) cycle
      
      ! three-state term
      do k=1,size(estuff)
        if (k==i) cycle
        e3(i) = e3(i) + m(i,j)*m(j,k)*m(k,i)*d(j,i)*d(k,i)
      enddo
      
      ! two-state term
      e3(i) = e3(i) - m(i,i)*m(i,j)*m(i,j)*d(i,j)*d(i,j)
      
    enddo
  enddo
  
  deallocate(d)
  deallocate(m)
end function

end module
