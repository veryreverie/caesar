! ======================================================================
! A class to help with calculating <bra|ket> integrals.
! ======================================================================
! Takes a bra and ket (and optionally a subspace and monomial),
!    and converts to a representation where all modes are represented,
!    where previously implicit modes are now explicit with power zero.
module state_helper_module
  use common_module
  implicit none
  
  type, extends(NoDefaultConstructor) :: StateHelper
    type(ComplexUnivariate), allocatable :: bra(:)
    type(ComplexUnivariate), allocatable :: ket(:)
    type(ComplexUnivariate), allocatable :: monomial(:)
    
    logical, allocatable :: bra_in_subspace(:)
    logical, allocatable :: ket_in_subspace(:)
    logical, allocatable :: monomial_in_subspace(:)
  end type
  
  interface StateHelper
    module procedure new_StateHelper_bra_ket
    module procedure new_StateHelper_bra_ket_monomial
    module procedure new_StateHelper_bra_ket_subspace
  end interface
  
  interface size
    module procedure size_StateHelper
  end interface
contains

! Constructors and size function.
function new_StateHelper_bra_ket(bra,ket) result(this)
  implicit none
  
  type(ComplexMonomial), intent(in) :: bra
  type(ComplexMonomial), intent(in) :: ket
  type(StateHelper)                 :: this
  
  integer :: i_bra,i_ket
  integer :: id_bra,id_ket
  integer :: n,ialloc
  
  n = size(bra)+size(ket)
  allocate( this%bra(n), &
          & this%ket(n), &
          & stat=ialloc); call err(ialloc)
  i_bra = 1
  i_ket = 1
  n = 0
  do while (i_bra<=size(bra) .or. i_ket<=size(ket))
    n = n+1
    
    if (i_bra>size(bra)) then
      id_bra = huge(0)
    else
      id_bra = bra%id(i_bra)
    endif
    
    if (i_ket>size(ket)) then
      id_ket = huge(0)
    else
      id_ket = ket%id(i_ket)
    endif
    
    if (id_bra==id_ket) then
      this%bra(n) = bra%mode(i_bra)
      this%ket(n) = ket%mode(i_ket)
      i_bra = i_bra + 1
      i_ket = i_ket + 1
    elseif (id_bra<id_ket) then
      this%bra(n) = bra%mode(i_bra)
      this%ket(n) = ComplexUnivariate(          &
         & id           = bra%id(i_bra),        &
         & paired_id    = bra%paired_id(i_bra), &
         & power        = 0,                    &
         & paired_power = 0                     )
      i_bra = i_bra + 1
    elseif (id_bra>id_ket) then
      this%bra(n) = ComplexUnivariate(          &
         & id           = ket%id(i_ket),        &
         & paired_id    = ket%paired_id(i_ket), &
         & power        = 0,                    &
         & paired_power = 0                     )
      this%ket(n) = ket%mode(i_ket)
      i_ket = i_ket + 1
    else
      call err()
    endif
  enddo
  
  this%bra = this%bra(:n)
  this%ket = this%ket(:n)
end function

function new_StateHelper_bra_ket_monomial(bra,ket,monomial) result(this)
  implicit none
  
  type(ComplexMonomial), intent(in) :: bra
  type(ComplexMonomial), intent(in) :: ket
  type(ComplexMonomial), intent(in) :: monomial
  type(StateHelper)                 :: this
  
  integer :: i_bra,i_ket,i_monomial
  integer :: id_bra,id_ket,id_monomial
  integer :: n,ialloc
  
  n = size(bra)+size(ket)+size(monomial)
  allocate( this%bra(n),      &
          & this%ket(n),      &
          & this%monomial(n), &
          & stat=ialloc); call err(ialloc)
  i_bra = 1
  i_ket = 1
  i_monomial = 1
  n = 0
  do while ( i_bra<=size(bra) .or.      &
           & i_ket<=size(ket) .or.      &
           & i_monomial<=size(monomial) )
    n = n+1
    
    if (i_bra>size(bra)) then
      id_bra = huge(0)
    else
      id_bra = bra%id(i_bra)
    endif
    
    if (i_ket>size(ket)) then
      id_ket = huge(0)
    else
      id_ket = ket%id(i_ket)
    endif
    
    if (i_monomial>size(monomial)) then
      id_monomial = huge(0)
    else
      id_monomial = monomial%id(i_monomial)
    endif
    
    if (id_bra<id_ket .and. id_bra<id_monomial) then
      this%bra(n) = bra%mode(i_bra)
      this%ket(n) = ComplexUnivariate(          &
         & id           = bra%id(i_bra),        &
         & paired_id    = bra%paired_id(i_bra), &
         & power        = 0,                    &
         & paired_power = 0                     )
      this%monomial(n) = ComplexUnivariate(     &
         & id           = bra%id(i_bra),        &
         & paired_id    = bra%paired_id(i_bra), &
         & power        = 0,                    &
         & paired_power = 0                     )
      i_bra = i_bra + 1
    elseif (id_ket<id_bra .and. id_ket<id_monomial) then
      this%bra(n) = ComplexUnivariate(          &
         & id           = ket%id(i_ket),        &
         & paired_id    = ket%paired_id(i_ket), &
         & power        = 0,                    &
         & paired_power = 0                     )
      this%ket(n) = ket%mode(i_ket)
      this%monomial(n) = ComplexUnivariate(     &
         & id           = ket%id(i_ket),        &
         & paired_id    = ket%paired_id(i_ket), &
         & power        = 0,                    &
         & paired_power = 0                     )
      i_ket = i_ket + 1
    elseif (id_monomial<id_bra .and. id_monomial<id_ket) then
      this%bra(n) = ComplexUnivariate(                    &
         & id           = monomial%id(i_monomial),        &
         & paired_id    = monomial%paired_id(i_monomial), &
         & power        = 0,                              &
         & paired_power = 0                               )
      this%ket(n) = ComplexUnivariate(                    &
         & id           = monomial%id(i_monomial),        &
         & paired_id    = monomial%paired_id(i_monomial), &
         & power        = 0,                              &
         & paired_power = 0                               )
      this%monomial(n) = monomial%mode(i_monomial)
      i_monomial = i_monomial + 1
    elseif (id_bra==id_ket .and. id_bra<id_monomial) then
      this%bra(n) = bra%mode(i_bra)
      this%ket(n) = ket%mode(i_ket)
      this%monomial(n) = ComplexUnivariate(     &
         & id           = bra%id(i_bra),        &
         & paired_id    = bra%paired_id(i_bra), &
         & power        = 0,                    &
         & paired_power = 0                     )
      i_bra = i_bra + 1
      i_ket = i_ket + 1
    elseif (id_bra==id_monomial .and. id_bra<id_ket) then
      this%bra(n) = bra%mode(i_bra)
      this%ket(n) = ComplexUnivariate(          &
         & id           = bra%id(i_bra),        &
         & paired_id    = bra%paired_id(i_bra), &
         & power        = 0,                    &
         & paired_power = 0                     )
      this%monomial(n) = monomial%mode(i_monomial)
      i_bra = i_bra + 1
      i_monomial = i_monomial + 1
    elseif (id_ket==id_monomial .and. id_ket<id_bra) then
      this%bra(n) = ComplexUnivariate(          &
         & id           = ket%id(i_bra),        &
         & paired_id    = ket%paired_id(i_bra), &
         & power        = 0,                    &
         & paired_power = 0                     )
      this%ket(n) = ket%mode(i_ket)
      this%monomial(n) = monomial%mode(i_monomial)
      i_ket = i_ket + 1
      i_monomial = i_monomial + 1
    elseif (id_bra==id_ket .and. id_bra==id_monomial) then
      this%bra(n) = bra%mode(i_bra)
      this%ket(n) = ket%mode(i_ket)
      this%monomial(n) = monomial%mode(i_monomial)
      i_bra = i_bra + 1
      i_ket = i_ket + 1
      i_monomial = i_monomial + 1
    else
      call err()
    endif
  enddo
  
  this%bra = this%bra(:n)
  this%ket = this%ket(:n)
  this%monomial = this%monomial(:n)
end function

function new_StateHelper_bra_ket_subspace(bra,ket,monomial,subspace) &
   & result(this)
  implicit none
  
  type(ComplexMonomial),    intent(in), optional :: bra
  type(ComplexMonomial),    intent(in), optional :: ket
  type(ComplexMonomial),    intent(in), optional :: monomial
  type(DegenerateSubspace), intent(in)           :: subspace
  type(StateHelper)                              :: this
  
  integer,                 allocatable :: unique_modes(:)
  type(ComplexUnivariate), allocatable :: bra_modes(:)
  type(ComplexUnivariate), allocatable :: ket_modes(:)
  type(ComplexUnivariate), allocatable :: monomial_modes(:)
  integer                              :: i_bra,i_ket,i_monomial,i_subspace
  integer                              :: id,paired_id
  type(ComplexUnivariate)              :: constant_mode
  
  integer :: i,ialloc
  
  unique_modes = filter(subspace%mode_ids<=subspace%paired_ids)
  
  if (present(bra)) then
    this%bra_in_subspace = [(                             &
       & any(bra%id(i)==subspace%mode_ids(unique_modes)), &
       & i=1,                                             &
       & size(bra)                                        )]
    bra_modes = bra%modes(filter(this%bra_in_subspace))
    allocate(this%bra(size(unique_modes)), stat=ialloc); call err(ialloc)
  endif
  
  if (present(ket)) then
    this%ket_in_subspace = [(                             &
       & any(ket%id(i)==subspace%mode_ids(unique_modes)), &
       & i=1,                                             &
       & size(ket)                                        )]
    ket_modes = ket%modes(filter(this%ket_in_subspace))
    allocate(this%ket(size(unique_modes)), stat=ialloc); call err(ialloc)
  endif
  
  if (present(monomial)) then
    this%monomial_in_subspace = [(                             &
       & any(monomial%id(i)==subspace%mode_ids(unique_modes)), &
       & i=1,                                                  &
       & size(monomial)                                        )]
    monomial_modes = monomial%modes(filter(this%monomial_in_subspace))
    allocate(this%monomial(size(unique_modes)), stat=ialloc); call err(ialloc)
  endif
  
  i_bra = 1
  i_ket = 1
  i_monomial = 1
  do i_subspace=1,size(unique_modes)
    id = subspace%mode_ids(unique_modes(i_subspace))
    paired_id = subspace%paired_ids(unique_modes(i_subspace))
    constant_mode = ComplexUnivariate(id, paired_id, 0, 0)
    
    if (present(bra)) then
      if (i_bra<=size(bra_modes)) then
        if (bra_modes(i_bra)%id==id) then
          this%bra(i_subspace) = bra_modes(i_bra)
          i_bra = i_bra + 1
        else
          this%bra(i_subspace) = constant_mode
        endif
      else
        this%bra(i_subspace) = constant_mode
      endif
    endif
    
    if (present(ket)) then
      if (i_ket<=size(ket_modes)) then
        if (ket_modes(i_ket)%id==id) then
          this%ket(i_subspace) = ket_modes(i_ket)
          i_ket = i_ket + 1
        else
          this%ket(i_subspace) = constant_mode
        endif
      else
        this%ket(i_subspace) = constant_mode
      endif
    endif
    
    if (present(monomial)) then
      if (i_monomial<=size(monomial_modes)) then
        if (monomial_modes(i_monomial)%id==id) then
          this%monomial(i_subspace) = monomial_modes(i_monomial)
          i_monomial = i_monomial + 1
        else
          this%monomial(i_subspace) = constant_mode
        endif
      else
        this%monomial(i_subspace) = constant_mode
      endif
    endif
  enddo
end function

function size_StateHelper(this) result(output)
  implicit none
  
  type(StateHelper), intent(in) :: this
  integer                       :: output
  
  output = size(this%bra)
end function
end module
