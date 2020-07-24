! ======================================================================
! Calculates the contribution to a set of dynamical matrices from
!    a given monomial or polynomial.
! ======================================================================
module polynomial_dynamical_matrices_module
  use common_module
  
  use anharmonic_common_module
  implicit none
  
  private
  
  public :: calculate_dynamical_matrices
  public :: calculate_correction
  
  interface calculate_dynamical_matrices
    module procedure calculate_dynamical_matrices_ComplexMonomial
    module procedure calculate_dynamical_matrices_ComplexPolynomial
  end interface
  
  interface calculate_correction
    module procedure calculate_correction_ComplexMonomial
    module procedure calculate_correction_ComplexPolynomial
  end interface
contains

! Calculate the monomial or polynomial's contribution to the effective
!    dynamical matrix from which the potential can be interpolated in the
!    large-supercell limit.
function calculate_dynamical_matrices_ComplexMonomial(term,qpoints, &
   & thermal_energy,subspaces,subspace_bases,subspace_states,       &
   & subspaces_in_coupling,anharmonic_data) result(output) 
  implicit none
  
  type(ComplexMonomial),    intent(in)    :: term
  type(QpointData),         intent(in)    :: qpoints(:)
  real(dp),                 intent(in)    :: thermal_energy
  type(DegenerateSubspace), intent(in)    :: subspaces(:)
  class(SubspaceBasis),     intent(in)    :: subspace_bases(:)
  class(BasisStates),       intent(inout) :: subspace_states(:)
  integer,                  intent(in)    :: subspaces_in_coupling(:)
  type(AnharmonicData),     intent(in)    :: anharmonic_data
  type(DynamicalMatrix), allocatable      :: output(:)
  
  integer                              :: prefactor
  type(ComplexUnivariate), allocatable :: modes(:)
  type(ComplexUnivariate), allocatable :: univariates(:)
  type(ComplexMode),       allocatable :: subspace_modes(:)
  type(ComplexMonomial)                :: monomial
  
  real(dp), allocatable :: expectations(:)
  real(dp)              :: coefficient
  
  type(ComplexMode) :: mode_k
  type(ComplexMode) :: mode_l
  
  integer :: q
  
  integer :: i,j,k,l,ialloc
  
  allocate( expectations(size(subspaces_in_coupling)), &
          & stat=ialloc); call err(ialloc)
  output = [( DynamicalMatrix(anharmonic_data%structure%no_atoms), &
            & i=1,                                                 &
            & size(qpoints)                                        )]
  prefactor = (term%total_power()*(term%total_power()-1))/2
  modes = term%modes()
  
  ! Calculate the expectation of the part of term in each subspace.
  do i=1,size(subspaces_in_coupling)
    j = subspaces_in_coupling(i)
    univariates = modes(filter(modes%id .in. subspaces(j)%mode_ids))
    
    monomial = ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp), &
                              & modes       = univariates     )
    call integrate( monomial,           &
                  & subspace_states(j), &
                  & subspaces(j),       &
                  & subspace_bases(j),  &
                  & anharmonic_data     )
    expectations(i) = real(monomial%coefficient)
  enddo
  
  ! Calculate the contribution to the dynamical matrices
  !    from each subspace.
  ! If a term is prod_{i=1}^n [u_i], then the u_i u_j* coefficient is
  !   binom(n,2) * <term / (u_i u_j*)>.
  do i=1,size(subspaces_in_coupling)
    j = subspaces_in_coupling(i)
    univariates = modes(filter(     &
       & modes%id .in. subspaces(j)%mode_ids ))
    subspace_modes = anharmonic_data%complex_modes([( &
         & first( univariates(k)%id          &
         &     == anharmonic_data%complex_modes%id ), &
         & k=1,                                       &
         & size(univariates)                 )])
    
    do k=1,size(univariates)
      ! If (u_k)^{n_k} has n_k==0, a factor of u_k cannot be removed.
      if (univariates(k)%power==0) then
        cycle
      endif
      
      mode_k = subspace_modes(k)
      
      q = first(qpoints%id==mode_k%qpoint_id)
      do l=1,size(univariates)
        mode_l = subspace_modes(l)
        
        ! If u_k = u_k*, then only loop over l<=k.
        if (univariates(k)%id == univariates(k)%paired_id) then
          if (l>k) then
            exit
          endif
        endif
        
        ! If q_l/=q_k, <u_k u_l*> must be zero by translational symmetry.
        ! If (u_l*)^{n_l} has n_l==0, a factor of u_l cannot be removed.
        if (mode_l%qpoint_id/=mode_k%qpoint_id) then
          cycle
        elseif (univariates(l)%paired_power==0) then
          cycle
        endif
        
        ! Construct the monomial corresponding to the part of the term in
        !    subspace j, but with (u_k u_l*) removed.
        if (univariates(k)%id==univariates(k)%paired_id) then
          if (k==l .and. univariates(k)%power<2) then
            cycle
          endif
          univariates(k)%power = univariates(k)%power-1
          univariates(k)%paired_power = univariates(k)%paired_power-1
          univariates(l)%power = univariates(l)%power-1
          univariates(l)%paired_power = univariates(l)%paired_power-1
          monomial = ComplexMonomial(                 &
             & coefficient = cmplx(1.0_dp,0.0_dp,dp), &
             & modes       = univariates     )
          univariates(k)%power = univariates(k)%power+1
          univariates(k)%paired_power = univariates(k)%paired_power+1
          univariates(l)%power = univariates(l)%power+1
          univariates(l)%paired_power = univariates(l)%paired_power+1
        else
          univariates(k)%power = univariates(k)%power-1
          univariates(l)%paired_power = univariates(l)%paired_power-1
          monomial = ComplexMonomial(                 &
             & coefficient = cmplx(1.0_dp,0.0_dp,dp), &
             & modes       = univariates           )
          univariates(k)%power = univariates(k)%power+1
          univariates(l)%paired_power = univariates(l)%paired_power+1
        endif
        
        ! Integrate this monomial, and multiply the result by the integrated
        !    parts of the term which are not in subspace j.
        ! This gives <term / u_ku_l*>.
        call integrate( monomial,           &
                      & subspace_states(j), &
                      & subspaces(j),       &
                      & subspace_bases(j),  &
                      & anharmonic_data     )
        coefficient = product(expectations(:i-1)) &
                  & * real(monomial%coefficient)  &
                  & * product(expectations(i+1:)) &
                  & * prefactor
        
        ! Construct the u_k u_l* matrix, multiply by the coefficient,
        !    and add it to the relevant q-point's dynamical matrix.
        output(q) = output(q) + DynamicalMatrix(mode_k,mode_l,coefficient)
        
        ! If u_k/=u_k* then the u_l u_k* term is just the conjugate of the
        !    u_k u_l* term.
        if (mode_l%id/=mode_k%id) then
          output(q) = output(q) + DynamicalMatrix(mode_l,mode_k,coefficient)
        endif
      enddo
    enddo
  enddo
  
  ! Multiply by the term's coefficient, then convert from a term coefficient to
  !    a dynamical matrix.
  ! The dynamical matrix elements are -2/N times the term coefficients.
  output = output           &
       & * term%coefficient &
       & / (-0.5_dp*anharmonic_data%anharmonic_supercell%sc_size)
end function

function calculate_dynamical_matrices_ComplexPolynomial(polynomial,qpoints, &
   & thermal_energy,subspaces,subspace_bases,subspace_states,               &
   & subspaces_in_coupling,anharmonic_data) result(output) 
  implicit none
  
  type(ComplexPolynomial),  intent(in)    :: polynomial
  type(QpointData),         intent(in)    :: qpoints(:)
  real(dp),                 intent(in)    :: thermal_energy
  type(DegenerateSubspace), intent(in)    :: subspaces(:)
  class(SubspaceBasis),     intent(in)    :: subspace_bases(:)
  class(BasisStates),       intent(inout) :: subspace_states(:)
  integer,                  intent(in)    :: subspaces_in_coupling(:)
  type(AnharmonicData),     intent(in)    :: anharmonic_data
  type(DynamicalMatrix), allocatable      :: output(:)
  
  integer :: i
  
  output = [( DynamicalMatrix(anharmonic_data%structure%no_atoms), &
            & i=1,                                                 &
            & size(qpoints)                                        )]
  do i=1,size(polynomial%terms)
    output = output                        &
         & + calculate_dynamical_matrices( &
         &          polynomial%terms(i),   &
         &          qpoints,               &
         &          thermal_energy,        &
         &          subspaces,             &
         &          subspace_bases,        &
         &          subspace_states,       &
         &          subspaces_in_coupling, &
         &          anharmonic_data        )
  enddo
end function

! Calculate the correction due to double-counting when interpolating
!    polynomials.
! It is assumed that a monomial of order n will appear in:
!    n/2     terms if n is even
!    (n-1)/2 terms if n is odd.
! The double counting of a term is (1-x)<term>, where x is the number of terms
!    where the term appears.
function calculate_correction_ComplexMonomial(monomial,subspaces, &
   & subspace_bases,subspace_states,anharmonic_data) result(output) 
  implicit none
  
  class(ComplexMonomial),   intent(in)    :: monomial
  type(DegenerateSubspace), intent(in)    :: subspaces(:)
  class(SubspaceBasis),     intent(in)    :: subspace_bases(:)
  class(BasisStates),       intent(inout) :: subspace_states(:)
  type(AnharmonicData),     intent(in)    :: anharmonic_data
  real(dp)                                :: output
  
  type(ComplexMonomial) :: term
  
  integer :: i
  
  ! Calculate <term>.
  term = monomial
  do i=1,size(subspaces)
    call integrate( term,               &
                  & subspace_states(i), &
                  & subspaces(i),       &
                  & subspace_bases(i),  &
                  & anharmonic_data     )
  enddo
  
  ! N.B. this intentionally uses integer floor division.
  output = real(term%coefficient * (1-monomial%total_power()/2))
end function

function calculate_correction_ComplexPolynomial(polynomial,subspaces, &
   & subspace_bases,subspace_states,anharmonic_data) result(output) 
  implicit none
  
  class(ComplexPolynomial), intent(in)    :: polynomial
  type(DegenerateSubspace), intent(in)    :: subspaces(:)
  class(SubspaceBasis),     intent(in)    :: subspace_bases(:)
  class(BasisStates),       intent(inout) :: subspace_states(:)
  type(AnharmonicData),     intent(in)    :: anharmonic_data
  real(dp)                                :: output
  
  integer :: i
  
  output = 0
  do i=1,size(polynomial%terms)
    output = output + calculate_correction( polynomial%terms(i), &
                                          & subspaces,           &
                                          & subspace_bases,      &
                                          & subspace_states,     &
                                          & anharmonic_data      )
  enddo
end function
end module
