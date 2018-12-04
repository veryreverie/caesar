! ======================================================================
! Calculates the VSCF correction to the thermodynamic data.
! As usual, the output is normalised to be per primitive cell.
! ======================================================================
module vscf_thermodynamics_module
  use common_module
  
  use states_module
  use anharmonic_common_module
  use potentials_module
  implicit none
  
  private
  
  public :: EnergySpectrum
  public :: calculate_vscha_thermodynamics
  public :: calculate_vscf_spectrum
  public :: VscfThermodynamics
  !public :: calculate_vscf_thermodynamics
  
  type, extends(NoDefaultConstructor) :: EnergySpectrum
    real(dp), allocatable :: vscf_energies(:)
    real(dp), allocatable :: vscha_energies(:)
    integer,  allocatable :: vscha_occupations(:)
  end type
  
  ! Return type.
  type, extends(NoDefaultConstructor) :: VscfThermodynamics
    type(ThermodynamicData) :: vscha
    type(ThermodynamicData) :: vscf
  end type
  
  interface EnergySpectrum
    module procedure new_EnergySpectrum
  end interface
contains

function new_EnergySpectrum(vscf_energies,vscha_energies,vscha_occupations) &
   & result(this)
  implicit none
  
  real(dp), allocatable :: vscf_energies(:)
  real(dp), allocatable :: vscha_energies(:)
  integer,  allocatable :: vscha_occupations(:)
  type(EnergySpectrum)  :: this
  
  this%vscf_energies = vscf_energies
  this%vscha_energies = vscha_energies
  this%vscha_occupations = vscha_occupations
end function

! Calculates VSCHA thermodynamics.
function calculate_vscha_thermodynamics(vscha_frequency,potential,subspace, &
   & anharmonic_data,no_vscha_basis_states,thermal_energy) result(output)
  implicit none
  
  real(dp),                 intent(in) :: vscha_frequency
  class(PotentialData),     intent(in) :: potential
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  integer,                  intent(in) :: no_vscha_basis_states
  real(dp),                 intent(in) :: thermal_energy
  type(ThermodynamicData)              :: output(2)
  
  real(dp) :: harmonic_potential_expectation
  real(dp) :: anharmonic_potential_expectation
  
  output(1) = ThermodynamicData(thermal_energy, vscha_frequency) &
          & * size(subspace)                                     &
          & / real(anharmonic_data%anharmonic_supercell%sc_size,dp)
  harmonic_potential_expectation = output(1)%energy / 2
  anharmonic_potential_expectation = potential%harmonic_expectation( &
                                            & vscha_frequency,       &
                                            & thermal_energy,        &
                                            & no_vscha_basis_states, &
                                            & subspace,              &
                                            & anharmonic_data        )
  output(2) = output(1)
  output(2)%energy = output(2)%energy               &
                 & - harmonic_potential_expectation &
                 & + anharmonic_potential_expectation
  output(2)%free_energy = output(2)%free_energy          &
                      & - harmonic_potential_expectation &
                      & + anharmonic_potential_expectation
end function

! Calculates the energy spectrum of the VSCF states.
function calculate_vscf_spectrum(vscha_frequency,potential,subspace, &
   & anharmonic_data,no_basis_states) result(output)
  implicit none
  
  real(dp),                 intent(in) :: vscha_frequency
  class(PotentialData),     intent(in) :: potential
  type(DegenerateSubspace), intent(in) :: subspace
  type(AnharmonicData),     intent(in) :: anharmonic_data
  integer,                  intent(in) :: no_basis_states
  type(EnergySpectrum)                 :: output
  
  ! Variables for constructing the basis.
  type(StructureData) :: supercell
  type(SubspaceBasis) :: basis
  
  ! Variables for calculating and diagonalising the VSCF Hamiltonian.
  integer                                :: no_states
  type(HarmonicState)                    :: bra
  type(HarmonicState)                    :: ket
  real(dp),                  allocatable :: hamiltonian(:,:)
  type(SymmetricEigenstuff), allocatable :: estuff(:)
  
  real(dp), allocatable :: vscf_energies(:)
  real(dp), allocatable :: vscha_energies(:)
  integer,  allocatable :: vscha_occupations(:)
  
  integer :: i,j,k,l,ialloc
  
  supercell = anharmonic_data%anharmonic_supercell
  
  basis = SubspaceBasis( subspace,                                  &
                       & vscha_frequency,                           &
                       & anharmonic_data%complex_modes,             &
                       & anharmonic_data%qpoints,                   &
                       & anharmonic_data%anharmonic_supercell,      &
                       & no_basis_states-1,                         &
                       & anharmonic_data%potential_expansion_order, &
                       & anharmonic_data%structure%symmetries       )
  
  vscf_energies = [real::]
  vscha_energies = [real::]
  vscha_occupations = [integer::]
  do i=1,size(basis)
    no_states = size(basis%wavevectors(i))
    allocate( hamiltonian(no_states,no_states), &
            & stat=ialloc); call err(ialloc)
    hamiltonian = 0
    do j=1,no_states
      bra = basis%wavevectors(i)%harmonic_states(j)
      do k=1,size(basis%wavevectors(i)%harmonic_couplings(j))
        l = basis%wavevectors(i)%harmonic_couplings(j)%id(k)
        ket = basis%wavevectors(i)%harmonic_states(l)
        
        hamiltonian(j,l) = kinetic_energy(bra, ket, subspace, supercell) &
                       & + potential%potential_energy(bra,ket,anharmonic_data)
        if (j==l) then
          vscha_energies = [vscha_energies, hamiltonian(j,l)]
          vscha_occupations = [vscha_occupations, bra%total_occupation()]
        endif
      enddo
    enddo
    
    estuff = diagonalise_symmetric(hamiltonian)
    do j=1,basis%wavevectors(i)%degeneracy
      vscf_energies = [vscf_energies, estuff%eval]
    enddo
    deallocate(hamiltonian, stat=ialloc); call err(ialloc)
  enddo
  
  vscf_energies = vscf_energies * supercell%sc_size
  vscha_energies = vscha_energies * supercell%sc_size
  
  output = EnergySpectrum(vscf_energies, vscha_energies, vscha_occupations)
end function

!function calculate_vscf_thermodynamics(vscf_frequency,potential,subspace, &
!   & anharmonic_data,no_vscha_basis_states,thermal_energy,vscf_spectrum)  &
!   & result(output)
!  implicit none
!  
!  real(dp)                 :: vscf_frequency
!  class(PotentialData)     :: potential
!  type(DegenerateSubspace) :: subspace
!  type(AnharmonicData)     :: anharmonic_data
!  integer                  :: no_vscha_basis_states
!  real(dp)                 :: thermal_energy
!  type(EnergySpectrum)     :: vscf_spectrum
!  type(ThermodynamicData)  :: output
!  
!  type(ThermodynamicData) :: vscha_thermodynamics(2)
!  
!  real(dp) :: vscha_zpe
!  real(dp) :: vscf_zpe
!  
!  real(dp) :: z1
!  real(dp) :: dz1
!  real(dp) :: z2
!  real(dp) :: dz2
!  
!  real(dp), allocatable :: boltzmann_factors(:)
!  
!  vscha_zpe = minval(vscf_spectrum%vscha_energies)
!  vscf_zpe = minval(vscf_spectrum%vscf_energies)
!  
!  if (vscf_frequency > 690*thermal_energy) then
!    ! Very low-temperature regime. exp(-690)<1e-300, below dp floating point.
!    energy = vscf_zpe
!    free_energy = vscf_zpe
!    entropy = 0
!  else
!    vscha_thermodynamics = calculate_vscha_thermodynamics( &
!                                  & vscf_frequency,        &
!                                  & potential,             &
!                                  & subspace,              &
!                                  & anharmonic_data,       &
!                                  & no_vscha_basis_states, &
!                                  & thermal_energy         )
!    
!    ! z is the partition function.
!    ! dz is -dz/d(1/T), the derivative of z w/r/t 1/T.
!    
!    ! Calculate z and dz for the VSCF potential with VSCHA states
!    !    and weightings.
!    !    F = zpe - T*ln(z)
!    ! -> z = ln((zpe-F)/T)
!    !
!    !     U = zpe + dz/z
!    ! -> dz = z(U-zpe)
!    z1 = exp((vscha_zpe-vscha_thermodynamics(2)%free_energy)/thermal_energy)
!    dz1 = z1*(vscha_thermodynamics(2)%energy-vscha_zpe)
!    
!    ! p_i = exp((zpe-E_i)/T)
!    ! z = sum_i p_i
!    ! dz = sum_i (E_i-zpe) p_i
!    boltzmann_factors = [                                           &
!       exp((vscha_zpe-vscf_spectrum%vscha_energies)/thermal_energy) ]
!    z2 = sum(boltzmann_factors)
!    dz2 = sum((vscf_spectrum%vscha_energies-vscha_zpe)*boltzmann_factors)
!  enddo
!  
!  !type(ThermodynamicData) :: thermodynamic_data(4)
!  !
!  !real(dp) :: partition_functions(4)
!  !real(dp) :: eps(4)
!  !
!  !real(dp) :: energy
!  !real(dp) :: free_energy
!  !real(dp) :: entropy
!  !
!  !integer :: i
!  !
!  !thermodynamic_data(1:2) = calculate_vscha_thermodynamics( &
!  !                                 & vscf_frequency,        &
!  !                                 & potential,             &
!  !                                 & subspace,              &
!  !                                 & anharmonic_data,       &
!  !                                 & no_vscha_basis_states, &
!  !                                 & thermal_energy         )
!  !
!  !thermodynamic_data(3) = ThermodynamicData( thermal_energy,           &
!  !                                         & vscf_spectrum(1)%energies )
!  !thermodynamic_data(4) = ThermodynamicData( thermal_energy,           &
!  !                                         & vscf_spectrum(2)%energies )
!  !do i=1,4
!  !  ! TODO: zero_energy(i) = ...
!  !  partitions_functions(i) = exp( -thermodynamic_data(i)%free_energy &
!  !                             & /  thermal_energy                    )
!  !  eps(i) = thermodynamic_data(i)%energy * partition_functions(i)
!  !enddo
!  !energy = 
!  !free_energy = -thermal_energy*log(                                  &
!  !   &   exp(-vscha_thermodynamic_data(2)%free_energy/thermal_energy) &
!  !   & + exp(-vscf_thermodynamic_data(1)%free_energy/thermal_energy)  &
!  !   & - exp(-vscf_thermodynamic_data(2)%free_energy/thermal_energy)  )
!end function

!function calculate_vscf_thermodynamics(vscha_frequency,potential,subspace, &
!   & anharmonic_data,thermal_energy,no_basis_states) result(output)
!  implicit none
!  
!  real(dp),                 intent(in) :: vscha_frequency
!  class(PotentialData),     intent(in) :: potential
!  type(DegenerateSubspace), intent(in) :: subspace
!  type(AnharmonicData),     intent(in) :: anharmonic_data
!  real(dp),                 intent(in) :: thermal_energy
!  integer,                  intent(in) :: no_basis_states
!  type(VscfThermodynamics)             :: output
!  
!  type(StructureData) :: supercell
!  type(SubspaceBasis) :: basis
!  
!  real(dp)                :: harmonic_zpe
!  type(ThermodynamicData) :: vscha
!  
!  integer             :: no_states
!  type(HarmonicState) :: bra
!  type(HarmonicState) :: ket
!  
!  real(dp) :: kinetic
!  
!  real(dp), allocatable :: hamiltonian(:,:)
!  real(dp), allocatable :: energies(:)
!  real(dp), allocatable :: harmonic_energies(:)
!  
!  real(dp) :: exp_minus
!  
!  ! Partition functions.
!  real(dp) :: z_vscha
!  real(dp) :: z_vscf
!  real(dp) :: z_old
!  real(dp) :: z_new
!  
!  ! sum_i [E_i P_i], where P_i = exp(-E_i/T) is the bose probability before
!  !    division by the partition function.
!  real(dp) :: ep_vscha
!  real(dp) :: ep_vscf
!  real(dp) :: ep_old
!  real(dp) :: ep_new
!  
!  real(dp) :: energy
!  real(dp) :: free_energy
!  real(dp) :: entropy
!  
!  type(SymmetricEigenstuff), allocatable :: estuff(:)
!  
!  integer :: i,j,k,l,ialloc
!  
!  logical :: low_temperature
!  
!  ! z(w,T) = 1/(1-exp(-w/T))
!  ! ep(w,T) = exp(-w/T) / (1-exp(-w/T))^2
!  if (vscha_frequency > 690*thermal_energy) then
!    ! Very low-temperature regime. exp(-690)<1e-300, below dp floating point.
!    z_old = 1.0_dp
!    ep_old = 0.0_dp
!  elseif (vscha_frequency > 23*thermal_energy) then
!    ! Low-temperature regime. exp(-23)<1e-9,
!    !    so O(exp(-2w/T)) < 1e-20 is below numerical error.
!    ! z(w,T)  = 1 + exp(-w/T)            + O(exp(-2w/T))
!    ! ep(w,T) = exp(-w/T) (1+2exp(-w/T)) + O(exp(-2w/T))
!    exp_minus = exp(-vscha_frequency/thermal_energy)
!    z_old = 1.0_dp + exp_minus
!    ep_old = exp_minus*(1+2*exp_minus)
!  elseif (frequency*1e10_dp > thermal_energy) then
!    ! Usual regieme. Neither high nor low temperature limits.
!    exp_minus = exp(-vscha_frequency/thermal_energy)
!    z_old = 1/(1-exp_minus)
!    ep_old = exp_minus/(1-exp_minus)**2
!  else
!    ! High-temperature regime. O((w/T)^2) < 1e-20, below numerical error.
!    ! z(w,T) = 1/(1 - (1 - w/T + 0.5(w/T)^2 + O((w/T)^3)))
!    !        = 1/(w/T - 0.5(w/T)^2 + O((w/T)^3))
!    !        = T/w (1 + 0.5(w/T) + O((w/T)^2))
!    !        = T/w + 0.5
!    ! ep(w,T) = (1 + w/T + O((w/T)^2)) / (1- (1-w/T + 0.5(w/T)^2 + O(w/T)^3))^2
!    !         = (1 + w/T + O((w/T)^2)) / (w/T - 0.5(w/T)^2 + O((w/T)^3))^2
!    !         = (T/w + 1 + O(w/T)) / (1 - 0.5(w/T) + O((w/T)^2))^2
!    !         = (T/w + 1 + O(w/T)) / (1 - w/T + O((w/T)^2))
!    !         = (T/w + 1 + O(w/T)) * (1 + w/T + O((w/T)^2))
!    !         = T/w + 2
!    z_old = thermal_energy/vscha_frequency + 0.5_dp
!    ep_old = thermal_energy/vscha_frequency + 2.0_dp
!  endif
!  
!  z_old = z_old**size(subspace)
!  ! TODO
!  
!  
!  
!  if (output%vscha%energy-output%vscha%free_energy<thermal_energy*1e20_dp) then
!    low_temperature = .false.
!  else
!    low_temperature = .true.
!  endif
!  
!  supercell = anharmonic_data%anharmonic_supercell
!  
!  harmonic_zpe = (vscha_frequency/2.0_dp)**size(subspace)
!  output%vscha = ThermodynamicData(thermal_energy, vscha_frequency) &
!             & * size(subspace)
!  
!  basis = SubspaceBasis( subspace,                                  &
!                       & vscha_frequency,                           &
!                       & anharmonic_data%complex_modes,             &
!                       & anharmonic_data%qpoints,                   &
!                       & anharmonic_data%anharmonic_supercell,      &
!                       & no_basis_states-1,                         &
!                       & anharmonic_data%potential_expansion_order, &
!                       & anharmonic_data%structure%symmetries       )
!  
!  ! Calculate corrections, subtracting the harmonic zero point energy.
!  ep_vscha = 0.0_dp
!  ep_vscf = 0.0_dp
!  z_vscha = 0.0_dp
!  z_vscf = 0.0_dp
!  do i=1,size(basis)
!    no_states = size(basis%wavevectors(i))
!    allocate( hamiltonian(no_states,no_states), &
!            & harmonic_energies(no_states),     &
!            & stat=ialloc); call err(ialloc)
!    hamiltonian = 0
!    do j=1,no_states
!      bra = basis%wavevectors(i)%harmonic_states(j)
!      do k=1,size(basis%wavevectors(i)%harmonic_couplings(j))
!        l = basis%wavevectors(i)%harmonic_couplings(j)%id(k)
!        ket = basis%wavevectors(i)%harmonic_states(l)
!        
!        kinetic = kinetic_energy(bra, ket, subspace, supercell)
!        hamiltonian(j,l) = kinetic &
!                       & + potential%potential_energy(bra,ket,anharmonic_data)
!        if (l==j) then
!          harmonic_energies(j) = kinetic                              &
!                             & + harmonic_potential_energy( bra,      &
!                             &                              ket,      &
!                             &                              subspace, &
!                             &                              supercell )
!        endif
!      enddo
!    enddo
!    
!    estuff = diagonalise_symmetric(hamiltonian)
!    energies = estuff%eval
!    
!    ! Subtract harmonic z.p.e.
!    harmonic_energies = harmonic_energies - harmonic_zpe
!    energies = energies - harmonic_zpe
!    
!    if (low_temperature) then
!      if (i==0) then
!        energy = minval(energies)
!      else
!        energy = min(energy, minval(energies))
!      endif
!    else
!      z_vscha = z_vscha                                     &
!            & + sum(exp(-harmonic_energies/thermal_energy)) &
!            & * basis%wavevectors(i)%degeneracy
!      z_vscf = z_vscf                             &
!           & + sum(exp(-energies/thermal_energy)) &
!           & * basis%wavevectors(i)%degeneracy
!      ep_vscha = ep_vscha                                      &
!             & + sum( harmonic_energies                        &
!             &      * exp(-harmonic_energies/thermal_energy) ) &
!             & * basis%wavevectors(i)%degeneracy
!      ep_vscf = ep_vscf                              &
!            & + sum( energies                        &
!            &      * exp(-energies/thermal_energy) ) &
!            & * basis%wavevectors(i)%degeneracy
!    endif
!    deallocate( hamiltonian,       &
!              & harmonic_energies, &
!              & stat=ialloc); call err(ialloc)
!  enddo
!  
!  if (low_temperature) then
!    energy = harmonic_zpe + energy
!    free_energy = energy
!    entropy = 0
!  else
!    ! Calculate z=exp(-F/T) and ep=sum_i[E_iP_i]=U*z for the VSCHA solution.
!    ep_old = (output%vscha%energy-harmonic_zpe) * z_old
!    
!    ! Calculate the corrections to z and ep from allowing the calculated states
!    !    to relax.
!    z_new  = z_old  + z_vscf  - z_vscha
!    ep_new = ep_old + ep_vscf - ep_vscha
!    
!    ! U=sum_i[E_iP_i]/z, F=-Tln(z), S=(U-F)/T.
!    energy      = harmonic_zpe + ep_new / z_new
!    free_energy = harmonic_zpe - thermal_energy*log(z_new)
!    entropy     = thermal_energy*(energy-free_energy)
!  endif
!  
!  
!  output%vscha = ThermodynamicData(thermal_energy, vscha_frequency) &
!             & * size(subspace)
!  output%vscf = ThermodynamicData(thermal_energy, energy, free_energy, entropy)
!  
!  ! Normalise to be per primitive cell.
!  output%vscha = output%vscha / supercell%sc_size
!  output%vscf  = output%vscf  / supercell%sc_size
!end function
end module
