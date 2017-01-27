module calculate_bs_module
  implicit none
contains

subroutine calculate_bs(args)
  use constants, only : dp, kB
  use file_module
  use string_module
  implicit none
  
  ! input variables
  type(String), intent(in) :: args(:)
  
  ! Parameters
  real(dp), parameter :: dtemperature = 50.0d0
  
  ! Input variables
  integer :: no_kpoints,no_modes,degeneracy
  integer, allocatable :: multiplicity(:)
  real(dp),allocatable :: bands(:,:),frequency(:,:),deformation(:,:)
  real(dp) :: mapping_amplitude
  
  ! Working variables
  real(dp) :: temperature
  integer  :: i,j,k,total_kpoints
  real(dp) :: renormalised_band,renormalised_band_kpoint
  real(dp) :: db
  real(dp) :: dump_r
  
  ! filenames
  type(String) :: ibz_filename
  type(String) :: in_dir
  type(String) :: file_end
  type(String) :: band_gap_correction_filename
  type(String) :: bg_correction_kp_filename
  
  ! file units
  integer :: ibz_file
  integer :: freq_file
  integer :: bs_file
  integer :: bgc_file
  integer :: bck_file
  
  ! Process inputs
  no_modes = int(args(1))
  degeneracy = int(args(2))
  mapping_amplitude = dble(args(3))
  ibz_filename = args(4)
  in_dir = args(5)
  band_gap_correction_filename = args(6)
  bg_correction_kp_filename = args(7)
  
  no_kpoints = count_lines(ibz_filename)
  
  ! allocate and initialise arrays
  allocate(bands(no_kpoints,no_modes))
  allocate(multiplicity(no_kpoints))
  allocate(frequency(no_kpoints,no_modes))
  allocate(deformation(no_kpoints,no_modes))
  bands=0.0
  frequency=1.0
  deformation=0.0
  
  ! read ibz file
  ibz_file = open_read_file(ibz_filename)
  do i=1,no_kpoints 
    read(ibz_file,*)dump_r,dump_r,dump_r,multiplicity(i)
  enddo
  close(ibz_file)
  
  total_kpoints = sum(multiplicity)

  ! Read in bands
  do i=1,no_kpoints
    do j=1,no_modes
      if (i==1 .and. j<4) cycle ! skip acoustic modes
      
      file_end = str('.')//i//'.'//j//'.dat'
      
      freq_file = open_read_file(in_dir//'/frequency'//file_end)
      read(freq_file,*) frequency(i,j)
      close(freq_file)
      
      bands(i,j) = 0
      bs_file = open_read_file(in_dir//'/bs'//file_end)
      
      ! b1
      do k=1,degeneracy
        read(bs_file,*) db
        bands(i,j) = bands(i,j)+db/(2.d0*degeneracy)
      enddo
      
      ! b2
      do k=1,degeneracy
        ! TODO: for the i/=1 case, this used the db value above. Bug?
        read(bs_file,*) db 
        bands(i,j) = bands(i,j)+db/(2.d0*degeneracy)
      enddo
      
      ! b0
      read(bs_file,*) db
      bands(i,j) = bands(i,j)-db
      
      close(bs_file)
    enddo
  enddo

  ! Calculate deformation potential
  deformation = bands/mapping_amplitude**2
  
  ! Calculate quadratic vibrational correction
  bgc_file = open_write_file(band_gap_correction_filename)
  bck_file = open_write_file(bg_correction_kp_filename)
  do k=1,21  ! loop over temperature
    renormalised_band=0.0
    temperature = (k-1)*dtemperature
    if(temperature<1.d-5)then
      do i=1,no_kpoints
        write(bck_file,*)'k-point',i
        renormalised_band_kpoint=0.0
        do j=1,no_modes
          renormalised_band = renormalised_band &
                          & + deformation(i,j)*multiplicity(i)/total_kpoints
          renormalised_band_kpoint = renormalised_band_kpoint+deformation(i,j)
          write(bck_file,*)i,j,deformation(i,j)
        enddo
        write(bck_file,*)i,renormalised_band_kpoint
      enddo
    else
      do i=1,no_kpoints
        do j=1,no_modes
          renormalised_band = renormalised_band                               &
                          & + deformation(i,j)                                &
                          & * (1.0                                            &
                          &   +2.0/(dexp(frequency(i,j)/(temperature*kB))-1)) &
                          & * multiplicity(i)                                 &
                          & / total_kpoints
        enddo
      enddo
    endif ! temperature
    write(bgc_file,*)temperature,renormalised_band 
  enddo
  close(bgc_file)
  close(bck_file)
end subroutine
end module
