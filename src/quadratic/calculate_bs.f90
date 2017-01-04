module calculate_bs_module
  implicit none
contains

subroutine calculate_bs(args)
  use constants, only : dp, kB
  use utils,     only : i2s
  use file_io,   only : open_read_file, open_write_file
  implicit none
  
  ! input variables
  character(100), intent(in) :: args(:)
  
  ! Parameters
  real(dp), parameter :: dtemperature = 50.0d0
  
  ! Input variables
  integer :: no_kpoints,no_modes,degeneracy
  integer, allocatable :: multiplicity(:)
  real(dp),allocatable :: bands(:,:),frequency(:,:),deformation(:,:)
  real(dp) :: mapping_amplitude
  
  ! Working variables
  real(dp) :: temperature
  integer :: i,j,k,total_kpoints
  real(dp) :: renormalised_band,renormalised_band_kpoint
  real(dp), allocatable :: db1(:),db2(:)
  real(dp) :: b1,b2,b0,dump_r
  
  ! filenames
  character(100) :: ibz_filename
  character(100) :: in_dir
  character(100) :: file_end
  character(100) :: band_gap_correction_filename
  character(100) :: bg_correction_kp_filename
  
  ! file units
  integer :: ibz_file
  integer :: freq_file
  integer :: bs_file
  integer :: bgc_file
  integer :: bck_file
  
  ! Process inputs
  read(args(1),*) no_kpoints
  read(args(2),*) no_modes
  read(args(3),*) degeneracy
  read(args(4),*) mapping_amplitude
  ibz_filename = args(5)
  in_dir = args(6)
  band_gap_correction_filename = args(7)
  bg_correction_kp_filename = args(8)
  
  ! allocate and initialise arrays
  allocate(bands(no_kpoints,no_modes))
  allocate(multiplicity(no_kpoints))
  allocate(frequency(no_kpoints,no_modes))
  allocate(deformation(no_kpoints,no_modes))
  bands=0.0
  frequency=1.0
  deformation=0.0
  allocate(db1(degeneracy))
  allocate(db2(degeneracy))
  
  ! read ibz file
  ibz_file = open_read_file(ibz_filename)
  do i=1,no_kpoints 
    read(ibz_file,*)dump_r,dump_r,dump_r,multiplicity(i)
  enddo
  close(ibz_file)
  
  total_kpoints = sum(multiplicity)

  ! Read in bands
  do i=1,no_kpoints
    if(i==1)then
      do j=4,no_modes
        file_end = '.'//trim(i2s(i))//'.'//trim(i2s(j))//'.dat'
        
        freq_file = open_read_file(trim(in_dir)//'/frequency'//trim(file_end))
        read(freq_file,*) frequency(i,j)
        close(freq_file)
        
        bs_file = open_read_file(trim(in_dir)//'/bs'//trim(file_end))
        b1=0.0
        do k=1,degeneracy
          read(bs_file,*)db1(k) 
          b1=b1+(db1(k)/degeneracy)
        enddo
        b2=0.0
        do k=1,degeneracy
          read(bs_file,*)db2(k) 
          b2=b2+(db2(k)/degeneracy)
        enddo
        read(bs_file,*)b0 
        close(bs_file)
        
        bands(i,j)=(b1+b2)/2.0-b0      
      enddo
    else
      do j=1,no_modes
        file_end = '.'//trim(i2s(i))//'.'//trim(i2s(j))//'.dat'
        
        freq_file = open_read_file(trim(in_dir)//'/frequency'//trim(file_end))
        read(freq_file,*) frequency(i,j)
        close(freq_file)
        
        bs_file = open_read_file(trim(in_dir)//'/bs'//trim(file_end))
        b1=0.0
        do k=1,degeneracy
          read(bs_file,*)db1(k) 
          b1=b1+(db1(k)/degeneracy)
        enddo
        b2=0.0
        do k=1,degeneracy
          read(bs_file,*)db2(k) 
          b2=b2+(db1(k)/degeneracy)
        enddo
        read(bs_file,*)b0 
        close(bs_file)
        
        bands(i,j)=(b1+b2)/2.0-b0      
      enddo
    endif ! skip acoustic modes
  enddo

  ! Calculate deformation potential
  do i=1,no_kpoints
    do j=1,no_modes
      deformation(i,j) = bands(i,j)                        &
                     & / ( mapping_amplitude               &
                     &   / dsqrt(2.0*dabs(frequency(i,j))) &
                     &   )**2                              &
                     & / (2.0*dabs(frequency(i,j)))
    enddo
  enddo
  
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
