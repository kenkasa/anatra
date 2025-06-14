!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_ctrl
  use mod_bootstrap
  use mod_random

  ! constants
  !

  ! structures
  !
  type :: s_booteach
    integer, allocatable :: rand(:, :)
    real(8), allocatable :: func(:, :)

  end type s_booteach

  type :: s_bootave
    real(8), allocatable :: ave(:)
    real(8), allocatable :: err(:) 
  end type s_bootave

  ! subroutines
  !
  public  :: analyze

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),   intent(in)    :: input
      type(s_output),  intent(in)    :: output
      type(s_option),  intent(in)    :: option

      ! I/O
      !
      integer                :: io
      character(len=MaxChar) :: fname 

      ! Local
      !
      integer :: ndim, nstep, nfile
      real(8) :: x, val, ave, dev, sterr

      ! Dummy
      !
      integer :: ifile, istep

      ! Arrays
      !
      type(s_cv), allocatable :: cv(:)
      real(8),    allocatable :: func_ave(:), func_stdev(:), func_sterr(:) 


      ! Setup
      !
      ndim  = 1
      nfile = input%ncv
      allocate(cv(nfile))

      ! Read CV files
      !
      write(iw,*)
      write(iw,'("Analyze> Read CV file")')
      do ifile = 1, nfile 
        call read_cv(input%fcv(ifile), ndim, cv(ifile))
      end do

      nstep = cv(1)%nstep
      allocate(func_ave(nstep), func_stdev(nstep), func_sterr(nstep))

      func_ave   = 0.0d0
      func_stdev = 0.0d0
      func_sterr = 0.0d0
      
      do istep = 1, nstep
        ave   = 0.0d0
        dev   = 0.0d0
        sterr = 0.0d0

        do ifile = 1, nfile
          ave = ave + cv(ifile)%data(1, istep)  
        end do
        ave = ave / dble(nfile)

        do ifile = 1, nfile
          dev = dev + (cv(ifile)%data(1, istep) - ave)**2
        end do
        dev   = sqrt(dev / (nfile - 1))
        sterr = dev / sqrt(dble(nfile))

        func_ave(istep)   = ave
        func_stdev(istep) = dev
        func_sterr(istep) = sterr
      end do

      ! Output
      !
      write(fname, '(a,".ave")') trim(output%fhead)
      call open_file(fname, io)

      do istep = 1, nstep 
        x = option%xsta + option%dx * (istep - 1)
        write(io,'(e15.7,2x)',advance='no') x
        write(io,'(e15.7,2x)',advance='no') func_ave(istep)
        if (nfile > 2) then
          !write(io,'(e15.7,2x)',advance='no') func_stdev(istep)
          write(io,'(e15.7,2x)',advance='no') func_sterr(istep)
        end if
        write(io,*)
      end do

      close(io)

      ! Deallocate
      !
      deallocate(func_ave, func_stdev, func_sterr)


    end subroutine analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine analyze_bootstrap(input, output, option, bootopt)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),   intent(in)    :: input
      type(s_output),  intent(in)    :: output
      type(s_option),  intent(in)    :: option
      type(s_bootopt), intent(inout) :: bootopt 

      type(s_booteach) :: beach
      type(s_bootave)  :: bave

      ! I/O
      !
      integer                :: io
      character(len=MaxChar) :: fname

      ! Local
      !
      integer :: ndim, ntrial, nsample
      integer :: nstep, nfile
      real(8) :: x, val, ave, dev, sterr

      ! Dummy
      !
      integer :: ifile, isample, istep, ig, itrial

      ! Arrays
      !
      type(s_cv), allocatable :: cv(:)
      real(8),    allocatable :: func_ave(:)
      !, func_stdev(:), func_sterr(:) 


      ! Setup
      !
      ndim  = 1
      nfile = input%ncv
      allocate(cv(nfile))

      ! Read CV files
      !
      write(iw,*)
      write(iw,'("Analyze> Read CV file")')
      do ifile = 1, nfile
        call read_cv(input%fcv(ifile), ndim, cv(ifile))
      end do

      nstep = cv(1)%nstep

      ! Setup Bootstrap
      !
      ntrial  = bootopt%ntrial
      nsample = bootopt%nsample

      allocate(beach%func(nstep, ntrial))
      allocate(bave%ave(nstep), bave%err(nstep))

      ! Generate random numbers 
      !
      call get_seed(bootopt%iseed)
      call initialize_random(bootopt%iseed)

      allocate(beach%rand(nsample, ntrial))

      do itrial = 1, ntrial
        call get_random_integer(nsample, 1, nfile, bootopt%duplicate, beach%rand(1, itrial)) 
      end do

      ! Run Bootstrap
      !
      allocate(func_ave(nstep))

      beach%func = 0.0d0

      !$omp parallel private(itrial, isample, istep, val, func_ave) &
      !$omp          default(shared)
      !$omp do
      do itrial = 1, ntrial

        if (mod(itrial, 10) == 0) then
          write(iw,'("trial : ", i0)') itrial
        end if

        func_ave = 0.0d0

        do isample = 1, nsample
          ifile = beach%rand(isample, itrial)

          do istep = 1, nstep
            val             = cv(ifile)%data(1, istep) 
            func_ave(istep) = func_ave(istep) + val 
          end do

        end do

        func_ave                    = func_ave / dble(nsample)
        beach%func(1:nstep, itrial) = func_ave(1:nstep) 

      end do

      !$omp end do
      !$omp end parallel

      ! Calculate average & error
      !
      bave%ave = 0.0d0
      bave%err = 0.0d0

      do istep = 1, nstep
        bave%ave(istep) = sum(beach%func(istep, 1:ntrial)) / dble(ntrial)
      end do

      do istep = 1, nstep
        dev = 0.0d0
        ave = bave%ave(istep)
        do itrial = 1, ntrial
          val = beach%func(istep, itrial)
          dev = dev + (val - ave)**2 
        end do
        bave%err(istep) = sqrt(dev / dble(ntrial - 1))
      end do 

      ! Output
      !
      write(fname, '(a,".ave")') trim(output%fhead)
      call open_file(fname, io)

      do istep = 1, nstep 
        x = option%xsta + option%dx * (istep - 1)
        write(io,'(e15.7,2x)',advance='no') x
        write(io,'(e15.7,2x)',advance='no') bave%ave(istep)
        write(io,'(e15.7,2x)')              bave%err(istep)
      end do

      close(io)

      ! Deallocate
      !
      deallocate(cv)
      deallocate(func_ave)
      deallocate(beach%rand, beach%func)
      deallocate(bave%ave, bave%err)


    end subroutine analyze_bootstrap
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
