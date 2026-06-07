!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_cv
  use mod_ctrl
  use mod_random
  use mod_bootstrap

  ! subroutines
  !
  public :: analyze 

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option, bootopt)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),   intent(in)    :: input
      type(s_output),  intent(in)    :: output
      type(s_option),  intent(in)    :: option
      type(s_bootopt), intent(inout) :: bootopt

      ! IO
      !
      integer                :: io
      character(len=MaxChar) :: fname

      ! Local
      !
      integer :: nfile, nsample, ntrial

      real(8) :: xsta, dx 

      ! Dummy
      !
      integer                :: ifile, jfile, itrial
      character(len=MaxChar) :: line

      ! Arrays
      !
      integer, allocatable :: rand(:, :)


      ! Setup 
      !
      nfile   = input%ncv
      ntrial  = bootopt%ntrial
      nsample = bootopt%nsample
      if (nsample == 0) &
        nsample = nfile

      allocate(rand(nsample, ntrial))

      ! Generate random seed
      !
      call get_seed         (bootopt%iseed)
      call initialize_random(bootopt%iseed)

      ! Generate random numbers 
      !
      do itrial = 1, ntrial
        call get_random_integer(nsample,            &
                                1,                  &
                                nfile,              &
                                bootopt%duplicate,  &
                                rand(1, itrial))
      end do

      ! Make file lists
      !
      do itrial = 1, ntrial
        write(fname,'(a, i4.4, ".", a)') &
          trim(output%fhead), itrial, trim(output%file_extension)

        call open_file(fname, io)
        do ifile = 1, nsample
          jfile = rand(ifile, itrial)
          write(io,'(a)') trim(input%fcv(jfile))
        end do
        close(io)
      end do
     

    end subroutine analyze
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
