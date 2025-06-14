!=======================================================================
module mod_ctrl
!=======================================================================
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_cv
  use mod_bootstrap
  implicit none

  ! constants
  !

  ! structures
  !
  type :: s_option
    logical :: use_bootstrap = .false.
    real(8) :: xsta          = 0.0d0
    real(8) :: dx            = 0.1d0
  end type s_option

  ! subroutines
  !
  public  :: read_ctrl
  private :: read_ctrl_option
  private :: show_input
  private :: show_output

  contains

!-----------------------------------------------------------------------
    subroutine read_ctrl(input, output, option, bootopt)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),   intent(out) :: input
      type(s_output),  intent(out) :: output
      type(s_option),  intent(out) :: option 
      type(s_bootopt), intent(out) :: bootopt 

      ! I/O
      !
      integer                :: io
      character(len=MaxChar) :: f_ctrl


      ! get control file name
      !
      call getarg(1, f_ctrl)

      write(iw,*)
      write(iw,'("Read_Ctrl> Reading parameters from ", a)') trim(f_ctrl)
      call open_file(f_ctrl, io)
      call read_ctrl_input  (io, input)
      call show_input(input)
      call read_ctrl_output (io, output)
      call show_output(output)
      call read_ctrl_option (io, option)

      if (option%use_bootstrap) then
        call read_ctrl_bootstrap(io, bootopt)
      end if

      close(io)

      !open(iunit, file=trim(input%flist))
      !  call read_ctrl_cvinfo  (iunit, cvinfo)
      !close(io) 

    end subroutine read_ctrl
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine read_ctrl_option(io, option)
!-----------------------------------------------------------------------
      implicit none
!
      integer, parameter :: ndim_max = 3 
!
      integer,        intent(in)  :: io
      type(s_option), intent(out) :: option 

      logical :: use_bootstrap = .false.
      real(8) :: xsta          = 0.0d0
      real(8) :: dx            = 0.1d0

      integer :: i, j
      integer :: iopt, ierr

      namelist /option_param/ &
              use_bootstrap,  &
              xsta,           &
              dx


      rewind io 
      read(io, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("use_bootstrap = ", a)') get_tof(use_bootstrap)
      write(iw,'("xsta          = ", f20.10)') xsta 
      write(iw,'("dx            = ", f20.10)') dx

      option%use_bootstrap = use_bootstrap
      option%xsta          = xsta
      option%dx            = dx

    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine show_input(input)
!-----------------------------------------------------------------------
      implicit none

      type(s_input), intent(in) :: input

      ! Dummy
      !
      integer :: i

      ! Check
      !
      if (input%ncv == 0) then
        write(iw,'("Error. fcv should be specified.")')
        stop
      end if

      ! Print
      !
      write(iw,*)
      write(iw,'(">> Input section parameters")')
      do i = 1, input%ncv
        write(iw,'("fcv", 3x, i0, 3x, " = ", a)') i, trim(input%fcv(i)) 
      end do


    end subroutine show_input
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine show_output(output)
!-----------------------------------------------------------------------
      implicit none

      type(s_output), intent(in) :: output 


      ! Check
      !
      if (trim(output%fhead) == '') then
        write(iw,'("Error. fhead should be specified.")')
        stop
      end if

      ! Print
      !
      write(iw,*)
      write(iw,'(">> Output section parameters")')
      write(iw,'("fhead          = ", a)') trim(output%fhead)


    end subroutine show_output
!-----------------------------------------------------------------------

end module
!=======================================================================
