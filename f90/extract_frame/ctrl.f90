!=======================================================================
module mod_ctrl
!=======================================================================
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  implicit none

  ! constants
  !
  integer, parameter, public :: Nstate   = 1 
  integer, parameter, public :: REACTIVE = 1 
  integer, parameter, public :: OTHERS   = 2 
  integer, parameter, public :: StateInfo(Nstate) = (/1/)

  ! structures
  !
  type :: s_option
    integer              :: ndim         = 1
    real(8)              :: dt           = 1.0d0
    real(8), allocatable :: state_def(:, :, :)
    real(8), allocatable :: react_range(:, :)
  end type s_option

  ! subroutines
  !
  public  :: read_ctrl
  private :: read_ctrl_option
  private :: show_input
  private :: show_output

  contains

!-----------------------------------------------------------------------
    subroutine read_ctrl(input, output, option)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),   intent(out) :: input
      type(s_output),  intent(out) :: output
      type(s_option),  intent(out) :: option

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
      close(io)

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

      integer :: ndim                       = 1
      real(8) :: dt                         = 1.0d0
      real(8) :: react_range(2, ndim_max)   = 0.0d0

      ! Dummy
      !
      integer :: i, j

      namelist /option_param/ &
        ndim,                 &
        dt,                   &
        react_range


      rewind io
      read(io, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("ndim          = ", i0)') ndim
      write(iw,*)
      write(iw,'("react_range   =")')
      do i = 1, ndim
        write(iw,'(es15.7, " <= component ",i0," < ",es15.7)') &
          react_range(1, i), i, react_range(2, i) 
      end do

      ! Send
      !
      option%ndim          = ndim
      option%dt            = dt

      allocate(option%react_range(1:2, ndim))
      allocate(option%state_def(2, ndim, Nstate))

      do i = 1, ndim
        option%react_range(1:2, i)   = react_range(1:2, i)
      end do

      do i = 1, ndim
        option%state_def(1:2, i, 1) = react_range(1:2, i)
      end do

      ! Combination check
      !


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
      if (input%ntraj == 0) then
        write(iw,'("Error. ftraj should be specified.")')
        stop
      end if

      ! Print
      !
      write(iw,*)
      write(iw,'(">> Input section parameters")')
      do i = 1, input%ntraj
        write(iw,'("ftraj", 3x, i0, 3x, " = ", a)') i, trim(input%ftraj(i))
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
