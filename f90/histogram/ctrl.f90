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

  ! structures
  !
  type :: s_option
    real(8) :: dx       = 0.1d0
    integer :: xsta     = 0.0d0 
    integer :: nx       = 100
    integer :: ncol     = 1
    integer :: data_sta = 1 
    integer :: data_end = 0
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
      write(iw,'("Read_Ctrl> Read parameters from ", a)') trim(f_ctrl)

      call open_file(f_ctrl, io, stat = 'old')

      call read_ctrl_input  (io, input)
      call show_input       (input)

      call read_ctrl_output (io, output)
      call show_output      (output)

      call read_ctrl_option (io, option)

      close(io)

    end subroutine read_ctrl
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine read_ctrl_option(iunit, option)
!-----------------------------------------------------------------------
      implicit none
!
      integer,        intent(in)  :: iunit
      type(s_option), intent(out) :: option 

      real(8)                :: dx       = 0.1d0
      real(8)                :: xsta     = 0.0d0 
      integer                :: nx       = 100
      integer                :: ncol     = 1
      integer                :: data_sta = 1 
      integer                :: data_end = 0


      namelist /option_param/ &
        dx,                   &
        xsta,                 &
        nx,                   &
        ncol,                 &
        data_sta,             &
        data_end

      rewind iunit
      read(iunit, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("dx       = ", f15.7)')   dx
      write(iw,'("xsta     = ", f15.7)')   xsta 
      write(iw,'("nx       = ", i0)')      nx
      write(iw,'("ncol     = ", i0)')      ncol
      write(iw,'("data_sta = ", i0)')      data_sta 
      write(iw,'("data_end = ", i0)')      data_end

      option%dx       = dx
      option%xsta     = xsta
      option%nx       = nx
      option%ncol     = ncol
      option%data_sta = data_sta
      option%data_end = data_end


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
      write(iw,'("file_extension = ", a)') trim(output%file_extension)


    end subroutine show_output
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end module
!=======================================================================
