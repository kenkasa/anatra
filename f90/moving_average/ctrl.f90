!=======================================================================
module mod_movave_ctrl
!=======================================================================
  use mod_util
  use mod_const
  use mod_cv
  use mod_input
  use mod_output
  implicit none

  ! constants
  !
  integer, parameter :: MaxSep = 10

  ! structures
  !
  type :: s_movave_option
    real(8) :: dt                    = 0.1d0
    real(8) :: t_sta                 = 0.0d0
    real(8) :: t_sep(MaxSep - 1)     = 0.0d0
    integer :: nregion               = 1
    integer :: npoint(MaxSep)        = 5
    logical :: include_zero = .false.
  end type s_movave_option

  type :: s_movave
    integer              :: ngrid
    real(8), allocatable :: grid(:)
    real(8), allocatable :: data(:)
    real(8), allocatable :: deriv(:)
  end type s_movave

  ! subroutines
  !
  public :: read_movave_ctrl
  public :: read_movave_ctrl_option

  contains

!-----------------------------------------------------------------------
    subroutine read_movave_ctrl(input, output, option, movave)
!-----------------------------------------------------------------------
      implicit none

      integer, parameter                   :: iunit = 10

      type(s_input),               intent(out) :: input
      type(s_output),              intent(out) :: output
      type(s_movave_option),       intent(out) :: option
      type(s_movave), allocatable, intent(out) :: movave(:)

      character(len=MaxChar)       :: f_ctrl


      ! get control file name
      !
      call getarg(1, f_ctrl)

      write(iw,*)
      write(iw,'("Read_Movave_Ctrl> Reading parameters from ", a)') trim(f_ctrl)
      open(iunit, file=trim(f_ctrl), status='old')
        call read_ctrl_input          (iunit, input)
        call read_ctrl_output         (iunit, output)
        call read_movave_ctrl_option  (iunit, option)
        allocate(movave(input%ncv))
      close(iunit)

    end subroutine read_movave_ctrl
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine read_movave_ctrl_option(iunit, option, external_use)
!-----------------------------------------------------------------------
      implicit none
!
      integer,                intent(in)  :: iunit
      type(s_movave_option),  intent(out) :: option 
      logical, optional,      intent(in)  :: external_use

      ! IO
      !

      ! Local
      !
      real(8) :: dt               = 0.1d0
      real(8) :: t_sta            = 0.0d0 
      real(8) :: t_sep(MaxSep -1) = 0.0d0
      integer :: nregion          = 1
      integer :: npoint(MaxSep)   = 5
      logical :: include_zero     = .true.

      logical :: extr

      ! Dummy
      !
      integer :: i

      namelist /movave_option_param/ &
        dt,                          &
        t_sta,                       &
        t_sep,                       &
        nregion,                     &
        npoint,                      &
        include_zero

      namelist /option_param/        &
        dt,                          &
        t_sta,                       &
        t_sep,                       &
        nregion,                     &
        npoint,                      &
        include_zero

      rewind iunit

      extr = .false.
      if (present(external_use)) & 
        extr = external_use

      if (extr) then
        read(iunit, movave_option_param)
      else
        read(iunit, option_param)
      end if

      write(iw,*)
      write(iw,'(">> Option section parameters")')

      if (.not. extr) then 
        write(iw,'("dt           = ", f15.7)')     dt
        write(iw,'("t_sta        = ", f15.7)')     t_sta
      end if

      write(iw,'("t_sep        = ", 10(f15.7))') t_sep(1:2)
      write(iw,'("nregion      = ", i0)')        nregion
      write(iw,'("npoint       = ", 10(i0,2x))') (npoint(i), i = 1, nregion)
      write(iw,'("include_zero = ", a)')         get_tof(include_zero)

      option%dt           = dt
      option%t_sta        = t_sta
      option%t_sep        = t_sep
      option%nregion      = nregion
      option%npoint       = npoint
      option%include_zero = include_zero 

    end subroutine read_movave_ctrl_option
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
      !write(iw,'("file_extension = ", a)') trim(output%file_extension)


    end subroutine show_output
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end module
!=======================================================================
