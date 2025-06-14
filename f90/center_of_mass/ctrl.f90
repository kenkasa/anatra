!=======================================================================
module mod_ctrl
!=======================================================================
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_traj
  use mod_com
  implicit none

  ! constants
  !
  integer,      parameter, public :: CoMFormatTypeXYZ        = 1
  integer,      parameter, public :: CoMFormatTypeTIMESERIES = 2 
  character(*), parameter, public :: CoMFormatType(2)        = (/'XYZ       ',&
                                                                 'TIMESERIES'/)


  ! structures
  !
  type :: s_option
    integer :: mode      = CoMModeRESIDUE
    integer :: comformat = CoMFormatTypeXYZ 
    integer :: msddim    = 3
    real(8) :: dt        = 1.0d0
    real(8) :: t_sparse  = -1.0d0
    real(8) :: t_range   = -1.0d0
    logical :: out_com   = .false. 
    logical :: out_msd   = .false.
    logical :: onlyz     = .false.
    logical :: unwrap    = .false.

    ! prepared after reading namelists
    !
    real(8) :: dt_out    = 0.0d0
    integer :: nstep     = 0
    integer :: nt_shift  = 0
    integer :: nt_range  = 0
  end type s_option

  type :: s_timegrid
    integer :: ng
    real(8), allocatable :: val(:)
    integer, allocatable :: ind(:)
  end type s_timegrid

  ! subroutines
  !
  public  :: read_ctrl
  private :: read_ctrl_option
  private :: show_input
  private :: show_output

  contains

!-----------------------------------------------------------------------
    subroutine read_ctrl(input, output, option, trajopt, timegrid)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),    intent(out) :: input
      type(s_output),   intent(out) :: output
      type(s_option),   intent(out) :: option
      type(s_trajopt),  intent(out) :: trajopt
      type(s_timegrid), intent(out) :: timegrid 

      ! I/O
      !
      integer                      :: io
      character(len=MaxChar)       :: f_ctrl

      ! Get control file name
      !
      call getarg(1, f_ctrl)

      write(iw,*)
      write(iw,'("Read_Ctrl> Reading parameters from ", a)') trim(f_ctrl)

      call open_file(f_ctrl, io, stat = 'old')

      call read_ctrl_input  (io, input)
      call show_input       (input)

      call read_ctrl_output (io, output)
      call show_output      (output)

      call read_ctrl_option (io, option, timegrid)
      call read_ctrl_trajopt(io, trajopt)

      close(io)

    end subroutine read_ctrl
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine read_ctrl_option(io, option, timegrid)
!-----------------------------------------------------------------------
      implicit none
!
      integer,          intent(in)  :: io
      type(s_option),   intent(out) :: option 
      type(s_timegrid), intent(out) :: timegrid 

      ! Local
      !
      character(len=MaxChar) :: mode      = "RESIDUE"
      character(len=MaxChar) :: comformat = "XYZ"
      logical                :: unwrap    = .false. 
      logical                :: out_com   = .false.
      logical                :: out_msd   = .false.
      logical                :: onlyz     = .false.
      integer                :: msddim    = 3 
      real(8)                :: dt        = 1.0d0
      real(8)                :: t_sparse  = -1.0d0
      real(8)                :: t_range   = -1.0d0

      ! Parser 
      !
      integer :: iopt, ierr

      ! Dummy
      !
      integer :: it
      real(8) :: t 

      namelist /option_param/ mode,        &
                              comformat,   &
                              unwrap,      &
                              out_msd,     &
                              out_com,     &
                              onlyz,       &
                              msddim,      &
                              dt,          &
                              t_sparse,    &
                              t_range


      rewind io
      read(io, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("mode      = ", a)')     trim(mode)
      write(iw,'("comformat = ", a)')     trim(comformat)
      write(iw,'("unwrap    = ", a)')     get_tof(unwrap)
      write(iw,'("out_com   = ", a)')     get_tof(out_com)
      write(iw,'("out_msd   = ", a)')     get_tof(out_msd)
      write(iw,'("onlyz     = ", a)')     get_tof(onlyz)
      write(iw,'("msddim    = ", i0)')    msddim 
      write(iw,'("dt        = ", f15.7)') dt
      write(iw,'("t_sparse  = ", f15.7)') t_sparse 
      write(iw,'("t_range   = ", f15.7)') t_range

      ! Get mode
      !
      iopt = get_opt(mode, CoMMode, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("mode = ",a," is not available.")') trim(mode)
        stop
      end if
      option%mode      = iopt

      ! Get comformat
      !
      iopt = get_opt(comformat, CoMFormatType, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("comformat = ",a," is not available.")') trim(comformat)
        stop
      end if
      option%comformat = iopt

      ! Namelist => Option
      !
      option%msddim    = msddim
      option%out_com   = out_com 
      option%out_msd   = out_msd
      option%onlyz     = onlyz
      option%unwrap    = unwrap
      option%dt        = dt
      option%t_sparse  = t_sparse
      option%t_range   = t_range

      ! Check combinations
      !

      ! Convert from real to integer
      !
      if (option%t_sparse < 0.0d0) then
        option%t_sparse = option%dt
      end if

      option%nt_shift  = nint(option%t_sparse / option%dt)
      option%dt_out    = option%dt * option%nt_shift

      if (option%t_range > 0.0d0) then
        option%nt_range = nint(option%t_range / option%dt)
      else
        write(iw,'("Read_Ctrl_Option> Error. t_range should be specified.")')
        stop
      end if

      if (option%nt_shift > 1) then
        option%nt_range = int(dble(option%nt_range) / dble(option%nt_shift))
      end if

      write(iw,'("dt_out    = ", f15.7)') option%dt_out
      write(iw,'("nt_shift  = ", i0)')    option%nt_shift

      ! Setup Time grids
      !
      allocate(timegrid%val(0:option%nt_range))
      allocate(timegrid%ind(0:option%nt_range))

      timegrid%ng = option%nt_range

      do it = 0, timegrid%ng
        timegrid%ind(it) = option%nt_shift * it
        timegrid%val(it) = option%dt_out   * it
      end do


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

!-----------------------------------------------------------------------
end module
!=======================================================================
