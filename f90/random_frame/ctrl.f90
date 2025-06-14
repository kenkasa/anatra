!=======================================================================
module mod_ctrl
!=======================================================================
  use mod_util
  use mod_const
  use mod_traj
  use mod_input
  use mod_output
  implicit none

  ! constants
  !

  ! structures
  !
  type :: s_option
    integer              :: nsample         =  1
    integer              :: iseed           = -1
    integer              :: output_trajtype = TrajTypeDCD 
    logical              :: duplicate       = .false.
    logical              :: out_rst7        = .false.
    logical              :: use_allsnap     = .false.
    logical              :: shuffle         = .false.
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
    subroutine read_ctrl_option(io, option)
!-----------------------------------------------------------------------
      implicit none
!
      integer, parameter :: ndim_max = 3 
!
      integer,        intent(in)  :: io 
      type(s_option), intent(out) :: option 

      integer                :: nsample         = 100
      integer                :: iseed           = -1
      character(len=MaxChar) :: output_trajtype = 'DCD' 
      logical                :: duplicate       = .false.
      logical                :: out_rst7        = .false.
      logical                :: use_allsnap     = .false.
      logical                :: shuffle         = .false.

      ! Parser
      !
      integer :: iopt, ierr

      ! Dummy
      !
      integer :: i, j

      namelist /option_param/ nsample,         &
                              iseed,           &
                              output_trajtype, &
                              duplicate,       &
                              out_rst7,        &
                              use_allsnap,     &
                              shuffle

      rewind io 
      read(io, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("nsample         = ", i0)') nsample 
      write(iw,'("iseed           = ", i0)') iseed 
      write(iw,'("output_trajtype = ", a)')  trim(output_trajtype)
      write(iw,'("duplicate       = ", a)')  get_tof(duplicate) 
      write(iw,'("out_rst7        = ", a)')  get_tof(out_rst7) 
      write(iw,'("use_allsnap     = ", a)')  get_tof(use_allsnap) 
      write(iw,'("shuffle         = ", a)')  get_tof(shuffle) 
      write(iw,*)

      iopt = get_opt(output_trajtype, TrajTypes, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("output_trajtype = ",a," is not available.")') trim(output_trajtype)
        stop
      end if
      option%output_trajtype = iopt
                           
      option%nsample       = nsample 
      option%iseed         = iseed
      option%duplicate     = duplicate
      option%use_allsnap   = use_allsnap
      option%out_rst7      = out_rst7
      option%shuffle       = shuffle

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
