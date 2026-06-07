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


  ! structures
  !
  type :: s_option
    integer :: mode      = CoMModeRESIDUE
    integer :: nparallel = 1
    real(8) :: dt        = 1.0d0
    real(8) :: rcut      = 3.5d0
  end type s_option

  ! subroutines
  !
  public  :: read_ctrl
  private :: read_ctrl_option
  private :: show_input
  private :: show_output

  contains

!-----------------------------------------------------------------------
    subroutine read_ctrl(input, output, option, trajopt)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),    intent(out) :: input
      type(s_output),   intent(out) :: output
      type(s_option),   intent(out) :: option
      type(s_trajopt),  intent(out) :: trajopt

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

      call read_ctrl_option (io, option)
      call read_ctrl_trajopt(io, trajopt)

      close(io)

    end subroutine read_ctrl
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine read_ctrl_option(io, option)
!-----------------------------------------------------------------------
      implicit none
!
      integer,          intent(in)  :: io
      type(s_option),   intent(out) :: option 

      ! Local
      !
      character(len=MaxChar) :: mode      = "RESIDUE"
      integer                :: nparallel = 1
      real(8)                :: dt        = 1.0d0
      real(8)                :: rcut      = 3.5d0

      ! Parser 
      !
      integer :: iopt, ierr

      ! Dummy
      !
      integer :: it
      real(8) :: t 

      namelist /option_param/ mode,        &
                              nparallel,   &
                              dt,          &
                              rcut


      rewind io
      read(io, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("mode      = ", a)')     trim(mode)
      write(iw,'("dt        = ", f15.7)') dt
      write(iw,'("rcut      = ", f15.7)') rcut 

      ! Get mode
      !
      iopt = get_opt(mode, CoMMode, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("mode = ",a," is not available.")') trim(mode)
        stop
      end if
      option%mode      = iopt

      ! Namelist => Option
      !
      option%nparallel = nparallel
      option%dt        = dt
      option%rcut      = rcut


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
