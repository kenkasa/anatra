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
  integer,      parameter, public :: OrientModeRESIDUE = 1
  integer,      parameter, public :: OrientModeWHOLE   = 2 
  character(*), parameter, public :: OrientMode(2) = (/'RESIDUE   ',&
                                                       'WHOLE     '/)

  ! structures
  !
  type :: s_option
    integer :: mode     = CoMModeRESIDUE 
    integer :: tcfrange = 1000
  end type s_option

  ! subroutines
  !
  public  :: read_ctrl
  private :: read_ctrl_option

  contains

!-----------------------------------------------------------------------
    subroutine read_ctrl(input, output, option, trajopt)
!-----------------------------------------------------------------------
      implicit none

      integer, parameter           :: iunit = 10

      type(s_input),   intent(out) :: input
      type(s_output),  intent(out) :: output
      type(s_option),  intent(out) :: option
      type(s_trajopt), intent(out) :: trajopt

      character(len=MaxChar)       :: f_ctrl

      ! get control file name
      !
      call getarg(1, f_ctrl)

      write(iw,*)
      write(iw,'("Read_Ctrl> Reading parameters from ", a)') trim(f_ctrl)
      open(iunit, file=trim(f_ctrl), status='old')
        call read_ctrl_input  (iunit, input)
        call read_ctrl_output (iunit, output)
        call read_ctrl_option (iunit, option)
        call read_ctrl_trajopt(iunit, trajopt)
      close(iunit)

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

      character(len=MaxChar) :: mode     = "RESIDUE"
      integer                :: tcfrange = 1000 

      integer                :: iopt, ierr

      namelist /option_param/ mode, tcfrange 

      rewind iunit
      read(iunit, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("mode     = ", a)')      trim(mode)
      write(iw,'("tcfrange = ", i0)')     tcfrange 


      iopt = get_opt(mode, CoMMode, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("mode = ",a," is not available.")') trim(mode)
        stop
      end if
      option%mode     = iopt

      option%tcfrange = tcfrange 

    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end module
!=======================================================================
