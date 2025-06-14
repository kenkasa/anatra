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
  integer,      parameter, public :: OrientTypeCOSTHETA = 1
  integer,      parameter, public :: OrientTypeTHETA    = 2 
  character(*), parameter, public :: OrientTypes(2) = (/'COSTHETA  ',&
                                                        'THETA     '/)

  integer,      parameter, public :: ZdefTypeCOORD      = 1
  integer,      parameter, public :: ZdefTypeOUTERPROD  = 2 
  character(*), parameter, public :: ZdefTypes(2) = (/'COORD     ',&
                                                      'OUTERPROD '/)


  integer,      parameter, public :: JudgeUpModeNMOLUP = 1
  integer,      parameter, public :: JudgeUpModeCOORD  = 2 
  integer,      parameter, public :: JudgeUpModeNONE   = 3 
  character(*), parameter, public :: JudgeUpMode(3) = (/'NMOLUP   ',&
                                                        'COORD    ',&
                                                        'NONE     '/)

  ! structures
  !
  type :: s_option
    integer :: judgeup      = JudgeUpModeCOORD
    integer :: zdef_type    = ZdefTypeCOORD 
    integer :: mode(7)      = (/CoMModeRESIDUE, CoMModeRESIDUE,  CoMModeWHOLE, &
                                CoMModeWHOLE,   ComModeWHOLE,    CoMModeWHOLE, &
                                CoMModeWHOLE/)
    integer :: orient_type  = OrientTypeCOSTHETA 
    real(8) :: dx           = 0.1d0
    integer :: nmolup       = 0

    integer :: nsel         = 7 
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

      type(s_input),   intent(out) :: input
      type(s_output),  intent(out) :: output
      type(s_option),  intent(out) :: option
      type(s_trajopt), intent(out) :: trajopt

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
      call read_ctrl_trajopt(io, trajopt)

      close(io)

    end subroutine read_ctrl
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine read_ctrl_option(io, option)
!-----------------------------------------------------------------------
      implicit none
!
      integer,        intent(in)  :: io 
      type(s_option), intent(out) :: option 

      character(len=MaxChar) :: judgeup      = "COORD"
      character(len=MaxChar) :: zdef_type    = "COORD" 
      character(len=MaxChar) :: mode(3)      = (/"RESIDUE", "RESIDUE", "RESIDUE"/) 
      character(len=MaxChar) :: orient_type  = "COSTHETA"
      real(8)                :: dx           = 0.1d0
      integer                :: nmolup       = 0

      ! Parser
      !
      integer                :: iopt, ierr

      ! Dummy
      !
      integer                :: itraj

      namelist /option_param/  &
        judgeup,               &
        zdef_type,             &
        nmolup,                &
        mode,                  & 
        orient_type,           &
        dx


      rewind io 
      read(io, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("judgeup     = ", a)')       trim(judgeup)
      write(iw,'("zdef_type   = ", a)')       trim(zdef_type)
      write(iw,'("nmolup      = ", i0)')      nmolup 
      write(iw,'("mode        = ", 3(a,2x))') trim(mode(1)), trim(mode(2)), trim(mode(3))
      write(iw,'("orient_type = ", a)')       trim(orient_type)
      write(iw,'("dx          = ", f15.7)')   dx 

      ! Parse
      !
      iopt = get_opt(judgeup, JudgeUpMode, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("judgeup = ",a," is not available.")') trim(judgeup)
        stop
      end if
      option%judgeup = iopt

      iopt = get_opt(zdef_type, ZdefTypes, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("zdef_type = ",a," is not available.")') trim(zdef_type)
        stop
      end if
      option%zdef_type = iopt

      do itraj = 1, 3 
        iopt = get_opt(mode(itraj), CoMMode, ierr)
        if (ierr /= 0) then
          write(iw,'("Read_Ctrl_Option> Error.")')
          write(iw,'("mode = ",a," is not available.")') trim(mode(itraj))
          stop
        end if
        option%mode(itraj) = iopt
      end do

      option%mode(4) = CoMModeWHOLE
      option%mode(5) = CoMModeWHOLE
      option%mode(6) = CoMModeWHOLE
      option%mode(7) = CoMModeWHOLE

      iopt = get_opt(orient_type, OrientTypes, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("xcoord = ",a," is not available.")') trim(orient_type)
        stop
      end if
      option%orient_type = iopt

      ! Send
      !
      option%dx       = dx
      option%nmolup   = nmolup

      if (option%zdef_type == ZdefTypeOUTERPROD) then
       option%nsel = 7 
      else
       option%nsel = 3 
      end if


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
