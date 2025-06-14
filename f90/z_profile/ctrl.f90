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
  integer,      parameter, public :: DensityTypeNUMBER   = 1
  integer,      parameter, public :: DensityTypeELECTRON = 2 
  character(*), parameter, public :: DensityType(2)      = (/'NUMBER   ',&
                                                             'ELECTRON '/)

  ! structures
  !
  type :: s_option
    real(8) :: dz         = 0.1d0
    integer :: mode(2)    = (/CoMModeRESIDUE, CoMModeRESIDUE/)
    integer :: denstype   = DensityTypeNUMBER
    integer :: symmetrize = .false. 
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

      real(8)                :: dz         = 0.1d0
      character(len=MaxChar) :: mode(2)    = (/'RESIDUE', 'RESIDUE'/)
      character(len=MaxChar) :: denstype   = 'NUMBER'
      character(len=MaxChar) :: centertype = 'ZERO'
      logical                :: symmetrize = .false.

      ! Parser
      !
      integer :: iopt, ierr

      ! Dummy
      !
      integer :: itraj

      namelist /option_param/  &
        dz,                    &
        mode,                  &
        denstype,              &
        symmetrize


      rewind io 
      read(io, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("dz         = ", f15.7)') dz 
      write(iw,'("mode       = ", 2(a,2x))') trim(mode(1)), trim(mode(2))
      write(iw,'("denstype   = ", a)')       trim(denstype)
      write(iw,'("symmetrize = ", a)')       get_tof(symmetrize)

      ! Parse
      !
      do itraj = 1, 2
        iopt = get_opt(mode(itraj), CoMMode, ierr)
        if (ierr /= 0) then
          write(iw,'("Read_Ctrl_Option> Error.")')
          write(iw,'("mode = ",a," is not available.")') trim(mode(itraj))
          stop
        end if
        option%mode(itraj) = iopt
      end do

      iopt = get_opt(denstype, DensityType, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("denstype = ",a," is not available.")') trim(denstype)
        stop
      end if
      option%denstype = iopt

      option%dz         = dz
      option%symmetrize = symmetrize

      ! Combination check
      !
      if (option%mode(1) /= ComModeATOM .and. &
          option%denstype == DensityTypeELECTRON ) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("Combination mode /= ATOM and denstype = ELECTIONS")')
        write(iw,'("is not available.")')
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
