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
  integer,      parameter :: ParmFormatANAPARM   = 1
  integer,      parameter :: ParmFormatPRMTOP    = 2
  character(*), parameter :: ParmFormatTypes(2)  = (/'ANAPARM   ',& 
                                                     'PRMTOP    '/)

  integer,      parameter :: CoordTypeZ        = 1
  character(*), parameter :: CoordTypes(1)     = (/'Z       '/)


  ! structures
  !
  type :: s_option
    integer :: parmformat       = ParmFormatPRMTOP
    integer :: mode(2)          = (/ComModeWHOLE, ComModeWHOLE/)
    integer :: coord_type       = CoordTypeZ 
    real(8) :: xsta             = 0.0d0
    real(8) :: dx               = 1.0d0
    integer :: ngrid            = 10
    integer :: nins             = 100000

    ! Hidden options
    !
    logical :: output_param     = .false.

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

      character(len=MaxChar) :: parmformat = "PRMTOP" 
      character(len=MaxChar) :: mode(2)    = (/"WHOLE", "WHOLE"/) ! Not in Namelist
      character(len=MaxChar) :: coord_type = "Z"
      real(8)                :: xsta       = 0.0d0
      real(8)                :: dx         = 1.0d0
      integer                :: ngrid      = 10
      integer                :: nins       = 1000

      ! Hidden options
      logical                :: output_param     = .false.

      integer                :: i
      integer                :: iopt, ierr

      namelist /option_param/ parmformat, &
                              coord_type, &
                              xsta,       &
                              dx,         &
                              ngrid,      &
                              nins


      rewind iunit
      read(iunit, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      !write(iw,'("parmformat       = ", a)')       trim(parmformat) 
      write(iw,'("coord_type       = ", a)')       trim(coord_type)
      write(iw,'("xsta             = ",f20.10)')   xsta 
      write(iw,'("dx               = ",f20.10)')   dx 
      write(iw,'("ngrid            = ",i0)')       ngrid
      write(iw,'("nins             = ",i0)')       nins 

      iopt = get_opt(parmformat, ParmFormatTypes, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("parmformat = ",a," is not available.")') trim(parmformat)
        stop
      end if
      option%parmformat = iopt

      iopt = get_opt(coord_type, CoordTypes, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("coord_type = ",a," is not available.")') trim(coord_type)
        stop
      end if
      option%coord_type = iopt

      option%xsta         = xsta 
      option%dx           = dx
      option%ngrid        = ngrid
      option%nins         = nins
      option%output_param = output_param
      option%mode(1:2)    = (/ComModeWHOLE, ComModeWHOLE/)

      ! Combination check
      !

    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end module
!=======================================================================
