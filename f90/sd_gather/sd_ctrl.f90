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
  type :: s_sdinfo
    integer                             :: nfile     = 1
    !integer                             :: nstep_tot = 0
    real(8)                             :: norm_tot  = 0.0d0
    character(len=MaxChar), allocatable :: fsdbin(:)
  end type s_sdinfo

  type :: s_option
    real(8) :: norm_tot  = 0.0d0
    !integer :: nstep_tot = 0
  end type s_option

  ! subroutines
  !
  public  :: read_ctrl
  private :: read_ctrl_option
  private :: read_ctrl_sdinfo

  contains

!-----------------------------------------------------------------------
    subroutine read_ctrl(input, output, option, sdinfo)
!-----------------------------------------------------------------------
      implicit none

      integer, parameter           :: iunit = 10

      type(s_input),   intent(out) :: input
      type(s_output),  intent(out) :: output
      type(s_option),  intent(out) :: option
      type(s_sdinfo),  intent(out) :: sdinfo

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
      close(iunit)

      open(iunit, file=trim(input%flist))
        call read_ctrl_sdinfo  (iunit, sdinfo)
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

      real(8)                :: norm_tot
      !integer                :: nstep_tot 
      integer                :: iopt, ierr

      namelist /option_param/ norm_tot !, nstep_tot 


      norm_tot = -1.0d0 
      !nstep_tot = 0

      rewind iunit
      read(iunit, option_param, end=101)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("norm_tot      = ", i0)')  norm_tot
      !write(iw,'("nstep_tot     = ", i0)')  nstep_tot 

      !do itraj = 1, 2
      !  iopt = get_opt(mode(itraj), CoMMode, ierr)
      !  if (ierr /= 0) then
      !    write(iw,'("Read_Ctrl_Option> Error.")')
      !    write(iw,'("mode = ",a," is not available.")') trim(mode(itraj))
      !    stop
      !  end if
      !  option%mode(itraj) = iopt
      !end do

101   continue 
      option%norm_tot = norm_tot
      !option%nstep_tot  = nstep_tot

      if (option%norm_tot < 0.0d0) then
      !if (option%nstep_tot == 0) then
        write(iw,'("Read_Ctrl_Option> norm_tot is not specified.")')
        write(iw,'(">> total norm will be read from sdbin files")')
        write(iw,*)
      end if

      ! Combination check
      !

    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

!
!-----------------------------------------------------------------------
    subroutine read_ctrl_sdinfo(iunit, sdinfo)
!-----------------------------------------------------------------------
      implicit none
!
      integer,        intent(in)  :: iunit
      type(s_sdinfo), intent(out) :: sdinfo


      integer :: ifile, nfile

      rewind iunit

      ifile = 0
      do while(.true.)
        read(iunit,*,end=100)
        ifile = ifile + 1
      end do
100   continue
      rewind iunit
      nfile = ifile
      allocate(sdinfo%fsdbin(ifile))
      
      do ifile = 1, nfile
        read(iunit,'(a)') sdinfo%fsdbin(ifile)
      end do

      write(iw,*)
      write(iw,'(">> RDINFO parameters")')
      write(iw,'("nfile    = ", i0)')  nfile
      do ifile = 1, nfile
        write(iw,'(a)') trim(sdinfo%fsdbin(ifile))
      end do

      sdinfo%nfile = nfile

    end subroutine read_ctrl_sdinfo
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end module
!=======================================================================
