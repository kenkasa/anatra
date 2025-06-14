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
  type :: s_rdinfo
    integer                             :: nfile     = 1
    character(len=MaxChar), allocatable :: frdbin(:)
  end type s_rdinfo

  ! subroutines
  !
  public  :: read_ctrl
  private :: read_ctrl_rdinfo

  contains

!-----------------------------------------------------------------------
    subroutine read_ctrl(input, output, rdinfo)
!-----------------------------------------------------------------------
      implicit none

      integer, parameter           :: iunit = 10

      type(s_input),   intent(out) :: input
      type(s_output),  intent(out) :: output
      type(s_rdinfo),  intent(out) :: rdinfo

      character(len=MaxChar)       :: f_ctrl

      ! get control file name
      !
      call getarg(1, f_ctrl)

      write(iw,*)
      write(iw,'("Read_Ctrl> Reading parameters from ", a)') trim(f_ctrl)
      open(iunit, file=trim(f_ctrl), status='old')
        call read_ctrl_input  (iunit, input)
        call read_ctrl_output (iunit, output)
      close(iunit)

      open(iunit, file=trim(input%flist))
        call read_ctrl_rdinfo  (iunit, rdinfo)
      close(iunit) 

    end subroutine read_ctrl
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine read_ctrl_rdinfo(iunit, rdinfo)
!-----------------------------------------------------------------------
      implicit none
!
      integer,        intent(in)  :: iunit
      type(s_rdinfo), intent(out) :: rdinfo


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
      allocate(rdinfo%frdbin(ifile))
      
      do ifile = 1, nfile
        read(iunit,'(a)') rdinfo%frdbin(ifile)
      end do

      write(iw,*)
      write(iw,'(">> RDINFO parameters")')
      write(iw,'("nfile    = ", i0)')  nfile
      do ifile = 1, nfile
        write(iw,'(a)') trim(rdinfo%frdbin(ifile))
      end do

      rdinfo%nfile = nfile

    end subroutine read_ctrl_rdinfo
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end module
!=======================================================================
