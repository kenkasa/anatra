!=======================================================================
module mod_output
!=======================================================================
  use mod_util
  use mod_const

  implicit none

  ! structures
  !
  type :: s_output
    integer                :: ndigit = 4
    character(len=MaxChar) :: fhead  
    character(len=MaxChar) :: file_extension
  end type s_output

  ! subroutines
  !
  public :: read_ctrl_output

  contains
!-----------------------------------------------------------------------
    subroutine read_ctrl_output(iunit, output, myrank)
!-----------------------------------------------------------------------
      implicit none
!
      integer,           intent(in)  :: iunit
      type(s_output),    intent(out) :: output
      integer, optional, intent(in)  :: myrank 

      integer                     :: ndigit         = 4 
      character(len=MaxChar)      :: fhead          = "run"
      character(len=MaxChar)      :: file_extension = "dat"

      integer :: irank

      namelist /output_param/ &
        ndigit,               &
        fhead,                &
        file_extension


      if (present(myrank)) then
        irank = myrank
      else
        irank = 0
      end if

      ! Initialize
      !
      ndigit         = 4 
      fhead          = "run"
      file_extension = "dat"

      ! Read Output namelist
      !
      rewind iunit
      read(iunit, output_param)

      ! Send
      !
      output%ndigit         = ndigit
      output%fhead          = fhead
      output%file_extension = file_extension


    end subroutine read_ctrl_output
!-----------------------------------------------------------------------
!

end module mod_output
!=======================================================================
