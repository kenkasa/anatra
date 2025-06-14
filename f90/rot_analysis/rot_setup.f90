!=======================================================================
module mod_setup
!=======================================================================
  use mod_const
  use mod_ctrl
  use mod_dcdio
  use mod_traj
  implicit none


  ! subroutines
  !
  public :: setup

  contains
!-----------------------------------------------------------------------
    subroutine setup(input, trajopt, traj)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),   intent(in)  :: input
      type(s_trajopt), intent(in)  :: trajopt
      type(s_traj),    intent(out) :: traj(2)
  
      type(s_dcd)                  :: dcd

      integer :: itraj

      do itraj = 1, 2
        call read_dcd(input%fdcd(itraj), 0, dcd)
        call setup_traj(trajopt, dcd, traj(itraj), itraj)
        call dealloc_dcd(dcd)
      end do

    end subroutine setup

!-----------------------------------------------------------------------

end module mod_setup
!=======================================================================
