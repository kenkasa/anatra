!=======================================================================
module mod_setup
!=======================================================================
  use mod_const
  use mod_ctrl
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

      ! Dummy
      !
      integer :: i

      if (input%ftraj(1) == "") then
        write(iw,'("Setup> Error.")')
        write(iw,'("ftraj in input_param is empty...")')
        stop
      end if

      do i = 1, 2
        call setup_traj_from_args(trajopt, 1, traj(i), trajid = i)
      end do


    end subroutine setup
!-----------------------------------------------------------------------

end module mod_setup
!=======================================================================
