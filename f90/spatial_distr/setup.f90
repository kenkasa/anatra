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
      type(s_traj),    intent(out) :: traj(3)

      integer :: itraj, nmolinfo  

      if (input%ftraj(1) == "") then
        write(iw,'("Setup> Error.")')
        write(iw,'("ftraj in input_param is empty...")')
        stop
      end if

      nmolinfo = trajopt%nmolinfo

      do itraj = 1, nmolinfo
        call setup_traj_from_args(trajopt,      &
                                  1,            &
                                  traj(itraj),  &
                                  trajid = itraj)
      end do

      if (nmolinfo == 2) then
        call setup_traj_from_args(trajopt,      &
                                  1,            &
                                  traj(3),      &
                                  trajid = 2)
      end if
      
    end subroutine setup

!-----------------------------------------------------------------------

end module mod_setup
!=======================================================================
