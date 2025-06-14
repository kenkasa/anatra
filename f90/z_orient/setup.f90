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
    subroutine setup(input, option, trajopt, traj)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),   intent(in)  :: input
      type(s_option),  intent(in)  :: option 
      type(s_trajopt), intent(in)  :: trajopt
      type(s_traj),    intent(out) :: traj(7)
  
      type(s_dcd)                  :: dcd

      integer :: itraj


      if (input%ftraj(1) == "") then
        write(iw,'("Setup> Error.")')
        write(iw,'("ftraj in input_param is empty...")')
        stop
      end if

      if (option%zdef_type == ZdefTypeCoord) then
        do itraj = 1, 3 
          call setup_traj_from_args(trajopt,      &
                                    1,            &
                                    traj(itraj),  &
                                    trajid = itraj)
        end do
      else ! OUTERPROD
        do itraj = 1, 7 
          call setup_traj_from_args(trajopt,      &
                                    1,            &
                                    traj(itraj),  &
                                    trajid = itraj)
        end do
      end if

    end subroutine setup

!-----------------------------------------------------------------------

end module mod_setup
!=======================================================================
