!=======================================================================
program main 
!=======================================================================
  use mod_const
  use mod_ctrl
  use mod_setup
  use mod_dcdio
  use mod_traj
  use mod_analyze

  type(s_input)    :: input
  type(s_output)   :: output
  type(s_option)   :: option
  type(s_timegrid) :: timegrid 
  type(s_trajopt)  :: trajopt
  type(s_traj)     :: traj

  call show_title
  call read_ctrl(input, output, option, trajopt, timegrid)
  call setup(input, trajopt, traj)
  call analyze(input, output, option, timegrid, traj)
  call dealloc_traj(traj)
  call termination("Center of Mass coordinate analysis")

end program main 
!=======================================================================

!-----------------------------------------------------------------------
subroutine show_title
!-----------------------------------------------------------------------
  implicit none

  write(6,'("==================================================")')
  write(6,*)
  write(6,'("     Center-of-Mass Coordinate Analysis")')
  write(6,*)
  write(6,'("==================================================")') 

end subroutine show_title
!-----------------------------------------------------------------------
