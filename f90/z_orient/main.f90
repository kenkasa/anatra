!=======================================================================
program main 
!=======================================================================
  use mod_const
  use mod_ctrl
  use mod_setup
  use mod_dcdio
  use mod_traj
  use mod_analyze

  type(s_input)   :: input
  type(s_output)  :: output
  type(s_option)  :: option
  type(s_trajopt) :: trajopt
  type(s_traj)    :: traj(7)

  call show_title
  call read_ctrl(input, output, option, trajopt)
  call setup(input, option, trajopt, traj)
  call analyze(input, output, option, traj)
  call dealloc_traj(traj(1))
  call dealloc_traj(traj(2))
  call dealloc_traj(traj(3))
  call termination("Z-orient analysis")

end program main 
!=======================================================================

!-----------------------------------------------------------------------
subroutine show_title
!-----------------------------------------------------------------------
  implicit none

  write(6,'("==================================================")')
  write(6,*)
  write(6,'("               Z-Orient Analysis")')
  write(6,*)
  write(6,'("==================================================")') 

end subroutine show_title
!-----------------------------------------------------------------------
