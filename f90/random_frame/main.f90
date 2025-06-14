!=======================================================================
program main 
!=======================================================================
  use mod_const
  use mod_ctrl
  use mod_analyze

  type(s_input)   :: input
  type(s_output)  :: output
  type(s_option)  :: option 


  call show_title
  call show_usage
  call read_ctrl(input, output, option)

  if (option%shuffle) then
    call analyze_shuffle(input, output, option)
  else
    call analyze(input, output, option)
  end if
  call termination("Random-frame analysis")


end program main 
!=======================================================================

!-----------------------------------------------------------------------
subroutine show_title
!-----------------------------------------------------------------------
  implicit none

  write(6,'("==================================================")')
  write(6,*)
  write(6,'("              Random-Frame Analysis")')
  write(6,*)
  write(6,'("==================================================")') 

end subroutine show_title
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine show_usage
!-----------------------------------------------------------------------
  use mod_const

  implicit none

  character(len=MaxChar) :: f_ctrl

  
  call getarg(1, f_ctrl)

  if (trim(f_ctrl) == "-h") then
    write(iw,'("&input_param")')
    write(iw,'(" ftraj      = ""inp.dcd"" ! input trajectories")')
    write(iw,'(" flist_traj = ""flist""   ! trajectory list neccesary if ftraj is not specified")')
    write(iw,'("/")')
    write(iw,*)
    write(iw,'("&output_param")')
    write(iw,'(" fhead = ""out"" ! header of output filename")')
    write(iw,'("/")')
    write(iw,*)
    write(iw,'("&option_param")')
    write(iw,'(" nsample         = 100       ! # of frames randomly extracted")')
    write(iw,'(" iseed           = 3141592   ! random seed")')
    write(iw,'(" output_trajtype = 3141592   ! output trajectory format (dcd or xtc)")')
    write(iw,'(" duplicate       = .false.   ! whether to allow the duplicated extraction")')
    write(iw,'(" out_rst7        = .false.   ! whether to output AMBER rst7 format file")')
    write(iw,'(" shuffle         = .false.   ! whether to shuffle frames (if true, all the frames are sampled")')
    write(iw,'("/")')
    stop
  end if


end subroutine show_usage
!-----------------------------------------------------------------------

