!=======================================================================
program main 
!=======================================================================
  use mod_const
  use mod_ctrl
  use mod_setup
  use mod_dcdio
  use mod_traj
  use mod_setup
  use mod_cv
  use mod_analyze

  type(s_input)   :: input
  type(s_output)  :: output
  type(s_option)  :: option
  type(s_trajopt) :: trajopt
  type(s_traj)    :: traj(3)
  type(s_cvinfo)  :: cvinfo

  call show_title
  call show_usage
  call read_ctrl(input, output, option, trajopt, cvinfo)
  call setup(input, trajopt, traj)
  call analyze(input, output, option, traj, cvinfo)
  call dealloc_traj(traj(1))
  call dealloc_traj(traj(2))
  call dealloc_traj(traj(3))

end program main 
!=======================================================================

!-----------------------------------------------------------------------
subroutine show_title
!-----------------------------------------------------------------------
  implicit none

  write(6,'("==================================================")')
  write(6,*)
  write(6,'("         Spatial-Distribution Analysis")')
  write(6,*)
  write(6,'("==================================================")') 

end subroutine show_title
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine show_usage
!-----------------------------------------------------------------------
  use mod_const

  implicit none

  character(len=MaxChar) :: f_ctrl


  call getarg(1, f_ctrl)

  if (trim(f_ctrl) == "-h") then
    write(iw,'("&input_param")')
    write(iw,'(" ftraj   = ""inp.dcd""     ! input trajectory file")')
    write(iw,'(" fcv     = ""inp.cv""      ! input CV file (neccessary if use_conditional = .true.)")')
    write(iw,'(" fweight = ""inp.weight""  ! input weight file (weight is considered if this file is specified))")')
    write(iw,'(" !!! Alternative !!!")')
    write(iw,'(" flist_traj   = ""flist.traj""   ")')
    write(iw,'(" flist_cv     = ""flist.cv""     ")')
    write(iw,'(" flist_weight = ""flist.weight"" ")')
    write(iw,'("/")')

    write(iw,'("&output_param")')
    write(iw,'(" fhead = ""out""    ! Header of Output iles")')
    write(iw,'("/")')

    write(iw,'("&option_param")')
    write(iw,'("  ndim        = 2")')
    write(iw,'("  xyzcol      = 1 2  ! column ID corresponding to x, y, z components in cv file")')
    write(iw,'("                     ! (please ignore the first column in cv file)")')
    write(iw,'("  ng3         = 10 20       ! # of grids for each cv")')
    write(iw,'("  del         = 0.1d0 0.2d0 ! grid spacing for each cv")')
    write(iw,'("  origin      = 0.0d0 1.0d0 ! origin for each cv")')
    write(iw,'("  temperature = 298.0d0     ! temperature [K]")')
    write(iw,'("  box_system  = 50.0d0 50.0d0 50.0d0 ! system box size (used if ndim = 3)")')
    write(iw,'("  nstep_tot   = 500000      ! totan number of steps (used if ndim = 3)")')
    write(iw,'()')
    write(iw,'("  use_spline  = .false.     ! whether spline is used or not")')
    write(iw,'("  spline_resolution = 4     ! make spline fine as increasing this parameter ")')
    write(iw,'("/")')
    write(iw,*)
    stop
  end if

end subroutine show_usage
!-----------------------------------------------------------------------
