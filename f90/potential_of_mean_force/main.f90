!=======================================================================
program main 
!=======================================================================
  use mod_const
  use mod_grid3d
  use mod_ctrl
  use mod_analyze

  type(s_input)   :: input
  type(s_output)  :: output
  type(s_option)  :: option
  type(s_bootopt) :: bootopt
  type(s_mpl2d)   :: mpl2d

  call show_title
  call show_usage
  call read_ctrl(input, output, option, mpl2d, bootopt)

  if (option%use_bootstrap) then
    !call analyze_bootstrap(input, output, option, mpl2d, bootopt)
  else
    call analyze(input, output, option, mpl2d)
  end if
  call termination('free_energy')

end program main 
!=======================================================================

!-----------------------------------------------------------------------
subroutine show_title
!-----------------------------------------------------------------------
  implicit none

  write(6,'("==================================================")')
  write(6,*)
  write(6,'("              Free-Energy Analysis")')
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
    write(iw,'(" flist = ""cvlist"" ! File that contains List of CV")')
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
    write(iw,'("  gen_script_mpl2d = .false.   ! whether generate python script or not (only ndim = 2 is supported)")')
    write(iw,'("  skip_calc = .false.       ! if .true., calculation is skipped")')
    write(iw,'("  ")')
    write(iw,'("  spline_resolution = 4     ! make spline fine as increasing this parameter ")')
    write(iw,'("  use_voronoi       = .false.     ! use voronoi tessellation")')
    write(iw,'("  use_cellsearch    = .false.     ! search voronoi cells ")')
    write(iw,'("                                  ! # of cells is specified by vr_ncell")')
    write(iw,'("  vr_ncell          = 3           ! # of voronoi cells ")')
    write(iw,'("  vr_separation     = 1.0d0 2.0d0 ! separation between cells")')
    write(iw,'("  vr_normvec        = 4.0d0 5.0d0 ! normalization vector of coordinates")')
    write(iw,'("  vr_cellpos        = 1.0d0 2.0d0")')
    write(iw,'("                      5.0d0 7.0d0")')
    write(iw,'("                      9.0d0 1.0d0")')
    write(iw,'("  ! cell positions.&
                  & if use_cellsearch = .false.,&
                  & user should specify the positions with vr_cellpos")')
    write(iw,'("/")')
    write(iw,*)
    write(iw,'("! if gen_script_mpl2d = .true., following section should be specified.")')
    write(iw,'("&matplotlib_param")')
    write(iw,'("  fpdf   = ""pmf.pdf""")')
    write(iw,'("  labels = ""$x$"" ""$y$"" ""PMF""     ! labels for x, y, and z axes")')
    write(iw,'("  ranges = -10.0 10.0 -5.0 5.0 0.0 10  ! plot ranges for x, y, and z axes")')
    write(iw,'("  tics   = 1.0 1.0 0.5                 ! tics for x, y, zna z axes")')
    write(iw,'("  scales = 1.0 1.0 1.0 ! unit conversion constants for x, y, and z axes")')
    write(iw,'("/")')
    stop
  end if

end subroutine show_usage
!-----------------------------------------------------------------------
