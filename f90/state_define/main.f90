!=======================================================================
program main 
!=======================================================================
  use mod_const
  use mod_grid3d
  use mod_ctrl
  use mod_cv
  use mod_analyze

  type(s_input)   :: input
  type(s_output)  :: output
  type(s_option)  :: option

  call show_title
  call show_usage
  call read_ctrl(input, output, option)

  if (option%use_voronoi) then
    call determine_voronoi_state(input, output, option)
  end if

  if (option%use_grids) then
    call determine_spatial_state(input, output, option)
  end if

  call termination('st_analysis')

end program main 
!=======================================================================

!-----------------------------------------------------------------------
subroutine show_title
!-----------------------------------------------------------------------
  implicit none

  write(6,'("==================================================")')
  write(6,*)
  write(6,'("                 State Define")')
  write(6,*)
  write(6,'("==================================================")')
  write(6,'("Remark: Currently, only Voronoi analysis is supported.")')

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
    write(iw,'(" flist = ""cvlist""   ! File that contains List of CV")')
    write(iw,'(" fdx   = ""input.dx"" ! used if use_grids = .true.")')
    write(iw,'("/")')

    write(iw,'("&output_param")')
    write(iw,'(" fhead = ""out""    ! Header of Output iles")')
    write(iw,'("/")')

    write(iw,'("&option_param")')
    write(iw,'("  ndim        = 2")')
    write(iw,'("  xyzcol      = 1 2  !&
              &column ID corresponding to x, y, z components in cv file")')
    write(iw,'("                            &
              &! (please ignore the first column in cv file)")')
    write(iw,'("  use_voronoi = .true.")')
    write(iw,'("  vr_ncell    = 3           ! # of voronoi cells")')
    write(iw,'("  vr_normvec  = 1.0d0 1.0d0 ! normalization vector")')
    write(iw,'("  vr_cellpos  = ")')
    write(iw,'("    1.0d0 2.0d0")')
    write(iw,'("    5.0d0 7.0d0")')
    write(iw,'("   10.0d0 2.0d0")')
    write(iw,'("  use_grids   = .false.")')
    write(iw,'("  pbc         = .false. ! pbc adopted (active only if use_grids = .true.")')
    write(iw,'("  box         = 50.0 50.0 50.0 ! box size (active only if use_grids = .true.")')
    write(iw,'("/")')
    write(iw,*)
    stop
  end if

end subroutine show_usage
!-----------------------------------------------------------------------
