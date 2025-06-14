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
  call analyze(input, output, option)
  call termination("Histogram analysis")

end program main 
!=======================================================================

!-----------------------------------------------------------------------
subroutine show_title
!-----------------------------------------------------------------------
  implicit none

  write(6,'("==================================================")')
  write(6,*)
  write(6,'("              Histogram Analysis")')
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
    write(iw,'(" fcv      = ""cvdata1"" ""cvdata2"" ! Time-series CV data")')
    write(iw,'(" flist_cv = ""list"" ! CV-file list (Either fcv or flist_cv should be specified)")')
    write(iw,'("/")')

    write(iw,'("&output_param")')
    write(iw,'(" fhead = ""run"" ! header of output filename")')
    write(iw,'("/")')

    write(iw,'("&option_param")')
    write(iw,'("  dx       = 0.1d0  ! grid spacing")')
    write(iw,'("  xsta     = 0.0d0  ! minimum value of the coordinate")')
    write(iw,'("  nx       = 100    ! number of grids")')
    write(iw,'("  ncol     = 1      ! number of data column")')
    write(iw,'("/")')
    write(iw,*)

    write(iw,*)
    write(iw,'("::: format of time-series data :::")')
    write(iw,'("1st column: time (arbitrary), Nth column (N>1): data")')
    write(iw,'("example)")')
    write(iw,*)
    write(iw,'("0.00   0.02120")')
    write(iw,'("0.01   0.03140")')
    write(iw,'("0.02   0.22310")')
    write(iw,'("....   .......")')
    stop
  end if

end subroutine show_usage
!-----------------------------------------------------------------------
