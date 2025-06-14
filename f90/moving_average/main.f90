!=======================================================================
program main 
!=======================================================================
  use mod_const
  use mod_movave_ctrl
  use mod_movave_analyze

  type(s_input)               :: input
  type(s_output)              :: output
  type(s_movave_option)       :: option
  type(s_movave), allocatable :: movave(:)

  call show_title
  call show_usage
  call read_movave_ctrl(input, output, option, movave)
  call movave_analyze(input, output, option, movave)
  call movave_write(output, movave)

end program main 
!=======================================================================

!-----------------------------------------------------------------------
subroutine show_title
!-----------------------------------------------------------------------
  implicit none

  write(6,'("==================================================")')
  write(6,*)
  write(6,'("             Moving Average Analysis")')
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
    write(iw,'(" fhead = ""out""    ! Output file header")')
    write(iw,'("/")')

    write(iw,'("&option_param")')
    write(iw,'("  dt           = 0.1d0   ! time grid")')
    write(iw,'("  t_sta        = 0.0d0   ! minimum value of time")')
    write(iw,'("  t_sep        = 1.0d0   ! boundary of time domains (if no boundary (nregion = 1), &
              &                          ! specify t_sep = 0.0d0)")')
    write(iw,'("  nregion      = 2       ! # of diffrent regions")')
    write(iw,'("  npoint       = 3 5     ! # of data points used for each region (odd number)")')
    write(iw,'("  include_zero = .true.  ! whether origin is changed or not by the average")')
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
