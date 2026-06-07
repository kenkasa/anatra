!=======================================================================
program main 
!=======================================================================
  use mod_const
  use mod_util
  use mod_ctrl
  use mod_analyze

  type(s_option)  :: option 

  call show_title
  call show_usage
  call read_ctrl(option)
  call analyze(option)
  call termination("rprate_analysis")

end program main 
!=======================================================================

!-----------------------------------------------------------------------
subroutine show_title
!-----------------------------------------------------------------------
  implicit none

  write(6,'("==================================================")')
  write(6,*)
  write(6,'("            Rate Constant Evaluation")')
  write(6,'("       Based on Returning Probability Theory")')
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
    write(iw,'("&option_param")')
    write(iw,'(" rporder       = 0")')
    write(iw,'(" use_sfe       = .false.")')
    write(iw,'(" timeunit      = 1.0d-9")')
    write(iw,'(" kins          = 0.2334")')
    write(iw,'(" kins_err      = 0.0200")')
    write(iw,'(" taud          = 3.9800")')
    write(iw,'(" taud_err      = 0.4000")')
    write(iw,*)
    write(iw,'(" Kstar         = 10.000    ! required if use_sfe = .false.")')
    write(iw,'(" Kstar_err     = 0.2300    ! required if use_sfe = .false.")')
    write(iw,*)
    write(iw,'(" temperature   = 300.0     ! required if use_sfe = .true.")')
    write(iw,'(" sfeb          = -11.03    ! required if use_sfe = .true.")')
    write(iw,'(" sfeb_err      =   0.10    ! required if use_sfe = .true.")')
    write(iw,'(" sfed          =  -9.20    ! required if use_sfe = .true.")')
    write(iw,'(" sfed_err      =   0.05    ! required if use_sfe = .true.")')
    write(iw,'("/")')
    stop
  end if

end subroutine show_usage
!-----------------------------------------------------------------------

