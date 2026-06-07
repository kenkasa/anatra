!=======================================================================
program main 
!=======================================================================
  use mod_const
  use mod_ctrl
  use mod_analyze

  type(s_input)   :: input
  type(s_output)  :: output
  type(s_option)  :: option
  type(s_bootopt) :: bootopt

  call show_title
  call show_usage
  call read_ctrl(input, output, option, bootopt)
  call analyze(input, output, option, bootopt)
  call termination("Bootstrap prepper")

end program main 
!=======================================================================

!-----------------------------------------------------------------------
subroutine show_title
!-----------------------------------------------------------------------
  implicit none

  write(6,'("==================================================")')
  write(6,*)
  write(6,'("               Bootstrap Prepper")')
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

    !write(iw,'("&option_param")')
    !write(iw,'("/")')
    !write(iw,*)

    write(iw,'("&bootopt_param")')
    write(iw,'(" duplicate     = .true.        ! Whether to allow duplication when generating random numbers")')
    write(iw,'("                               ! default: .true.,")')
    write(iw,'(" iseed         =  3141592      ! input seed_number (if iseed <= 0, seed is generated)")')
    write(iw,'(" nsample       =  1000         ! Number of samples to be selected for each trial")')
    write(iw,'(" ntrial        =  100          ! Number of bootstrap trials")')
    write(iw,'(" ")')
    write(iw,'("/")')

    stop
  end if

end subroutine show_usage
!-----------------------------------------------------------------------
