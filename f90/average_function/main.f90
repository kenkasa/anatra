!=======================================================================
program main 
!=======================================================================
  use mod_const
  use mod_ctrl
  use mod_analyze
  use mod_bootstrap

  type(s_input)   :: input
  type(s_output)  :: output
  type(s_option)  :: option 
  type(s_bootopt) :: bootopt

  call show_title
  call show_usage
  call read_ctrl(input, output, option, bootopt)

  if (option%use_bootstrap) then
    call analyze_bootstrap(input, output, option, bootopt)
  else
    call analyze(input, output, option)
  end if

  call termination("Average function analysis")

end program main 
!=======================================================================

!-----------------------------------------------------------------------
subroutine show_title
!-----------------------------------------------------------------------
  implicit none

  write(6,'("==================================================")')
  write(6,*)
  write(6,'("          Average Function Analysis")')
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
    write(iw,'(" fhead         = ""filehead""  ! header of output file")')
    write(iw,'("/")')
    write(iw,'("&option_param")')
    write(iw,'(" use_bootstrap =  .false.      ! use bootstrap or not")')
    write(iw,'(" xsta          =  0.0d0        ! origin of grid")')
    write(iw,'(" dx            =  0.1d0        ! grid spacing")')
    write(iw,'("/")')
    write(iw,'(" ")')
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

