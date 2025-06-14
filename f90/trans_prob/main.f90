!=======================================================================
program main 
!=======================================================================
  use mod_const
  use mod_ctrl
  use mod_analyze
  use mod_bootstrap

  type(s_input)    :: input
  type(s_output)   :: output
  type(s_option)   :: option 
  type(s_bootopt)  :: bootopt
  type(s_timegrid) :: timegrid 

  call show_title
  call show_usage
  call read_ctrl(input, output, option, bootopt, timegrid)

  if (option%use_bootstrap) then
    !call analyze_bootstrap(input, output, option, bootopt)
  else
    call analyze(input, output, option, timegrid)
  end if

  call termination("Transition Probability Analysis")

end program main 
!=======================================================================

!-----------------------------------------------------------------------
subroutine show_title
!-----------------------------------------------------------------------
  implicit none

  write(6,'("==================================================")')
  write(6,*)
  write(6,'("        Transition Probability Analysis")')
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
    write(iw,'(" use_bootstrap        =  .false.       ! use bootstrap or not")')
    write(iw,'(" use_reflection_state =  .false.       ! not implemented")')
    write(iw,'(" use_product_state    =  .false.       ! not implemented")')
    write(iw,'(" kinetic_mode         = ""TRANSITION""   ! TRANSITION or REACTION")')
    write(iw,'(" nstate               =  4             ! # of states")')
    write(iw,'(" ndim                 =  1             ! # of dimensions")')
    write(iw,'(" nmol                 =  1             ! # of target molecules")')
    write(iw,'(" dt                   =  0.1d0         ! Time grid")')
    write(iw,'(" t_sparse             =  0.1d0         ! Sparse time-grid (used for output)")')
    write(iw,'(" t_range              =  10.0          ! Timescale for outputting TCFs")')
    write(iw,'(" initial_state_ids    =  2 3           ! Mandatory")')
    write(iw,'(" reflection_state_ids =  4             ! Optional")')
    write(iw,'(" product_state_ids    =  1             ! Optional")')
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
    write(iw,'("&state")')
    write(iw,'("-500.0  -16.0  0.5  ! State 1 (weight of 0.5)")')
    write(iw,'(" -16.0  -15.0  0.5  ! State 2 (weight of 0.5)")')
    write(iw,'(" -15.0   -2.0  0.0  ! State 3")')
    write(iw,'("  -2.0  500.0  0.0  ! State 4")')
    write(iw,'("/")')
    stop
  end if


end subroutine show_usage
!-----------------------------------------------------------------------

