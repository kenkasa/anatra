!=======================================================================
program main 
!=======================================================================
  use mod_const
  use mod_util
  use mod_ctrl
  use mod_analyze
  use mod_bootstrap

  type(s_input)   :: input
  type(s_output)  :: output
  type(s_option)  :: option 
  type(s_cvinfo)  :: cvinfo
  type(s_bootopt) :: bootopt

  call show_title
  call show_usage
  call read_ctrl(input, output, option, cvinfo, bootopt)
  call analyze(input, output, option, cvinfo, bootopt)
  call termination("rp_analysis")

end program main 
!=======================================================================

!-----------------------------------------------------------------------
subroutine show_title
!-----------------------------------------------------------------------
  implicit none

  write(6,'("==================================================")')
  write(6,*)
  write(6,'("          Returning Probability Analysis")')
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
    write(iw,'(" flist         = ""cvlist""    ! File that contains List of CV")')
    write(iw,'("/")')
    write(iw,'("&output_param")')
    write(iw,'(" fhead         = ""rpfilehead""! returning probability file")')
    write(iw,'("/")')
    write(iw,'("&option_param")')
    write(iw,'(" calc_rp       = .true.        ! calculate RP or not (default: .true.")')
    write(iw,'(" check_tscale  = .false.       ! check transient time scale (default: .true.)")')
    write(iw,'("                               ! if true, parameters ntr_sta, ntr_interval,")') 
    write(iw,'("                               ! and nntr should be specified.")')
    write(iw,'(" ntr_sta       =  10           ! minimum transient time scale size")')
    write(iw,'(" ntr_interval  =  100          ! interval of transient time scale checked for the convergence")')
    write(iw,'(" nntr          =  10           ! number of transient time scales checked for the convergence")')
    write(iw,'(" ndim          =  1            ! CV dimensions")')
    write(iw,'(" nt_range      =  1000         ! time range of RP")')
    write(iw,'(" nt_transient  =  100          ! transient time scale")')
    write(iw,'(" dt            =  0.001d0      ! time grid spacing")')
    write(iw,'(" use_moving    = .false.       ! use moving average scheme for detecting transition events")')
    write(iw,'(" use_bootstrap = .false.       ! use boot-strap analyze (defaule:.false.)")')
    write(iw,'(" bound_range   = -10.0   0.0   ! range of bound state (1st: minimum, 2nd: maximum)")')
    write(iw,'(" react_range   =  0.0   10.0   ! range of reactive state (1st: minimum, 2nd: maximum)")')
    write(iw,'(" unbound_range =  10.0  20.0   ! range of unbound state (1st: minimum, 2nd: maximum")')
    write(iw,'(" ")')
    write(iw,'("/")')
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

