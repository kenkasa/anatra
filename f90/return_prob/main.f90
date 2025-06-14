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
  type(s_bootopt) :: bootopt

  call show_title
  call show_usage
  call read_ctrl(input, output, option, bootopt)
  call analyze(input, output, option, bootopt)
  call termination("sr_analysis")

end program main 
!=======================================================================

!-----------------------------------------------------------------------
subroutine show_title
!-----------------------------------------------------------------------
  implicit none

  write(6,'("==================================================")')
  write(6,*)
  write(6,'("           State Relaxation Analysis")')
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
    write(iw,'(" flist                = ""cvlist""  ! File that contains List of CV")')
    write(iw,'("/")')
    write(iw,'("&output_param")')
    write(iw,'(" fhead                = ""out""    ! header of output files")')
    write(iw,'("/")')
    write(iw,'("&option_param")')
    write(iw,'(" use_bootstrap        = .false.  ! use boot-strap analyze (defaule:.false.)")')
    write(iw,'(" check_timescale      = .false.  ! check transient time scale (default: .true.)")')
    write(iw,'("                                 ! if true, parameters t_transient_sta, t_transient_interval,")') 
    write(iw,'("                                 ! and n_check_transient should be specified.")')
    write(iw,'(" use_quench           = .false.  ! whether the snapshots appearing after entering quench zone is discard")')
    write(iw,'("                                 ! for transition calc. (default: .false.)")')
    write(iw,'(" use_reactraj         = .false.  ! whether reactive trajectory is used for transition calc. (default: .false.)")')
    write(iw,'(" out_normfactor       = .true.   ! whether normalized factor of Pret is outputted or not. ")')
    write(iw,'(" t_transient_sta      = 10.0     ! minimum transient time scale")')
    write(iw,'(" t_transient_interval =  100     ! interval of transient time scale checked for the convergence")')
    write(iw,'(" n_check_transient    =  100     ! number of transient time scales checked for the convergence")')
    write(iw,'(" nmol                 =  1       ! number of ligand molecules")')
    write(iw,'(" ndim                 =  1       ! CV dimensions")')
    write(iw,'(" nstate               =  4       ! # of states")')
    write(iw,'(" bound_id             =  1       ! ID of bound state (each state is defined in &state section")')
    write(iw,'(" reaczone_id          =  1       ! ID of reactive state (each state is defined in &state section)")')
    write(iw,'(" t_transient          =  0.0     ! transient time scale")')
    write(iw,'(" t_cut                =  0.0     ! trajectory length used for transition calc.")')
    write(iw,'(" t_range              =  1.0     ! time range of transition calc.")')
    write(iw,'(" dt                   =  0.001d0 ! time grid spacing")')
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
    write(iw,*)
    write(iw,'("&state")')
    write(iw,'("-500.0d0  -16.0d0   ! State 1 (In this sample, bound state (bound_id))")')
    write(iw,'(" -16.0d0  -15.0d0   ! State 2 (In this sample, reaction zone (reaczone_id))")')
    write(iw,'(" -15.0d0   -2.0d0   ! State 3")')
    write(iw,'("  -2.0d0  500.0d0   ! State 4")')
    write(iw,'("  -2.0d0  500.0d0   ! Quench zone")')
    write(iw,'("/")')
    stop
  end if

end subroutine show_usage
!-----------------------------------------------------------------------

