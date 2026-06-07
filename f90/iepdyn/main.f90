!=======================================================================
program main 
!=======================================================================
  use mod_const
  use mod_ctrl
  use mod_analyze

  type(s_input)       :: input
  type(s_extra_input) :: einput
  type(s_output)      :: output
  type(s_option)      :: option 
  type(s_timegrid)    :: timegrid 

  call show_title
  call show_usage
  call read_ctrl  (input, einput, output, option, timegrid)
  call analyze    (input, einput, output, option, timegrid)
  call termination("IEPDYN Analysis")

end program main 
!=======================================================================

!-----------------------------------------------------------------------
subroutine show_title
!-----------------------------------------------------------------------
  implicit none

  write(6,'("========================================================")')
  write(6,*)
  write(6,'("                    I  E  P  D  Y  N")')
  write(6,*)
  write(6,'("   Integral-Equation formalism of Population DYNamics")')
  write(6,*)
  write(6,'("========================================================")')
  write(6,'("[Developer]")')
  write(6,'("Kento Kasahara (The Univ. of Osaka)")')
  write(6,*)
  write(6,'("[Reference]")')
  write(6,'("K. Kasahara, R. Okabe, C. A. Chang, T. Mori, and N. Matubayasi,")')
  write(6,'("IEPDYN: Integral-equation formalism of population dynamics,")')
  write(6,'("J. Chem. Phys., 164, 124112 (2026).")')
  write(6,*)

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
    write(iw,'(" fhead         = ""output""  ! header of output file")')
    write(iw,'("/")')
    write(iw,'("&option_param")')
    write(iw,'(" input_type             =  ""TIMESERIES"" ! Input fcv format (Default: TIMESERIES)")')
    write(iw,'(" ! TIMESERIES or HISTOGRAM")')
    write(iw,'(" ! TIMESERIES format is ANATRA standard")')
    write(iw,'(" ! HISTOGRAM-formatted is restart file created with this program when output_histgram = .true.")')
    write(iw,*)
    write(iw,'(" output_histogram       =  .false.       ! create restart files (krhist for Kijk and Rij)")')
    write(iw,'(" use_dissociate_state   =  .false.       ! define dissociate state or not")')
    write(iw,'(" use_reflection_state   =  .false.       ! define reflection state or not")')
    write(iw,'(" use_product_state      =  .false.       ! define product (absorbing) state or not")')
    write(iw,'(" calc_steady            =  .false.       ! calculate steady-state (equilibrium) populations or not")')
    write(iw,'(" calc_Pint              =  .false.       ! calculate time integral of Pj analytically or not")')
    write(iw,'(" extrapolate            =  .false.       ! evaluate the time development of Pj(t) based on the integral equations")')
    write(iw,'(" nstate                 =  4             ! # of states")')
    write(iw,'(" ndim                   =  1             ! # of dimensions.")')
    write(iw,'(" nmol                   =  1             ! # of target molecules (typically 1)")')
    write(iw,'(" dt                     =  0.1d0         ! Time grid")')
    write(iw,'(" t_sparse               =  0.1d0         ! Sparse time-grid (used for computing TCF)")')
    write(iw,'(" t_range                =  10.0          ! Timescale for computing K-, M-, R-, and P0-functions")')
    write(iw,'(" t_extend               = 100.0          ! Extended timescale for outputting P- and Q-functions")')
    write(iw,'(" dt_tcfout              =  1.0           ! Time grid for outputting P- and Q-functions")')
    write(iw,'(" initial_state_ids      =  2 3           ! Intial state IDs")')
    write(iw,'(" reflection_state_ids   =  4             ! Reflection state IDs")')
    write(iw,'(" product_state_ids      =  1             ! Product (absorbing) state IDs")')
    write(iw,'("/")')
    write(iw,*)
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

