!=======================================================================
program main 
!=======================================================================
  use mod_const
  use mod_ctrl
  use mod_analyze

  type(s_input)   :: input
  type(s_output)  :: output
  type(s_option)  :: option 
  type(s_cvinfo)  :: cvinfo
  type(s_bootopt) :: bootopt 

  call show_title
  call show_usage
  call read_ctrl(input, output, option, cvinfo, bootopt)
  call analyze(input, output, option, cvinfo, bootopt)
  call termination("rr_analysis")

end program main 
!=======================================================================

!-----------------------------------------------------------------------
subroutine show_title
!-----------------------------------------------------------------------
  implicit none

  write(6,'("==================================================")')
  write(6,*)
  write(6,'("          Restricted RDF Analysis")')
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
    write(iw,'(" flist_cv      = ""cvlist""    ! File that contains List of CV")')
    write(iw,'(" flist_weight  = ""cvlist""    ! File that contains List of Weight (if necessary)")')
    write(iw,'("/")')
    write(iw,'("&output_param")')
    write(iw,'(" fhead         = ""out"" ! Header of Output files")')
    write(iw,'("/")')
    write(iw,'("&option_param")')
    write(iw,'(" calcfe        = .true.        ! calculate standard free energy or not (default: .true.")')
    write(iw,'(" ndim          =  1            ! CV dimensions")')
    write(iw,'(" temperature   =  298.0        ! temperature (K)")')
    write(iw,'(" dr            =  0.2d0        ! distance grid spacing (A)")')
    write(iw,'(" ngrid         = 1000          ! number of distance grids")')
    write(iw,'(" nsta          =  1            ! first frame for analysis")')
    write(iw,'(" use_bootstrap = .false.       ! use bootstrap analysis (default:.false.)")')
    write(iw,'(" vol0          = 1661.0d0      ! standard volume (A^3) (default: 1661.0)")')
    write(iw,'(" urange        = 25.0   30.0   ! distance range of unbound state (1st: minimum, 2nd: maximum)")')
    write(iw,'(" react_range   =  0.0   10.0   ! range of reactive state (1st: minimum, 2nd: maximum)")')
    write(iw,'("                  5.0   20.0   ! range of reactive state (1st: minimum, 2nd: maximum)")')
    write(iw,'(" ")')
    write(iw,'("/")')
    write(iw,'("&bootopt_param")')
    write(iw,'(" iseed         =  3141592      ! input seed_number (if iseed <= 0, seed is generated)")')
    write(iw,'(" nsample       =  1000         ! Number of samples to be selected for each trial")')
    write(iw,'(" duplicate     = .true.        ! Whether to allow duplicated selections for each trial")')
    write(iw,'(" ntrial        =  100          ! Number of bootstrap trials")')
    write(iw,'(" ")')
    write(iw,'("/")')

    stop
  end if


end subroutine show_usage
!-----------------------------------------------------------------------

