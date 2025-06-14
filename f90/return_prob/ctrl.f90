!=======================================================================
module mod_ctrl
!=======================================================================
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_cv
  use mod_bootstrap
  implicit none

  ! constants
  !
  integer, parameter, public :: ON  =  1
  integer, parameter, public :: OFF = -1

  integer,      parameter, public :: KinsModeCoarse = 1
  integer,      parameter, public :: KinsModeQuench = 2
  character(*), parameter, public :: KinsModes(2)   = (/'COARSE', &
                                                        'QUENCH'/)

  integer,      parameter, public :: TcfModePji = 1
  integer,      parameter, public :: TcfModeRji = 2
  character(*), parameter, public :: TcfModes(2)   = (/'PJI   ', &
                                                       'RJI   '/)

  ! structures
  !
  type :: s_option
    logical              :: check_timescale  = .false.
    logical              :: use_bootstrap    = .false.
    logical              :: use_quench       = .false.
    logical              :: use_zeropadding  = .false.
    logical              :: use_reactraj     = .false.
    logical              :: allow_state_jump = .true.    ! used only if use_quench = .true., allow the reaction 
                                                         ! state i (not reaczone) => state quench 
                                                         ! is not counted as quench event 
    logical              :: calc_pret        = .true.
    logical              :: calc_kins        = .true.
    logical              :: calc_2nd         = .false.
    logical              :: out_normfactor   = .false.
    logical              :: out_tpm          = .false.
    integer              :: kins_mode        = KinsModeQuench
    integer              :: tcf_mode         = TcfModePji

    ! for state definition
    !
    integer              :: nmol             = 1
    integer              :: ndim             = 1
    integer              :: nstate           = 3 
    integer              :: bound_id         = 1
    integer              :: reaczone_id      = 2

    ! for time scale definition
    !
    real(8)              :: dt               = - 1.0d0 
    real(8)              :: dt_out           = - 1.0d0 
    real(8)              :: t_range          = - 1.0d0
    real(8)              :: t_transient      = - 1.0d0
    real(8)              :: t_cut            = - 1.0d0

    ! for checking transient time scale
    !
    real(8)              :: t_transient_sta      = - 1.0d0
    real(8)              :: t_transient_interval = - 1.0d0
    real(8)              :: n_check_transient    = 0

    ! prepared after reading namelists
    !
    integer              :: nstep        = 0
    integer              :: nt_sparse    = 0
    integer              :: nt_range     = 0
    integer              :: nt_transient = 0
    integer              :: nt_cut       = 0

    integer              :: ntr_sta      = 0
    integer              :: ntr_interval = 0
    integer              :: nntr

    real(8), allocatable :: state_def(:, :, :)

  end type s_option

  ! subroutines
  !
  public  :: read_ctrl
  private :: read_ctrl_option

  contains

!-----------------------------------------------------------------------
    subroutine read_ctrl(input, output, option, bootopt)
!-----------------------------------------------------------------------
      implicit none

      !integer, parameter           :: iunit = 10

      type(s_input),   intent(out) :: input
      type(s_output),  intent(out) :: output
      type(s_option),  intent(out) :: option 
      type(s_bootopt), intent(out) :: bootopt 

      character(len=MaxChar)       :: f_ctrl
      integer                      :: iunit

      ! get control file name
      !
      call getarg(1, f_ctrl)

      write(iw,*)
      write(iw,'("Read_Ctrl> Reading parameters from ", a)') trim(f_ctrl)

      call open_file(trim(f_ctrl), iunit)
      call read_ctrl_input  (iunit, input)
      call read_ctrl_output (iunit, output)
      call read_ctrl_option (iunit, option)
      call read_ctrl_state (iunit, option)

      if (option%use_bootstrap) then
        call read_ctrl_bootstrap(iunit, bootopt)
      end if
      close(iunit)

    end subroutine read_ctrl
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine read_ctrl_option(iunit, option)
!-----------------------------------------------------------------------
      implicit none
!
      integer, parameter :: ndim_max = 3 
!
      integer,        intent(in)  :: iunit
      type(s_option), intent(out) :: option 

      logical                :: check_timescale = .false.
      logical                :: use_bootstrap   = .false.
      logical                :: use_quench      = .false.
      logical                :: use_zeropadding = .false.
      logical                :: use_reactraj    = .false.
      logical                :: allow_state_jump= .false.
      logical                :: calc_pret       = .true.
      logical                :: calc_kins       = .true.
      logical                :: calc_2nd        = .false.
      logical                :: out_normfactor  = .false.
      logical                :: out_tpm         = .false.
      character(len=MaxChar) :: kins_mode       = "COARSE"
      character(len=MaxChar) :: tcf_mode        = "PJI"

      ! for state definition
      !
      integer                :: nmol        = 1
      integer                :: ndim        = 1
      integer                :: nstate      = 3
      integer                :: bound_id    = 1
      integer                :: reaczone_id = 2

      ! for time scale 
      !
      real(8)                :: dt          = - 1.0d0
      real(8)                :: t_sparse    = - 1.0d0
      real(8)                :: t_range     = - 1.0d0 
      real(8)                :: t_transient = - 1.0d0 
      real(8)                :: t_cut       = - 1.0d0

      ! for checking time scale
      !
      real(8)                :: t_transient_sta      = - 1.0d0
      real(8)                :: t_transient_interval = - 1.0d0
      integer                :: n_check_transient    = 0

      integer :: i, j, nrange
      integer :: n_return_cut
      integer :: iopt, ierr

      namelist /option_param/ &
        check_timescale,      &
        use_bootstrap,        &
        use_quench,           &
        use_zeropadding,      & 
        use_reactraj,         &
        allow_state_jump,     &
        calc_pret,            &
        calc_kins,            &
        calc_2nd,             &
        out_normfactor,       &
        out_tpm,              &
        kins_mode,            &
        tcf_mode,             &
        nmol,                 &
        ndim,                 &
        nstate,               &
        bound_id,             &
        reaczone_id,          &
        dt,                   &
        t_sparse,             &
        t_range,              &
        t_transient,          &
        t_cut,                &
        t_transient_sta,      &
        t_transient_interval, &
        n_check_transient

      rewind iunit
      read(iunit, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("check_timescale      = ", a)')      get_tof(check_timescale)
      write(iw,'("use_bootstrap        = ", a)')      get_tof(use_bootstrap)
      write(iw,'("use_quench           = ", a)')      get_tof(use_quench)
      write(iw,'("use_zeropadding      = ", a)')      get_tof(use_zeropadding)
      write(iw,'("use_reactraj         = ", a)')      get_tof(use_reactraj)
      write(iw,'("allow_state_jump     = ", a)')      get_tof(allow_state_jump)
      write(iw,'("calc_pret            = ", a)')      get_tof(calc_pret)
      write(iw,'("calc_kins            = ", a)')      get_tof(calc_kins)
      write(iw,'("calc_2nd             = ", a)')      get_tof(calc_2nd)
      write(iw,'("out_normfactor       = ", a)')      get_tof(out_normfactor)
      write(iw,'("out_tpm              = ", a)')      get_tof(out_tpm)
      write(iw,'("kins_mode            = ", a)')      trim(kins_mode)
      write(iw,'("tcf_mode             = ", a)')      trim(tcf_mode)
      write(iw,*)
      write(iw,'("nmol                 = ", i0)')     nmol
      write(iw,'("ndim                 = ", i0)')     ndim
      write(iw,'("nstate               = ", i0)')     nstate
      write(iw,'("bound_id             = ", i0)')     bound_id
      write(iw,'("reaczone_id          = ", i0)')     reaczone_id
      write(iw,*)
      write(iw,'("dt                   = ", f15.7)')  dt
      write(iw,'("t_sparse             = ", f15.7)')  t_sparse
      write(iw,'("t_range              = ", f15.7)')  t_range
      write(iw,'("t_transient          = ", f15.7)')  t_transient
      write(iw,'("t_cut                = ", f15.7)')  t_cut
      write(iw,*)

      if (check_timescale) then
        write(iw,'("t_transient_sta      = ", f15.7)') t_transient_sta
        write(iw,'("t_transient_interval = ", f15.7)') t_transient_interval
        write(iw,'("n_check_transient    = ", i0)')    n_check_transient
      end if


      iopt = get_opt(kins_mode, KinsModes, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("kins_mode = ",a," is not available.")') trim(kins_mode)
        stop
      end if
      option%kins_mode = iopt

      iopt = get_opt(tcf_mode, TcfModes, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("tcf_mode = ",a," is not available.")') trim(tcf_mode)
        stop
      end if
      option%tcf_mode = iopt

      option%check_timescale      = check_timescale
      option%use_bootstrap        = use_bootstrap
      option%use_quench           = use_quench
      option%use_zeropadding      = use_zeropadding
      option%use_reactraj         = use_reactraj
      option%allow_state_jump     = allow_state_jump 
      option%calc_pret            = calc_pret
      option%calc_kins            = calc_kins
      option%calc_2nd             = calc_2nd
      option%out_normfactor       = out_normfactor
      option%out_tpm              = out_tpm

      option%nmol                 = nmol
      option%ndim                 = ndim
      option%nstate               = nstate
      option%bound_id             = bound_id
      option%reaczone_id          = reaczone_id

      option%dt                   = dt
      option%t_range              = t_range
      option%t_transient          = t_transient
      option%t_cut                = t_cut

      option%t_transient_sta      = t_transient_sta
      option%t_transient_interval = t_transient_interval
      option%n_check_transient    = n_check_transient


      ! Check input parameter
      !
      call check_input_parameter_real8(dt,          0.0d0, "up", "dt",          .true.)
      call check_input_parameter_real8(t_range,     0.0d0, "up", "t_range",     .true.)
      if (.not. check_timescale) then 
        if (calc_kins) &
          call check_input_parameter_real8(t_transient, 0.0d0, "up", "t_transient", .true.)
      !call check_input_parameter_real8(t_cut,       0.0d0, "up", "t_cut",       .true.)
      end if

      if (option%check_timescale) then
        call check_input_parameter_real8(t_transient_sta,      0.0d0, "up", "t_transient_sta",      .true.)
        call check_input_parameter_real8(t_transient_interval, 0.0d0, "up", "t_transient_interval", .true.)
        call check_input_parameter_integer(n_check_transient,  1,     "up", "n_check_transient",    .true.)
      end if

      ! Convert time constants from real to integer
      !
      if (option%t_cut > 0.0d0) then
        option%nstep = nint(option%t_cut / option%dt) 
      else
        option%nstep = 0
      end if

      option%nt_range     = nint(option%t_range / option%dt)
      option%nt_transient = nint(option%t_transient / option%dt)

      if (option%check_timescale) then
        option%ntr_sta      = nint(option%t_transient_sta / option%dt)
        option%ntr_interval = nint(option%t_transient_interval / option%dt)
        option%nntr         = option%n_check_transient
      end if


      ! Introduce sparse time grid (output time grid) 
      ! (If t_sparse < 0.0, no changes) 
      if (t_sparse < 0.0d0) then
        t_sparse = dt
      end if

      option%nt_sparse = nint(t_sparse / dt)
      option%dt_out    = dt * option%nt_sparse
      if (option%nt_sparse > 1) then
        option%nt_range = nint(dble(option%nt_range) / dble(option%nt_sparse))
      end if

      ! Combination check
      !
      if (use_reactraj) then
        if (.not. use_quench) then
          write(iw,'("Read_Ctrl_Option> Error.")')
          write(iw,'("use_reactraj = .true. & use_quench = .false. &
                     &should give unreasonable results.")')
          stop
        end if
      end if

      if (.not. use_quench) then
        if (use_zeropadding) then
          write(iw,'("Read_Ctrl_Option> Error.")')
          write(iw,'("use_zeropading = .true. is not available if use_quench = .false.")')
          stop
        end if

        if (option%kins_mode == KinsModeQuench) then
          write(iw,'("Read_Ctrl_Option> Error.")')
          write(iw,'("kins_mode = QUENCH is not available if use_quench = .false.")')
          stop
        end if
      end if

      if (use_quench .and. option%kins_mode == KinsModeCoarse) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("use_quench = .true. and kins_mode = COARSE can not be used at the same time")')
        stop
      end if

      if (.not. calc_kins .and. check_timescale) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("check_timescale = .true. is not available if calc_kins = .false.")')
        stop
      end if

      if (.not. calc_pret .and. calc_2nd) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("calc_2nd = .true. is not available if calc_pret = .false.")')
        stop
      end if

      ! Memory allocation
      !
      allocate(option%state_def(2, option%ndim, option%nstate + 1))

    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine read_ctrl_state(iunit, option)
!-----------------------------------------------------------------------
      implicit none

      integer,        intent(in)    :: iunit
      type(s_option), intent(inout) :: option

      integer                :: istate, idm, imm, ierr, add
      character(len=MaxChar) :: line


      add = 0
      if (option%use_quench) then
        add = 1 
      end if

      call seek_line(iunit, '&state', ierr)

      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_State> Error.")')
        write(iw,'("&state section is not found.")')
        stop
      end if

      istate           = 0
      option%state_def = 0.0d0
      do while (.true.)
        read(iunit,'(a)', end = 100) line

        if (line(1:1) /= "#") then
          istate = istate + 1

          read(line,*) ((option%state_def(imm, idm, istate), imm = 1, 2), &
                         idm = 1, option%ndim) 
        end if

        if (istate == option%nstate + add) exit

      end do

100   continue

      write(iw,*)
      write(iw,'("Read_Ctrl_State>  States are defined as")')
      do istate = 1, option%nstate + add
        write(iw,'(i5)', advance = 'no') istate
        if (istate == option%bound_id) then
          write(iw,'("  BOUNDED STATE")')
        else if (istate == option%reaczone_id) then
          write(iw,'("  REACTION ZONE")')
        else if (istate == option%nstate + 1) then
          write(iw,'("  QUENCH ZONE")')
        else
          write(iw,*)
        end if

        do idm = 1, option%ndim
          write(iw,'("  component ", i5, " : ", f15.7, " <===> ", f15.7)') &
            idm,                              &
            option%state_def(1, idm, istate), &
            option%state_def(2, idm, istate)
        end do

        write(iw,*)

      end do

    end subroutine read_ctrl_state
!-----------------------------------------------------------------------

end module
!=======================================================================
