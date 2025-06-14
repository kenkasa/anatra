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
  integer, parameter, public :: Nstate   = 3 
  integer, parameter, public :: REACTIVE = 1 
  integer, parameter, public :: BOUND    = 2 
  integer, parameter, public :: UNBOUND  = 3 
  integer, parameter, public :: OTHERS   = 4
  integer, parameter, public :: StateInfo(Nstate) = (/1, 2, 3/)

  integer, parameter, public :: ON  =  1
  integer, parameter, public :: OFF = -1
  integer, parameter, public :: StateConv(Nstate+1) = (/OFF, ON, OFF, OFF/)

  integer,      parameter, public :: TubTypeWHOLE    = 1
  integer,      parameter, public :: TubTypeCUTFINAL = 2
  character(*), parameter, public :: TubTypes(2)     = (/'WHOLE     ',&
                                                         'CUTFINAL  '/)

  ! structures
  !
  type :: s_option
    logical              :: calc_rp       = .true.
    logical              :: check_tscale  = .false.
    logical              :: use_moving    = .true.
    logical              :: use_window    = .false.
    logical              :: use_bootstrap = .false.
    integer              :: nstep         = 0
    integer              :: ntr_sta       = 10 
    integer              :: ntr_interval  = 10
    integer              :: nntr          = 200
    integer              :: ndim          = 1
    integer              :: nt_range      = 1000
    integer              :: nt_transient  = 100
    integer              :: nsta          = 1
    integer              :: tubtype       = TubTypeWHOLE
    real(8)              :: dt            = 1.0d0
    real(8), allocatable :: state_def(:, :, :)
    real(8), allocatable :: react_range(:, :)
    real(8), allocatable :: bound_range(:, :)
    real(8), allocatable :: unbound_range(:, :)

    real(8)              :: t_transient   = 0.0d0  ! Hidden option.
    real(8)              :: t_rpcut       = 0.0d0  ! Hidden option. trajectory length (time unit) 
                                                   ! used for P(t) calc.
    real(8)              :: t_range       = 0.0d0  ! Hidden option. calculation time range of P(t)
                                                   ! (time unit). If non-zero value is specified,
                                                   ! the setting specified with nt_range is replaced 
                                                   ! by # of points determined from t_range
    real(8)              :: t_return_cut  = 0.0d0  ! Hidden option. If non-zero value is specified,
                                                   ! returned configurations after t_return_cut dissociation
                                                   ! is neglected
    integer              :: n_return_cut  = 0
    !logical              :: duplicate     = .true.
    !logical              :: seed_input    = .false.
    !character(len=MaxChar) :: fseed       = ""
    !integer              :: gen_num
    !integer              :: samples
  end type s_option

  ! subroutines
  !
  public  :: read_ctrl
  private :: read_ctrl_option
  !private :: read_ctrl_btstrpoption

  contains

!-----------------------------------------------------------------------
    subroutine read_ctrl(input, output, option, cvinfo, bootopt)
!-----------------------------------------------------------------------
      implicit none

      integer, parameter           :: iunit = 10

      type(s_input),   intent(out) :: input
      type(s_output),  intent(out) :: output
      type(s_option),  intent(out) :: option 
      type(s_cvinfo),  intent(out) :: cvinfo
      type(s_bootopt), intent(out) :: bootopt 

      character(len=MaxChar)       :: f_ctrl

      ! get control file name
      !
      call getarg(1, f_ctrl)

      write(iw,*)
      write(iw,'("Read_Ctrl> Reading parameters from ", a)') trim(f_ctrl)
      open(iunit, file=trim(f_ctrl), status='old')
        call read_ctrl_input  (iunit, input)
        call read_ctrl_output (iunit, output)
        call read_ctrl_option (iunit, option)

        if (option%use_bootstrap) then
          call read_ctrl_bootstrap(iunit, bootopt)
        end if
      close(iunit)

      open(iunit, file=trim(input%flist))
        call read_ctrl_cvinfo  (iunit, cvinfo)
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

      logical                :: calc_rp                    = .true.
      logical                :: check_tscale               = .false.
      logical                :: use_moving                 = .false.
      logical                :: use_window                 = .false.
      logical                :: use_bootstrap              = .false.
      integer                :: nstep                      = 0 
      integer                :: ntr_sta                    = 10
      integer                :: ntr_interval               = 10
      integer                :: nntr                       = 200
      integer                :: ndim                       = 1
      integer                :: nt_range                   = 100
      integer                :: nt_transient               = 100
      integer                :: nsta                       = 1
      character(len=MaxChar) :: tubtype                    = "WHOLE"
      real(8)                :: dt                         = 1.0d0
      real(8)                :: t_transient                = 0.0d0  ! Hidden option
      real(8)                :: t_rpcut                    = 0.0d0  ! Hidden option
      real(8)                :: t_range                    = 0.0d0  ! Hidden option
      real(8)                :: t_return_cut               = 0.0d0  ! Hidden option
      real(8)                :: react_range(2, ndim_max)   = 0.0d0
      real(8)                :: bound_range(2, ndim_max)   = 0.0d0
      real(8)                :: unbound_range(2, ndim_max) = 0.0d0

      integer :: i, j, nrange
      integer :: n_return_cut
      integer :: iopt, ierr

      namelist /option_param/ calc_rp, check_tscale, use_moving,                 &
                              use_window, use_bootstrap, nstep,                  &
                              ntr_sta, ntr_interval,                             &
                              nntr, ndim, nt_range, nt_transient,                &
                              nsta, dt, t_transient, t_rpcut, t_range,           &
                              t_return_cut, tubtype, react_range,                &
                              bound_range, unbound_range

      rewind iunit
      read(iunit, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("calc_rp          = ", a)')  get_tof(calc_rp)
      write(iw,'("use_moving       = ", a)')  get_tof(use_moving)
      write(iw,'("use_window       = ", a)')  get_tof(use_window)
      write(iw,'("use_bootstrap    = ", a)')  get_tof(use_bootstrap)
      write(iw,'("ndim             = ", i0)') ndim
      write(iw,'("step             = ", i0)') nstep 
      write(iw,'("nt_range         = ", i0)') nt_range
      write(iw,'("nt_transient     = ", i0)') nt_transient
      write(iw,'("nsta             = ", i0)') nsta
      write(iw,'("tubtype          = ", a)')  trim(tubtype)
      write(iw,*)
      write(iw,'("react_range      =")')
      do i = 1, ndim
        write(iw,'(es15.7, " <= component ",i0," < ",es15.7)') &
          react_range(1, i), i, react_range(2, i) 
      end do
      write(iw,*)
      write(iw,'("bound_range      =")')
      do i = 1, ndim
        write(iw,'(es15.7, " <= component ",i0," < ",es15.7)') &
          bound_range(1, i), i, bound_range(2, i) 
      end do
      write(iw,*)
      write(iw,'("unbound_range    =")')
      do i = 1, ndim
        write(iw,'(es15.7, " <= component ",i0," < ",es15.7)') &
          unbound_range(1, i), i, unbound_range(2, i) 
      end do

      iopt = get_opt(tubtype, TubTypes, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("tubtype = ",a," is not available.")') trim(tubtype)
        stop
      end if
      option%tubtype       = iopt
                          
      option%calc_rp       = calc_rp
      option%use_moving    = use_moving
      option%use_window    = use_window
      option%use_bootstrap = use_bootstrap 
      option%check_tscale  = check_tscale
      option%nstep         = nstep
      option%ntr_sta       = ntr_sta
      option%ntr_interval  = ntr_interval
      option%nntr          = nntr 
      option%ndim          = ndim
      option%nt_range      = nt_range
      option%nt_transient  = nt_transient
      option%nsta          = nsta
      option%dt            = dt
      option%t_transient   = t_transient
      option%t_rpcut       = t_rpcut
      option%t_range       = t_range
      option%t_return_cut  = t_return_cut

      nrange = nint(t_range / option%dt)
      if (nrange /= 0) then
        write(iw,*)
        write(iw,'("Read_Ctrl_Option> Remark: Detect non-zero t_range value.")')
        write(iw,'("nt_range is changed from ",i0," to ",i0)') nt_range, nrange
        option%nt_range = nrange
      end if

      n_return_cut        = nint(t_return_cut / option%dt)
      option%n_return_cut = 0
      if (n_return_cut /= 0) then
        write(iw,*)
        write(iw,'("Read_Ctrl_Option> Remark: Detect non-zero t_return_cut value.")')
        write(iw,'("")')
        option%n_return_cut = n_return_cut
      end if

      nrange = nint(t_transient / option%dt)
      if (nrange /= 0) then
        write(iw,*)
        write(iw,'("Read_Ctrl_Option> Remark: Detect non-zero t_transient value.")')
        write(iw,'("nt_transient is changed from ",i0," to ",i0)') &
          nt_transient, nrange

        option%nt_transient = nrange
      end if

      allocate(option%react_range(1:2, ndim))
      allocate(option%bound_range(1:2, ndim))
      allocate(option%unbound_range(1:2, ndim))
      allocate(option%state_def(2, ndim, Nstate))

      do i = 1, ndim
        option%react_range(1:2, i)   = react_range(1:2, i)
        option%bound_range(1:2, i)   = bound_range(1:2, i)
        option%unbound_range(1:2, i) = unbound_range(1:2, i)
      end do

      do i = 1, ndim
        option%state_def(1:2, i, 1) = react_range(1:2, i)
      end do

      do i = 1, ndim
        option%state_def(1:2, i, 2) = bound_range(1:2, i)
      end do

      do i = 1, ndim
        option%state_def(1:2, i, 3) = unbound_range(1:2, i)
      end do

      !do i = 1, 3
      !  do j = 1, ndim
      !    write(iw,'(2f15.7)') option%state_def(1, j, i), &
      !                         option%state_def(2, j, i)
      !  end do
      !end do

      ! Combination check
      !
      if (option%tubtype == TubTypeCUTFINAL &
          .and. .not. option%use_moving) then

        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("tubtype = CUTFINAL is not available &
                   &when use_moving = .true.")')
        stop
      end if 

    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

end module
!=======================================================================
