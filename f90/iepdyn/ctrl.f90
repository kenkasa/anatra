!=======================================================================
module mod_ctrl
!=======================================================================
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_cv
  implicit none

  ! constants
  !
  integer,      parameter, public :: MaxStates    =  100 
  integer,      parameter, public :: NotSpecified = -100

  integer,      parameter, public :: InputTypeTimeSeries   = 1
  integer,      parameter, public :: InputTypeHistogram    = 2
  character(*), parameter, public :: InputTypes(2) = (/'TIMESERIES', &
                                                       'HISTOGRAM '/)

  integer,      parameter, public :: ErrexTypeExclude = 1
  integer,      parameter, public :: ErrexTypeInclude = 2
  character(*), parameter, public :: ErrexTypes(2)    = (/'EXCLUDE', &
                                                          'INCLUDE'/)

  integer,      parameter, public :: CumDirecIncrease = 1
  integer,      parameter, public :: CumDirecDecrease = 2 
  character(*), parameter, public :: CumDirecTypes(2) = (/'INCREASE', &
                                                          'DECREASE'/)

  ! structures
  !
  type :: s_option
    logical :: use_perturbed_traj   = .false.
    logical :: use_reflection_state = .false.
    logical :: use_product_state    = .false.
    logical :: use_dissociate_state = .false.
    logical :: use_constant_Qij     = .false.
    logical :: use_rsto_Rij         = .true.
    logical :: output_histogram     = .false.
    logical :: extrapolate          = .false.
    logical :: calc_Pint            = .false.
    logical :: calc_Steady          = .false.
    logical :: check_Kijk           = .false.
    logical :: check_senserr        = .false.
    logical :: check_blockave       = .false.
    logical :: check_cumulative     = .false.

    integer :: input_type           = InputTypeTimeSeries
    integer :: errex_type           = ErrexTypeExclude 
    integer :: cumdirec             = CumDirecIncrease

    integer :: nmol                               = NotSpecified 
    integer :: ndim                               = NotSpecified 
    integer :: nstate                             = NotSpecified
    integer :: reflection_state_ids   (MaxStates) = NotSpecified 
    integer :: product_state_ids      (MaxStates) = NotSpecified
    integer :: dissociate_state_ids   (MaxStates) = NotSpecified
    integer :: initial_state_ids      (MaxStates) = NotSpecified
    integer :: err_exception_state_ids(MaxStates) = NotSpecified
    integer :: nkmax                              = 1000
    integer :: nblock                             = 5
    integer :: ncum                               = 10

    ! File names 
    !
    character(len=MaxChar) :: f_unperturbed_id = ''
    character(len=MaxChar) :: f_cQij           = ''

    ! for free-energy calculation 
    !
    real(8)              :: temperature  = 298.0d0

    ! for time scale definition
    !
    real(8)              :: dt           = - 1.0d0 
    real(8)              :: t_sparse     = - 1.0d0
    real(8)              :: t_range      = - 1.0d0
    real(8)              :: t_extend     = - 1.0d0
    real(8)              :: dt_tcfout    = - 1.0d0

    ! prepared after reading namelists
    !
    integer              :: nreflect     = 0
    integer              :: nproduct     = 0
    integer              :: ndissoc      = 0
    integer              :: ninitial     = 0
    integer              :: nselect      = 0

    real(8)              :: dt_out
    integer              :: nstep        = 0
    integer              :: nt_sparse    = 0
    integer              :: nt_range     = 0
    integer              :: nt_extend    = 0
    integer              :: nt_tcfout    = 0 

    logical, allocatable :: is_initial  (:)
    logical, allocatable :: is_product  (:)
    logical, allocatable :: is_reflect  (:)
    logical, allocatable :: is_dissoc   (:)
    logical, allocatable :: is_errex    (:) 
    real(8), allocatable :: state_def   (:, :, :)
    real(8), allocatable :: state_weight(:)
    real(8), allocatable :: state_weight_unnorm(:)

  end type s_option

  type s_extra_input
    character(len=MaxChar) :: f_connect   = ''
    character(len=MaxChar) :: f_init_d    = ''
    character(len=MaxChar) :: f_rbin_head = ''
    character(len=MaxChar) :: f_kbin_head = ''
  end type s_extra_input

  type :: s_timegrid
    integer :: ng
    real(8), allocatable :: val(:)
    integer, allocatable :: ind(:)
  end type s_timegrid

  ! subroutines
  !
  public  :: read_ctrl
  private :: read_ctrl_option
  private :: read_ctrl_state
  private :: show_input
  private :: show_output

  contains

!-----------------------------------------------------------------------
    subroutine read_ctrl(input, einput, output, option, timegrid)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),       intent(out) :: input
      type(s_extra_input), intent(out) :: einput
      type(s_output),      intent(out) :: output
      type(s_option),      intent(out) :: option 
      type(s_timegrid),    intent(out) :: timegrid 

      ! I/O
      !
      integer                :: io
      character(len=MaxChar) :: f_ctrl


      ! get control file name
      !
      call getarg(1, f_ctrl)

      write(iw,*)
      write(iw,'("Read_Ctrl> Reading parameters from ", a)') trim(f_ctrl)
      call open_file        (f_ctrl, io)
      call read_ctrl_input  (io, input)
      call show_input       (input)
      call read_ctrl_output (io, output)
      call show_output      (output)
      call read_ctrl_option (io, option, timegrid)
      call read_ctrl_state  (io, option)
      close(io)

    end subroutine read_ctrl
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine read_ctrl_option(io, option, timegrid)
!-----------------------------------------------------------------------
      implicit none
!
      integer, parameter :: ndim_max = 3 
!
      integer,          intent(in)  :: io
      type(s_option),   intent(out) :: option
      type(s_timegrid), intent(out) :: timegrid

      logical :: use_perturbed_traj   = .false.
      logical :: use_reflection_state = .false.
      logical :: use_product_state    = .false.
      logical :: use_dissociate_state = .false.
      logical :: use_constant_Qij     = .false.
      logical :: use_rsto_Rij         = .true.
      logical :: output_histogram     = .false.
      logical :: extrapolate          = .false.
      logical :: calc_Pint            = .false.
      logical :: calc_Steady          = .false.
      logical :: check_Kijk           = .false.
      logical :: check_senserr        = .false.
      logical :: check_blockave       = .false.
      logical :: check_cumulative     = .false.

      character(len=MaxChar) :: input_type       = 'TIMESERIES'
      character(len=MaxChar) :: errex_type       = 'EXCLUDE'
      character(len=MaxChar) :: cumdirec         = 'INCREASE'
      character(len=MaxChar) :: f_unperturbed_id = ''
      character(len=MaxChar) :: f_cQij           = '' 
      
      integer :: nmol                               = NotSpecified
      integer :: ndim                               = NotSpecified
      integer :: nstate                             = NotSpecified
      integer :: reflection_state_ids   (MaxStates) = NotSpecified
      integer :: product_state_ids      (MaxStates) = NotSpecified
      integer :: dissociate_state_ids   (MaxStates) = NotSpecified
      integer :: initial_state_ids      (MaxStates) = NotSpecified
      integer :: err_exception_state_ids(MaxStates) = NotSpecified
      integer :: nkmax                              = 1000
      integer :: nblock                             = 5
      integer :: ncum                               = 10
      real(8) :: temperature                        = 300.0d0
      real(8) :: dt
      real(8) :: t_sparse
      real(8) :: t_range
      real(8) :: t_extend
      real(8) :: dt_tcfout

      ! Local
      !
      integer :: nreflect = 0, nproduct = 0, ndissoc = 0, ninitial = 0, nselect = 0
      integer :: nerrex   = 0

      ! Dummy
      !
      integer :: i, j, it, id
      integer :: iopt, ierr


      namelist /option_param/   &
        use_perturbed_traj,     &
        use_reflection_state,   &
        use_product_state,      &
        use_dissociate_state,   &
        use_constant_Qij,       &
        use_rsto_Rij,           &
        output_histogram,       &
        extrapolate,            &
        check_Kijk,             &
        check_senserr,          &
        check_blockave,         &
        check_cumulative,       &
        calc_Pint,              &
        calc_Steady,            &
        input_type,             &
        errex_type,             &
        cumdirec,               &
        f_unperturbed_id,       &
        f_cQij,                 &
        nmol,                   &
        ndim,                   &
        nstate,                 &
        reflection_state_ids,   &
        product_state_ids,      &
        dissociate_state_ids,   &
        initial_state_ids,      &
        err_exception_state_ids,&
        nkmax,                  &
        nblock,                 &
        ncum,                   &
        temperature,            &
        dt,                     &
        t_sparse,               &
        t_range,                &
        t_extend,               &
        dt_tcfout

      rewind io 
      read(io, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("input_type           = ", a)')     trim(input_type)

      write(iw,'("use_perturbed_traj   = ", a)')   get_tof(use_perturbed_traj)
      write(iw,'("use_reflection_state = ", a)')   get_tof(use_reflection_state)
      write(iw,'("use_product_state    = ", a)')   get_tof(use_product_state)
      write(iw,'("use_dissociate_state = ", a)')   get_tof(use_dissociate_state)
      write(iw,'("use_constant_Qij     = ", a)')   get_tof(use_constant_Qij)
      write(iw,'("use_rsto_Rij         = ", a)')   get_tof(use_rsto_Rij)
      write(iw,'("output_histogram     = ", a)')   get_tof(output_histogram)
      write(iw,'("check_Kijk           = ", a)')   get_tof(check_Kijk)
      write(iw,'("check_senserr        = ", a)')   get_tof(check_senserr)
      write(iw,'("check_blockave       = ", a)')   get_tof(check_blockave)
      write(iw,'("calc_Pint            = ", a)')   get_tof(calc_Pint)
      write(iw,'("calc_Steady          = ", a)')   get_tof(calc_Steady)
      write(iw,'("f_unperturbed_id     = ", a)')   trim(f_unperturbed_id)
      write(iw,'("f_cQij               = ", a)')   trim(f_cQij)

      write(iw,'("extrapolate          = ", a)')      get_tof(extrapolate)
      write(iw,'("nmol                 = ", i0)')     nmol
      write(iw,'("ndim                 = ", i0)')     ndim
      write(iw,'("nstate               = ", i0)')     nstate
      write(iw,'("temperature          = ", f20.10)') temperature
      write(iw,'("nblock               = ", i0)')     nblock 

      write(iw,'("dt                   = ", f20.10)') dt
      write(iw,'("t_sparse             = ", f20.10)') t_sparse
      write(iw,'("t_range              = ", f20.10)') t_range

      if (extrapolate) then
        write(iw,'("t_extend             = ", f20.10)') t_extend
      end if

      if (nstate == NotSpecified) then
        write(iw,*)
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("nstate should be specified.")')
        stop
      end if


      ! Reflection
      !
      allocate(option%is_reflect(nstate))
      option%is_reflect = .false.

      if (use_reflection_state) then
        nreflect = 0
        write(iw,*)
        do i = 1, MaxStates
          if (reflection_state_ids(i) /= NotSpecified)  then
            write(iw,'("reflection_state_ids  ", i0, " : ", i0)') i, reflection_state_ids(i)
            option%is_reflect(reflection_state_ids(i)) = .true.
            nreflect                                   = nreflect + 1
          else
            exit
          end if
        end do

      end if

      ! Product
      !
      allocate(option%is_product(nstate))
      option%is_product = .false.

      if (use_product_state) then
        nproduct = 0
        write(iw,*)
        do i = 1, MaxStates
          if (product_state_ids(i) /= NotSpecified)  then
            write(iw,'("product_state_ids  ", i0, " : ", i0)') i, product_state_ids(i)
            option%is_product(product_state_ids(i)) = .true.
            nproduct                                = nproduct + 1
          else
            exit
          end if
        end do
      end if

      ! Dissociate
      !
      allocate(option%is_dissoc(nstate))
      option%is_dissoc = .false.

      if (use_dissociate_state) then
        ndissoc = 0
        write(iw,*)
        do i = 1, MaxStates
          if (dissociate_state_ids(i) /= NotSpecified)  then
            write(iw,'("dissociate_state_ids  ", i0, " : ", i0)') i,  dissociate_state_ids(i)
            option%is_dissoc(dissociate_state_ids(i)) = .true.
            ndissoc                                   = ndissoc + 1
          else
            exit
          end if
        end do
      end if

      ! Initial
      !
      ninitial = 0
      write(iw,*)
      do i = 1, MaxStates
        if (initial_state_ids(i) /= NotSpecified) then
          write(iw,'("initial_state_ids  ", i0, " : ", i0)') i, initial_state_ids(i)
          ninitial = ninitial + 1
        else
          exit
        end if
      end do

      allocate(option%is_initial(nstate))

      option%is_initial = .false.
      do i = 1, ninitial 
        option%is_initial(initial_state_ids(i)) = .true. 
      end do


      iopt = get_opt(input_type, InputTypes, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("input_type = ",a," is not available.")') trim(input_type)
      end if
      option%input_type = iopt

      iopt = get_opt(errex_type, ErrexTypes, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("errex_type = ",a," is not available.")') trim(errex_type)
      end if
      option%errex_type = iopt

      iopt = get_opt(cumdirec, CumDirecTypes, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("cumdirec = ",a," is not available.")') trim(cumdirec)
      end if
      option%cumdirec = iopt

      option%use_perturbed_traj   = use_perturbed_traj
      option%use_reflection_state = use_reflection_state
      option%use_product_state    = use_product_state
      option%use_dissociate_state = use_dissociate_state
      option%use_constant_Qij     = use_constant_Qij
      option%use_rsto_Rij         = use_rsto_Rij
      option%output_histogram     = output_histogram
      option%extrapolate          = extrapolate
      option%check_Kijk           = check_Kijk
      option%check_senserr        = check_senserr
      option%check_blockave       = check_blockave
      option%check_cumulative     = check_cumulative
      option%calc_Pint            = calc_Pint
      option%calc_Steady          = calc_Steady

      option%f_unperturbed_id     = f_unperturbed_id
      option%f_cQij               = f_cQij

      option%nmol                 = nmol
      option%ndim                 = ndim
      option%nstate               = nstate

      option%nreflect             = nreflect
      option%reflection_state_ids = reflection_state_ids

      option%nproduct             = nproduct
      option%product_state_ids    = product_state_ids

      option%ndissoc              = ndissoc
      option%dissociate_state_ids = dissociate_state_ids

      option%ninitial             = ninitial
      option%initial_state_ids    = initial_state_ids

      option%nkmax                = nkmax
      option%nblock               = nblock
      option%ncum                 = ncum
      option%temperature          = temperature
      option%dt                   = dt
      option%t_sparse             = t_sparse
      option%t_range              = t_range
      option%t_extend             = t_extend
      option%dt_tcfout            = dt_tcfout

      ! Error exception states (hidden options)
      ! is_errex(is) = .true.  ---> exclude state is for error estimation
      ! is_errex(is) = .false. ---> include

      ! if errex_type = EXCLUDE, specify states to be excluded
      !                 INCLUDE, specify states to be included
      allocate(option%is_errex(nstate))
      if (option%errex_type == ErrexTypeExclude) option%is_errex = .false.
      if (option%errex_type == ErrexTypeInclude) option%is_errex = .true.

      nerrex = 0
      do i = 1, MaxStates
        if (err_exception_state_ids(i) /= NotSpecified) then
          id = err_exception_state_ids(i)
          write(iw,'("err_exception_state_ids  ", i0, " : ", i0)') i, id 
          if (option%errex_type == ErrexTypeExclude) option%is_errex(id) = .true.
          if (option%errex_type == ErrexTypeInclude) option%is_errex(id) = .false.  
        end if
      end do

      ! Setup time grids
      !
      if (option%t_sparse <= 0.0d0) then
        option%t_sparse = option%dt
      end if

      if (option%dt_tcfout <= 0.0d0) then
        option%dt_tcfout = option%dt
      end if

      option%nt_sparse = nint(option%t_sparse / option%dt)
      option%dt_out    = option%dt * option%nt_sparse

      !option%nt_range  = nint(option%t_range   / option%dt)
      !option%nt_extend = nint(option%t_extend  / option%dt)
      !option%nt_tcfout = nint(option%dt_tcfout / option%dt) 
      option%nt_range  = nint(option%t_range   / option%dt_out)
      option%nt_extend = nint(option%t_extend  / option%dt_out)
      option%nt_tcfout = nint(option%dt_tcfout / option%dt_out) 

      write(iw,*)
      write(iw,'("dt_tcfout            = ", f20.10)') option%dt_tcfout 
      write(iw,'("dt_out               = ", f20.10)') option%dt_out 
      write(iw,'("nt_sparse            = ", i0)')     option%nt_sparse
      write(iw,'("nt_tcfout            = ", i0)')     option%nt_tcfout

      allocate(timegrid%val(0:option%nt_range)) 
      allocate(timegrid%ind(0:option%nt_range))

      timegrid%ng = option%nt_range

      do it = 0, option%nt_range
        timegrid%ind(it) = option%nt_sparse * it
        timegrid%val(it) = option%dt_out    * it
      end do 
      
      ! Check
      !
      if (nmol == NotSpecified) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("nmol shouldb be specified.")')
        stop
      end if

      if (ndim == NotSpecified) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("ndim shouldb be specified.")')
        stop
      end if

      if (nstate == NotSpecified) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("nstate shouldb be specified.")')
        stop
      end if

      if (use_perturbed_traj .and. nmol > 1) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("nmol should be 1 &
                   &if use_perturbed_traj = .true.")')
        stop
      end if

      if (use_reflection_state .and. nreflect < 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("reflect_state_ids should be specified &
                   &if use_reflection_state = .true.")')
        stop
      end if

      if (use_product_state .and. nproduct < 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("product_state_ids should be specified &
                   &if use_product_state = .true.")')
        stop
      end if

      if (output_histogram) then
        if (check_senserr .or. check_blockave .or. check_cumulative) then
          write(iw,'("Read_Ctrl_Option> Error.")')
          write(iw,'("output_histogram should be turned off when&
                    & check_senserr, check_blockave, or check_cumulative is true")')
          stop
        end if
      end if

      i = 0
      if (check_senserr)    i = i + 1
      if (check_blockave)   i = i + 1
      if (check_cumulative) i = i + 1

      if (i == 1) then
        if (.not. use_perturbed_traj) then
          write(iw,'("Read_Ctrl_Option> Error.")')
          write(iw,'("use_perturbed_traj should be turned on when &
                     &check_senserr, check_blockave, and check_cumulative is true.")')
          stop
        end if
      else if (i > 1) then
         write(iw,'("Read_Ctrl_Option> Error.")')
         write(iw,'("check_blockave, check_senserr, and check_cumulative &
                    &can not be performed at the same time")')
         stop
      end if

      if (calc_Pint .and. calc_Steady) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("Calc_Pint and Calc_Steady can not be used at the same time.")')
        stop
      end if

      if (use_constant_Qij) then
        if (calc_Pint .or. calc_Steady) then
          write(iw,'("Read_Ctrl_Option> Error.")')
          write(iw,'("use_constant_Qij can not be used &
                     &when calc_Pint = .true. or calc_Steady = .true.")')
          stop
        end if
      end if

      if (use_rsto_Rij) then
        if (option%input_type == InputTypeHistogram) then
          !write(iw,'("Read_Ctrl_Option> Warning.")')
          !write(iw,'("use_rsto_Rij = .true. is effective only if input_type = TIMESERIES.")')
          !write(iw,'(">> changed to use_rsto_Rij = .false.")')
          option%use_rsto_Rij = .false.
        end if

        !if (.not. option%use_perturbed_traj) then
        !  write(iw,'("Read_Ctrl_Option> Error.")')
        !  write(iw,'("use_rsto_Rij = .true. should  be used with use_perturbed_traj = .true.")')
        !  stop
        !end if 
      end if

      ! Memory allocation
      !
      allocate(option%state_def(2, option%ndim, option%nstate))
      allocate(option%state_weight(option%nstate))
      allocate(option%state_weight_unnorm(option%nstate))


    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine read_ctrl_state(io, option)
!-----------------------------------------------------------------------
      implicit none

      integer,        intent(in)    :: io
      type(s_option), intent(inout) :: option

      ! Local
      !
      character(len=MaxChar) :: line
      integer                :: nstate
      real(8)                :: weight_sum

      ! Dummy
      !
      integer :: istate, jstate, idm, imm, ierr
      integer :: init


      ! Setup
      !
      nstate   = option%nstate

      call seek_line(io, '&state', ierr)

      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_State> Error.")')
        write(iw,'("&state section is not found.")')
        stop
      end if

      istate                     = 0
      option%state_def           = 0.0d0
      option%state_weight        = 0.0d0
      option%state_weight_unnorm = 0.0d0
      do while (.true.)
        read(io,'(a)', end = 100) line

        if (line(1:1) /= "#") then
          istate = istate + 1

          read(line,*) ((option%state_def(imm, idm, istate), imm = 1, 2), &
                         idm = 1, option%ndim), option%state_weight(istate)
        end if

        if (istate == option%nstate) exit

      end do

100   continue

      write(iw,*)
      write(iw,'("Read_Ctrl_State>  States are defined as")')
      do istate = 1, option%nstate
        write(iw,'("State ", i5)') istate

        do idm = 1, option%ndim
          write(iw,'("  component ", i5, " : ", f15.7, " <===> ", f15.7)') &
            idm,                              &
            option%state_def(1, idm, istate), &
            option%state_def(2, idm, istate)
        end do

        write(iw,*)

      end do

      ! Weight
      !
      weight_sum = 0.0d0
      do istate = 1, nstate
        if (.not. option%is_initial(istate)) then
          option%state_weight(istate) = 0.0d0
          cycle
        end if
        weight_sum = weight_sum + option%state_weight(istate)
      end do
      option%state_weight_unnorm = option%state_weight
      option%state_weight        = option%state_weight / weight_sum 

      write(iw,*)
      write(iw,'("Read_Ctrl_State> Weight for each state")')
      do istate = 1, option%ninitial
        write(iw,'(i5,f15.7)') option%initial_state_ids(istate), &
                               option%state_weight(option%initial_state_ids(istate))
      end do
      write(iw,*)

    end subroutine read_ctrl_state
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine show_input(input)
!-----------------------------------------------------------------------
      implicit none

      type(s_input), intent(in) :: input

      ! Dummy
      !
      integer :: i

      ! Check
      !
      if (input%ncv == 0) then
        write(iw,'("Error. fcv should be specified.")')
        stop
      end if

      ! Print
      !
      write(iw,*)
      write(iw,'(">> Input section parameters")')
      do i = 1, input%ncv
        write(iw,'("fcv", 3x, i0, 3x, " = ", a)') i, trim(input%fcv(i)) 
      end do


    end subroutine show_input
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine show_output(output)
!-----------------------------------------------------------------------
      implicit none

      type(s_output), intent(in) :: output 


      ! Check
      !
      if (trim(output%fhead) == '') then
        write(iw,'("Error. fhead should be specified.")')
        stop
      end if

      ! Print
      !
      write(iw,*)
      write(iw,'(">> Output section parameters")')
      write(iw,'("fhead          = ", a)') trim(output%fhead)


    end subroutine show_output
!-----------------------------------------------------------------------

end module
!=======================================================================
