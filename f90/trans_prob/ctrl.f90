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
  integer,      parameter, public :: MaxStates    =  100 
  integer,      parameter, public :: NotSpecified = -100

  integer,      parameter, public :: KineticModeTransition = 1
  integer,      parameter, public :: KineticModeReaction   = 2
  character(*), parameter, public :: KineticModes(2) = (/'TRANSITION', &
                                                         'REACTION  '/) 

  ! structures
  !
  type :: s_option
    logical :: use_bootstrap        = .false.
    logical :: use_reflection_state = .false.
    logical :: use_product_state    = .false.
    logical :: use_dissociate_state = .false.
    logical :: write_Rij_bin        = .false. 
    logical :: read_Rij_bin         = .false.
    logical :: write_Kijk_bin       = .false.
    logical :: read_Kijk_bin        = .false. 
    logical :: extrapolate          = .false.
    logical :: check_Kijk           = .false.

    integer :: kinetic_mode  = KineticModeTransition

    integer :: nmol                               = NotSpecified 
    integer :: ndim                               = NotSpecified 
    integer :: nstate                             = NotSpecified
    integer :: out_id_rij                         = NotSpecified
    integer :: out_id_kijk                        = NotSpecified 
    integer :: reflection_state_ids(MaxStates)    = NotSpecified 
    integer :: product_state_ids   (MaxStates)    = NotSpecified
    integer :: dissociate_state_ids(MaxStates)    = NotSpecified
    integer :: initial_state_ids   (MaxStates)    = NotSpecified
    integer :: sel_ijk             (3, MaxStates) = NotSpecified

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
    logical, allocatable :: is_dissoc   (:) 
    real(8), allocatable :: state_def   (:, :, :)
    real(8), allocatable :: state_weight(:)
    real(8), allocatable :: state_weight_unnorm(:)

  end type s_option

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
    subroutine read_ctrl(input, output, option, bootopt, timegrid)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),    intent(out) :: input
      type(s_output),   intent(out) :: output
      type(s_option),   intent(out) :: option 
      type(s_bootopt),  intent(out) :: bootopt 
      type(s_timegrid), intent(out) :: timegrid 

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

      if (option%use_bootstrap) then
        call read_ctrl_bootstrap(io, bootopt)
      end if

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

      logical :: use_bootstrap        = .false.
      logical :: use_reflection_state = .false.
      logical :: use_product_state    = .false.
      logical :: use_dissociate_state = .false.
      logical :: write_Kijk_bin       = .false.
      logical :: read_Kijk_bin        = .false.
      logical :: write_Rij_bin        = .false. 
      logical :: read_Rij_bin         = .false. 
      logical :: extrapolate          = .false.
      logical :: check_Kijk           = .false.

      character(len=MaxChar) :: kinetic_mode     = 'TRANSITION' 
      
      integer :: nmol                            = NotSpecified
      integer :: ndim                            = NotSpecified
      integer :: nstate                          = NotSpecified
      integer :: out_id_rij                      = NotSpecified
      integer :: out_id_kijk                     = NotSpecified
      integer :: reflection_state_ids(MaxStates) = NotSpecified
      integer :: product_state_ids(MaxStates)    = NotSpecified
      integer :: dissociate_state_ids(MaxStates) = NotSpecified
      integer :: initial_state_ids(MaxStates)    = NotSpecified
      integer :: sel_ijk(3, MaxStates)           = NotSpecified
      real(8) :: dt
      real(8) :: t_sparse
      real(8) :: t_range
      real(8) :: t_extend
      real(8) :: dt_tcfout

      ! Local
      !
      integer :: nreflect = 0, nproduct = 0, ndissoc = 0, ninitial = 0, nselect = 0

      ! Dummy
      !
      integer :: i, j, it
      integer :: iopt, ierr


      namelist /option_param/ &
        use_bootstrap,        &
        use_reflection_state, &
        use_product_state,    &
        use_dissociate_state, &
        write_Kijk_bin,       &
        read_Kijk_bin,        &
        write_Rij_bin,        &
        read_Rij_bin,         &
        extrapolate,          &
        check_Kijk,           &
        kinetic_mode,         &
        nmol,                 &
        ndim,                 &
        nstate,               &
        out_id_rij,           &
        out_id_kijk,          &
        reflection_state_ids, &
        product_state_ids,    &
        dissociate_state_ids, &
        initial_state_ids,    &
        sel_ijk,              &
        dt,                   &
        t_sparse,             &
        t_range,              &
        t_extend,             &
        dt_tcfout

      rewind io 
      read(io, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("use_bootstrap        = ", a)')     get_tof(use_bootstrap)
      write(iw,'("kinetic_mode         = ", a)')     trim(kinetic_mode)

      if (trim(kinetic_mode) == 'REACTION') then
        write(iw,'("use_reflection_state = ", a)')   get_tof(use_reflection_state)
        write(iw,'("use_product_state    = ", a)')   get_tof(use_product_state)
        write(iw,'("use_dissociate_state = ", a)')   get_tof(use_dissociate_state)
        write(iw,'("write_Kijk_bin       = ", a)')   get_tof(write_Kijk_bin)
        write(iw,'("read_Kijk_bin        = ", a)')   get_tof(read_Kijk_bin)
        write(iw,'("write_Rij_bin        = ", a)')   get_tof(write_Rij_bin)
        write(iw,'("read_Rij_bin         = ", a)')   get_tof(read_Rij_bin)
        write(iw,'("check_Kijk           = ", a)')   get_tof(check_Kijk)
      end if

      write(iw,'("extrapolate          = ", a)')     get_tof(extrapolate)
      write(iw,'("nmol                 = ", i0)')    nmol
      write(iw,'("ndim                 = ", i0)')    ndim
      write(iw,'("nstate               = ", i0)')    nstate

      if (write_Kijk_bin) then
        if (out_id_kijk == NotSpecified) then
          out_id_kijk = 1
        end if
        write(iw,'("out_id_kijk          = ", i0)')  out_id_kijk 
      end if

      if (write_Rij_bin) then
        if (out_id_rij == NotSpecified) then
          out_id_rij = 1
        end if
        write(iw,'("out_id_rij           = ", i0)')    out_id_rij 
      end if

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
      if (use_reflection_state) then
        nreflect = 0
        write(iw,*)
        do i = 1, MaxStates
          if (reflection_state_ids(i) /= NotSpecified)  then
            write(iw,'("reflection_state_ids  ", i0, " : ", i0)') reflection_state_ids(i)
            nreflect = nreflect + 1
          else
            exit
          end if
        end do
      end if

      ! Product
      !
      allocate(option%is_product(option%nstate))
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
      allocate(option%is_dissoc(option%nstate))
      option%is_dissoc = .false.

      if (use_dissociate_state) then
        ndissoc = 0
        write(iw,*)
        do i = 1, MaxStates
          if (dissociate_state_ids(i) /= NotSpecified)  then
            write(iw,'("dissociate_state_ids  ", i0, " : ", i0)')i,  dissociate_state_ids(i)
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

      ! Specific element
      !
      if (write_Kijk_bin .or. read_Kijk_bin) then
        nselect = 0
        do i = 1, MaxStates
          if (sel_ijk(1, i) /= NotSpecified) then
            nselect = nselect + 1
            write(iw,'("sel_ijk      = ", 3(i0,2x))') sel_ijk(1:3, i)
          end if
        end do
        option%nselect = nselect 
      end if

      iopt = get_opt(kinetic_mode, KineticModes, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("kinetic_mode = ",a," is not available.")') trim(kinetic_mode)
      end if
      option%kinetic_mode = iopt

      option%use_bootstrap        = use_bootstrap
      option%use_reflection_state = use_reflection_state
      option%use_product_state    = use_product_state
      option%use_dissociate_state = use_dissociate_state
      option%write_Kijk_bin       = write_Kijk_bin
      option%read_Kijk_bin        = read_Kijk_bin
      option%write_Rij_bin        = write_Rij_bin
      option%read_Rij_bin         = read_Rij_bin
      option%extrapolate          = extrapolate
      option%check_Kijk           = check_Kijk
      option%nmol                 = nmol
      option%ndim                 = ndim
      option%nstate               = nstate
      option%out_id_rij           = out_id_rij
      option%out_id_kijk          = out_id_kijk 

      option%nreflect             = nreflect
      option%reflection_state_ids = reflection_state_ids

      option%nproduct             = nproduct
      option%product_state_ids    = product_state_ids

      option%ndissoc              = ndissoc
      option%dissociate_state_ids = dissociate_state_ids

      option%ninitial             = ninitial
      option%initial_state_ids    = initial_state_ids

      option%sel_ijk              = sel_ijk

      option%dt                   = dt
      option%t_sparse             = t_sparse
      option%t_range              = t_range
      option%t_extend             = t_extend
      option%dt_tcfout            = dt_tcfout

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
