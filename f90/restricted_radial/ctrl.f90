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
  integer, parameter, public :: Nstate   = 1 
  integer, parameter, public :: REACTIVE = 1 
  integer, parameter, public :: OTHERS   = 2 
  integer, parameter, public :: StateInfo(Nstate) = (/1/)

  ! structures
  !
  !type :: s_cvinfo
  !  integer                             :: nfile     = 1
  !  character(len=MaxChar), allocatable :: fcv(:)
  !end type s_cvinfo

  type :: s_option
    logical              :: calcfe        = .false.
    logical              :: use_bootstrap = .false.
    logical              :: refs_system   = .false.
    logical              :: only_r        = .false.
    integer              :: ndim          = 1
    integer              :: ngrid         = 1000
    integer              :: nsta          = 1
    real(8)              :: dr            = 0.1d0
    real(8)              :: temperature   = 298.0d0
    real(8)              :: vol0          = 1661.0d0
    real(8)              :: box_ref(3)    = 0.0d0
    real(8)              :: urange(2)     = 0.0d0 
    real(8), allocatable :: state_def(:, :, :)
    real(8), allocatable :: react_range(:, :)
  end type s_option

  ! subroutines
  !
  public  :: read_ctrl
  private :: read_ctrl_option
  !private :: read_ctrl_cvinfo

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

      !open(iunit, file=trim(input%flist))
      !  call read_ctrl_cvinfo  (iunit, cvinfo)
      !close(iunit) 

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

      logical :: calcfe                     = .false.
      logical :: use_bootstrap              = .false.
      logical :: refs_system                = .false.
      logical :: only_r                     = .false.
      integer :: ndim                       = 1
      integer :: ngrid                      = 1000
      integer :: nsta                       = 1 
      real(8) :: dr                         = 0.1d0
      real(8) :: temperature                = 298.0d0
      real(8) :: vol0                       = 1661.0d0
      real(8) :: box_ref(3)                 = 0.0d0
      real(8) :: urange(2)                  = 0.0d0
      real(8) :: react_range(2, ndim_max)   = 0.0d0

      integer :: i, j
      integer :: iopt, ierr

      namelist /option_param/ calcfe, use_bootstrap, refs_system, only_r, ndim, ngrid, nsta, dr, &
                              temperature, vol0, box_ref, urange, react_range

      rewind iunit
      read(iunit, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("calcfe        = ", a)')  get_tof(calcfe)
      write(iw,'("use_bootstrap = ", a)')  get_tof(use_bootstrap)
      write(iw,'("refs_system   = ", a)')  get_tof(refs_system)
      write(iw,'("only_r        = ", a)')  get_tof(only_r)
      write(iw,'("ndim          = ", i0)') ndim
      write(iw,'("ngrid         = ", i0)') ngrid
      write(iw,'("nsta          = ", i0)') nsta
      write(iw,'("dr            = ", f20.10)') dr
      write(iw,'("temperature   = ", f20.10)') temperature 
      write(iw,'("vol0          = ", f20.10)') vol0
      write(iw,'("box_ref       = ", 3f20.10)') box_ref(1:3)
      write(iw,'("urange        = ", 2f20.10)') urange(1), urange(2) 
      write(iw,*)
      write(iw,'("react_range   =")')
      do i = 1, ndim
        write(iw,'(es15.7, " <= component ",i0," < ",es15.7)') &
          react_range(1, i), i, react_range(2, i) 
      end do

      !iopt = get_opt(mode, CoMMode, ierr)
      !if (ierr /= 0) then
      !  write(iw,'("Read_Ctrl_Option> Error.")')
      !  write(iw,'("mode = ",a," is not available.")') trim(mode)
      !  stop
      !end if
      !option%mode          = iopt
                          
      option%calcfe        = calcfe
      option%use_bootstrap = use_bootstrap
      option%refs_system   = refs_system
      option%only_r        = only_r 
      option%ndim          = ndim
      option%ngrid         = ngrid
      option%nsta          = nsta
      option%dr            = dr 
      option%temperature   = temperature 
      option%vol0          = vol0
      option%box_ref       = box_ref
      option%urange        = urange

      allocate(option%react_range(1:2, ndim))
      allocate(option%state_def(2, ndim, Nstate))

      do i = 1, ndim
        option%react_range(1:2, i)   = react_range(1:2, i)
      end do

      do i = 1, ndim
        option%state_def(1:2, i, 1) = react_range(1:2, i)
      end do

      !do i = 1, 3
      !  do j = 1, ndim
      !    write(iw,'(2f15.7)') option%state_def(1, j, i), &
      !                         option%state_def(2, j, i)
      !  end do
      !end do
      ! Combination check
      !

    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

end module
!=======================================================================
