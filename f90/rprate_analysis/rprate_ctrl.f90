!=======================================================================
module mod_ctrl
!=======================================================================
  use mod_util
  use mod_const
  implicit none

  ! constants
  !

  ! structures
  !
  type :: s_option
    integer :: rporder     = 0
    logical :: use_sfe     = .false.
    logical :: use_shift   = .false.

    real(8) :: timeunit    = 0.0d0
    real(8) :: kins        = 0.0d0
    real(8) :: kins_err    = 0.0d0
    real(8) :: taud        = 0.0d0
    real(8) :: taud_err    = 0.0d0
    real(8) :: Kstar       = 0.0d0
    real(8) :: Kstar_err   = 0.0d0
    real(8) :: taupa3      = 0.0d0
    real(8) :: taupa3_err  = 0.0d0

    real(8) :: temperature = 0.0d0
    real(8) :: sfeb        = 0.0d0
    real(8) :: sfeb_err    = 0.0d0
    real(8) :: sfed        = 0.0d0
    real(8) :: sfed_err    = 0.0d0
    real(8) :: dGcorr      = 0.0d0

    real(8) :: shift_sol     = 0.0d0
    real(8) :: shift_sol_err = 0.0d0 
    real(8) :: shift_ref     = 0.0d0
    real(8) :: shift_ref_err = 0.0d0

  end type s_option

  ! subroutines
  !
  public  :: read_ctrl
  private :: read_ctrl_option

  contains

!-----------------------------------------------------------------------
    subroutine read_ctrl(option)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),  intent(out) :: option 

      integer                :: iunit
      character(len=MaxChar) :: f_ctrl

      ! get control file name
      !
      call getarg(1, f_ctrl)

      call open_file(trim(f_ctrl), iunit)
      write(iw,*)
      write(iw,'("Read_Ctrl> Reading parameters from ", a)') trim(f_ctrl)
      open(iunit, file=trim(f_ctrl), status='old')
        call read_ctrl_option (iunit, option)
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

      integer :: rporder       = 0
      logical :: use_sfe       = .false.
      logical :: use_shift     = .false.
      real(8) :: timeunit      = 0.0d0
      real(8) :: kins          = 0.0d0
      real(8) :: kins_err      = 0.0d0
      real(8) :: taud          = 0.0d0
      real(8) :: taud_err      = 0.0d0
      real(8) :: taupa3        = 0.0d0
      real(8) :: taupa3_err    = 0.0d0
      real(8) :: Kstar         = 0.0d0
      real(8) :: Kstar_err     = 0.0d0
      real(8) :: temperature   = 0.0d0
      real(8) :: sfeb          = 0.0d0
      real(8) :: sfeb_err      = 0.0d0
      real(8) :: sfed          = 0.0d0
      real(8) :: sfed_err      = 0.0d0
      real(8) :: dGcorr        = 0.0d0
      real(8) :: shift_sol     = 0.0d0
      real(8) :: shift_sol_err = 0.0d0
      real(8) :: shift_ref     = 0.0d0
      real(8) :: shift_ref_err = 0.0d0

      namelist /option_param/ rporder,       &
                              use_sfe,       &
                              use_shift,     &
                              timeunit,      &
                              kins,          &
                              kins_err,      &
                              taud,          &
                              taud_err,      &
                              Kstar,         &
                              Kstar_err,     &
                              taupa3,        &
                              taupa3_err,    &
                              temperature,   &
                              sfeb,          &
                              sfeb_err,      &
                              sfed,          &
                              sfed_err,      &
                              dGcorr,        &
                              shift_sol,     &
                              shift_sol_err, &
                              shift_ref,     &
                              shift_ref_err


      ! Initialize
      !
      use_sfe = .false. 

      ! Read
      !
      rewind iunit
      read(iunit, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("rporder         = ", i0)')     rporder
      write(iw,'("use_sfe         = ", a)')      get_tof(use_sfe)
      write(iw,'("use_shift       = ", a)')      get_tof(use_shift)
      write(iw,'("timeunit        = ", es15.7)') timeunit 
      write(iw,'("kins            = ", f15.7)')  kins
      write(iw,'("kins_err        = ", f15.7)')  kins_err
      write(iw,'("taud            = ", f15.7)')  taud
      write(iw,'("taud_err        = ", f15.7)')  taud_err

      if (rporder == 2) then
        write(iw,*)
        write(iw,'("taupa3          = ", f15.7)')  taupa3
        write(iw,'("taupa3_err      = ", f15.7)')  taupa3_err
      end if

      if (use_sfe) then
        write(iw,*)
        write(iw,'("temperature     = ", f15.7)')  temperature
        write(iw,'("sfeb            = ", f15.7)')  sfeb
        write(iw,'("sfeb_err        = ", f15.7)')  sfeb_err
        write(iw,'("sfed            = ", f15.7)')  sfed
        write(iw,'("sfed_err        = ", f15.7)')  sfed_err
        write(iw,'("dGcorr          = ", f15.7)')  dGcorr
        if (use_shift) then
          write(iw,'("shift_sol       = ", f15.7)')  shift_sol 
          write(iw,'("shift_sol_err   = ", f15.7)')  shift_sol_err
          write(iw,'("shift_ref       = ", f15.7)')  shift_ref
          write(iw,'("shift_ref_err   = ", f15.7)')  shift_ref_err
        end if
      else
        write(iw,*)
        write(iw,'("Kstar           = ", f15.7)')  Kstar 
        write(iw,'("Kstar_err       = ", f15.7)')  Kstar_err
      end if


      option%rporder       = rporder
      option%use_sfe       = use_sfe
      option%use_shift     = use_shift
      option%timeunit      = timeunit
      option%kins          = kins
      option%kins_err      = kins_err
      option%taud          = taud
      option%taud_err      = taud_err
      option%taupa3        = taupa3
      option%taupa3_err    = taupa3_err

      option%Kstar         = Kstar
      option%Kstar_err     = Kstar_err

      option%temperature   = temperature
      option%sfeb          = sfeb
      option%sfeb_err      = sfeb_err
      option%sfed          = sfed
      option%sfed_err      = sfed_err
      option%dGcorr        = dGcorr
      option%shift_sol     = shift_sol
      option%shift_sol_err = shift_sol_err
      option%shift_ref     = shift_ref
      option%shift_ref_err = shift_ref_err

    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

end module
!=======================================================================
