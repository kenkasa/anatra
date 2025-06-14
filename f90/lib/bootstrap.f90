!=======================================================================
! Module: Mod_Bootstrap
!
!   Module for using bootstrap structure 
!   and read parameters from ctrl file
!
!   (c) Copyright 2024 Osaka Univ. All rights reserved.
!
!=======================================================================
module mod_bootstrap

  use mod_util
  use mod_const

  implicit none

  ! parameters
  !
  integer,      parameter, public :: BootTypeTRAJ = 1
  integer,      parameter, public :: BootTypeSNAP = 2
  character(*), parameter, public :: BootTypes(2) = (/'TRAJ  ', &
                                                      'SNAP  '/)

  ! structures
  !
  type :: s_bootopt
    logical :: duplicate  = .true.
    integer :: boottype   = BootTypeSNAP
    integer :: iseed      = 0
    integer :: nsample    = 1000
    integer :: ntrial     = 1000
  end type s_bootopt

  ! subroutines
  !
  public :: read_ctrl_bootstrap

  contains
!-----------------------------------------------------------------------
    subroutine read_ctrl_bootstrap(iunit, bootopt)
!-----------------------------------------------------------------------
      implicit none
!
      integer,           intent(in)  :: iunit
      type(s_bootopt),   intent(out) :: bootopt

      integer                :: iseed       = -1
      character(len=MaxChar) :: boottype    = 'TRAJ'
      logical                :: duplicate   = .true.
      integer                :: nsample     = 1000
      integer                :: ntrial      = 1000

      integer :: i, j
      integer :: iopt, ierr

      namelist /bootopt_param/ &
        iseed,                 &
        boottype,              &
        duplicate,             &
        nsample,               &
        ntrial


      rewind iunit
      read(iunit,bootopt_param)

      write(iw,*)
      write(iw,'(">> Bootstrap section parameters")')
      write(iw,'("iseed     = ",i0)') iseed
      write(iw,'("boottype  = ",a)')  trim(boottype)
      write(iw,'("duplicate = ",a)')  get_tof(duplicate)
      write(iw,'("nsample   = ",i0)') nsample
      write(iw,'("ntrial    = ",i0)') ntrial


      iopt = get_opt(boottype, BootTypes, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Bootstrap> Error.")')
        write(iw,'("boottype = ", a , " is not available.")') trim(boottype)
        stop
      end if
      bootopt%boottype  = iopt

      bootopt%iseed     = iseed
      bootopt%duplicate = duplicate
      bootopt%nsample   = nsample
      bootopt%ntrial    = ntrial

    end subroutine read_ctrl_bootstrap
!-----------------------------------------------------------------------
!

end module mod_bootstrap
!=======================================================================
