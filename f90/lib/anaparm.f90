!=======================================================================
! Module: mod_anaparm
!
!   Module for defining ANATRA formatted potential parameter (anaparm) file 
!
!   (c) Copyright 2024 Osaka Univ. All rights reserved.
!
!=======================================================================
module mod_anaparm

  use mod_util
  use mod_const

  implicit none

  ! parameters
  !

  ! structures
  !
  type :: s_anaparm
    integer :: natm
    real(8), allocatable :: sgm(:)
    real(8), allocatable :: eps(:)
    real(8), allocatable :: charge(:) 
  end type s_anaparm

  ! subroutines
  !
  public   :: read_anaparm 

  contains

!-----------------------------------------------------------------------
    subroutine read_anaparm(fanaparm, anaparm) 
!-----------------------------------------------------------------------
      implicit none

      character(len=MaxChar), intent(in)  :: fanaparm
      type(s_anaparm),        intent(out) :: anaparm

      integer                :: iunit
      integer                :: iatm, natm
      character(len=MaxChar) :: cdum


      iunit = UnitANAPARM
      open(iunit, file=trim(fanaparm))

      ! Get # of atoms (lines) 
      !
      natm = 0
      do while (.true.)
        read(iunit,*,end=100)
        natm = natm + 1 
      end do

100   rewind iunit

      ! Allocate memories 
      !
      anaparm%natm = natm
      allocate(anaparm%sgm(natm))
      allocate(anaparm%eps(natm))
      allocate(anaparm%charge(natm))

      ! Read parameters
      !
      do iatm = 1, natm
        read(iunit,*) cdum, anaparm%sgm(iatm), anaparm%eps(iatm), &
                      anaparm%charge(iatm)
      end do

      close(iunit)

    end subroutine read_anaparm
!-----------------------------------------------------------------------

end module mod_anaparm
!=======================================================================
