!-------------------------------------------------------------------------------
!
!  Module   mod_xyz
!
!    Module for reading xyz format files 
!
!    (c) Copyright 2024 Osaka Univ. All rights reserved.
!
!-------------------------------------------------------------------------------

module mod_xyz

  use mod_util
  use mod_const

  implicit none

  ! parameters
  !

  ! structures
  !
  type :: s_xyz
    integer                       :: natm
    character(len=4), allocatable :: atmname(:)
    real(8),          allocatable :: coord(:, :) 
  end type s_xyz
  
  ! subroutines
  !
  public :: read_xyz 
  public :: write_xyz 

  contains

    !---------------------------------------------------------------------------
    !
    !  Subroutine  read_xyz
    !  @brief      read xyz format file 
    !  @authors    KK
    !  @param[in]  iunit   : file unit number for ctrl file
    !  @param[out] xyz     : structure of xyz structure 
    !
    !---------------------------------------------------------------------------

    subroutine read_xyz(iunit, xyz)
      implicit none

      integer,           intent(in)  :: iunit
      type(s_xyz),       intent(out) :: xyz 

      integer :: i, j, irank


      rewind iunit
      read(iunit,*) xyz%natm
      read(iunit,*)

      allocate(xyz%coord(3, xyz%natm), xyz%atmname(xyz%natm))

      do i = 1, xyz%natm
        read(iunit,*) xyz%atmname(i), (xyz%coord(j, i), j = 1, 3)
      end do

    end subroutine read_xyz

    !---------------------------------------------------------------------------
    !
    !  Subroutine  write_xyz
    !  @brief      read xyz format file 
    !  @authors    KK
    !  @param[in]  iunit   : file unit number for ctrl file
    !  @param[in]  xyz     : structure of xyz structure 
    !
    !---------------------------------------------------------------------------

    subroutine write_xyz(iunit, xyz)
      implicit none

      integer,           intent(in) :: iunit
      type(s_xyz),       intent(in) :: xyz 

      integer :: i, j, irank

      write(iunit,'(i0)') xyz%natm
      write(iunit,*)

      do i = 1, xyz%natm
        write(iunit,'(a4,3(e20.10,2x))') xyz%atmname(i), (xyz%coord(j, i), j = 1, 3)
      end do

    end subroutine write_xyz

end module mod_xyz
!=======================================================================
