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
  integer, parameter, public :: MaxDim   = 3
  integer, parameter, public :: MaxNcell = 100 

  ! structures
  !

  type :: s_option
    integer :: ndim                   = MaxDim
    integer :: xyzcol(MaxDim)         = (/2, 3, 4/) 

    ! use voronoi tessellation (ugly implementation) 
    logical :: use_voronoi            = .true.
    integer :: vr_ncell               = 5
    real(8) :: vr_normvec(3)          = (/1.0d0, 1.0d0, 1.0d0/)
    real(8), allocatable :: vr_cellpos(:, :)

    ! use grid information
    logical :: use_grids              = .false.
    logical :: pbc                    = .false.
    real(8) :: box(3)                 = 0.0d0

  end type s_option

  ! subroutines
  !
  public  :: read_ctrl
  private :: read_ctrl_option

  contains

!-----------------------------------------------------------------------
    subroutine read_ctrl(input, output, option)
!-----------------------------------------------------------------------
      implicit none

      integer, parameter           :: iunit = 10

      type(s_input),   intent(out) :: input
      type(s_output),  intent(out) :: output
      type(s_option),  intent(out) :: option

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

      close(iunit)

    end subroutine read_ctrl
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine read_ctrl_option(iunit, option)
!-----------------------------------------------------------------------
      implicit none
!
      integer,        intent(in)  :: iunit
      type(s_option), intent(out) :: option 

      integer :: ndim                    = 2
      integer :: xyzcol(MaxDim)          = (/1, 2, 3/)
      ! voronoi
      logical :: use_voronoi             = .true.
      integer :: vr_ncell                = 0
      real(8) :: vr_normvec(3)           = 1.0d0
      real(8) :: vr_cellpos(3*MaxNcell)  = 0.0d0
      ! use_grids
      logical :: use_grids               = .false.
      logical :: pbc                     = .false.
      real(8) :: box(3)                  = 0.0d0

      integer :: i, j, k 


      namelist /option_param/ ndim,        &
                              xyzcol,      &
                              use_voronoi, &
                              vr_ncell,    &
                              vr_normvec,  &
                              vr_cellpos,  &
                              use_grids,   &
                              pbc,         &
                              box

      rewind iunit
      read(iunit, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("ndim              = ", i0)')              ndim 
      write(iw,'("xyzcol            = ", 3(i0,2x))')        (xyzcol(i), i = 1, ndim) 

      write(iw,'("use_voronoi       = ", a)')               get_tof(use_voronoi)
      if (use_voronoi) then
        write(iw,'("vr_ncell          = ", i0)')            vr_ncell
        write(iw,'("vr_normvec        = ", 3(f15.7,2x))')   (vr_normvec(i),    i = 1, ndim)

        write(iw,'("vr_cellpos : ")')
        j = 0
        do i = 1, vr_ncell
          write(iw,'(i3,3f20.10)') i, (vr_cellpos(j + k), k = 1, ndim)
          j = j + ndim
        end do
        write(iw,*)
      end if
      write(iw,'("use_grids         = ", a)')               get_tof(use_grids)
      write(iw,'("pbc               = ", a)')               get_tof(pbc)
      write(iw,'("box               = ", 3(f15.7))')        box(1:3) 

      ! memory allocation
      !
      if (use_voronoi) then
        allocate(option%vr_cellpos(ndim, vr_ncell))
      end if

      ! setup option variables
      !
      option%xyzcol = 0

      option%ndim                        = ndim
      !option%nstart                      = nstart
      option%xyzcol(1:ndim)              = xyzcol(1:ndim)
      option%use_voronoi                 = use_voronoi
      option%use_grids                   = use_grids
      option%pbc                         = pbc
      option%box                         = box

      if (use_voronoi) then
        option%vr_ncell                    = vr_ncell
        option%vr_normvec                  = vr_normvec

        k = 0
        do i = 1, vr_ncell
          do j = 1, ndim
            k = k + 1
            option%vr_cellpos(j, i) = vr_cellpos(k)
          end do
        end do
        !option%vr_cellpos(1:3, 1:vr_ncell) = vr_cellpos(1:3, 1:vr_ncell)
      end if
      
      ! combination check
      !
      if (option%use_voronoi) then
        !if (option%ndim > 2) then
        !  write(iw,'("Read_Ctrl_Option> Error.")')
        !  write(iw,'("ndim > 2 is currently not supported.")')
        !  write(iw,'("sorry for the inconvenience.")')
        !  stop
        !end if

        !if (option%ndim /= 2) then
        !  write(iw,'("Read_Ctrl_Option> Error.")')
        !  write(iw,'("voronoi tessellation is available &
        !             &only if ndim = 2")')
        !  stop
        !end if
      end if

      if (option%use_voronoi .and. option%use_grids) then
        write(iw,'("Read_Ctrol_Option>")')
        write(iw,'("use_voronoi = .true. and use_grids = .true. &
                   &can not be used at the same time.")')
        stop
      end if

      if (option%use_grids) then
        if (option%ndim /= 3) then
          write(iw,'("Read_Ctrl_Option>")')
          write(iw,'("use_grids = .true. is available &
                     &only if ndim = 3")')
          stop
        end if
      end if

    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end module
!=======================================================================
