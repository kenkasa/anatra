!-------------------------------------------------------------------------------
!
!  Module   mod_traj  
!    Module for defining trajectory parameters and read DCD trajectory  
!
!    (c) Copyright 2024 Osaka Univ. All rights reserved.
!
!-------------------------------------------------------------------------------

module mod_traj

  use mod_util
  use mod_const
  use mod_netcdfio
  use mod_dcdio
  use mod_xtcio
  use xdr, only: xtcfile

  implicit none

  ! parameters
  !
  integer, parameter, public      :: TrajTypeDCD  = 1
  integer, parameter, public      :: TrajTypeXTC  = 2
  integer, parameter, public      :: TrajTypeNCD  = 3
  character(*), parameter, public :: TrajTypes(3) = (/'DCD  ',&
                                                      'XTC  ',&
                                                      'NCD  '/)

  ! structures
  !
  type :: s_trajopt
    real(8)                       :: dt
    integer                       :: trajtype_in       = TrajTypeDCD
    integer                       :: trajtype_out      = TrajTypeXTC
    character(len=MaxChar)        :: molinfo(MaxTraj)  = "molinfo"
    character(len=MaxChar)        :: molinfo_refu      = "molinfo"
    character(len=MaxChar)        :: molinfo_refv      = "molinfo"

    integer                       :: nmolinfo 
  end type s_trajopt
  
  type :: s_traj
    integer(4)                    :: dcdinfo(20)
    integer                       :: natm
    integer                       :: nstep
    real(8)                       :: dt
    character(len=4), allocatable :: atmname(:)
    character(len=4), allocatable :: resname(:)
    character(len=4), allocatable :: segname(:)
    integer,          allocatable :: resid(:)
    real(8),          allocatable :: box(:, :)
    real(8),          allocatable :: coord(:,:,:)
    real(8),          allocatable :: mass(:)
    real(8),          allocatable :: charge(:) 
    integer,          allocatable :: ind(:)
    integer,          allocatable :: ind_ref(:)
  end type s_traj

  ! subroutines
  !
  public  :: read_ctrl_trajopt
  public  :: setup_traj
  public  :: setup_traj_from_args
  public  :: get_trajtype
  private :: alloc_traj
  public  :: dealloc_traj

  contains

    !---------------------------------------------------------------------------
    !
    !  Subroutine  read_ctrl_trajopt
    !  @brief      read trajopt variables from ctrl file
    !  @authors    KK
    !  @param[in]  iunit   : file unit number for ctrl file
    !  @param[out] trajopt : structure of trajectory option information 
    !
    !---------------------------------------------------------------------------

    subroutine read_ctrl_trajopt(iunit, trajopt, myrank)
      implicit none

      integer,           intent(in)  :: iunit
      type(s_trajopt),   intent(out) :: trajopt
      integer, optional, intent(in)  :: myrank

      real(8)                      :: dt               = 0.0d0
      character(len=MaxChar)       :: molinfo(MaxTraj) = ""
      character(len=MaxChar)       :: molinfo_refu     = ""
      character(len=MaxChar)       :: molinfo_refv     = ""

      integer :: i, irank
      integer :: nmolinfo

      namelist /trajopt_param/ dt, molinfo, molinfo_refu, molinfo_refv 

      molinfo = ""
      rewind iunit
      read(iunit, trajopt_param)

      if (present(myrank)) then
        irank = myrank
      else
        irank = 0
      end if

      nmolinfo = 0
      do i = 1, MaxTraj
        if (trim(molinfo(i)) /= "") then
          nmolinfo = nmolinfo + 1
        end if
      end do

      if (irank == 0) then
        write(iw,*)
        write(iw,'(">> Trajopt section parameters")')
        write(iw,'("dt          = ", e15.7)') dt
        do i = 1, MaxTraj
          if (trim(molinfo(i)) /= "") then
            write(iw,'("molinfo     = ", a)') trim(molinfo(i))
          end if
        end do
      end if
      if (trim(molinfo_refu) /= "") then
        write(iw,'("molinfo_refu = ", a)') trim(molinfo_refu)
        write(iw,'("molinfo_refv = ", a)') trim(molinfo_refv)
      end if

      trajopt%dt           = dt
      trajopt%molinfo      = molinfo 
      trajopt%molinfo_refu = molinfo_refu
      trajopt%molinfo_refv = molinfo_refv

      trajopt%nmolinfo     = nmolinfo


    end subroutine read_ctrl_trajopt

    !---------------------------------------------------------------------------
    !
    !  Subroutine  setup_traj 
    !  @brief      setup ANATRA trajectory from DCD trajectory 
    !  @authors    KK
    !  @param[in]  trajopt : structure of trajectory option information 
    !  @param[in]  dcd     : structure of dcd trajectory 
    !                                 (already defined)
    !  @param[out] traj    : structure of ANATRA trajectory
    !  @param[in, optional]  trajid : trajectory ID 
    !
    !---------------------------------------------------------------------------

    subroutine setup_traj(trajopt, dcd, traj, trajid)

      implicit none

      integer, parameter :: ncolmax = 7 

      ! formal arguments
      type(s_trajopt),            intent(in)  :: trajopt
      type(s_dcd),                intent(in)  :: dcd
      type(s_traj),               intent(out) :: traj
      integer,          optional, intent(in)  :: trajid

      ! local variables
      integer                :: iatm, icol, ncol, id, ierr
      character(len=MaxChar) :: col(ncolmax)
      character(len=MaxChar) :: line


      ! specify trajectory id if trajid is present
      ! (default: id = 1)
      !
      if (present(trajid)) then
        id = trajid 
      else
        id = 1
      end if

      ! allocate traj structure
      !
      call alloc_traj(dcd%natm, dcd%nstep, traj)

      ! setup traj variables
      !
      traj%dcdinfo = dcd%dcdinfo
      traj%dt      = trajopt%dt
      traj%nstep   = dcd%nstep
      traj%natm    = dcd%natm
      traj%box     = dcd%box
      traj%coord   = dcd%coord

      ! read molinfo
      !
      open(11,file=trim(trajopt%molinfo(id)))
        ! check # of columns
        !
        ncol = ncolmax + 1
100     ncol = ncol - 1
        if (ncol == 0) then
          write(iw,'("Setup_Traj> Error.")')
          write(iw,'("molinfo file is empty.")')
          stop
        end if

        rewind 11
        read(11,'(a)') line
        read(line,*,iostat=ierr) (col(icol), icol = 1, ncol)
        if (ierr /= 0) then
          go to 100
        end if

        rewind 11
        if (ncol == 4) then
          do iatm = 1, dcd%natm
            read(11,*) traj%resid(iatm), traj%resname(iatm), &
                       traj%atmname(iatm), traj%mass(iatm)
          end do
        else if (ncol == 5) then
          do iatm = 1, dcd%natm
            read(11,*) traj%resid(iatm),   traj%resname(iatm), &
                       traj%atmname(iatm), traj%mass(iatm),    &
                       traj%charge(iatm)
          end do
        else if (ncol == 6) then
          do iatm = 1, dcd%natm
            read(11,*) traj%resid(iatm),   traj%resname(iatm), &
                       traj%atmname(iatm), traj%mass(iatm),    &
                       traj%charge(iatm),  traj%ind(iatm)
          end do
        else if (ncol == 7) then
          do iatm = 1, dcd%natm
            read(11,*) traj%resid(iatm),   traj%resname(iatm), &
                       traj%atmname(iatm), traj%mass(iatm),    &
                       traj%charge(iatm),  traj%ind(iatm),     &
                       traj%segname(iatm)
          end do
        end if
      close(11)

    end subroutine setup_traj

    !---------------------------------------------------------------------------
    !
    !  Subroutine  setup_traj_from_args 
    !  @brief      setup ANATRA trajectory from arguments
    !  @authors    KK
    !  @param[in]  trajopt : structure of trajectory option information 
    !  @param[in]  nstep   : # of time steps 
    !  @param[out] traj    : structure of ANATRA trajectory
    !  @param[in, optional]  trajid : trajectory ID 
    !
    !---------------------------------------------------------------------------

    subroutine setup_traj_from_args(trajopt, nstep, traj, trajid, refinfo)

      implicit none

      integer, parameter :: ncolmax = 7 

      ! formal arguments
      type(s_trajopt),            intent(in)  :: trajopt
      integer,                    intent(in)  :: nstep
      type(s_traj),               intent(out) :: traj
      integer,          optional, intent(in)  :: trajid
      character(len=4), optional, intent(in)  :: refinfo

      ! local variables
      integer                :: iatm, icol, ncol, id, ierr
      character(len=MaxChar) :: col(ncolmax)
      character(len=MaxChar) :: line


      ! specify trajectory id if trajid is present
      ! (default: id = 1)
      !
      if (present(trajid)) then
        id = trajid 
      else
        id = 1
      end if

      ! get # of atoms
      !
      open(11,file=trim(trajopt%molinfo(id)))
        iatm = 0
        do while(.true.) 
          read(11, *, end = 99)
          iatm = iatm + 1
        end do
99      rewind 11
        traj%natm = iatm
      close(11)

      ! allocate traj structure
      !
      call alloc_traj(traj%natm, nstep, traj)

      ! setup traj variables
      !
      traj%dt      = trajopt%dt
      traj%nstep   = nstep
      traj%box     = 0.0d0 
      traj%coord   = 0.0d0 


      ! read molinfo
      !
      open(11,file=trim(trajopt%molinfo(id)))
        ! check # of columns
        !
        ncol = ncolmax + 1
100     ncol = ncol - 1
        if (ncol == 0) then
          write(iw,'("Setup_Traj> Error.")')
          write(iw,'("molinfo file is empty.")')
          stop
        end if

        rewind 11
        read(11,'(a)') line
        read(line,*,iostat=ierr) (col(icol), icol = 1, ncol)
        if (ierr /= 0) then
          go to 100
        end if

        rewind 11
        if (ncol == 4) then
          do iatm = 1, traj%natm
            read(11,*) traj%resid(iatm), traj%resname(iatm), &
                       traj%atmname(iatm), traj%mass(iatm)
          end do
        else if (ncol == 5) then
          do iatm = 1, traj%natm
            read(11,*) traj%resid(iatm),   traj%resname(iatm), &
                       traj%atmname(iatm), traj%mass(iatm),    &
                       traj%charge(iatm)
          end do
        else if (ncol == 6) then
          do iatm = 1, traj%natm
            read(11,*) traj%resid(iatm),   traj%resname(iatm), &
                       traj%atmname(iatm), traj%mass(iatm),    &
                       traj%charge(iatm),  traj%ind(iatm)
          end do
        else if (ncol == 7) then
          do iatm = 1, traj%natm
            read(11,*) traj%resid(iatm),   traj%resname(iatm), &
                       traj%atmname(iatm), traj%mass(iatm),    &
                       traj%charge(iatm),  traj%ind(iatm),     &
                       traj%segname(iatm)
          end do
        end if
      close(11)

      if (present(refinfo)) then
        if (refinfo == "refu") then
          open(11,file=trim(trajopt%molinfo_refu))
            do iatm = 1, traj%natm
              read(11,*) col(1), col(2),               &
                         col(3), col(4),               &
                         col(5), traj%ind_ref(iatm),   &
                         col(7) 
            end do
          close(11)
        else if (refinfo == "refv") then
          open(11,file=trim(trajopt%molinfo_refv))
            do iatm = 1, traj%natm
              read(11,*) col(1), col(2),               &
                         col(3), col(4),               &
                         col(5), traj%ind_ref(iatm),   &
                         col(7) 
            end do

          close(11)
        end if
      end if


    end subroutine setup_traj_from_args

!-------------------------------------------------------------------------------
    subroutine get_trajtype(fname, trajtype)
!-------------------------------------------------------------------------------
      implicit none

      character(*), intent(in)  :: fname
      integer,      intent(out) :: trajtype

      character(len=MaxChar)    :: ext
      integer                   :: ierr


      call get_file_extention(fname, ext, ierr)

      if (ierr /= 0) then
        write(iw,'("Get_Trajtype> Error.")')
        write(iw,'("File extention of trajectory is not found.")')
        write(iw,'(a)') trim(ext)
        stop
      end if

      if (trim(ext) == "dcd") then
        trajtype = TrajTypeDCD
      else if (trim(ext) == "xtc") then
        trajtype = TrajTypeXTC
      else if (trim(ext) == "nc" .or. trim(ext) == "rst") then
        trajtype = TrajTypeNCD
      else
        write(iw,'("Get_TrajType> Error.")')
        write(iw,'("Unknown trajectory type: ", a)') trim(ext)
      end if

    end subroutine get_trajtype
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
    subroutine open_trajfile(fname, trajtype, iunit, dcd, xtc, nc, mode)
!-------------------------------------------------------------------------------
      implicit none

      character(*),               intent(in)    :: fname
      integer,                    intent(in)    :: trajtype
      integer,                    intent(inout) :: iunit
      type(s_dcd),                intent(inout) :: dcd
      type(xtcfile),              intent(inout) :: xtc
      type(s_netcdf),   optional, intent(inout) :: nc
      character(len=1), optional, intent(in)    :: mode

      character(len=1) :: rw


      if (trajtype == TrajTypeDCD) then

        call dcd_open(fname, iunit)

      else if (trajtype == TrajTypeXTC) then

        rw = 'r'
        if (present(mode)) &
          rw = mode
          
        call xtc%init(fname, rw)

      else if (trajtype == TrajTypeNCD) then

        call netcdf_open(fname, iunit)

      end if

    end subroutine open_trajfile
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
    subroutine init_trajfile(trajtype, iunit, dcd, xtc, nc, natm, init)
!-------------------------------------------------------------------------------
      implicit none

      integer,                  intent(in)    :: trajtype
      integer,                  intent(in)    :: iunit
      type(s_dcd),              intent(inout) :: dcd
      type(xtcfile),            intent(inout) :: xtc
      type(s_netcdf), optional, intent(inout) :: nc 
      integer,                  intent(out)   :: natm
      logical,        optional, intent(in)    :: init

      logical :: initialize


      initialize = .false.
      if (present(init)) &
        initialize = init

      if (trajtype == TrajTypeDCD) then
        call dcd_read_header(iunit, dcd)
        natm = dcd%natm

        if (.not. allocated(dcd%coord)) then
          call alloc_dcd(natm, 1, dcd) 
        else
          if (initialize) then
            call dealloc_dcd(dcd)
            call alloc_dcd(natm, 1, dcd)
          end if
        end if

      else if (trajtype == TrajTypeXTC) then
        natm = xtc%natoms

      else if (trajtype == TrajTypeNCD) then

        call netcdf_read_dimension(iunit, nc)

      end if

!
    end subroutine init_trajfile 
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
    subroutine read_trajfile_oneframe(trajtype, iunit, istep, dcd, xtc, nc, is_end)
!-------------------------------------------------------------------------------
      implicit none

      integer,        intent(in)    :: trajtype
      integer,        intent(in)    :: iunit
      integer,        intent(in)    :: istep
      type(s_dcd),    intent(inout) :: dcd
      type(xtcfile),  intent(inout) :: xtc
      type(s_netcdf), intent(inout) :: nc 
      logical,        intent(out)   :: is_end


      is_end = .false.

      if (trajtype == TrajTypeDCD) then

        if (istep <= dcd%nstep) then
          call read_dcd_oneframe(iunit, dcd)
        else
          is_end = .true.
        end if

      else if (trajtype == TrajTypeXTC) then
       
        call xtc%read

        if (xtc%STAT /= 0) then
          is_end = .true.
        end if

      else if (trajtype == TrajTypeNCD) then

        if (istep <= nc%nstep) then
          call netcdf_read_oneframe(iunit, istep, nc)
        else
          is_end = .true.
        end if

      end if


    end subroutine read_trajfile_oneframe
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
    subroutine close_trajfile(trajtype, iunit, dcd, xtc, nc)
!-------------------------------------------------------------------------------
      implicit none

      integer,        intent(in)    :: trajtype
      integer,        intent(in)    :: iunit
      type(s_dcd),    intent(inout) :: dcd
      type(xtcfile),  intent(inout) :: xtc
      type(s_netcdf), intent(inout) :: nc 


      if (trajtype == TrajTypeDCD) then
        call dcd_close(iunit)
      else if (trajtype == TrajTypeXTC) then
        call xtc%close
      else if (trajtype == TrajTypeNCD) then
        call netcdf_close(iunit)
      end if

    end subroutine close_trajfile 
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    subroutine send_coord_to_traj(istep, trajtype, dcd, xtc, nc, traj)
!-------------------------------------------------------------------------------
      implicit none

      integer,        intent(in)     :: istep
      integer,        intent(in)     :: trajtype
      type(s_dcd),    intent(in)     :: dcd
      type(xtcfile),  intent(in)     :: xtc
      type(s_netcdf), intent(in)     :: nc 
      type(s_traj),   intent(inout)  :: traj


      if (trajtype == TrajTypeDCD) then

        call get_coord_from_dcd(istep, dcd, traj)

      else if (trajtype == TrajTypeXTC) then

        call get_coord_from_xtc(istep, xtc, traj) 

      else if (trajtype == TrajTypeNCD) then

        call get_coord_from_netcdf(istep, nc, traj)

      end if


    end subroutine send_coord_to_traj
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
    subroutine get_coord_from_dcd(istep, dcd, traj)
!-------------------------------------------------------------------------------
      implicit none

      integer,      intent(in)    :: istep
      type(s_dcd),  intent(in)    :: dcd
      type(s_traj), intent(inout) :: traj

      integer :: iatm, id


      do iatm = 1, traj%natm
        id                           = traj%ind(iatm)
        traj%coord(1:3, iatm, istep) = dcd%coord(1:3, id, istep) 
      end do

      traj%box(1, istep) = dcd%box(1, istep)
      traj%box(2, istep) = dcd%box(2, istep)
      traj%box(3, istep) = dcd%box(3, istep)


    end subroutine get_coord_from_dcd
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
    subroutine get_coord_from_xtc(istep, xtc, traj)
!-------------------------------------------------------------------------------
      implicit none

      integer,       intent(in)    :: istep
      type(xtcfile), intent(in)    :: xtc 
      type(s_traj),  intent(inout) :: traj

      real(8), parameter :: nm2ang = 10.0d0

      integer :: iatm, id


      do iatm = 1, traj%natm
        id                           = traj%ind(iatm)
        traj%coord(1:3, iatm, istep) = xtc%pos(1:3, id) * nm2ang
      end do

      traj%box(1, istep) = xtc%box(1, 1) * nm2ang
      traj%box(2, istep) = xtc%box(2, 2) * nm2ang
      traj%box(3, istep) = xtc%box(3, 3) * nm2ang


    end subroutine get_coord_from_xtc
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
    subroutine get_coord_from_netcdf(istep, nc, traj)
!-------------------------------------------------------------------------------
      implicit none

      integer,        intent(in)    :: istep
      type(s_netcdf), intent(in)    :: nc 
      type(s_traj),   intent(inout) :: traj

      integer :: iatm, id


      do iatm = 1, traj%natm
        id                           = traj%ind(iatm)
        traj%coord(1:3, iatm, istep) = nc%coord(1:3, id, istep) 
      end do

      traj%box(1, istep) = nc%box(1, istep)
      traj%box(2, istep) = nc%box(2, istep)
      traj%box(3, istep) = nc%box(3, istep)


    end subroutine get_coord_from_netcdf
!-------------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !
    !  Subroutine    alloc_traj 
    !  @brief        allocate memories for traj structure 
    !  @authors      KK
    !  @param[in]    natm  : # of atoms 
    !  @param[in]    nstep : # of steps 
    !  @param[inout] traj  : structure of ANATRA trajectory
    !
    !---------------------------------------------------------------------------

    subroutine alloc_traj(natm, nstep, traj)

      implicit none

      integer,      intent(in)    :: natm
      integer,      intent(in)    :: nstep
      type(s_traj), intent(inout) :: traj 


      allocate(traj%atmname(natm))
      allocate(traj%resname(natm))
      allocate(traj%segname(natm))
      allocate(traj%resid(natm))
      allocate(traj%box(1:3, nstep))
      allocate(traj%coord(1:3, natm, nstep))
      allocate(traj%mass(natm))
      allocate(traj%charge(natm))
      allocate(traj%ind(natm))
      allocate(traj%ind_ref(natm))

    end subroutine alloc_traj

    !---------------------------------------------------------------------------
    !
    !  Subroutine    dealloc_traj 
    !  @brief        deallocate memories in traj structure 
    !  @authors      KK
    !  @param[inout] traj : structure of ANATRA trajecotry 
    !
    !---------------------------------------------------------------------------

    subroutine dealloc_traj(traj)

      implicit none

      type(s_traj), intent(inout) :: traj


      deallocate(traj%atmname)
      deallocate(traj%resname)
      deallocate(traj%segname)
      deallocate(traj%resid)
      deallocate(traj%box)
      deallocate(traj%coord)
      deallocate(traj%mass)
      deallocate(traj%charge)
      deallocate(traj%ind)
      deallocate(traj%ind_ref)


    end subroutine dealloc_traj

end module mod_traj
!=======================================================================
