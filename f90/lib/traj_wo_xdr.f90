!-------------------------------------------------------------------------------
!
!  Module   mod_traj  
!
!    Module for defining trajectory parameters and read DCD trajectory  
!
!    (c) Copyright 2024 Osaka Univ. All rights reserved.
!
!-------------------------------------------------------------------------------

module mod_traj

  use mod_util
  use mod_const
  use mod_dcdio

  implicit none

  ! parameters
  !
  integer, parameter :: TrajTypeDCD = 1
  integer, parameter :: TrajTypeXTC = 2

  ! structures
  !
  type :: s_trajopt
    real(8)                       :: dt
    integer                       :: trajtype_in           = TrajTypeDCD
    integer                       :: trajtype_out          = TrajTypeXTC
    character(len=MaxChar)        :: molinfo(MaxTraj)      = ""
    character(len=MaxChar)        :: molinfo_ref(MaxTraj)  = ""

    integer                       :: nmolinfo
    integer                       :: nmolinfo_ref
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

      real(8)                      :: dt                    = 0.0d0
      character(len=MaxChar)       :: molinfo(MaxTraj)      = ""
      character(len=MaxChar)       :: molinfo_ref(MaxTraj)  = ""

      integer :: i, irank
      integer :: nmolinfo, nmolinfo_ref

      namelist /trajopt_param/ dt, molinfo, molinfo_ref

      molinfo      = ""
      molinfo_ref  = ""

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

      nmolinfo_ref = 0
      do i = 1, MaxTraj
        if (trim(molinfo_ref(i)) /= "") nmolinfo_ref = nmolinfo_ref + 1
      end do

      if (irank == 0) then
        do i = 1, nmolinfo_ref
          write(iw,'("molinfo_ref ", i0, " : ", a)') &
            i, trim(molinfo_ref(i))
        end do
      end if

      trajopt%dt           = dt
      trajopt%molinfo      = molinfo 
      trajopt%molinfo_ref  = molinfo_ref

      trajopt%nmolinfo     = nmolinfo
      trajopt%nmolinfo_ref = nmolinfo_ref


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
      character(len=3), optional, intent(in)  :: refinfo

      ! local variables
      integer                :: iatm, icol, ncol, id, ierr
      integer                :: imol
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
        if (refinfo == "ref") then
          open(11,file=trim(trajopt%molinfo_ref(id)))
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
        stop
      end if

      if (trim(ext) == "dcd") then
        trajtype = TrajTypeDCD
      else if (trim(ext) == "xtc") then
        trajtype = TrajTypeXTC
      else
        write(iw,'("Get_TrajType> Error.")')
        write(iw,'("Unknown trajectory type: ", a)') trim(ext)
      end if

    end subroutine get_trajtype
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
