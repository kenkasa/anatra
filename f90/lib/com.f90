!===============================================================================
! Module: Mod_CoM
!
!   Module for defining molecules and calculating Center-of-Mass coordinates
!   according to Molinfo and Trajectory structure 
!
!   (c) Copyright 2024 Osaka Univ. All rights reserved.
!
!===============================================================================
module mod_com

  use mod_util
  use mod_const
  use mod_traj

  implicit none

  ! Constants
  !
  integer,      parameter, public :: CoMModeRESIDUE = 1
  integer,      parameter, public :: CoMModeWHOLE   = 2
  integer,      parameter, public :: ComModeATOM    = 3 
  integer,      parameter, public :: ComModeSEGMENT = 4 
  character(*), parameter, public :: CoMMode(4)     = (/'RESIDUE   ',&
                                                        'WHOLE     ',&
                                                        'ATOM      ',&
                                                        'SEGMENT   '/)

  ! Structures
  !
  type :: s_com
    integer              :: nstep = 0
    integer              :: nmol  = 0
    integer, allocatable :: molid(:)
    integer, allocatable :: molsta(:) 
    integer, allocatable :: molend(:)
    real(8), allocatable :: coord(:,:,:) ! (3, imol, istep)
    real(8), allocatable :: mtot(:)      ! (imol) 
  end type s_com

  ! Subroutines
  !
  public :: get_com 

  contains

!-------------------------------------------------------------------------------
!   Subroutine Get_CoM
!
!     Define molecules and calculate CoM
!
!     - INTEGER mode [IN]      
!       Mode for how molecules are defined
!       (see ComMode parameters)
!
!     - s_traj  traj [IN]
!       Trajectory info. (see traj.f90) 
!
!     - s_com   com  [IN] 
!       CoM info. (defined in this module file)
!
!     - LOGICAL unwrap [IN, OPTIONAL]
!       Unwrap procedure is performed or not
!       If this variable is not specified, unwrap is NOT performed
!
!     - LOGICAL setup [IN, OPTIONAL]
!       Defining molecule procedure is performed or not
!       If this variable is not specified, this procedure is performed
!        
!     - LOGICAL calc_coord [IN, OPTIONAL]
!       Calculating CoM coordinate procedure is performed or not
!       If this variable is not specified, this procedure is performed
!
!     - INTEGER myrank [IN, OPTIONAL]
!       MPI process id or OpenMP thread id
!       if myrank = 0, several messeges will be outputted 
!
!-------------------------------------------------------------------------------
    subroutine get_com(mode, traj, com, unwrap, setup, calc_coord, myrank)

      implicit none

      integer,           intent(in)    :: mode
      type(s_traj),      intent(in)    :: traj
      type(s_com),       intent(inout) :: com
      logical, optional, intent(in)    :: unwrap
      logical, optional, intent(in)    :: setup
      logical, optional, intent(in)    :: calc_coord
      integer, optional, intent(in)    :: myrank

      integer              :: istep, ixyz, ires, ires_prev, id
      integer              :: iatm, imol, jmol, nmol, sgn
      integer              :: irank
      real(8)              :: m
      real(8)              :: d(3), boxh(3)
      character(len=4)     :: seg_prev, seg
      logical              :: su, cc

      integer, allocatable :: njump(:, :)
      real(8), allocatable :: crd_prev(:, :), crd_now(:, :)


      ! Setup Procedures 
      !
      irank = 0
      if (present(myrank)) &
        irank = myrank

      su = .true.
      if (present(setup)) &
        su = setup

      cc = .true.
      if (present(calc_coord)) &
        cc = calc_coord

      ! Setup Molecule
      !
      if (su) then
        allocate(com%molid(traj%natm))
       
        ! Get number of molecules
        !
        if (mode == ComModeRESIDUE) then

          ires_prev = traj%resid(1)
          imol = 1

          do iatm = 1, traj%natm
            ires = traj%resid(iatm)
            if (ires == ires_prev) then
              com%molid(iatm) = imol 
            else
              imol             = imol + 1
              com%molid(iatm)  = imol 
              ires_prev        = ires 
            end if
          end do
          nmol = imol

        else if (mode == ComModeWHOLE) then

          com%molid = 1
          nmol      = 1

        else if (mode == ComModeATOM) then

          do iatm = 1, traj%natm
            com%molid(iatm) = iatm
          end do
          nmol = traj%natm

        else if (mode == ComModeSEGMENT) then

          seg_prev = traj%segname(1)
          imol     = 1

          do iatm = 1, traj%natm
            seg = traj%segname(iatm)
            if (trim(seg) == trim(seg_prev)) then
              com%molid(iatm) = imol 
            else
              imol             = imol + 1
              com%molid(iatm)  = imol 
              seg_prev         = seg 
            end if
          end do
          nmol = imol

        end if
      
        com%nmol  = nmol

        ! Get first and final atom indice for each molecule
        !
        allocate(com%molsta(nmol), com%molend(nmol))
       
        jmol = 0 
        do iatm = 1, traj%natm
          imol = com%molid(iatm)
          if (imol /= jmol) then
            com%molsta(imol) = iatm
            if (jmol /= 0) then
              !com%molend(jmol) = iatm
              com%molend(jmol) = iatm - 1
            end if
          end if
          jmol = imol
        end do

        com%molend(nmol) = traj%natm
       
        ! Get total mass for each molecule
        !
        allocate(com%mtot(nmol))
       
        com%mtot = 0.0d0
        do iatm = 1, traj%natm
          m              = traj%mass(iatm)
          imol           = com%molid(iatm)
          com%mtot(imol) = com%mtot(imol) + m 
        end do
        com%nstep = 1 ! changed at later step
      end if

      ! Get CoM
      !
      if (.not. cc) &
        return

      nmol      = com%nmol
      com%nstep = traj%nstep 

      if (.not. allocated(com%coord)) &
        allocate(com%coord(3, nmol, traj%nstep))

      com%coord = 0.0d0
      do istep = 1, traj%nstep
        do iatm = 1, traj%natm
          m    = traj%mass(iatm)
          imol = com%molid(iatm)
          com%coord(:, imol, istep) = &
            com%coord(:, imol, istep) + m * traj%coord(:, iatm, istep)
        end do
      end do


      do istep = 1, traj%nstep
        do imol = 1, nmol
          do ixyz = 1, 3
           com%coord(ixyz, imol, istep) = &
             com%coord(ixyz, imol, istep) / com%mtot(imol)
          end do
        end do
      end do

      if (irank == 0) then
        write(iw,*)
        write(iw,'("Get_CoM>")')
        write(iw,'(" # of molecules = ", i0)') nmol
      end if

      ! Perform unwrap
      !
      if (present(unwrap)) then

        if (unwrap) then

          write(iw,*)
          write(iw,'(" unwrap         = .true.")')
          write(iw,'(" >> perform unwrap")')
          write(iw,'(" Note that reliability of unwrap depends on ")')
          write(iw,'(" the time interval between snapshots.")')
          write(iw,'(" Please carefully check the results")')
          write(iw,*)

          ! allocate memory for unwrap process
          !
          allocate(njump(1:3, nmol))
          allocate(crd_prev(1:3, nmol), crd_now(1:3,nmol))

          crd_prev(1:3, 1:nmol) = com%coord(1:3, 1:nmol, 1)
          njump                 = 0

          do istep = 2, traj%nstep

            do imol = 1, nmol
              do ixyz = 1, 3
                com%coord(ixyz, imol, istep) =   &
                    com%coord(ixyz, imol, istep) &
                  + njump(ixyz, imol) * traj%box(ixyz, istep) 
              end do
            end do

            boxh(1:3)            = traj%box(1:3, istep) * 0.5d0
            crd_now(1:3, 1:nmol) = com%coord(1:3, 1:nmol, istep)

            do imol = 1, nmol
              d(1:3) = crd_now(1:3, imol) - crd_prev(1:3, imol)

              do ixyz = 1, 3
                if (abs(d(ixyz)) > boxh(ixyz)) then
                  sgn                 = nint(sign(1.0d0, d(ixyz)))
                  njump(ixyz, imol)   = njump(ixyz, imol) - sgn 
                  crd_now(ixyz, imol) = &
                    crd_now(ixyz, imol) - sgn * traj%box(ixyz, istep)

                  com%coord(ixyz, imol, istep) = crd_now(ixyz, imol)  
                end if
              end do

            end do

            crd_prev(1:3, 1:nmol) = crd_now(1:3, 1:nmol)

          end do

          ! Deallocate memory
          !
          deallocate(njump)
          deallocate(crd_prev, crd_now)

        end if 
      end if

    end subroutine get_com

!-------------------------------------------------------------------------------
!   Subroutine Get_CoM_Again
!
!     Calculate CoM. Unlike subroutine Get_CoM, defining molecules procedure is
!     not performed
!
!     - s_traj  traj   [IN]
!       Trajectory info. (see traj.f90) 
!
!     - s_com   com    [IN] 
!       CoM info. (defined in this module file)
!
!     - LOGICAL unwrap [IN, OPTIONAL]
!       Unwrap procedure is performed or not
!       If this variable is not specified, unwrap is NOT performed
!
!-------------------------------------------------------------------------------
    subroutine get_com_again(traj, com, unwrap)

      implicit none

      type(s_traj),      intent(in)    :: traj
      type(s_com),       intent(inout) :: com
      logical, optional, intent(in)    :: unwrap

      integer              :: istep, ixyz, ires, ires_prev, id
      integer              :: iatm, imol, jmol, nmol, sgn
      real(8)              :: m
      real(8)              :: d(3), boxh(3)
      character(len=4)     :: seg_prev, seg

      integer, allocatable :: njump(:, :)
      real(8), allocatable :: crd_prev(:, :), crd_now(:, :)


      ! get CoM
      !
      nmol      = com%nmol
      com%coord = 0.0d0
      do istep = 1, traj%nstep
        do iatm = 1, traj%natm
          m    = traj%mass(iatm)
          imol = com%molid(iatm)
          com%coord(:, imol, istep) = &
            com%coord(:, imol, istep) + m * traj%coord(:, iatm, istep)
        end do
      end do

      do istep = 1, traj%nstep
        do imol = 1, nmol
          do ixyz = 1, 3
           com%coord(ixyz, imol, istep) = &
             com%coord(ixyz, imol, istep) / com%mtot(imol)
          end do
        end do
      end do

      if (present(unwrap)) then

        if (unwrap) then

          ! allocate memory for unwrap process
          !
          allocate(njump(1:3, nmol))
          allocate(crd_prev(1:3, nmol), crd_now(1:3,nmol))

          ! perform unwrap
          !
          crd_prev(1:3, 1:nmol) = com%coord(1:3, 1:nmol, 1)
          njump                 = 0

          do istep = 2, traj%nstep

            do imol = 1, com%nmol
              do ixyz = 1, 3
                com%coord(ixyz, imol, istep) =   &
                    com%coord(ixyz, imol, istep) &
                  + njump(ixyz, imol) * traj%box(ixyz, istep) 
              end do
            end do

            boxh(1:3)            = traj%box(1:3, istep) * 0.5d0
            crd_now(1:3, 1:nmol) = com%coord(1:3, 1:nmol, istep)

            do imol = 1, nmol
              d(1:3) = crd_now(1:3, imol) - crd_prev(1:3, imol)

              do ixyz = 1, 3
                if (abs(d(ixyz)) > boxh(ixyz)) then
                  sgn                 = nint(sign(1.0d0, d(ixyz)))
                  njump(ixyz, imol)   = njump(ixyz, imol) - sgn 
                  crd_now(ixyz, imol) = &
                    crd_now(ixyz, imol) - sgn * traj%box(ixyz, istep)

                  com%coord(ixyz, imol, istep) = crd_now(ixyz, imol)  
                end if
              end do

            end do

            crd_prev(1:3, 1:nmol) = crd_now(1:3, 1:nmol)

          end do

          ! deallocate memory
          !
          deallocate(njump)
          deallocate(crd_prev, crd_now)

        end if 
      end if

    end subroutine get_com_again 

end module
!=======================================================================
