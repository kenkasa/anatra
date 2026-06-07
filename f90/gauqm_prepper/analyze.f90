!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_ctrl
  use mod_traj
  use mod_com

  ! constants
  !
  integer, parameter :: ID_T = 1
  integer, parameter :: ID_C = 2

  ! subroutines
  !
  public :: analyze 

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option, traj)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in)    :: input
      type(s_output), intent(in)    :: output
      type(s_option), intent(in)    :: option
      type(s_traj),   intent(inout) :: traj(2)

      type(s_com)    :: com(2)
      type(s_dcd)    :: dcd
      type(xtcfile)  :: xtc
      type(s_netcdf) :: nc

      ! I/O
      !
      integer                :: io, io_g
      character(len=MaxChar) :: fname

      ! Local
      !
      integer :: natm, natm_u, natm_v, nmol(2), nstep_tot
      real(8) :: box(3), tmp_crd(3), tvec(3) 
      real(8) :: m, d(3)
      integer :: trajtype
      logical :: is_end

      ! Dummy
      !
      integer :: i, itraj, istep, istep_tot
      integer :: iatm, imol, ixyz, ierr

      ! Arrays
      !
      integer, allocatable :: atmnum(:)


      ! Get trajectory type
      !
      call get_trajtype(input%ftraj(1), trajtype)

      ! Setup molecule
      !
      do i = 1, 2
        call get_com(option%mode(i),       &
                     traj(i),              &
                     com(i),               &
                     setup = .true.,       &
                     calc_coord = .false., &
                     myrank = 0)
      end do

      nmol(1) = com(1)%nmol
      nmol(2) = com(2)%nmol

      natm_u  = traj(1)%natm
      natm_v  = traj(2)%natm

      if (nmol(1) /= 1) then
        write(iw,'("Analyze> Error. # of molecules in solute should be 1")')
        stop
      end if
      
      allocate(atmnum(natm_u))

      atmnum = 0
      do iatm = 1, natm_u
        atmnum(iatm) = get_atomicnum(traj(1)%mass(iatm), ierr)
        if (ierr /= 0) then
          write(iw,'("Analyze> Error.")')
          write(iw,'("Failed to assign atomic number of ",a)') &
                  trim(traj(1)%atmname(iatm))
          stop
        end if
      end do

      ! Analyze
      !
      is_end    = .false.
      istep_tot = 0
      do itraj = 1, input%ntraj
        call open_trajfile(input%ftraj(itraj), trajtype, io, dcd, xtc, nc)
        call init_trajfile(trajtype, io, dcd, xtc, nc, natm)

        istep  = 0
        is_end = .false.
        do while (.not. is_end)
          istep     = istep     + 1
          istep_tot = istep_tot + 1

          call read_trajfile_oneframe(trajtype, io, istep, dcd, xtc, nc, is_end)

          if (is_end) exit

          if (mod(istep_tot, 100) == 0) then
            write(iw,'("Step ",i0)') istep_tot
          end if

          do i = 1, 2
            call send_coord_to_traj(1, trajtype, dcd, xtc, nc, traj(i))
            call get_com(option%mode(i),      &
                         traj(i),             &
                         com(i),              &
                         setup = .false.,     &
                         calc_coord = .true., & 
                         myrank = 1)
          end do

          tvec(1:3) = com(1)%coord(1:3, 1, 1)
          box(1:3)  = traj(1)%box(1:3, 1)

          write(iw,'("Coord ", f20.13)') traj(1)%coord(1, 1, 1)
          write(fname,'(a,i4.4,".dat")') trim(output%fhead), istep_tot
          call open_file(fname, io_g)

          ! Output Solute Coordinate
          !
          do iatm = 1, natm_u
            write(io_g, '(i5,3f20.13)') & 
              atmnum(iatm),            &
              (traj(1)%coord(ixyz, iatm, 1) - tvec(ixyz), ixyz = 1, 3) 
          end do
          write(io_g,*)

          ! Output Solvent Coordinate with charge
          !
          do iatm = 1, natm_v
            tmp_crd(1:3) = traj(2)%coord(1:3, iatm, 1) - tvec(1:3)
            tmp_crd(1:3) = tmp_crd(1:3) - box(1:3) * nint(tmp_crd(1:3) / box(1:3))
            write(io_g, '(4f20.13)') (tmp_crd(ixyz), ixyz = 1, 3), traj(2)%charge(iatm)
          end do 
          write(io_g,*)

          ! Output Solvent Coordinate 
          !
          do iatm = 1, natm_v
            tmp_crd(1:3) = traj(2)%coord(1:3, iatm, 1) - tvec(1:3)
            tmp_crd(1:3) = tmp_crd(1:3) - box(1:3) * nint(tmp_crd(1:3) / box(1:3))
            write(io_g, '(3f20.13)') (tmp_crd(ixyz), ixyz = 1, 3)
          end do 
          write(io_g,*)
          close(io_g)

        end do ! step

        istep_tot = istep_tot - 1
        call close_trajfile(trajtype, io, dcd, xtc, nc)

      end do   ! traj

      nstep_tot = istep_tot

      ! Deallocate
      !

    end subroutine analyze
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
