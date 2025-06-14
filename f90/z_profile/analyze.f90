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
      integer                :: io
      character(len=MaxChar) :: fname 

      ! Local
      !
      integer :: natm, nmol, nmol_c, nz, nstep_tot
      real(8) :: boxave(3), vol, dv, zsta, zmin, zmax, zij 
      real(8) :: m, d(3), d2, z, wd
      integer :: trajtype
      logical :: is_end

      ! Dummy
      !
      integer :: i, itraj, istep, istep_tot
      integer :: iatm, imol, imol_c, ixyz, iz, ierr

      ! Arrays
      !
      real(8), allocatable :: weight(:) 
      real(8), allocatable :: gz(:) 


      allocate(weight(traj(1)%natm))

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

      nmol   = com(ID_T)%nmol
      nmol_c = com(ID_C)%nmol

      if (nmol_c /= 1) then

        if (nmol /= nmol_c) then
          write(iw,'("Analyze> Error.")')  
          write(iw,'("# of molecules in target (nmol) should be the same as that in center (nmol_c)")')
          write(iw,'("when nmol_c /= 1 (mode(2) /= WHOLE)")')  
          stop
        end if

      end if 

      ! Get weight
      !
      if (option%denstype == DensityTypeNUMBER) then
        weight = 1.0d0
      else if (option%denstype == DensityTypeELECTRON) then
        do iatm = 1, traj(ID_T)%natm
          weight(iatm) = - traj(ID_T)%charge(iatm) &
                         + get_atomicnum(traj(ID_T)%mass(iatm), ierr)
          if (ierr /= 0) then
            write(iw,'("Analyze> Error.")')
            write(iw,'("Failed to assign atomic number of ",a)') &
                    trim(traj(ID_T)%atmname(iatm))
          end if
        end do
      end if

      ! Get boxsize at first step for determining # of grids
      !
      call open_trajfile(input%ftraj(1), trajtype, io, dcd, xtc, nc)
      call init_trajfile(trajtype, io, dcd, xtc, nc, natm)
      call read_trajfile_oneframe(trajtype, io, 1, dcd, xtc, nc, is_end)
      call send_coord_to_traj(1, trajtype, dcd, xtc, nc, traj(1))

      boxave(:) = traj(1)%box(:, 1)

      zmax      = boxave(3)

      if (option%symmetrize) then
        zsta      = 0.0d0 
        nz        = 0.5d0 * zmax / option%dz + 1
        wd        = 0.5d0
      else
        zsta      = - 0.5d0 * zmax
        nz        = zmax / option%dz + 1
        wd        = 1.0d0
      end if

      call close_trajfile(trajtype, io, dcd, xtc, nc)
      boxave(:) = 0.0d0

      ! Allocate memory
      !
      allocate(gz(0:nz))

      ! Analyze
      !
      gz        = 0.0d0
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

          boxave(:) = boxave(:) + traj(ID_T)%box(:, 1)

          do imol = 1, nmol

            if (option%mode(ID_C) == CoMModeWHOLE) then
              imol_c = 1
            else
              imol_c = imol
            end if

            zij = com(ID_T)%coord(3, imol, 1) - com(ID_C)%coord(3, imol_c, 1)
            zij = zij - traj(ID_T)%box(3, 1) * nint(zij / traj(ID_T)%box(3, 1))
            iz  = nint((zij - zsta) / option%dz) 
            !iz  = nint((com(ID_T)%coord(3, imol, 1) - com(ID_C)%coord(3, 1, 1) - zsta) / option%dz)

            if (option%symmetrize .and. iz < 0) then
              iz = - iz 
            end if

            !if (iz > 0 .and. iz < nz) then 
            !  gz(iz) = gz(iz) + weight(imol)
            !else if (iz == 0 .or. iz == nz) then
            !  gz(iz) = gz(iz) + weight(imol) * 0.5d0
            !end if

            if (iz >= 0 .and. iz <= nz) then 
              gz(iz) = gz(iz) + wd * weight(imol)
            end if

          end do
        end do ! step

        istep_tot = istep_tot - 1
        call close_trajfile(trajtype, io, dcd, xtc, nc)

      end do   ! traj

      nstep_tot = istep_tot

      boxave(:) = boxave(:) / dble(nstep_tot)
      dv        = option%dz * boxave(1) * boxave(2)
      vol       = boxave(1) * boxave(2) * boxave(3)

      gz = gz / (dv * nstep_tot)

      ! Generate files
      !
      write(fname,'(a,".zprof")') trim(output%fhead)
      call open_file(fname, io)
      do iz = 0, nz
        z = zsta + option%dz * iz 
        write(io,'(2f20.10)') z, gz(iz) 
      end do
      close(io)

      ! Deallocate
      !
      deallocate(weight, gz)


    end subroutine analyze
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
