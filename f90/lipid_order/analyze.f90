!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_ctrl
  use mod_traj
  use mod_netcdfio
  use mod_com

  ! subroutines
  !
  public :: analyze 

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option, traj)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in)     :: input
      type(s_output), intent(in)     :: output
      type(s_option), intent(in)     :: option
      type(s_traj),   intent(inout)  :: traj(:)

      type(s_dcd)              :: dcd
      type(xtcfile)            :: xtc
      type(s_netcdf)           :: nc

      ! I/O
      !
      integer                :: io, io_t
      character(len=MaxChar) :: fname

      ! Local
      !
      integer :: natm, nstep_tot, nmolinfo
      integer :: trajtype
      real(8) :: c0(3), c1(3), d(3), d2
      real(8) :: cossq, dev
      logical :: is_end

      ! Dummy
      !
      integer :: i, j
      integer :: istep, istep_tot, itraj
      integer :: iatm, jatm, imol, id

      ! Arrays
      !
      type(s_com), allocatable :: com(:)
      integer,     allocatable :: nmol(:), nsite(:)
      real(8),     allocatable :: scd(:,:), scdave(:), sterr(:)


      write(iw,*)
      write(iw,'("Analyze> Start the analysis")')

      ! Get trajectory types
      !
      call get_trajtype(input%ftraj(1), trajtype)


      ! Allocate memory
      !
      nmolinfo = size(traj(:))
      allocate(nmol(nmolinfo), nsite(nmolinfo))
      allocate(com(nmolinfo))

      ! Setup molecule
      !
      nsite = 0
      do i = 1, nmolinfo 
        call get_com(ComModeRESIDUE,       &
                     traj(i),              &
                     com(i),               &
                     setup = .true.,       &
                     calc_coord = .false., &
                     myrank = 0)

        nmol(i) = com(i)%nmol

        do j = 1, traj(i)%natm
          id = com(i)%molid(j)
          if (id == 1) then 
            nsite(i) = nsite(i) + 1
          end if
        end do
      end do

      ! Allocate memory
      !
      allocate(scd(MaxTraj, nmolinfo))
      allocate(scdave(nmolinfo), sterr(nmolinfo))

      scd       = 0.0d0
      istep_tot = 0
      do itraj = 1, input%ntraj
        call open_trajfile(input%ftraj(itraj), trajtype, io, dcd, xtc, nc)
        call init_trajfile(trajtype, io, dcd, xtc, nc, natm)

        is_end = .false.
        istep  = 0

        do while (.not. is_end)
          istep     = istep     + 1
          istep_tot = istep_tot + 1

          call read_trajfile_oneframe(trajtype, io, istep, dcd, xtc, nc, is_end)

          if (is_end) exit

          if (mod(istep_tot, 100) == 0) then
            write(iw,'("Step ",i0)') istep_tot
          end if

          do i = 1, nmolinfo
            call send_coord_to_traj(1, trajtype, dcd, xtc, nc, traj(i))
            call get_com(ComModeRESIDUE,       &
                         traj(i),              &
                         com(i),               &
                         setup = .false.,      &
                         calc_coord = .true.,  &
                         myrank = 1)
          end do

          do i = 1, nmolinfo
          
            jatm = 0
            do imol = 1, nmol(i)
              jatm  = jatm + 1
              c0(:) = traj(i)%coord(:, jatm, 1) 

              cossq = 0.0d0
              do iatm = 2, nsite(i)
                jatm  = jatm + 1
                c1(:) = traj(i)%coord(:, jatm, 1)
                d(:)  = c1(:) - c0(:)
                d2    = dot_product(d, d)

                cossq = cossq + d(3) * d(3) / d2
              end do

              scd(imol, i) = scd(imol, i) + 1.5d0 * cossq / dble(nsite(i) - 1) - 0.5d0 

            end do  ! imol

          end do    ! molinfo

        end do      ! step 

        istep_tot = istep_tot - 1

        call close_trajfile(trajtype, io, dcd, xtc, nc)

      end do

      nstep_tot = istep_tot

      ! averaging for each molecule
      !
      do i = 1, nmolinfo
        do imol = 1, nmol(i)
          scd(imol, i) = scd(imol, i) / dble(nstep_tot)
        end do
      end do

      do i = 1, nmolinfo
        scdave(i) = sum(scd(1:nmol(i), i)) / dble(nmol(i))
      end do

      sterr = 0.0d0
      do i = 1, nmolinfo

        dev = 0.0d0

        do imol = 1, nmol(i)
          dev      = dev + (scd(imol, i) - scdave(i))**2
        end do

        dev      = sqrt(dev / dble(nmol(i) - 1))
        sterr(i) = dev / sqrt(dble(nmol(i))) 

      end do

      ! Output
      !
      write(fname, '(a,".scd")') trim(output%fhead)

      call open_file(fname, io)
        do i = 1, nmolinfo
          write(io,'(i10, 2f15.7)') option%carbon_id(i), abs(scdave(i)), sterr(i) 
        end do
      close(io)


      ! Deallocate memory
      !
      deallocate(nmol, nsite)
      deallocate(scd, scdave, sterr)

    end subroutine analyze
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
