!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_ctrl
  use mod_traj
  use mod_com

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

      type(s_dcd)            :: dcd
      type(xtcfile)          :: xtc
      type(s_netcdf)         :: nc 
      type(s_com)            :: com(2)

      ! I/O
      !
      integer                :: io
      character(len=MaxChar) :: fname

      ! Local
      ! 
      integer                :: nend
      integer                :: nmol, natm, nstep, nstep_tot
      integer                :: trajtype
      real(8)                :: cost, vec(3), veclength
      logical                :: is_end

      ! Dummy
      !
      integer                :: i, ixyz, itraj
      integer                :: istep, istep_tot, jstep, ij, imol

      ! Arrays
      !
      real(8), allocatable :: orient_vec(:, :, :)
      real(8), allocatable :: tcf(:, :), tcfave(:)


      ! Get trajectory type
      !
      call get_trajtype(input%ftraj(1), trajtype)

      ! Get Total number of steps
      !
      if (trajtype == TrajTypeDCD) then
        call get_total_step_from_dcd(input%ftraj, nstep_tot)
      else if (trajtype == TrajTypeXTC) then
        call get_total_step_from_xtc(input%ftraj, nstep_tot)
      else if (trajtype == TrajTypeNCD) then
        call get_total_step_from_netcdf(input%ftraj, nstep_tot)
      end if 

      do i = 1, 2
        call get_com(option%mode,          &
                     traj(i),              &
                     com(i),               &
                     setup = .true.,       &
                     calc_coord = .false., &
                     myrank = 0)
      end do

      nmol         = com(1)%nmol
      com(1)%nstep = nstep_tot
      nstep        = nstep_tot 

      if (com(1)%nmol /= com(2)%nmol) then
        write(iw,'("Analyze> Error.")')
        write(iw,'("Number of molecules in two dcd should be the same.")')
        stop
      end if

      ! Allocate memory
      !
      allocate(orient_vec(3, nstep, nmol))
      allocate(tcf(0:option%tcfrange, nmol))
      allocate(tcfave(0:option%tcfrange))

      orient_vec = 0.0d0
      tcf        = 0.0d0
      tcfave     = 0.0d0

      ! Get CoM
      !
      write(iw,*)
      write(iw,'("Analyze> Get CoM coordinates")')
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
            write(iw,'("Step ", i0)') istep_tot
          end if

          do i = 1, 2
            call send_coord_to_traj(1, trajtype, dcd, xtc, nc, traj(i))
            call get_com(option%mode,         &
                         traj(i),             &
                         com(i),              &
                         setup = .false.,     &
                         calc_coord = .true., & 
                         myrank = 1)
          end do

          do imol = 1, nmol
            vec(1:3)  = com(2)%coord(1:3, imol, 1) &
                      - com(1)%coord(1:3, imol, 1)
            veclength = sqrt(dot_product(vec(1:3), vec(1:3)))
            orient_vec(1:3, istep, imol) = vec(1:3) / veclength
          end do 
        end do
        istep_tot = istep_tot - 1
        call close_trajfile(trajtype, io, dcd, xtc, nc)
      end do

      ! calculate tcf
      !
      tcf    = 0.0d0
      tcfave = 0.0d0

!$omp parallel private(imol, istep, jstep, nend, ij, cost) &
!$omp        & default(shared)
!$omp do
      do imol = 1, nmol
        do istep = 1, nstep - 1

          nend = istep + option%tcfrange + 1
          if (nend > nstep) &
            nend = nstep

          do jstep = istep, nend
            ij            = jstep - istep
            cost          = dot_product(orient_vec(1:3, jstep, imol), &
                                        orient_vec(1:3, istep, imol))
            tcf(ij, imol) = tcf(ij, imol) + cost * cost 
                            
          end do
        end do 
      end do
!$omp end do
!$omp end parallel

      do imol = 1, nmol
        tcf(0, imol) = 1.0d0
        do istep = 1, option%tcfrange
          tcf(istep, imol) = tcf(istep, imol) / (nstep - istep)
          tcf(istep, imol) = 1.5d0 * tcf(istep, imol) - 0.5d0
        end do
      end do

      ! calculate averaged tcf
      !
      do istep = 0, option%tcfrange - 1
        tcfave(istep) = sum(tcf(istep, 1:nmol)) / dble(nmol)
      end do

      ! generate files
      !
      write(fname,'(a,".tcf")') trim(output%fhead)
      call open_file(fname, io)
      open(io, file=trim(fname))
      do istep = 0, option%tcfrange - 1
        write(io, '(f20.10)', advance='no') traj(1)%dt * istep
        do imol = 1, nmol
          write(io,'(f20.10)', advance='no') tcf(istep, imol)
        end do
        write(io, *)
      end do
      close(io)

      write(fname,'(a,".tcfave")') trim(output%fhead)
      call open_file(fname, io)
      open(io, file=trim(fname))
      do istep = 0, option%tcfrange - 1
        write(io, '(2f20.10)') traj(1)%dt * istep, tcfave(istep)
      end do
      close(io)

      ! Deallocate memory
      !
      deallocate(orient_vec)
      deallocate(tcf, tcfave)

    end subroutine analyze
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
