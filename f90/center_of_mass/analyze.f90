!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_input
  use mod_ctrl
  use mod_xtcio
  use mod_dcdio
  use mod_netcdfio
  use mod_traj
  use mod_com
  use xdr, only: xtcfile

  ! subroutines
  !
  public  :: analyze
  private :: calc_msd 

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option, timegrid, traj)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),    intent(in)    :: input 
      type(s_output),   intent(in)    :: output
      type(s_option),   intent(in)    :: option
      type(s_timegrid), intent(in)    :: timegrid
      type(s_traj),     intent(inout) :: traj

      ! I/O
      !
      integer                :: io

      ! Local
      ! 
      type(s_dcd)            :: dcd
      type(xtcfile)          :: xtc
      type(s_netcdf)         :: nc 
      type(s_com)            :: com, com_all

      integer                :: nmol, natm, nstep_tot
      integer                :: trajtype
      logical                :: is_end

      ! Dummy
      !
      integer                :: itraj
      integer                :: istep, istep_tot

      ! Arrays
      !
      real(8), allocatable   :: msd(:, :), msdave(:)
      real(8), allocatable   :: msd_sterr(:)


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

      call get_com(option%mode,          &
                   traj,                 &
                   com,                  &
                   setup = .true.,       &
                   calc_coord = .false., &
                   myrank = 0)

      call get_com(option%mode,          &
                   traj,                 &
                   com_all,              &
                   setup = .true.,       &
                   calc_coord = .false., &
                   myrank = 1)

      nmol          = com_all%nmol
      com_all%nstep = nstep_tot

      allocate(com_all%coord(1:3, nmol, nstep_tot))

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

          call send_coord_to_traj(1, trajtype, dcd, xtc, nc, traj)

          call get_com(option%mode,         &
                       traj,                &
                       com,                 &
                       setup = .false.,     &
                       calc_coord = .true., & 
                       myrank = 1)

          com_all%coord(1:3, 1:nmol, istep_tot) &
            = com%coord(1:3, 1:nmol, 1)

        end do

        istep_tot = istep_tot - 1

        call close_trajfile(trajtype, io, dcd, xtc, nc)

      end do

      ! Calculate MSD
      !
      if (option%out_msd) then
        allocate(msd(0:option%nt_range, nmol))
        allocate(msdave(0:option%nt_range))
        allocate(msd_sterr(0:option%nt_range))

        write(iw,*)
        write(iw,'("Analyze> Calculate MSD")')

        call calc_msd(com_all, option, timegrid, msd, msdave, msd_sterr)

        write(iw,'("Finished")')
        write(iw,*)
      end if

      ! Output
      !
      call generate_comfile(output,   &
                            option,   &
                            traj,     &
                            com_all)

      call generate_msdfile(output,   &
                            option,   &
                            timegrid, &
                            traj,     &
                            com_all,  &
                            msd,      &
                            msdave,   &
                            msd_sterr)

      ! Deallocate 
      !
      if (allocated(msd)) &
        deallocate(msd, msdave, msd_sterr)


    end subroutine analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine calc_msd(com, option, timegrid, msd, msdave, msd_sterr)
!-----------------------------------------------------------------------
      implicit none

      type(s_com),      intent(in)    :: com
      type(s_option),   intent(in)    :: option
      type(s_timegrid), intent(in)    :: timegrid 
      real(8),          intent(inout) :: msd(0:option%nt_range, com%nmol)
      real(8),          intent(inout) :: msdave(0:option%nt_range)
      real(8),          intent(inout) :: msd_sterr(0:option%nt_range)

      integer :: imol, istep, jstep, it, ij, nend
      integer :: nmol, nstep, nt_range, msddim 
      real(8) :: d(3), d2, d4, dev

      real(8), allocatable :: crd(:, :, :)
      integer, allocatable :: msd_count(:, :)


      ! setup variables
      !
      nmol     = com%nmol
      nstep    = com%nstep
      nt_range = option%nt_range
      msddim   = option%msddim

      ! allocate memory 
      !
      allocate(crd(3, nstep, nmol))
      allocate(msd_count(0:nt_range, nmol))

      ! prepare modified array of CoM 
      !
      do istep = 1, nstep
        do imol = 1, nmol
          crd(1:3, istep, imol) = com%coord(1:3, imol, istep)
        end do
      end do

      ! calculate MSD
      !
      msd = 0.0d0

!$omp parallel private(imol, istep, jstep, it, d, d2), &
!$omp        & default(shared)
!$omp do

      ! Calculate MSD for each
      !
      do imol = 1, nmol

        do istep = 1, nstep - 1

          do it = 0, timegrid%ng
            jstep = istep + timegrid%ind(it)
            if (jstep > nstep) exit

            d(1:3)        = crd(1:3, jstep, imol) - crd(1:3, istep, imol)
            d2            = dot_product(d(1:msddim), d(1:msddim))

            msd(it, imol)       = msd(it, imol)       + d2
            msd_count(it, imol) = msd_count(it, imol) + 1 

          end do

        end do

        do istep = 0, timegrid%ng 
          msd(istep, imol) = msd(istep, imol) / msd_count(istep, imol) 
        end do

      end do
!$omp end do
!$omp end parallel

      ! Average over molecules
      !
      do istep = 0, timegrid%ng 
        msdave(istep) = sum(msd(istep, 1:nmol)) / dble(nmol)
      end do

      do istep = 0, timegrid%ng
        dev = 0.0d0
        do imol = 1, nmol
          dev = dev + (msd(istep, imol) - msdave(istep)) ** 2
        end do
        dev              = sqrt(dev / dble(nmol - 1))
        msd_sterr(istep) = dev / sqrt(dble(nmol))
      end do

      ! Deallocate
      !
      deallocate(crd)
      deallocate(msd_count)

    end subroutine calc_msd
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine generate_comfile(output, option, traj, com)
!-----------------------------------------------------------------------
      implicit none

      type(s_output), intent(in) :: output
      type(s_option), intent(in) :: option
      type(s_traj),   intent(in) :: traj
      type(s_com),    intent(in) :: com

      ! I/O
      !
      integer                :: io
      character(len=MaxChar) :: fname

      ! Local
      !
      integer :: nmol 

      ! Dummy
      ! 
      integer :: istep, imol, ixyz


      if (.not. option%out_com) &
        return

      nmol = com%nmol

      write(fname,'(a,".com")') trim(output%fhead)

      call open_file(fname, io)

      if (option%onlyz) then

        do istep = 1, com%nstep
          write(io,'(f20.10)',advance='no') traj%dt * istep
          do imol = 1, nmol
            write(io,'(f20.10)',advance='no') com%coord(3, imol, istep)
          end do
          write(io,*)
        end do

      else
        if (option%comformat == CoMFormatTypeXYZ) then

          do istep = 1, com%nstep
            write(io,'(i0)') nmol 
            write(io,*)
            do imol = 1, nmol 
              write(io,'("Ar ",3f20.10)') &
                (com%coord(ixyz, imol, istep), ixyz = 1, 3)
            end do
          end do

        else if (option%comformat == CoMFormatTypeTIMESERIES) then

          do istep = 1, com%nstep
            write(io,'(f20.10)',advance='no') traj%dt * istep
            do imol = 1, nmol
              write(io,'(3f20.10)',advance='no') & 
                (com%coord(ixyz, imol, istep), ixyz = 1, 3)
            end do
            write(io,*)
          end do

        end if 
      end if

      close(io)

    end subroutine generate_comfile
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine generate_msdfile(output,    &
                                option,    &
                                timegrid,  &
                                traj,      &
                                com,       &
                                msd,       &
                                msdave,    &
                                msd_sterr)
!-----------------------------------------------------------------------
      implicit none

      type(s_output),   intent(in) :: output
      type(s_option),   intent(in) :: option
      type(s_timegrid), intent(in) :: timegrid 
      type(s_traj),     intent(in) :: traj
      type(s_com),      intent(in) :: com 
      real(8),          intent(in) :: msd(0:option%nt_range, com%nmol)
      real(8),          intent(in) :: msdave(0:option%nt_range)
      real(8),          intent(in) :: msd_sterr(0:option%nt_range)

      ! I/O
      !
      integer                :: io
      character(len=MaxChar) :: fname

      ! Local
      !
      integer :: nmol

      ! Dummy
      !
      integer :: istep, imol, ixyz


      nmol = com%nmol

      if (.not. option%out_msd) &
        return

      ! MSD (each)
      !
      write(fname,'(a,".msd")') trim(output%fhead)

      call open_file(fname, io)
      do istep = 0, timegrid%ng
        write(io,'(f20.10)', advance='no') timegrid%val(istep) 
        do imol = 1, nmol
          write(io,'(f20.10)', advance='no') msd(istep, imol)
        end do
        write(io,*) 
      end do
      close(io)

      ! MSD (average)
      !
      write(fname,'(a,".msdave")') trim(output%fhead)

      call open_file(fname, io)
      do istep = 0, timegrid%ng 
        write(io,'(3e20.10)') timegrid%val(istep),   &
                              msdave(istep),         &
                              msd_sterr(istep)
      end do
      close(io)


    end subroutine generate_msdfile
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
