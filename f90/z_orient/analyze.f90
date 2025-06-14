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
  integer, parameter :: ID_B   = 1 ! Bottom 
  integer, parameter :: ID_H   = 2 ! Head 
  integer, parameter :: ID_C   = 3 ! Center
  integer, parameter :: ID_V0B = 4
  integer, parameter :: ID_V0H = 5 
  integer, parameter :: ID_V1B = 6 
  integer, parameter :: ID_V1H = 7 

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
      type(s_traj),   intent(inout) :: traj(7)

      type(s_com)    :: com(7)
      type(s_dcd)    :: dcd
      type(xtcfile)  :: xtc
      type(s_netcdf) :: nc 

      ! I/O
      !
      integer                :: io, io_t
      character(len=MaxChar) :: fname

      ! Local
      !
      integer :: natm, ngrid, nmol(7), nstep_tot, nsel
      real(8) :: sumd, zave, dlen, prod, ez(3), vec0(3), vec1(3)
      real(8) :: dg, xsta, cossq, sorder_ave
      integer :: trajtype
      logical :: is_end
      
      ! Dummy
      !
      integer                :: itraj
      integer                :: i, ig, igc, imol
      integer                :: istep, istep_tot

      ! Arrays
      !
      real(8), allocatable :: dp(:, :), angle(:), cost(:), theta(:)
      real(8), allocatable :: sorder(:)
      real(8), allocatable :: distr(:)


      if (option%orient_type == OrientTypeCOSTHETA) then
        dg    = option%dx
        ngrid = 2.0d0 / dble(option%dx)
        xsta  = -1.0d0
      else if (option%orient_type == OrientTypeTHETA) then
        dg    = option%dx
        ngrid = 180.0d0 / dble(option%dx)
        xsta  = 0.0d0
      end if

      ! Get trajectory types
      !
      call get_trajtype(input%ftraj(1), trajtype)

      ! Setup molecule
      !
      nsel = option%nsel
      do i = 1, nsel 
        call get_com(option%mode(i),       &
                     traj(i),              &
                     com(i),               &
                     setup = .true.,       &
                     calc_coord = .false., &
                     myrank = 0)

        nmol(i) = com(i)%nmol
      end do

      if (nmol(ID_B) /= nmol(ID_H)) then
        write(iw,'("Analyze> Error.")')
        write(iw,'("Number of molecules in two selections should be the same.")')
        stop
      end if

      ! allocate memory
      !

      allocate(dp(3, nmol(ID_B)), angle(nmol(ID_B)), cost(nmol(ID_B)))
      allocate(theta(nmol(ID_B)), sorder(nmol(ID_B)))
      allocate(distr(0:ngrid))

      ! calculate histogram
      !
      write(iw,*)
      write(iw,'("Analyze> Start orientation analysis")')
      !
      distr      = 0.0d0
      sorder_ave = 0.0d0
      istep_tot  = 0

      if (option%orient_type == OrientTypeCOSTHETA) then
        write(fname,'(a,".costheta")') trim(output%fhead)
      else if  (option%orient_type == OrientTypeTHETA) then
        write(fname,'(a,".theta")')    trim(output%fhead)
      end if

      call open_file(fname, io)

      do itraj = 1, input%ntraj
        call open_trajfile(input%ftraj(itraj), trajtype, io_t, dcd, xtc, nc)
        call init_trajfile(trajtype, io_t, dcd, xtc, nc, natm)

        is_end = .false.
        istep  = 0

        do while (.not. is_end)
          istep     = istep     + 1
          istep_tot = istep_tot + 1

          call read_trajfile_oneframe(trajtype, io_t, istep, dcd, xtc, nc, is_end)

          if (is_end) exit

          if (mod(istep_tot, 100) == 0) then
            write(iw,'("Step ",i0)') istep_tot
          end if

          do i = 1, nsel 
            call send_coord_to_traj(1, trajtype, dcd, xtc, nc, traj(i))
            call get_com(option%mode(i),       &
                         traj(i),              &
                         com(i),               &
                         setup = .false.,      &
                         calc_coord = .true.,  &
                         myrank = 1)
          end do

          if (option%zdef_type == ZdefTypeOUTERPROD) then
            vec0(:) = com(ID_V0H)%coord(:, 1, 1) - com(ID_V0B)%coord(:, 1, 1)
            vec1(:) = com(ID_V1H)%coord(:, 1, 1) - com(ID_V1B)%coord(:, 1, 1)

            ez(1) = vec0(2) * vec1(3) - vec0(3) * vec1(2) 
            ez(2) = vec0(3) * vec1(1) - vec0(1) * vec1(3) 
            ez(3) = vec0(1) * vec1(2) - vec0(2) * vec1(1)
            ez(:) = ez(:) / sqrt(dot_product(ez(1:3), ez(1:3))) 
          else
            ez(:) = 0.0d0
            ez(3) = 1.0d0
          end if

!$omp parallel private(imol, dlen, zave, cossq), &
!$omp        & shared(dp, cost, theta, angle),   &
!$omp        & default(shared),                  &
!$omp        & reduction(+:sorder_ave)
!$omp do
!
          do imol = 1, nmol(ID_B)
            dp(:, imol) = com(ID_H)%coord(:, imol, 1) &
                        - com(ID_B)%coord(:, imol, 1)

            dlen        = sqrt(dot_product(dp(1:3, imol), dp(1:3, imol)))
            dp(:, imol) = dp(:, imol) / dlen

            prod        = dot_product(dp(:, imol), ez(:))

            ! cos(theta) & theta
            !
            if (option%judgeup == JudgeUpModeNMOLUP) then

              if (imol <= option%nmolup) then
                cost(imol) =  prod ! dp(3, imol) / dlen 
              else
                cost(imol) = -prod ! dp(3, imol) / dlen
              end if
            else if (option%judgeup == JudgeUpModeCOORD) then
              zave = 0.5d0 * (com(1)%coord(3, imol, 1) + com(2)%coord(3, imol, 1))
              zave = zave - com(ID_C)%coord(3, 1, 1) 

              if (zave >= 0.0d0) then
                cost(imol) =  prod ! dp(3, imol) / dlen
              else
                cost(imol) = -prod ! dp(3, imol) / dlen
              end if
            else
              cost(imol) =  prod   !dp(3, imol) / dlen 
            endif

            theta(imol) = acos(cost(imol)) / PI * 180.0d0

            ! order parameter (3/2) * cos^2 - 1/2  
            !
            cossq = cost(imol) ** 2
            sorder(imol) = 1.5d0 * cossq - 0.5d0

            sorder_ave = sorder_ave + sorder(imol)

            ! Selected orient-coordinate is saved in angle
            !
            if (option%orient_type == OrientTypeCOSTHETA) then
              angle(imol) = cost(imol)
            else if (option%orient_type == OrientTypeTHETA) then
              angle(imol) = theta(imol)
            end if
          end do
!$omp     end do
!$omp     end parallel

          do imol = 1, nmol(ID_B)
            ig  = nint((angle(imol) - xsta)/ dg)
            igc = (angle(imol) - xsta)/ dg

            if (ig < 0 .or. ig > ngrid) then
              write(iw,'("istep = ", i0)') istep_tot
              write(iw,'("imol  = ", i0)') imol
              write(iw,'("ig = ",i0)') ig
              write(iw,'("cost  = ",f15.7)') cost(imol)
            else
              if (igc == ngrid - 1) then
                distr(igc)     = distr(igc)     + 0.5d0
                distr(igc + 1) = distr(igc + 1) + 0.5d0
              else
                distr(ig) = distr(ig) + 1.0d0
              end if
            end if
            
          end do

          ! Print out Time-series
          !
          write(io,'(f20.10)', advance='no') traj(1)%dt * istep_tot
          do imol = 1, nmol(ID_B) 
            write(io,'(f20.10)',advance='no') angle(imol)
          end do
          write(io,*)

        end do

        istep_tot = istep_tot - 1
        call close_trajfile(trajtype, io_t, dcd, xtc, nc)

      end do
      close(io)

      nstep_tot = istep_tot

      do ig = 0, ngrid
        distr(ig) = distr(ig) / (nstep_tot * dg * nmol(ID_B))
      end do

      sorder_ave = sorder_ave / (nstep_tot * nmol(ID_B))

      ! Normalize
      !
      sumd = 0.0d0
      do ig = 0, ngrid
        sumd = sumd + distr(ig) * dg
      end do
      distr = distr / sumd
      
      ! Output 
      !
      write(fname,'(a,".distr")') trim(output%fhead) 
      call open_file(fname, io)
      do ig = 0, ngrid
        write(io,'(2(e20.10,2x))') xsta + dble(ig) * dg, distr(ig)
      end do
      close(io)

      ! Deallocate
      !
      deallocate(dp, angle, cost, theta, distr)
      deallocate(sorder)

    end subroutine analyze
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
