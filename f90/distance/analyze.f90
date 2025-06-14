!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_ctrl
  use mod_traj
  use mod_xtcio
  use mod_dcdio
  use mod_netcdfio
  use mod_com

  ! subroutines
  !
  public  :: analyze
  private :: get_dist
  private :: get_mindist

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option, traj)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in)    :: input
      type(s_output), intent(in)    :: output
      type(s_option), intent(in)    :: option
      type(s_traj),   intent(inout) :: traj(2)

      type(s_dcd)    :: dcd 
      type(xtcfile)  :: xtc 
      type(s_netcdf) :: nc 
      type(s_com)    :: com(2)

      ! I/O
      !
      integer                :: io_d, io_cd, io_t
      character(len=MaxChar) :: fname_d, fname_cd

      ! Local
      !
      integer :: natm, nmol(2)
      integer :: trajtype
      real(8) :: m, d2, d(3), ti, tj
      logical :: is_end

      ! Dummy
      !
      integer :: i, istep, istep_tot, istep_use
      integer :: itraj, imol, jmol

      ! Arrays
      !
      real(8), allocatable :: dist(:, :)
      real(8), allocatable :: dist_closest_pair(:)


      do itraj = 1, 2
        call get_com(option%mode(itraj),   &
                     traj(itraj),          &
                     com(itraj),           &
                     setup = .true.,       &
                     calc_coord = .false., &
                     myrank = 0)
        nmol(itraj) = com(itraj)%nmol
      end do

      ! Allocate memory
      !
      allocate(dist(nmol(2), nmol(1)))
      if (option%distance_type == DistanceTypeMINIMUM) then
        allocate(dist_closest_pair(nmol(1)))
      end if

      ! Calculate distance 
      !
      write(fname_d,'(a,".dis")') trim(output%fhead)
      call open_file(fname_d, io_d)

      if (option%distance_type == DistanceTypeMINIMUM) then
        write(fname_d,'(a,".dis_closest")') trim(output%fhead)
        call open_file(fname_cd, io_cd)
      end if

      ! Get trajectory types
      !
      call get_trajtype(input%ftraj(1), trajtype)

      write(iw,*)
      write(iw,'("Analyze> Start")')

      dist = 0.0d0
      istep_tot = 0
      istep_use = 0
      do itraj = 1, input%ntraj
        call open_trajfile(input%ftraj(itraj), trajtype, io_t, dcd, xtc, nc)
        call init_trajfile(trajtype, io_t, dcd, xtc, nc, natm)

        is_end = .false.
        istep  = 0
        do while (.not. is_end)
          istep     = istep     + 1
          istep_tot = istep_tot + 1

          if (mod(istep_tot, 1) == 100) then
            write(iw,'("Progress: ",i0)') istep_tot
          end if

          call read_trajfile_oneframe(trajtype, io_t, istep, dcd, xtc, nc, is_end)

          if (is_end) exit

          ti = traj(1)%dt * istep_tot

          if (ti > option%t_end) then
            is_end = .true.
            exit
          end if

          if (ti >= option%t_sta .and. ti <= option%t_end) then
            istep_use = istep_use + 1
          else
            cycle  
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

          if (option%distance_type == DistanceTypeSTANDARD) then
         
            ! Calculate pair distance
            !
            call get_dist(1, option%pbc, traj(1)%box(1:3, 1), &
                          nmol, com, dist)
         
          else if (option%distance_type == DistanceTypeMINIMUM) then
         
            ! Calculate pair minimum distance
            !
            call get_mindist(1, option, traj(1)%box(1:3, 1),   &
                             nmol, com, traj, dist, dist_closest_pair)
             
          else if (option%distance_type == DistanceTypeINTRA) then
         
            ! Calculate intramolecular distance
            !
            call get_intradist(1, option%pbc, traj(1)%box(1:3, 1), &
                               nmol, com, dist)
         
          end if
         
          ! Write pair distance at istep  
          !
         
          if (     option%distance_type == DistanceTypeSTANDARD      &
              .or. option%distance_type == DistanceTypeMINIMUM) then
            write(io_d,'(f20.10)',advance="no") traj(1)%dt * istep_use
            do imol = 1, nmol(1)
              do jmol = 1, nmol(2)
                write(io_d,'(f20.10)',advance="no") dist(jmol, imol)
              end do
            end do
            write(io_d,*)
         
            if (option%distance_type == DistanceTypeMINIMUM) then
              write(io_cd,'(f20.10)',advance="no") traj(1)%dt * istep_use
              do imol = 1, nmol(1)
                write(io_cd,'(f20.10)',advance="no") &
                  dist_closest_pair(imol)
              end do
              write(io_cd,*)
            end if
         
          else if (option%distance_type == DistanceTypeINTRA) then
            write(io_d,'(f20.10)',advance="no") traj(1)%dt * istep_use
            do imol = 1, nmol(1)
              write(io_d,'(f20.10)',advance="no") dist(imol, imol)
            end do
            write(io_d,*)
          end if

        end do

        istep_tot = istep_tot - 1

      end do

      close(io_d)
      if (option%distance_type == DistanceTypeMINIMUM) &
        close(io_cd)

      ! Deallocate
      !
      deallocate(dist)

      if (option%distance_type == DistanceTypeMINIMUM) &
        deallocate(dist_closest_pair)


    end subroutine analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_dist(istep, pbc, box, nmol, com, dist)
!-----------------------------------------------------------------------
      implicit none

      integer,      intent(in)  :: istep
      logical,      intent(in)  :: pbc
      real(8),      intent(in)  :: box(3)
      integer,      intent(in)  :: nmol(2)
      type(s_com),  intent(in)  :: com(2)
      real(8),      intent(out) :: dist(:,:)

      ! Local
      !
      real(8) :: d(3)

      ! Dummy
      ! 
      integer :: imol, jmol

     
      dist = 0.0d0

      if(pbc) then
        do imol = 1, nmol(1)
          do jmol = 1, nmol(2)
            d(1:3) = com(2)%coord(1:3, jmol, istep) &
                   - com(1)%coord(1:3, imol, istep)
            d(1:3) = d(1:3) - box(1:3) * nint(d(1:3) / box(1:3))

            dist(jmol, imol)  = sqrt(dot_product(d, d))
          end do
        end do
      else
        do imol = 1, nmol(1)
          do jmol = 1, nmol(2)
            d(1:3) = com(2)%coord(1:3, jmol, istep) &
                   - com(1)%coord(1:3, imol, istep)

            dist(jmol, imol) = sqrt(dot_product(d, d))
          end do
        end do
      end if


    end subroutine get_dist
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_mindist(istep, option, box, nmol, com, traj, dist, &
                           dist_closest_pair)
!-----------------------------------------------------------------------
      implicit none

      integer,        intent(in)  :: istep
      type(s_option), intent(in)  :: option
      real(8),        intent(in)  :: box(3)
      integer,        intent(in)  :: nmol(2)
      type(s_com),    intent(in)  :: com(2) 
      type(s_traj),   intent(in)  :: traj(2)
      real(8),        intent(out) :: dist(:,:)
      real(8),        intent(out) :: dist_closest_pair(:)

      ! Local
      !
      real(8) :: d(3), r, dist_min

      ! Dummy
      !
      integer :: iatm, jatm
      integer :: imol, jmol

     
      dist = 0.0d0

      if (option%pbc) then
        if (option%mindist_type(1) == MinDistTypeSITE &
            .and. option%mindist_type(2) == MinDistTypeSITE) then
          do imol = 1, nmol(1)
            do jmol = 1, nmol(2)
              dist_min = 1.0d10 
              do iatm = com(1)%molsta(imol), com(1)%molend(imol)
                do jatm = com(2)%molsta(jmol), com(2)%molend(jmol)
                  d(1:3) = traj(2)%coord(1:3, jatm, istep) &
                         - traj(1)%coord(1:3, iatm, istep)
                  d(1:3) = d(1:3) - box(1:3) * nint(d(1:3) / box(1:3))
                  r      = sqrt(dot_product(d, d))
         
                  if (r <= dist_min) &
                    dist_min = r
         
                end do
              end do
              dist(jmol, imol) = dist_min
         
            end do
          end do
        else if (option%mindist_type(1) == MinDistTypeSITE &
            .and. option%mindist_type(2) == MinDistTypeCOM) then
          do imol = 1, nmol(1)
            do jmol = 1, nmol(2)
              dist_min = 1.0d10 
              do iatm = com(1)%molsta(imol), com(1)%molend(imol)
                d(1:3) = com(2)%coord(1:3, jmol, istep) &
                       - traj(1)%coord(1:3, iatm, istep)
                d(1:3) = d(1:3) - box(1:3) * nint(d(1:3) / box(1:3))
                r      = sqrt(dot_product(d, d))
         
                if (r <= dist_min) &
                  dist_min = r
         
              end do
              dist(jmol, imol) = dist_min
         
            end do
          end do

        else if (option%mindist_type(1) == MinDistTypeCOM   &
          .and.  option%mindist_type(2) == MinDistTypeSITE) then 
          do imol = 1, nmol(1)
            do jmol = 1, nmol(2)
              dist_min = 1.0d10 
              do jatm = com(2)%molsta(jmol), com(2)%molend(jmol)
                d(1:3) = traj(2)%coord(1:3, jatm, istep) &
                       - com(1)%coord(1:3, imol, istep) 
                d(1:3) = d(1:3) - box(1:3) * nint(d(1:3) / box(1:3))
                r      = sqrt(dot_product(d, d))
         
                if (r <= dist_min) &
                  dist_min = r
         
              end do
              dist(jmol, imol) = dist_min
         
            end do
          end do
        end if
      else
        if (option%mindist_type(1) == MinDistTypeSITE &
            .and. option%mindist_type(2) == MinDistTypeSITE) then
          do imol = 1, nmol(1)
            do jmol = 1, nmol(2)
              dist_min = 1.0d10 
              do iatm = com(1)%molsta(imol), com(1)%molend(imol)
                do jatm = com(2)%molsta(jmol), com(2)%molend(jmol)
                  d(1:3) = traj(2)%coord(1:3, jatm, istep) &
                         - traj(1)%coord(1:3, iatm, istep)
                  r      = sqrt(dot_product(d, d))
         
                  if (r <= dist_min) &
                    dist_min = r
         
                end do
              end do
              dist(jmol, imol) = dist_min
            end do
          end do
        else if (option%mindist_type(1) == MinDistTypeSITE &
            .and. option%mindist_type(2) == MinDistTypeCOM) then
          do imol = 1, nmol(1)
            do jmol = 1, nmol(2)
              dist_min = 1.0d10 
              do iatm = com(1)%molsta(imol), com(1)%molend(imol)
                d(1:3) = com(2)%coord(1:3, jmol, istep) &
                       - traj(1)%coord(1:3, iatm, istep)
                r      = sqrt(dot_product(d, d))
         
                if (r <= dist_min) &
                  dist_min = r
         
              end do
              dist(jmol, imol) = dist_min
            end do
          end do
        else if (option%mindist_type(1) == MinDistTypeCOM) then 
          write(iw,'("Get_Mindist> Error.")')
          write(iw,'("mindist_type(1) = COM is not supported")')
          stop
        end if
      end if

      do imol = 1, nmol(1)
        dist_closest_pair(imol) = minval(dist(:, imol))
      end do


    end subroutine get_mindist
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_intradist(istep, pbc, box, nmol, com, dist)
!-----------------------------------------------------------------------
      implicit none

      integer,      intent(in)  :: istep
      logical,      intent(in)  :: pbc
      real(8),      intent(in)  :: box(3)
      integer,      intent(in)  :: nmol(2)
      type(s_com),  intent(in)  :: com(2)
      real(8),      intent(out) :: dist(:, :)

      ! Local
      !
      real(8) :: d(3)

      ! Dummy
      ! 
      integer :: imol

     
      dist = 0.0d0

      if(pbc) then
        do imol = 1, nmol(1)
          d(1:3) = com(2)%coord(1:3, imol, istep) &
                 - com(1)%coord(1:3, imol, istep)
          d(1:3) = d(1:3) - box(1:3) * nint(d(1:3) / box(1:3))

          dist(imol, imol)  = sqrt(dot_product(d, d))
        end do
      else
        do imol = 1, nmol(1)
          d(1:3) = com(2)%coord(1:3, imol, istep) &
                 - com(1)%coord(1:3, imol, istep)

          dist(imol, imol) = sqrt(dot_product(d, d))
        end do
      end if


    end subroutine get_intradist
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
