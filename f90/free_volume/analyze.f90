!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_ctrl
  use mod_traj
  use mod_anaparm
  use mod_prmtop
  use mod_xtcio
  use mod_dcdio
  use mod_com
  use mod_potential
  use mod_random
  !use mod_fftw3i

  ! constants
  !
  integer, parameter :: ID_T = 1
  integer, parameter :: ID_C = 2
  real(8), parameter :: maxrad = 5.0

  ! structures
  !
  type :: s_neighbor
    integer :: natm
    integer, allocatable :: ind(:)
    real(8), allocatable :: rad(:) 
    real(8), allocatable :: pos(:, :) 
  end type s_neighbor

  ! subroutines
  !
  public  :: analyze

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option, traj)
!-----------------------------------------------------------------------
      implicit none

      ! parameters
      real(8), parameter :: charge_const = 332.05221729d0

      ! formal arguments
      type(s_input),  intent(in)    :: input
      type(s_output), intent(in)    :: output
      type(s_option), intent(in)    :: option
      type(s_traj),   intent(inout) :: traj(2)

      ! I/O
      !
      integer                :: io
      character(len=MaxChar) :: fname

      type(s_dcd)            :: dcd
      type(xtcfile)          :: xtc
      type(s_netcdf)         :: nc 
      type(s_com)            :: com(2)
      type(s_prmtop)         :: prmtop
      type(s_anaparm)        :: anaparm
      type(s_pot)            :: pot
      type(s_neighbor)       :: nb

      integer                :: iseed
      integer                :: i, j, ix, ixyz, ins
      integer                :: itraj, istep, istep_tot
      integer                :: iatm, jatm
      integer                :: natm, nmol(2)

      integer                :: ngrid, nins
      real(8)                :: xsta, dx, x, diff, d2, zc, L, rand
      real(8)                :: dvec(3), ins_pos(3)
      logical                :: do_overlap

      integer                :: trajtype
      logical                :: is_end
      logical                :: is_first

      ! Arrays
      !
      real(8), allocatable   :: gx(:), cnt(:)
      real(8), allocatable   :: gx_p(:), cnt_p(:)



      ! Read parameter file
      !
      write(iw,'("")')
      write(iw,'("Analyze> Read potential parameter files")')
      if (option%parmformat == ParmFormatANAPARM) then
        call read_anaparm(input%fanaparm, anaparm)
      else if (option%parmformat == ParmFormatPRMTOP) then
        call read_prmtop(input%fprmtop, prmtop)
      end if

      ! Setup
      !
      iseed = - 1
      ngrid = option%ngrid
      dx    = option%dx
      xsta  = option%xsta
      nins  = option%nins

      ! Generate random seed
      !
      call get_seed         (iseed)
      call initialize_random(iseed)

      write(iw,'("")')
      write(iw,'("Analyze> Setup molecules")')
      call get_com(option%mode(1),       &
                   traj(1),              &
                   com(1),               &
                   setup = .true.,       &
                   calc_coord = .false., &
                   myrank = 0)

      call get_com(option%mode(2),       &
                   traj(2),              &
                   com(2),               &
                   setup = .true.,       &
                   calc_coord = .false., &
                   myrank = 0)

      nmol(1) = com(1)%nmol
      nmol(2) = com(2)%nmol

      ! Prepare LJ coefficients & point charges
      !
      write(iw,'("")')
      write(iw,'("Analyze> Setup potential")')
      
      call alloc_pot((/traj(ID_T)%natm, traj(ID_T)%natm/), pot, diagonal = .true.)
      call setup_pot((/traj(ID_T), traj(ID_T)/), pot, prmtop = prmtop, diagonal = .true.)

      ! Hidden: output potential parameter (stop after outputting)
      !
      !if (option%output_param) then
      !  call output_pot_parameter(output, option, com, traj, pot)
      !end if

      if (option%coord_type == CoordTypeZ) then
        allocate(gx(0:ngrid), cnt(0:ngrid))
        allocate(gx_p(0:ngrid), cnt_p(0:ngrid))
        gx  = 0.0d0
        cnt = 0.0d0
      end if

      ! Get trajectory types
      !
      call get_trajtype(input%ftraj(1), trajtype)

      !
      write(iw,*)
      write(iw,'("Analyze> Start the Free-Volume calculation")')
      !
      istep_tot = 0
      do itraj = 1, input%ntraj
        call open_trajfile(input%ftraj(itraj), trajtype, io, dcd, xtc, nc)
        call init_trajfile(trajtype, io, dcd, xtc, nc, natm)

        is_end = .false.
        istep  = 0
        do while (.not. is_end)
          istep     = istep     + 1
          istep_tot = istep_tot + 1

          if (mod(istep_tot, 100) == 0) then
            write(iw,'("Step ",i0)') istep_tot
          end if

          call read_trajfile_oneframe(trajtype, io, istep, dcd, xtc, nc, is_end)

          if (is_end) exit

          do i = 1, 2
            call send_coord_to_traj(1, trajtype, dcd, xtc, nc, traj(i))
            call get_com(option%mode(i),      &
                         traj(i),             &
                         com(i),              &
                         setup = .false.,     &
                         calc_coord = .true., & 
                         myrank = 1)
          end do

          if (option%coord_type == CoordTypeZ) then
            zc = com(ID_C)%coord(3, 1, 1)

            do ix = 0, ngrid

              x       = xsta + ix * dx

              ! Initialize
              !
              if (allocated(nb%pos)) then
                deallocate(nb%pos, nb%rad, nb%ind)
              end if
              nb%natm = 0

              ! Get # of neighbor atoms 
              !
              do iatm = 1, traj(ID_T)%natm
                diff = abs(x - traj(ID_T)%coord(3, iatm, 1) + zc)
                if (diff < maxrad) then
                  nb%natm = nb%natm + 1
                end if
              end do

              ! Allocate
              !
              allocate(nb%pos(3, nb%natm), nb%rad(nb%natm), nb%ind(nb%natm))

              ! Store atom info.
              !
              jatm = 0
              do iatm = 1, traj(ID_T)%natm
                diff = abs(x - traj(ID_T)%coord(3, iatm, 1) + zc)
                if (diff < maxrad) then
                  jatm              = jatm + 1
                  nb%ind(jatm)      = iatm
                  !nb%rad(jatm)      = pot%ljsgm(iatm, iatm) * 0.5d0
                  nb%rad(jatm)      = pot%ljsgm(iatm, 1) * 0.5d0
                  nb%pos(1:3, jatm) = traj(ID_T)%coord(1:3, iatm, 1)
                end if
              end do

              ! Particle insertion
              !
              gx_p  = 0.0d0
              cnt_p = 0.0d0

              !$omp parallel private(ins, ixyz, iatm, rand, L, ins_pos, do_overlap, dvec, d2) &
              !$omp default(shared), reduction(+:cnt_p), reduction(+:gx_p)
              !$omp do 
              do ins = 1, nins

                ! Generate insertion pos.
                !
                do ixyz = 1, 3
                  call get_random(rand)
                  if (ixyz <= 2) then
                    L             = traj(ID_T)%box(ixyz, 1) 
                    ins_pos(ixyz) = - 0.5d0 * L + rand * L
                  else
                    ins_pos(ixyz) = x + rand * dx
                  end if 
                end do

                cnt_p(ix) = cnt_p(ix) + 1.0d0

                ! Judge
                !
                do_overlap = .false.
                do iatm = 1, nb%natm
                  dvec(1:3) = nb%pos(1:3, iatm) - zc - ins_pos(1:3)
                  dvec(1:3) = dvec(1:3) &
                    - traj(ID_T)%box(1:3, 1) * nint(dvec(1:3) / traj(ID_T)%box(1:3, 1))
                  d2 = dot_product(dvec(1:3), dvec(1:3))

                  if (d2 < nb%rad(iatm)**2) then
                    do_overlap = .true.
                    exit
                  end if 
                end do

                if (.not. do_overlap) then
                  !gx(ix) = gx(ix) + 1.0d0
                  gx_p(ix) = gx_p(ix) + 1.0d0 
                end if

              end do
              !$omp end do
              !$omp end parallel

              gx  = gx  + gx_p
              cnt = cnt + cnt_p

            end do

          end if

        end do ! step

        istep_tot = istep_tot - 1
        call close_trajfile(trajtype, io, dcd, xtc, nc)

      end do   ! traj

      if (option%coord_type == CoordTypeZ) then

        ! Normalize
        !
        do ix = 0, ngrid
          gx(ix) = gx(ix) / cnt(ix) 
        end do

        ! Output
        !
        write(fname, '(a, ".fvz")') trim(output%fhead)
        call open_file(fname, io)
        do ix = 0, ngrid
          x = xsta + ix * dx
          write(io, '(2f20.10)') x, gx(ix)
        end do 
        close(io)
      end if

      call dealloc_pot(pot)

    end subroutine analyze
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
