!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_voronoi
  use mod_input
  use mod_output
  use mod_ctrl
  use mod_cv
  use mod_grid3d

  ! constants
  !
  real(8), parameter, private :: EPS  = 1.0d-8 
  real(8), parameter, private :: EPS2 = 1.0d-5 

  ! structures
  !
  type :: s_states
    integer              :: nstep
    integer, allocatable :: data(:)
  end type

  ! subroutines
  !
  public :: determine_voronoi_state
  public :: determine_spatial_state

  contains
!-----------------------------------------------------------------------
    subroutine determine_voronoi_state(input, output, option)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in)   :: input
      type(s_output), intent(in)   :: output
      type(s_option), intent(in)   :: option

      type(s_voronoi)             :: voronoi
      type(s_cv),     allocatable :: cv(:)
      type(s_states), allocatable :: states(:)

      integer                 :: colmax
      integer                 :: ifile, istep, icv, nfile
      integer                 :: ib, ic, jdim 
      real(8)                 :: crd(3)
      character(len=MaxChar)  :: fhead_out  


      ! setup
      !
      nfile  = input%ncv
      colmax = maxval(option%xyzcol(1:option%ndim))
      allocate(cv(nfile), states(nfile))

      ! setup voronoi
      !
      write(iw,*)
      write(iw,'("Determine_Voronoi_State> Setup Voronoi")')

      call setup_voronoi(option%ndim, option%vr_ncell, option%vr_normvec, voronoi)
      voronoi%cellpos(1:3, 1:voronoi%ncell) &
        = option%vr_cellpos(1:3, 1:voronoi%ncell)
      write(iw,*)
      write(iw,'("Input minimum positions")')
      do ic = 1, voronoi%ncell
          write(iw,'(i3,2x,3f20.10)') &
            ic, (option%vr_cellpos(jdim, ic), jdim = 1, option%ndim)
            !ic, (voronoi%cellpos(jdim, ic), jdim = 1, option%ndim)
      end do
      write(iw,*)
      write(iw,'(">> Finished")')

      ! read cv files
      !
      write(iw,*)
      write(iw,'("> Read CV file")')
      do ifile = 1, nfile
        call read_cv(input%fcv(ifile), colmax, cv(ifile))
      end do
      write(iw,'(">> Finished")')

      ! memory allocation
      !
      do ifile = 1, nfile
        states(ifile)%nstep = cv(ifile)%nstep
        allocate(states(ifile)%data(cv(ifile)%nstep))
      end do

      write(iw,*)
      write(iw,'("> Determine states")')

      ! determine 
      !
      do ifile = 1, nfile
        !do istep = option%nstart, cv(ifile)%nstep
        do istep = 1, cv(ifile)%nstep

          do icv = 1, option%ndim
            crd(icv) = cv(ifile)%data(option%xyzcol(icv), istep)
          end do

          if (option%ndim == 3) then
            states(ifile)%data(istep) = voronoi_state(voronoi, crd(1), crd(2), z = crd(3))
          else
            states(ifile)%data(istep) = voronoi_state(voronoi, crd(1), crd(2))
          end if
        end do
      end do

      write(iw,'(">> Finished")')


      do ifile = 1, nfile
        write(fhead_out, '(a,i4.4,".dat")') trim(output%fhead), ifile
        open(UnitOUT, file=trim(fhead_out))
        do istep = 1, states(ifile)%nstep
          write(UnitOUT,'(2i10)') istep, states(ifile)%data(istep) 
        end do
        close(UnitOUT)
      end do

      deallocate(cv, states)

    end subroutine determine_voronoi_state 
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine determine_spatial_state(input, output, option)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in)   :: input
      type(s_output), intent(in)   :: output
      type(s_option), intent(in)   :: option

      type(s_func3d)    :: sdf

      type(s_cv),     allocatable :: cv(:)
      type(s_states), allocatable :: states(:)

      integer                 :: colmax
      integer                 :: ifile, istep, icv, nfile
      integer                 :: ib, ic
      integer                 :: igx, igy, igz, gind(3)
      real(8)                 :: crd(3), val
      logical                 :: is_inside
      character(len=MaxChar)  :: fhead_out  


      ! setup
      !
      nfile  = input%ncv
      colmax = maxval(option%xyzcol(1:option%ndim))
      allocate(cv(nfile), states(nfile))

      ! read cv files
      !
      write(iw,*)
      write(iw,'("> Read CV file")')
      do ifile = 1, nfile
        call read_cv(input%fcv(ifile), colmax, cv(ifile))
      end do
      write(iw,'(">> Finished")')

      ! read dx file
      !
      write(iw,*)
      write(iw,'("> Read DX file")')
      call read_dx(input%fdx, sdf)
      write(iw,'(">> Finished")')

      ! memory allocation
      !
      do ifile = 1, nfile
        states(ifile)%nstep = cv(ifile)%nstep
        allocate(states(ifile)%data(cv(ifile)%nstep))
      end do

      write(iw,*)
      write(iw,'("> Determine states")')

      ! determine 
      !
      do ifile = 1, nfile
        !do istep = option%nstart, cv(ifile)%nstep
        do istep = 1, cv(ifile)%nstep

          do icv = 1, option%ndim
            crd(icv) = cv(ifile)%data(option%xyzcol(icv), istep)
          end do

          if (option%pbc) then
            do icv = 1, option%ndim
              crd(icv) = crd(icv) &
                       - option%box(icv) * anint(crd(icv) / option%box(icv))
            end do
          end if

          is_inside = .true.
          do icv = 1, option%ndim
            gind(icv) = (crd(icv) - sdf%origin(icv)) / sdf%del(icv) + 1

            if (gind(icv) < 1 .or. gind(icv) > sdf%ng3(icv)) then
              is_inside = .false.
            end if

          end do

          if (is_inside) then
            val = sdf%data(gind(1), gind(2), gind(3))

            if (abs(val) >= EPS2) then
              states(ifile)%data(istep) = 1
            else
              states(ifile)%data(istep) = 0
            end if
          else
            states(ifile)%data(istep)   = 0
          end if

        end do
      end do

      write(iw,'(">> Finished")')


      do ifile = 1, nfile
        write(fhead_out, '(a,".dat")') trim(output%fhead)
        open(UnitOUT, file=trim(fhead_out))
        do istep = 1, states(ifile)%nstep
          write(UnitOUT,'(2i10)') istep, states(ifile)%data(istep) 
        end do
        close(UnitOUT)
      end do

      deallocate(cv, states)

    end subroutine determine_spatial_state
!-----------------------------------------------------------------------
!
end module mod_analyze
!=======================================================================
