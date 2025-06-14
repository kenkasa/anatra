!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_akima
  use mod_grid3d
  use mod_output
  use mod_ctrl

  ! constants
  !
  integer, parameter, public :: Ndim = 3

  ! structures
  !
  type :: s_sdbin
    real(8)        :: norm
    !integer        :: nstep_tot
    integer        :: nmol
    real(8)        :: vol
    real(8)        :: dv
    logical        :: use_spline
    integer        :: spline_resolution
    type(s_func3d) :: func3d 
  end type s_sdbin

  ! subroutines
  !
  public :: analyze 

  contains
!-----------------------------------------------------------------------
    subroutine analyze(output, option, sdinfo)
!-----------------------------------------------------------------------
      implicit none

      type(s_output), intent(in)   :: output
      type(s_option), intent(in)   :: option 
      type(s_sdinfo), intent(in)   :: sdinfo

      type(s_sdbin), allocatable   :: sdbin(:)
      type(s_func3d)               :: sdf0, sdf1

      integer :: ifile

      integer :: nstep_tot, nmol
      real(8) :: vol, dv, norm_tot
      logical :: use_spline
      integer :: spline_resolution
      integer :: ngrids_nonzero
      integer :: igx, igy, igz, ng3(3)
      real(8) :: del(3), origin(3), box(3)
      real(8) :: val, vol_nonzero

      integer :: ng3_spl(3)
      real(8) :: del_spl(3)

      character(len=MaxChar) :: fhead_out


      allocate(sdbin(sdinfo%nfile))
      do ifile = 1, sdinfo%nfile
        open(UnitOut, file=trim(sdinfo%fsdbin(ifile)), form='unformatted')
          read(UnitOut) sdbin(ifile)%norm
          !read(UnitOut) sdbin(ifile)%nstep_tot
          read(UnitOut) nmol
          read(UnitOut) sdbin(ifile)%vol
          read(UnitOut) dv
          read(UnitOut) use_spline
          read(UnitOut) spline_resolution
          read(UnitOut) ng3
          read(UnitOut) del 
          read(UnitOut) origin
          read(UnitOut) box

          !sdbin(ifile)%norm              = norm
          sdbin(ifile)%nmol              = nmol
          sdbin(ifile)%dv                = dv
          sdbin(ifile)%use_spline        = use_spline
          sdbin(ifile)%spline_resolution = spline_resolution
         
          sdbin(ifile)%func3d%ng3        = ng3   
          sdbin(ifile)%func3d%del        = del  
          sdbin(ifile)%func3d%origin     = origin
          sdbin(ifile)%func3d%box        = box  
          allocate(sdbin(ifile)%func3d%data(ng3(1), ng3(2), ng3(3)))
          read(UnitOut) sdbin(ifile)%func3d%data
        close(UnitOut)
      end do

      norm_tot          = sum(sdbin(:)%norm)
      !nstep_tot         = sum(sdbin(:)%nstep_tot)
      vol               = sum(sdbin(:)%vol) / sdinfo%nfile

      !if (option%nstep_tot /= 0) then
      if (option%norm_tot > 0.0d0) then
        norm_tot = option%norm_tot
        !nstep_tot = option%nstep_tot
      end if
      !write(iw,*)
      !write(iw,'(" Analyze> Show parameters read from sdbin files")')
      !write(iw,'(" vol           = ",f20.10)') vol
      !write(iw,*)

      allocate(sdf0%data(ng3(1), ng3(2), ng3(3)))
      sdf0%ng3    = ng3
      sdf0%del    = del
      sdf0%origin = origin
      sdf0%box    = box 

      sdf0%data = 0.0d0
      do ifile = 1, sdinfo%nfile
        sdf0%data(:, :, :) = sdf0%data(:, :, :) + sdbin(ifile)%func3d%data 
      end do

      write(iw,'("vol       = ",f20.10)')  vol
      write(iw,'("dv        = ",f20.10)')  dv
      write(iw,'("nmol      = ",i0)')      nmol 
      write(iw,'("norm_tot  = ",e15.7)')   norm_tot
      !write(iw,'("nstep_tot = ",i0)')      nstep_tot 
      !sdf0%data(:, :, :) = vol * sdf0%data(:, :, :) / (nmol * nstep_tot * dv)
      !sdf0%data(:, :, :) = vol * sdf0%data(:, :, :)  / dble(nmol * dv * nstep_tot)
      sdf0%data(:, :, :) = vol * sdf0%data(:, :, :)  / dble(norm_tot * dv)
      !/ (nmol * nstep_tot * dv)

      if (use_spline) then
        ng3_spl = 1
        del_spl = 1.0d0
        ng3_spl(1:Ndim) = ng3(1:Ndim) * spline_resolution
        del_spl(1:Ndim) = del(1:Ndim) / spline_resolution 
        call setup_func3d(ng3_spl, del_spl, origin, sdf1)
        call akima3d(sdf0%ng3, sdf1%ng3, sdf0%box, 0.0d0, 1.0d20, &
                     sdf0%data, sdf1%data)  
      end if

      ! Calculate the volume in which g is larger than zero
      !
      ngrids_nonzero = 0
      do igz = 1, ng3(3)
        do igy = 1, ng3(2)
          do igx = 1, ng3(1)
            val = sdf0%data(igx, igy, igz)  
            if (abs(val) >= 1.0d-5) &
              ngrids_nonzero = ngrids_nonzero + 1
          end do
        end do
      end do

      vol_nonzero = ngrids_nonzero * dv

      write(iw,'(" Volume where g(r) is larger than zero / A^3 = ", f20.10)') &
        vol_nonzero

      write(fhead_out,'(a,"_g_original")') trim(output%fhead)
      call write_dx(fhead_out, sdf0)
      if(use_spline) then
        write(fhead_out,'(a,"_g_spline")')   trim(output%fhead)
        call write_dx(fhead_out, sdf1)
      end if


    end subroutine analyze
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
