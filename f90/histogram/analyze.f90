!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_cv
  use mod_ctrl

  ! subroutines
  !
  public :: analyze 

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in)   :: input
      type(s_output), intent(in)   :: output
      type(s_option), intent(in)   :: option

      ! IO
      !
      integer                :: io
      character(len=MaxChar) :: fname

      ! Local
      !
      integer :: nfile, ncol, ista, iend, nstep
      real(8) :: val, ave, dev 

      integer :: nx
      real(8) :: xsta, dx 

      ! Dummy
      !
      integer                :: ifile, istep, ix, icol
      real(8)                :: x
      character(len=MaxChar) :: line

      ! Arrays
      !
      type(s_cv), allocatable :: cv(:) 
      real(8),    allocatable :: distr(:), distr_each(:, :)
      real(8),    allocatable :: sterr(:) 


      ! Read CV files
      !
      nfile = input%ncv
      ncol  = option%ncol
      ista  = option%data_sta
      iend  = option%data_end

      allocate(cv(nfile))
     
      if (iend /= 0) then 

        do ifile = 1, nfile
          call read_cv(input%fcv(ifile),             &
                       ncol,                         &
                       cv(ifile),                    &
                       nsta = ista,                  &
                       nstep_read = iend - ista + 1) 
        end do

      else

        do ifile = 1, nfile
          call read_cv(input%fcv(ifile),             &
                       ncol,                         &
                       cv(ifile),                    &
                       nsta = ista)
        end do

      end if

      ! Calculate distribution 
      !
      nx   = option%nx
      xsta = option%xsta
      dx   = option%dx

      allocate(distr(0:nx), distr_each(0:nx, nfile))
      allocate(sterr(0:nx))

      distr_each = 0.0d0

      do ifile = 1, nfile
        nstep = cv(ifile)%nstep

        do istep = 1, nstep
          do icol = 1, ncol
            val = cv(ifile)%data(icol, istep)
            ix = nint((val - xsta) / dx)

            if (ix > 0 .and. ix < nx) then
              distr_each(ix, ifile) = distr_each(ix, ifile) + 1.0d0 
            else if (ix == 0 .or. ix == option%nx) then
              distr_each(ix, ifile) = distr_each(ix, ifile) + 2.0d0 
            end if

          end do
        end do

        distr_each(:, ifile) = distr_each(:, ifile) / (dx * nstep * ncol)

      end do

      ! Average
      !
      do ix = 0, nx
        distr(ix) = sum(distr_each(ix, :)) / dble(nfile)
      end do


      do ix = 0, nx

        ave = distr(ix)
        dev = 0.0d0

        do ifile = 1, nfile
          val = distr_each(ix, ifile)
          dev = dev + (val - ave)**2
        end do

        dev       = sqrt(dev / (nfile - 1))
        sterr(ix) = dev / sqrt(dble(nfile))

      end do 
      
      ! Output
      !
      write(fname,'(a,".",a)') trim(output%fhead), trim(output%file_extension)
      call open_file(fname, io)

      if (nfile > 3) then
        do ix = 0, option%nx
          x = option%xsta + ix * option%dx
          write(io,'(3(e20.10,2x))') x, distr(ix), sterr(ix)
        end do
      else
        do ix = 0, option%nx
          x = option%xsta + ix * option%dx
          write(io,'(2(e20.10,2x))') x, distr(ix)
        end do
      end if
      close(io)

      ! Deallocate
      !
      deallocate(cv)
      deallocate(distr, distr_each, sterr)


    end subroutine analyze
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
