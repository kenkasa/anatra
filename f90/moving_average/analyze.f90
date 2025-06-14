!=======================================================================
module mod_movave_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_cv
  use mod_input
  use mod_output
  use mod_movave_ctrl

  ! structures
  !

  ! subroutines
  !
  public :: movave_analyze
  public :: movave_write

  contains
!-----------------------------------------------------------------------
    subroutine movave_analyze(input, output, option, movave, finp)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),                    intent(in)    :: input
      type(s_output),                   intent(in)    :: output
      type(s_movave_option),            intent(in)    :: option
      type(s_movave),                   intent(inout) :: movave(:)
      character(len=MaxChar), optional, intent(in)    :: finp


      ! IO
      !
      integer                :: io
      character(len=MaxChar) :: fname

      ! Local
      !
      integer :: nfile, nstep
      integer :: nsplit(0:MaxSep)

      ! Dummy
      !
      integer :: i, ifile, istep, isep
      integer :: ista, iend, iorg
      integer :: reg_sta, reg_end
      real(8) :: rdum, diff

      ! Arrays
      !
      type(s_cv), allocatable :: cv(:)
      integer,    allocatable :: range(:, :)
      real(8),    allocatable :: input_data(:)


      ! Read CV files
      !
      nfile = input%ncv

      allocate(cv(nfile))

      do ifile = 1, nfile
        call read_cv(input%fcv(ifile),              &
                     1,                             &
                     cv(ifile))
      end do

      ! Setup 
      !
      !   Read CV file
      !
      do ifile = 1, nfile
        nstep               = cv(ifile)%nstep
        movave(ifile)%ngrid = nstep

        allocate(movave(ifile)%grid (0:nstep - 1), &
                 movave(ifile)%data (0:nstep - 1), &
                 movave(ifile)%deriv(0:nstep - 1))

        do istep = 1, nstep
          read(cv(ifile)%x(istep), *) movave(ifile)%grid(istep - 1)
        end do 

      end do

      !   Setup starting point for averaging
      !
      iorg = 1
      if (option%include_zero) then
        iorg = 0
      end if

      !   Setup boundary point
      !
      nstep = cv(1)%nstep

      if (option%nregion == 1) then
        nsplit(1) = nstep - 1
        !if (option%include_zero) then
          nsplit(0) = iorg
        !end if
      else
        do i = 1, option%nregion - 1
          nsplit(i) = nint((option%t_sep(i) - option%t_sta) / option%dt) 
        end do
        nsplit(0)              = iorg
        nsplit(option%nregion) = nstep - 1 
      end if

      !   Setup average range for each point
      !
      allocate(range(1:3, 0:nstep - 1))
      range = 0

      !   ... Initialize
      do istep = 0, nstep - 1
        range(1:2, istep) = istep  ! 1: start , 2: end
        range(3,   istep) = 1      ! # of points used for averaging
      end do

      !   ... Determine the average range for each point
      do isep = 0, option%nregion - 1 
        reg_sta = nsplit(isep) 
        reg_end = nsplit(isep + 1) 

        do istep = reg_sta, reg_end 
          ista = istep - option%npoint(isep + 1) / 2
          iend = istep + option%npoint(isep + 1) / 2
       
          if (ista < 0) then
            ista = 0
          end if
       
          if (iend > nstep - 1) then
            iend = nstep - 1
          end if
       
          range(1, istep) = ista
          range(2, istep) = iend
          range(3, istep) = iend - ista + 1
       
        end do
      end do

      ! Perform averaging
      !
      do ifile = 1, nfile
        movave(ifile)%data(:)  = 0.0d0
        movave(ifile)%deriv(:) = 0.0d0
      end do

      do ifile = 1, nfile

        nstep = cv(ifile)%nstep

        if (.not. option%include_zero) then
          movave(ifile)%data(0) = cv(ifile)%data(1, 1)
        end if

        if (.not. allocated(input_data)) then
          allocate(input_data(0:nstep - 1))
        end if

        input_data(0:nstep - 1) = cv(ifile)%data(1, 1:nstep)

        do istep = iorg, nstep - 1
          ista = range(1, istep)
          iend = range(2, istep)
          movave(ifile)%data(istep) &
            = sum(input_data(ista:iend)) / range(3, istep)
        end do

        do istep = 0, nstep - 2
          diff                       =   movave(ifile)%data(istep + 1) &
                                       - movave(ifile)%data(istep)
          movave(ifile)%deriv(istep) =   diff / option%dt 
        end do

        diff                           =   movave(ifile)%data(nstep) &
                                         - movave(ifile)%data(nstep - 1)
        movave(ifile)%deriv(nstep - 1) =   diff / option%dt

        deallocate(input_data)
      end do 

    end subroutine movave_analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine movave_write(output, movave)
!-----------------------------------------------------------------------
      implicit none

      type(s_output), intent(in) :: output
      type(s_movave), intent(in) :: movave(:)

      ! IO
      !
      integer :: io

      ! Local
      !
      integer                :: nfile
      character(len=MaxChar) :: fwrite, fmt_str, ext

      ! Dummy
      !
      integer                :: ifile, istep


      nfile = size(movave(:))
      ext   = output%file_extension

      if (trim(ext) == '') then
        ext = 'dat'
      end if 

      if (nfile == 1) then
        write(fwrite,'(a,".",a)') trim(output%fhead), trim(ext) 
      else
        write(fmt_str,'("(a,i",i0,".",i0,",""."",a)")') output%ndigit, output%ndigit
      end if

      do ifile = 1, nfile

        if (nfile == 1) then
          call open_file(fwrite, io)
        else
          write(fwrite, fmt=trim(fmt_str)) &
            trim(output%fhead), ifile, trim(ext)
          call open_file(fwrite, io)
        end if

        do istep = 0, movave(ifile)%ngrid - 1
          write(io,'(3(e15.7,2x))')     &
            movave(ifile)%grid (istep), &
            movave(ifile)%data (istep), &
            movave(ifile)%deriv(istep)
        end do

        close(io)
      end do

    end subroutine movave_write
!-----------------------------------------------------------------------

end module mod_movave_analyze
!=======================================================================
