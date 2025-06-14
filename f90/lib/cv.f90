!=======================================================================
module mod_cv
!=======================================================================
  use mod_util
  use mod_const
  implicit none

  ! constants
  !

  ! structures
  !
  type :: s_cvinfo
    integer                             :: nfile     = 1
    character(len=MaxChar), allocatable :: fcv(:)
    character(len=MaxChar), allocatable :: fweight(:)
  end type s_cvinfo

  type :: s_cv
    integer                             :: nstep
    integer                             :: ndim
    character(len=MaxChar), allocatable :: x(:)
    real(8),                allocatable :: data(:,:)
    real(8),                allocatable :: weight(:)
  end type s_cv

  ! subroutines
  !
  public :: read_ctrl_cvinfo
  public :: read_ctrl_cvinfo_weight
  public :: read_cv

  contains
!-----------------------------------------------------------------------
    subroutine read_ctrl_cvinfo(iunit, cvinfo)
!-----------------------------------------------------------------------
      implicit none
!
      integer,        intent(in)  :: iunit
      type(s_cvinfo), intent(out) :: cvinfo


      integer :: ifile, nfile

      rewind iunit

      ifile = 0
      do while(.true.)
        read(iunit,*,end=100)
        ifile = ifile + 1
      end do
100   continue
      rewind iunit
      nfile = ifile
      allocate(cvinfo%fcv(ifile))
      
      do ifile = 1, nfile
        read(iunit,'(a)') cvinfo%fcv(ifile)
      end do

      write(iw,*)
      write(iw,'(">> CVINFO parameters")')
      write(iw,'("nfile    = ", i0)')  nfile
      if (ifile >= 20) then
        do ifile = 1, 10
          write(iw,'(a)') trim(cvinfo%fcv(ifile))
        end do
        write(iw,'(" ... skipped ...")')
        do ifile = nfile - 10, nfile
          write(iw,'(a)') trim(cvinfo%fcv(ifile))
        end do
      else
        do ifile = 1, nfile
          write(iw,'(a)') trim(cvinfo%fcv(ifile))
        end do
      end if

      cvinfo%nfile = nfile

    end subroutine read_ctrl_cvinfo
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine read_ctrl_cvinfo_weight(iunit, cvinfo)
!-----------------------------------------------------------------------
      implicit none
!
      integer,        intent(in)    :: iunit
      type(s_cvinfo), intent(inout) :: cvinfo

      integer :: ifile, nfile

      
      if (.not. allocated(cvinfo%fcv)) then
        write(iw,'("Read_Ctrl_Cvinfo_Weight> Error.")')
        write(iw,'("This subroutine should be called &
                   &after calling read_ctrl_cvinfo.")')
        stop
      end if

      nfile = cvinfo%nfile
      allocate(cvinfo%fweight(nfile))

      do ifile = 1, cvinfo%nfile
        read(iunit,'(a)') cvinfo%fweight(ifile)
      end do

      write(iw,*)
      write(iw,'(">> CVINFO Weight parameter parameters")')
      write(iw,'("nfile    = ", i0)')  nfile
      if (ifile >= 20) then
        do ifile = 1, 10
          write(iw,'(a)') trim(cvinfo%fweight(ifile))
        end do
        write(iw,'(" ... skipped ...")')
        do ifile = nfile - 10, nfile
          write(iw,'(a)') trim(cvinfo%fweight(ifile))
        end do
      else
        do ifile = 1, nfile
          write(iw,'(a)') trim(cvinfo%fweight(ifile))
        end do
      end if


    end subroutine read_ctrl_cvinfo_weight
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine read_cv(fcv, ndim, cv, plus, nsta, nstep_read)
!-----------------------------------------------------------------------
      implicit none

      character(len=MaxChar), intent(in)  :: fcv
      integer,                intent(in)  :: ndim
      type(s_cv),             intent(out) :: cv
      integer, optional,      intent(in)  :: plus
      integer, optional,      intent(in)  :: nsta
      integer, optional,      intent(in)  :: nstep_read

      integer                :: istep, jstep, icv
      integer                :: nstep, np, ns, ncomm
      character(len=MaxChar) :: cdum, cdum2


      np = 0
      if (present(plus)) then
        np = plus
      end if

      ns = 0
      if (present(nsta)) then
        ns = nsta - 1
      end if

      open(UnitCV,file=trim(fcv))
      nstep = 0
      ncomm = 0
      do while(.true.)
        read(UnitCV,'(a)',end=100) cdum
        cdum2 = trim(adjustl(cdum))
        if (cdum2(1:1) /= "#") then
          nstep = nstep + 1
        else
          ncomm = ncomm + 1
        end if
      end do

100   rewind UnitCV

      cv%nstep = nstep -  ns
      if (present(nstep_read)) then
        if (nstep_read /= 0) then
          if (nstep_read > cv%nstep) then
            write(iw,'("Read_Cv> Warning: nstep_read > nstep in CV file")')
          else
            cv%nstep = nstep_read
          end if
        end if
      end if

      allocate(cv%x(cv%nstep))
      allocate(cv%data(ndim + np, cv%nstep))

      cv%x  = '' 
      jstep = 0
      !do istep = 1, nstep
      do istep = 1, ncomm + ns + cv%nstep
        if (istep <= ncomm + ns) then
          read(UnitCV,*) cdum
        else
          jstep = jstep + 1
          read(UnitCV,*)                              &
            cv%x(jstep),                              &
            (cv%data(icv, jstep), icv = 1, ndim + np)
        end if
      end do
      close(UnitCV)

      cv%ndim = ndim
      

    end subroutine read_cv
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine read_cv_weight(fweight, cv, nsta, nstep_read)
!-----------------------------------------------------------------------
      implicit none

      character(len=MaxChar), intent(in)    :: fweight
      type(s_cv),             intent(inout) :: cv
      integer, optional,      intent(in)    :: nsta
      integer, optional,      intent(in)    :: nstep_read

      ! I/O
      !
      integer :: io

      ! Dummy
      !
      character(len=MaxChar) :: cdum, cdum2, cdum3 

      integer                :: istep, jstep, icv
      integer                :: nstep, np, ns, nstep_chk
      real(8)                :: w


      ns = 0
      if (present(nsta)) then
        ns = nsta - 1
      end if

      call open_file(fweight, io, stat = 'old')

      nstep = 0
      if (.not. allocated(cv%data)) then

        do while (.true.)

          read(io,'(a)', end = 100) cdum
          cdum2 = trim(adjustl(cdum))

          if (cdum2(1:1) /= '#') then
            nstep = nstep + 1 
          end if

        end do

100     rewind io

        cv%nstep = nstep - ns
        if (present(nstep_read)) then
          if (nstep_read /= 0) then
            if (nstep_read > cv%nstep) then
              write(iw,'("Read_Cv> Warning: nstep_read > nstep in CV file")')
            else
              cv%nstep = nstep_read
            end if
          end if
        end if

      end if

      allocate(cv%weight(cv%nstep))

      istep = 0
      jstep = 0
      do while (.true.)
        read(io,'(a)', end = 101) cdum

        cdum2 = trim(adjustl(cdum))

        if (cdum2(1:1) /= '#') then
          istep = istep + 1
          if (istep > ns .and. istep <= ns + cv%nstep) then
            jstep = jstep + 1
            read(cdum, *) cdum3, cv%weight(jstep)  
          end if
        end if

      end do

101   close(io)


    end subroutine read_cv_weight
!-----------------------------------------------------------------------
end module
!=======================================================================
