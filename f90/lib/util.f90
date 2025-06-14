!=======================================================================
module mod_util
!=======================================================================
  implicit none

  integer, parameter, private :: MaxChar = 10000
  integer, parameter, private :: iw      = 6

  public :: open_file
  public :: is_specified
  public :: tolower
  public :: toupper
  public :: get_tof
  public :: get_opt
  public :: get_atomicnum
  public :: get_file_extention
  public :: gen_file_name
  public :: seek_line
  public :: check_input_parameter_real8
  public :: check_input_parameter_integer
  public :: trape_integral
  public :: standard_deviation

  contains
!-----------------------------------------------------------------------
    subroutine open_file(fname, iunit, frmt, stat)
!-----------------------------------------------------------------------
      implicit none

      character(*),           intent(in)  :: fname
      integer,                intent(out) :: iunit
      character(*), optional, intent(in)  :: frmt
      character(*), optional, intent(in)  :: stat 

      integer :: i
      logical :: chk

      iunit = 10
      do i = 1, 99
        inquire(unit=iunit, opened=chk)
        if (chk) then
          iunit = iunit + 1
        else

          if (present(frmt)) then
            if (present(stat)) then
              open(iunit, file=trim(fname), form=trim(frmt), status=trim(stat))
            else
              open(iunit, file=trim(fname), form=trim(frmt))
            end if
          else
            if (present(stat)) then
              open(iunit, file=trim(fname), status=trim(stat))
            else
              open(iunit, file=trim(fname))
            end if
          end if

          exit
        end if
      end do

    end subroutine open_file
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    function is_specified(string)
!-----------------------------------------------------------------------
      implicit none

      character(*), intent(in) :: string

      logical :: is_specified


      is_specified = .false.
      if (trim(string) /= "") then
        is_specified = .true.
      end if

    end function is_specified
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    function tolower(string)
!-----------------------------------------------------------------------
      implicit none

      character(*), intent(in) :: string 
      character(len=MaxChar)   :: tolower 

      integer :: i
      integer :: lenstr

      tolower = "          "
      lenstr  = len(string)
      do i = 1, lenstr
        if (string(i:i) >= 'A' .and. string(i:i) <= 'Z') then
          tolower(i:i) = char(ichar(string(i:i)) + 32)
        else
          tolower(i:i) = string(i:i)
        end if
      end do

    end function tolower 
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    function toupper(string)
!-----------------------------------------------------------------------
      implicit none

      character(*), intent(in) :: string 
      character(len=MaxChar)   :: toupper

      integer :: i
      integer :: lenstr

      toupper = "          "
      lenstr = len(string)
      do i = 1, lenstr
        if (string(i:i) >= 'a' .and. string(i:i) <= 'z') then
          toupper(i:i) = char(ichar(string(i:i)) - 32)
        else
          toupper(i:i) = string(i:i)
        end if
      end do

    end function toupper 
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    function get_tof(tof)
!-----------------------------------------------------------------------
      implicit none

      logical, intent(in) :: tof
      character(len=7)    :: get_tof


      if (tof) then
        get_tof = ".true. " 
      else
        get_tof = ".false."
      end if

    end function get_tof
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    function get_opt(opt, optlist, ierr)
!-----------------------------------------------------------------------
      implicit none
      
      character(*), intent(in)  :: opt
      character(*), intent(in)  :: optlist(:)
      integer,      intent(out) :: ierr

      integer                   :: get_opt

      integer                   :: iopt, nopt
      character(len=MaxChar)    :: optconv
      logical                   :: match


      optconv = toupper(opt) 

      nopt    = size(optlist)
      get_opt = 0
      match   = .false.
      ierr    = 0 
      do iopt = 1, nopt
        if (trim(adjustl(optconv)) == trim(adjustl(optlist(iopt)))) then
          get_opt = iopt
          match   = .true. 
        end if 
      end do

      if (.not. match) then
        ierr = 1
      end if

    end function get_opt 
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    function get_atomicnum(mass, ierr)
!-----------------------------------------------------------------------
      implicit none

      character(*), parameter   :: alist(11)  = (/'H ','C ','N ','O ', &
                                                  'NA','MG','P ','S ', &
                                                  'CL','K ','CA'/)
      integer,      parameter   :: anlist(11) = (/  1,   6,   7,   8, &
                                                  11,  12,  15,  16, &
                                                  17,  19,  20/)
      real(8),      parameter   :: mlist(11)  = (/ 1.00798d0, 12.0106d0, 14.0069d0, 15.994d0,&
                                                  22.9898d0,  24.306d0,  30.9738d0, 32.068d0,&
                                                  35.452d0,   39.0983d0, 40.078d0/) 

      real(8),      intent(in)  :: mass 
      integer,      intent(out) :: ierr

      integer                   :: get_atomicnum

      integer                   :: iopt, nopt
      real(8)                   :: dev
      logical                   :: match


      nopt     = size(alist)
      match    = .false.
      ierr     = 0

      get_atomicnum = 0
      do iopt = 1, nopt
        dev = abs(mass - mlist(iopt)) / mlist(iopt) * 100.0d0
        if (dev < 1.0d0) then 
          get_atomicnum = anlist(iopt) 
          match         = .true. 
        end if 
      end do

      if (.not. match) then
        ierr = 1
      end if

    end function get_atomicnum
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_file_extention(fname, ext, ierr)
!-----------------------------------------------------------------------
      implicit none

      character(*), intent(in)  :: fname
      character(*), intent(out) :: ext
      integer,      intent(out) :: ierr

      integer :: dotloc, length


      ierr = 0

      dotloc = index(fname, ".", back = .true.)
      if (dotloc == 0) then
        ierr = 1
      else
        length = len_trim(fname)
        ext    = fname(dotloc+1:length)
      end if
      

    end subroutine get_file_extention
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine gen_file_name(fhead, ext, fname, fspecific)
!-----------------------------------------------------------------------
      implicit none

      character(*),           intent(in)  :: fhead
      character(*),           intent(in)  :: ext 
      character(*),           intent(out) :: fname
      character(*), optional, intent(in)  :: fspecific


      if (present(fspecific)) then
        if (trim(fspecific) /= "") then
          write(fname,'(a)') trim(fspecific)
        else
          write(fname,'(a,".",a)') trim(fhead), trim(ext)
        end if
      else
        write(fname,'(a,".",a)') trim(fhead), trim(ext)
      end if


    end subroutine gen_file_name
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine seek_line(iunit, keyword, ierr)
!-----------------------------------------------------------------------
      implicit none

      integer,      intent(in)    :: iunit
      character(*), intent(in)    :: keyword
      integer,      intent(inout) :: ierr

      character(len=MaxChar) :: line


      rewind iunit

      ierr = -1 
      do while(.true.)
        read(iunit,'(a)',end=100) line
        if (trim(keyword) == trim(adjustl(line))) then
          ierr = 0
          exit
        end if
      end do

100   continue

    end subroutine seek_line
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine termination(program_name) 
!-----------------------------------------------------------------------
      implicit none

      character(*), intent(in) :: program_name


      write(6,'("============================================================")')
      write(6,'(a," terminated normally")') trim(program_name)
      stop


    end subroutine termination 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine check_input_parameter_real8(val,           &
                                           criteria,      &
                                           should_updown, &
                                           valname,       &
                                           abort)
!-----------------------------------------------------------------------
      implicit none

      real(8),      intent(in) :: val
      real(8),      intent(in) :: criteria
      character(*), intent(in) :: should_updown
      character(*), intent(in) :: valname
      logical,      intent(in) :: abort


      if (trim(should_updown) == "up") then
        if (val < criteria) then
          write(iw,*)
          write(iw,'("Check_Input_Parameter_Real8> ")')
          write(iw,'(a, " is smaller than ", f15.7)') trim(valname), criteria

          if (abort) then
            stop
          end if

        end if
      else if (trim(should_updown) == "down") then
        if (val > criteria) then
          write(iw,*)
          write(iw,'("Check_Input_Parameter_Real8> ")')
          write(iw,'(a, " is larger than ", f15.7)')  trim(valname), criteria

          if (abort) then
            stop
          end if

        end if
      end if

    end subroutine check_input_parameter_real8
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine check_input_parameter_integer(val,           &
                                             criteria,      &
                                             should_updown, &
                                             valname,       &
                                             abort)
!-----------------------------------------------------------------------
      implicit none

      integer,      intent(in) :: val
      integer,      intent(in) :: criteria
      character(*), intent(in) :: should_updown
      character(*), intent(in) :: valname
      logical,      intent(in) :: abort


      if (trim(should_updown) == "up") then
        if (val < criteria) then
          write(iw,*)
          write(iw,'("Check_Input_Parameter_Integer> ")')
          write(iw,'(a, " is smaller than ", i0)') trim(valname), criteria

          if (abort) then
            stop
          end if

        end if
      else if (trim(should_updown) == "down") then
        if (val > criteria) then
          write(iw,*)
          write(iw,'("Check_Input_Parameter_Integer> ")')
          write(iw,'(a, " is larger than ", f15.7)')  trim(valname), criteria

          if (abort) then
            stop
          end if

        end if
      end if

    end subroutine check_input_parameter_integer
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine trape_integral(nsta, nend, dx, func, val, cumm)
!-----------------------------------------------------------------------
      implicit none 

      integer,           intent(in)  :: nsta
      integer,           intent(in)  :: nend
      real(8),           intent(in)  :: dx
      real(8),           intent(in)  :: func(nsta:nend)
      real(8),           intent(out) :: val
      real(8), optional, intent(out) :: cumm(nsta:nend)  
     
      integer :: i
      real(8) :: f, p1, p2
     
      f = 0.5d0 * dx

      if (present(cumm)) then
        cumm = 0.0d0
      end if

      val = 0.0d0
      do i = nsta, nend - 1
        p1  = func(i)
        p2  = func(i+1)
        val = val + f * (p1 + p2)

        if (present(cumm)) then
          cumm(i) = val
        end if
      end do
      val        = val + f * p2

      if (present(cumm)) then
        cumm(nend) = val
      end if


    end subroutine trape_integral
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine standard_deviation(nsample, ave, flc, stdev)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in)  :: nsample
      real(8), intent(in)  :: ave
      real(8), intent(in)  :: flc(1:nsample) 
      real(8), intent(out) :: stdev

      integer :: i
      real(8) :: val, dev


      dev = 0.0d0
      do i = 1, nsample
        val = flc(i) - ave
        dev = dev + val * val
      end do
      stdev = sqrt(dev / dble(nsample - 1))

    end subroutine standard_deviation
!-----------------------------------------------------------------------


end module mod_util
!=======================================================================

