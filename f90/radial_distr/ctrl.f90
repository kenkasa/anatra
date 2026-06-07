!=======================================================================
module mod_ctrl
!=======================================================================
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_traj
  use mod_com

  implicit none

  ! constants
  !
  integer,      parameter, public :: DensityTypeNUMBER   = 1
  integer,      parameter, public :: DensityTypeELECTRON = 2 
  character(*), parameter, public :: DensityType(2)      = (/'NUMBER   ',&
                                                             'ELECTRON '/)

  ! structures
  !
  type :: s_option
    real(8) :: dr            = 0.1d0
    integer :: mode(2)       = (/CoMModeRESIDUE, CoMModeRESIDUE/)
    logical :: identical     = .false.
    logical :: normalize     = .false.
    logical :: separate_self = .false.

    real(8) :: dt            = 1.0d0
    real(8) :: t_sta         = -1.0d0
    real(8) :: t_end         = -1.0d0 

    ! prepared after reading namelists
    !
    integer :: nt_sta        = 0
    integer :: nt_end        = 0
  end type s_option

  ! subroutines
  !
  public  :: read_ctrl
  private :: read_ctrl_option
  private :: show_input
  private :: show_output

  contains

!-----------------------------------------------------------------------
    subroutine read_ctrl(input, output, option, trajopt)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),   intent(out) :: input
      type(s_output),  intent(out) :: output
      type(s_option),  intent(out) :: option
      type(s_trajopt), intent(out) :: trajopt

      ! I/O
      !
      integer                :: io
      character(len=MaxChar) :: f_ctrl

      ! get control file name
      !
      call getarg(1, f_ctrl)

      write(iw,*)
      write(iw,'("Read_Ctrl> Reading parameters from ", a)') trim(f_ctrl)

      call open_file(f_ctrl, io, stat = 'old')

      call read_ctrl_input  (io, input)
      call show_input       (input)

      call read_ctrl_output (io, output)
      call show_output      (output)

      call read_ctrl_option (io, option)
      call read_ctrl_trajopt(io, trajopt)

      close(io)

    end subroutine read_ctrl
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine read_ctrl_option(io, option)
!-----------------------------------------------------------------------
      implicit none
!
      integer,        intent(in)  :: io
      type(s_option), intent(out) :: option 

      real(8)                :: dr            = 0.1d0
      character(len=MaxChar) :: mode(2)       = (/"RESIDUE", "RESIDUE"/)
      logical                :: identical     = .false.
      logical                :: normalize     = .false.
      logical                :: separate_self = .false.

      real(8)                :: dt            =  1.0d0
      real(8)                :: t_sta         = -1.0d0
      real(8)                :: t_end         = -1.0d0

      ! Parser
      !
      integer :: iopt, ierr

      ! Dummy
      !
      integer :: itraj

      namelist /option_param/ &
        dr,                   &
        mode,                 &
        identical,            &
        normalize,            &
        separate_self,        &
        dt,                   &
        t_sta,                &
        t_end


      rewind io
      read(io, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("dr            = ", f15.7)')   dr 
      write(iw,'("mode          = ", a,2x,a)')  trim(mode(1)), trim(mode(2))
      write(iw,'("identical     = ", a)')       get_tof(identical)
      write(iw,'("normalize     = ", a)')       get_tof(normalize)
      write(iw,'("separate_self = ", a)')       get_tof(separate_self)
      write(iw,'("dt            = ", f15.7)')   dt
      write(iw,'("t_sta         = ", f15.7)')   t_sta
      write(iw,'("t_end         = ", f15.7)')   t_end

      do itraj = 1, 2
        iopt = get_opt(mode(itraj), CoMMode, ierr)
        if (ierr /= 0) then
          write(iw,'("Read_Ctrl_Option> Error.")')
          write(iw,'("mode = ",a," is not available.")') trim(mode(itraj))
          stop
        end if
        option%mode(itraj) = iopt
      end do
                           
      option%dr            = dr
      option%identical     = identical
      option%normalize     = normalize
      option%separate_self = separate_self
      option%dt            = dt
      option%t_sta         = t_sta
      option%t_end         = t_end

      ! Combination check
      !

      ! Convert from real to integer
      !
      option%nt_sta = 0
      option%nt_end = 0
      if (option%t_sta >= -1e-5 .and. option%t_end > 0.0d0) then
        option%nt_sta = nint(option%t_sta / option%dt)
        option%nt_end = nint(option%t_end / option%dt)
        if (option%nt_sta == 0) option%nt_sta = 1
      end if

    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine show_input(input)
!-----------------------------------------------------------------------
      implicit none

      type(s_input), intent(in) :: input

      ! Dummy
      !
      integer :: i

      ! Check
      !
      if (input%ntraj == 0) then
        write(iw,'("Error. ftraj should be specified.")')
        stop
      end if

      ! Print
      !
      write(iw,*)
      write(iw,'(">> Input section parameters")')
      do i = 1, input%ntraj
        write(iw,'("ftraj", 3x, i0, 3x, " = ", a)') i, trim(input%ftraj(i))
      end do


    end subroutine show_input
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine show_output(output)
!-----------------------------------------------------------------------
      implicit none

      type(s_output), intent(in) :: output 


      ! Check
      !
      if (trim(output%fhead) == '') then
        write(iw,'("Error. fhead should be specified.")')
        stop
      end if

      ! Print
      !
      write(iw,*)
      write(iw,'(">> Output section parameters")')
      write(iw,'("fhead          = ", a)') trim(output%fhead)


    end subroutine show_output
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
end module
!=======================================================================
