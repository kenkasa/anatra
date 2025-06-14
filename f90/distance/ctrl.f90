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
  integer,      parameter, public :: DistanceTypeSTANDARD = 1 
  integer,      parameter, public :: DistanceTypeMINIMUM  = 2
  integer,      parameter, public :: DistanceTypeINTRA    = 3
  character(*), parameter, public :: DistanceTypes(3)     = (/'STANDARD',&
                                                              'MINIMUM ',&
                                                              'INTRA   '/)
   
  integer,      parameter, public :: MinDistTypeSITE = 1
  integer,      parameter, public :: MinDistTypeCOM  = 2
  character(*), parameter, public :: MinDistTypes(2) = (/'SITE', 'COM '/)


  ! structures
  !
  type :: s_option
    logical :: pbc             = .false.
    integer :: mode(2)         = (/ComModeRESIDUE, ComModeRESIDUE/)
    integer :: distance_type   = DistanceTypeSTANDARD
    integer :: mindist_type(2) = (/MinDistTypeSITE, MinDistTypeSITE/)
    real(8) :: t_sta           = 0.0d0
    real(8) :: t_end           = 0.0d0
  end type s_option

  ! subroutines
  !
  public  :: read_ctrl
  private :: read_ctrl_option

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
      integer                 :: io
      character(len=MaxChar)  :: f_ctrl

      ! Get control file name
      !
      call getarg(1, f_ctrl)

      write(iw,*)
      write(iw,'("Read_Ctrl> Reading parameters from ", a)') trim(f_ctrl)
      call open_file(f_ctrl, io)
      call read_ctrl_input  (io, input)
      call read_ctrl_output (io, output)
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

      logical                :: pbc             = .false.
      character(len=MaxChar) :: mode(2)         = (/"RESIDUE", "RESIDUE"/)
      character(len=MaxChar) :: distance_type   = "STANDARD"
      character(len=MaxChar) :: mindist_type(2) = (/"SITE", "SITE"/)
      real(8)                :: t_sta
      real(8)                :: t_end

      ! Parser
      !
      integer                :: iopt, ierr

      ! Dummy
      !
      integer                :: i

      namelist /option_param/ &
        pbc,                  &
        mode,                 &
        distance_type,        &
        mindist_type,         &
        t_sta,                &
        t_end 


      ! Initialize
      !
      pbc           = .false.
      mode(1:2)     = "RESIDUE"
      distance_type = "STANDARD"
      mindist_type  = "SITE"
      t_sta         = 0.0d0
      t_end         = 0.0d0


      ! Read namelist
      !
      rewind io
      read(io, option_param)

      ! Show parameters
      !
      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("pbc           = ", a)') get_tof(pbc)
      write(iw,'("mode          = ",2(a,2x))') (trim(mode(i)), i = 1, 2)
      write(iw,'("distance_type = ", a)') trim(distance_type)
      write(iw,'("t_sta         = ", f15.7)') t_sta
      write(iw,'("t_end         = ", f15.7)') t_end

     
      if (t_sta <= 1.0d-8) then
        t_sta = -1.0d20
      end if

      if (t_end <= 1.0d-8) then
        write(iw,*)
        write(iw,'("Remark: Detect t_end = 0.0")')
        write(iw,'(">> input trajectories are read till the end of the records")')
        write(iw,'("(This is not error)")')
        t_end = 1.0d20
      end if

      do i = 1, 2
        iopt = get_opt(mode(i), CoMMode, ierr)
        if (ierr /= 0) then
          write(iw,'("Read_Ctrl_Option> Error.")')
          write(iw,'("mode = ",a," is not available.")') trim(mode(i))
          stop
        end if
        option%mode(i)     = iopt
      end do

      iopt = get_opt(distance_type, DistanceTypes, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("distance_type = ",a," is not available.")') trim(distance_type)
        stop
      end if
      option%distance_type = iopt

      do i = 1, 2
        iopt = get_opt(mindist_type(i), MinDistTypes, ierr)
        if (ierr /= 0) then
          write(iw,'("Read_Ctrl_Option> Error.")')
          write(iw,'("mode = ",a," is not available.")') trim(mode(i))
          stop
        end if
        option%mindist_type(i) = iopt
      end do
                          
      option%pbc = pbc

      if (option%distance_type == DistanceTypeMINIMUM) then
        write(iw,'("mindist_type  = ",2(a,2x))') &
          (trim(mindist_type(i)), i = 1, 2)
      end if

      option%t_sta = t_sta
      option%t_end = t_end

      ! Combination check
      !


    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end module
!=======================================================================
