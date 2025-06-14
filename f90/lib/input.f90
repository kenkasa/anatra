!=======================================================================
module mod_input
!=======================================================================
  use mod_util
  use mod_const

  implicit none

  ! parameters
  !

  ! structures
  !
  type :: s_input

    integer                :: ntraj
    character(len=MaxChar) :: ftraj(MaxTraj)
    character(len=MaxChar) :: flist_traj

    integer                :: ncv
    character(len=MaxChar) :: fcv(MaxTraj)
    character(len=MaxChar) :: flist_cv

    integer                :: nweight
    character(len=MaxChar) :: fweight(MaxTraj)
    character(len=MaxChar) :: flist_weight

    integer                :: nprof
    character(len=MaxChar) :: fprof(MaxTraj)
    character(len=MaxChar) :: flist_prof

    character(len=MaxChar) :: fdx

    character(len=MaxChar) :: fanaparm, fanaparm2  
    character(len=MaxChar) :: fprmtop,  fprmtop2

  end type s_input

  ! subroutines
  !
  public  :: read_ctrl_input
  private :: setup_filelist 

  contains
!-----------------------------------------------------------------------
    subroutine read_ctrl_input(iunit, input, myrank)
!-----------------------------------------------------------------------
      implicit none

      integer,           intent(in)  :: iunit
      type(s_input),     intent(out) :: input
      integer, optional, intent(in)  :: myrank

      ! Local
      !
      character(len=MaxChar) :: ftraj(1:MaxTraj)   = ''
      character(len=MaxChar) :: flist_traj         = ''

      character(len=MaxChar) :: fcv(1:MaxTraj)     = ''
      character(len=MaxChar) :: flist_cv           = ''

      character(len=MaxChar) :: fweight(1:MaxTraj) = ''
      character(len=MaxChar) :: flist_weight       = ''

      character(len=MaxChar) :: fprof(1:MaxTraj)   = ''
      character(len=MaxChar) :: flist_prof         = ''

      character(len=MaxChar) :: fdx                = ''

      character(len=MaxChar) :: fanaparm           = ''
      character(len=MaxChar) :: fanaparm2          = ''

      character(len=MaxChar) :: fprmtop            = ''
      character(len=MaxChar) :: fprmtop2           = ''

      integer :: ntraj, ncv, nweight, nprof
      integer :: irank

      ! I/O
      !
      integer :: io

      ! Dummy
      !
      integer :: i

      namelist /input_param/  &
        ftraj,                &
        flist_traj,           &
        fcv,                  &
        flist_cv,             &
        fweight,              &
        flist_weight,         &
        fprof,                &
        flist_prof,           &
        fdx,                  &
        fanaparm,             &
        fanaparm2,            &
        fprmtop,              &
        fprmtop2


      if (present(myrank)) then
        irank = myrank
      else
        irank = 0
      end if

      ! Initialize
      !
      ntraj        = 0
      ftraj        = ''
      flist_traj   = ''

      ncv          = 0
      fcv          = ''
      flist_cv     = ''

      nweight      = 0
      fweight      = ''
      flist_weight = ''

      nprof        = 0
      fprof        = ''
      flist_prof   = ''

      fdx          = ''

      ! Read INPUT_PARAM namelist
      !
      rewind iunit
      read(iunit,input_param)

      call setup_filelist(ftraj,   flist_traj,   ntraj)
      call setup_filelist(fcv,     flist_cv,     ncv)
      call setup_filelist(fweight, flist_weight, nweight)
      call setup_filelist(fprof,   flist_prof,   nprof)

      ! Send values in input parameters to input structure
      !
      input%ntraj        = ntraj
      input%ftraj        = ftraj
      input%flist_traj   = flist_traj

      input%ncv          = ncv
      input%fcv          = fcv
      input%flist_cv     = flist_cv

      input%nweight      = nweight
      input%fweight      = fweight
      input%flist_weight = flist_weight

      input%nprof        = nprof
      input%fprof        = fprof
      input%flist_prof   = flist_prof

      input%fdx          = fdx

      input%fprmtop      = fprmtop
      input%fprmtop2     = fprmtop2

      input%fanaparm     = fanaparm
      input%fanaparm2    = fanaparm2

    end subroutine read_ctrl_input
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine setup_filelist(fs, flist, n) 
!-----------------------------------------------------------------------
      implicit none

      character(len=MaxChar), intent(inout) :: fs(1:MaxTraj)
      character(len=MaxChar), intent(inout) :: flist
      integer,                intent(inout) :: n

      ! I/O
      !
      integer :: io

      ! Dummy 
      !
      integer :: i


      n = 0
      do i = 1, MaxTraj
        if (trim(fs(i)) /= '') then
          n = n + 1
        end if
      end do

      if (trim(flist) /= '') then

        n = 0
        call open_file(flist, io)
        
        do while (.true.)
          read(io, *, end = 100)
          n = n + 1
        end do

100     rewind io

        fs = ''
        do i = 1, n
          read(io, '(a)') fs(i)
        end do

        close(io)

      end if
      

    end subroutine setup_filelist
!-----------------------------------------------------------------------

end module mod_input
!=======================================================================
