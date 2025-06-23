!===============================================================================
module mod_parse_arguments
!===============================================================================
  use mod_const
  use mod_util

  implicit none

  ! constants
  !
  integer, parameter :: MaxArgs = 100000

  ! structures
  !
  type :: s_parg
    integer                             :: num_args      = 0
    integer                             :: num_opts      = 0
    integer                             :: num_opts_used = 0
    character(len=MaxChar), allocatable :: args(:)
    integer,                allocatable :: opt_pos(:)
    integer,                allocatable :: opt_len(:)
    integer,                allocatable :: opt_ind(:)
    integer,                allocatable :: opt_pos_align(:)
  end type s_parg

  ! subroutines
  !
  public :: get_arguments
  public :: assign_option_names

  contains
!-------------------------------------------------------------------------------
    subroutine get_arguments(parg)
!-------------------------------------------------------------------------------
      implicit none

      type(s_parg), intent(inout) :: parg

      ! function
      integer :: iargc

      integer :: i, j
      integer :: narg, nopts, nlen


      ! get number of arguments
      !
      parg%num_args = iargc()
      narg          = parg%num_args

      allocate(parg%args(narg))

      ! get all the arguments
      !
      parg%args = ""
      do i = 1, narg
        call getarg(i, parg%args(i))
      end do

      ! get number of options
      !
      nopts = 0
      do i = 1, narg
        if (parg%args(i)(1:1) == "-") then
          nopts = nopts + 1
        end if
      end do
      parg%num_opts = nopts

      allocate(parg%opt_pos(0:nopts), parg%opt_len(nopts))

      ! get option positions (opt_pos) 
      !
      j            = 0
      parg%opt_pos = 0
      do i = 1, narg
        if (parg%args(i)(1:1) == "-") then
          j = j + 1
          parg%opt_pos(j) = i
        end if
      end do

      ! get option length (opt_len)
      !
      parg%opt_len = 0
      do i = 1, nopts - 1
        nlen            = parg%opt_pos(i + 1) - parg%opt_pos(i) - 1
        parg%opt_len(i) = nlen 
      end do

      nlen                = narg - parg%opt_pos(nopts) 
      parg%opt_len(nopts) = nlen


    end subroutine get_arguments
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
    subroutine assign_option_names(indices, names, parg) 
!-------------------------------------------------------------------------------
      implicit none

      integer,         intent(in)    :: indices(:)
      character(*),    intent(in)    :: names(:)
      type(s_parg),    intent(inout) :: parg

      integer                :: i, j, ipos
      integer                :: nopts         = 0
      integer                :: nopts_in_args = 0
      character(len=MaxChar) :: opt, opt2 


      nopts         = size(names)
      nopts_in_args = parg%num_opts


      allocate(parg%opt_pos_align(nopts), parg%opt_ind(nopts))

      parg%num_opts_used = nopts
      parg%opt_pos_align = 0
      parg%opt_ind       = 0
      do i = 1, nopts
        opt = toupper(names(i))
        do j = 1, nopts_in_args
          ipos = parg%opt_pos(j)
          opt2 = toupper(parg%args(ipos))
          if (trim(opt) == trim(opt2)) then
            parg%opt_pos_align(i) = ipos
            parg%opt_ind(i)       = j
          end if
        end do
      end do

    end subroutine assign_option_names
!-------------------------------------------------------------------------------

end module mod_parse_arguments
!===============================================================================
