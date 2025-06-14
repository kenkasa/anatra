!==============================================================================
module mod_random
!==============================================================================
  use mod_const
  use mtmod
  implicit none
 
  ! constants
  !

  ! structures 
  !
  
  !subroutines
  public :: get_seed
  public :: initialize_random
  public :: get_random
  public :: get_random_integer  

  contains
!------------------------------------------------------------------------------
    subroutine get_seed(iseed)
!------------------------------------------------------------------------------
      implicit none

      integer, intent(inout) :: iseed

      integer :: iclock 


      if (iseed <= 0) then
        call system_clock(count=iclock)
        iseed = iclock
      end if

      write(iw,*)
      write(iw,'("Get_Seed> generatege random seed : ",i0)') iseed


    end subroutine get_seed
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
    subroutine initialize_random(iseed)
!------------------------------------------------------------------------------
      implicit none
!
      integer, intent(in) :: iseed
!
      call sgrnd(iseed)

    end subroutine initialize_random
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
    subroutine get_random(rand)
!------------------------------------------------------------------------------
      implicit none

      real(8), intent(out) :: rand


      rand = grnd()

    end subroutine get_random
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
    subroutine get_random_integer(nsample, rand_min, rand_max, duplicate, rand)
!------------------------------------------------------------------------------
      implicit none

      integer, intent(in)  :: nsample
      integer, intent(in)  :: rand_min
      integer, intent(in)  :: rand_max
      logical, intent(in)  :: duplicate
      integer, intent(out) :: rand(nsample) 

      integer :: i, j
      integer :: curr, prev, nget
      logical :: is_finished, is_duplicated


      nget = 0

      is_finished = .false.
      rand        = 0
      prev        = 0

      do while (.not. is_finished)
        curr = nint(grnd() * (rand_max - rand_min))
        curr = curr + rand_min
        if (curr < rand_min) then
          curr = rand_min 
        end if

        if (curr > rand_max) then
          curr = rand_max
        end if

        if (duplicate) then
          nget       = nget + 1
          rand(nget) = curr
        else
          !if (nget /= 0) then
            is_duplicated = .false.
            do i = 1, nget
              prev = rand(i)
              if (curr == prev) then
                is_duplicated = .true.    
              end if 
            end do

            if (.not. is_duplicated) then
              nget       = nget + 1
              rand(nget) = curr 
            end if 
          !else
          !  nget       = nget + 1
          !  rand(nget) = curr
          !end if
        end if
        
        if (nget == nsample) then 
          is_finished = .true.
        end if
      end do


    end subroutine get_random_integer
!------------------------------------------------------------------------------

end module mod_random
!==============================================================================
