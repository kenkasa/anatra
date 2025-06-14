!=======================================================================
module mod_xtcio
!=======================================================================
  use mod_const
  use xdr, only : xtcfile

  implicit none

  ! structures
  !


  ! subroutines
  !
  public  :: get_natm_from_xtc
  public  :: get_total_step_from_xtc

  contains

!-----------------------------------------------------------------------
    subroutine get_natm_from_xtc(xtcfilename, natm)
!-----------------------------------------------------------------------
      implicit none

      character(len=MaxChar), intent(in)  :: xtcfilename
      integer,                intent(out) :: natm

      integer                             :: i
      integer                             :: iunit

      type(xtcfile) :: xtc


      iunit = UnitXTC

      write(iw,*)
      write(iw,'("Get_Natm_From_XTC> Reading XTC file ", a)') trim(xtcfilename)

      call xtc%init(trim(xtcfilename))
      call xtc%read
      
      natm = xtc%natoms

      call xtc%close

!
    end subroutine get_natm_from_xtc
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine get_total_step_from_xtc(fset, nstep_tot)
!-----------------------------------------------------------------------
      implicit none

      character(len=MaxChar), intent(in)  :: fset(:)
      integer,                intent(out) :: nstep_tot

      type(xtcfile)          :: xtc
      integer                :: itraj, ntraj
      integer                :: iunit
      logical                :: is_end
      character(len=MaxChar) :: f


      ntraj     = size(fset(:))
      nstep_tot = 0
      do itraj = 1, ntraj
        f = fset(itraj)

        if (trim(f) /= "") then
          call xtc%init(trim(fset(itraj)))
          is_end = .false.

          do while(.not. is_end)
            call xtc%read
            if (xtc%STAT /= 0) then
              is_end = .true.
              exit
            end if

            nstep_tot = nstep_tot + 1

            if (mod(nstep_tot, 100) == 0) then
              write(iw,'("Read ",i0)') nstep_tot
            end if

            if (xtc%STAT /= 0) then
              is_end = .true.
            end if
          end do

          call xtc%close

        end if
      end do

    end subroutine get_total_step_from_xtc
!-----------------------------------------------------------------------
!
end module mod_xtcio
!=======================================================================
