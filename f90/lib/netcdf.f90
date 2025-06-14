!=======================================================================
module mod_netcdfio
!=======================================================================
  use mod_const
  use mod_util
  use netcdf
  implicit none

  ! structures
  !
  type :: s_netcdf
    integer :: nstep
    integer :: natm
    integer :: coord_id
    integer :: box_id

    real(8), allocatable :: coord(:, :, :)
    real(8), allocatable :: box  (:, :)
  end type s_netcdf


  ! subroutines
  !
  public :: netcdf_open
  public :: netcdf_close
  public :: netcdf_read_dimension
  public :: netcdf_read_oneframe
  public :: get_total_step_from_netcdf

  contains
!-----------------------------------------------------------------------
    subroutine netcdf_open(netcdffile, iunit)
!-----------------------------------------------------------------------
      implicit none

      character(len=MaxChar), intent(in)  :: netcdffile
      integer,                intent(out) :: iunit

      integer :: retval 


      retval = nf90_open(netcdffile, NF90_NOWRITE, iunit)

      if (retval /= nf90_noerr) then
        write(iw,'("NetCDF_Open> Error.")')
        write(iw,'("Failed in opening ",a,". Stop.")') trim(netcdffile)
        stop
      end if

    end subroutine netcdf_open
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine netcdf_close(iunit, netcdffile)
!-----------------------------------------------------------------------
      implicit none

      integer,                intent(in)            :: iunit
      character(len=MaxChar), optional, intent(in)  :: netcdffile 


      close(iunit)
!
    end subroutine netcdf_close
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine netcdf_read_dimension(io, nc, nstore)
!-----------------------------------------------------------------------
      implicit none

      integer,        intent(in)    :: io
      type(s_netcdf), intent(inout) :: nc
      
      integer, optional, intent(in) :: nstore

      integer :: retval, dimid_frame, dimid_atom
      integer :: nst


      nst = 1 

      ! # of frames
      !
      retval = nf90_inq_dimid        (io, "frame", dimid_frame)
      retval = nf90_inquire_dimension(io, dimid_frame, len=nc%nstep)

      ! # of atoms
      !
      retval = nf90_inq_dimid        (io, "atom", dimid_atom) 
      retval = nf90_inquire_dimension(io, dimid_atom, len=nc%natm)

      ! Coord & Box IDs 
      !
      retval = nf90_inq_varid        (io, "coordinates" , nc%coord_id)
      retval = nf90_inq_varid        (io, "cell_lengths", nc%box_id)
      
      ! Memory allocation  
      !
      if (.not. allocated(nc%coord)) then

        if (present(nstore)) then
          if (nstore == -1) then
            nst = nc%nstep
          else
            nst = nstore 
          end if
        else
          nst = 1
        end if

        allocate(nc%coord(1:3, nc%natm, nst))
        allocate(nc%box(1:3, nst))

      end if
       
!
    end subroutine netcdf_read_dimension
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine netcdf_read_oneframe(io, istep, nc)
!-----------------------------------------------------------------------
      implicit none

      integer,           intent(in)    :: io
      integer,           intent(in)    :: istep
      type(s_netcdf),    intent(inout) :: nc 

      integer :: rstval
     

      rstval = nf90_get_var(                               &
                 io,                                       &
                 nc%coord_id, nc%coord(1:3, 1:nc%natm, 1), &
                 start = (/1,       1, istep/),            &
                 count = (/3, nc%natm,     1/))

      rstval = nf90_get_var(                               &
                 io,                                       &
                 nc%box_id,  nc%box(1:3, 1),               &
                 start = (/1, istep/),                     &
                 count = (/3,     1/))

!
    end subroutine netcdf_read_oneframe 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_total_step_from_netcdf(fset, nstep_tot)
!-----------------------------------------------------------------------
      implicit none

      character(len=MaxChar), intent(in)  :: fset(:)
      integer,                intent(out) :: nstep_tot

      type(s_netcdf)         :: nc
      integer                :: itraj, ntraj
      integer                :: iunit
      character(len=MaxChar) :: f


      ntraj = size(fset(:))

      nstep_tot = 0
      do itraj = 1, ntraj
        f = fset(itraj)

        if (trim(f) /= "") then
          call netcdf_open(fset(itraj), iunit)
          call netcdf_read_dimension(iunit, nc)
          nstep_tot = nstep_tot + nc%nstep
          call netcdf_close(iunit)
        end if
      end do

!
    end subroutine get_total_step_from_netcdf
!-----------------------------------------------------------------------

end module mod_netcdfio
!=======================================================================
