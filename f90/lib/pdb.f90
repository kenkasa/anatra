!=======================================================================
module mod_pdb
!=======================================================================
  implicit none

  ! structures
  !
  type :: s_pdb
    integer                       :: natm
    integer,          allocatable :: resid(:)
    integer,          allocatable :: atomid(:)
    character(len=5), allocatable :: atomname(:)
    character(len=5), allocatable :: resname(:)
    real(8),          allocatable :: coord(:, :)
    real(8),          allocatable :: box(:)
    character(len=4), allocatable :: segid(:)
  end type s_pdb

  ! subroutines
  !
  public :: read_pdb
  public :: write_pdb
  public :: write_pdb_oneline
  public :: alloc_pdb
  public :: dealloc_pdb

  contains

!-----------------------------------------------------------------------
    subroutine read_pdb(iunit, pdb)
!-----------------------------------------------------------------------
      implicit none

      integer,     intent(in)  :: iunit
      type(s_pdb), intent(out) :: pdb

      integer           :: iatom, ixyz
      character(len=80) :: line


      rewind iunit

      ! get number of atoms
      !
      iatom = 0
      do while (.true.)
        read(iunit,'(a)', end=100) line
        if (line(1:4) == "ATOM") then
          iatom = iatom + 1
        end if
      end do

100   continue
      rewind iunit

      pdb%natm = iatom
      call alloc_pdb(pdb%natm, pdb)

      ! read pdb information
      !
      iatom = 0
      do while (.true.)
        read(iunit,'(a)',end=101) line
        if (line(1:4) == "ATOM") then
          iatom = iatom + 1
          pdb%atomid(iatom) = iatom
          read(line(13:16),*) pdb%atomname(iatom)
          read(line(18:21),*) pdb%resname(iatom)
          read(line(23:26),*) pdb%resid(iatom)
          read(line(31:38),'(f8.3)') pdb%coord(1, iatom) 
          read(line(39:46),'(f8.3)') pdb%coord(2, iatom)  
          read(line(47:54),'(f8.3)') pdb%coord(3, iatom)  
        end if
      end do

101   continue


    end subroutine read_pdb
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine write_pdb(iunit, pdb)
!-----------------------------------------------------------------------
      implicit none

      integer,     intent(in) :: iunit
      type(s_pdb), intent(in) :: pdb

      integer          :: iatom, ixyz
      character(len=5) :: iatom_c


      do iatom = 1, pdb%natm
        if (iatom < 99999) then
          write(iatom_c,'(i5)') iatom
        else
          write(iatom_c,'(a5)') '*****'
        end if
        write(iunit,'(a6,a5,1x,a4,1x,a4,1x,i4,4x,3f8.3,2f6.2,6x,a4)') &
          "ATOM  ",                              &
          iatom_c,                               &
          pdb%atomname(iatom),                   &
          pdb%resname(iatom),                    &
          pdb%resid(iatom),                      &
          (pdb%coord(ixyz, iatom), ixyz = 1, 3), &
          0.0d0,                                 &
          0.0d0,                                 &
          pdb%segid(iatom) 
      end do

    end subroutine write_pdb 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine write_pdb_oneline(iunit,    &
                                 iatm,     &
                                 atomname, &
                                 resname,  &
                                 resid,    &
                                 coord)
!-----------------------------------------------------------------------
      implicit none

      integer,      intent(in) :: iunit
      integer,      intent(in) :: iatm
      character(*), intent(in) :: atomname
      character(*), intent(in) :: resname
      integer,      intent(in) :: resid
      real(8),      intent(in) :: coord(3)

      character(len=5) :: iatmc


      if (iatm < 99999) then
        write(iatmc, '(i5)') iatm
      else
        write(iatmc, '(a5)') '*****'
      end if

      write(iunit,'(a6)',    advance = 'no') "ATOM  " 
      write(iunit,'(a5)',    advance = 'no') iatmc
      write(iunit,'(1x)',    advance = 'no')
      write(iunit,'(a4)',    advance = 'no') trim(atomname)
      write(iunit,'(1x)',    advance = 'no')  
      write(iunit,'(a4)',    advance = 'no') trim(resname)
      write(iunit,'(1x)',    advance = 'no')  
      write(iunit,'(i4)',    advance = 'no') resid
      write(iunit,'(4x)',    advance = 'no') 
      write(iunit,'(3f8.3)', advance = 'no') coord(1:3)
      write(iunit,*)


    end subroutine write_pdb_oneline
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine alloc_pdb(natm, pdb)
!-----------------------------------------------------------------------
      implicit none

      integer,     intent(in)    :: natm
      type(s_pdb), intent(inout) :: pdb


      pdb%natm = natm
      allocate(pdb%resid(natm))
      allocate(pdb%atomid(natm))
      allocate(pdb%atomname(natm))
      allocate(pdb%resname(natm))
      allocate(pdb%coord(3, natm))
      allocate(pdb%box(3))
      allocate(pdb%segid(natm))

    end subroutine alloc_pdb
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine dealloc_pdb(pdb)
!-----------------------------------------------------------------------
      implicit none

      type(s_pdb), intent(inout) :: pdb


      deallocate(pdb%resid)
      deallocate(pdb%atomid)
      deallocate(pdb%atomname)
      deallocate(pdb%resname)
      deallocate(pdb%coord)
      deallocate(pdb%box)

    end subroutine dealloc_pdb
!-----------------------------------------------------------------------

end module mod_pdb
!=======================================================================
