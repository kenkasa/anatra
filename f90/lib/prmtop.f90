!=======================================================================
module mod_prmtop
!=======================================================================
  use mod_util
  use mod_const

  implicit none

  ! parameters
  !

  ! structures
  !
  type :: s_prmtop
    ! FLAG POINTER
    integer :: natom
    integer :: ntypes
    integer :: nbonh
    integer :: mbona
    integer :: ntheth
    integer :: mtheta
    integer :: nphih
    integer :: mphia
    integer :: nhparm
    integer :: nparm
    integer :: nnb
    integer :: nres
    integer :: nbona
    integer :: ntheta
    integer :: nphia
    integer :: numbnd
    integer :: numang
    integer :: nptra
    integer :: natyp
    integer :: nphb
    integer :: ifpert
    integer :: nbper
    integer :: ngper
    integer :: ndper
    integer :: mbper
    integer :: mgper
    integer :: mdper
    integer :: ifbox
    integer :: nmxrs
    integer :: ifcap
    integer :: numextra
    integer :: ncopy

    ! FLAG ATOM_TYPE_INDEX
    integer, allocatable          :: iac(:)
    ! FLAG NONBONDED_PARM_INDEX
    integer, allocatable          :: ico(:)
    ! FLAG LENNARD_JONES_ACOEF
    real(8), allocatable          :: cn1(:)
    ! FLAG LENNARD_JONES_BCOEF
    real(8), allocatable          :: cn2(:)
    ! FLAG CHARGE
    real(8), allocatable          :: charge(:)
    ! FLAG NUMBER_EXCLUDED_ATOMS
    integer, allocatable          :: num_excl_atm(:)
    ! FLAG EXCLUDED_ATOMS_LIST
    integer, allocatable          :: excl_atm_list(:)
    ! FLAG SCEE_SCALE_FACTOR
    real(8), allocatable          :: scee_scale_fact(:)
    ! FLAG SCNB_SCALE_FACTOR
    real(8), allocatable          :: scnb_scale_fact(:)
    ! FLAG DIHEDRALS_INC_HYDROGEN
    integer, allocatable          :: dihed_inc_hyd(:)
    ! FLAG DIHEDRALS_WITHOUT_HYDROGEN
    integer, allocatable          :: dihed_wo_hyd(:)

  end type s_prmtop

  ! subroutines
  !
  public   :: read_prmtop
  public   :: seek_prmtop_next_flag
  private  :: seek_prmtop_flag
  private  :: read_prmtop_integer 
  private  :: read_prmtop_real

  contains

!-----------------------------------------------------------------------
    subroutine read_prmtop(fprmtop, prmtop) 
!-----------------------------------------------------------------------
      implicit none

      integer, parameter :: npsize = 31 

      character(len=MaxChar), intent(in)  :: fprmtop
      type(s_prmtop),         intent(out) :: prmtop 

      integer :: i, j, k
      integer :: i10, ires, nsize 
      integer :: iunit, ierr
      integer :: pts(npsize)


      iunit = UnitPRMTOP
      open(iunit, file=trim(fprmtop))

      ! read POINTERS
      !
      nsize = npsize 
      call seek_prmtop_flag(iunit, 'POINTERS', ierr)
      call read_prmtop_integer(iunit, nsize, pts)

      prmtop%natom    = pts(1)
      prmtop%ntypes   = pts(2)
      prmtop%nbonh    = pts(3)
      prmtop%mbona    = pts(4)
      prmtop%ntheth   = pts(5)
      prmtop%mtheta   = pts(6)
      prmtop%nphih    = pts(7)
      prmtop%mphia    = pts(8)
      prmtop%nhparm   = pts(9)
      prmtop%nparm    = pts(10)
      prmtop%nnb      = pts(11)
      prmtop%nres     = pts(12)
      prmtop%nbona    = pts(13)
      prmtop%ntheta   = pts(14)
      prmtop%nphia    = pts(15)
      prmtop%numbnd   = pts(16)
      prmtop%numang   = pts(17)
      prmtop%nptra    = pts(18)
      prmtop%natyp    = pts(19)
      prmtop%nphb     = pts(20)
      prmtop%ifpert   = pts(21)
      prmtop%nbper    = pts(22)
      prmtop%ngper    = pts(23)
      prmtop%ndper    = pts(24)
      prmtop%mbper    = pts(25)
      prmtop%mgper    = pts(26)
      prmtop%mdper    = pts(27)
      prmtop%ifbox    = pts(28)
      prmtop%nmxrs    = pts(29)
      prmtop%ifcap    = pts(30)
      prmtop%numextra = pts(31)

      ! read ATOM_TYPE_INDEX 
      !
      allocate(prmtop%iac(prmtop%natom))
      !write(iw,'("ATOM_TYPE_INDEX")')
      call seek_prmtop_flag(iunit, 'ATOM_TYPE_INDEX', ierr)
      call read_prmtop_integer(iunit, prmtop%natom, prmtop%iac)

      ! read NONBONDED_PARM_INDEX
      !
      nsize = prmtop%ntypes * prmtop%ntypes
      allocate(prmtop%ico(nsize))
      !write(iw,'("NONBONDED_PARM_INDEX")')
      call seek_prmtop_flag(iunit, 'NONBONDED_PARM_INDEX', ierr)
      call read_prmtop_integer(iunit, nsize, prmtop%ico)

      ! read LENNARD_JONES_ACOEF
      !
      nsize = prmtop%ntypes * (prmtop%ntypes + 1) / 2
      allocate(prmtop%cn1(nsize))
      !write(iw,'("LENNARD_JONES_ACOEF")')
      call seek_prmtop_flag(iunit, 'LENNARD_JONES_ACOEF', ierr)
      call read_prmtop_real(iunit, nsize, prmtop%cn1)

      ! read LENNARD_JONES_BCOEF
      !
      nsize = prmtop%ntypes * (prmtop%ntypes + 1) / 2
      allocate(prmtop%cn2(nsize))
      !write(iw,'("LENNARD_JONES_BCOEF")')
      call seek_prmtop_flag(iunit, 'LENNARD_JONES_BCOEF', ierr)
      call read_prmtop_real(iunit, nsize, prmtop%cn2)

      ! read CHARGE 
      !
      nsize = prmtop%natom
      allocate(prmtop%charge(nsize))
      !write(iw,'("CHARGE")')
      call seek_prmtop_flag(iunit, 'CHARGE', ierr)
      call read_prmtop_real(iunit, nsize, prmtop%charge)

      ! read NUMBER_EXCLUDED_ATOMS
      !
      nsize = prmtop%natom
      allocate(prmtop%num_excl_atm(nsize))
      !write(iw,'("NUMBER_EXCLUDED_ATOMS")')
      call seek_prmtop_flag(iunit, 'NUMBER_EXCLUDED_ATOMS', ierr)
      call read_prmtop_integer(iunit, nsize, prmtop%num_excl_atm)

      ! read EXCLUDED_ATOMS_LIST
      !
      nsize = prmtop%nnb
      allocate(prmtop%excl_atm_list(nsize))
      !write(iw,'("EXCLUDED_ATOMS_LIST")')
      call seek_prmtop_flag(iunit, 'EXCLUDED_ATOMS_LIST', ierr)
      call read_prmtop_integer(iunit, nsize, prmtop%excl_atm_list)

      ! read SCEE_SCALE_FACTOR 
      !
      nsize = prmtop%nptra
      allocate(prmtop%scee_scale_fact(nsize))
      !write(iw,'("SCEE_SCALE_FACTOR")')
      call seek_prmtop_flag(iunit, 'SCEE_SCALE_FACTOR', ierr)
      call read_prmtop_real(iunit, nsize, prmtop%scee_scale_fact)

      ! read SCNB_SCALE_FACTOR 
      !
      nsize = prmtop%nptra
      allocate(prmtop%scnb_scale_fact(nsize))
      !write(iw,'("SCNB_SCALE_FACTOR")')
      call seek_prmtop_flag(iunit, 'SCNB_SCALE_FACTOR', ierr)
      call read_prmtop_real(iunit, nsize, prmtop%scnb_scale_fact)

      ! read DIHEDRALS_INC_HYDROGEN 
      !
      nsize = prmtop%nphih * 5
      allocate(prmtop%dihed_inc_hyd(nsize))
      !write(iw,'("DIHEDRALS_INC_HYDROGEN")')
      call seek_prmtop_flag(iunit, 'DIHEDRALS_INC_HYDROGEN', ierr)
      call read_prmtop_integer(iunit, nsize, prmtop%dihed_inc_hyd)

      ! read DIHEDRALS_WITHOUHT_HYDROGEN 
      !
      nsize = prmtop%nphia * 5
      allocate(prmtop%dihed_wo_hyd(nsize))
      !write(iw,'("DIHEDRALS_WITHOUT_HYDROGEN")')
      call seek_prmtop_flag(iunit, 'DIHEDRALS_WITHOUT_HYDROGEN', ierr)
      call read_prmtop_integer(iunit, nsize, prmtop%dihed_wo_hyd)

      close(iunit)


    end subroutine read_prmtop 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine seek_prmtop_next_flag(iunit, line, ierr) 
!-----------------------------------------------------------------------
      implicit none

      integer,                intent(in)  :: iunit
      character(*),           intent(out) :: line 
      integer,                intent(out) :: ierr

      character(len=MaxChar) :: c1, c2 
      logical                :: is_found


      is_found = .false.
      ierr = 0
      do while (.not. is_found)

        read(iunit,'(a)',end=100) line

        if (line(1:5) == '%FLAG') then
          return
        end if

      end do

      100 continue

      if (.not. is_found) then
        ierr = 1 
      end if


    end subroutine seek_prmtop_next_flag
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine seek_prmtop_flag(iunit, flag, ierr) 
!-----------------------------------------------------------------------
      implicit none

      integer,                intent(in)  :: iunit
      character(*),           intent(in)  :: flag
      integer,                intent(out) :: ierr

      character(len=MaxChar) :: line, c1, c2 
      logical                :: is_found


      rewind iunit

      is_found = .false.
      ierr = 0
      do while (.not. is_found)

        read(iunit,'(a)',end=100) line

        if (line(1:5) == '%FLAG') then
          read(line,*) c1, c2
          if (trim(c2) == trim(flag)) then
            is_found = .true.
            ! Skip FORMAT line
            read(iunit,*) 
          end if
        end if

      end do

      100 continue

      if (.not. is_found) then
        ierr = 1 
      end if


    end subroutine seek_prmtop_flag 
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine read_prmtop_integer(iunit, nsize, arr) 
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in)  :: iunit
      integer, intent(in)  :: nsize
      integer, intent(out) :: arr(nsize) 

      integer :: i, j, k
      integer :: i10, ires, nl


      nl   = 10

      i10  = nsize / nl 
      ires = nsize - i10 * nl 
      arr  = 0

      j    = 0
      do i = 1, i10
        read(iunit,*) (arr(j+k), k = 1, nl)
        j = j + nl 
      end do
      if (ires /= 0) then
        read(iunit,*)   (arr(j+k), k = 1, ires)
      end if

    end subroutine read_prmtop_integer
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine read_prmtop_real(iunit, nsize, arr) 
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in)  :: iunit
      integer, intent(in)  :: nsize
      real(8), intent(out) :: arr(nsize) 

      integer, parameter   :: ncolmax = 5 

      integer                :: i, j, k
      integer                :: il, ires, nl
      integer                :: ierr, icol, ncol
      integer                :: col(ncolmax+1) 
      character(len=MaxChar) :: line

      ! check # of columns
      !
      read(iunit,'(a)') line
      backspace(iunit)

      ncol = ncolmax + 1
100   ncol = ncol - 1
      if (ncol == 0) then
        write(iw,'("Read_Prmtop_Real> Error.")')
        write(iw,'("No data on FLAG.")')
        stop
      end if

      read(line,*,iostat=ierr) (col(icol), icol = 1, ncol)
      if (ierr /= 0) then
        go to 100
      end if

      nl  = ncol
      !nl = 5

      il   = nsize / nl 
      ires = nsize - il * nl 
      arr  = 0

      j    = 0
      do i = 1, il
        read(iunit,*) (arr(j+k), k = 1, nl)
        j = j + nl 
      end do
      if (ires /= 0) then
        read(iunit,*)   (arr(j+k), k = 1, ires)
      end if

    end subroutine read_prmtop_real
!-----------------------------------------------------------------------

end module mod_prmtop
!=======================================================================
