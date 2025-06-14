!=======================================================================
module mod_fitting
!=======================================================================
  use mod_util
  use mod_const
  use mod_traj

  type :: s_fit
    integer              :: natm
    real(8), allocatable :: mass(:)
    real(8), allocatable :: refcoord(:, :)
    real(8), allocatable :: movcoord(:, :)
    real(8)              :: refcom(3)
    real(8)              :: movcom(3)
    real(8)              :: rot_matrix(3, 3)
    integer, allocatable :: ind(:)
  end type s_fit

  ! subroutines
  !
  public :: setup_fit
  public :: get_com_refmov
  public :: get_trrot
  public :: operate_trrot 

  contains

!-----------------------------------------------------------------------
    subroutine setup_fit(traj, fit, refcoord) 
!-----------------------------------------------------------------------
      implicit none

      type(s_traj),               intent(in)    :: traj
      type(s_fit),                intent(inout) :: fit
      real(8),          optional, intent(in)    :: refcoord(:, :)

      integer :: i
      integer :: natm


      natm = traj%natm

      allocate(fit%mass(natm), fit%refcoord(3, natm), fit%movcoord(3, natm))
      allocate(fit%ind(natm))
      do i = 1, natm
        fit%mass(i) = traj%mass(i)
        !write(6,'("mass ",f20.10)') fit%mass(i)
      end do

      fit%refcoord = 0.0d0
      fit%movcoord = 0.0d0

      if (present(refcoord)) then
        fit%refcoord(:, :) = refcoord(:, :)
      end if

      fit%natm   = traj%natm
      fit%ind(:) = traj%ind(:)
!
    end subroutine setup_fit
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_com_refmov(fit) 
!-----------------------------------------------------------------------
      implicit none

      integer,           parameter     :: lwork = max(1,4*3-1)

      type(s_fit),       intent(inout) :: fit

      integer :: i, natom, ierr
      real(8) :: total_mass

      natom      = fit%natm
      fit%refcom = 0.0d0
      fit%movcom = 0.0d0
      total_mass = 0.0d0 

      do i = 1, natom
        fit%refcom(1:3) = fit%refcom(1:3) + fit%mass(i) * fit%refcoord(1:3, i) 
        fit%movcom(1:3) = fit%movcom(1:3) + fit%mass(i) * fit%movcoord(1:3, i) 
        total_mass      = total_mass + fit%mass(i)
      end do

      total_mass = 1.0d0 / total_mass
      fit%refcom(1:3) = fit%refcom(1:3) * total_mass
      fit%movcom(1:3) = fit%movcom(1:3) * total_mass

!
    end subroutine get_com_refmov
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_trrot(fit) 
!-----------------------------------------------------------------------
      implicit none

      integer,           parameter     :: lwork = max(1,4*3-1)

      type(s_fit),       intent(inout) :: fit

      real(8) :: vec(4), eigen_work(4, 4) 

      integer :: i, natom, ierr
      real(8) :: total_mass
      real(8) :: sym_matrix(4, 4), eval(1:4), work(1:lwork), evec(1:4)
      real(8) :: dref(3), dmov(3), dsub(3), dadd(3)


      natom      = fit%natm
      fit%refcom = 0.0d0
      fit%movcom = 0.0d0
      total_mass = 0.0d0 

      do i = 1, natom
        fit%refcom(1:3) = fit%refcom(1:3) + fit%mass(i) * fit%refcoord(1:3, i) 
        fit%movcom(1:3) = fit%movcom(1:3) + fit%mass(i) * fit%movcoord(1:3, i) 
        total_mass      = total_mass + fit%mass(i)

        !write(6,'("refcoord ",3f20.10)') fit%refcoord(1:3, i)
        !write(6,'("movcoord ",3f20.10)') fit%movcoord(1:3, i)
      end do

      total_mass = 1.0d0 / total_mass
      fit%refcom(1:3) = fit%refcom(1:3) * total_mass
      fit%movcom(1:3) = fit%movcom(1:3) * total_mass

      sym_matrix(1:4, 1:4) = 0.0d0

      do i = 1, natom
        dref(1:3) = fit%refcoord(1:3, i) - fit%refcom(1:3)
        dmov(1:3) = fit%movcoord(1:3, i) - fit%movcom(1:3)
        dsub(1:3) = fit%mass(i) * (dref(1:3) - dmov(1:3))
        dadd(1:3) = fit%mass(i) * (dref(1:3) + dmov(1:3))

        sym_matrix(1,1) = sym_matrix(1,1) + dsub(1)*dsub(1) &
                                          + dsub(2)*dsub(2) &
                                          + dsub(3)*dsub(3)  
        sym_matrix(1,2) = sym_matrix(1,2) + dadd(2)*dsub(3) - dsub(2)*dadd(3)
        sym_matrix(1,3) = sym_matrix(1,3) + dsub(1)*dadd(3) - dadd(1)*dsub(3)
        sym_matrix(1,4) = sym_matrix(1,4) + dadd(1)*dsub(2) - dsub(1)*dadd(2)
        sym_matrix(2,2) = sym_matrix(2,2) + dsub(1)*dsub(1) &
                                          + dadd(2)*dadd(2) &
                                          + dadd(3)*dadd(3)  
        sym_matrix(2,3) = sym_matrix(2,3) + dsub(1)*dsub(2) - dadd(1)*dadd(2)
        sym_matrix(2,4) = sym_matrix(2,4) + dsub(1)*dsub(3) - dadd(1)*dadd(3)
        sym_matrix(3,3) = sym_matrix(3,3) + dadd(1)*dadd(1) &
                                          + dsub(2)*dsub(2) &
                                          + dadd(3)*dadd(3)  
        sym_matrix(3,4) = sym_matrix(3,4) + dsub(2)*dsub(3) - dadd(2)*dadd(3)
        sym_matrix(4,4) = sym_matrix(4,4) + dadd(1)*dadd(1) &
                                          + dadd(2)*dadd(2) &
                                          + dsub(3)*dsub(3)  
      end do

      sym_matrix(2,1) = sym_matrix(1,2)
      sym_matrix(3,1) = sym_matrix(1,3)
      sym_matrix(3,2) = sym_matrix(2,3)
      sym_matrix(4,1) = sym_matrix(1,4)
      sym_matrix(4,2) = sym_matrix(2,4)
      sym_matrix(4,3) = sym_matrix(3,4)

      call dsyev('V', 'U', 4, sym_matrix, 4, eval, work, lwork, ierr)
      !call eigen(4, 4, sym_matrix, vec, eval, eigen_work) 

      evec(1:4) = sym_matrix(1:4, 1)

      fit%rot_matrix = 0.0d0

      fit%rot_matrix(1,1) =           evec(1)*evec(1) + evec(2)*evec(2)  &
                                     -evec(3)*evec(3) - evec(4)*evec(4)
      fit%rot_matrix(1,2) = 2.0d0  * (evec(2)*evec(3) + evec(1)*evec(4))
      fit%rot_matrix(1,3) = 2.0d0  * (evec(2)*evec(4) - evec(1)*evec(3))
      fit%rot_matrix(2,1) = 2.0d0  * (evec(2)*evec(3) - evec(1)*evec(4))
      fit%rot_matrix(2,2) =           evec(1)*evec(1) - evec(2)*evec(2)  &
                                     +evec(3)*evec(3) - evec(4)*evec(4)
      fit%rot_matrix(2,3) = 2.0d0  * (evec(3)*evec(4) + evec(1)*evec(2))
      fit%rot_matrix(3,1) = 2.0d0  * (evec(2)*evec(4) + evec(1)*evec(3))
      fit%rot_matrix(3,2) = 2.0d0  * (evec(3)*evec(4) - evec(1)*evec(2))
      fit%rot_matrix(3,3) =           evec(1)*evec(1) - evec(2)*evec(2)  &
                                     -evec(3)*evec(3) + evec(4)*evec(4)

      !fit%rot_matrix(:,:) = 0.0d0
      !fit%rot_matrix(1,1) = 1.0d0
      !fit%rot_matrix(2,2) = 1.0d0
      !fit%rot_matrix(3,3) = 1.0d0

!
    end subroutine get_trrot 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine operate_trrot(fit, natm, coord) 
!-----------------------------------------------------------------------
      implicit none

      type(s_fit),       intent(inout)   :: fit
      integer,           intent(in)      :: natm
      real(8),           intent(inout)   :: coord(3, natm)

      integer :: iatm, natom
      real(8) :: dmov(1:3), rot_mov(1:3)


      natom      = natm

      do iatm = 1, natom
      
        dmov(1:3)       = coord(1:3,iatm) - fit%movcom(1:3)
        rot_mov(1:3)    = fit%rot_matrix(1:3,1) * dmov(1) + &
                          fit%rot_matrix(1:3,2) * dmov(2) + &
                          fit%rot_matrix(1:3,3) * dmov(3)
        coord(1:3,iatm) = rot_mov(1:3) + fit%refcom(1:3)
      
      end do
!
    end subroutine operate_trrot 
!-----------------------------------------------------------------------

end module mod_fitting
!=======================================================================
