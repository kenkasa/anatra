!=======================================================================
module mod_pme_str
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_grid3d

  ! constants
  !
  integer, parameter :: ngmax = 1000 

  ! structures
  !
  type :: s_pmevar
    integer :: ng3(3)
    real(8) :: box(3)
    real(8) :: del(3)
    real(8) :: dv(3)

    real(8),    allocatable :: dfunc (:, :, :)
    real(8),    allocatable :: qfunck(:, :, :)
    complex(kind(0d0)), allocatable :: qfuncm(:, :, :)
    real(8),    allocatable :: gfunck(:, :, :)
    complex(kind(0d0)), allocatable :: gfuncm(:, :, :)
    real(8),    allocatable :: splfc (:, :)

    real(8),    allocatable :: qfunck_dual(:, :, :)
    complex(kind(0d0)), allocatable :: qfuncm_dual(:, :, :)
    real(8),    allocatable :: gfunck_dual(:, :, :)
    complex(kind(0d0)), allocatable :: gfuncm_dual(:, :, :)

    real(8),    allocatable :: qderiv(:, :, :, :)
    real(8),    allocatable :: qderiv_dual(:, :, :, :)
    !real,    allocatable :: dfunc (:, :, :)
    !real,    allocatable :: qfunck(:, :, :)
    !complex, allocatable :: qfuncm(:, :, :)
    !real,    allocatable :: gfunck(:, :, :)
    !complex, allocatable :: gfuncm(:, :, :)
    !real,    allocatable :: splfc (:, :)
    !
    !real,    allocatable :: qfunck_dual(:, :, :)
    !complex, allocatable :: qfuncm_dual(:, :, :)
    !real,    allocatable :: gfunck_dual(:, :, :)
    !complex, allocatable :: gfuncm_dual(:, :, :)
    !
    !real,    allocatable :: qderiv(:, :, :, :)
    !real,    allocatable :: qderiv_dual(:, :, :, :)
  end type s_pmevar

  ! subroutines
  !
  public :: dealloc_pmevar

  contains
!-----------------------------------------------------------------------
    subroutine dealloc_pmevar(pmevar)
!-----------------------------------------------------------------------
      implicit none

      type(s_pmevar), intent(inout) :: pmevar


      deallocate(pmevar%dfunc)
      deallocate(pmevar%qfunck)
      deallocate(pmevar%qfuncm)
      deallocate(pmevar%gfunck)
      deallocate(pmevar%gfuncm)
      deallocate(pmevar%splfc)
      deallocate(pmevar%qderiv)

      deallocate(pmevar%qfunck_dual)
      deallocate(pmevar%qfuncm_dual)
      deallocate(pmevar%gfunck_dual)
      deallocate(pmevar%gfuncm_dual)
      deallocate(pmevar%qderiv_dual)

    end subroutine dealloc_pmevar
!-----------------------------------------------------------------------
!
end module mod_pme_str
!=======================================================================
