!-------------------------------------------------------------------------------
!
!  Module   mod_fftmkl
!
!    Module for using MKL-FFT library   
!
!    (c) Copyright 2024 Osaka Univ. All rights reserved.
!
!-------------------------------------------------------------------------------

include 'mkl_dfti.f90'

module mod_fftmkl
  use MKL_DFTI 

  ! constants
  !
  integer, parameter :: FFT_SignForward  = 1
  integer, parameter :: FFT_SignBackward = 2 

  ! structures
  !
  type s_fftinfo
    type(Dfti_Descriptor), pointer :: desc_forward
    type(Dfti_Descriptor), pointer :: desc_backward
    integer    :: ng3(3)
  end type s_fftinfo

  ! variables
  !
  integer :: fftsize(3)

  ! subroutines
  !
  public :: fftmkl_init
  public :: fftmkl_r2c
  public :: fftmkl_c2r

  contains
!-----------------------------------------------------------------------
      subroutine fftmkl_init(fftinfo, funck, funcm)
!-----------------------------------------------------------------------
      implicit none
!
      type(s_fftinfo),     intent(inout) :: fftinfo
      !real,                intent(inout) :: funck(0:fftinfo%ng3(1)-1,&
      !                                            0:fftinfo%ng3(2)-1,&
      !                                            0:fftinfo%ng3(3)-1)
      !complex,             intent(inout) :: funcm(0:fftinfo%ng3(1)/2,&
      !                                            0:fftinfo%ng3(2)-1,&
      !                                            0:fftinfo%ng3(3)-1)
      real(8),             intent(inout) :: funck(fftinfo%ng3(1),       &
                                                  fftinfo%ng3(2),       &
                                                  fftinfo%ng3(3))
      !real,                intent(inout) :: funck(fftinfo%ng3(1),       &
      !                                            fftinfo%ng3(2),       &
      !                                            fftinfo%ng3(3))
      complex(kind(0d0)),  intent(inout) :: funcm(fftinfo%ng3(1)/2 + 1, &
                                                  fftinfo%ng3(2),       &
                                                  fftinfo%ng3(3))
      !complex,             intent(inout) :: funcm(fftinfo%ng3(1)/2 + 1, &
      !                                            fftinfo%ng3(2),       &
      !                                            fftinfo%ng3(3))

      integer :: igx, igy, igz, igr, igk
      integer :: ngx, ngy, ngz, ngr, ngk
      integer :: statf, statb
      integer :: strides(4)

      
      ! setup module variables
      !
      fftsize(1:3) = fftinfo%ng3(1:3)

      
      ngx = fftinfo%ng3(1)
      ngy = fftinfo%ng3(2)
      ngz = fftinfo%ng3(3)

      strides(1) = 0
      strides(2) = 1
      strides(3) = ngx / 2 + 1
      strides(4) = strides(3) * ngy

      ! Forward
      !
      statf = DftiCreateDescriptor(fftinfo%desc_forward,  &
                                   DFTI_DOUBLE,           &
                                   DFTI_REAL,             &
                                   3,                     &
                                   fftinfo%ng3)
      !statf = DftiCreateDescriptor(fftinfo%desc_forward,  &
      !                             DFTI_SINGLE,           &
      !                             DFTI_REAL,             &
      !                             3,                     &
      !                             fftinfo%ng3)

      statf = DftiSetValue(fftinfo%desc_forward,          &
                           DFTI_PLACEMENT,                &
                           DFTI_NOT_INPLACE)
   
      statf = DftiSetValue(fftinfo%desc_forward,          &
                           DFTI_CONJUGATE_EVEN_STORAGE,   &
                           DFTI_COMPLEX_COMPLEX)

      statf = DftiSetValue(fftinfo%desc_forward,          &
                           DFTI_OUTPUT_STRIDES,           &
                           strides)

      statf = DftiSetValue(fftinfo%desc_forward,          &
                          DFTI_PACKED_FORMAT,             &
                          DFTI_CCE_FORMAT)

      statf = DftiCommitDescriptor(fftinfo%desc_forward)

      ! Backward
      !
      statb = DftiCreateDescriptor(fftinfo%desc_backward, &
                                   DFTI_DOUBLE,           &
                                   DFTI_REAL,             &
                                   3,                     &
                                   fftinfo%ng3)
      !statb = DftiCreateDescriptor(fftinfo%desc_backward, &
      !                             DFTI_SINGLE,           &
      !                             DFTI_REAL,             &
      !                             3,                     &
      !                             fftinfo%ng3)

      statb = DftiSetValue(fftinfo%desc_backward,         &
                           DFTI_PLACEMENT,                &
                           DFTI_NOT_INPLACE)
   
      statb  = DftiSetValue(fftinfo%desc_backward,        &
                           DFTI_CONJUGATE_EVEN_STORAGE,   &
                           DFTI_COMPLEX_COMPLEX)

      statb = DftiSetValue(fftinfo%desc_backward,         &
                           DFTI_INPUT_STRIDES,            &
                           strides)

      statb = DftiSetValue(fftinfo%desc_backward,         &
                          DFTI_PACKED_FORMAT,             &
                          DFTI_CCE_FORMAT)

      statb = DftiCommitDescriptor(fftinfo%desc_backward)

      end subroutine fftmkl_init
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine fftmkl_r2c(fftinfo, in, out)
!-----------------------------------------------------------------------
      use MKL_DFTI 
      implicit none 

      type(s_fftinfo),     intent(in)  :: fftinfo
      real(8),             intent(in)  :: in(*)
      complex(kind(0d0)),  intent(out) :: out(*)

      integer :: stat
      integer :: ng

      stat = DftiComputeForward(fftinfo%desc_forward, in, out)

      !ng = (fftsize(1) / 2 + 1) * fftsize(2) * fftsize(3) 
      !out(1:ng) = out(1:ng) &
      !            / dble(fftsize(1) * fftsize(2) * fftsize(3))

      end subroutine fftmkl_r2c
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine fftmkl_c2r(fftinfo, in, out)
!-----------------------------------------------------------------------
      use MKL_DFTI 
      implicit none 

      type(s_fftinfo),     intent(in)  :: fftinfo
      complex(kind(0d0)),  intent(in)  :: in(*)
      real(8),             intent(out) :: out(*)

      integer :: stat

      stat = DftiComputeBackward(fftinfo%desc_backward, in, out) 
!
      end subroutine fftmkl_c2r
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine fftmkl_cleanup(fftinfo)
!-----------------------------------------------------------------------
      implicit none 

      type(s_fftinfo),     intent(in)  :: fftinfo

      integer :: stat

      stat = DftiFreeDescriptor(fftinfo%desc_forward)
      stat = DftiFreeDescriptor(fftinfo%desc_backward)
!
      end subroutine fftmkl_cleanup
!----------------------------------------------------------------------

end module mod_fftmkl
