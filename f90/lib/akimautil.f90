!======================================================================
! Module: mod_akima
!
!   Module for Akima interpolation scheme
!   (Ref. H. Akima, J. ACM, 17, 589-602 (1970)) 
!
!   (c) Copyright 2023 Osaka Univ. All rights reserved.
!
!   TODO: Modified Akima should be implemented 
!
!======================================================================
module mod_akima


  ! Subroutines
  public :: akima1d
  public :: akima2d
  public :: akima3d
  public :: akima
  public :: akima_prepare
  public :: akima_interpolate
  public :: get_lag2nd_intpl
  public :: extend_data2d
  public :: extend_data3d

  contains
!----------------------------------------------------------------------
!   Subroutine Akima1d
!
!   Perform Akima interpolation in 1D-space
!
!   - integer ng3t(1:3)     [IN]  :
!   - integer ng3(1:3)      [IN]  : 
!   - real(8) box(1:3)      [IN]  :
!   - real(8) lcut          [IN]  :
!   - real(8) ucut          [IN]  :
!   - real(8) fold(ng3t(1)) [IN]  : 
!   - real(8) fnew(ng3(1))  [OUT] : 
!     
    subroutine akima1d(ng3t,ng3,box,lcut,ucut,fold,fnew)
!----------------------------------------------------------------------
!$    use omp_lib
!
      implicit none
!
      integer, intent(in)              :: ng3t(3),ng3(3)
      real(8), intent(in)              :: box(3), ucut, lcut, fold(ng3t(1))
      real(8), intent(out)             :: fnew(ng3(1))  
!
      real(8)                          :: del3t(3),del3(3)
      real(8), allocatable             :: ftemp(:)
!
      integer                          :: ig,igx
      integer                          :: i1
      real(8)                          :: x
      real(8)                          :: fout


      allocate(ftemp(ng3(1)))
!
      del3t(:)= box(:)/dfloat(ng3t(:))
      del3(:) = box(:)/dfloat(ng3(:))
!
      do igx=1,ng3t(1)
        ftemp(igx)= fold(igx) 
      end do
!
      do i1=0,ng3(1)-1    
        x= del3(1)* dfloat(i1)
        call akima(ng3t(1),0.0d0,del3t(1),x,ftemp(1),fout)
!
        if( fout < lcut ) fout= lcut
        if( fout > ucut ) fout= ucut 
!
        fnew(i1+1)= fout
      end do           ! i1
!
    end subroutine akima1d
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
    subroutine akima2d(ng3t,ng3,box,lcut,ucut,fold,fnew)
!----------------------------------------------------------------------
!$    use omp_lib
!
      implicit none
!
      integer, intent(in)              :: ng3t(3),ng3(3)
      real(8), intent(in)              :: box(3), ucut, lcut, fold(ng3t(1),ng3t(2))
      real(8), intent(out)             :: fnew(ng3(1),ng3(2))  
!
      real(8)                          :: del3t(3),del3(3)
      real(8), allocatable             :: ftemp(:,:), f2d(:,:),f1d(:)
!
      integer                          :: ig,igx,igy,igz
      integer                          :: i1,i2,i3
      real(8)                          :: x,y,z
      real(8)                          :: fout
!
!$    integer :: nchnk
!
      allocate(ftemp(ng3(1), ng3(2)), f1d(ng3t(2)))
!
      del3t(:)= box(:)/dfloat(ng3t(:))
      del3(:) = box(:)/dfloat(ng3(:))
!
      do igy=1,ng3t(2)
        do igx=1,ng3t(1)
          ftemp(igx, igy)= fold(igx, igy) 
        end do
      end do
!
      !call extend_data2d(ng3t(1),ng3t(2),ftemp) 
!
!!$    nchnk= ng3(1)/omp_get_max_threads()
!!$    if( nchnk == 0 ) nchnk= 1
!!$omp parallel
!!$omp do schedule(static,nchnk) private(i1,i2,x,y,ig,igy,f1d,fout)
      do i1=0,ng3(1)-1    
        x= del3(1)* dfloat(i1)
        do igy=1,ng3t(2)
          call akima(ng3t(1),0.0d0,del3t(1),x,ftemp(1,igy),f1d(igy))
        end do
!
        do i2=0,ng3(2)-1
          y= del3(2)* dfloat(i2)
!          do igz=1,ng3t(3)+1
            call akima(ng3t(2),0.0d0,del3t(2),y,f1d,fout)
!          end do
!
!
          if( fout < lcut ) fout= lcut
          if( fout > ucut ) fout= ucut 
!
          fnew(i1+1,i2+1)= fout
!
        end do         ! i2
      end do           ! i1
!
!!$omp end do
!!$omp end parallel
!
    end subroutine akima2d
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
    subroutine akima3d(ng3t,ng3,box,lcut,ucut,fold,fnew)
!----------------------------------------------------------------------
!$    use omp_lib
!
      implicit none
!
      integer, intent(in)              :: ng3t(3),ng3(3)
      real(8), intent(in)              :: box(3), ucut, lcut, fold(ng3t(1),ng3t(2),ng3t(3))
      real(8), intent(out)             :: fnew(ng3(1),ng3(2),ng3(3))  
!
      real(8)                          :: del3t(3),del3(3)
      real(8), allocatable             :: ftemp(:,:,:),f2d(:,:),f1d(:)
!
      integer                          :: ig,igx,igy,igz
      integer                          :: i1,i2,i3
      real(8)                          :: x,y,z
      real(8)                          :: fout
!
!$    integer :: nchnk
!
      allocate(ftemp(ng3(1),ng3(2),ng3(3)),f2d(ng3t(2),ng3t(3)),f1d(ng3t(3)))
!
      del3t(:)= box(:)/dfloat(ng3t(:))
      del3(:)= box(:)/dfloat(ng3(:))
!
      do igz=1,ng3t(3)
        do igy=1,ng3t(2)
          do igx=1,ng3t(1)
            ftemp(igx,igy,igz)= fold(igx,igy,igz) 
          end do
        end do
      end do
!
      !call extend_data3d(ng3t(1),ng3t(2),ng3t(3),ftemp) 
!
!$    nchnk= ng3(1)/omp_get_max_threads()
!$    if( nchnk == 0 ) nchnk= 1
!$omp parallel
!$omp do schedule(static,nchnk) private(i1,i2,i3,x,y,z,ig,igy,igz,f2d,f1d,fout)
      do i1=0,ng3(1)-1    
        x= del3(1)* dfloat(i1)
        do igz=1,ng3t(3)
          do igy=1,ng3t(2)
            call akima(ng3t(1),0.0d0,del3t(1),x,ftemp(1,igy,igz),f2d(igy,igz))
          end do
        end do
!
        do i2=0,ng3(2)-1
          y= del3(2)* dfloat(i2)
          do igz=1,ng3t(3)
            call akima(ng3t(2),0.0d0,del3t(2),y,f2d(1,igz),f1d(igz))
          end do
!
          do i3=0,ng3(3)-1
            z= del3(3) * dfloat(i3)
            call akima(ng3t(3),0.0d0,del3t(3),z,f1d,fout)
!
            if( fout < lcut ) fout= lcut
            if( fout > ucut ) fout= ucut 
!
            ig= 1 + i1 + i2*ng3(1) + i3*ng3(1)*ng3(2)
            fnew(i1+1,i2+1,i3+1)= fout
!
          end do       ! i3
        end do         ! i2
      end do           ! i1
!
!$omp end do
!$omp end parallel
!
    end subroutine akima3d
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
    subroutine akima(ngrid,shift,del,val,func,fout)
!----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
!
      real*8 func(ngrid),fout
      real*8 rtab(ngrid)
      real*8 rtab_ex(ngrid+4),func_ex(ngrid+4)
!
      do i=0,ngrid-1
        rtab(i+1) = shift + dfloat(i)*del
      end do
!      
      call akima_prepare(ngrid,rtab,func,rtab_ex,func_ex)
      call akima_scheme(ngrid,rtab_ex,func_ex,val,fout)
!
    end subroutine akima
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
    subroutine akima_prepare(ngrid,rtab,func,rtab_ex,func_ex)
!----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
!
      real*8 rtab(ngrid),func(ngrid)
      real*8 rtab_ex(ngrid+4),func_ex(ngrid+4)
!
      do i=3,ngrid+2
        rtab_ex(i) = rtab(i-2)
        func_ex(i) = func(i-2)
      end do
!
      del = rtab(2) - rtab(1)
!
      rtab_ex(1) = rtab(1)-del*2.0d0
      rtab_ex(2) = rtab(1)-del
!
      rtab_ex(ngrid+3) = rtab(ngrid) + del
      rtab_ex(ngrid+4) = rtab(ngrid) + del*2.0d0
!
! estimate func at i=1,2 with 2nd order lagrange interpolation formula
!
      r1=rtab(1) ; f1=func(1)
      r2=rtab(2) ; f2=func(2)
      r3=rtab(3) ; f3=func(3)
!
      r=rtab_ex(1)
      func_ex(1)=get_lag2nd_intpl(r1,r2,r3,f1,f2,f3,r) 
      r=rtab_ex(2)
      func_ex(2)=get_lag2nd_intpl(r1,r2,r3,f1,f2,f3,r)
!
! estimate func at i=ngrid+1,ngrid+2 with 2nd order lagrange interpolation formula
!
      r1=rtab(ngrid-2) ; f1=func(ngrid-2)
      r2=rtab(ngrid-1) ; f2=func(ngrid-1)
      r3=rtab(ngrid)   ; f3=func(ngrid)
!
      r=rtab_ex(ngrid+3)
      func_ex(ngrid+3)=get_lag2nd_intpl(r1,r2,r3,f1,f2,f3,r)
      r=rtab_ex(ngrid+4)
      func_ex(ngrid+4)=get_lag2nd_intpl(r1,r2,r3,f1,f2,f3,r)
!
!  do i=1,ngrid+4
!    write(*,'(1x,e15.7,5x,e15.7)') rtab_ex(i),func_ex(i)
!  end do
!
    end subroutine akima_prepare
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
    subroutine akima_scheme(ngrid,rtab_ex,func_ex,r,fout)
!----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
!
      real*8 rtab_ex(ngrid+4),func_ex(ngrid+4) 
      real*8 rval(6),fval(6)
!
      ngrid_ex = ngrid+4
!
      do i=1,ngrid_ex
        r0 = rtab_ex(i)
        r1 = rtab_ex(i+1)
        if(r.ge.r0.and.r.lt.r1) then
          ipos=i
          exit
        end if 
      end do
!
      if(ipos.eq.ngrid+2) then
        fout=func_ex(ngrid+2)
      else
!
        do i=1,6
          rval(i) = rtab_ex(ipos+i-3)
          fval(i) = func_ex(ipos+i-3)
        end do 
!
        call akima_interpolate(rval,fval,r,fout)
!
      end if
!
    end subroutine
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
    subroutine akima_interpolate(rval,fval,r,fout)
!----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
!
      real*8 rval(6),fval(6)
!
      f_m2=fval(1)
      f_m1=fval(2)
      f_0 =fval(3)
      f_p1=fval(4)
      f_p2=fval(5)
      f_p3=fval(6)
!
      r_m2=rval(1)
      r_m1=rval(2)
      r_0 =rval(3)
      r_p1=rval(4)
      r_p2=rval(5)
      r_p3=rval(6)
!
      grad_m2=(f_m1-f_m2)/(r_m1-r_m2)
      grad_m1=(f_0 -f_m1)/(r_0 -r_m1)
      grad_0 =(f_p1-f_0 )/(r_p1-r_0 )
      grad_p1=(f_p2-f_p1)/(r_p2-r_p1) 
      grad_p2=(f_p3-f_p2)/(r_p3-r_p2)
!
      denom_rt0=dabs(grad_p1-grad_0)+dabs(grad_m1-grad_m2)
      r_num_rt0=dabs(grad_p1-grad_0)*grad_m1+dabs(grad_m1-grad_m2)*grad_0
!
      denom_rt1=dabs(grad_p2-grad_p1)+dabs(grad_0-grad_m1)
      r_num_rt1=dabs(grad_p2-grad_p1)*grad_0+dabs(grad_0-grad_m1)*grad_p1
!
      if(denom_rt0.le.1.0d-10) then
        rt0 = 0.5d0*(grad_m1+grad_0) 
      else
        rt0 = r_num_rt0/denom_rt0
      end if
!
      if(denom_rt1.le.1.0d-10) then
        rt1 = 0.5d0*(grad_0+grad_p1) 
      else
        rt1 = r_num_rt1/denom_rt1
      end if
!
      a0 = f_0
      a1 = rt0
      a2 = (3.0d0*grad_0-2.0d0*rt0-rt1)/(r_p1-r_0)
      a3 = (rt0+rt1-2.0d0*grad_0)/((r_p1-r_0)*(r_p1-r_0))
!
      x  = r-r_0
      x2 = x*x
      x3 = x2*x
!
      fout = a0 + a1*x + a2*x2 + a3*x3
!
    end subroutine akima_interpolate
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
    function get_lag2nd_intpl(r1,r2,r3,f1,f2,f3,r)
!----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
!
      real*8 get_lag2nd_intpl
!
      get_lag2nd_intpl = ( (r-r2)*(r-r3)/((r1-r2)*(r1-r3)) ) * f1 &
                   &+( (r-r1)*(r-r3)/((r2-r1)*(r2-r3)) ) * f2 &
                   &+( (r-r1)*(r-r2)/((r3-r1)*(r3-r2)) ) * f3
!
    end function get_lag2nd_intpl
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
    subroutine extend_data2d(ngx,ngy,wrk)
!----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
!
      real*8 wrk(ngx+1,ngy+1)
!
      do igz=1,ngz
        do igy=1,ngy
          wrk(ngx+1,igy) = wrk(ngx,igy)
        end do
      end do
!
      do igz=1,ngz
        do igx=1,ngx+1
          wrk(igx,ngy+1) = wrk(igx,ngy)
        end do
      end do
!
    end subroutine extend_data2d
!----------------------------------------------------------------------

!----------------------------------------------------------------------
    subroutine extend_data3d(ngx,ngy,ngz,wrk)
!----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
!
      real*8 wrk(ngx+1,ngy+1,ngz+1)
!
      do igz=1,ngz
        do igy=1,ngy
          wrk(ngx+1,igy,igz) = wrk(ngx,igy,igz)
        end do
      end do
!
      do igz=1,ngz
        do igx=1,ngx+1
          wrk(igx,ngy+1,igz) = wrk(igx,ngy,igz)
        end do
      end do
!
      do igy=1,ngy+1
        do igx=1,ngx+1
          wrk(igx,igy,ngz+1) = wrk(igx,igy,ngz)
        end do
      end do
!
    end subroutine extend_data3d
!----------------------------------------------------------------------

end module mod_akima
!======================================================================
