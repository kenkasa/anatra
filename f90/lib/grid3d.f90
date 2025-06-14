!=======================================================================
module mod_grid3d
!=======================================================================
  use mod_const

  implicit none

  ! structures 
  !
  type s_func3d
    integer                  :: ng3(3)
    real(8)                  :: del(3)
    real(8)                  :: origin(3)
    real(8)                  :: box(3)
    real(8), allocatable     :: data(:, :, :) 
  end type s_func3d

  type s_mesh3d
    integer                  :: ng3(3)
    real(8)                  :: del(3)
    real(8)                  :: origin(3)
    real(8)                  :: box(3)
    real(8), allocatable     :: data(:, :, :, :) 
  end type s_mesh3d

  type s_spregion
    integer                  :: nspr
    real(8)                  :: del(3)
    real(8)                  :: origin(3)
    real(8)                  :: box(3)
    real(8), allocatable     :: data(:, :)
  end type s_spregion

  type s_mpl2d
    character(len=MaxChar)   :: fpdf          = "pmf.pdf"
    character(len=MaxChar)   :: labels(3)     = (/"$x$", "$y$", "$f$"/)
    real(8)                  :: ranges(2, 3)  = reshape((/1,2,3,4,5,6/), (/2,3/)) 
    real(8)                  :: tics(3)       = 0.0d0 
    real(8)                  :: scales(3)     = 0.0d0 
  end type s_mpl2d


  ! subroutines
  !
  public :: read_ctrl_matplotlib2d
  public :: setup_func3d
  public :: setup_mesh3d
  public :: generate_script_matplotlib2d

  contains
!-----------------------------------------------------------------------
    subroutine read_ctrl_matplotlib2d(iunit, mpl2d)
!-----------------------------------------------------------------------
      implicit none
!
      integer,        intent(in)  :: iunit
      type(s_mpl2d),  intent(out) :: mpl2d 

      character(len=MaxChar) :: fpdf          = "pmf.pdf"
      character(len=MaxChar) :: labels(3)     = (/"$x$", "$y$", "$f$"/)
      real(8)                :: ranges(2, 3)  = reshape((/1,2,3,4,5,6/), (/2,3/)) 
      real(8)                :: tics(3)       = (/2.0d0, 0.5d0, 0.1d0/)
      real(8)                :: scales(3)     = (/1.0d0, 1.0d0, 0.5d0/)

      integer :: i, j
      integer :: iopt, ierr

      namelist /matplotlib_param/ fpdf, labels, ranges, tics, scales 

      rewind iunit
      read(iunit, matplotlib_param)

      write(iw,*)
      write(iw,'(">> Matplotlib section parameters")')
      write(iw,'("fpdf             = ", a)')        trim(fpdf) 
      write(iw,*)
      write(iw,'("labels")')
      write(iw,'("  x : ",a)') trim(labels(1))
      write(iw,'("  y : ",a)') trim(labels(2))
      write(iw,'("  f : ",a)') trim(labels(3))
      write(iw,*)
      write(iw,'("ranges")')
      write(iw,'(es15.7, " <= x <= ",es15.7)') ranges(1, 1), ranges(2, 1)
      write(iw,'(es15.7, " <= y <= ",es15.7)') ranges(1, 2), ranges(2, 2)
      write(iw,'(es15.7, " <= z <= ",es15.7)') ranges(1, 3), ranges(2, 3)
      write(iw,*)
      write(iw,'("tics")')
      write(iw,'("  x : ",f15.7)') tics(1)
      write(iw,'("  y : ",f15.7)') tics(2)
      write(iw,'("  f : ",f15.7)') tics(3)
      write(iw,*)

      !iopt = get_opt(mode, CoMMode, ierr)
      !if (ierr /= 0) then
      !  write(iw,'("Read_Ctrl_Option> Error.")')
      !  write(iw,'("mode = ",a," is not available.")') trim(mode)
      !  stop
      !end if
      !option%mode          = iopt
      
      mpl2d%fpdf   = fpdf
      mpl2d%labels = labels
      mpl2d%ranges = ranges
      mpl2d%tics   = tics
      mpl2d%scales = scales

    end subroutine read_ctrl_matplotlib2d
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine setup_func3d(ng3, del, origin, func3d)
!-----------------------------------------------------------------------
      implicit none
 
      integer,        intent(in)  :: ng3(3)
      real(8),        intent(in)  :: del(3)
      real(8),        intent(in)  :: origin(3)
      type(s_func3d), intent(out) :: func3d
 
 
      func3d%ng3    = ng3
      func3d%del    = del
      func3d%origin = origin
      func3d%box(:) = ng3(:) * del(:)
      allocate(func3d%data(ng3(1), ng3(2), ng3(3)))
 
      func3d%data = 0.0d0
 
    end subroutine setup_func3d
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine setup_mesh3d(ng3, del, origin, mesh3d)
!-----------------------------------------------------------------------
      implicit none
 
      integer,        intent(in)  :: ng3(3)
      real(8),        intent(in)  :: del(3)
      real(8),        intent(in)  :: origin(3)
      type(s_mesh3d), intent(out) :: mesh3d
 
      integer :: igx, igy, igz
      real(8) :: x, y, z
 
 
      mesh3d%ng3    = ng3
      mesh3d%del    = del
      mesh3d%origin = origin
      mesh3d%box(:) = ng3(:) * del(:)
 
      allocate(mesh3d%data(3, ng3(1), ng3(2), ng3(3)))
 
      do igz = 1, ng3(3)
        z = origin(3) + (igz - 1) * del(3)
        do igy = 1, ng3(2)
          y = origin(2) + (igy - 1) * del(2)
          do igx = 1, ng3(1)
            x = origin(1) + (igx - 1) * del(1)
            mesh3d%data(1, igx, igy, igz) = x 
            mesh3d%data(2, igx, igy, igz) = y 
            mesh3d%data(3, igx, igy, igz) = z 
          end do
        end do
      end do
 
    end subroutine setup_mesh3d
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine search_nonzero_region(f, spreg, threshold)
!-----------------------------------------------------------------------
      implicit none

      type(s_func3d),    intent(in)    :: f
      type(s_spregion),  intent(inout) :: spreg
      real(8), optional, intent(in)    :: threshold

      integer :: ix, iy, iz, ispr
      integer :: ngx, ngy, ngz, nspr
      real(8) :: val, x, y, z
      real(8) :: thr


      if (present(threshold)) then
        thr = threshold
      else
        thr = 1.0d-5
      end if

      if (allocated(spreg%data)) &
        deallocate(spreg%data)

      ngx = f%ng3(1)
      ngy = f%ng3(2)
      ngz = f%ng3(3)

      spreg%del    = f%del
      spreg%origin = f%origin
      spreg%box    = f%box

      nspr = 0
      do iz = 1, ngz
        do iy = 1, ngy
          do ix = 1, ngx
            val = f%data(ix, iy, iz)

            if (abs(val) >= thr) & 
              nspr = nspr + 1

          end do
        end do
      end do
   
      spreg%nspr = nspr
      allocate(spreg%data(1:3, nspr))

      ispr = 0
      do iz = 1, ngz
        do iy = 1, ngy
          do ix = 1, ngx
            val = f%data(ix, iy, iz)

            if (abs(val) >= 1.0d-5) then 
              ispr                = ispr + 1
              x                   = f%del(1) * (ix - 1) + f%origin(1)
              y                   = f%del(2) * (iy - 1) + f%origin(2)
              z                   = f%del(3) * (iz - 1) + f%origin(3)
              spreg%data(1, ispr) = x
              spreg%data(2, ispr) = y
              spreg%data(3, ispr) = z
            end if

          end do
        end do
      end do


    end subroutine search_nonzero_region
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine generate_1dmap_gnuplot(fhead, func3d, func_err)
!-----------------------------------------------------------------------
      implicit none
     
      character(len=MaxChar), intent(in) :: fhead
      type(s_func3d),         intent(in) :: func3d
      real(8), optional,      intent(in) :: func_err(:)
     
      character(len=MaxChar) :: fg
      integer                :: ix
      real(8)                :: x
     
     
      write(fg,'(a,".gnplt")') trim(fhead)
     
      open(UnitGNPLT,file=trim(fg))

        if (present(func_err)) then
          do ix = 1, func3d%ng3(1)
            x = func3d%origin(1) + (ix - 1) * func3d%del(1)
            write(UnitGNPLT,'(3(es20.10,2x))') x, func3d%data(ix, 1, 1), func_err(ix) 
          end do
        else
          do ix = 1, func3d%ng3(1)
            x = func3d%origin(1) + (ix - 1) * func3d%del(1)
            write(UnitGNPLT,'(2(es20.10,2x))') x, func3d%data(ix, 1, 1) 
          end do
        end if
      close(UnitGNPLT)

   end subroutine generate_1dmap_gnuplot
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine generate_2dmap_matplotlib(fhead, func3d)
!-----------------------------------------------------------------------
      implicit none
 
      character(len=MaxChar), intent(in) :: fhead
      type(s_func3d),         intent(in) :: func3d
 
      character(len=MaxChar) :: fxi, fyi, fg
      integer                :: ix, iy
      real(8)                :: x, y
 
      write(fxi,'(a,".xi")') trim(fhead)
      write(fyi,'(a,".yi")') trim(fhead)
      write(fg,'(a,".mpl")') trim(fhead)
 
      open(UnitMPL,file=trim(fxi))
        do ix = 1, func3d%ng3(1)
          x = func3d%origin(1) + (ix - 1) * func3d%del(1)
          write(UnitMPL,'(f20.10)',advance='no') x 
        end do
        write(UnitMPL,*) 
      close(UnitMPL)
 
      open(UnitMPL,file=trim(fyi))
        do iy = 1, func3d%ng3(2)
          y = func3d%origin(2) + (iy - 1) * func3d%del(2)
          write(UnitMPL,'(f20.10)',advance='no') y 
        end do
        write(UnitMPL,*) 
      close(UnitMPL)
 
      open(UnitMPL,file=trim(fg))
        do iy = 1, func3d%ng3(2)
          do ix = 1, func3d%ng3(1)
             write(UnitMPL,'(f20.10)',advance='no') func3d%data(ix, iy, 1)
          end do
          write(UnitMPL,*) 
        end do
      close(UnitMPL)
 
    end subroutine generate_2dmap_matplotlib
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine generate_2dmap_gnuplot(fhead, func3d)
!-----------------------------------------------------------------------
      implicit none
     
      character(len=MaxChar), intent(in) :: fhead
      type(s_func3d),         intent(in) :: func3d
     
      character(len=MaxChar) :: fg
      integer                :: ix, iy
      real(8)                :: x, y
     
     
      write(fg,'(a,".gnplt")') trim(fhead)
     
      open(UnitGNPLT,file=trim(fg))
        do ix = 1, func3d%ng3(1)
          x = func3d%origin(1) + (ix - 1) * func3d%del(1)
          do iy = 1, func3d%ng3(2)
            y = func3d%origin(2) + (iy - 1) * func3d%del(2)
             write(UnitGNPLT,'(3(es20.10,2x))') x, y, func3d%data(ix, iy, 1) 
          end do
          write(UnitGNPLT,*) 
        end do
      close(UnitGNPLT)

   end subroutine generate_2dmap_gnuplot
!-----------------------------------------------------------------------
!
!------------------------------------------------------------------------
    subroutine write_dx(fhead, func3d)
!------------------------------------------------------------------------
      implicit none

      character(len=MaxChar), intent(in) :: fhead
      type(s_func3d),         intent(in) :: func3d 

      character(len=MaxChar) :: fg 
      integer                :: n, ngr 
      integer                :: cyl, rem, lin, ati, i, j
      integer                :: igx, igy, igz, ig
      integer                :: ngx, ngy, ngz

      real(8), allocatable   :: f3d(:)


      write(fg,'(a,".dx")') trim(fhead)

      ngr = func3d%ng3(1) * func3d%ng3(2) * func3d%ng3(3)
      ngx = func3d%ng3(1)
      ngy = func3d%ng3(2)
      ngz = func3d%ng3(3)

      allocate(f3d(ngr))

      do igz = 1, ngz 
        do igy = 1, ngy 
          do igx = 1, ngx
            !ig      = igx + (igy -1) * ngx + (igz - 1) * ngx * ngy
            ig      = igz + (igy -1) * ngz + (igx - 1) * ngy * ngz
            f3d(ig) = func3d%data(igx, igy, igz)
          end do
        end do
      end do

      open(UnitDX, file = trim(fg))
        write(UnitDX, "(A,3(i0,2x))") "object 1 class gridpositions counts ", func3d%ng3(1), func3d%ng3(2), func3d%ng3(3)
        write(UnitDX, "(A,3(e13.6,3x))") "origin ", func3d%origin(1), func3d%origin(2), func3d%origin(3)
        write(UnitDX, "(A,3(e13.6,3x))") "delta ", func3d%del(1), 0.0d0, 0.0d0
        write(UnitDX, "(A,3(e13.6,3x))") "delta ", 0.0d0, func3d%del(2), 0.0d0
        write(UnitDX, "(A,3(e13.6,3x))") "delta ", 0.0d0, 0.0d0, func3d%del(3)
        write(UnitDX, "(A,3(i0,2x))") "object 2 class gridpositions counts ", func3d%ng3(1), func3d%ng3(2), func3d%ng3(3)
        write(UnitDX, "(A,i0,A)") "object 3 class array type double rank 0 items ", ngr, " data follows"
        n = 3
        cyl = ngr / n
        rem = ngr - cyl * n
        lin = 0
        do ati = 1, cyl
           write(UnitDX, 753) (f3d(lin + j), j = 1, n)
           lin = lin + n
        end do
        if(rem /= 0) then
           write(UnitDX, 753) (f3d(lin + j), j = 1, rem)
        endif
753     format(e13.6,2e14.6)

        write(UnitDX, '(A)') 'attribute "dep" string "positions"'
        write(UnitDX, '(A)') 'object "regular positions regular connections" class field'
        write(UnitDX, '(A)') 'component "positions" value 1'
        write(UnitDX, '(A)') 'component "connections" value 2'
        write(UnitDX, '(A)') 'component "data" value 3'
      close(UnitDX)

      deallocate(f3d)

    end subroutine write_dx
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
    subroutine read_dx(fname, f)
!------------------------------------------------------------------------
      implicit none

      character(len=MaxChar), intent(in)  :: fname
      type(s_func3d),         intent(out) :: f

      integer           :: iunit
      integer           :: n, ngr
      integer           :: i, j, cyl, rem, lin, ati
      integer           :: igx, igy, igz, ig
      integer           :: ngx, ngy, ngz
      character(len=10) :: cdum1, cdum2, cdum3, cdum4, cdum5
      logical           :: chk

      real(8), allocatable :: wrk(:)


      iunit = 10
      do i = 1, 99
        inquire(unit=iunit, opened=chk)
        if (chk) then
          iunit = iunit + 1
        else
          exit
        end if
      end do

      open(unit = iunit, file = trim(fname))
        read(iunit,*) cdum1, cdum2, cdum3, cdum4, cdum5, f%ng3(1), f%ng3(2), f%ng3(3) 
        read(iunit,*) cdum1, f%origin(1), f%origin(2), f%origin(3)
        read(iunit,*) cdum1, f%del(1), cdum2,     cdum3
        read(iunit,*) cdum1, cdum2,     f%del(2), cdum3
        read(iunit,*) cdum1, cdum2,     cdum3,     f%del(3)
        read(iunit,*) cdum1, cdum2, cdum3, cdum4, cdum5, f%ng3(1), f%ng3(2), f%ng3(3)
        read(iunit,*)

        ngx      = f%ng3(1)
        ngy      = f%ng3(2)
        ngz      = f%ng3(3)
        ngr      = f%ng3(1) * f%ng3(2) * f%ng3(3)

        f%box(:) = f%ng3(:) * f%del(:)

        ! allocate memory
        ! - global
        allocate(f%data(f%ng3(1), f%ng3(2), f%ng3(3)))
        ! - local
        allocate(wrk(ngr))

        n = 3
        cyl = ngr / n
        rem = ngr - cyl * n
        lin = 0
        do ati = 1, cyl
          read(iunit,*) (wrk(lin + j), j = 1, n)
          lin = lin + n
        end do
        if(rem /= 0) then
          read(iunit,*) (wrk(lin + j), j = 1, rem)
        endif
      close(iunit)


      do igz = 1, ngz 
        do igy = 1, ngy 
          do igx = 1, ngx
            ig      = igz + (igy -1) * ngz + (igx - 1) * ngy * ngz
            f%data(igx, igy, igz) = wrk(ig)
          end do
        end do
      end do

      deallocate(wrk)
!
    end subroutine read_dx
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
    subroutine generate_script_matplotlib2d(fhead, mpl2d) 
!------------------------------------------------------------------------
      implicit none

      character(len=MaxChar), intent(in) :: fhead
      type(s_mpl2d),          intent(in) :: mpl2d

      character(len=MaxChar) :: fs, fx, fy, ff
      character(len=MaxChar) :: cmd


      write(fx,'(a,".xi")')  trim(fhead)
      write(fy,'(a,".yi")')  trim(fhead)
      write(ff,'(a,".mpl")') trim(fhead)

      write(fs,'(a,"_mpl2d.py")') trim(fhead)
      open(UnitMPL, file=trim(fs))
        write(UnitMPL, '("#!/usr/bin/env python3")')
        write(UnitMPL,*)
        write(UnitMPL, '("from pylab import meshgrid, cm, imshow, contour,\")')
        write(UnitMPL, '("                  clabel, colorbar, axis, title,\")')
        write(UnitMPL, '("                  show, pcolor")')
        write(UnitMPL, '("import numpy as np")')
        write(UnitMPL, '("import matplotlib.pyplot as plt")')
        write(UnitMPL,*)

        write(UnitMPL, '("plt.switch_backend(''agg'')")')
        write(UnitMPL, '("plt.style.use(''seaborn-bright'')")')
        write(UnitMPL, '("")')
        write(UnitMPL, '("plt.rcParams[''font.family'']         = ''serif''")')
        write(UnitMPL, '("plt.rcParams[''font.serif'']          = ''Times New Roman''")')
        write(UnitMPL, '("plt.rcParams[''font.size'']           = 12")')
        write(UnitMPL, '("plt.rcParams[''axes.labelsize'']      = 14")')
        write(UnitMPL, '("plt.rcParams[''mathtext.cal'']        = ''serif''")')
        write(UnitMPL, '("plt.rcParams[''mathtext.rm'']         = ''serif''")')
        write(UnitMPL, '("plt.rcParams[''mathtext.it'']         = ''serif:italic''")')
        write(UnitMPL, '("plt.rcParams[''mathtext.bf'']         = ''serif:bold''")')
        write(UnitMPL, '("plt.rcParams[''mathtext.fontset'']    = ''cm''")')
        write(UnitMPL, '("")')
        write(UnitMPL, '("plt.rcParams[''pdf.fonttype'']        = 42")')
        write(UnitMPL, '("plt.rcParams[''ps.fonttype'']         = 42")')
        write(UnitMPL, '("")')
        write(UnitMPL, '("plt.rcParams[''axes.axisbelow'']      = True")')
        write(UnitMPL, '("")')
        write(UnitMPL, '("plt.rcParams[''xtick.direction'']     = ''in''")')
        write(UnitMPL, '("plt.rcParams[''ytick.direction'']     = ''in''")')
        write(UnitMPL, '("plt.rcParams[''xtick.minor.visible''] = True")')
        write(UnitMPL, '("plt.rcParams[''ytick.minor.visible''] = True")')
        write(UnitMPL, '("")')
        write(UnitMPL, '("plt.figure(figsize=(4,3))")')
        write(UnitMPL, '("plt.grid(color=''gray'',linestyle=''dotted'', linewidth=1)")')
        write(UnitMPL, '("x = np.loadtxt(''",a,"'')")') trim(fx)
        write(UnitMPL, '("y = np.loadtxt(''",a,"'')")') trim(fy)
        write(UnitMPL, '("z = np.loadtxt(''",a,"'')")') trim(ff)
        write(UnitMPL, '("X, Y = meshgrid(x, y)")')
        write(UnitMPL, '("interval = np.arange(",f15.7,",",f15.7,",",f15.7,")")') &
          mpl2d%ranges(1,3), mpl2d%ranges(2, 3), mpl2d%tics(3)
        write(UnitMPL, '("c1 = plt.contourf(X, Y, z, interval, cmap=cm.jet)")')
        write(UnitMPL, '("c2 = plt.contour(c1, interval, colors=''k'', linewidths=0.3)")')
        write(UnitMPL, '("")')
        write(UnitMPL, '("cbar = plt.colorbar(c1)")')
        write(UnitMPL, '("cbar.set_label(''",a,"'')")') trim(mpl2d%labels(3))
        write(UnitMPL, '("")')
        write(UnitMPL, '("plt.xlim([", f15.7, ",", f15.7, "])")') &
          mpl2d%ranges(1, 1), mpl2d%ranges(2, 1)
        write(UnitMPL, '("plt.ylim([", f15.7, ",", f15.7, "])")') &
          mpl2d%ranges(1, 2), mpl2d%ranges(2, 2)
        write(UnitMPL, '("")')
        write(UnitMPL, '("plt.xticks(np.arange(", f15.7, ",", f15.7, ",", f15.7,"))")') &
          mpl2d%ranges(1, 1), mpl2d%ranges(2, 1) + mpl2d%tics(1), mpl2d%tics(1)
        write(UnitMPL, '("plt.yticks(np.arange(", f15.7, ",", f15.7, ",", f15.7,"))")') &
          mpl2d%ranges(1, 2), mpl2d%ranges(2, 2) + mpl2d%tics(2), mpl2d%tics(2)
        write(UnitMPL, '("plt.xlabel(''",a,"'')")') trim(mpl2d%labels(1))
        write(UnitMPL, '("plt.ylabel(''",a,"'')")') trim(mpl2d%labels(2))

        write(UnitMPL, '("")')
        write(UnitMPL, '("plt.tick_params(left=True, right=True,")')
        write(UnitMPL, '("                top=True,  bottom=True,")')
        write(UnitMPL, '("                labelleft=True,")')
        write(UnitMPL, '("                labelright=False)")')
        write(UnitMPL, '("plt.tight_layout()")')
        write(UnitMPL, '("plt.savefig(''",a,"'')")') trim(mpl2d%fpdf)

      close(UnitMPL)


      write(cmd, '("chmod +x ",a)') trim(fs)
      call system(trim(cmd))

    end subroutine generate_script_matplotlib2d
!------------------------------------------------------------------------

end module mod_grid3d
!=======================================================================

