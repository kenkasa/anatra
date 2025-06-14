!-------------------------------------------------------------------------------
!
!  Module   mod_voronoi  
!
!    Module for voronoi tessellation 
!
!    (c) Copyright 2021 Osaka Univ. All rights reserved.
!
!-------------------------------------------------------------------------------

module mod_voronoi

  use mod_util
  use mod_const
  use mod_grid3d

  implicit none

  ! constants
  !
  integer,      parameter, public :: VrCellTypeMIN  = 1
  integer,      parameter, public :: VrCellTypeMAX  = 2

  integer :: MaxNcell = 100

  ! structures
  !
  type :: s_voronoi
    integer              :: ndim          = 2
    integer              :: ncell         = 1.0d0
    real(8)              :: normvec(3)    = 1.0d0
    real(8)              :: separation(3) = 0.0d0
    real(8), allocatable :: cellpos(:, :)
    real(8), allocatable :: cellval(:)
    type(s_func3d)       :: states

    integer              :: nbound_points = 0
    real(8), allocatable :: boundary(:, :)
  end type s_voronoi

  ! subroutines
  !
  public :: setup_voronoi 
  public :: search_cellpos
  public :: voronoi_boundary
  public :: voronoi_state

  contains

    !---------------------------------------------------------------------------
    !
    !  Subroutine    setup_voronoi 
    !  @brief        setup voronoi tessellation  
    !  @authors      KK
    !  @param[in]    ndim       : number of dimensions 
    !  @param[in]    ncell      : number of cells
    !  @param[in]    normvec    : normalization vector 
    !  @param[inout] voronoi    : structure of voronoi 
    !  @param[in]    funcin     : input 2d-function (OPTIONAL) 
    !  @param[in]    separation : separation between cells (OPTIONAL)
    !
    !---------------------------------------------------------------------------

    subroutine setup_voronoi(ndim, ncell, normvec, voronoi, funcin, separation)

      implicit none

      integer, parameter                      :: up_resol = 5

      integer,                  intent(in)    :: ndim
      integer,                  intent(in)    :: ncell
      real(8),                  intent(in)    :: normvec(3)
      type(s_voronoi),          intent(inout) :: voronoi
      type(s_func3d), optional, intent(in)    :: funcin
      real(8),        optional, intent(in)    :: separation(3)


      voronoi%ndim       = ndim
      voronoi%ncell      = ncell
      voronoi%normvec    = normvec

      if (present(separation)) then
        voronoi%separation = separation
      end if

      if (present(funcin)) then
        call setup_func3d(funcin%ng3 * up_resol,        &
                          funcin%del / float(up_resol), &
                          funcin%origin,                &
                          voronoi%states)
      end if

      allocate(voronoi%cellpos(3, ncell), voronoi%cellval(ncell))
      voronoi%cellpos = 0.0d0
      voronoi%cellval = 0.0d0

    end subroutine setup_voronoi 

    !---------------------------------------------------------------------------
    !
    !  Subroutine  serach_cellpos 
    !  @brief      search cell position  
    !  @authors    KK
    !  @param[in]    celltype : definition of cell (See VrCellTypeXX) 
    !  @param[in]    func     : input of 2d-function 
    !  @param[inout] voronoi  : structure of voronoi 
    !
    !---------------------------------------------------------------------------

    subroutine search_cellpos(celltype, func, voronoi)

      implicit none

      integer,           intent(in)    :: celltype 
      type(s_func3d),    intent(in)    :: func
      type(s_voronoi),   intent(inout) :: voronoi

      integer              :: ic, jc, igx, igy, igz
      integer              :: ndim, ncell, ngx, ngy, ngz, ndetect
      real(8)              :: fmin, fval, fn, fs, fe, fw, fu, fd
      real(8)              :: x, y, z, dx, dy, dz, xmin, ymin, zmin
      real(8)              :: origin(3), del(3), sep(3)
      logical              :: is_within
   
      real(8), allocatable :: f(:, :, :)
      real(8), allocatable :: xyzmin(:, :)


      ! setup variables
      !
      ndim        = voronoi%ndim
      ncell       = voronoi%ncell
      ngx         = func%ng3(1)
      ngy         = func%ng3(2)
      ngz         = func%ng3(3)
      origin(1:3) = func%origin(1:3)
      del(1:3)    = func%del(1:3)
      sep(1:3)    = voronoi%separation(1:3)

      ! allocate memory
      !
      allocate(f(ngx, ngy, ngz), xyzmin(3, ncell))

     
      ! change sign of function according to celltype
      !
      if (celltype == VrCellTypeMIN) then
        f(1:ngx, 1:ngy, 1:ngz) =   func%data(1:ngx, 1:ngy, 1:ngz)
      else if (celltype == VrCellTypeMAX) then
        f(1:ngx, 1:ngy, 1:ngz) = - func%data(1:ngx, 1:ngy, 1:ngz)
      end if

      if (ndim == 2) then

        ! search global minimum
        !
        fmin  = 1.0d30
        xyzmin = 0.0d0
        do igy = 2, ngy - 1
          do igx = 2, ngx - 1
            x = origin(1) + del(1) * (igx - 1)
            y = origin(2) + del(2) * (igy - 1)
            
            fval = f(igx, igy, 1)
            if (fmin > fval) then
              fmin = fval
              xmin = x
              ymin = y
            end if

          end do
        end do
       
        xyzmin(1, 1)        = xmin
        xyzmin(2, 1)        = ymin
        voronoi%cellval(1)  = fmin
       
        ! search local minimum
        !
        ndetect = 1
        do ic = 2, ncell
          fmin = 1.0d30
          do igy = 2, ngy - 1
            do igx = 2, ngx - 1
              x    = origin(1) + del(1) * (igx - 1)
              y    = origin(2) + del(2) * (igy - 1)
              fval = f(igx, igy, 1) 
              
              is_within = .false.
              do jc = 1, ndetect
                dx = abs(x - xyzmin(1, jc))
                dy = abs(y - xyzmin(2, jc))
              
                if (dx < sep(1) .and. dy < sep(2)) then
                  is_within = .true.
                end if
              
              end do
              
              if (.not. is_within) then
                if (fmin > fval) then
                  fn   = f(igx    , igy + 1, igz)
                  fs   = f(igx    , igy - 1, igz)
                  fe   = f(igx + 1, igy    , igz)
                  fw   = f(igx - 1, igy    , igz)
              
                  if (fval <= fn .and. fval <= fs .and. fval <= fe &
                      .and. fval <= fw) then 
                    fmin = fval
                    xmin = x
                    ymin = y
                  end if
                end if
              end if
       
       
            end do
          end do
       
          xyzmin(1, ic)        = xmin
          xyzmin(2, ic)        = ymin
          voronoi%cellval(ic)  = fmin
       
          ndetect             = ndetect + 1
        end do
       
        voronoi%cellpos = xyzmin
       
        ! print out results
        !
        write(iw,'("Search_CellPos> Detected minimum (x, y, fmin)")')
        do ic = 1, ncell
          write(iw,'(i3,2x,3f20.10)')                           &
            ic, voronoi%cellpos(1, ic), voronoi%cellpos(2, ic), &
            voronoi%cellval(ic)
        end do

      else if (ndim == 3) then
        fmin  = 1.0d30
        xyzmin = 0.0d0
        do igz = 2, ngz - 1
          do igy = 2, ngy - 1
            do igx = 2, ngx - 1
              x = origin(1) + del(1) * (igx - 1)
              y = origin(2) + del(2) * (igy - 1)
              z = origin(3) + del(3) * (igz - 1)
             
              fval = f(igx, igy, igz)
              if (fmin > fval) then
                fmin = fval
                xmin = x
                ymin = y
                zmin = z
              end if
       
            end do
          end do
        end do
       
        xyzmin(1, 1)        = xmin
        xyzmin(2, 1)        = ymin
        xyzmin(3, 1)        = zmin
        voronoi%cellval(1)  = fmin
       
        ! search local minimum
        !
        ndetect = 1
        do ic = 2, ncell
          fmin = 1.0d30
          do igz = 2, ngz - 1
            do igy = 2, ngy - 1
              do igx = 2, ngx - 1
                x    = origin(1) + del(1) * (igx - 1)
                y    = origin(2) + del(2) * (igy - 1)
                z    = origin(3) + del(3) * (igz - 1)
                fval = f(igx, igy, igz) 
               
                is_within = .false.
                do jc = 1, ndetect
                  dx = abs(x - xyzmin(1, jc))
                  dy = abs(y - xyzmin(2, jc))
                  dz = abs(z - xyzmin(3, jc))
               
                  if (dx < sep(1) .and. dy < sep(2) .and. dz < sep(3)) then
                    is_within = .true.
                  end if
               
                end do
               
                if (.not. is_within) then
                  if (fmin > fval) then
                    fn   = f(igx    , igy + 1, igz)
                    fs   = f(igx    , igy - 1, igz)
                    fe   = f(igx + 1, igy    , igz)
                    fw   = f(igx - 1, igy    , igz)
                    fu   = f(igx,     igy    , igz + 1)
                    fd   = f(igx,     igy    , igz - 1)
               
                    if (fval <= fn .and. fval <= fs .and. fval <= fe &
                        .and. fval <= fw .and. fval <= fu .and. fval <= fd) then 
                      fmin = fval
                      xmin = x
                      ymin = y
                      zmin = z
                    end if
                  end if
                end if
       
              end do
       
            end do
          end do
       
          xyzmin(1, ic)        = xmin
          xyzmin(2, ic)        = ymin
          xyzmin(3, ic)        = zmin
          voronoi%cellval(ic)  = fmin
       
          ndetect             = ndetect + 1
        end do
       
        voronoi%cellpos = xyzmin
       
        ! print out results
        !
        write(iw,'("Search_CellPos> Detected minimum (x, y, z, fmin)")')
        do ic = 1, ncell
          write(iw,'(i3,2x,4f20.10)')                                                   &
            ic, voronoi%cellpos(1, ic), voronoi%cellpos(2, ic), voronoi%cellpos(3, ic), &
            voronoi%cellval(ic)
        end do
      end if

      ! deallocate memory
      !
      deallocate(f, xyzmin)

    end subroutine search_cellpos 

    !---------------------------------------------------------------------------
    !
    !  Subroutine  voronoi_boundary 
    !  @brief      deterimine voronoi boundary 
    !  @authors    KK
    !  @param[inout] voronoi  : structure of voronoi 
    !
    !---------------------------------------------------------------------------

    subroutine voronoi_boundary(voronoi)

      implicit none

      type(s_voronoi),   intent(inout) :: voronoi

      integer              :: ib, ic, jc, igx, igy
      integer              :: ncell, ngx, ngy, curr, prev 
      real(8)              :: x, y, dx, dy, xmin, ymin
      real(8)              :: origin(2), del(2)
      logical              :: is_within

      integer, allocatable :: bound_label(:, :)
   

      ! setup variables
      !
      ncell       = voronoi%ncell
      ngx         = voronoi%states%ng3(1)     
      ngy         = voronoi%states%ng3(2)     
      origin(1:2) = voronoi%states%origin(1:2)
      del(1:2)    = voronoi%states%del(1:2)   

      ! allocate memory
      !
      allocate(bound_label(ngx, ngy))

      ! determine state of each grid
      !
      voronoi%states%data = 0.0d0
      do igy = 1, ngy
        do igx = 1, ngx
          x = origin(1) + del(1) * (igx - 1)
          y = origin(2) + del(2) * (igy - 1)
          
          voronoi%states%data(igx, igy, 1) &
            = float(voronoi_state(voronoi, x, y))
          !write(iw,'("state : ", 3f20.10)') x, y, voronoi%states%data(igx, igy, 1)
        end do
      end do

      ! determine boundary
      !
      bound_label = 0

      ! sweep along x direction
      do igy = 1, ngy

        prev = nint(voronoi%states%data(1, igy, 1))
        do igx = 2, ngx
          curr = nint(voronoi%states%data(igx, igy, 1))

          if (curr /= prev) then
            bound_label(igx, igy) = 1
          end if

          prev = curr
        end do
      end do

      ! sweep along y direction
      do igx = 1, ngx

        prev = nint(voronoi%states%data(igx, 1, 1))
        do igy = 2, ngy
          curr = nint(voronoi%states%data(igx, igy, 1))

          if (curr /= prev) then
            bound_label(igx, igy) = 1
          end if

          prev = curr
        end do
      end do

      voronoi%nbound_points = sum(bound_label)
      allocate(voronoi%boundary(2, voronoi%nbound_points))

      ib = 0
      do igy = 1, ngy
        do igx = 1, ngx
          x = origin(1) + del(1) * (igx - 1)
          y = origin(2) + del(2) * (igy - 1)
          
          if (bound_label(igx, igy) == 1) then
            ib = ib + 1
            voronoi%boundary(1, ib) = x
            voronoi%boundary(2, ib) = y
          end if 
        end do
      end do


      ! print out results
      !
      !write(iw,'("Voronoi_Boundary> Boundary (x, y)")')
      !do igy = 1, ngy
      !  do igx = 1, ngx
      !    x = origin(1) + del(1) * (igx - 1)
      !    y = origin(2) + del(2) * (igy - 1)
      !    
      !    if (bound_label(igx, igy) == 1) then
      !      write(iw,'(2f20.10)') x, y
      !    end if 
      !  end do
      !end do

      ! deallocate memory
      !
      deallocate(bound_label)

    end subroutine voronoi_boundary 

    !---------------------------------------------------------------------------
    !
    !  Function    voronoi_state
    !  @brief      return voronoi state 
    !  @authors    KK
    !  @param[in]  voronoi  : structure of voronoi 
    !  @param[in]  x        : x coordinate 
    !  @param[in]  y        : y coordinate
    !  @param[in]  z        : z coordinate
    !
    !---------------------------------------------------------------------------

    function voronoi_state(voronoi, x, y, z)

      implicit none

      type(s_voronoi), intent(in) :: voronoi
      real(8),         intent(in) :: x
      real(8),         intent(in) :: y
      real(8),         intent(in), optional :: z 

      ! return value
      integer :: voronoi_state 

      ! local variables
      integer :: ic
      real(8) :: dx, dy, dz, rsq(MaxNcell)


      if (present(z)) then

        rsq = 0.0d0
        do ic = 1, voronoi%ncell
          dx      = (x - voronoi%cellpos(1, ic)) / voronoi%normvec(1)
          dy      = (y - voronoi%cellpos(2, ic)) / voronoi%normvec(2)
          dz      = (z - voronoi%cellpos(3, ic)) / voronoi%normvec(3)
          rsq(ic) = dx * dx + dy * dy + dz * dz
        end do

      else
        rsq = 0.0d0
        do ic = 1, voronoi%ncell
          dx      = (x - voronoi%cellpos(1, ic)) / voronoi%normvec(1)
          dy      = (y - voronoi%cellpos(2, ic)) / voronoi%normvec(2)
       
          rsq(ic) = dx * dx + dy * dy
         end do
      end if

      voronoi_state = minloc(rsq(1:voronoi%ncell), 1)

    end function voronoi_state

end module
!=======================================================================
