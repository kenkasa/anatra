!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_akima
  use mod_grid3d
  use mod_voronoi
  use mod_input
  use mod_output
  use mod_ctrl
  use mod_cv
  use mod_bootstrap
  use mod_random

  ! constants
  !
  real(8), parameter, private :: EPS = 1.0d-8 

  ! structures
  !
  type :: s_booteach
    integer :: ntrial
    integer :: nsample

    integer,        allocatable :: rand(:, :)
    type(s_func3d), allocatable :: g0(:)
    integer,        allocatable :: state_ratio(:)

  end type s_booteach

  type :: s_bootave
    type(s_func3d) :: g0
    real(8), allocatable :: yerr(:)
    real(8) :: Gint
    real(8) :: Gint_err
  end type s_bootave

  ! subroutines
  !
  public :: analyze 

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option, cvinfo, mpl2d, bootopt)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),   intent(in)    :: input
      type(s_output),  intent(in)    :: output
      type(s_option),  intent(in)    :: option
      type(s_cvinfo),  intent(in)    :: cvinfo
      type(s_mpl2d),   intent(in)    :: mpl2d
      type(s_bootopt), intent(inout) :: bootopt

      type(s_func3d)          :: g0,  g1 
      type(s_func3d)          :: pmf0, pmf1
      type(s_voronoi)         :: voronoi

      type(s_cv), allocatable :: cv(:)
      real(8),    allocatable :: state_count(:), fe_state(:)
      real(8),    allocatable :: fe_ave(:), fe_err(:)
      real(8),    allocatable :: pmf_err(:)

      integer                 :: colmax, nstep_tot
      integer                 :: ng3_spl(3)
      real(8)                 :: del_spl(3)
                              
      integer                 :: ifile, iunit, istep, icv
      integer                 :: ib, ic, jdim, igx, igy, igz, igr, istate, nstate, active
      integer                 :: igx_low, igx_up, igx_low_ref, igx_up_ref
      integer                 :: gind(MaxDim) 
      real(8)                 :: dv, dr, r, val, vol, gmax, gmax0, gmax1
      real(8)                 :: gval, crd(3), rcmin, rcmax, sref, sval 
      real(8)                 :: pmin, kT, fourpi, fact, weight_sum
      real(8)                 :: Gint, Gint_err, Gint_ref
      logical                 :: is_inside
      character(len=MaxChar)  :: fhead_out  


      ! generate matplotlib python script 
      !
      if (option%ndim == 2 .and. option%gen_script_mpl2d) then

        write(fhead_out,'(a,"_original")') trim(output%fhead)
        call generate_script_matplotlib2d(output%fhead, mpl2d)

        if (option%use_spline) then
          write(fhead_out,'(a,"_spline")') trim(output%fhead)
          call generate_script_matplotlib2d(fhead_out, mpl2d)
        end if

      end if

      if (option%skip_calc) & 
        return

      ! setup
      !
      allocate(cv(cvinfo%nfile))
      allocate(pmf_err(1:option%ng3(1)))

      pmf_err = 0.0d0

      dv = 1.0d0
      do icv = 1, option%ndim
        dv = dv * option%del(icv) 
      end do

      colmax = maxval(option%xyzcol(1:option%ndim))

      kT = Boltz * option%temperature

      if (option%use_spline) then
        ng3_spl = 1
        del_spl = 1.0d0
        ng3_spl(1:option%ndim) = option%ng3(1:option%ndim) &
                                 * option%spline_resolution
        del_spl(1:option%ndim) = option%del(1:option%ndim) &
                                 / option%spline_resolution 
      end if 

      if (option%calc_statepop) then
        nstate = option%nstate
        allocate(state_count(nstate))
        allocate(fe_state(nstate))
        allocate(fe_ave(nstate), fe_err(nstate))
        state_count = 0
        fe_state    = 0.0d0
        fe_ave      = 0.0d0
        fe_err      = 0.0d0
      end if

      ! read cv files
      !
      write(iw,*)
      write(iw,'("Analyze> Read CV file")')
      do ifile = 1, cvinfo%nfile
        call read_cv(cvinfo%fcv(ifile), colmax, cv(ifile))
      end do
      write(iw,'(">> Finished")')

      if (trim(input%flist_weight) /= "") then
        do ifile = 1, cvinfo%nfile
          call read_cv_weight(cvinfo%fweight(ifile), cv(ifile))
        end do
      else
        do ifile = 1, cvinfo%nfile
          istep = cv(ifile)%nstep
          allocate(cv(ifile)%weight(istep))
          cv(ifile)%weight = 1.0d0
        end do
      end if

      nstep_tot = 0
      do ifile = 1, cvinfo%nfile
        nstep_tot = nstep_tot + (cv(ifile)%nstep - option%nstart + 1) 
      end do

      write(iw,*)
      write(iw,'("Analyze> Calculate Histogram")')

      ! setup 3d-functions
      !
      call setup_func3d(option%ng3, option%del, option%origin, g0) 
      call setup_func3d(option%ng3, option%del, option%origin, pmf0)
      if (option%use_spline) then
        call setup_func3d(ng3_spl, del_spl, option%origin, g1) 
        call setup_func3d(ng3_spl, del_spl, option%origin, pmf1)
      end if

      ! generate histogram
      !
      if (option%use_bootstrap) then
        call analyze_bootstrap(option, output, bootopt, cv, cvinfo%nfile, g0, pmf_err, fe_ave, fe_err, Gint, Gint_err)
      else

        weight_sum = 0.0d0
        do ifile = 1, cvinfo%nfile
          do istep = option%nstart, cv(ifile)%nstep
            is_inside = .true.
            do icv = 1, option%ndim
              val       = cv(ifile)%data(option%xyzcol(icv),istep)
              gind(icv) = (val - option%origin(icv)) / option%del(icv) + 1
              if (gind(icv) < 1 .or. gind(icv) > option%ng3(icv)) then
                is_inside = .false.
              end if
            end do

            if (option%calc_statepop) then
              do istate = 1, nstate
                active = 0
                do icv = 1, option%ndim
                  rcmin = option%staterange(1, icv, istate)
                  rcmax = option%staterange(2, icv, istate)
                  val = cv(ifile)%data(option%xyzcol(icv),istep)

                  if (val >= rcmin .and. val < rcmax) then
                    active = active + 1
                  end if
                end do

                if (active == option%ndim) then
                  state_count(istate) = state_count(istate) + cv(ifile)%weight(istep)
                end if

              end do
            end if
       
            if (option%ndim == 1) then
              gind(2:3) = 1
            else if (option%ndim == 2) then
              gind(3)   = 1
            end if
       
            if (is_inside) then
              !g0%data(gind(1), gind(2), gind(3)) = g0%data(gind(1), gind(2), gind(3)) + 1.0d0
              g0%data(gind(1), gind(2), gind(3)) = g0%data(gind(1), gind(2), gind(3)) &
                + cv(ifile)%weight(istep)
              weight_sum = weight_sum + cv(ifile)%weight(istep)
            end if
       
          end do 
        end do

        if (option%ndim < 3) then
          !g0%data = g0%data / (nstep_tot * dv)
          g0%data = g0%data / (weight_sum * dv)
        else
          if (option%space3D) then
            vol     = option%box_system(1) * option%box_system(2) * option%box_system(3)
            g0%data = vol * g0%data / (option%nstep_system * dv) 
          else
            !g0%data = g0%data / (nstep_tot * dv)
            g0%data = g0%data / (weight_sum * dv)
          end if
        end if


        if (option%calc_gintegral .and. option%ndim == 1) then
          igx_low = (option%gint_range(1) - option%origin(1)) / option%del(1) + 1 
          igx_up  = (option%gint_range(2) - option%origin(1)) / option%del(1) + 1

          igx_low_ref = (option%gint_refrange(1) - option%origin(1)) / option%del(1) + 1 
          igx_up_ref  = (option%gint_refrange(2) - option%origin(1)) / option%del(1) + 1 

          ! Note: rectangular integration scheme is used
          Gint     = 0.0d0
          Gint_err = 0.0d0
          do igx = igx_low, igx_up - 1
            Gint = Gint + g0%data(igx, 1, 1) * option%del(1)
          end do

          if (igx_low_ref /= igx_up_ref) then
            Gint_ref = 0.0d0
            do igx = igx_low_ref, igx_up_ref - 1
              Gint_ref = Gint_ref + g0%data(igx, 1, 1) * option%del(1)
            end do

            Gint = (Gint / Gint_ref) * option%del(1)
          end if

        end if

      end if ! use_bootstrap

      ! normalization for radial coordinate
      !
      if (.not. option%use_bootstrap) then
        if (option%use_radial) then
          dr     = option%del(option%radcol)
          fourpi = 4.0d0 * PI
       
          if (option%radcol == 1) then
            do igr = 1, option%ng3(1)
              r                   = option%origin(1) + (igr - 1) * dr
              fact                = fourpi * r * r
              if (r < EPS) then
                fact              = fourpi / 3.0d0 * dr * dr 
              end if
              g0%data(igr, :, :)  = g0%data(igr, :, :) / fact 
            end do 
          else if (option%radcol == 2) then
            do igr = 1, option%ng3(2)
              r                   = option%origin(2) + (igr - 1) * dr
              fact                = fourpi * r * r
              if (r < EPS) then
                fact              = fourpi / 3.0d0 * dr * dr 
              end if
              g0%data(:, igr, :)  = g0%data(:, igr, :) / fact 
            end do 
          else if (option%radcol == 3) then
            do igr = 1, option%ng3(3)
              r                   = option%origin(3) + (igr - 1) * dr
              fact                = fourpi * r * r
              if (r < EPS) then
                fact              = fourpi / 3.0d0 * dr * dr 
              end if
              g0%data(:, :, igr)  = g0%data(:, :, igr) / fact 
            end do 
          end if
        end if
      end if

      if (option%calc_statepop) then

        if (option%use_bootstrap) then
          fe_state(:) = fe_ave(:)
        else
          sref = state_count(1)
          do istate = 1, nstate
            sval = state_count(istate)
            fe_state(istate) = - kT * log(sval/sref)
          end do
        end if
      end if

      write(iw,'(">> Finished")')

      ! spline interpolation
      !
      if (option%use_spline) then

        write(iw,*)
        write(iw,'("Analyze> Spline interpolation")')

        if (option%ndim == 1) then
          call akima1d(g0%ng3, g1%ng3, g0%box, 0.0d0, 1.0d20, g0%data, g1%data) 
        else if (option%ndim == 2) then
          call akima2d(g0%ng3, g1%ng3, g0%box, 0.0d0, 1.0d20, g0%data, g1%data) 
        else if (option%ndim == 3) then
          call akima3d(g0%ng3, g1%ng3, g0%box, 0.0d0, 1.0d20, g0%data, g1%data) 
        end if

        write(iw,'(">> Finished")')
      end if

      ! Erase the values outside the region of interest
      !
      if (option%use_cutrange) then
        call cut_gval(option, g0, EPS)
        call cut_gval(option, g1, EPS)
      end if

      ! Normalize G 
      !
      !if (.not. option%space3D) then

      if (option%ndim < 3) then
        call get_gmax(option, g0, gmax0)
        g0%data = g0%data / gmax0

        if (option%use_spline) then
          call get_gmax(option, g1, gmax1)
          g1%data = g1%data / gmax1
        end if
      else
        if (.not. option%space3D) then
          call get_gmax(option, g0, gmax0)
          g0%data = g0%data / gmax0

          if (option%use_spline) then
            call get_gmax(option, g1, gmax1)
            g1%data = g1%data / gmax1
          end if
        end if
      end if

      ! PMF calculation 
      !
      call generate_pmf(kT, g0, pmf0)
      if (option%use_spline) then
        call generate_pmf(kT, g1, pmf1)
      end if

      ! Voronoi tessellation
      !
      if (option%use_voronoi) then
        write(iw,*)
        write(iw,'("Analyze> Start Voronoi analysis")')
        call setup_voronoi(option%ndim,                       &
                           option%vr_ncell,                   &
                           option%vr_normvec,                 &
                           voronoi,                           &
                           funcin     = pmf1,                 &
                           separation = option%vr_separation)

        ! determine cell positions
        !
        if (option%use_cellsearch) then
          call search_cellpos(VrCellTypeMIN, pmf1, voronoi)
        else

          voronoi%cellpos(1:3, 1:voronoi%ncell) &
            = option%vr_cellpos(1:3, 1:voronoi%ncell)

          write(iw,*)
          write(iw,'("Input minimum (x, y, z)")')
          do ic = 1, voronoi%ncell
            write(iw,'(i3,2x,3f20.10)') &
              ic, (voronoi%cellpos(jdim, ic), jdim = 1, option%ndim)
          end do
          write(iw,*)
        end if

        ! determine boundary
        !
        if (option%ndim) then
          call voronoi_boundary(voronoi)
        end if

        write(iw,'(">> Finished")')
      end if

      if (option%use_voronoi) then
        write(fhead_out,'(a,".cellpos")') trim(output%fhead)
        open(UnitOUT, file=trim(fhead_out))
        do ic = 1, voronoi%ncell
          write(UnitOUT,'(3f20.10)') &
            (voronoi%cellpos(jdim, ic), jdim = 1, option%ndim) 
        end do
        close(UnitOUT)

        if (option%ndim) then
          write(fhead_out,'(a,".boundary")') trim(output%fhead)
          open(UnitOUT, file=trim(fhead_out))
          do ib = 1, voronoi%nbound_points
            write(UnitOUT,'(2f20.10)') &
              voronoi%boundary(1, ib), voronoi%boundary(2, ib)
          end do
          close(UnitOUT)
        end if
      end if

      if (option%ndim == 1) then

        if (option%calc_gintegral .and. option%ndim == 1) then
          write(fhead_out,'(a,"_Gint.dat")') trim(output%fhead)
          call open_file(fhead_out, iunit)
          write(iunit,'(2f15.7)') Gint, Gint_err
        end if

        write(fhead_out,'(a,"_g_original")') trim(output%fhead)
        call generate_1dmap_gnuplot(fhead_out, g0)
        write(fhead_out,'(a,"_original")') trim(output%fhead)
        call generate_1dmap_gnuplot(fhead_out, pmf0, pmf_err)
        if (option%use_spline) then
          write(fhead_out,'(a,"_g_spline")') trim(output%fhead)
          call generate_1dmap_gnuplot(fhead_out, g1)
          write(fhead_out,'(a,"_spline")') trim(output%fhead)
          call generate_1dmap_gnuplot(fhead_out, pmf1)
        end if
      else if (option%ndim == 2) then
        write(fhead_out,'(a,"_g_original")') trim(output%fhead)
        call generate_2dmap_matplotlib(fhead_out, pmf0) 
        call generate_2dmap_gnuplot(fhead_out, pmf0)
        write(fhead_out,'(a,"_original")') trim(output%fhead)
        call generate_2dmap_matplotlib(fhead_out, g0) 
        call generate_2dmap_gnuplot(fhead_out, g0)

        if (option%use_spline) then
          write(fhead_out,'(a,"_g_spline")') trim(output%fhead)
          call generate_2dmap_matplotlib(fhead_out, g1) 
          call generate_2dmap_gnuplot(fhead_out, g1)
          write(fhead_out,'(a,"_spline")') trim(output%fhead)
          call generate_2dmap_matplotlib(fhead_out, pmf1) 
          call generate_2dmap_gnuplot(fhead_out, pmf1)
        end if
      else if (option%ndim == 3) then
        write(fhead_out,'(a,"_g_original")') trim(output%fhead)
        call write_dx(fhead_out, g0) 
        write(fhead_out,'(a,"_original")')   trim(output%fhead)
        call write_dx(fhead_out, pmf0)
        if (option%use_spline) then
          write(fhead_out,'(a,"_g_spline")') trim(output%fhead)
          call write_dx(fhead_out, g1) 
          write(fhead_out,'(a,"_spline")')   trim(output%fhead)
          call write_dx(fhead_out, pmf1) 
        end if 
      end if

      if (option%calc_statepop) then
        write(fhead_out,'(a,".festate")') trim(output%fhead)
        call open_file(fhead_out, iunit)
        do istate = 1, nstate
          write(iunit,'(2f15.7)') fe_state(istate), fe_err(istate)
        end do
        close(iunit)
      end if

      deallocate(cv)

    end subroutine analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine analyze_bootstrap(option, output, bootopt, cv, nfile, g0, pmf_err, fe_ave, fe_err, Gint, Gint_err)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),  intent(in)    :: option
      type(s_output),  intent(in)    :: output
      type(s_bootopt), intent(inout) :: bootopt
      type(s_cv),      intent(in)    :: cv(:)
      integer,         intent(in)    :: nfile
      type(s_func3d),  intent(inout) :: g0
      real(8),         intent(inout) :: pmf_err(:)
      real(8),         intent(inout) :: fe_ave(:)
      real(8),         intent(inout) :: fe_err(:)
      real(8),         intent(inout) :: Gint
      real(8),         intent(inout) :: Gint_err

      ! local variables
      !
      type(s_booteach) :: beach
      type(s_bootave)  :: bave

      real(8), allocatable :: state_count(:)
      real(8), allocatable :: fe_state(:, :)

      integer :: ig, istep, isample, itrial, ifile, icv, istate, nstep, active
      integer :: igx, igy, igz, igx_low, igx_up, igx_low_ref, igx_up_ref
      integer :: ntrial, ngrid, nsample, nstate
      integer :: gind(MaxDim)
      real(8) :: val, weight_sum, dv, vol, rcmin, rcmax, kT, fave, ferr
      real(8) :: gmax0, gmaxe, gerr, v1, v2, Gint_ref
      logical :: is_inside


      ! Prepare parameters
      !
      ntrial  = bootopt%ntrial
      nsample = bootopt%nsample

      beach%ntrial  = ntrial
      beach%nsample = nsample

      dv = 1.0d0
      do icv = 1, option%ndim
        dv = dv * option%del(icv)
      end do

      kT = Boltz * option%temperature

      if (option%space3D) then
        vol = option%box_system(1) * option%box_system(2) * option%box_system(3)
      end if
  
      ! Allocate memory
      !
      if (bootopt%boottype == BootTypeTRAJ) then
        allocate(beach%rand(nsample, ntrial))
      else if (bootopt%boottype == BootTypeSNAP) then
        allocate(beach%rand(nsample, 1))
      end if

      allocate(beach%g0(ntrial))
      allocate(bave%yerr(option%ng3(1)))

      if (option%calc_statepop) then
        nstate = option%nstate
        allocate(state_count(nstate))
        allocate(fe_state(nstate, ntrial))
      end if

      do itrial = 1, ntrial
        call setup_func3d(option%ng3, option%del, option%origin, beach%g0(itrial))
      end do
      call setup_func3d(option%ng3, option%del, option%origin, bave%g0)

      ! Generate random seed
      !
      call get_seed(bootopt%iseed)
      call initialize_random(bootopt%iseed)

      ! Generate random numbers
      !
      if (bootopt%boottype == BootTypeTRAJ) then
        do itrial = 1, ntrial
          call get_random_integer(nsample, 1, nfile, bootopt%duplicate, &
                                  beach%rand(1, itrial))
        end do
      end if 

      write(iw,*)
      write(iw,'("Analyze_Bootstrap> Start MC-Bootstrap")')

      if (bootopt%boottype == BootTypeTRAJ) then

        call bootstrap_traj(option,              &
                            bootopt,             &
                            cv,                  &
                            nfile,               &
                            g0,                  &
                            beach,               &
                            bave,                &
                            state_count,         &
                            fe_state)

      else if (bootopt%boottype == BootTypeSNAP) then

        call bootstrap_snap(option,              &
                            bootopt,             &
                            cv,                  &
                            nfile,               &
                            g0,                  &
                            beach,               &
                            bave,                &
                            state_count,         &
                            fe_state)

      end if


      if (option%calc_statepop) then
        do istate = 1, nstate
          fave           = sum(fe_state(istate, 1:ntrial)) / dble(ntrial)
          fe_ave(istate) = fave 

          ferr = 0.0d0
          do itrial = 1, ntrial
            ferr = ferr + (fe_state(istate, itrial) - fave)**2
          end do
          ferr = sqrt(ferr / dble(ntrial - 1))

          fe_err(istate) = ferr
        end do
      end if

      ! Average
      !
      do igz = 1, option%ng3(3)
        do igy = 1, option%ng3(2)
          do igx = 1, option%ng3(1)
            do itrial = 1, ntrial
              bave%g0%data(igx, igy, igz) = bave%g0%data(igx, igy, igz) &
                + beach%g0(itrial)%data(igx, igy, igz) 
            end do
          end do
        end do
      end do
      g0%data(:, :, :) = bave%g0%data(:, :, :)

      if (option%use_cutrange) then

        call cut_gval(option, g0, EPS)

        do itrial = 1, ntrial
          call cut_gval(option, beach%g0(itrial), EPS)
        end do

      end if

      if (option%ndim < 3) then
        call get_gmax(option, g0, gmax0)
        g0%data = g0%data / gmax0

        do itrial = 1, ntrial
          call get_gmax(option, beach%g0(itrial), gmaxe)
          beach%g0(itrial)%data = beach%g0(itrial)%data / gmaxe
        end do
      end if

      if (option%ndim == 1) then
        do igx = 1, option%ng3(1)
          gerr = 0.0d0
          do itrial = 1, ntrial
            v1   = g0%data(igx, 1, 1)
            v2   = beach%g0(itrial)%data(igx, 1, 1)

            if (v1 <= EPS) then
              v1 = - kT * log(EPS)
            else
              v1 = - kT * log(v1)
            end if

            if (v2 <= EPS) then
              v2 = - kT * log(EPS)
            else
              v2 = - kT * log(v2)
            end if

            gerr = gerr + (v2 - v1)**2 
          end do
          gerr   = sqrt(gerr / dble(ntrial - 1)) 
          bave%yerr(igx) = gerr
        end do

        pmf_err(:) = bave%yerr(:)

      end if

      ! gintegral
      !
      if (option%calc_gintegral .and. option%ndim == 1) then
        igx_low     = (option%gint_range(1) - option%origin(1)) / option%del(1) + 1 
        igx_up      = (option%gint_range(2) - option%origin(1)) / option%del(1) + 1
        igx_low_ref = (option%gint_refrange(1) - option%origin(1)) / option%del(1) + 1 
        igx_up_ref  = (option%gint_refrange(2) - option%origin(1)) / option%del(1) + 1

        ! Note: rectangular integration scheme is used
        
        bave%Gint = 0.0d0
        do igx = igx_low, igx_up - 1
          bave%Gint = bave%Gint + g0%data(igx, 1, 1) * option%del(1)
        end do

        if (igx_low_ref /= igx_up_ref) then
          Gint_ref = 0.0d0
          do igx = igx_low_ref, igx_up_ref - 1
            Gint_ref = Gint_ref + g0%data(igx, 1, 1) * option%del(1)
          end do

          bave%Gint = (bave%Gint / Gint_ref) * option%del(1) 
        end if

        gerr = 0.0d0
        do itrial = 1, ntrial
          Gint = 0.0d0
          do igx = igx_low, igx_up - 1
            Gint = Gint + beach%g0(itrial)%data(igx, 1, 1) * option%del(1)
          end do

          if (igx_low_ref /= igx_up_ref) then
            Gint_ref = 0.0d0
            do igx = igx_low_ref, igx_up_ref - 1
              Gint_ref = Gint_ref + beach%g0(itrial)%data(igx, 1, 1) * option%del(1)
            end do
          
            Gint = (Gint / Gint_ref) * option%del(1) 
          end if

          gerr = gerr + (Gint - bave%Gint)**2

        end do

        gerr          = sqrt(gerr / dble(ntrial - 1))
        bave%Gint_err = gerr

        Gint     = bave%Gint
        Gint_err = bave%Gint_err

      end if


      ! For using main routine 
      !
      if (option%ndim < 3) then
        g0%data = g0%data * gmax0
      end if

      ! Deallocate memory
      !
      deallocate(beach%rand)
      if (option%calc_statepop) then
        deallocate(state_count)
      end if

    end subroutine analyze_bootstrap
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine bootstrap_traj(option,                                  &
                              bootopt,                                 &
                              cv,                                      &
                              nfile,                                   &
                              g0,                                      &
                              beach,                                   &
                              bave,                                    &
                              state_count,                             &
                              fe_state)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_bootopt),  intent(inout) :: bootopt
      type(s_cv),       intent(in)    :: cv(:)
      integer,          intent(in)    :: nfile
      type(s_func3d),   intent(inout) :: g0
      type(s_booteach), intent(inout) :: beach
      type(s_bootave),  intent(inout) :: bave
      real(8),          intent(inout) :: state_count(:)
      real(8),          intent(inout) :: fe_state(:, :)
     
      integer :: icv, itrial, isample, istep, ifile, istate, igr, nstep, active
      integer :: gind(MaxDim)
      real(8) :: weight_sum, rcmin, rcmax, val, dv 
      logical :: is_inside

      integer :: ntrial, nsample, nstate
      real(8) :: kT, vol, r, dr, fourpi, fact


      ! Prepare parameters
      !
      nstate  = option%nstate
      ntrial  = bootopt%ntrial
      nsample = bootopt%nsample

      beach%ntrial  = ntrial
      beach%nsample = nsample

      dv = 1.0d0
      do icv = 1, option%ndim
        dv = dv * option%del(icv)
      end do

      kT     = Boltz * option%temperature
      dr     = option%del(option%radcol)
      fourpi = 4.0d0 * PI

      if (option%space3D) then
        vol = option%box_system(1) * option%box_system(2) * option%box_system(3)
      end if

      !$omp parallel private(itrial, isample, ifile, istep, istate, icv,    &
      !$omp                  weight_sum, nstep, is_inside,                  &
      !$omp                  active, rcmin, rcmax, val, gind, state_count), &
      !$omp          default(shared)
      !$omp do
      
      do itrial = 1, ntrial

        if (mod(itrial, 10) == 0) then
          write(iw,'("trial : ", i0)') itrial
        end if

        weight_sum = 0.0d0

        if (option%calc_statepop) then
          state_count = 0
        end if

        do isample = 1, nsample
          ifile = beach%rand(isample, itrial)
          nstep = cv(ifile)%nstep
          do istep = option%nstart, cv(ifile)%nstep
            is_inside = .true.
            do icv = 1, option%ndim
              val       = cv(ifile)%data(option%xyzcol(icv), istep)
              gind(icv) = (val - option%origin(icv)) / option%del(icv) + 1

              if (gind(icv) < 1 .or. gind(icv) > option%ng3(icv)) then
                is_inside = .false.
              end if
            end do

            if (option%calc_statepop) then
              do istate = 1, nstate
                active = 0
                do icv = 1, option%ndim
                  rcmin = option%staterange(1, icv, istate)
                  rcmax = option%staterange(2, icv, istate)
                  val   = cv(ifile)%data(option%xyzcol(icv),istep)

                  if (val >= rcmin .and. val < rcmax) then
                    active = active + 1 
                  end if

                end do

                if (active == option%ndim) then
                  state_count(istate) = state_count(istate) + cv(ifile)%weight(istep) 
                end if

              end do
            end if

            if (option%ndim == 1) then
              gind(2:3) = 1
            else if (option%ndim == 2) then
              gind(3)   = 1
            end if

            if (is_inside) then
              beach%g0(itrial)%data(gind(1), gind(2), gind(3)) = beach%g0(itrial)%data(gind(1), gind(2), gind(3)) &
                + cv(ifile)%weight(istep)
              weight_sum = weight_sum + cv(ifile)%weight(istep)
            end if
          end do

        end do ! isample

        if (option%ndim < 3) then
          beach%g0(itrial)%data = beach%g0(itrial)%data / (weight_sum * dv) 
        else
          if (option%space3D) then
            beach%g0(itrial)%data = vol * g0%data / (option%nstep_system * dv)
          else
            beach%g0(itrial)%data = beach%g0(itrial)%data / (weight_sum * dv)
          end if

        end if

        if (option%calc_statepop) then
          do istate = 1, nstate
            fe_state(istate, itrial) = - kT * log(state_count(istate)/state_count(1))
          end do
        end if

      end do   ! itrial

      !$omp end do
      !$omp end parallel

      ! Consider Jacobian for radial coordinate
      !
      if (option%use_radial) then
        do itrial = 1, ntrial
        
          dr     = option%del(option%radcol)
          fourpi = 4.0d0 * PI
          
          if (option%radcol == 1) then
            do igr = 1, option%ng3(1)
              r                   = option%origin(1) + (igr - 1) * dr
              fact                = fourpi * r * r
              if (r < EPS) then
                fact              = fourpi / 3.0d0 * dr * dr 
              end if
              beach%g0(itrial)%data(igr, :, :) &
                = beach%g0(itrial)%data(igr, :, :) / fact 
            end do 
          else if (option%radcol == 2) then
            do igr = 1, option%ng3(2)
              r                   = option%origin(2) + (igr - 1) * dr
              fact                = fourpi * r * r
              if (r < EPS) then
                fact              = fourpi / 3.0d0 * dr * dr 
              end if
              beach%g0(itrial)%data(:, igr, :)  &
                = beach%g0(itrial)%data(:, igr, :) / fact 
            end do 
          else if (option%radcol == 3) then
            do igr = 1, option%ng3(3)
              r                   = option%origin(3) + (igr - 1) * dr
              fact                = fourpi * r * r
              if (r < EPS) then
                fact              = fourpi / 3.0d0 * dr * dr 
              end if
              beach%g0(itrial)%data(:, :, igr) &
                = beach%g0(itrial)%data(:, :, igr) / fact 
            end do 
          end if
        
        end do
      end if

    end subroutine bootstrap_traj
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine bootstrap_snap(option,                                  &
                              bootopt,                                 &
                              cv,                                      &
                              nfile,                                   &
                              g0,                                      &
                              beach,                                   &
                              bave,                                    &
                              state_count,                             &
                              fe_state)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_bootopt),  intent(inout) :: bootopt
      type(s_cv),       intent(in)    :: cv(:)
      integer,          intent(in)    :: nfile
      type(s_func3d),   intent(inout) :: g0
      type(s_booteach), intent(inout) :: beach
      type(s_bootave),  intent(inout) :: bave
      real(8),          intent(inout) :: state_count(:)
      real(8),          intent(inout) :: fe_state(:, :)
     
      integer :: icv, itrial, irand, isample, istep, ifile, istate, igr, nstep, active
      integer :: gind(MaxDim)
      real(8) :: weight_sum, rcmin, rcmax, val, dv 
      logical :: is_inside

      integer :: nstep_tot, ntrial, nsample, nstate
      real(8) :: kT, vol, r, dr, fourpi, fact


      ! Prepare parameters
      !
      nstep   = cv(1)%nstep
      nstate  = option%nstate
      ntrial  = bootopt%ntrial
      nsample = bootopt%nsample

      beach%ntrial  = ntrial
      beach%nsample = nsample

      dv = 1.0d0
      do icv = 1, option%ndim
        dv = dv * option%del(icv)
      end do

      nstep_tot = 0
      do ifile = 1, nfile
        nstep_tot = nstep_tot + cv(ifile)%nstep
      end do

      kT     = Boltz * option%temperature
      dr     = option%del(option%radcol)
      fourpi = 4.0d0 * PI

      if (option%space3D) then
        vol = option%box_system(1) * option%box_system(2) * option%box_system(3)
      end if

      ! Check
      !
      do ifile = 1, nfile
        if (cv(ifile)%nstep /= cv(1)%nstep) then
          write(iw,'("Bootstrap_Snap> Error.")')
          write(iw,'("Time steps in file (", i0,") is different from that in in file 1.")')
          stop
        end if
      end do

      ! $omp parallel private(itrial, irand, isample, ifile, istep, istate, icv,    &
      ! $omp                  weight_sum, is_inside,                                &
      ! $omp                  active, rcmin, rcmax, val, gind, state_count),        &
      ! $omp          default(shared)
      ! $omp do
      
      do itrial = 1, ntrial

        if (mod(itrial, 10) == 0) then
          write(iw,'("trial : ", i0)') itrial
        end if

        ! Generate random numbers
        !
        call get_random_integer(nsample,           &
                                1,                 &
                                nstep_tot,         &
                                bootopt%duplicate, &
                                beach%rand(1, 1))

        ! Initialize
        !
        weight_sum = 0.0d0

        if (option%calc_statepop) then
          state_count = 0
        end if

        ! Analyze
        !
        do isample = 1, nsample
          irand = beach%rand(isample, 1)
          ifile = (irand - 1) / nstep + 1 
          istep = irand - nstep * (ifile - 1)

          is_inside = .true.
          do icv = 1, option%ndim
            val       = cv(ifile)%data(option%xyzcol(icv), istep)
            gind(icv) = (val - option%origin(icv)) / option%del(icv) + 1

            if (gind(icv) < 1 .or. gind(icv) > option%ng3(icv)) then
              is_inside = .false.
            end if
          end do

          if (option%calc_statepop) then
            do istate = 1, nstate
              active = 0
              do icv = 1, option%ndim
                rcmin = option%staterange(1, icv, istate)
                rcmax = option%staterange(2, icv, istate)
                val   = cv(ifile)%data(option%xyzcol(icv),istep)

                if (val >= rcmin .and. val < rcmax) then
                  active = active + 1 
                end if

              end do

              if (active == option%ndim) then
                state_count(istate) = state_count(istate) + cv(ifile)%weight(istep) 
              end if

            end do
          end if

          if (option%ndim == 1) then
            gind(2:3) = 1
          else if (option%ndim == 2) then
            gind(3)   = 1
          end if

          if (is_inside) then
            beach%g0(itrial)%data(gind(1), gind(2), gind(3))     &
              = beach%g0(itrial)%data(gind(1), gind(2), gind(3)) &
                + cv(ifile)%weight(istep)
            weight_sum = weight_sum + cv(ifile)%weight(istep)
          end if

        end do ! isample

        if (option%ndim < 3) then
          beach%g0(itrial)%data = beach%g0(itrial)%data / (weight_sum * dv) 
        else
          if (option%space3D) then
            beach%g0(itrial)%data = vol * g0%data / (option%nstep_system * dv)
          else
            beach%g0(itrial)%data = beach%g0(itrial)%data / (weight_sum * dv)
          end if

        end if

        if (option%calc_statepop) then
          do istate = 1, nstate
            fe_state(istate, itrial) = - kT * log(state_count(istate)/state_count(1))
          end do
        end if

        ! Consider Jacobian for radial coordinate
        !
        if (option%use_radial) then

          dr     = option%del(option%radcol)
          fourpi = 4.0d0 * PI
         
          if (option%radcol == 1) then
            do igr = 1, option%ng3(1)
              r                   = option%origin(1) + (igr - 1) * dr
              fact                = fourpi * r * r
              if (r < EPS) then
                fact              = fourpi / 3.0d0 * dr * dr 
              end if
              beach%g0(itrial)%data(igr, :, :) &
                = beach%g0(itrial)%data(igr, :, :) / fact 
            end do 
          else if (option%radcol == 2) then
            do igr = 1, option%ng3(2)
              r                   = option%origin(2) + (igr - 1) * dr
              fact                = fourpi * r * r
              if (r < EPS) then
                fact              = fourpi / 3.0d0 * dr * dr 
              end if
              beach%g0(itrial)%data(:, igr, :)  &
                = beach%g0(itrial)%data(:, igr, :) / fact 
            end do 
          else if (option%radcol == 3) then
            do igr = 1, option%ng3(3)
              r                   = option%origin(3) + (igr - 1) * dr
              fact                = fourpi * r * r
              if (r < EPS) then
                fact              = fourpi / 3.0d0 * dr * dr 
              end if
              beach%g0(itrial)%data(:, :, igr) &
                = beach%g0(itrial)%data(:, :, igr) / fact 
            end do 
          end if

        end if

      end do   ! itrial

      ! $omp end do
      ! $omp end parallel


    end subroutine bootstrap_snap
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine generate_pmf(kT, g, pmf)
!-----------------------------------------------------------------------
      implicit none

      real(8),        intent(in)    :: kT
      type(s_func3d), intent(in)    :: g
      type(s_func3d), intent(inout) :: pmf

      integer :: igx, igy, igz
      real(8) :: val, gmax


      !gmax    = maxval(g%data)
      do igz = 1, g%ng3(3)
        do igy = 1, g%ng3(2)
          do igx = 1, g%ng3(1)
            !val = g%data(igx, igy, igz) / gmax
            val = g%data(igx, igy, igz)
            if (val <= EPS) then
              pmf%data(igx, igy, igz) = -kT*log(EPS) 
            else
              pmf%data(igx, igy, igz) = -kT*log(val)
            end if
          end do
        end do
      end do

    end subroutine generate_pmf
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine get_gmax(option, g, gmax)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)  :: option
      type(s_func3d), intent(in)  :: g
      real(8),        intent(out) :: gmax 

      integer :: igx, igy, igz, icv
      real(8) :: gval, crd(3)
      logical :: is_inside

      
      gmax = -1.0d10

      if (option%use_refrange) then
        if (option%ndim == 1) then
          do igx = 1, g%ng3(1)
            crd(1) = g%origin(1) + (igx - 1) * g%del(1)
            if (crd(1) >= option%refrange(1, 1) .and. &
                crd(1) <= option%refrange(2, 1)) then
              gval = g%data(igx, 1, 1)
              if (gmax < gval) then
                gmax = gval 
              end if
            end if
          end do
        else if (option%ndim == 2) then
          do igy = 1, g%ng3(2)
          do igx = 1, g%ng3(1)
       
            crd(1) = g%origin(1) + (igx - 1) * g%del(1)
            crd(2) = g%origin(2) + (igy - 1) * g%del(2)
       
            is_inside = .true.
            do icv = 1, 2
              if (crd(icv) < option%refrange(1, icv) .and. &
                  crd(icv) > option%refrange(2, icv)) then
                is_inside = .false. 
              end if
            end do
       
            if (is_inside) then
              gval = g%data(igx, igy, 1)
              if (gmax < gval) then
                gmax = gval 
              end if
            end if
          end do
          end do
       
        end if

      else

        gmax = maxval(g%data)

      end if

    end subroutine get_gmax 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine cut_gval(option, g, lowval)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)     :: option
      type(s_func3d), intent(inout)  :: g
      real(8),        intent(in)     :: lowval 

      integer :: igx, igy, igz, icv
      real(8) :: crd(3)
      logical :: is_inside

      
      if (option%ndim == 1) then
        do igx = 1, g%ng3(1)
          crd(1) = g%origin(1) + (igx - 1) * g%del(1)
          is_inside = .true.
          if (crd(1) < option%cutrange(1, 1) .or. &
              crd(1) > option%cutrange(2, 1)) then
            is_inside = .false.
          end if

          if (.not. is_inside) then
            g%data(igx, 1, 1) = lowval
          end if

        end do
      else if (option%ndim == 2) then
        do igy = 1, g%ng3(2)
        do igx = 1, g%ng3(1)
      
          crd(1) = g%origin(1) + (igx - 1) * g%del(1)
          crd(2) = g%origin(2) + (igy - 1) * g%del(2)
      
          is_inside = .true.
          do icv = 1, 2
            if (crd(icv) < option%cutrange(1, icv) .or. &
                crd(icv) > option%cutrange(2, icv)) then
              is_inside = .false. 
            end if
          end do
      
          if (.not. is_inside) then
            g%data(igx, igy, 1) = lowval
          end if

        end do
        end do
      
      end if

    end subroutine cut_gval 
!-----------------------------------------------------------------------


end module mod_analyze
!=======================================================================
