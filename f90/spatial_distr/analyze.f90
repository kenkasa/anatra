!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_xyz
  use mod_akima
  use mod_grid3d
  use mod_output
  use mod_ctrl
  use mod_traj
  use mod_fitting
  use mod_cv
  use mod_com

  ! constants
  !
  real(8), parameter, private :: EPS = 1.0d-8 

  ! structures
  !
  type s_state
    integer :: nmol
    integer :: nstep 
    integer, allocatable :: data(:, :)
  end type s_state

  ! subroutines
  !
  public :: analyze 

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option, traj, cvinfo)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in)    :: input
      type(s_output), intent(in)    :: output
      type(s_option), intent(in)    :: option
      type(s_traj),   intent(inout) :: traj(3)
      type(s_cvinfo), intent(in)    :: cvinfo 

      type(s_dcd)    :: dcd
      type(xtcfile)  :: xtc
      type(s_netcdf) :: nc 
      type(s_com)    :: com
      type(s_func3d) :: g0,  g1 

      ! I/O
      !
      integer                :: io
      character(len=MaxChar) :: fname 

      ! Local
      !
      integer                :: nfile, ndim
      integer                :: natm, nmol, nstep_tot
      real(8)                :: weight, norm
      real(8)                :: weight_sum_reac, weight_sum
      integer                :: ng3_spl(3), gind(3)
      real(8)                :: del_spl(3)
      integer                :: ngrids_count
      real(8)                :: dv, val, vol, vol_count, boxave(3)
      real(8)                :: d(3), d2, shift(3)
      logical                :: is_inside

      integer                :: trajtype
      logical                :: is_end

      integer                :: itraj, istep, istep_tot
      integer                :: igx, igy, igz, ixyz
      integer                :: iatm, imol

      character(len=MaxChar) :: fxyz 

      type(s_xyz)            :: fit_xyz
      type(s_fit)            :: fitu 
     
      ! Dummy
      !
      integer :: ifile, idim

      ! Arrays
      !
      type(s_cv),    allocatable :: cv
      type(s_state), allocatable :: state

      integer,    allocatable :: state_tr(:, :)


      ! Setup parameters
      !
      dv = option%del(1) * option%del(2) * option%del(3)

      if (option%use_spline) then
        ng3_spl = 1
        del_spl = 1.0d0
        ng3_spl(1:3) = option%ng3(1:3) &
                                 * option%spline_resolution
        del_spl(1:3) = option%del(1:3) &
                                 / option%spline_resolution 
      end if 

      ! Setup 3d-functions
      !
      call setup_func3d(option%ng3, option%del, option%origin, g0)

      if (option%use_spline) then
        call setup_func3d(ng3_spl, del_spl, option%origin, g1) 
      end if

      ! Load xyz (if needed)
      !
      if (option%fit) then
        write(fxyz,'(a,".sd.xyz")') trim(output%fhead)
        call open_file(fxyz, io)
        call read_xyz (io, fit_xyz)
        close(io)

        traj(3)%coord(:, :, 1) = fit_xyz%coord(:, :)
        call setup_fit(traj(2), fitu)
      end if

      ! Get trajectory type (dcd or xtc)
      !
      call get_trajtype(input%ftraj(1), trajtype)

      ! Setup CoM structure
      !
      call get_com(option%mode,         &
                  traj(1),              &
                  com,                  &
                  setup = .true.,       &
                  calc_coord = .false., &
                  myrank = 0)

      nmol = com%nmol

      ! Start main calculation
      !
      write(iw,*)
      write(iw,'("Analyze> Start SDF calculation")')

      boxave          = 0.0d0
      istep_tot       = 0
      weight_sum      = 0.0d0
      weight_sum_reac = 0.0d0

      do itraj = 1, input%ntraj

        ! Open trajectory
        !
        call open_trajfile(input%ftraj(itraj), trajtype, io, dcd, xtc, nc)

        ! Initialize trajectory
        !
        call init_trajfile(trajtype, io, dcd, xtc, nc, natm)

        ! Get restriction info.
        !
        call get_restriction(input, option, nmol, itraj, cv, state, state_tr)

        ! Get weight info. 
        !
        call get_weight(input, option, itraj, cv)

        is_end     = .false.
        istep      = 0
        do while (.not. is_end)

          istep     = istep     + 1
          istep_tot = istep_tot + 1

          ! Read frame
          !
          call read_trajfile_oneframe(trajtype, io, istep, dcd, xtc, nc, is_end)

          if (is_end) exit

          if (mod(istep_tot, 1000) == 0) then
            write(iw,'("Step ",i0)') istep_tot
          end if

          ! Transfer from dcd or xtc to traj
          !
          call send_coord_to_traj(1, trajtype, dcd, xtc, nc, traj(1))

          ! Get CoM info.
          !
          call get_com(option%mode,          &
                       traj(1),              &
                       com,                  &
                       setup = .false.,      &
                       calc_coord = .true.,  &
                       myrank = 1)

          if (option%fit) then
            call send_coord_to_traj(1, trajtype, dcd, xtc, nc, traj(2))

            do iatm = 1, fitu%natm
              fitu%movcoord(1:3, iatm) = traj(2)%coord(1:3, iatm, 1)
              fitu%refcoord(1:3, iatm) = fit_xyz%coord(1:3, iatm)
            end do

            call get_com_refmov(fitu)
            shift(1:3) = - fitu%movcom(1:3)

            do iatm = 1, fitu%natm
              fitu%movcoord(1:3, iatm) = fitu%movcoord(1:3, iatm) - shift(1:3)
            end do

            call get_trrot(fitu)

            nmol = com%nmol
            do imol = 1, nmol

              ! Shift
              !
              com%coord(1:3, imol, 1) = com%coord(1:3, imol, 1) - shift(1:3)

              ! PBC
              !
              do ixyz = 1, 3
                com%coord(ixyz, imol, 1) = com%coord(ixyz, imol, 1) &
                  - nint((com%coord(ixyz, imol, 1)) / traj(1)%box(ixyz, 1)) * traj(1)%box(ixyz, 1) 
              end do 

              ! Fitting
              ! 
              call operate_trrot(fitu, 1, com%coord(1:3, imol, 1))

            end do 

          end if

          ! Set weight value
          !
          if (option%use_weight) then
            weight = cv%weight(istep)
          else
            weight = 1.0d0
          end if

          ! Wrap (if use_pbcwrap = .true.)
          !
          call pbcwrap(option, traj(1), com)

          ! Update histogram
          !
          call update_hist(option,          &
                           com,             &
                           state_tr,        &
                           istep,           &
                           nmol,            &
                           weight,          &
                           g0,              &
                           weight_sum,      &
                           weight_sum_reac)

          boxave(:) = boxave(:) + traj(1)%box(:, 1)
        end do

        istep_tot = istep_tot - 1

        call close_trajfile(trajtype, io, dcd, xtc, nc)

      end do

      ! Determine normalization const.
      !
      !if (option%is_whole_snap) then

      !  if (option%normalize_reacsnap) then
      !    norm = weight_sum_reac
      !  else
      !    norm = weight_sum
      !  end if

      !else

      !norm = option%nstep_whole * nmol

      !end if

      if (option%use_conditional) then 
        norm = weight_sum_reac
      else
        norm = weight_sum
      end if

      boxave(:) = boxave(:) / dble(istep_tot)
      vol       = boxave(1) * boxave(2) * boxave(3)

      write(iw,'("dv = ", f15.7)') dv
      write(iw,'("box ave = ", 3f15.7)') boxave(1), boxave(2), boxave(3)


      write(fname,'(a,".sdbin")') trim(output%fhead)
      call open_file(fname, io, frmt = 'unformatted')
        write(io) norm
        write(io) nmol
        write(io) vol
        write(io) dv
        write(io) option%use_spline
        write(io) option%spline_resolution
        write(io) g0%ng3
        write(io) g0%del
        write(io) g0%origin
        write(io) g0%box
        write(io) g0%data
      close(io)

      ! Normalize G
      !
      g0%data = vol * g0%data / (norm * dv)

      ! Spline 
      !
      if (option%use_spline) then
        call akima3d(g0%ng3, g1%ng3, g0%box, 0.0d0, 1.0d20, g0%data, g1%data)
      end if

      ! Calculate the volume in which g is larger than zero
      !
      ngrids_count = 0
      do igz = 1, option%ng3(3)
        do igy = 1, option%ng3(2)
          do igx = 1, option%ng3(1)
            val = g0%data(igx, igy, igz)  
            if (abs(val) >= option%count_threshold) &
              ngrids_count = ngrids_count + 1
          end do
        end do
      end do

      vol_count = ngrids_count * dv

      write(iw,*)
      write(iw,'(" Volume where g(r) is larger than ", &
                  e15.7, " / A^3 = ", f20.10)')        &
                  option%count_threshold, vol_count
      write(iw,'(" Percentage = ", f20.10)') vol_count / vol * 100.0d0

      ! generate dx file
      !
      write(fname,'(a,"_g_original")') trim(output%fhead)
      call write_dx(fname, g0) 
      if (option%use_spline) then
        write(fname,'(a,"_g_spline")') trim(output%fhead)
        call write_dx(fname, g1) 
      end if 

      !deallocate(molid, mtot, com)

    end subroutine analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine pbcwrap(option, traj, com)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      type(s_traj),   intent(in)    :: traj
      type(s_com),    intent(inout) :: com

      integer :: id, ixyz, nmol
      real(8) :: shift(3)


      shift = 0.0d0
      nmol  = com%nmol

      if (.not. option%use_pbcwrap) &
        return

      if (option%centertype == CenterTypeZERO) then
        shift(1:3) = 0.0d0
      else if (option%centertype == CenterTypeHALF) then
        shift(1:3) = traj%box(1:3, 1) * 0.5d0
      end if

      do id = 1, nmol
        do ixyz = 1, 3
          com%coord(ixyz, id, 1) = com%coord(ixyz, id, 1) &
            - nint((com%coord(ixyz, id, 1) - shift(ixyz)) &
              / traj%box(ixyz, 1)) * traj%box(ixyz, 1) 
        end do
      end do

    end subroutine pbcwrap
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_restriction(input, option, nmol, itraj, cv, state, state_tr)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),        intent(in)    :: input
      type(s_option),       intent(in)    :: option
      integer,              intent(in)    :: nmol 
      integer,              intent(in)    :: itraj
      type(s_cv),           intent(inout) :: cv
      type(s_state),        intent(inout) :: state
      integer, allocatable, intent(inout) :: state_tr(:, :)


      if (.not. option%use_conditional) then
        if (.not. allocated(state_tr)) then
          allocate(state_tr(1, 1))
        end if
        return
      end if

      if (allocated(cv%data)) then

        deallocate(cv%data)
        deallocate(state%data)
        deallocate(state_tr)
        cv%nstep = 0
        cv%ndim  = 0

      end if

      call read_cv(input%fcv(itraj), option%ndim * nmol, cv)
      call get_state(option%ndim, nmol, cv%nstep, option%state_def, cv, state)
      allocate(state_tr(cv%nstep, nmol))
      state_tr = transpose(state%data)

    end subroutine get_restriction
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_weight(input, option, itraj, cv)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in)    :: input 
      type(s_option), intent(in)    :: option
      integer,        intent(in)    :: itraj
      type(s_cv),     intent(inout) :: cv


      if (.not. option%use_weight) &
        return

      if (allocated(cv%weight)) then
        deallocate(cv%weight)
      end if

      call read_cv_weight(input%fweight(itraj), cv)

    end subroutine get_weight
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine update_hist            (option,                         &
                                       com,                            &
                                       state_tr,                       &
                                       istep,                          &
                                       nmol,                           &
                                       weight,                         &
                                       g0,                             &
                                       weight_sum,                     &
                                       weight_sum_reac)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),       intent(in)    :: option
      type(s_com),          intent(in)    :: com
      integer, allocatable, intent(inout) :: state_tr(:, :)
      integer,              intent(in)    :: istep
      integer,              intent(in)    :: nmol
      real(8),              intent(in)    :: weight
      type(s_func3d),       intent(inout) :: g0
      real(8),              intent(inout) :: weight_sum
      real(8),              intent(inout) :: weight_sum_reac


      ! Restricted Sampling
      !
      if (option%use_conditional) then

        call update_hist_restriction(option,          &
                                     com,             &
                                     state_tr,        &
                                     istep,           &
                                     nmol,            &
                                     weight,          &
                                     g0,              &
                                     weight_sum,      &
                                     weight_sum_reac)


      ! Normal Sampling
      !
      else
        call update_hist_normal     (option,          &
                                     com,             &
                                     state_tr,        &
                                     istep,           &
                                     nmol,            &
                                     weight,          &
                                     g0,              &
                                     weight_sum,      &
                                     weight_sum_reac)

      end if

    end subroutine update_hist
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine update_hist_normal     (option,                         &
                                       com,                            &
                                       state_tr,                       &
                                       istep,                          &
                                       nmol,                           &
                                       weight,                         &
                                       g0,                             &
                                       weight_sum,                     &
                                       weight_sum_reac)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),       intent(in)    :: option
      type(s_com),          intent(in)    :: com
      integer, allocatable, intent(inout) :: state_tr(:, :)
      integer,              intent(in)    :: istep
      integer,              intent(in)    :: nmol
      real(8),              intent(in)    :: weight
      type(s_func3d),       intent(inout) :: g0
      real(8),              intent(inout) :: weight_sum
      real(8),              intent(inout) :: weight_sum_reac

      integer :: id, ixyz, gind(1:3)
      real(8) :: val
      logical :: is_inside

      do id = 1, nmol 
          weight_sum = weight_sum + weight
          is_inside = .true.
          do ixyz = 1, 3 
            val        = com%coord(ixyz, id, 1) 
            gind(ixyz) = (val - option%origin(ixyz)) / option%del(ixyz) + 1
            if (gind(ixyz) < 1 .or. gind(ixyz) > option%ng3(ixyz)) then
              is_inside = .false.
            end if
          end do
      
          if (is_inside) then
            g0%data(gind(1), gind(2), gind(3)) = g0%data(gind(1), gind(2), gind(3)) + weight 
          end if

      end do


    end subroutine update_hist_normal
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine update_hist_restriction(option,                         &
                                       com,                            &
                                       state_tr,                       &
                                       istep,                          &
                                       nmol,                           &
                                       weight,                         &
                                       g0,                             &
                                       weight_sum,                     &
                                       weight_sum_reac)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),       intent(in)    :: option
      type(s_com),          intent(in)    :: com
      integer, allocatable, intent(inout) :: state_tr(:, :)
      integer,              intent(in)    :: istep
      integer,              intent(in)    :: nmol
      real(8),              intent(in)    :: weight
      type(s_func3d),       intent(inout) :: g0
      real(8),              intent(inout) :: weight_sum
      real(8),              intent(inout) :: weight_sum_reac

      integer :: id, ixyz, gind(1:3)
      real(8) :: val
      logical :: is_inside


      do id = 1, nmol

        weight_sum = weight_sum + weight

        if (state_tr(istep, id) == REACTIVE) then

          is_inside       = .true.
          weight_sum_reac = weight_sum_reac + weight

          do ixyz = 1, 3 
            val        = com%coord(ixyz, id, 1)
            gind(ixyz) = (val - option%origin(ixyz)) / option%del(ixyz) + 1

            if (gind(ixyz) < 1 .or. gind(ixyz) > option%ng3(ixyz)) then
              is_inside = .false.
            end if
          end do
         
          if (is_inside) then
            g0%data(gind(1), gind(2), gind(3)) &
              = g0%data(gind(1), gind(2), gind(3)) + weight 
          end if
        end if

      end do

    end subroutine update_hist_restriction
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


      gmax    = maxval(g%data)
      do igz = 1, g%ng3(3)
        do igy = 1, g%ng3(2)
          do igx = 1, g%ng3(1)
            val = g%data(igx, igy, igz) / gmax
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
    subroutine get_state(ndim, nmol, nstep, state_def, cv, state)
!-----------------------------------------------------------------------
      implicit none

      integer,                intent(in)   :: ndim
      integer,                intent(in)   :: nmol
      integer,                intent(in)   :: nstep
      real(8),                intent(in)   :: state_def(2, ndim)
      type(s_cv),             intent(in)   :: cv
      type(s_state),          intent(out)  :: state

      integer :: istep, icv, ic, istate, ia, imol, ncv 
      logical :: is_assigned
      real(8) :: val 


      state%nmol  = nmol
      state%nstep = nstep

      ncv         = ndim

      ! for global
      !if (allocated(state%data)) then
      !  deallocate(state%data)
        allocate(state%data(nmol, nstep))
      !end if


      do istep = 1, nstep
        ic = 0
        do imol = 1, nmol
          ia = 0
          do icv = 1, ncv 
            ic = ic + 1 
            val = cv%data(ic, istep)
            if (val >= state_def(1, icv) .and. val < state_def(2, icv)) then
              ia = ia + 1
            end if 
          end do

          if (ia == ncv) then
            is_assigned = .true.
            state%data(imol, istep) = StateInfo(1) 
          else
            is_assigned = .false.
            state%data(imol, istep) = OTHERS 
          end if 
        end do
      end do


    end subroutine get_state 
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
