!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_ctrl
  use mod_cv
  use mod_bootstrap
  use mod_random

  ! constants
  !

  ! structures
  !
  type :: s_state
    integer :: nstep
    integer :: step_final
    integer :: ndim
    integer :: nBevents   = 0
    integer :: nUBevents  = 0
    logical :: is_reacted = .false.
    logical :: is_reacted_final = .false.
    integer, allocatable :: data(:)
    real(8), allocatable :: data_ts(:) 
    real(8), allocatable :: avedata(:) 
    real(8), allocatable :: hist(:)
    real(8), allocatable :: norm(:)
  end type s_state

  type :: s_booteach
    integer              :: ntrial
    integer              :: nsample
    integer              :: nt_range
    integer, allocatable :: rand(:, :)
    real(8), allocatable :: taud(:), kins(:)
    real(8), allocatable :: hist(:, :)
  end type s_booteach

  type :: s_bootave
    integer              :: nt_range
    real(8)              :: taud, taud_stdev, taud_err
    real(8)              :: kins, kins_stdev, kins_err
    real(8), allocatable :: hist(:), hist_err(:), hist_stdev(:)
    real(8), allocatable :: cumm(:)
  end type s_bootave

  ! subroutines
  !
  public  :: analyze 
  private :: analyze_bootstrap
  private :: calc_kine_bootstrap
  private :: calc_kineave_bootstrap

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option, cvinfo, bootopt)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),   intent(in)    :: input
      type(s_output),  intent(in)    :: output
      type(s_option),  intent(in)    :: option
      type(s_cvinfo),  intent(in)    :: cvinfo
      type(s_bootopt), intent(inout) :: bootopt

      type(s_cv),    allocatable :: cv(:)
      type(s_state), allocatable :: state(:)


      integer :: ifile, istep, itr, ntr, ncut
      integer :: nreact, nreact_final, nstay, nstay_bound
      integer :: nbound, nubound
      real(8) :: rpinit, kr, t, dth
      real(8) :: tau_r, tau_d

      integer, allocatable :: nbcheck(:, :), nubcheck(:, :)
      real(8), allocatable :: rp(:), norm_tot(:) 
      character(len=MaxChar) :: frp


!$    write(iw,*)
!$    write(iw,'("Analyze> OpenMP parallelization is activated")')
      allocate(cv(cvinfo%nfile), state(cvinfo%nfile))

      ! read cv files
      !
      write(iw,*)
      write(iw,'("Analyze> Read CV file")')
      do ifile = 1, cvinfo%nfile
        call read_cv(cvinfo%fcv(ifile),         &
                     option%ndim,               &
                     cv(ifile),                 &
                     nsta = option%nsta,        &
                     nstep_read = option%nstep)
      end do

      ! get state
      !
      write(iw,*)
      write(iw,'("Analyze> Get State")')
!$omp parallel private(ifile), default(shared)
!$omp do
      do ifile = 1, cvinfo%nfile
          call get_state(option%ndim, option%state_def, &
                       cv(ifile), state(ifile))
      end do
!$omp end do
!$omp end parallel

      ! check transient time scale if check_tscale = .true. 
      !
      if (option%check_tscale) then
        write(iw,*)
        write(iw,'("Analyze> Check transient time-scale")')
        allocate(nbcheck(0:option%nntr, cvinfo%nfile))
        allocate(nubcheck(0:option%nntr, cvinfo%nfile))

        nbcheck  = 0
        nubcheck = 0
!$omp parallel private(ifile), default(shared)
!$omp do
        do ifile = 1, cvinfo%nfile
          call check_transient(option%use_moving,    &
                  option%ntr_sta,                    &
                  option%ntr_interval,               &
                  option%nntr,                       &
                  state(ifile),                      &
                  nbcheck(0, ifile),                 &
                  nubcheck(0, ifile))
        end do
!$omp end do
!$omp end parallel

        ! check the convergence of # of events
        !
        do itr = 0, option%nntr
          ntr = option%ntr_sta + option%ntr_interval * itr
          write(iw,'("ntr = ",i10," nBevents = ",i10,&
                    " nUBevents = ",i10)') &
                    ntr, sum(nbcheck(itr,  1:cvinfo%nfile)), &
                    sum(nubcheck(itr, 1:cvinfo%nfile))  
        end do

        ! calculate rate coefficients 
        !
        nstay = 0
        do ifile = 1, cvinfo%nfile 
          do istep = 1, state(ifile)%nstep
            if (state(ifile)%data(istep) == REACTIVE) &
              nstay = nstay + 1
          end do
        end do

        nstay_bound = 0
        do ifile = 1, cvinfo%nfile 
          do istep = 1, state(ifile)%nstep
            if (state(ifile)%data(istep) == BOUND) &
              nstay_bound = nstay_bound + 1
          end do
        end do


        write(frp, '(a,".kcheck")') trim(output%fhead)
        open(UnitOUT, file=trim(frp))
          write(UnitOUT,'("# Delta T  k_insertion  k_exclude")')
          do itr = 0, option%nntr
            ntr = option%ntr_sta + itr * option%ntr_interval
            t    = ntr * option%dt
            write(UnitOUT,'(3(e20.10))') &
              t, sum(nbcheck(itr,  1:cvinfo%nfile))  / (dble(nstay)       * option%dt), &
                 sum(nubcheck(itr,  1:cvinfo%nfile)) / (dble(nstay_bound) * option%dt)
          end do
        close(UnitOUT)

        deallocate(nbcheck, nubcheck)
        stop
      end if 

      ! search reaction
      !
      write(iw,*)
      write(iw,'("Analyze> Search Reaction")')

      if (option%use_moving) then
!$omp parallel private(ifile), default(shared)
!$omp do
        do ifile = 1, cvinfo%nfile
          call search_reaction_moving(&
                 option%tubtype,      &
                 option%nt_transient, &
                 state(ifile))
        end do
!$omp end do
!$omp end parallel
        do ifile = 1, cvinfo%nfile
          write(iw,'("Traj. ",i10," : # of Binds   = ",i10,   &
                                  " , # of UnBinds = ",i10)') &
            ifile, state(ifile)%nBevents, state(ifile)%nUBevents 
        end do

      else
!$omp parallel private(ifile), default(shared)
!$omp do
        do ifile = 1, cvinfo%nfile
          call search_reaction(ifile, option%nt_transient, state(ifile))
        end do
!$omp end do
!$omp end parallel
      end if

      ! FOR CHECK
      do ifile = 1, cvinfo%nfile
        write(iw,'(" Step_Final of file ", i0 , " = ", i0)') &
          ifile, state(ifile)%step_final
      end do

      ! Analyze RP
      !
      if (option%calc_rp) then
        write(iw,*)
        write(iw,'("Analyze> Calculate RP")')
        ncut = nint(option%t_rpcut / option%dt)
        if (ncut /= 0) then
          write(iw,'("Analyze> Detect non-zero t_rpcut")')
          write(iw,'("First ", f15.7," [time] trajectories are used for RP calculation.")') &
            option%t_rpcut
        end if

        if (ncut /= 0) then
          do ifile = 1, cvinfo%nfile
            write(iw,'("Progress ",i0," / ",i0)') ifile, cvinfo%nfile
            call get_rphist(option%nt_range,     &
                            option%use_window,   &
                            state(ifile),        &
                            option%n_return_cut, &
                            ncut = ncut)
          end do
        else
          do ifile = 1, cvinfo%nfile
            write(iw,'("Progress ",i0," / ",i0)') ifile, cvinfo%nfile
            call get_rphist(option%nt_range,     &
                            option%use_window,   &
                            state(ifile),        &
                            option%n_return_cut)
          end do
        end if
      end if
      
      if (option%use_bootstrap) then 
        call analyze_bootstrap(option, output, bootopt, state, cvinfo%nfile)
      else

        if (option%calc_rp) then
          allocate(rp(0:option%nt_range), norm_tot(0:option%nt_range))
          rp       = 0.0d0
          norm_tot = 0.0d0
          do istep = 0, option%nt_range
            do ifile = 1, cvinfo%nfile
              rp(istep)       = rp(istep)       + state(ifile)%hist(istep)
              norm_tot(istep) = norm_tot(istep) + state(ifile)%norm(istep)
            end do
          end do
       
          if (option%use_window) then
            rpinit = rp(0)
            rp(:)  = rp(:) / rpinit
          else
            rp(:)  = rp(:) / norm_tot(:) 
          end if
        
          tau_d = 0.0d0
          dth   = 0.5d0 * option%dt
          do istep = 0, option%nt_range - 1
            tau_d = tau_d + dth * (rp(istep) + rp(istep + 1)) 
          end do
        end if

        ! Analyze Reaction
        !
        if (option%use_moving) then 
          nreact = 0
          do ifile = 1, cvinfo%nfile
            nreact = nreact + state(ifile)%nBevents
          end do
        else
          nreact = 0
          do ifile = 1, cvinfo%nfile
            if (state(ifile)%is_reacted) then
              nreact = nreact + 1 
            end if
          end do
        end if

        !write(iw,'("Binding Trajectory ID: ")')
        !do ifile = 1, cvinfo%nfile
        !  if (state(ifile)%is_reacted) then
        !    write(iw,'(i5)') ifile
        !  end if
        !end do

        !nstay = 0
        !do ifile = 1, cvinfo%nfile
        !  if (state(ifile)%is_reacted) then 
        !    do istep = 1, state(ifile)%nstep
        !      if (state(ifile)%data(istep) == REACTIVE) &
        !        nstay = nstay + 1
        !    end do
        !  else
        !    do istep = 1, state(ifile)%nstep
        !      if (state(ifile)%data(istep) /= UNBOUND) &
        !        nstay = nstay + 1
        !    end do
        !  end if
        !end do

        nstay = 0
        do ifile = 1, cvinfo%nfile 
          do istep = 1, state(ifile)%step_final
            if (state(ifile)%data(istep) == REACTIVE) &
              nstay = nstay + 1
          end do
        end do

        ! NEW
        kr    = dble(nreact) / (dble(nstay) * option%dt)
        tau_r = 1.0d0 / kr 

        write(iw,*)
        write(iw,'("----- Results -----")')
        !write(iw,'("# of reaction events : ", i0)') nreact
        write(iw,'("# of reaction events : ", i0)') nreact
        write(iw,'("k_r    (time^-1)     : ", es15.7)') kr
        write(iw,'("tau_r  (time)        : ", es15.7)') tau_r

        if (option%calc_rp) then
          write(iw,'("tau_d  (time)        : ", es15.7)') tau_d

          write(frp, '(a,".rp")') trim(output%fhead)
          open(10,file=trim(frp))
          do istep = 0, option%nt_range
            t = istep * option%dt
            write(10,'(2f20.10)') t, rp(istep)
          end do
          close(10)

          deallocate(rp, norm_tot)
        end if

      end if

      !deallocate(rp)

    end subroutine analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_state(ndim, state_def, cv, state)
!-----------------------------------------------------------------------
      implicit none

      integer,                intent(in)  :: ndim
      real(8),                intent(in)  :: state_def(2, ndim, Nstate)
      type(s_cv),             intent(in)  :: cv
      type(s_state),          intent(out) :: state

      integer :: istep, istate, icv, ia
      integer :: nstep
      logical :: is_assigned
      real(8) :: wrk(ndim) 


      nstep       = cv%nstep

      ! for global
      allocate(state%data(nstep), state%data_ts(nstep), &
               state%avedata(nstep))

      do istep = 1, cv%nstep
        wrk(:) = cv%data(:, istep) 
        is_assigned = .false.
        do istate = 1, nstate
          if (is_assigned) then
            exit
          else
            ia = 0
            do icv = 1, ndim 
              if (wrk(icv) >= state_def(1, icv, istate) &
                .and. wrk(icv) < state_def(2, icv, istate)) then
                ia = ia + 1
              end if
            end do

            if (ia == ndim) then 
              is_assigned       = .true.
              state%data(istep) = StateInfo(istate)
            end if

          end if
        end do

        if (.not. is_assigned) then
          state%data(istep) = OTHERS 
        end if
      end do

      do istep = 1, nstep
        state%data_ts(istep) = StateConv(state%data(istep))
      end do
      
      state%nstep = nstep
      state%ndim  = ndim

    end subroutine get_state 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine check_transient(use_moving, ntr_sta, ntr_interval, nntr, state, &
                               nbcheck, nubcheck)
!-----------------------------------------------------------------------
      implicit none

      logical,                intent(in)    :: use_moving 
      integer,                intent(in)    :: ntr_sta, ntr_interval, nntr
      type(s_state),          intent(inout) :: state
      integer,                intent(out)   :: nbcheck(0:nntr)
      integer,                intent(out)   :: nubcheck(0:nntr)

      integer :: istep, jstep, icv, iseq, current, next
      integer :: itr, ntr, nstep, nstep_sparse, ns
      integer :: nBevents, nUBevents
      real(8) :: r, rcurr, rnext

      logical :: is_reacted, is_changed


      nstep = state%nstep

      if (use_moving) then
        do itr = 0, nntr
          ntr = ntr_sta + itr * ntr_interval

          state%avedata = 0.0d0
          ns = 2*ntr + 1
          jstep = 0 
          do istep = ntr + 1, nstep - ntr, ntr
            jstep = jstep + 1
            r = sum(state%data_ts(istep-ntr/2:istep+ntr/2)) / dble(ns) 
            state%avedata(jstep) = r 
            !write(iw,'("step ",i10," : ",f15.7)') istep, r 
          end do
          !nstep_sparse = jstep
          nstep_sparse = jstep + 1
          
          r = sum(state%data_ts(nstep-ntr/2:nstep)) / dble(ns)
          state%avedata(nstep_sparse) = r

          nBevents  = 0
          nUBevents = 0

          rcurr = state%data_ts(1) 
          istep = 0 
          is_changed = .true. 
          do while(is_changed)
            is_changed = .false.
            do jstep = istep + 1, nstep_sparse
              rnext = state%avedata(jstep)
              if (rcurr <  0.0d0 .and. rnext  >=  0.0d0) then
                nBevents = nBevents + 1
                rcurr    = rnext
                istep    = jstep
                is_changed = .true.
                exit
              else if (rcurr >= 0.0d0 .and. rnext < 0.0d0) then
                nUBevents = nUBevents + 1
                rcurr     = rnext
                istep     = jstep
                is_changed = .true.
                exit
              end if
            end do 

          end do 
          !do istep = 2 * (ntr + 1), nstep - ntr, ntr
          !do istep = 2, nstep_sparse
          !  next  = state%avedata(istep)

          !  if (current < -0.25d0) then
          !    if (next >= +0.25d0) then
          !      nBevents = nBevents + 1
          !    end if 
          !  else if (current >= 0.25d0) then
          !    if (next < -0.25d0) then
          !      nUBevents = nUBevents + 1
          !    end if 
          !  end if
          !
          !  current = next
          !end do

          nbcheck(itr)  = nBevents
          nubcheck(itr) = nUBevents

        end do

      else
        do itr = 0, nntr
          ntr = ntr_sta + itr * ntr_interval
       
          nBevents  = 0
          nUBevents = 0 
          !current = state%data(1)
          current = UNBOUND 
          do istep = 2, nstep - ntr + 1
            next = state%data(istep)
            if (current /= BOUND) then
         
              if (next == BOUND) then
                is_changed = .true.
                do jstep = istep + 1, istep + ntr - 1
                  if (state%data(jstep) /= BOUND) then
                    is_changed = .false.
                  end if
                end do
         
                if (is_changed) then
                  nBevents = nBevents + 1
                  current  = next 
                end if
              end if
         
            else if (current == BOUND) then
         
              if (next /= BOUND) then
                is_changed = .true.
                do jstep = istep + 1, istep + ntr - 1
                  if (state%data(jstep) == BOUND) then
                    is_changed = .false.
                  end if
                end do
         
                if (is_changed) then
                  nUBevents = nUBevents + 1
                  current   = next
                end if
              end if
         
            end if
          end do
       
          nbcheck(itr)  = nBevents
          nubcheck(itr) = nUBevents
       
        end do
      end if

    end subroutine check_transient 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine search_reaction_moving(tubtype, nt_transient, state)
!-----------------------------------------------------------------------
      implicit none

      integer,                intent(in)    :: tubtype
      integer,                intent(in)    :: nt_transient
      type(s_state),          intent(inout) :: state

      integer :: istep, jstep, icv, iseq, current, next
      integer :: itr, ntr, nstep, nstep_sparse, ns
      integer :: nBevents, nUBevents
      real(8) :: r, rcurr, rnext
      real(8) :: state_final 

      logical :: is_reacted, is_changed


      nstep = state%nstep

      state%avedata = 0.0d0
      ntr = nt_transient 
      ns  = 2*ntr + 1
      jstep = 0 
      do istep = ntr + 1, nstep - ntr, ntr
        jstep = jstep + 1
        r = sum(state%data_ts(istep-ntr/2:istep+ntr/2)) / dble(ns)
        state%avedata(jstep) = r
      end do
      !nstep_sparse = jstep
      nstep_sparse = jstep + 1
      
      r = sum(state%data_ts(nstep-ntr/2:nstep)) / dble(ns)
      state%avedata(nstep_sparse) = r

      nBevents  = 0
      nUBevents = 0

      rcurr = state%data_ts(1) 
      istep = 0 
      is_changed = .true. 
      do while(is_changed)
        is_changed = .false.
        do jstep = istep + 1, nstep_sparse
          rnext = state%avedata(jstep)
          if (rcurr <  0.0d0 .and. rnext  >=  0.0d0) then
              nBevents = nBevents + 1
              rcurr    = 1.0d0 
              istep    = jstep
              is_changed = .true.
              exit
          else if (rcurr >= 0.0d0 .and. rnext < 0.0d0) then
              nUBevents = nUBevents + 1
              rcurr     = -1.0d0 
              istep     = jstep
              is_changed = .true.
              exit
          end if
        end do 

      end do

      state%nBevents  = nBevents
      state%nUBevents = nUBevents

      state%is_reacted = .false.
      if (nBevents > 0) then
        state%is_reacted = .true.
      end if

      state_final = state%avedata(nstep_sparse)

      state%step_final  = nstep
      if (tubtype == TubTypeCUTFINAL) then
        if (state_final < 0.0d0) then ! if unbound state at final step
          do istep = 1, nstep_sparse
            jstep = nstep_sparse - istep
            rcurr = state%avedata(jstep)
 
            if (rcurr > 0.0d0) then
              !state%step_final = ntr * (istep - 1) + ntr / 2  
              state%step_final = ntr * jstep + 1 
              exit 
            end if

          end do 
        end if
      else
        state%step_final = nstep
      end if 

    end subroutine search_reaction_moving 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine search_reaction(ifile, nt_transient, state)
!-----------------------------------------------------------------------
      implicit none

      integer,                intent(in)    :: ifile
      integer,                intent(in)    :: nt_transient
      type(s_state),          intent(inout) :: state

      integer :: istep, jstep, icv, iseq, current, next
      integer :: nstep
      integer :: nBevents, nUBevents

      logical :: is_reacted, is_changed


      nstep = state%nstep

      is_reacted = .false.
      do istep = 1, nstep - nt_transient + 1

        if (is_reacted) exit

        if (state%data(istep) == BOUND) then
          is_reacted = .true.
          do jstep = istep + 1, istep + nt_transient - 1
            !if (state%data(jstep) /= BOUND) then
            !  write(iw,'(i0)') state%data(jstep)
            !end if
            if (state%data(jstep) /= BOUND) then
              is_reacted = .false.
            end if
          end do

        end if
      end do

      nBevents  = 0
      nUBevents = 0 
      current = state%data(1)   
      do istep = 2, nstep - nt_transient + 1
        next = state%data(istep)
        if (current /= BOUND) then

          if (next == BOUND) then
            is_changed = .true.
            do jstep = istep + 1, istep + nt_transient - 1
              if (state%data(jstep) /= BOUND) then
                is_changed = .false.
              end if
            end do

            if (is_changed) then
              nBevents = nBevents + 1
              current  = next 
            end if
          end if

        else ! current == BOUND

          if (next /= BOUND) then
            is_changed = .true.
            do jstep = istep + 1, istep + nt_transient - 1
              if (state%data(jstep) == BOUND) then
                is_changed = .false.
              end if
            end do

            if (is_changed) then
              nUBevents = nUBevents + 1
              current   = next
            end if
          end if

        end if
      end do 


      if (state%data(nstep) == BOUND) then
        state%is_reacted_final = .true. 
      end if

      state%is_reacted = is_reacted
      state%nBevents   = nBevents
      state%nUBevents  = nUBevents
      state%step_final = nstep 


    end subroutine search_reaction 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_rphist(nt_range, use_window, state, n_return_cut, ncut)
!-----------------------------------------------------------------------
      implicit none

      integer,                intent(in)    :: nt_range
      logical,                intent(in)    :: use_window
      type(s_state),          intent(inout) :: state
      integer,                intent(in)    :: n_return_cut
      integer, optional,      intent(in)    :: ncut 

      integer :: istep, jstep, kstep, it
      integer :: nstep, nt, nrc
      logical :: chk, count_return

      real(8), allocatable :: hist(:), norm(:)


      nstep = state%nstep
      if (present(ncut)) &
        nstep = ncut

      nrc = n_return_cut

      ! for global
      if (.not. allocated(state%hist)) & 
        allocate(state%hist(0:nt_range))
      if (.not. allocated(state%norm)) &
        allocate(state%norm(0:nt_range))
      ! for local
      if (.not. allocated(hist)) & 
        allocate(hist(0:nt_range))
      if (.not. allocated(norm)) & 
        allocate(norm(0:nt_range))

      state%hist = 0.0d0
      state%norm = 0.0d0

      if (state%is_reacted) then
        deallocate(hist)
        return
      end if

      hist = 0.0d0
      norm = 0.0d0
      if (nrc == 0) then
        if (use_window) then
          do istep = 1, nstep - nt_range - 1
            if (state%data(istep) == REACTIVE) then
            !if (state%data(istep) /= UNBOUND) then
              it = 0
              do jstep = istep, istep + nt_range
               
                if (state%data(jstep) == REACTIVE) then
                !if (state%data(jstep) /= UNBOUND) then
                  hist(it) = hist(it) + 1.0d0
                end if
                it = it + 1
              end do
            end if
          end do
        else if (.not. use_window) then
          do istep = 1, nstep - 1
            if (state%data(istep) == REACTIVE) then
              it = 0
       
              nt = istep + nt_range
              if (nt > nstep) then
                nt = nstep
              end if
       
              do jstep = istep, nt
                norm(it) = norm(it) + 1.0d0 
                if (state%data(jstep) == REACTIVE) then
                  hist(it) = hist(it) + 1.0d0
                end if
                it       = it + 1
              end do
            end if
          end do
        end if
      else if (nrc /= 0) then
        if (use_window) then
          do istep = 1, nstep - nt_range - 1
            if (state%data(istep) == REACTIVE) then
              it           = 0
              count_return = .true.
              do jstep = istep, istep + nt_range
                if (count_return) then 
                  if (state%data(jstep) == REACTIVE) then
                 
                    if (it <= nrc) then
                      hist(it) = hist(it) + 1.0d0
                    else
                      chk = .false.
                      do kstep = 1, nrc - 1
                        if (state%data(jstep - kstep) == REACTIVE) then
                          chk = .true.
                          exit
                        end if
                      end do
                 
                      if (chk) then
                        hist(it) = hist(it) + 1.0d0
                      else
                        count_return = .false. 
                      end if
                 
                    end if
                  end if

                end if
                it = it + 1
              end do

            end if
          end do
        else if (.not. use_window) then
          do istep = 1, nstep - 1
            if (state%data(istep) == REACTIVE) then
            !if (state%data(istep) /= UNBOUND) then
              it = 0
       
              nt = istep + nt_range
              if (nt > nstep) then
                nt = nstep
              end if
      
              count_return = .true.
              do jstep = istep, nt
                norm(it) = norm(it) + 1.0d0
                if (count_return) then
                  if (state%data(jstep) == REACTIVE) then
                    if (it <= nrc) then
                      hist(it) = hist(it) + 1.0d0
                    else
                      chk = .false.
                      do kstep = 1, nrc - 1
                        if (state%data(jstep - kstep) == REACTIVE) then
                          chk          = .true.
                          exit
                        end if
                      end do
                 
                      if (chk) then
                        hist(it) = hist(it) + 1.0d0
                      else
                        count_return = .false.
                      end if
                    end if
                 
                  end if
                end if
                it       = it + 1
              end do
            end if
          end do
        end if


      end if

      state%hist = hist
      state%norm = norm
      deallocate(hist, norm)

    end subroutine get_rphist 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine analyze_bootstrap(option, output, bootopt, state, nfile)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),         intent(in)    :: option
      type(s_output),         intent(in)    :: output 
      type(s_bootopt),        intent(inout) :: bootopt
      type(s_state),          intent(in)    :: state(:)
      integer,                intent(in)    :: nfile

      ! structures
      !

      ! local variables
      !
      type(s_booteach)       :: beach
      type(s_bootave)        :: bave

      integer, allocatable   :: nstay(:)

      integer                :: i, itrial, isample, istep
      integer                :: j, k
      integer                :: ntrial, nsample
      integer                :: nt_range
      character(len=MaxChar) :: fout


      ! Prepare parameters
      !
      ntrial   = bootopt%ntrial
      nsample  = bootopt%nsample
      nt_range = option%nt_range

      beach%ntrial   = ntrial
      beach%nsample  = nsample
      beach%nt_range = nt_range 

    
      ! Allocate memory
      !
      allocate(nstay(1:nfile))
      allocate(beach%rand(nsample, ntrial))
      allocate(beach%taud(ntrial))
      allocate(beach%kins(ntrial))
      allocate(beach%hist(0:nt_range, ntrial))
      allocate(bave%hist(0:nt_range))
      allocate(bave%hist_err(0:nt_range), bave%hist_stdev(0:nt_range))
      allocate(bave%cumm(0:nt_range))

      ! Generate random seed
      !
      call get_seed(bootopt%iseed)
      call initialize_random(bootopt%iseed)

      ! Generate random numbers
      !
      do itrial = 1, ntrial
        call get_random_integer(nsample, 1, nfile, &
                                bootopt%duplicate, beach%rand(1, itrial))
      end do

      ! Calculate time spent in reactive region for each traj.
      !
      nstay = 0 
      !$omp parallel private(i, j), default(shared)
      !$omp do
      do i = 1, nfile
        do j = 1, state(i)%step_final
          if (state(i)%data(j) == REACTIVE) &
            nstay(i) = nstay(i) + 1
        end do
      end do
      !$omp end do
      !$omp end parallel

      ! Calculate quantities for each trial
      !
      call calc_kine_bootstrap(option, state, nstay, beach)

      ! Calculate average
      !
      call calc_kineave_bootstrap(option, beach, bave)

      ! Output
      !
      write(fout,'(a,".rp")') trim(output%fhead)
      open(UnitOUT, file=trim(fout))
      do i = 0, option%nt_range
        write(UnitOUT, '(4f20.10)') &
          i*option%dt, bave%hist(i), bave%hist_err(i), bave%cumm(i)
      end do
      close(UnitOUT)

      write(fout,'(a,".prop")') trim(output%fhead)
      open(UnitOUT, file=trim(fout))
      !write(UnitOUT,'(a20,2x,a15,2x,a15,2x,a15)') &
      !        "Property", "Value", "SD", "SE" 
      write(UnitOUT,'(a20,2x,a15,2x,a15)') &
              "Property", "Value", "SE"
      write(UnitOUT,'(a20,2x,2(es15.7,2x))') &
              "tau_ret  (time)    :", bave%taud, bave%taud_stdev !, bave%taud_err
      write(UnitOUT,'(a20,2x,2(es15.7,2x))') &
              "k_ins    (time^1)  :", bave%kins, bave%kins_stdev !, bave%kins_err

      close(UnitOUT)

      deallocate(nstay)
      deallocate(beach%rand)
      deallocate(beach%taud)
      deallocate(beach%kins)
      deallocate(beach%hist)
      deallocate(bave%hist)
      deallocate(bave%hist_err)
      deallocate(bave%hist_stdev)
      deallocate(bave%cumm)

    end subroutine analyze_bootstrap
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine calc_kine_bootstrap(option, state, nstay, beach)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_state),    intent(in)    :: state(:)
      integer,          intent(in)    :: nstay(:)
      type(s_booteach), intent(inout) :: beach

      integer :: istep, itraj, itrial, isample
      integer :: ntrial, nsample, nt_range
      integer :: nev, nuev

      real(8), allocatable :: norm_tot(:)  


      ntrial   = beach%ntrial
      nsample  = beach%nsample
      nt_range = beach%nt_range

      allocate(norm_tot(0:nt_range))

      beach%hist = 0.0d0
      beach%taud = 0.0d0
      beach%kins = 0.0d0
      norm_tot   = 0.0d0

      write(iw,'("Calc_Kine_Bootstrap> Start MC-Bootstrap")')

      !$omp parallel private(itrial, isample, itraj, istep, &
      !$omp                  nev, nuev, norm_tot), &
      !$omp          shared(nstay)
      !$omp do 
      do itrial = 1, ntrial

        if (mod(itrial, 10) == 0) then
          write(iw,'("trial : ",i0)') itrial
        end if

        nev      = 0
        nuev     = 0
        norm_tot = 0.0d0
        if (option%use_moving) then
          do isample = 1, nsample
            itraj = beach%rand(isample, itrial)
            do istep = 0, nt_range
              beach%hist(istep, itrial) = beach%hist(istep, itrial) + state(itraj)%hist(istep)
              norm_tot(istep)           = norm_tot(istep) + state(itraj)%norm(istep)
            end do
            nev  = nev  + state(itraj)%nBevents
            nuev = nuev + nstay(itraj) 
          end do
        else
          do isample = 1, nsample
            itraj = beach%rand(isample, itrial)
            if (state(itraj)%is_reacted) then
              nev = nev + 1
            end if
            nuev = nuev + nstay(itraj)
          end do
        end if

        if (option%use_window) then
          beach%hist(:, itrial) = beach%hist(:, itrial) &
                                / beach%hist(0, itrial)
        else
          beach%hist(:, itrial) = beach%hist(:, itrial) &
                                / norm_tot(:) 

        end if

        beach%kins(itrial) = dble(nev) / (dble(nuev) * option%dt)
      end do
      !$omp end do
      !$omp end parallel

      beach%taud = 0.0d0
      do itrial = 1, ntrial
        do istep = 0, nt_range - 1
          beach%taud(itrial) = beach%taud(itrial) &
                             + 0.5d0 * option%dt &
                               * (beach%hist(istep, itrial) &
                                + beach%hist(istep + 1, itrial))
        end do
      end do

      deallocate(norm_tot)

    end subroutine calc_kine_bootstrap
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine calc_kineave_bootstrap(option, beach, bave)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_booteach), intent(in)    :: beach
      type(s_bootave),  intent(inout) :: bave

      integer :: itrial, istep
      integer :: ntrial, nt_range
      real(8) :: val, dev, sterr 


      ntrial        = beach%ntrial
      nt_range      = beach%nt_range
      bave%nt_range = nt_range

      ! average
      !
      bave%hist = 0.0d0
      bave%taud = 0.0d0
      bave%kins = 0.0d0

      do itrial = 1, ntrial
        bave%hist(:) = bave%hist(:) + beach%hist(:, itrial)
        bave%taud    = bave%taud    + beach%taud(itrial)
        bave%kins    = bave%kins    + beach%kins(itrial)
      end do

      bave%hist(:) = bave%hist(:) / ntrial
      bave%taud    = bave%taud    / ntrial
      bave%kins    = bave%kins    / ntrial

      bave%cumm    = 0.0d0

      bave%cumm(1) = 0.5d0 * option%dt &
                   * (bave%hist(1) + bave%hist(2))
      do istep = 2, nt_range - 1
        bave%cumm(istep) = bave%cumm(istep - 1)         &
          + 0.5d0 * option%dt                           &
            * (bave%hist(istep) + bave%hist(istep + 1)) 
      end do
      bave%cumm(nt_range) = bave%cumm(nt_range - 1)     &
        + option%dt * bave%hist(nt_range)

      ! error
      !
      ! - taud
      dev = 0.0d0
      do itrial = 1, ntrial
        val = beach%taud(itrial) - bave%taud
        dev = dev + val * val
      end do
      sterr    = dev / ((ntrial-1)*ntrial)
      sterr    = sqrt(sterr)
      bave%taud_err   = sterr
      bave%taud_stdev = sqrt(dev / (ntrial - 1))

      ! - kins 
      dev = 0.0d0
      do itrial = 1, ntrial
        val = beach%kins(itrial) - bave%kins
        dev = dev + val * val
      end do
      sterr    = dev / ((ntrial-1)*ntrial)
      sterr    = sqrt(sterr)
      bave%kins_err   = sterr
      bave%kins_stdev = sqrt(dev / (ntrial - 1)) 

      ! - hist 
      do istep = 1, nt_range
        dev = 0.0d0 
        do itrial = 1, ntrial
          val = beach%hist(istep, itrial) - bave%hist(istep) 
          dev = dev + val * val
        end do
        sterr = dev / (ntrial - 1)
        sterr = sqrt(sterr)
        bave%hist_err(istep)   = sterr
        bave%hist_stdev(istep) = sterr 
        !bave%hist_stdev(istep) = sqrt(dev / (ntrial - 1))
      end do

    end subroutine calc_kineave_bootstrap
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
