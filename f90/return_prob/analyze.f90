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
  use mod_analyze_str
  use mod_reac
  use mod_tcf

  ! constants
  !

  ! structures
  !

  ! subroutines
  !
  public  :: analyze 
  private :: analyze_bootstrap
  private :: calc_kine_bootstrap
  private :: calc_kineave_bootstrap

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option, bootopt)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),   intent(in)    :: input
      type(s_output),  intent(in)    :: output
      type(s_option),  intent(in)    :: option
      type(s_bootopt), intent(inout) :: bootopt

      type(s_cv),    allocatable :: cv(:)
      type(s_state), allocatable :: state(:)


      integer :: ifile, istep, itr, ntr, ncut
      integer :: nBevents, nUBevents
      integer :: nreact, nreact_final, nstay, nstay_bound
      integer :: nstep_sparse
      integer :: nbound, nubound
      integer :: nthreads, thread_id, nwork, ifinish, iper
      real(8) :: rpinit, kr, t, dth
      real(8) :: tau_r, tau_d
      real(8) :: progress

      integer, allocatable :: nbcheck(:, :), nubcheck(:, :)
      real(8), allocatable :: rp(:), norm_tot(:) 
      character(len=MaxChar) :: frp


      ! Memory allocation
      !
      allocate(cv(input%ncv), state(input%ncv))

      ! Read cv files
      !
      write(iw,*)
      write(iw,'("Analyze> Read CV file")')

      iper = 1
      do ifile = 1, input%ncv

        ! Check progress
        !
        progress = ifile / dble(input%ncv)
        if (progress >= iper * 0.1d0) then
          write(iw,'("  Progress : ", f6.2, "%")') progress * 100.0d0
          iper = iper + 1
        end if

        ! Read
        !
        call read_cv(input%fcv(ifile),          &
                     option%ndim * option%nmol, &
                     cv(ifile),                 &
                     nstep_read = option%nstep)
        cv(ifile)%ndim = option%ndim
      end do

      ! Get state info
      !
      write(iw,*)
      write(iw,'("Analyze> Get State")')

!$omp parallel private(ifile, nthreads, thread_id, nwork, progress, ifinish), & 
!$omp          default(shared)
      nthreads  = omp_get_num_threads()
      thread_id = omp_get_thread_num()
      nwork     = input%ncv / nthreads + 1
      ifinish   = 0
!$omp do
      do ifile = 1, input%ncv

        call get_state(option, cv(ifile), state(ifile))

        if (thread_id == 0) then
          ifinish = ifinish + 1
          progress = ifinish / dble(nwork) * 100.0d0
          write(iw,'("  Progress : ", f6.2, "%")') progress
        end if

      end do
!$omp end do
!$omp end parallel

      ! Check transient time scale
      !
      if (option%check_timescale) then
        write(iw,*)
        write(iw,'("Analyze> Check transient time-scale")')
        call check_transient(output, option, state)
        call termination("sr_analysis")
      end if

      ! Count # of binding/unbinding events
      ! (stored in state%nBevents and state%nUBevents)
      !
      write(iw,*)
      write(iw,'("Analyze> Analyze reaction events")')
      call get_reaction_info(option, state)

      ! Calculate transition time-correlation functions 
      !
      write(iw,*)
      write(iw,'("Analyze> Get transition TCF")')
      call get_transtcf(option, state)

      write(iw,*)
      write(iw,'("Analyze> Get kinetic properties")')
      call get_kineprop(option, output, bootopt, state)

      !deallocate(rp)

    end subroutine analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_state(option, cv, state)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),         intent(in)  :: option
      type(s_cv),             intent(in)  :: cv
      type(s_state),          intent(out) :: state 

      integer                :: istep, jstep, kstep, istate
      integer                :: imol, icv, imc, ia, sid 
      integer                :: nstate, nmol, ndim, nstep
      real(8)                :: wrk(option%ndim)
      logical                :: end_judge  
      character(len=MaxChar) :: fname

      logical, allocatable   :: is_assigned(:), is_quenched(:)

      nstate = option%nstate
      nmol   = option%nmol
      ndim   = option%ndim
      nstep  = cv%nstep

      ! for global
      allocate(state%quench_step(nmol), state%read_step(nmol))
      allocate(state%data(nstep, nmol), state%data_ts(nstep, nmol), &
               state%avedata(nstep, nmol))
      allocate(state%is_reacted(nmol))

      ! for local
      if (.not. allocated(is_assigned)) &
        allocate(is_assigned(nmol))
      if (.not. allocated(is_quenched)) then 
        allocate(is_quenched(nmol))
      end if

      state%quench_step(:) = nstep
      state%read_step(:)   = nstep
      is_quenched(:)       = .false.

      do istep = 1, cv%nstep
        is_assigned(:) = .false.
        do imol = 1, nmol
          imc = (imol - 1) * ndim
          wrk(1:ndim) = cv%data(imc+1:imc+ndim, istep)
          do istate = 1, nstate
            if (is_assigned(imol)) & 
              exit

            ia  = 0
            do icv = 1, ndim
              if (wrk(icv) >= option%state_def(1, icv, istate) &
                .and. wrk(icv) < option%state_def(2, icv, istate)) then
                ia = ia + 1
              end if

            end do

            if (ia == ndim) then
              is_assigned(imol)       = .true.
              state%data(istep, imol) = istate
            end if

          end do  ! istate

          if (.not. is_assigned(imol)) then
            state%data(istep, imol) = 0
          end if

          ! Check quenching
          !
          if (option%use_quench) then
            if (.not. is_quenched(imol)) then
              ia  = 0
              do icv = 1, ndim
                if (wrk(icv) >= option%state_def(1, icv, nstate + 1) &
                  .and. wrk(icv) < option%state_def(2, icv, nstate + 1)) then
                  ia = ia + 1
                end if
              end do

              if (ia == ndim) then

                if (option%allow_state_jump) then
                  is_quenched(imol)       = .true.
                  state%quench_step(imol) = istep
                  state%read_step(imol)   = istep
                else

                  end_judge = .false.
                  do jstep = 1, istep
                    kstep = istep - jstep
                    sid   = state%data(kstep, imol)

                    if (sid == option%reaczone_id) then
                      is_quenched(imol)       = .true.
                      state%quench_step(imol) = istep
                      state%read_step(imol)   = istep
                      end_judge               = .true. 
                    else

                      if (sid >= 1 .and. sid <= nstate) then
                        end_judge = .true. 
                      end if

                    end if

                    if (end_judge) exit
                  end do

                end if   ! allow_state_jump
              end if     ! ia ? ndim
            end if       ! is_quenched
          end if         ! use_quench

        end do
      end do

      ! Zero padding 
      !
      if (option%use_quench) then
        if (option%use_zeropadding) then
          do imol = 1, nmol
            if (is_quenched(imol)) then
              state%data(state%quench_step(imol):cv%nstep, imol) = 0 
              !state%quench_step(imol) = cv%nstep
              state%read_step(imol) = cv%nstep
            end if
          end do
        end if
      end if
      

      state%data_ts = OFF
      do imol = 1, nmol
        do istep = 1, nstep
          if (state%data(istep, imol) == option%bound_id) then
            state%data_ts(istep, imol) = ON
          end if
        end do
      end do

      state%nstep = nstep
      state%nmol  = nmol
      state%ndim  = ndim

    end subroutine get_state 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_kineprop(option, output, bootopt, state)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),  intent(in)    :: option
      type(s_output),  intent(in)    :: output
      type(s_bootopt), intent(inout) :: bootopt
      type(s_state),   intent(inout) :: state(:)

      ! Local variables
      !
      type(s_booteach) :: beach
      type(s_bootave)  :: bave

      integer :: nstate, nstep, ntrial, nsample, nt_range, nfile, nmol
      integer :: ifile, istep, itrial, imol, is, iunit
      character(len=MaxChar) :: fname

      integer, allocatable :: nstay(:, :)


      ! Prepare parameters 
      !
      nstate   = option%nstate
      nstep    = state(1)%nstep
      ntrial   = bootopt%ntrial
      nsample  = bootopt%nsample
      nt_range = option%nt_range
      nfile    = size(state(:))
      nmol     = option%nmol

      if (.not. option%use_bootstrap) then
        ntrial  = 1
        nsample = nfile 
      end if

      ! Allocate memory
      !
      allocate(nstay(0:nstate, 1:nfile))

      allocate(beach%rand(nsample, ntrial))
      allocate(beach%taud(ntrial), beach%kd(ntrial))
      allocate(beach%tauins(ntrial), beach%kins(ntrial))
      allocate(beach%hist(0:nt_range, 0:nstate, 0:nstate, ntrial))
      allocate(beach%tau_tpm(0:nstate, 0:nstate, ntrial))
      allocate(beach%k_tpm(0:nstate, 0:nstate, ntrial))

      allocate(bave%hist      (0:nt_range, 0:nstate, 0:nstate))
      allocate(bave%hist_err  (0:nt_range, 0:nstate, 0:nstate))
      allocate(bave%hist_stdev(0:nt_range, 0:nstate, 0:nstate))
      allocate(bave%cumm      (0:nt_range, 0:nstate, 0:nstate))
      allocate(bave%tau_tpm   (0:nstate, 0:nstate))
      allocate(bave%tau_tpm_err(0:nstate, 0:nstate))
      allocate(bave%k_tpm     (0:nstate, 0:nstate))
      allocate(bave%k_tpm_err (0:nstate, 0:nstate))

      if (option%calc_2nd) then
        allocate(beach%hist2(0:nt_range, 0:nstate, ntrial))
        allocate(beach%pa3(0:nt_range, ntrial))
        allocate(beach%taud2(0:nstate, ntrial))
        allocate(beach%taupa3(ntrial))

        allocate(bave%hist2(0:nt_range, 0:nstate))
        allocate(bave%pa3(0:nt_range))
        allocate(bave%taud2(0:nstate))
        allocate(bave%taud2_err(0:nstate))
        allocate(bave%cumm2(0:nt_range, 0:nstate))
        allocate(bave%cummpa3(0:nt_range))
      end if

      ! Setup Bootstrap (if use_bootstrap = .true.) 
      !
      if (option%use_bootstrap) then

        write(iw,*)
        write(iw,'("Get_Kineprop> Setup Bootstrap")')

        ! Generate random seed
        !
        call get_seed(bootopt%iseed)
        call initialize_random(bootopt%iseed)

        ! Generate random numbers
        !
        do itrial = 1, ntrial
          call get_random_integer(nsample,              &
                                  1,                    &
                                  nfile,                &
                                  bootopt%duplicate,    &
                                  beach%rand(1, itrial))
        end do

      end if

      ! Calculate time spent in each region for each traj  
      !
      nstay = 0
!$omp parallel private(ifile, istep, imol, is), shared(nstay), &
!$omp          default(shared)
!$omp do
      do ifile = 1, nfile
        do imol = 1, nmol
          !do istep = 1, nstep
          do istep = 1, state(ifile)%read_step(imol)
            is = state(ifile)%data(istep, imol)
            nstay(is, ifile) = nstay(is, ifile) + 1
          end do
        end do
      end do
!$omp end do
!$omp end parallel

      ! Calculate kinetic properties
      !
      write(iw,*)
      write(iw,'("Get_Kineprop> Calculate average values")')
      call get_kinetics_average(option, output, bootopt, state, nstay, nfile, beach, bave)

      call output_kinetics(option, output, bave)

    end subroutine get_kineprop
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_kinetics_average(option, output, bootopt, state, nstay, nfile, beach, bave)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_output),   intent(in)    :: output
      type(s_bootopt),  intent(in)    :: bootopt
      type(s_state),    intent(in)    :: state(:)
      integer,          intent(in)    :: nstay(0:option%nstate, nfile)
      integer,          intent(in)    :: nfile
      type(s_booteach), intent(inout) :: beach
      type(s_bootave),  intent(inout) :: bave

      integer :: nstate, ntrial, nsample, nt_range, nmol
      integer :: itrial, isample, itraj, istep, imol, is, js, rid
      integer :: nthreads, thread_id, nwork, ifinish
      integer :: nev, nuev
      real(8) :: t, p1, p2, f, rval
      real(8) :: dev, val, sterr
      real(8) :: progress
      real(8), allocatable :: norm_tot(:, :)
      real(8), allocatable :: norm_tot_save(:, :)

      integer                :: iunit
      character(len=MaxChar) :: fname



      ! Setup parameters
      !
      nstate   = option%nstate
      nt_range = option%nt_range
      ntrial   = bootopt%ntrial
      nsample  = bootopt%nsample
      rid      = option%reaczone_id
      nmol     = option%nmol

      if (.not. option%use_bootstrap) then
        ntrial  = 1
        nsample = nfile 
      end if

      ! Memory allocation
      !
      allocate(norm_tot(0:nt_range, 0:nstate))
      allocate(norm_tot_save(0:nt_range, 0:nstate))
      norm_tot = 0.0d0


      beach%hist = 0.0d0
      beach%taud = 0.0d0
      beach%kd   = 0.0d0
      beach%kins = 0.0d0

      beach%hist2 = 0.0d0

!$omp parallel private(itrial, isample, itraj, istep, is, js, nev, nuev, norm_tot, t, &
!$omp                  nthreads, thread_id, nwork, progress, ifinish),   &
!$omp          shared(norm_tot_save, nstay), default(shared)
      nthreads  = omp_get_num_threads()
      thread_id = omp_get_thread_num()
      nwork     = ntrial / nthreads + 1
      ifinish   = 0
!$omp do 
      do itrial = 1, ntrial
        
        if (option%use_bootstrap .and. thread_id == 0) then
          ifinish  = ifinish + 1
          progress = ifinish / dble(nwork) * 100.0d0
          write(iw,'("  Progress : ", f6.2, "%")') progress
        end if

        nev      = 0
        nuev     = 0
        norm_tot = 0.0d0

        do isample = 1, nsample

          if (option%use_bootstrap) then
            itraj = beach%rand(isample, itrial)
          else
            itraj = isample
          end if

          ! Pret part
          !
          if (option%calc_pret) then
            do is = 0, nstate
           
              do istep = 0, nt_range
                norm_tot(istep, is) &
                  = norm_tot(istep, is) + state(itraj)%norm(istep, is)
              end do
           
              do js = 0, nstate
                do istep = 0, nt_range
                  beach%hist(istep, js, is, itrial)       &
                    = beach%hist(istep, js, is, itrial)   &
                    + state(itraj)%hist(istep, js, is)
                end do
              end do
           
            end do
           
            ! Integrate TCFs
            !
            call trape_integral(0,                               &
                                nt_range,                        &
                                option%dt_out,                   &
                                beach%hist(:, rid, rid, itrial), &
                                beach%taud(itrial))

            ! 2nd-order transition
            !
           
            if (option%calc_2nd) then
              do is = 0, nstate
                do istep = 0, nt_range
                  beach%hist2(istep, is, itrial)     &
                    = beach%hist2(istep, is, itrial) &
                    + state(itraj)%hist2(istep, is)
                end do
              end do
            end if

          end if

          ! Kins part
          !
          if (option%calc_kins) then
            nev  = nev  + state(itraj)%nBevents
            nuev = nuev + nstay(option%reaczone_id, itraj)
          end if


        end do

        ! Pret part
        !
        if (option%calc_pret) then
          do is = 0, nstate
            if (norm_tot(0, is) /= 0.0d0) then
              do js = 0, nstate
                do istep = 0, nt_range
                  beach%hist(istep, js, is, itrial)      &
                    = beach%hist(istep, js, is, itrial)  &
                      / norm_tot(istep, is)
                end do
              end do
            end if
          end do

          if (.not. option%use_bootstrap) then
            norm_tot_save(:, :) = norm_tot(:, :) 
          end if

          if (option%calc_2nd) then
         
            ! P_A<=i<=A (t)
            !
            do is = 0, nstate
              do istep = 0, nt_range
                beach%hist2(istep, is, itrial) &
                  = beach%hist2(istep, is, itrial) / norm_tot(istep, option%reaczone_id)
              end do
            end do
         
            ! P_A<=A<=A (t)
            !
            beach%pa3(:, itrial) &
              = beach%hist2(:, option%reaczone_id, itrial)
         
            ! P_A<=i<=A (t) / t
            !
            beach%hist2(0, option%reaczone_id, itrial) &
              = beach%hist2(0, option%reaczone_id, itrial) / option%dt_out
         
            do is = 0, nstate
              do istep = 1, nt_range
                t = option%dt_out * istep
                beach%hist2(istep, is, itrial) = beach%hist2(istep, is, itrial) / t
              end do
            end do
         
          end if

        end if

        ! Kins part
        !
        if (option%calc_kins) then
          beach%kins(itrial)   = dble(nev) / (dble(nuev) * option%dt)
          beach%tauins(itrial) = 1.0d0 / beach%kins(itrial)
        end if


      end do
!$omp end do
!$omp end parallel

!
      if (.not. option%use_bootstrap) then
        if (option%out_normfactor) then
          is = option%reaczone_id 
          write(fname,'(a,".unnorm")') trim(output%fhead)
          call open_file(fname, iunit)
          do istep = 0, nt_range
            t = option%dt_out * istep
            write(iunit,'(f15.7)',advance='no') t
            write(iunit,'(f20.7)',advance='no') beach%hist(istep, is, is, 1) * norm_tot_save(istep, is)
            write(iunit,'(f20.7)',advance='no') norm_tot_save(istep, is)
            write(iunit,*)
          end do
          close(iunit)
        end if
      end if
!

      ! Integrate TCFs for each trial
      !
      if (option%calc_pret) then
        do itrial = 1, ntrial
          call trape_integral(0, nt_range, option%dt_out, beach%hist(:, rid, rid, itrial), beach%taud(itrial))
          beach%kd(itrial) = 1.0d0 / beach%taud(itrial)
        end do


        do itrial = 1, ntrial
          do is = 1, option%nstate
            do js = 1, option%nstate
              call trape_integral(0, nt_range, option%dt_out, beach%hist(:, js, is, itrial), beach%tau_tpm(js, is, itrial))
              beach%k_tpm(js, is, itrial) = 1.0d0 / beach%tau_tpm(js, is, itrial)
              
            end do
          end do
        end do
      end if

      if (option%calc_2nd) then

        do itrial = 1, ntrial

          ! P_A<i<A (t) / t
          do is = 1, nstate
            call trape_integral(0, nt_range, option%dt_out, beach%hist2(:, is, itrial), beach%taud2(is, itrial))
          end do

          ! P_A<A<A (t)
          call trape_integral(0, nt_range, option%dt_out, beach%pa3(:, itrial), beach%taupa3(itrial))

        end do

      end if

      ! MC boostrap average 
      !
      if (option%use_bootstrap) then

        bave%hist    = 0.0d0
        bave%taud    = 0.0d0
        bave%kd      = 0.0d0
        bave%tauins  = 0.0d0
        bave%kins    = 0.0d0
        bave%tau_tpm = 0.0d0
        bave%k_tpm   = 0.0d0

        if (option%calc_2nd) then
          bave%hist2     = 0.0d0
          bave%taud2     = 0.0d0
          bave%taud2_err = 0.0d0
          bave%pa3       = 0.0d0
          bave%taupa3    = 0.0d0
        end if

        do itrial = 1, ntrial
          bave%hist(:, :, :) = bave%hist(:, :, :) + beach%hist(:, :, :, itrial)
          bave%taud          = bave%taud          + beach%taud(itrial)
          !bave%tauins        = bave%tauins        + beach%tauins(itrial)
          bave%kins          = bave%kins          + beach%kins(itrial)

          bave%tau_tpm(:, :) = bave%tau_tpm(:, :) + beach%tau_tpm(:, :, itrial)
          !bave%k_tpm(:, :)   = bave%k_tpm(:, :)   + beach%k_tpm  (:, :, itrial)

          if (option%calc_2nd) then
            bave%hist2(:, :) = bave%hist2(:, :)   + beach%hist2(:, :, itrial)
            bave%taud2(:)    = bave%taud2(:)      + beach%taud2(:, itrial)
            bave%pa3(:)      = bave%pa3(:)        + beach%pa3(:, itrial)
            bave%taupa3      = bave%taupa3        + beach%taupa3(itrial)
          end if

        end do

        bave%hist(:, :, :) = bave%hist(:, :, :) / ntrial
        bave%taud          = bave%taud          / ntrial
        bave%kins          = bave%kins          / ntrial

        bave%kd            = 1.0d0 / bave%taud
        bave%tauins        = 1.0d0 / bave%kins

        bave%tau_tpm(:, :) = bave%tau_tpm(:, :) / ntrial
        bave%k_tpm(:, :)   = 1.0d0 / bave%tau_tpm(:, :)

        if (option%calc_2nd) then
          bave%hist2       = bave%hist2  / ntrial
          bave%taud2       = bave%taud2  / ntrial
          bave%pa3         = bave%pa3    / ntrial
          bave%taupa3      = bave%taupa3 / ntrial
        end if

        ! Running integral of averaged TCFs
        !
        ! P_A<A (t)
        if (option%calc_pret) then

          do is = 0, option%nstate
            do js = 0, option%nstate
              call trape_integral(0, nt_range, option%dt_out, bave%hist(:, js, is), rval, cumm = bave%cumm(0, js, is))
            end do
          end do 
        end if

        if (option%calc_2nd) then

          ! P_A<i<A (t) / t
          do is = 0, nstate 
            call trape_integral(0, nt_range, option%dt_out, bave%hist2(:, is), rval, cumm = bave%cumm2(:, is))
          end do

          ! P_A<A<A (t)
          call trape_integral(0, nt_range, option%dt_out, bave%pa3, rval, cumm = bave%cummpa3)

        end if

        ! Error analysis
        !   Note: StdDev of Bootstrap is StdErr

        if (option%calc_pret) then

          ! - taud
          !
          call standard_deviation(ntrial, bave%taud, beach%taud, bave%taud_stdev)
          bave%taud_err = bave%taud_stdev

          ! - kd
          !
          bave%kd_stdev = bave%taud_stdev / bave%taud**2 
          bave%kd_err   = bave%kd_stdev

          ! - tau_tpm
          !
          do is = 1, nstate
            do js = 1, nstate
              call standard_deviation(ntrial, bave%tau_tpm(js, is), &
                beach%tau_tpm(js, is, 1:ntrial), bave%tau_tpm_err(js, is))

              bave%k_tpm_err(js, is) = bave%tau_tpm_err(js, is) / bave%tau_tpm(js, is)**2 
            end do
          end do

          ! - hist
          !
          do is = 0, nstate
            do js = 0, nstate
              do istep = 0, nt_range
                call standard_deviation(ntrial,                         &
                                        bave%hist(istep, js, is),       &
                                        beach%hist(istep, js, is, :),   &
                                        bave%hist_stdev(istep, js, is))
                bave%hist_err = bave%hist_stdev
              end do
            end do
          end do

        end if


        if (option%calc_kins) then

          ! - kins
          !
          call standard_deviation(ntrial, bave%kins, beach%kins, bave%kins_stdev)
          bave%kins_err = bave%kins_stdev

          ! - tauins
          !
          bave%tauins_stdev = bave%kins_stdev / bave%kins**2 
          bave%tauins_err   = bave%tauins_stdev

        end if


        if (option%calc_2nd) then

          ! - taud2
          !
          do is = 1, nstate
            call standard_deviation(ntrial,             &
                                    bave%taud2(is),     &
                                    beach%taud2(is, :), &
                                    bave%taud2_err(is))
          end do

          ! - taupa3
          !
          call standard_deviation(ntrial,         &
                                  bave%taupa3,    &
                                  beach%taupa3,   &
                                  bave%taupa3_err)

        end if

      else
        bave%taud         = beach%taud(1)
        bave%taud_stdev   = 0.0d0

        bave%kd           = beach%kd(1)
        bave%taud_stdev   = 0.0d0

        bave%tauins       = beach%tauins(1)
        bave%tauins_stdev = 0.0d0

        bave%kins         = beach%kins(1)
        bave%kins_stdev   = 0.0d0

        bave%hist         = beach%hist(:, :, :, 1)
        bave%hist_stdev   = 0.0d0

        bave%tau_tpm(:, :) = beach%tau_tpm(:, :, 1)
        bave%tau_tpm_err   = 0.0d0

        bave%k_tpm(:, :)   = beach%k_tpm(:, :, 1)
        bave%k_tpm_err     = 0.0d0 

        if (option%calc_2nd) then
          bave%hist2   = beach%hist2(:, :, 1)
          bave%pa3     = beach%pa3(:, 1)
          bave%taud2   = beach%taud2(:, 1)
          bave%taupa3  = beach%taupa3(1)

          bave%taud_err   = 0.0d0
          bave%taupa3_err = 0.0d0

          bave%cumm2   = 0.0d0
          bave%cummpa3 = 0.0d0
        end if

        if (option%calc_pret) then

          bave%cumm    = 0.0d0

          do is = 0, option%nstate
            do js = 0, option%nstate

              bave%cumm(0, js, is) = 0.5d0 * option%dt_out &
                           * (bave%hist(0, js, is) + bave%hist(1, js, is))
              do istep = 1, nt_range - 1
                bave%cumm(istep, js, is) = bave%cumm(istep - 1, js, is) &
                  + 0.5d0 * option%dt_out &
                    * (bave%hist(istep, js, is) + bave%hist(istep + 1, js, is))
              end do
             
              bave%cumm(nt_range, js, is) = bave%cumm(nt_range - 1, js, is) &
                + option%dt_out * bave%hist(nt_range, js, is)
             
            end do
          end do

        end if

      end if


    end subroutine get_kinetics_average
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine output_kinetics(option, output, bave)
!-----------------------------------------------------------------------
      implicit none 

      type(s_option),  intent(in) :: option
      type(s_output),  intent(in) :: output
      type(s_bootave), intent(in) :: bave

      integer                :: iunit, istep, is, js, k
      real(8)                :: t
      character(len=MaxChar) :: fname
    

      ! PROP
      !
      write(fname,'(a,".prop")') trim(output%fhead)
      call open_file(fname, iunit)
      
      write(iunit,'("property      value          err")')
      if (option%calc_pret) then
        write(iunit,'("taud    ", 2f15.7)') bave%taud,   bave%taud_stdev
        write(iunit,'("kd      ", 2f15.7)') bave%kd,     bave%kd_stdev
      end if

      if (option%calc_kins) then 
        write(iunit,'("tauins  ", 2f15.7)') bave%tauins, bave%tauins_stdev
        write(iunit,'("kins    ", 2f15.7)') bave%kins,   bave%kins_stdev
      end if

      if (option%calc_2nd) then
        write(iunit,*)
        do is = 1, option%nstate
          if (is /= option%bound_id) then
            write(iunit,'("taud2 ",i0,1x,2f15.7)') is, bave%taud2(is), bave%taud2_err(is)
          end if
        end do

        write(iunit,'("taupa3  ", 2f15.7)') bave%taupa3, bave%taupa3_err
      end if

      close(iunit)

      ! Generate input for RP rate 
      !
      write(fname,'(a,".rprate_input")') trim(output%fhead)
      call open_file(fname, iunit)
    
      write(iunit,'("&option_param")')
      k = 1 
      if (option%calc_2nd) then
        k = 2 
      end if

      write(iunit,'("  rporder     = ", i0)') k
      write(iunit,'("  timeunit    = 1.0d-9 ! please modify it by yourself")')
      if (option%calc_kins) then
        write(iunit,'("  kins        = ", f15.7)') bave%kins
        write(iunit,'("  kins_err    = ", f15.7)') bave%kins_err
      else
        write(iunit,'("  kins        = VALUE1 ! please modify it by yourself")') 
        write(iunit,'("  kins_err    = VALUE2 ! please modify it by yourself")') 
      end if

      if (option%calc_pret) then
        write(iunit,'("  taud        = ", f15.7)') bave%taud
        write(iunit,'("  taud_err    = ", f15.7)') bave%taud_err
      else
        write(iunit,'("  taud        = VALUE1 ! please modify it by yourself")')
        write(iunit,'("  taud_err    = VALUE2 ! please modify it by yourself")')
      end if
      write(iunit,'("  Kstar       = VALUE1 ! please modify it by yourself")')
      write(iunit,'("  Kstar_err   = VALUE2 ! please modify it by yourself")')

      if (option%calc_2nd) then
        write(iunit,'("  taupa3      = ", f15.7)') bave%taupa3
        write(iunit,'("  taupa3_err  = ", f15.7)') bave%taupa3_err
      end if

      write(iunit,'("/")')

      close(iunit)

      if (.not. option%calc_pret) &
        return

      ! Pret(t)
      !
      write(fname,'(a,".pret")') trim(output%fhead)
      call open_file(fname, iunit)
      is = option%reaczone_id
      do istep = 0, option%nt_range
        t = istep * option%dt_out
        write(iunit,'(f15.7)', advance='no') t
        write(iunit,'(f15.7)', advance='no') bave%hist(istep, is, is)
        write(iunit,'(f15.7)', advance='no') bave%hist_err(istep, is, is)
        write(iunit,'(f15.7)', advance='no') bave%cumm(istep, is, is)

        if (option%calc_2nd) then
          do js = 1, option%nstate
            if (js /= option%bound_id) then
              write(iunit,'(f15.7)', advance='no') bave%hist2(istep, js)
              write(iunit,'(f15.7)', advance='no') bave%cumm2(istep, js)
            end if
          end do
        end if
        write(iunit,*)
      end do
      close(iunit)

      ! Transition probability matrix
      !
      if (option%out_tpm) then

        ! TCF
        !
        write(fname,'(a,".tpm")') trim(output%fhead)
        call open_file(fname, iunit)


        write(iunit,'("# col. ",i5," : time")') 1
        k = 1
        do is = 1, option%nstate
          do js = 1, option%nstate
            k = k + 1
            write(iunit,'("# col. ",i5," : state ",i5, " to ", i5)') k, is, js
          end do
        end do

        do istep = 0, option%nt_range
          write(iunit,'(f15.7)', advance='no') option%dt_out * istep
          do is = 1, option%nstate
            do js = 1, option%nstate
              write(iunit,'(f15.7)', advance='no') bave%hist(istep, js, is)
            end do 
          end do
          write(iunit,*)
        end do
        close(iunit)

        ! Running integral
        !
        write(fname,'(a,".tpm_cumm")') trim(output%fhead)
        call open_file(fname, iunit)

        write(iunit,'("# col. ",i5," : time")') 1
        k = 1
        do is = 1, option%nstate
          do js = 1, option%nstate
            k = k + 1
            write(iunit,'("# col. ",i5," : state ",i5, " <= ", i5)') k, js, is
          end do
        end do

        do istep = 0, option%nt_range
          write(iunit,'(f15.7)', advance='no') option%dt_out * istep
          do is = 1, option%nstate
            do js = 1, option%nstate
              write(iunit,'(f15.7)', advance='no') bave%cumm(istep, js, is)
            end do 
          end do
          write(iunit,*)
        end do

        close(iunit)

        ! Time constants
        !
        write(fname,'(a,".tpm_tau")') trim(output%fhead)
        call open_file(fname, iunit)

        write(iunit,'("# j i tau_{j<=i} Err")')
        do is = 1, option%nstate
          do js = 1, option%nstate
            write(iunit,'(i5,i5,2f15.7)') js, is, bave%tau_tpm(js, is), bave%tau_tpm_err(js, is) 
          end do
        end do

        close(iunit)

      end if


!      write(fname,'(a,".tpm_to_reac")') trim(output%fhead)
!      call open_file(fname, iunit)
!      is = option%reaczone_id
!      write(iunit,'("# col. ",i5," : time")') 1
!      k = 1 
!      do js = 0, option%nstate
!        if (js /= is) then
!          k = k + 1
!          write(iunit,'("# col. ",i5," : state ",i5, " to ", i5)') k, js, is
!          if (js == 0) then
!            write(iunit,'("#   Note: state 0 is a region which is not defined in fstate file")')
!          end if
!        end if
!      end do
!
!      do istep = 0, option%nt_range
!        write(iunit,'(f15.7)', advance='no') option%dt_out * istep
!        do js = 0, option%nstate
!          if (js /= is) then 
!            write(iunit,'(f15.7)', advance='no') bave%hist(istep, is, js)
!          end if
!        end do
!        write(iunit,*)
!      end do
!
!      close(iunit)
!
!      write(fname,'(a,".tpm_from_reac")') trim(output%fhead)
!      call open_file(fname, iunit)
!      is = option%reaczone_id
!      write(iunit,'("# col. ",i5," : time")') 1
!      k = 1 
!      do js = 0, option%nstate
!        if (js /= is) then
!          k = k + 1
!          write(iunit,'("# col. ",i5," : state ",i5, " to ", i5)') k, is, js
!          if (js == 0) then
!            write(iunit,'("#   Note: state 0 is a region which is not defined in fstate file")')
!          end if
!        end if
!      end do
!
!      do istep = 0, option%nt_range
!        write(iunit,'(f15.7)', advance='no') option%dt_out * istep
!        do js = 0, option%nstate
!          if (js /= is) then 
!            write(iunit,'(f15.7)', advance='no') bave%hist(istep, js, is)
!          end if
!        end do
!        write(iunit,*)
!      end do
!
!      close(iunit)

      ! 2nd-order transition
      !
      if (option%calc_2nd) then
        write(fname,'(a,".pret2nd")') trim(output%fhead)
        call open_file(fname, iunit)
        do istep = 0, option%nt_range
          write(iunit,'(f15.7)', advance='no') option%dt_out * istep
          write(iunit,'(f15.7)', advance='no') bave%pa3(istep)
          write(iunit,'(f15.7)', advance='no') bave%cummpa3(istep)
          write(iunit,*)
        end do
        close(iunit)
      end if

    end subroutine output_kinetics
!-----------------------------------------------------------------------


end module mod_analyze
!=======================================================================
