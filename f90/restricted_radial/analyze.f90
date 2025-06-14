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
  real(8), parameter :: c0    = 1.0d0   ! standard conc. [1 M] 
  real(8), parameter :: c0inv = 1.0d0   ! inverse of standard conc. [1 M^-1]
  real(8), parameter :: EPS   = 1.0d-20

  ! structures
  !
  type :: s_state
    integer :: nstep
    integer :: ndim
    logical :: is_reacted = .false.
    integer, allocatable :: data(:) 
    real(8), allocatable :: hist(:)
  end type s_state

  type :: s_booteach
    integer :: ntrial
    integer :: nsample
    integer :: ngrid 
    real(8), allocatable :: Keq(:)
    integer, allocatable :: rand(:, :)
    real(8), allocatable :: pmf(:, :), rpmf(:, :)
    real(8), allocatable :: rdf(:, :), rrdf(:, :)
  end type s_booteach

  type :: s_bootave
    integer :: ngrid
    real(8) :: Keq, Keq_stdev, Keq_err
    real(8),           allocatable :: rdf(:), rdf_err(:) 
    real(8),           allocatable :: rrdf(:),rrdf_err(:)  
    character(len=20), allocatable :: pmf(:), pmf_err(:) 
    character(len=20), allocatable :: rpmf(:), rpmf_err(:) 
  end type s_bootave

  ! subroutines
  !
  public  :: analyze
  private :: analyze_bootstrap 

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


      integer :: ifile, istep, ig, jg
      integer :: nreact, nstay
      real(8) :: rpinit, kr, t, dth
      real(8) :: tau_r, tau_d

      integer :: ngrid, ncount
      real(8) :: fourpi, fact, dr, r, gmax
      real(8) :: pop_f, pop_b

      real(8) :: vol, beta, kT
      real(8) :: vu, vb, wu, wb, dw
      real(8) :: dGpmf, dGv, dGzero, Keq0, Keq
      real(8) :: dG3d
      real(8) :: weight, weight_sum, weight_sum_bound
      real(8) :: p_shift, rp_shift, pmin
      logical :: use_weight

      real(8), allocatable :: rdf(:), rrdf(:)
      real(8), allocatable :: pmf(:), rpmf(:)
      character(len=MaxChar) :: frd


!$    write(iw,*)
!$    write(iw,'("Analyze> OpenMP parallelization is activated")')
      allocate(cv(input%ncv), state(input%ncv))
      allocate(rdf(0:option%ngrid), rrdf(0:option%ngrid)) 
      allocate(pmf(0:option%ngrid), rpmf(0:option%ngrid)) 

      ! read cv files
      !
      write(iw,*)
      write(iw,'("Analyze> Read CV file")')

      if (option%only_r) then
        do ifile = 1, input%ncv
          call read_cv(input%fcv(ifile), option%ndim, cv(ifile))
        end do
      else
        do ifile = 1, input%ncv
          call read_cv(input%fcv(ifile), option%ndim, cv(ifile), plus = 1)
        end do
      end if

      use_weight = .false.
      if (trim(input%fweight(1)) /= '') then
        use_weight = .true.
        do ifile = 1, input%ncv
          call get_weight(input, option, ifile, cv(ifile))
        end do
      else
        do ifile = 1, input%ncv
          allocate(cv(ifile)%weight(cv(ifile)%nstep))
          cv(ifile)%weight = 1.0d0
        end do
      end if 

      ! get state
      !
      write(iw,*)
      write(iw,'("Analyze> Get State")')
!$omp parallel private(ifile), default(shared)
!$omp do
      do ifile = 1, input%ncv
        call get_state(option%ndim, option%state_def, &
                       cv(ifile), state(ifile))
      end do
!$omp end do
!$omp end parallel 

      fourpi = 4.0d0 * Pi

      if (option%use_bootstrap) then
        if (use_weight) then
          write(iw,'("Analyze> Error.")')
          write(iw,'("Sorry, bootstrap is not supported when weight file is used.")')
          stop
        end if
        call analyze_bootstrap(option, output, bootopt, cv, state, input%ncv)
      else
        
        ngrid            = option%ngrid
        dr               = option%dr
                         
        rdf              = 0.0d0
        rrdf             = 0.0d0
        ncount           = 0
        pop_b            = 0.0d0
        pop_f            = 0.0d0
        weight           = 0.0d0
        weight_sum       = 0.0d0
        weight_sum_bound = 0.0d0
        do ifile = 1, input%ncv
          do istep = option%nsta, cv(ifile)%nstep

            if (option%only_r) then
              r = cv(ifile)%data(1, istep)
            else
              r = cv(ifile)%data(option%ndim + 1, istep)
            end if

            ig     = nint(r / option%dr)
            weight = cv(ifile)%weight(istep)
            ncount = ncount  + 1
            pop_f  = pop_f   + weight 
            weight_sum = weight_sum + weight 
            if (ig >= 0 .and. ig <= ngrid) then
              rdf(ig) = rdf(ig) + weight
              if (state(ifile)%data(istep) == REACTIVE) then
                rrdf(ig)         = rrdf(ig)         + weight 
                pop_b            = pop_b            + weight 
                weight_sum_bound = weight_sum_bound + weight 
              end if
            end if
          end do
        end do

        ! consider jacobian 
        !
        do ig = 1, ngrid
          r       = ig * dr
          fact    = 1.0d0 / (fourpi * r * r * dr)
          rdf(ig)  = rdf(ig)  * fact
          rrdf(ig) = rrdf(ig) * fact 
        end do

        fact = 1.0d0 / (fourpi/3.0d0*dr*dr*dr)
        rdf(0)  = rdf(0)  * fact
        rrdf(0) = rrdf(0) * fact 

        if (option%refs_system) then

          vol    = option%box_ref(1) * option%box_ref(2) &
                 * option%box_ref(3)
          vb     = vol * (pop_b / pop_f)
          kT     = boltz * option%temperature
          beta   = 1.0d0 / kT 
          dGzero = - kT * log(vb / option%vol0) 
!
          write(iw,'("... Thermodynamic analysis ...")')
          write(iw,'("vb      (A^3)      = ",f15.7)') vb 
          write(iw,'("dGzero  (kcal/mol) = ",f15.7)') dGzero

!          vol  = option%box_ref(1) * option%box_ref(2) * option%box_ref(3)
!          rdf  = vol * rdf  / dble(ncount)
!          rrdf = vol * rrdf / dble(ncount)
!
!          kT   = boltz * option%temperature
!          beta = 1.0d0 / kT 
!          do ig = 0, ngrid
!            pmf(ig)   = - kT * log(rdf(ig))
!            rpmf(ig)  = - kT * log(rrdf(ig))
!            if (rdf(ig) <= EPS) then
!              pmf(ig)  = 0.0d0 
!            end if
!            
!            if (rrdf(ig) <= EPS) then
!              rpmf(ig) = 0.0d0
!            end if 
!         
!          end do
!
!          if (option%calcfe) then
!            beta = 1.0d0 / (boltz * option%temperature)
!            ! vb
!            vb = 0.0d0
!            do ig = 0, ngrid
!              r  = ig * dr
!              vb = vb + fourpi * r * r * dr * rrdf(ig) 
!            end do
!
!            dGzero = - kT * log(vb / option%vol0) 
!
!            write(iw,'("... Thermodynamic analysis ...")')
!            write(iw,'("vb      (A^3)      = ",f15.7)') vb 
!            write(iw,'("dGzero  (kcal/mol) = ",f15.7)') dGzero
!
!          end if
!
        else

          ! search maximam and normalize
          !
          gmax = maxval(rrdf)
          rdf  = rdf  / gmax
          rrdf = rrdf / gmax
         
          kT   = boltz * option%temperature
          beta = 1.0d0 / kT 
          do ig = 0, ngrid
           pmf(ig)   = - kT * log(rdf(ig))
           rpmf(ig)  = - kT * log(rrdf(ig))
           if (rdf(ig) <= EPS) then
             pmf(ig)  = 0.0d0 
           end if
         
           if (rrdf(ig) <= EPS) then
             rpmf(ig) = 0.0d0
           end if 
         
          end do
         
          ! calculate standard free energy change (1 M) 
          !
          if (option%calcfe) then
            beta = 1.0d0 / (boltz * option%temperature)
            ! vb
            vb = 0.0d0
            do ig = 0, ngrid
              r  = ig * dr
              vb = vb + fourpi * r * r * dr * rrdf(ig) 
            end do
         
            ! wb
            jg = maxloc(rrdf,1) - 1
            wb = rpmf(jg)
         
            ! wu
            wu = 0.0d0
            jg = 0
            do ig = 0, ngrid
              r = ig * dr
              if (r >= option%urange(1) .and. r <= option%urange(2)) then
                jg = jg + 1
                wu = wu + pmf(ig)
              end if
            end do
            wu = wu / dble(jg)
         
            dw     =  wu - wb
            dGv    = -kT * log(vb/option%vol0)
            dGzero = -dw + dGv
            Keq    =  c0inv * exp(-beta*dGzero)
         
            write(iw,'("... Thermodynamic analysis ...")')
            write(iw,'("vb      (A^3)      = ",f15.7)') vb 
            write(iw,'("wu      (kcal/mol) = ",f15.7)') wu
            write(iw,'("wb      (kcal/mol) = ",f15.7)') wb
            write(iw,'("dw      (kcal/mol) = ",f15.7)') dw
            write(iw,'("dGv     (kcal/mol) = ",f15.7)') dGv
            write(iw,'("dGzero  (kcal/mol) = ",f15.7)') dGzero
            !write(iw,'("Keq0    (-)        = ",es15.7)') Keq0
            write(iw,'("Keq     (M^-1)     = ",es15.7)') Keq
          end if
         
          write(frd, '(a,".rd")') trim(output%fhead)
          
          pmin = minval(pmf)
          open(10,file=trim(frd))
            do ig = 0, ngrid 
              r = dr * ig
              p_shift  = pmf(ig)  - pmin
              rp_shift = rpmf(ig) - pmin 
              write(10,'(f20.10,a20,a20)') r, &
                conv_f2a(rdf(ig), p_shift, EPS), conv_f2a(rrdf(ig), rp_shift, EPS)
              !write(10,'(f20.10,a20,a20)') r, &
              !  conv_f2a(rdf(ig), pmf(ig), EPS), conv_f2a(rrdf(ig), rpmf(ig), EPS)
              !write(10,'(f20.10,f20.10,f20.10)') r, &
              !  pmf(ig), rpmf(ig)
            end do
          close(10)

        end if

      end if
      !if (option%btstrp) then
      !  write(bts%fbts_Keq, '(a,".Keq")') trim(output%fhead_btstrp)
      !  write(bts%fbts_pmf, '(a,".pmf")') trim(output%fhead_btstrp)
      !  call analyze_btstrp(option, state, cvinfo%nfile, bts)
      !end if
      deallocate(rdf, rrdf)

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
      allocate(state%data(nstep))

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
      
      state%nstep = nstep
      state%ndim  = ndim

    end subroutine get_state 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    function conv_f2a(gval, pval, thr)
!-----------------------------------------------------------------------
      implicit none

      real(8), intent(in) :: gval, pval, thr
      character(len=20)   :: conv_f2a

      character(len=20)   :: infi 


      infi = '            Infinity'

      if (gval <= thr) then
        conv_f2a = infi
      else
        write(conv_f2a,'(f20.10)') pval 
      end if

    end function conv_f2a
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine analyze_bootstrap(option, output, bootopt, cv, state, nfile)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),  intent(in)    :: option
      type(s_output),  intent(in)    :: output 
      type(s_bootopt), intent(inout) :: bootopt
      type(s_cv),      intent(in)    :: cv(:)
      type(s_state),   intent(in)    :: state(:)
      integer,         intent(in)    :: nfile

      ! local variables
      !
      type(s_booteach) :: beach
      type(s_bootave)  :: bave

      integer :: ig, itrial
      integer :: ntrial, ngrid, nsample

      character(len=MaxChar) :: fout


      ! Prepare parameters
      !
      ntrial  = bootopt%ntrial
      nsample = bootopt%nsample
      ngrid   = option%ngrid

      beach%ntrial  = ntrial
      beach%nsample = nsample
      beach%ngrid   = ngrid

      ! Allocate memory
      !
      allocate(beach%rand(nsample, ntrial))
      allocate(beach%Keq(ntrial))
      allocate(beach%pmf(0:ngrid, ntrial))
      allocate(beach%rdf(0:ngrid, ntrial))
      allocate(beach%rrdf(0:ngrid, ntrial))
      allocate(beach%rpmf(0:ngrid, ntrial))

      allocate(bave%rdf(0:ngrid), bave%rdf_err(0:ngrid))
      allocate(bave%rrdf(0:ngrid), bave%rrdf_err(0:ngrid))
      allocate(bave%pmf(0:ngrid), bave%pmf_err(0:ngrid))
      allocate(bave%rpmf(0:ngrid), bave%rpmf_err(0:ngrid))

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

      ! Calculate quantities for each trial
      !
      call calc_fe_bootstrap(option, state, cv, beach)

      ! Calculate average
      !
      call calc_feave_bootstrap(option, beach, bave)


      ! Output
      !
      write(fout, '(a,".pmf")') trim(output%fhead)
      open(UnitOUT, file=trim(fout))
        do ig = 0, option%ngrid
          write(UnitOUT,'(f20.10,4a20)')         &
                ig * option%dr,                  &
                bave%pmf(ig),  bave%pmf_err(ig), &
                bave%rpmf(ig), bave%rpmf_err(ig)
        end do  
      close(UnitOUT)

      if (option%calcfe) then
        write(fout, '(a,".Keq")') trim(output%fhead)
        open(UnitOUT, file=trim(fout))
      write(UnitOUT,'(a20,2x,a15,2x,a15,2x,a15)') &
              "Property", "Value", "SD", "SE" 
      write(UnitOUT,'(a20,2x,3(es15.7,2x))') &
              "Keq      (M^-1)    :", bave%Keq, bave%Keq_stdev, bave%Keq_err
        !write(UnitOUT,'("Keq    (M^-1)     : ", 2es15.7)') &
        !        bave%Keq, bave%Keq_err 
        close(UnitOUT)
      end if


    end subroutine analyze_bootstrap
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine calc_fe_bootstrap(option, state, cv, beach)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_state),    intent(in)    :: state(:)
      type(s_cv),       intent(in)    :: cv(:) 
      type(s_booteach), intent(inout) :: beach 

      integer :: istep, ig, jg, itraj, itrial, isample
      integer :: ntrial, nsample, ngrid, nstep
      real(8) :: kT, beta, dr
      real(8) :: r, fact, gmax
      real(8) :: wb, vb, wu
      real(8) :: dw, dGv, dGzero, Keq

      real(8), allocatable :: rdf(:), rrdf(:)
      real(8), allocatable :: pmf(:), rpmf(:)


      ! Prepare parameters
      !
      ntrial  = beach%ntrial
      nsample = beach%nsample
      ngrid   = beach%ngrid 

      dr      = option%dr
      kT      = boltz * option%temperature
      beta    = 1.0d0 / kT

      ! Allocate memory
      !
      allocate(rdf(0:ngrid), rrdf(0:ngrid))
      allocate(pmf(0:ngrid), rpmf(0:ngrid))

      beach%Keq  = 0.0d0
      beach%pmf  = 0.0d0 
      beach%rdf  = 0.0d0
      beach%rrdf = 0.0d0

      rdf        = 0.0d0
      rrdf       = 0.0d0
      pmf        = 0.0d0
      rpmf       = 0.0d0

      write(iw,'("Calc_Fe_Bootstrap> Start MC-Bootstrap")')

      !$omp parallel private(itrial, rdf, rrdf, pmf, rpmf, &
      !$omp                  isample, nstep, r, ig, fact,  &
      !$omp                  istep, itraj, gmax),          &
      !$omp          default(shared)
      !$omp do
      do itrial = 1, ntrial

        if (mod(itrial, 10) == 0) then
          write(iw,'("trial : ",i0)') itrial
        end if

        do isample = 1, nsample
          itraj = beach%rand(isample, itrial)
          nstep = cv(itraj)%nstep
          do istep = option%nsta, nstep

            if (option%only_r) then
              r  = cv(itraj)%data(1, istep)
            else
              r  = cv(itraj)%data(option%ndim + 1, istep)
            end if
            ig = nint(r / dr)

            if (ig >= 0 .and. ig <= ngrid) then
              rdf(ig) = rdf(ig) + 1.0d0
              if (state(itraj)%data(istep) == REACTIVE) then
                rrdf(ig) = rrdf(ig) + 1.0d0
              end if
            end if

          end do
        end do

        do ig = 1, ngrid
          r        = ig * dr
          fact     = 1.0d0 / (4.0d0 * PI * r * r * dr)
          rdf(ig)  = rdf(ig)  * fact
          rrdf(ig) = rrdf(ig) * fact
        end do

        fact    = 1.0d0 / ((4.0d0/3.0d0) * PI * dr * dr * dr)
        rdf(0)  = rdf(0)  * fact 
        rrdf(0) = rrdf(0) * fact

        gmax = maxval(rrdf)
        rdf  = rdf  / gmax
        rrdf = rrdf / gmax

        do ig = 0, ngrid
          if (rdf(ig) <= EPS) then
             pmf(ig) = 0.0d0
          else
             pmf(ig)  = - kT * log(rdf(ig))
          end if

          if (rrdf(ig) <= EPS) then
            rpmf(ig) = 0.0d0
          else
            rpmf(ig) = - kT * log(rrdf(ig))
          end if
        end do

        beach%rdf(:,  itrial) = rdf(:)
        beach%rrdf(:, itrial) = rrdf(:)
        beach%pmf(:,  itrial) = pmf(:)
        beach%rpmf(:, itrial) = rpmf(:) 
      end do
      !$omp end do
      !$omp end parallel

      if (option%calcfe) then
        do itrial = 1, ntrial
          ig = maxloc(beach%rrdf(:, itrial), 1) - 1

          ! wb
          wb = beach%rpmf(ig, itrial)

          ! vb, wu
          vb = 0.0d0
          wu = 0.0d0
          jg = 0
          do ig = 0, option%ngrid
            r  = ig * dr
            vb = vb + 4.0d0 * PI * r * r * dr * beach%rrdf(ig, itrial)

            if (r >= option%urange(1) .and. r <= option%urange(2)) then
              jg = jg + 1
              wu = wu + beach%pmf(ig, itrial) 
            end if 
          end do
          wu = wu / dble(jg)

          ! dw, dGv, dGzero, Keq
          dw          = wu - wb
          dGv         = -kT * log(vb / option%vol0)
          dGzero      = -dw + dGv

          beach%Keq(itrial) = exp(-beta * dGzero)

        end do

      end if

      !write(iw,'(" vb  = ", f15.7)') vb 
      !write(iw,'(" wu  = ", f15.7)') wu
      !write(iw,'(" dw  = ", f15.7)') dw 
      !write(iw,'(" dGv = ", f15.7)') dGv
      !write(iw,'(" dGz = ", f15.7)') dGzero
      !do ig = 0, ngrid
      !  write(99,'(2f20.10)') ig * dr, beach%rpmf(ig, 1) 
      !end do

      deallocate(rdf, rrdf)
      deallocate(pmf, rpmf)

    end subroutine calc_fe_bootstrap
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine calc_feave_bootstrap(option, beach, bave)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_booteach), intent(in)    :: beach
      type(s_bootave),  intent(inout) :: bave

      integer :: itrial, ntrial, ngrid
      integer :: ig, ncount
      real(8) :: kT, beta
      real(8) :: Keq
      real(8) :: val, dev, sterr

      real(8), allocatable :: rdf(:), rrdf(:)
      real(8), allocatable :: pmf(:), rpmf(:)
      real(8), allocatable :: pmf_err(:), rpmf_err(:)


      ! Prepare parameters
      !
      ntrial = beach%ntrial
      ngrid  = beach%ngrid

      kT   = boltz * option%temperature
      beta = 1.0d0 / kT 

      ! Allocate memory
      !
      allocate(rdf(0:ngrid))
      allocate(rrdf(0:ngrid))
      allocate(pmf(0:ngrid))
      allocate(rpmf(0:ngrid))
      allocate(pmf_err(0:ngrid))
      allocate(rpmf_err(0:ngrid))

      ! Calculate Average & Error
      !
      ! - rdf
      rdf  = 0.0d0
      rrdf = 0.0d0
      do ig = 0, ngrid
        do itrial = 1, ntrial 
          rdf(ig)  = rdf(ig)  + beach%rdf(ig, itrial)
          rrdf(ig) = rrdf(ig) + beach%rrdf(ig, itrial) 
        end do
        rdf(ig)  = rdf(ig)  / ntrial
        rrdf(ig) = rrdf(ig) / ntrial 
      end do

      pmf      = 0.0d0
      rpmf     = 0.0d0
      pmf_err  = 0.0d0
      rpmf_err = 0.0d0
      do ig = 0, ngrid
        ! - pmf (average)
        if (rdf(ig) <= EPS) then
          pmf(ig)      = 0.0d0
        else
          pmf(ig) = - kT * log(rdf(ig)) 
        end if
        bave%pmf(ig) = conv_f2a(rdf(ig), pmf(ig), EPS)

        ! - pmf (error)
        ncount = 0
        dev    = 0.0d0
        do itrial = 1, ntrial
          if (beach%rdf(ig, itrial) > EPS) then
            ncount = ncount + 1
            val    = -kT * log(beach%rdf(ig, itrial)) - pmf(ig)
            dev    = dev + val * val
          end if
        end do
        sterr = dev / ((ncount-1)*ncount)
        sterr = sqrt(sterr)

        bave%pmf_err(ig) = conv_f2a(rdf(ig), sterr, EPS) 

        ! - rpmf (average)
        if (rrdf(ig) <= EPS) then
          rpmf(ig)      = 0.0d0
        else
          rpmf(ig) = - kT * log(rrdf(ig)) 
        end if 
        bave%rpmf(ig) = conv_f2a(rrdf(ig), rpmf(ig), EPS)

        ! - rpmf (error)
        ncount = 0
        dev    = 0.0d0
        do itrial = 1, ntrial
          if (beach%rrdf(ig, itrial) > EPS) then
            ncount = ncount + 1
            val    = -kT * log(beach%rrdf(ig, itrial)) - rpmf(ig)
            dev    = dev + val * val
          end if
        end do
        sterr = dev / ((ncount-1)*ncount)
        sterr = sqrt(sterr)

        bave%rpmf_err(ig) = conv_f2a(rrdf(ig), sterr, EPS)

      end do

      ! - Keq
      if (option%calcfe) then
        Keq      = sum(beach%Keq(:)) / ntrial
        bave%Keq = Keq

        dev = 0.0d0
        do itrial = 1, ntrial
          val = beach%Keq(itrial) - Keq 
          dev = dev + val * val
        end do 
        sterr = dev / ((ntrial-1)*ntrial)
        sterr = sqrt(sterr)
        
        bave%Keq_err   = sterr
        bave%Keq_stdev = sqrt(dev / (ntrial - 1))
      end if 

      deallocate(rdf, rrdf)
      deallocate(pmf, rpmf)
      deallocate(pmf_err, rpmf_err)

    end subroutine calc_feave_bootstrap
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_weight(input, option, itraj, cv)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in)    :: input 
      type(s_option), intent(in)    :: option
      integer,        intent(in)    :: itraj
      type(s_cv),     intent(inout) :: cv


      if (allocated(cv%weight)) then
        deallocate(cv%weight)
      end if

      call read_cv_weight(input%fweight(itraj), cv)

    end subroutine get_weight
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
