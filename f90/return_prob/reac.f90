!=======================================================================
module mod_reac
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

  ! constants
  !

  ! structures
  !

  ! subroutines
  !
  public :: get_reaction_info
  public :: check_transient
  public :: state_average
  public :: search_reaction

  contains

!-----------------------------------------------------------------------
    subroutine get_reaction_info(option, state)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      type(s_state),  intent(inout) :: state(:)

      integer :: ifile, imol, nmol, nfile
      integer :: nBevents, nUBevents
      
      integer, allocatable :: nstep_sparse(:)


      if (.not. option%calc_kins) &
        return

      nmol  = option%nmol
      nfile = size(state(:))

      if (option%kins_mode == KinsModeCOARSE) then
        write(iw,'("Get_Reaction_Info> kins_mode = COARSE")')
        write(iw,'(">> reaction events are counted by introducing coarse-grained timescale")')
        ! for local 
        !
        allocate(nstep_sparse(nfile))

        ! Average state
        !
!$omp parallel private(ifile), default(shared)
!$omp do
        do ifile = 1, nfile
       
          call state_average(option%nt_transient, &
                             nmol,                &
                             state(ifile),        &
                             nstep_sparse(ifile))
       
        end do
!$omp end do
!$omp end parallel

        ! Search reaction
        !
!$omp parallel private(ifile, nBevents, nUBevents),  &
!$omp          shared(state),                        &
!$omp          default(shared)
!$omp do
        do ifile = 1, nfile
       
          call search_reaction(option,              &
                               nstep_sparse(ifile), &
                               state(ifile),        &
                               nBevents,            &
                               nUBevents)
       
          state(ifile)%nBevents   = nBevents
          state(ifile)%nUBevents  = nUBevents
       
        end do
!$omp end do
!$omp end parallel

       
        deallocate(nstep_sparse)
      else if (option%kins_mode == KinsModeQUENCH) then
        write(iw,'("Get_Reaction_Info> kins_mode = QUENCH")')
        write(iw,'(">> reaction events are counted as the number of &
                   &trajectories visiting quench state")')
        do ifile = 1, nfile
          state(ifile)%nBevents  = 0
          state(ifile)%nUBevents = 0
          do imol = 1, nmol
            if (state(ifile)%quench_step(imol) < state(ifile)%nstep) then
              state(ifile)%nBevents = state(ifile)%nBevents + 1
            end if
          end do
        end do
        
      end if

      write(iw,*)
      write(iw,'("Get_Reaction_Info> Detected reaction events")')
      write(iw,'("Traj.    Binds    Unbinds")')
      do ifile = 1, nfile
        write(iw,'(i6)', advance='no') ifile
        write(iw,'(i6)', advance='no') state(ifile)%nBevents
        write(iw,'(i6)')               state(ifile)%nUBevents
      end do
      write(iw,*)

    end subroutine get_reaction_info
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine check_transient(output, option, state)
!-----------------------------------------------------------------------
      implicit none

      type(s_output), intent(in)    :: output
      type(s_option), intent(in)    :: option
      type(s_state),  intent(inout) :: state(:)

      integer                :: ntr_sta, ntr_interval, nntr
      integer                :: nfile, iunit
      integer                :: ifile, itr, ntr, istep, imol, nmol
      integer                :: nstep_sparse
      real(8)                :: dt_sparse
      character(len=MaxChar) :: fname

      integer, allocatable :: nBevents(:, :), nUBevents(:, :)
      integer, allocatable :: nstay(:, :)
      real(8), allocatable :: krate(:)


      ! Setup parameters
      !
      ntr_sta      = option%ntr_sta
      ntr_interval = option%ntr_interval
      nntr         = option%nntr
      nmol         = option%nmol
      nfile        = size(state(:)) 

      ! Memory allocation
      !
      allocate(nBevents(0:nntr, nfile))
      allocate(nUBevents(0:nntr, nfile))
      allocate(nstay(0:nntr, nfile))
      allocate(krate(0:nntr))
      
      nBevents  = 0
      nUBevents = 0

!$omp parallel private(ifile, itr, ntr, nstep_sparse),          &
!$omp          shared(nmol, state, nBevents, nUBevents, nstay), &
!$omp          default(shared)
!$omp do
      do ifile = 1, nfile
        do itr = 0, nntr
          ntr = ntr_sta + itr * ntr_interval

          ! Average states using moving average scheme
          !
          call state_average(ntr, nmol, state(ifile), nstep_sparse)

          ! Search Binding/Unbinding
          !
          call search_reaction(option,                &
                               nstep_sparse,          &
                               state(ifile),          &
                               nBevents(itr, ifile),  &
                               nUBevents(itr, ifile))

          ! Calculate time in which complex stays in Reaction zone
          !
          call calculate_staytime(nmol, state(ifile),  &
                                  option%reaczone_id,  &
                                  nstay(itr, ifile))

        end do
      end do
!$omp end do
!$omp end parallel

      ! Calculate 1st-order rate constant (krate) 
      !
      do itr = 0, nntr
        krate(itr) = dble(sum(nBevents(itr, :)))        &
                     / (sum(nstay(itr, :)) * option%dt)
      end do

      ! Output Dt-dependency of krate
      !
      write(fname, '(a,".kcheck")') trim(output%fhead)
      call open_file(fname, iunit)
      write(iunit,'("# col-1 : time")')
      write(iunit,'("# col-2 : kins")')
      write(iunit,'("# col-3 : tauins")')
      do itr = 0, nntr
        ntr       = ntr_sta + itr * ntr_interval
        dt_sparse = ntr * option%dt

        write(iunit,'(f15.7,2x)', advance='no') dt_sparse
        write(iunit,'(f15.7)',    advance='no') krate(itr)
        write(iunit,'(f15.7)')          1.0d0 / krate(itr)

      end do
      close(iunit)

      ! Memory Deallocation
      !
      deallocate(nBevents)
      deallocate(nUBevents)
      deallocate(nstay)
      deallocate(krate)

    end subroutine check_transient
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine state_average(ntr, nmol, state, nstep_sparse)
!-----------------------------------------------------------------------
      implicit none

      integer,        intent(in)    :: ntr
      integer,        intent(in)    :: nmol
      type(s_state),  intent(inout) :: state
      integer,        intent(out)   :: nstep_sparse

      integer :: istep, jstep, nstep
      integer :: imin, imax, imol
      real(8) :: r


      state%avedata = 0.0d0
      nstep         = state%nstep

      do imol = 1, nmol
        jstep = 0 
        do istep = ntr, nstep - ntr, ntr
          jstep = jstep + 1
       
          imin  = istep - ntr / 2
          imax  = istep + ntr / 2
          r     = sum(state%data_ts(imin:imax, imol)) / dble(ntr + 1)
       
          state%avedata(jstep, imol) = r
        end do

        imin = nstep - ntr / 2
        imax = nstep
        r    = sum(state%data_ts(imin:imax, imol)) / dble(ntr /2 + 1)
       
        state%avedata(jstep + 1, imol) = r
        nstep_sparse                   = jstep + 1
      end do

    end subroutine state_average
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine search_reaction(option,       &
                               nstep_sparse, &
                               state,        &
                               nBevents,     &
                               nUBevents)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      integer,        intent(in)    :: nstep_sparse
      type(s_state),  intent(inout) :: state
      integer,        intent(out)   :: nBevents
      integer,        intent(out)   :: nUBevents

      integer :: istep, jstep, imol, nmol
      integer :: nBevents_prev, diff_nBevents
      real(8) :: rcurr, rnext
      logical :: is_changed
    

      nmol          = option%nmol
                    
      nBevents      = 0
      nUBevents     = 0

      nBevents_prev = 0
      do imol = 1, nmol 
        istep      = 0
        is_changed = .true.
        rcurr      = state%data_ts(1, imol)
        do while (is_changed)
          is_changed = .false.
          do jstep = istep + 1, nstep_sparse
            rnext = state%avedata(jstep, imol)
       
            if (rcurr < 0.0d0 .and. rnext >= 0.0d0) then
              nBevents   = nBevents + 1
              rcurr      = rnext
              istep      = jstep
              is_changed = .true.
              exit
            else if (rcurr >= 0.0d0 .and. rnext < 0.0d0) then
              nUBevents  = nUBevents + 1
              rcurr      = rnext
              istep      = jstep
              is_changed = .true.
              exit
            end if
          end do
       
        end do

        diff_nBevents = nBevents - nBevents_prev

        ! Judge whethre this trajectory includes the reactions
        ! (if use_reactraj = .true., the trajectory is 
        !  regarded as unreactive)
        !
        if (option%use_reactraj) then
          state%is_reacted(imol) = .false.
        else
          state%is_reacted(imol) = .false.
          if (diff_nBevents > 0) then
            state%is_reacted(imol) = .true.
          end if
        end if

        nBevents_prev = nBevents

      end do

    end subroutine search_reaction
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine calculate_staytime(nmol, state, target_id, nstay)
!-----------------------------------------------------------------------
      implicit none

      integer,       intent(in)  :: nmol
      type(s_state), intent(in)  :: state
      integer,       intent(in)  :: target_id
      integer,       intent(out) :: nstay

      integer :: istep, imol, nstep

      
      ! Setup parameters
      !
      nstep = state%nstep

      ! Calculate
      !
      nstay = 0
      do imol = 1, nmol
        do istep = 1, nstep
          if (state%data(istep, imol) == target_id) then
            nstay = nstay + 1
          end if
        end do
      end do

    end subroutine calculate_staytime
!-----------------------------------------------------------------------

end module mod_reac
!=======================================================================
