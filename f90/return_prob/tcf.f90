!=======================================================================
module mod_tcf
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
  public :: get_transtcf

  contains

!-----------------------------------------------------------------------
    subroutine get_transtcf(option, state) 
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      type(s_state),  intent(inout) :: state(:)

      integer :: ifile, istep, jstep, kstep, imol, it, jt
      integer :: is, js, ks
      integer :: nmol, nstep, nt_range, nt, nt2, nstate, nfile
      integer :: id_b, id_r
      integer :: nthreads, thread_id, nwork, ifinish, iper
      real(8) :: progress

      real(8), allocatable :: hist(:, :, :), norm(:, :)
      real(8), allocatable :: hist2(:, :)
      real(8), allocatable :: hsum(:)


      nmol     = option%nmol
      nstep    = state(1)%nstep
      nt_range = option%nt_range
      nstate   = option%nstate
      nfile    = size(state(:))

      id_b     = option%bound_id
      id_r     = option%reaczone_id

      ! Memory allocation (if needed) 
      !
      do ifile = 1, nfile

        ! for 1st-order calc.
        !
        if (.not. allocated(state(ifile)%hist)) then
          allocate(state(ifile)%hist(0:nt_range, 0:nstate, 0:nstate))
          allocate(state(ifile)%norm(0:nt_range, 0:nstate))
          state(ifile)%hist = 0.0d0
        end if

        ! for 2nd-order calc.
        !
        if (option%calc_2nd) then
          if (.not. allocated(state(ifile)%hist2)) then
            allocate(state(ifile)%hist2(0:nt_range, 0:nstate))
            state(ifile)%hist2 = 0.0d0
          end if
        end if
      end do

      if (.not. option%calc_pret) &
        return

      ! work space
      !
      allocate(hist(0:nt_range, 0:nstate, 0:nstate))
      allocate(norm(0:nt_range, 0:nstate))
      allocate(hsum(0:nstate))
      hist  = 0.0d0
      norm  = 0.0d0
      hsum  = 0.0d0


      ! Calculate tcf
      !   Remark: if use_reactraj = .true., is_reacted is always .false.
      !
      write(iw,*)
      write(iw,'("Get_TransTCF> Calculate 1st-order transition")')

!$omp parallel private(ifile, imol, istep, jstep, it, is, js, nt, hist, norm, &
!$omp                  nthreads, thread_id, nwork, progress, ifinish),        &
!$omp          default(shared)
      nthreads  = omp_get_num_threads()
      thread_id = omp_get_thread_num()
      nwork     = nfile / nthreads + 1
      ifinish   = 0
!$omp do schedule(dynamic)
      do ifile = 1, nfile
        
        if (thread_id == 0) then
          ifinish  = ifinish + 1
          progress = ifinish / dble(nwork) * 100.0d0
          write(iw,'("  Progress : ", f6.2, "%")') progress
        end if

        hist = 0.0d0
        norm = 0.0d0

        do imol = 1, nmol

          if (.not. state(ifile)%is_reacted(imol)) then
            do istep = 1, min(nstep, state(ifile)%read_step(imol))
           
              is = state(ifile)%data(istep, imol)
              !nt = min(istep + nt_range,               &
              !         state(ifile)%quench_step(imol), &
              !         nstep)
              nt = min(istep + nt_range * option%nt_sparse, &
                       state(ifile)%read_step(imol),        &
                       nstep)
           
              it = - 1

              if (option%tcf_mode == TcfModePji) then
                do jstep = istep, nt, option%nt_sparse 
                  it = it + 1
                  js = state(ifile)%data(jstep, imol)
                  hist(it, js, is) = hist(it, js, is) + 1.0d0
                  norm(it, is)     = norm(it, is)     + 1.0d0 
                end do
              else if (option%tcf_mode == TcfModeRji) then
                do jstep = istep, nt, option%nt_sparse 
                  it = it + 1
                  js = state(ifile)%data(jstep, imol)
                  
                  if (js /= is .and. js /= 0) then
                    hist(it:nt_range, js, is) = hist(it:nt_range, js, is) + 1.0d0
                    norm(it:nt_range, is)     = norm(it:nt_range, is)     + 1.0d0
                    exit
                  else
                    norm(it, is) = norm(it, is) + 1.0d0
                  end if

                end do
              end if
           
            end do
         
          end if

        end do

        state(ifile)%hist(:, :, :) = hist(:, :, :)
        state(ifile)%norm(:, :)    = norm(:, :)

        ! Trial implementation
        !
        !state(ifile)%hist(:, id_r, id_r) &
        !  = state(ifile)%hist(:, id_r, id_r) + state(ifile)%hist(:, id_b, id_r)
      end do
!$omp end do
!$omp end parallel

      ! Memory deallocation for 1st-order transition
      !
      deallocate(hist, norm)

      if (.not. option%calc_2nd) then
        return
      end if

      allocate(hist2(0:nt_range, 0:nstate))
      hist2 = 0.0d0

      write(iw,*)
      write(iw,'("Get_TransTCF> Calculate 2nd-order transition")')
!$omp parallel private(ifile, imol, istep, jstep, kstep, it, jt, is, js, ks, &
!$omp                  nt, nt2, hist2, hsum,                                 &
!$omp                  nthreads, thread_id, nwork, progress, ifinish),       &
!$omp          default(shared)
      nthreads  = omp_get_num_threads()
      thread_id = omp_get_thread_num()
      nwork     = nfile / nthreads + 1
      ifinish   = 0
!$omp do schedule(dynamic)
      do ifile = 1, nfile

        if (thread_id == 0) then
          ifinish  = ifinish + 1
          progress = ifinish / dble(nwork) * 100.0d0
          write(iw,'("  Progress : ", f6.2, "%")') progress
        end if

        hist2    = 0.0d0

        do imol = 1, nmol
          if (.not. state(ifile)%is_reacted(imol)) then
            do istep = 1, min(nstep, state(ifile)%read_step(imol))
           
              is = state(ifile)%data(istep, imol)
              hsum = 0.0d0
         
              if (is == id_r) then
                nt = min(istep + (nt_range - 1) * option%nt_sparse,  &
                         state(ifile)%read_step(imol),             &
                         nstep)
          
                hist2(0, is) = hist2(0, is) + option%dt_out 

                it = - 1
                do jstep = istep, nt, option%nt_sparse 
                  it = it + 1
                  js = state(ifile)%data(jstep, imol)
                  ks = state(ifile)%data(jstep + option%nt_sparse, imol)
        
                  hsum(js) = hsum(js) + option%dt_out
                  if (ks == id_r) then
                    hist2(it + 1, :) = hist2(it + 1, :) + hsum(:)
                  end if

                end do 

              end if
           
            end do
         
          end if
        end do ! imol

        state(ifile)%hist2(:, :) = hist2(:, :)

        ! Trial implementation
        !
        !state(ifile)%hist2(:, id_r) &
        !  = state(ifile)%hist2(:, id_r) + state(ifile)%hist2(:, id_b)
      end do
!$omp end do
!$omp end parallel

      write(iw,'("Finished!")')
      ! Memory deallocation
      !
      deallocate(hist2, hsum)

    end subroutine get_transtcf 
!-----------------------------------------------------------------------


end module mod_tcf
!=======================================================================
