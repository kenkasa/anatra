!-----------------------------------------------------------------------
    subroutine reacdyn_tcf(output, option, boundary, f)
!-----------------------------------------------------------------------
      implicit none

      type(s_output),   intent(in)    :: output
      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(in)    :: boundary
      type(s_func),     intent(inout) :: f

      ! I/O
      !
      integer :: io, io_q, io_p

      ! Local
      !
      integer                :: nstate, nt_range, nt_extend, nboundary
      integer                :: nt_tcfout
      real(8)                :: dt
      real(8)                :: calc_time_sta, calc_time_end
      character(len=MaxChar) :: fname

      integer                :: thread_id

      ! Dummy
      !
      integer :: is, js, is1, is2, js1, js2
      integer :: ib, jb, inflx, jnflx, id, istep, jstep, jsta
      integer :: it, it_kfr, it_mfr, it_rel
      integer :: ncount, nt_life
      real(8) :: psum, pint
      real(8) :: pval, qval, rval

      ! Arrays
      !
      real(8), allocatable :: Qij(:, :)
      real(8), allocatable :: Pi(:, :)
      real(8), allocatable :: Pold(:), Pnew(:)

      real(8), allocatable :: Qint(:), Mfinal(:)

      real(8), allocatable :: kfr(:, :, :), mfr(:, :)
      integer, allocatable :: tind_kfr(:, :, :), tind_mfr(:, :)
      integer, allocatable :: tind_kfr_final(:, :), tind_mfr_final(:)


      if (.not. option%extrapolate) return

      write(iw,*)
      write(iw,'("Analyze> Start propagation")')

      ! Setup
      !
      nstate    = option%nstate
      nt_range  = option%nt_range
      nt_extend = option%nt_extend
      nt_tcfout = option%nt_tcfout
      dt        = option%dt_out
      nboundary = boundary%nboundary

      ! Allocate
      !
      allocate(Pold(nstate), Pnew(nstate))

      nt_life = 0
      do ib = -nboundary, nboundary

        if (ib == 0) cycle

        is1 = boundary%b2p(1, ib)
        is2 = boundary%b2p(2, ib)

        do js = 1, nstate
          do istep = 0, nt_range
            if (f%K(istep, js, ib) > 1.0d-10 .and. nt_life < istep) then
              nt_life = istep 
            end if 
          end do
        end do

        if (option%is_initial(is1)) then
          do istep = 0, nt_range
            if (f%R(istep, is2, is1) > 1.0d-10 .and. nt_life < istep) then
              nt_life = istep
            end if
          end do
        end if
        
      end do

      ! Add margin
      !
      nt_life = nt_life + 10
      allocate(Qij            (0:nt_life, -nboundary:nboundary))
      allocate(Pi             (0:nt_life, nstate))
      allocate(Qint           (-nboundary:nboundary))
      allocate(Mfinal         (-nboundary:nboundary))

      allocate(kfr            (0:nt_range, nstate, -nboundary:nboundary))
      allocate(tind_kfr       (0:nt_range, nstate, -nboundary:nboundary))
      allocate(tind_kfr_final (nstate, -nboundary:nboundary))
      allocate(mfr            (0:nt_range, -nboundary:nboundary))
      allocate(tind_mfr       (0:nt_range, -nboundary:nboundary))
      allocate(tind_mfr_final (-nboundary:nboundary))

      kfr            = 0.0d0
      tind_kfr       = 0
      tind_kfr_final = 0
      mfr            = 0.0d0
      tind_mfr       = 0
      tind_mfr_final = 0

      ! Construct kfr
      !
      write(iw,'("ReacDyn_TCF> Construct K-matrix")')
      do ib = -nboundary, nboundary
        if (ib == 0) cycle
        is1 = boundary%b2p(1, ib)
        is2 = boundary%b2p(2, ib)
        do inflx = 1, boundary%n_influx_boundary(is1)
          jb    = boundary%influx_boundary(inflx, is1)
          js    = boundary%b2p(1, jb)
          jstep = 0
          do istep = 0, nt_life
            if (f%K(istep, is2, jb) > 1.0d-10) then
              jstep                    = jstep + 1
              tind_kfr(jstep, is2, jb) = istep
              kfr     (jstep, is2, jb) = f%K(istep, is2, jb)
            end if
          end do
          tind_kfr_final(is2, jb) = jstep
          write(iw,'(4i8)') is2, is1, js, jstep 
        end do 
      end do

      ! Construct mfr 
      !
      do ib = -nboundary, nboundary
        if (ib == 0) cycle
        jstep = 0
        do istep = 0, nt_life
          if (abs(f%M(istep, ib)) > 1.0d-10) then
            jstep               = jstep + 1
            tind_mfr(jstep, ib) = istep
            mfr     (jstep, ib) = f%M(istep, ib)
          end if
        end do
        tind_mfr_final(ib) = jstep
      end do

      ! Propagation
      !
      write(fname,'(a,".tcf")') trim(output%fhead)
      call open_file(fname, io)

      write(fname,'(a,".q")')   trim(output%fhead)
      call open_file(fname, io_q)

      write(fname,'(a,".pj")')  trim(output%fhead)
      call open_file(fname, io_p)

      Qij       = 0.0d0
      Pi        = 0.0d0
      pint      = 0.0d0
      Qint      = 0.0d0
      Mfinal(:) = f%M(nt_range, :) 

      calc_time_sta = omp_get_wtime()

      !!$omp parallel private(istep, jstep, it, ib, jb, is, is1, is2, js1, js2, &
      !!$omp                  inflx, jsta, it_kfr, it_mfr, it_rel, qval, pval, rval, psum)  &
      !!$omp          default(shared) 

      it = -1
      do istep = 0, nt_extend
        it = it + 1

        !!$omp single
        if (mod(istep, 1000) == 0) then
        !if (mod(istep, 1) == 0) then
          write(iw,'("Step: ", i20, " Time: ", f20.10)') istep, istep * dt
        end if
        !!$omp end single 

        ! Calc. Q
        !
        !!$omp do schedule(dynamic) 
        do ib = -nboundary, nboundary

          qval = 0.0d0
          if (ib == 0) cycle 

          if (option%use_reflection_state) then
            if (boundary%conv_direc(-ib)) cycle
          end if

          if (option%use_constant_Qij .and. boundary%is_cQij(ib)) then
            Qij(it, ib) = boundary%cQij(ib)
            !if (istep == 0) then
            !  write(iw,'("Constant-Qij ", i5,2x,e15.7)') ib, Qij(ib, it)
            !end if
            cycle
          end if

          is1  = boundary%b2p(1, ib)
          is2  = boundary%b2p(2, ib)

          if (istep <= nt_range .and. option%is_initial(is1)) then
            qval = qval + f%R(it, is2, is1)
          end if

          if (istep == 0) then
            Qij(it, ib) = qval
            if (option%use_reflection_state) then
              if (boundary%conv_direc(ib)) then
                Qij(it, -ib) = Qij(it, ib)
                Qij(it,  ib) = 0.0d0
              end if
            end if
            cycle
          end if

          !if (istep == 0) cycle
          do inflx = 1, boundary%n_influx_boundary(is1)
            jb  = boundary%influx_boundary(inflx, is1)
            js1 = boundary%b2p(1, jb)
            js2 = boundary%b2p(2, jb)

            jsta = max(0, istep - nt_life)
            rval = 0.0d0
            do it_kfr = 1, tind_kfr_final(is2, jb)
              jstep  = tind_kfr(it_kfr, is2, jb)
              it_rel = istep - jstep - jsta
              if (it_rel < 0) exit 
              rval = rval + dt * kfr(it_kfr, is2, jb) * Qij(it_rel, jb) 
            end do
            qval = qval + rval

          end do
          Qij(it, ib) = qval
          if (option%use_reflection_state) then
            if (boundary%conv_direc(ib)) then
              Qij(it, -ib) = Qij(it, ib)
              Qij(it,  ib) = 0.0d0
            end if
          end if

        end do  ! ib
        !!$omp end do
        !!$omp barrier

        ! Calc. P
        !
        do is = 1, nstate

          ! Initial
          !
          Pnew(is) = 0.0d0
          if (istep == 0) then
            Pold(is) = option%state_weight(is)
          end if

          ! Update
          !
          Pnew(is) = Pold(is)

          do inflx = 1, boundary%n_influx_boundary(is)
            ib       = boundary%influx_boundary(inflx, is)
            Pnew(is) = Pnew(is) + (Qij(it, ib) - Qij(it, -ib)) * dt
          end do

          Pi(it, is) = Pold(is)

        end do

!        !!$omp do schedule(dynamic)
!        do is = 1, nstate
!
!          pval = 0.0d0
!          if (istep <= nt_range .and. option%is_initial(is)) then
!            pval = pval + f%P0(it, is)
!          end if
!
!          if (istep == 0) then
!            Pi(it, is) = pval
!            cycle
!          end if
!
!          do inflx = 1, boundary%n_influx_boundary(is)
!            ib = boundary%influx_boundary(inflx, is)
!            jsta = max(0, istep - nt_life)
!            rval = 0.0d0
!            do it_mfr = 2, tind_mfr_final(ib)
!              jstep  = tind_mfr(it_mfr, ib)
!              it_rel = istep - jstep - jsta
!              if (it_rel < 0) cycle 
!              rval = rval + dt * mfr(it_mfr, ib) * Qij(it_rel, ib)
!            end do
!
!            pval = pval + rval
!            pval = pval + Mfinal(ib) * Qint(ib)
!          end do
!          Pi(it, is) = pval
!        end do
!        !!$omp end do
!        !!$omp barrier

        ! Output
        !
        !!$omp single
        psum = 0.0d0
        do is = 1, nstate
          if (option%is_initial(is)) then
            !psum = psum + Pi(it, is)
            psum = psum + Pold(is)
          end if
        end do
        pint = pint + psum * dt


        if (istep == 0) then
          jb = 0
          do ib = - nboundary, nboundary
            if (ib == 0) cycle
            jb  = jb + 1
            is1 = boundary%b2p(1, ib)
            is2 = boundary%b2p(2, ib) 
            write(io_q, '("# Col. ", i0, " : ", i0, " <--- ", i0)') jb + 1, is2, is1 
          end do
        end if

        if (mod(istep, nt_tcfout) == 0 .or. istep == 0) then
         
          write(io_q,'(e15.7,2x)', advance = 'no') dt * istep
          write(io_p,'(e15.7,2x)', advance = 'no') dt * istep
          write(io,  '(3(e15.7,2x))') dt * istep, truncate_tiny_val(psum), pint

          do ib = -nboundary, nboundary
            if (ib == 0) cycle 
            write(io_q,'(e15.7,2x)', advance = 'no') Qij(it, ib)
          end do
          write(io_q,*)

          do is = 1, nstate
            write(io_p,'(e15.7,2x)', advance = 'no') truncate_tiny_val(Pi(it, is))
          end do
          write(io_p,*)

        end if

        ! Shift time origin
        !
        if (it == nt_life) then

          Qint(:) = Qint(:) + Qij(0, :) * dt

          do ib = -nboundary, nboundary
            if (ib == 0) cycle
            Qij(0:nt_life - 1, ib) = Qij(1:nt_life, ib) 
          end do

          do is = 1, nstate
            Pi(0:nt_life - 1, is) = Pi(1:nt_life, is)
          end do

          Qij(nt_life, :) = 0.0d0
          Pi(nt_life, :)  = 0.0d0

          !it = it - 1

        end if

        Pold(:) = Pnew(:)

        !!$omp end single

        if (it == nt_life) then
          it = it - 1
        end if 

      end do    ! istep
      !!$omp end parallel

      calc_time_end = omp_get_wtime()
      write(iw,'("Calculation time (sec): ", f15.7)') calc_time_end - calc_time_sta

      close(io)
      close(io_q)
      close(io_p) 

      ! Deallocate memory
      !
      deallocate(Pi, Qij)
      deallocate(Qint, Mfinal)
      
    end subroutine reacdyn_tcf 
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    function truncate_tiny_val(val)
!-----------------------------------------------------------------------
      implicit none

      ! Constants
      !
      real(8) :: eps_tr = 1.0d-10

      real(8) :: val
      real(8) :: truncate_tiny_val

      
      if (abs(val) < eps_tr) then
        truncate_tiny_val = 0.0d0 
      else
        truncate_tiny_val = val
      end if

    end function truncate_tiny_val
!-----------------------------------------------------------------------
