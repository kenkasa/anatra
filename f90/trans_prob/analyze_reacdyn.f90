!-----------------------------------------------------------------------
    subroutine reacdyn_tcf(output, option, boundary, Rij, P0, Kijk, Mij)
!-----------------------------------------------------------------------
      implicit none

      type(s_output),   intent(in) :: output
      type(s_option),   intent(in) :: option
      type(s_boundary), intent(in) :: boundary
      real(8),          intent(in) :: Rij(0:option%nt_range, &
                                          option%nstate,     &
                                          option%nstate)
      real(8),          intent(in) :: P0(0:option%nt_range,  &
                                         option%nstate)
      real(8),          intent(in) :: Kijk(0:option%nt_range, &
                                          option%nstate,      &
                                          -boundary%nboundary:boundary%nboundary)
      real(8),          intent(in) :: Mij(0:option%nt_range,  &
                                         -boundary%nboundary:boundary%nboundary)

      ! I/O
      !
      integer :: io, io_q, io_p

      ! Local
      !
      integer                :: nstate, nt_range, nt_extend, nboundary
      integer                :: nt_tcfout
      real(8)                :: dt
      character(len=MaxChar) :: fname

      ! Dummy
      !
      integer :: is, js, is1, is2, ib, jb, id, istep, jstep, jsta
      integer :: it, it_rel
      integer :: js1, js2
      integer :: ncount, nt_life
      real(8) :: psum, pint

      ! Arrays
      !
      real(8), allocatable :: Qij(:, :)
      real(8), allocatable :: Pi(:, :)


      ! Setup
      !
      nstate    = option%nstate
      nt_range  = option%nt_range
      nt_extend = option%nt_extend
      nt_tcfout = option%nt_tcfout
      dt        = option%dt_out
      nboundary = boundary%nboundary

      nt_life = 0
      do ib = -nboundary, nboundary

        if (ib == 0) cycle

        is1 = boundary%b2p(1, ib)
        is2 = boundary%b2p(2, ib)

        do js = 1, nstate
          do istep = 0, nt_range
            if (Kijk(istep, js, ib) > 1.0d-10 .and. nt_life < istep) then
              nt_life = istep 
            end if 
          end do
        end do

        if (option%is_initial(is1)) then
          do istep = 0, nt_range
            if (Rij(istep, is2, is1) > 1.0d-10 .and. nt_life < istep) then
              nt_life = istep
            end if
          end do
        end if
        
      end do

      ! Add margin
      !
      nt_life = nt_life + 10
      allocate(Qij(0:nt_life, -nboundary:nboundary))
      allocate(Pi(0:nt_life, nstate))

      ! Propagation
      !
      write(fname,'(a,".tcf")') trim(output%fhead)
      call open_file(fname, io)

      write(fname,'(a,".q")') trim(output%fhead)
      call open_file(fname, io_q)

      write(fname,'(a,".pj")') trim(output%fhead)
      call open_file(fname, io_p)

      Qij   = 0.0d0
      Pi    = 0.0d0
      pint  = 0.0d0

      it = -1 
      do istep = 0, nt_extend
        it = it + 1

        if (mod(istep, 1000) == 0) then
          write(iw,'("Step: ", i0, "Time: ", f20.10)') istep, istep * dt
        end if 

        ! Calc. Q
        !
        do ib = -nboundary, nboundary

          if (ib == 0) cycle 

          is1 = boundary%b2p(1, ib)
          is2 = boundary%b2p(2, ib)

          if (istep <= nt_range .and. option%is_initial(is1)) then
            Qij(it, ib) = Qij(it, ib) + Rij(it, is2, is1)
          end if

          if (istep == 0) cycle 

          do jb = -nboundary, nboundary

            if (jb == 0) cycle

            js1 = boundary%b2p(1, jb)
            js2 = boundary%b2p(2, jb)

            if (js2 /= is1) cycle

            jsta   = istep - nt_life
            it_rel = - 1 
            do jstep = max(0, jsta), istep - 1
              it_rel     = it_rel + 1
              Qij(it, ib) = Qij(it, ib)  &
                + dt * Kijk(istep - jstep, is2, jb) * Qij(it_rel, jb)
            end do 

          end do

          if (option%use_reflection_state) then
            if (boundary%conv_direc(ib)) then
              Qij(it, -ib) = Qij(it, ib)
            end if
          end if

        end do  ! ib

        ! Calc. P
        !
        do is = 1, nstate

          if (istep <= nt_range .and. option%is_initial(is)) then
            Pi(it, is) = Pi(it, is) + P0(it, is)
          end if

          if (istep == 0) cycle

          do js = 1, nstate

            ib    = boundary%p2b(is, js)
            if (ib == 0) cycle

            jsta   = istep - nt_life !- 1
            it_rel = - 1
            do jstep = max(0, jsta), istep - 1
              it_rel = it_rel + 1
              Pi(it, is) = Pi(it, is)  &
                + dt * Mij(istep - jstep, ib) * Qij(it_rel, ib)
            end do

          end do 

        end do

        ! Output
        !
        psum = 0.0d0
        do is = 1, nstate
          if (option%is_initial(is)) then
            psum = psum + Pi(it, is)
          end if
        end do
        pint = pint + psum * dt

        if (mod(istep, nt_tcfout) == 0 .or. istep == 0) then
         
          write(io_q,'(f20.10)', advance = 'no') dt * istep
          write(io_p,'(f20.10)', advance = 'no') dt * istep
          write(io, '(3f20.10)') dt * istep, psum, pint

          do ib = -nboundary, nboundary
            if (ib == 0) cycle 
            write(io_q,'(f20.10)', advance = 'no') Qij(it, ib)
          end do
          write(io_q,*)

          do is = 1, nstate
            write(io_p,'(f20.10)', advance = 'no') Pi(it, is)
          end do
          write(io_p,*)

        end if

        ! Shift time origin
        !
        if (it == nt_life) then

          do jstep = 1, nt_life
            do ib = -nboundary, nboundary
              if (ib == 0) cycle
              Qij(jstep - 1, ib) = Qij(jstep, ib) 
            end do

            do is = 1, nstate
              Pi(jstep - 1, is) = Pi(jstep, is)
            end do
          end do

          Qij(nt_life, :) = 0.0d0
          Pi(nt_life, :) = 0.0d0

          it = it - 1

        end if 

      end do    ! istep

      close(io)
      close(io_q)
      close(io_p) 

!
!      ! Propagating P
!      !
!      Pi = 0.0d0
!      do is = 1, nstate
!        !if (.not. option%is_initial(is)) cycle
!
!        do istep = 0, nt_extend
!
!          if (option%is_initial(is)) then 
!            Pi(istep, is) = Pi(istep, is) + ZP(nt_range, istep, P0(:, is))
!          end if
!
!          if (istep == 0) cycle
!
!          do js = 1, nstate
!            ib = boundary%p2b(is, js)
!
!            if (ib == 0) cycle
!           
!            jsta = istep - nt_life - 1 
!            !do jstep = 0, istep - 1
!            do jstep = max(0, jsta), istep - 1
!              Pi(istep, is) = Pi(istep, is)                   &
!                + dt * ZP(nt_range, istep - jstep, Mij(:, ib)) &
!                     * Qb(jstep, ib)
!            end do 
!
!          end do
!            
!        end do
!      end do
!
!      !open(99, file = 'out.q')
!      !do istep = 0, nt_extend
!      !  write(99,'(f20.10)', advance = 'no') dt * istep
!
!      !  do ib = -nboundary, nboundary
!      !    if (ib == 0) cycle 
!      !    write(99,'(f20.10)', advance = 'no') Qb(istep, ib)
!      !  end do
!      !  write(99,*)
!      !end do
!      !close(99) 
!
!      open(99, file = 'out.tcf')
!
!      pint = 0.0d0
!      do istep = 0, nt_extend
!        !write(99,'(f20.10)', advance = 'no') dt * istep
!
!        psum = 0.d0
!        do is = 1, nstate
!          if (option%is_initial(is)) then
!            psum = psum + Pi(istep, is)
!          end if
!        end do
!        pint = pint + psum * dt
!
!        write(99,'(3f20.10)') dt * istep, psum, pint
!        !do is = 1, nstate
!        !  !if (option%is_initial(is)) then
!        !    write(99,'(f20.10)', advance = 'no') Pi(istep, is)
!        !  !end if
!        !end do
!        !write(99,*)
!      end do
!      close(99) 

      ! Deallocate memory
      !
      deallocate(Pi, Qij)
      
    end subroutine reacdyn_tcf 
!-----------------------------------------------------------------------
