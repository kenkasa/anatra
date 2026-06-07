!-----------------------------------------------------------------------
    subroutine update_rsto_Rij(option, state, f)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      type(s_state),  intent(in)    :: state
      type(s_func),   intent(inout) :: f

      ! Local
      !
      integer :: nmol, nstep, nt_range, nt_sparse
      integer :: nini
      integer :: use_for_Rij

      ! Dummy 
      !
      integer :: istep, jstep, ito, imol, iend
      integer :: is, js, it, it_reac, it_diff, nt
      integer :: unp_id
      real(8) :: dt

      ! Arrays
      !
      integer, allocatable :: init_set(:)
      integer, allocatable :: rand(:)
      

      !if (.not. option%use_perturbed_traj) then
      !  write(iw,'("Update_Rsto_Rij> Error.")')
      !  write(iw,'("use_perturbed_traj should be .true.")')
      !  stop
      !end if

      if (.not. (state%use_for_Rij == -1 .or. state%use_for_Rij == 1)) return

      ! Setup
      !
      nmol        = option%nmol
      nstep       = state%nstep
      unp_id      = state%unperturbed_id
      use_for_Rij = state%use_for_Rij
      nt_range    = option%nt_range
      nt_sparse   = option%nt_sparse
      dt          = option%dt_out

      ! Allocate
      !
      if (.not. allocated(init_set)) then
        allocate(init_set(nstep))
        allocate(rand(nstep))
      end if

      do imol = 1, nmol

        nini     =  0
        init_set = -1

        ! Check state id at each time step
        !
        do jstep = 1, nstep
          js = state%data(jstep, imol)
          if ((js == unp_id .or. unp_id == -1) .and. option%is_initial(js)) then
            nini           = nini + 1
            init_set(nini) = jstep 
          end if
        end do

        if (nini == 0) return

        ! Generate random numbers
        !
        call get_random_integer(nini, 1, nini, .true., rand(1:nini)) 

        ! Generate histogram
        !
        do ito = 1, nini * 0.5
          jstep = init_set(rand(ito))
          if (jstep + nt_sparse > nstep) cycle

          ! Search reaction time
          !
          js      = state%data(jstep, imol)
          it_reac = -1

          if (js /= unp_id .and. unp_id /= -1) then
            write(iw,'("Update_rsto_Rij> Error. Something went wrong")')
            stop
          end if 

          do istep = nt_sparse + jstep, nstep, nt_sparse
            is = state%data(istep, imol)
            if (is /= js) then
              iend    = istep - nt_sparse
              it_reac = iend  - jstep
              exit 
            end if
          end do

          !if (it_reac == -1) then
          !  write(iw,'("Calc_Rij_wo_normalize> ")')
          !  write(iw,'("Molecule ", i5, ": No reaction is observed.")') imol
          !  cycle 
          !end if

          do istep = 0, it_reac, nt_sparse
            it_diff              = it_reac - istep
            it_diff              = int(it_diff / dble(nt_sparse))
            f%R(it_diff, is, js) = f%R(it_diff, is, js) + 1.0d0 
          end do

        end do
      end do

    end subroutine update_rsto_Rij
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine update_Rij_wo_normalize(option, state, f)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      type(s_state),  intent(in)    :: state
      type(s_func),   intent(inout) :: f

      ! Local
      !
      integer :: nmol, nstep, nt_range, nt_sparse
      integer :: use_for_Rij

      ! Dummy 
      !
      integer :: istep, jstep, imol
      integer :: is, js, it, it_reac, it_diff, nt
      integer :: unp_id
      real(8) :: dt


      ! Setup
      !
      nmol        = option%nmol
      nstep       = state%nstep
      unp_id      = state%unperturbed_id
      use_for_Rij = state%use_for_Rij
      nt_range    = option%nt_range
      nt_sparse   = option%nt_sparse
      dt          = option%dt_out

      do imol = 1, nmol

        ! Search reaction time
        !
        js = state%data(1, imol)

        if (.not. option%is_initial(js)) cycle

        if (.not. (use_for_Rij == -1 .or. use_for_Rij == 1)) return

        it_reac = 0 
        do istep = nt_sparse + 1, nstep, nt_sparse
          is = state%data(istep, imol)
          if (is /= js) then
            it_reac = istep - nt_sparse 
            exit
          end if
        end do

        if (it_reac == 0) then
          write(iw,'("Calc_Rij_wo_normalize> ")')
          write(iw,'("Molecule ", i5, ": No reaction is observed.")') imol
          cycle 
        end if

        ! Calc.
        !
        do istep = 0, it_reac, nt_sparse
          it_diff              = it_reac - istep 
          it_diff              = int(it_diff / dble(nt_sparse))
          f%R(it_diff, is, js) = f%R(it_diff, is, js) + 1.0d0
        end do

      end do

    end subroutine update_Rij_wo_normalize
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine update_Rij_from_hist(io, option, f)
!-----------------------------------------------------------------------
      implicit none

      integer,        intent(in)    :: io
      type(s_option), intent(in)    :: option
      type(s_func),   intent(inout) :: f

      ! Local
      !

      ! Dummy 
      !
      character(len=MaxChar) :: line
      integer                :: istep
      integer                :: is, js
      real(8)                :: val

      ! Setup
      !

      do while (.true.)
       read(io,'(a)',end=100) line
       if (line(1:3) == 'END') go to 100
       
       read(line,*) is, js, istep, val
       f%R(istep, is, js) = f%R(istep, is, js) + val
      end do
 100  return 

    end subroutine update_Rij_from_hist
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine normalize_Rij(option, f)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      type(s_func),   intent(inout) :: f

      ! Local
      !
      integer :: nt_range, nstate
      real(8) :: dt

      ! Dummy
      !
      integer :: istep, is, js
      real(8) :: state_sum, weight


      ! Setup
      !
      nstate   = option%nstate
      nt_range = option%nt_range
      nstate   = option%nstate
      dt       = option%dt_out

      do js = 1, nstate

        if (.not. option%is_initial(js)) then
          f%R(:, :, js) = 0.0d0
          cycle
        end if 

        weight = option%state_weight(js)

        state_sum = 0.0d0
        do is = 1, nstate
          if (is /= js) then
            state_sum = state_sum + sum(f%R(:, is, js)) 
          end if
        end do

        do is = 1, nstate
          if (is /= js) then
            f%R(:, is, js) = weight * f%R(:, is, js) / (state_sum * dt) 
          end if
        end do
      end do 


    end subroutine normalize_Rij
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine running_integral_Rij(option, f)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      type(s_func),   intent(inout) :: f

      ! Local
      !
      integer :: nt_range, nstate
      real(8) :: dt

      ! Dummy
      !
      integer :: istep, is, js
      real(8) :: prev, now 


      ! Setup
      !
      nstate   = option%nstate
      nt_range = option%nt_range
      nstate   = option%nstate
      dt       = option%dt_out

      f%Rint  = 0.0d0
      do js = 1, nstate
        do is = 1, nstate
          if (is /= js) then
            prev = 0.0d0
            do istep = 0, nt_range - 1
              now                       = prev + f%R(istep, is, js) * dt
              prev                      = now
              f%Rint(istep + 1, is, js) = now 
            end do
          end if
        end do
      end do 


    end subroutine running_integral_Rij 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine write_Rij(output, option, boundary, f)
!-----------------------------------------------------------------------
      implicit none

      type(s_output),   intent(in) :: output
      type(s_option),   intent(in) :: option
      type(s_boundary), intent(in) :: boundary
      type(s_func),     intent(in) :: f

      ! I/O
      !
      integer :: io

      ! Local
      !
      integer                :: nstate, nt_range
      character(len=MaxChar) :: fname

      ! Dummy
      !
      integer :: is, js, id, istep


      ! Setup
      !
      nstate   = option%nstate
      nt_range = option%nt_range

      write(fname, '(a,".Rij")') trim(output%fhead)
      call open_file(fname, io)
      
      write(io,'("# J <--- I")')

      id = 1 
      do js = 1, nstate
        if (option%is_initial(js)) then
          do is = 1, nstate
            if (boundary%is_connected(is, js)) then
              write(io,'("# Col. ", i0, " Rij : ", i0, " - ", i0)') id + 1, is, js
              write(io,'("# Col. ", i0, " Int : ", i0, " - ", i0)') id + 2, is, js
              id = id + 2 
            end if
          end do
        end if
      end do

      do istep = 0, nt_range
        write(io,'(f20.10)', advance = 'no') option%dt_out * istep 
        do js = 1, nstate
          if (option%is_initial(js)) then
            do is = 1, nstate
              if (boundary%is_connected(is, js)) then
                write(io,'(e15.7,2x)', advance = 'no') f%R   (istep, is, js) 
                write(io,'(e15.7,2x)', advance = 'no') f%Rint(istep, is, js) 
              end if
            end do
          end if
        end do
        write(io,*)
      end do

      close(io) 

    end subroutine write_Rij
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine output_Rij_hist(option, output, f, is_gen)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      type(s_output), intent(in)    :: output
      type(s_func),   intent(inout) :: f
      logical,        intent(inout) :: is_gen 

      ! I/O
      !
      integer :: io

      ! Local
      !
      character(len=MaxChar) :: fname
      integer                :: nt_range, nstate

      ! Dummy
      !
      integer :: is, js, istep
      real(8) :: val


      ! Setup
      !
      nstate   = option%nstate
      nt_range = option%nt_range

      val = sum(f%R(:, :, :))
      if (val < 0.999d0) then
        return
      end if

      write(fname,'(a,".rhist")') trim(output%fhead)
      call open_file(fname, io)

      write(io,'("RIJ")')

      do js = 1, nstate
      do is = 1, nstate
      do istep = 0, nt_range
        val = f%R(istep, is, js)
        if (val >= 0.999d0) then
          write(io, '(i5, 2x, i5, 2x, i10, 2x, f20.10)') is, js, istep, val
        end if
      end do 
      end do
      end do 
      write(io, '("END")')

      close(io)


      write(fname, '(a,".krhist")') trim(output%fhead)
      if (is_gen) then
        call open_file(fname, io, pos = 'APPEND') 
      else
        call open_file(fname, io) 
      end if
      is_gen = .true.

      write(io,'("RIJ")')

      do js = 1, nstate
      do is = 1, nstate
      do istep = 0, nt_range
        val = f%R(istep, is, js)
        if (val >= 0.999d0) then
          write(io, '(i5, 2x, i5, 2x, i10, 2x, f20.10)') is, js, istep, val
        end if
      end do 
      end do
      end do 
      write(io, '("END")')

      close(io)

    end subroutine output_Rij_hist
!-----------------------------------------------------------------------
