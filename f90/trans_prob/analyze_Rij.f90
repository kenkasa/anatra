!-----------------------------------------------------------------------
    subroutine calc_Rij_wo_normalize(option, state, Rij)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      type(s_state),  intent(in)    :: state
      real(8),        intent(inout) :: Rij(0:option%nt_range, &
                                           option%nstate,     &
                                           option%nstate)

      ! Local
      !
      integer :: nmol, nstep, nt_range, nt_sparse

      ! Dummy 
      !
      integer :: istep, jstep, imol
      integer :: is, js
      integer :: it, it_reac, it_diff, nt
      real(8) :: dt


      ! Setup
      !
      nmol      = option%nmol
      nstep     = state%nstep
      nt_range  = option%nt_range
      nt_sparse = option%nt_sparse
      dt        = option%dt_out

      do imol = 1, nmol

        ! Search reaction time
        !
        js      = state%data(1, imol)
        it_reac = 0 
        do istep = nt_sparse + 1, nstep, nt_sparse
          is = state%data(istep, imol)
          if (is /= js) then
            it_reac = istep - nt_sparse 
            exit
          end if
        end do

        if (option%is_initial(js) .and. it_reac == 0) then
          write(iw,'("Calc_Rij_wo_normalize> Error.")')
          write(iw,'("No reaction is observed.")')
          stop
        end if

        ! Calc.
        !
        do istep = 0, it_reac, nt_sparse
          it_diff              = it_reac - istep 
          it_diff              = int(it_diff / dble(nt_sparse))
          Rij(it_diff, is, js) = Rij(it_diff, is, js) + 1.0d0
        end do

      end do

    end subroutine calc_Rij_wo_normalize
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine normalize_Rij(option, Rij)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      real(8),        intent(inout) :: Rij(0:option%nt_range, &
                                           option%nstate,     &
                                           option%nstate)


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
          Rij(:, :, js) = 0.0d0
          cycle
        end if 

        weight = option%state_weight(js)

        state_sum = 0.0d0
        do is = 1, nstate
          if (is /= js) then
            state_sum = state_sum + sum(Rij(:, is, js)) 
          end if
        end do

        do is = 1, nstate
          if (is /= js) then
            Rij(:, is, js) = weight * Rij(:, is, js) / (state_sum * dt) 
          end if
        end do
      end do 


    end subroutine normalize_Rij
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine running_integral_Rij(option, Rij, Rij_int)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      real(8),        intent(inout) :: Rij(0:option%nt_range, &
                                           option%nstate,     &
                                           option%nstate)
      real(8),        intent(inout) :: Rij_int(0:option%nt_range, &
                                               option%nstate,     &
                                               option%nstate)


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

      Rij_int  = 0.0d0
      do js = 1, nstate
        do is = 1, nstate
          if (is /= js) then
            prev = 0.0d0
            do istep = 0, nt_range - 1
              now                        = prev + Rij(istep, is, js) * dt
              prev                       = now
              Rij_int(istep + 1, is, js) = now 
            end do
          end if
        end do
      end do 


    end subroutine running_integral_Rij 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine write_Rij(output, option, boundary, Rij, Rij_int)
!-----------------------------------------------------------------------
      implicit none

      type(s_output),   intent(in) :: output
      type(s_option),   intent(in) :: option
      type(s_boundary), intent(in) :: boundary
      real(8),          intent(in) :: Rij(0:option%nt_range, &
                                          option%nstate,     &
                                          option%nstate)
      real(8),          intent(in) :: Rij_int(0:option%nt_range, &
                                              option%nstate,     &
                                              option%nstate)


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
                write(io,'(f20.10)', advance = 'no') Rij    (istep, is, js) 
                write(io,'(f20.10)', advance = 'no') Rij_int(istep, is, js) 
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
    subroutine Rij_bin(output, option, boundary, Rij)
!-----------------------------------------------------------------------
      implicit none

      type(s_output),   intent(in)    :: output 
      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(in)    :: boundary 
      real(8),          intent(inout) :: Rij(0:option%nt_range, &
                                             option%nstate,     &
                                             option%nstate)

      ! I/O
      !
      integer :: io
      logical :: io_status

      ! Local
      !
      integer                :: nt_range, nstate
      integer                :: nset 
      character(len=MaxChar) :: fname
      
      ! Dummy
      !
      integer :: i, j, k, ib, is, js, ifile, iset

      ! Arrays
      !
      real(8), allocatable :: wrk(:)


      ! Setup
      !
      nt_range = option%nt_range
      nstate   = option%nstate 

      allocate(wrk(0:nt_range))
  
      if (.not. option%write_Rij_bin &
    .and. .not. option%read_Rij_bin) return 

      ! Read
      !
      if (option%read_Rij_bin) then
        ifile = 0
        do while (.true.)
          write(fname,'(a,i4.4,".rbin")') trim(output%fhead), ifile + 1
          inquire(file=trim(fname), exist=io_status)
          if (io_status) then
            call open_file(fname, io, frmt = 'unformatted')
            read(io) nset
            do iset = 1, nset
              read(io) i, j, wrk
              Rij(:, i, j) = wrk(:) * option%state_weight(j)
              write(iw,'("Read Rij ",i0,2x,i0)') i, j
            end do
            close(io)
            ifile = ifile + 1 
          else
            return
          end if
        end do 
      end if

      ! Write 
      !
      if (option%write_Rij_bin) then

        nset = 0
        do j = 1, nstate
          if (.not. option%is_initial(j)) cycle
          do i = 1, nstate
            if (i == j .or. .not. boundary%is_connected(i, j)) cycle
              nset = nset + 1
          end do
        end do

        write(fname,'(a,i4.4,".rbin")') trim(output%fhead), option%out_id_rij
        write(iw,'("Write Rij")')
        call open_file(fname, io, frmt = 'unformatted')
        write(io) nset
        do j = 1, nstate
          if (.not. option%is_initial(j)) cycle
          do i = 1, nstate
            if (i == j .or. .not. boundary%is_connected(i, j)) cycle
            wrk(:) = Rij(:, i, j) / option%state_weight(j)
            write(iw,'(i0,2x,i0)') i, j
            write(io) i, j, wrk
          end do
        end do
        close(io)
      end if

      deallocate(wrk)

    end subroutine Rij_bin 
!-----------------------------------------------------------------------
