!-----------------------------------------------------------------------
    subroutine calc_Kijk_wo_normalize(option, boundary, state, Kijk, hit_count)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(in)    :: boundary 
      type(s_state),    intent(in)    :: state
      real(8),          intent(inout) :: Kijk(0:option%nt_range,        &
                                             option%nstate,             &
                                             -boundary%nboundary:       &
                                              boundary%nboundary)
      real(8),          intent(inout) :: hit_count(-boundary%nboundary: &
                                                    boundary%nboundary) 

      ! Local
      !
      integer :: nmol, nstep, nt_range, nt_sparse, nboundary
      real(8) :: dt
      logical :: is_final

      ! Dummy 
      !
      integer :: istep, jstep, imol
      integer :: is, js, ib
      integer :: it, it_reac, it_diff, nt
      integer :: ista, iend, ireac
      integer :: iprod, idissoc


      ! Setup
      !
      nmol      = option%nmol
      nstep     = state%nstep
      nt_range  = option%nt_range
      nt_sparse = option%nt_sparse
      nboundary = boundary%nboundary
      dt        = option%dt_out

      do imol = 1, nmol

        ! Search reaction time
        !
        js       = state%data(1, imol)
        ireac    = 0
        is_final = .false. 
        do istep = nt_sparse + 1, nstep, nt_sparse
          is = state%data(istep, imol)
          if (is /= js) then
            ireac  = ireac + 1
            
            if (ireac == 1) then
              ista          = istep - nt_sparse
              ib            = boundary%p2b(is, js)
              js            = is
              hit_count(ib) = hit_count(ib) + 1.0d0
            else
              iend                  = istep - nt_sparse 
              it_diff               = iend - ista
              it_diff               = int(it_diff / dble(nt_sparse))
              Kijk(it_diff, is, ib) = Kijk(it_diff, is, ib) + 1.0d0

              ib                    = boundary%p2b(is, js)
              hit_count(ib)         = hit_count(ib) + 1.0d0
                                    
              ista                  = istep - nt_sparse 
              js                    = is

              if (istep == nstep) then
                is_final = .true.
              end if

            end if

          end if
        end do

        if (ireac == 0 .or. ireac == 1) then
          if (option%use_product_state) then
            do iprod = 1, option%nproduct
              if (js == option%product_state_ids(iprod)) then
                return
              end if
            end do 
          else if (option%use_dissociate_state) then
            do idissoc = 1, option%ndissoc
              if (js == option%dissociate_state_ids(idissoc)) then
                return
              end if
            end do 
          else
            write(iw,'("Calc_Kijk_wo_normalize> Error.")')
            write(iw,'("No reaction is observed.")')
            stop
          end if

          if (ireac == 1) then
            if (.not. option%is_dissoc(js)) then 
              hit_count(ib) = hit_count(ib) - 1.0d0
            end if
          end if

        else

          if (.not. option%is_dissoc(js)) then
            hit_count(ib) = hit_count(ib) - 1.0d0
          end if

        end if

      end do

    end subroutine calc_Kijk_wo_normalize
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine normalize_Kijk(option, boundary, Kijk, hit_count)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(in)    :: boundary 
      real(8),          intent(inout) :: Kijk(0:option%nt_range,        &
                                             option%nstate,             &
                                             -boundary%nboundary:       &
                                              boundary%nboundary)
      real(8),          intent(inout) :: hit_count(-boundary%nboundary: &
                                                    boundary%nboundary)
                                           


      ! Local
      !
      integer :: nt_range, nstate, nboundary
      real(8) :: dt

      ! Dummy
      !
      integer :: istep, is, is1, is2, js, ib
      real(8) :: state_sum


      ! Setup
      !
      nstate    = option%nstate
      nt_range  = option%nt_range
      nstate    = option%nstate
      nboundary = boundary%nboundary 
      dt        = option%dt_out

      do ib = -nboundary, nboundary

        if (ib == 0) then
          cycle
        end if

        is1 = boundary%b2p(1, ib)
        is2 = boundary%b2p(2, ib)

        state_sum      = hit_count(ib)
        Kijk(:, :, ib) = Kijk(:, :, ib) / (state_sum * dt) 

      end do

    end subroutine normalize_Kijk
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine renormalize_Kijk(option, boundary, Kijk)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(in)    :: boundary 
      real(8),          intent(inout) :: Kijk(0:option%nt_range,        &
                                             option%nstate,             &
                                             -boundary%nboundary:       &
                                              boundary%nboundary)


      ! Local
      !
      integer :: nt_range, nstate, nboundary
      real(8) :: dt

      ! Dummy
      !
      integer :: istep, is, is1, is2, js, ib
      real(8) :: state_sum


      ! Setup
      !
      nstate    = option%nstate
      nt_range  = option%nt_range
      nstate    = option%nstate
      nboundary = boundary%nboundary 
      dt        = option%dt_out

      do ib = -nboundary, nboundary

        if (ib == 0) then
          cycle
        end if

        is1 = boundary%b2p(1, ib)
        is2 = boundary%b2p(2, ib)

        if (.not. option%is_dissoc(is2)) then 
          state_sum = 0.0d0
          do js = 1, nstate
            if (.not. boundary%is_connected(js, is2)) cycle 
            do istep = 0, nt_range
              state_sum = state_sum + Kijk(istep, js, ib)
            end do
          end do
          Kijk(:, :, ib) = Kijk(:, :, ib) / (state_sum * dt)
        end if

      end do

    end subroutine renormalize_Kijk
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine write_Kijk(output, option, boundary, Kijk)
!-----------------------------------------------------------------------
      implicit none

      type(s_output),   intent(in) :: output
      type(s_option),   intent(in) :: option
      type(s_boundary), intent(in) :: boundary
      real(8),          intent(in) :: Kijk(0:option%nt_range,  &
                                           option%nstate,      &
                                          -boundary%nboundary:boundary%nboundary)

      ! I/O
      !
      integer :: io

      ! Local
      !
      integer                :: nstate, nt_range, nboundary
      character(len=MaxChar) :: fname

      ! Dummy
      !
      integer :: is, js, is1, is2, ib, id, istep


      ! Setup
      !
      nstate    = option%nstate
      nt_range  = option%nt_range
      nboundary = boundary%nboundary


      write(fname, '(a,".Kijk")') trim(output%fhead)
      call open_file(fname, io)
      
      write(io,'("# J <--- B")')
      id = 1 
      do ib = -nboundary, nboundary

        if (ib == 0) then
          cycle
        end if

        is1 = boundary%b2p(1, ib)
        is2 = boundary%b2p(2, ib)

        do js = 1, nstate
          if (boundary%is_connected(js, is2)) then
            write(io,'(" # Col. ", i0, " Kijk : ", i0, " <= ", i0, " < ", i0)')  &
              id + 1, js, is2, is1

            id = id + 1
          end if

        end do

      end do

      do istep = 0, nt_range
        write(io,'(f20.10)', advance = 'no') option%dt_out * istep 
        do ib = -nboundary, nboundary

          if (ib == 0) then
            cycle
          end if

          is1 = boundary%b2p(1, ib)
          is2 = boundary%b2p(2, ib)

          do js = 1, nstate
            if (boundary%is_connected(js, is2)) then
              write(io,'(f20.10)', advance = 'no') Kijk(istep, js, ib)
            end if 
          end do
        end do
        write(io,*)
      end do
      close(io)

    end subroutine write_Kijk
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine check_Kijk(option, boundary, Kijk)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(in)    :: boundary 
      real(8),          intent(inout) :: Kijk(0:option%nt_range,        &
                                             option%nstate,             &
                                             -boundary%nboundary:       &
                                              boundary%nboundary)
                                           


      ! Local
      !
      integer :: nt_range, nstate, nboundary
      real(8) :: dt

      ! Dummy
      !
      integer :: istep, is, is1, is2, js, ib
      real(8) :: state_sum


      if (.not. option%check_Kijk) return

      ! Setup
      !
      nstate    = option%nstate
      nt_range  = option%nt_range
      nstate    = option%nstate
      nboundary = boundary%nboundary 
      dt        = option%dt_out

      do ib = -nboundary, nboundary

        if (ib == 0) then
          cycle
        end if

        is1 = boundary%b2p(1, ib)
        is2 = boundary%b2p(2, ib)

        do js = 1, nstate
          state_sum = 0.0d0
          if (boundary%is_connected(js, is2)) then 
            state_sum = sum(Kijk(:, js, ib))
            if (state_sum < 1.0d-6) then
              write(iw,'("Check_Kijk> Warning.")')
              write(iw,'("K(",i0,2x,i0,2x,i0,") is zero")') js, is2, is1
            end if 
          end if
        end do

      end do

    end subroutine check_Kijk
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine Kijk_bin(output, option, boundary, Kijk)
!-----------------------------------------------------------------------
      implicit none

      type(s_output),   intent(in)    :: output 
      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(in)    :: boundary
      real(8),          intent(inout) :: Kijk(0:option%nt_range,  &
                                              option%nstate,      &
                                             -boundary%nboundary: &
                                              boundary%nboundary)

      ! I/O
      !
      integer :: io
      logical :: io_status

      ! Local
      !
      integer                :: nfile, nselect, nt_range
      character(len=MaxChar) :: fname
      
      ! Dummy
      !
      integer :: ifile, i, j, k, ib, isel, jsel

      ! Arrays
      !
      real(8), allocatable :: wrk(:) 


      ! Setup
      ! 
      nselect  = option%nselect
      nt_range = option%nt_range

      allocate(wrk(0:nt_range))
  
      if (.not. option%read_Kijk_bin &
    .and. .not. option%write_Kijk_bin) return 

      ! Read
      !
      if (option%read_Kijk_bin) then 
        ifile = 0 
        do while (.true.)
          write(fname,'(a,i4.4,".kbin")') trim(output%fhead), ifile + 1
          inquire(file=trim(fname), exist=io_status)
          if (io_status) then
            call open_file(fname, io, frmt = 'unformatted')

            read(io) nselect

            do isel = 1, nselect 
              read(io) i, j, k, wrk
              ib             = boundary%p2b(j, k)
              Kijk(:, i, ib) = wrk(:)

              write(iw,'("Read Kijk ", 3(i3,2x))') i, j, k
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
      if (option%write_Kijk_bin) then

        write(fname,'(a,i4.4,".kbin")') trim(output%fhead), option%out_id_kijk

        call open_file(fname, io, frmt = 'unformatted')

        write(io) nselect

        do isel = 1, nselect
          i      = option%sel_ijk(1, isel)
          j      = option%sel_ijk(2, isel)
          k      = option%sel_ijk(3, isel)
          ib     = boundary%p2b(j, k)
          wrk(:) = Kijk(:, i, ib)

          write(io) i, j, k, wrk

        end do 

        close(io)
      end if

      deallocate(wrk) 

    end subroutine Kijk_bin
!-----------------------------------------------------------------------
