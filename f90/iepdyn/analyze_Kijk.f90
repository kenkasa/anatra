!-----------------------------------------------------------------------
    subroutine update_Kijk_wo_normalize(option, state, fwrk)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_state),    intent(in)    :: state
      type(s_fwrk),     intent(inout) :: fwrk

      ! Local
      !
      integer :: nmol, nstep, nt_range, nt_sparse
      integer :: unp_id
      real(8) :: dt
      logical :: is_final, is_prod, is_dissoc

      ! Dummy 
      !
      integer :: istep, jstep, imol
      integer :: is, js, ks, ls, ib, jb, ik
      integer :: it, it_reac, it_diff, nt
      integer :: ista, iend, ireac
      integer :: iprod, idissoc


      ! Setup
      !
      nmol             = option%nmol
      nstep            = state%nstep
      unp_id           = state%unperturbed_id
      nt_range         = option%nt_range
      nt_sparse        = option%nt_sparse
      dt               = option%dt_out

      do imol = 1, nmol

        ! Search reaction time
        !
        ks = state%data(1, imol)
        js = ks

        ireac    = 0
        is_final = .false. 
        do istep = nt_sparse + 1, nstep, nt_sparse
          is = state%data(istep, imol)
          if (is /= js) then
            ireac  = ireac + 1
            
            if (ireac == 1) then
              ista              = istep - nt_sparse
              js                = is
              if (unp_id == js .or. unp_id == -1) then 
                fwrk%h(js, ks) = fwrk%h(js, ks) + 1.0d0
              end if
            else
              iend    = istep - nt_sparse 
              it_diff = iend - ista
              it_diff = int(it_diff / dble(nt_sparse))

              if (unp_id == js .or. unp_id == -1) then
                ik = fwrk%kmesh(is, js, ks)
                if (ik == 0) then
                  fwrk%nk              = fwrk%nk + 1
                  ik                   = fwrk%nk 
                  fwrk%kmesh(is, js, ks) = ik

                  if (ik >= fwrk%nkmax) then
                    write(iw,'("Update_Kijk_wo_Normalize> Error.")')
                    write(iw,'("# of elements exceeded array size of K &
                               &specified by nkmax = ", i0)') fwrk%nkmax
                    write(iw,'("Please specify larger nkmax value in \&option_param section")')
                    stop
                  end if

                end if 
                fwrk%K(it_diff, ik) = fwrk%K(it_diff, ik) + 1.0d0
              end if

              if (unp_id == is .or. unp_id == -1) then
                fwrk%h(is, js) = fwrk%h(is, js) + 1.0d0
              end if

              ls      = ks
              ks      = js
              js      = is
              ista    = istep - nt_sparse 

              !if (istep == nstep) then
              if (istep == nstep .or. istep + nt_sparse > nstep) then
                is_final = .true.
              end if

            end if

          end if
        end do

        ! New implementation
        !
        if (ireac == 0) cycle

        if (option%use_product_state) then
          is_prod = .false. 
          do iprod = 1, option%nproduct
            if (js == option%product_state_ids(iprod)) then
              is_prod = .true.
              exit 
            end if
          end do
          if (is_prod) cycle
        end if
!
!        if (option%use_dissociate_state) then
!          is_dissoc = .false.
!          do idissoc = 1, option%ndissoc
!            if (js == option%dissociate_state_ids(idissoc)) then
!              is_dissoc = .true.
!              exit 
!            end if
!          end do
!          if (is_dissoc) cycle 
!        end if
!
!        if (.not. option%is_dissoc(js)) then
!          if (unp_id == js .or. unp_id == -1) then
!            hit_count(js, ks) = hit_count(js, ks) - 1.0d0
!          end if
!        end if

        if (option%is_dissoc(js)) then
          if (is_final .and. (unp_id == js .or. unp_id == -1)) then
            fwrk%h(js, ks) = fwrk%h(js, ks) - 1.0d0
          end if
        else
          if (unp_id == js .or. unp_id == -1) then
            fwrk%h(js, ks) = fwrk%h(js, ks) - 1.0d0
          end if
        end if

        ! End of New implementation
        !

!        if (ireac == 0 .or. ireac == 1) then
!
!          if (option%use_product_state) then
!            do iprod = 1, option%nproduct
!              if (js == option%product_state_ids(iprod)) then
!                return
!              end if
!            end do 
!          end if
!
!          if (option%use_dissociate_state) then
!            do idissoc = 1, option%ndissoc
!              if (js == option%dissociate_state_ids(idissoc)) then
!                return
!              end if
!            end do 
!          end if
!
!          if (ireac == 1) then
!            if (.not. option%is_dissoc(js)) then
!              if (unp_id == js .or. unp_id == -1) then
!                hit_count(js, ks) = hit_count(js, ks) - 1.0d0
!              end if
!            end if
!            return
!          end if
!
!          !write(iw,'("Update_Kijk_wo_normalize> Error.")')
!          !write(iw,'("No reaction is observed.")')
!          !stop
!
!          !if (option%use_product_state) then
!          !  do iprod = 1, option%nproduct
!          !    if (js == option%product_state_ids(iprod)) then
!          !      return
!          !    end if
!          !  end do 
!          !else if (option%use_dissociate_state) then
!          !  do idissoc = 1, option%ndissoc
!          !    if (js == option%dissociate_state_ids(idissoc)) then
!          !      return
!          !    end if
!          !  end do 
!          !else
!          !  write(iw,'("Update_Kijk_wo_normalize> Error.")')
!          !  write(iw,'("No reaction is observed.")')
!          !  stop
!          !end if
!
!          !if (ireac == 1) then
!          !  if (.not. option%is_dissoc(js)) then 
!          !    hit_count(js, ks) = hit_count(js, ks) - 1.0d0
!          !  end if
!          !end if
!
!          !if (ireac == 1) then
!          !  if (.not. option%is_dissoc(js)) then
!          !    if (unp_id == js .or. unp_id == -1) then
!          !      hit_count(js, ks) = hit_count(js, ks) - 1.0d0
!          !    end if
!          !  end if
!          !end if
!
!        else
!
!          if (.not. option%is_dissoc(js)) then
!            if (unp_id == js .or. unp_id == -1) then
!              hit_count(js, ks) = hit_count(js, ks) - 1.0d0
!            end if
!          end if
!
!        end if

      end do

    end subroutine update_Kijk_wo_normalize
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine update_Kijk_from_hist(io, option, fwrk)
!-----------------------------------------------------------------------
      implicit none

      integer,          intent(in)    :: io
      type(s_option),   intent(in)    :: option
      type(s_fwrk),     intent(inout) :: fwrk

      ! Local
      !
      integer :: nmol, nstep, nt_range, nt_sparse
      integer :: unp_id
      real(8) :: dt
      logical :: is_final, is_prod, is_dissoc

      ! Dummy 
      !
      character(len=MaxChar) :: line, typ
      integer                :: istep
      integer                :: is, js, ks, ik
      real(8)                :: val

      ! Setup
      !

      do while (.true.)
        read(io,'(a)',end=100) line
        if (line(1:3) == 'END') go to 100

        read(line,*) typ
        if (trim(typ) == 'H') then
          read(line,*) typ, js, ks, val
          fwrk%h(js, ks) = fwrk%h(js, ks) + val
        else if (trim(typ) == 'K') then
          read(line,*) typ, is, js, ks, istep, val
          ik = fwrk%kmesh(is, js, ks)
          if (ik == 0) then
            fwrk%nk              = fwrk%nk + 1
            ik                   = fwrk%nk
            fwrk%kmesh(is, js, ks) = ik

            if (ik >= fwrk%nkmax) then
              write(iw,'("Update_Kijk_wo_Normalize> Error.")')
              write(iw,'("# of elements exceeded array size of K &
                         &specified by nkmax = ", i0)') fwrk%nkmax
              write(iw,'("Please specify larger nkmax value in \&option_param section")')
              stop
            end if

          end if
          fwrk%K(istep, ik) = fwrk%K(istep, ik) + val 
        end if 
      end do
 100  return 

    end subroutine update_Kijk_from_hist
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine convert_Kijk_arrays(option, boundary, fwrk, f)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(in)    :: boundary
      type(s_fwrk),     intent(inout) :: fwrk
      type(s_func),     intent(inout) :: f

      ! Local
      integer :: nstate, nboundary, nt_range

      ! Dummy
      integer :: ib, jb, is1, is2, js1, js2, ik
      integer :: inflx


      ! Setup
      !
      nstate    = option%nstate
      nboundary = boundary%nboundary
      nt_range  = option%nt_range

      ! Allocate
      !
      if (.not. allocated(f%K)) then
        allocate(f%K(0:nt_range, nstate, -nboundary:nboundary))
        allocate(f%hit_count(-nboundary:nboundary))
      end if
      f%K         = 0.0d0
      f%hit_count = 0.0d0

      ! Convert K-arrays
      !
      do ib = -nboundary, nboundary
        if (ib == 0) cycle
        is1             = boundary%b2p(1, ib)
        is2             = boundary%b2p(2, ib)
        f%hit_count(ib) = fwrk%h(is2, is1)

        do inflx = 1, boundary%n_influx_boundary(is1)
          jb   = boundary%influx_boundary(inflx, is1)
          js1  = boundary%b2p(1, jb)
          js2  = boundary%b2p(2, jb)
          ik   = fwrk%kmesh(is2, is1, js1)
          if (ik == 0) cycle 
          f%K(:, is2, jb) = fwrk%K(:, ik)
        end do 
      end do

    end subroutine convert_Kijk_arrays
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine normalize_Kijk(option, boundary, f, verbose)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),    intent(in)    :: option
      type(s_boundary),  intent(in)    :: boundary
      type(s_func),      intent(inout) :: f
      logical, optional, intent(in)    :: verbose 


      ! Local
      !
      integer :: nt_range, nstate, nboundary
      real(8) :: dt
      logical :: vb

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

      vb = .true.
      if (present(verbose)) vb = verbose

      ! Normalize
      !
      do ib = -nboundary, nboundary

        if (ib == 0) then
          cycle
        end if

        is1       = boundary%b2p(1, ib)
        is2       = boundary%b2p(2, ib)
        state_sum = f%hit_count(ib)

        if (state_sum < 1.0d-10) then
          f%K(:, :, ib) = 0.0d0
        else 
          f%K(:, :, ib) = f%K(:, :, ib) / (state_sum * dt) 
        end if
      end do

      ! Report Hit counts
      !
      if (vb) then
        write(iw,'("Normalize_Kijk> Summary of Hit counts")')
        write(iw,'("I <--- J")')
        do ib = -nboundary, nboundary
          if (ib == 0) then
            cycle
          end if
        
          is1 = boundary%b2p(1, ib)
          is2 = boundary%b2p(2, ib)
        
          write(iw,'(i5,i5," : ", f20.10)') is2, is1, f%hit_count(ib) 
        end do
        write(iw,*)
      end if
      

    end subroutine normalize_Kijk
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine write_Kijk(output, option, boundary, f)
!-----------------------------------------------------------------------
      implicit none

      type(s_output),   intent(in)    :: output
      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(in)    :: boundary
      type(s_func),     intent(inout) :: f

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
              write(io,'(e15.7,2x)', advance = 'no') f%K(istep, js, ib)
            end if 
          end do
        end do
        write(io,*)
      end do
      close(io)

    end subroutine write_Kijk
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine check_Kijk(option, boundary, f)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(in)    :: boundary
      type(s_func),     intent(inout) :: f 

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

        state_sum = 0.0d0
        do js = 1, nstate
          if (boundary%is_connected(js, is2)) then
            state_sum = state_sum + sum(f%K(:, js, ib)) * dt
          end if
        end do
        write(iw,'("Kint value at Boundary ", i5, i5, " : ", f15.7)') &
              is2, is1, state_sum
      end do

    end subroutine check_Kijk
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine output_Kijk_hist(option, output, fwrk, is_gen)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      type(s_output), intent(in)    :: output
      type(s_fwrk),   intent(in)    :: fwrk
      logical,        intent(inout) :: is_gen

      ! I/O
      !
      integer :: io

      ! Local
      !
      character(len=MaxChar) :: fname
      integer                :: nt_range
      integer                :: nstate
      real(8)                :: hsum, val

      ! Dummy
      !
      integer :: is, js, ks, istep, ik


      ! Setup
      !
      nstate   = option%nstate
      nt_range = option%nt_range

      hsum = sum(fwrk%h(:, :))
      if (hsum < 0.999d0) then
        return
      end if

      write(fname,'(a,".khist")') trim(output%fhead)
      call open_file(fname, io)

      write(io,'("KIJK")')

      do ks = 1, nstate
      do js = 1, nstate
        val = fwrk%h(js, ks)
        if (val >= 0.999d0) then
          write(io, '("H", 2x, i5, 2x, i5, 2x, f20.10)') js, ks, val 
        end if
      end do
      end do

      do ks = 1, nstate
      do js = 1, nstate
      do is = 1, nstate
        ik = fwrk%kmesh(is, js, ks)
        if (ik == 0) cycle
        do istep = 0, nt_range
          val = fwrk%K(istep, ik) 
          if (val >= 0.999d0) then
            write(io, '("K", 2x, i5, 2x, i5, 2x, i5, 2x, i10, 2x, f20.10)') is, js, ks, istep, val
          end if
        end do 
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

      write(io,'("KIJK")')

      do ks = 1, nstate
      do js = 1, nstate
        val = fwrk%h(js, ks)
        if (val >= 0.999d0) then
          write(io, '("H", 2x, i5, 2x, i5, 2x, f20.10)') js, ks, val 
        end if
      end do
      end do

      do ks = 1, nstate
      do js = 1, nstate
      do is = 1, nstate
        ik = fwrk%kmesh(is, js, ks)
        if (ik == 0) cycle
        do istep = 0, nt_range
          val = fwrk%K(istep, ik) 
          if (val >= 0.999d0) then
            write(io, '("K", 2x, i5, 2x, i5, 2x, i5, 2x, i10, 2x, f20.10)') is, js, ks, istep, val
          end if
        end do 
      end do
      end do 
      end do
      write(io, '("END")')

      close(io)


    end subroutine output_Kijk_hist
!-----------------------------------------------------------------------
