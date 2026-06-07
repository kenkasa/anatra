!-----------------------------------------------------------------------
    subroutine calc_Mjk_from_Kijk(option, boundary, f)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(in)    :: boundary
      type(s_func),     intent(inout) :: f 

      ! Local
      !
      integer :: nmol, nstate, nt_range, nt_sparse, nboundary
      real(8) :: dt

      ! Dummy 
      !
      integer :: istep, jstep, imol
      integer :: is, js, is1, is2, ib
      integer :: it, it_reac, it_diff, nt
      integer :: ista, iend, ireac

      ! Arrays
      !
      real(8), allocatable :: ksum(:, :)


      ! Setup
      !
      nmol      = option%nmol
      nstate    = option%nstate
      nboundary = boundary%nboundary 
      nt_range  = option%nt_range
      nt_sparse = option%nt_sparse
      dt        = option%dt_out

      ! Allocate
      !
      if (.not. allocated(f%M)) then
        allocate(f%M(0:nt_range, -nboundary:nboundary))
        f%M = 0.0d0
      end if

      ! Allocate work space
      !
      if (.not. allocated(ksum)) then
        allocate(ksum(nstate, -nboundary:nboundary))
      end if
      ksum = 0.0d0

      ! Calculate M from K
      !
      f%M = 1.0d0 
      do istep = 0, nt_range

        do ib = -nboundary, nboundary
          if (ib == 0) cycle

          is1 = boundary%b2p(1, ib)
          is2 = boundary%b2p(2, ib)

          do js = 1, nstate
            if (.not. boundary%is_connected(js, is2)) cycle

            if (istep > 0) then
              ksum(js, ib)   = ksum(js, ib) - dt * f%K(istep - 1, js, ib)
              f%M(istep, ib) = f%M(istep, ib) + ksum(js, ib) 
            end if 

            if (f%M(istep, ib) < 0.0d0) then
              if (abs(f%M(istep, ib)) > 1.0d-3) then
                write(iw,'("Calc_Mjk_from_Kijk> Error.")')
                write(iw,'("Negative population has been detected. stop")')
                stop
              else
                f%M(istep, ib) = 0.0d0
              end if 
            end if
          end do

        end do

      end do

      deallocate(ksum)

    end subroutine calc_Mjk_from_Kijk
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine write_Mjk(output, option, boundary, f)
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

      write(fname, '(a,".Mij")') trim(output%fhead)
      call open_file(fname, io)
      write(io,'("# ib")')

      id = 1
      do ib = -nboundary, nboundary
        if (ib == 0) then
          cycle
        end if
        is1 = boundary%b2p(1, ib)
        is2 = boundary%b2p(2, ib)

        write(io,'("# Col. ", i0, " Mij : ", i0, &
                   " <= ", i0, " < ", i0)') id + 1, is2, is2, is1
        id = id + 1 

      end do

      do istep = 0, nt_range
        write(io,'(f20.10)', advance = 'no') option%dt_out * istep
        do ib = -nboundary, nboundary
          if (ib == 0) then
            cycle
          end if
          write(io,'(e15.7,2x)', advance = 'no') f%M(istep, ib)
        end do
        write(io,*)
      end do

      close(io) 

    end subroutine write_Mjk
!-----------------------------------------------------------------------
