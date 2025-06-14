!-----------------------------------------------------------------------
    subroutine calc_Mij_from_Kijk(option, boundary, Kijk, Mij)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(in)    :: boundary 
      real(8),          intent(in)    :: Kijk(0:option%nt_range, &
                                             option%nstate,     &
                                             -boundary%nboundary:boundary%nboundary)
      real(8),          intent(inout) :: Mij(0:option%nt_range, &
                                            -boundary%nboundary:boundary%nboundary) 

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


      ! Setup
      !
      nmol      = option%nmol
      nstate    = option%nstate
      nboundary = boundary%nboundary 
      nt_range  = option%nt_range
      nt_sparse = option%nt_sparse
      dt        = option%dt_out

      Mij = 1.0d0 
      do istep = 0, nt_range

        do ib = -nboundary, nboundary
          if (ib == 0) cycle

          is1 = boundary%b2p(1, ib)
          is2 = boundary%b2p(2, ib)

          do js = 1, nstate
            if (.not. boundary%is_connected(js, is2)) cycle
            do jstep = 0, istep - 1
              Mij(istep, ib) = Mij(istep, ib) - dt * Kijk(jstep, js, ib)  
            end do
          end do

        end do

      end do

    end subroutine calc_Mij_from_Kijk
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine write_Mij(output, option, boundary, Mij)
!-----------------------------------------------------------------------
      implicit none

      type(s_output),   intent(in) :: output
      type(s_option),   intent(in) :: option
      type(s_boundary), intent(in) :: boundary
      real(8),          intent(in) :: Mij(0:option%nt_range, &
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

          write(io,'(f20.10)', advance = 'no') Mij(istep, ib)

        end do
        write(io,*)
      end do

      close(io) 

    end subroutine write_Mij
!-----------------------------------------------------------------------
