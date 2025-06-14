!-----------------------------------------------------------------------
    subroutine calc_P0_from_Rij(option, boundary, Rij, P0)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(in)    :: boundary
      real(8),          intent(in)    :: Rij(0:option%nt_range, &
                                             option%nstate,     &
                                             option%nstate)
      real(8),          intent(inout) :: P0(0:option%nt_range,  &
                                          option%nstate)


      ! Local
      !
      integer :: nt_range, nstate
      real(8) :: dt

      ! Dummy
      !
      integer :: istep, is, js
      real(8) :: rval, weight 


      ! Setup
      !
      nstate   = option%nstate
      nt_range = option%nt_range
      nstate   = option%nstate
      dt       = option%dt_out

      P0 = 0.0d0
      do is = 1, nstate
        if (.not. option%is_initial(is)) cycle

        weight    = option%state_weight(is) 
        P0(:, is) = weight
        do js = 1, nstate
          if (.not. boundary%is_connected(js, is)) cycle

          do istep = 1, nt_range
            P0(istep, is) = P0(istep, is) - dt * sum(Rij(0:istep - 1, js, is))
          end do
        end do

      end do


    end subroutine calc_P0_from_Rij 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine write_P0(output, option, P0)
!-----------------------------------------------------------------------
      implicit none

      type(s_output),   intent(in) :: output
      type(s_option),   intent(in) :: option
      real(8),          intent(in) :: P0(0:option%nt_range, &
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

      write(fname, '(a,".P0")') trim(output%fhead)
      call open_file(fname, io)
      
      id = 1 
      do is = 1, nstate
        if (option%is_initial(is)) then
          write(io,'("# Col. ", i0, " P0 : ", i0)') id + 1, is 
          id = id + 1 
        end if
      end do

      do istep = 0, nt_range
        write(io,'(f20.10)', advance = 'no') option%dt_out * istep 
        do is = 1, nstate
          if (option%is_initial(is)) then
            write(io,'(f20.10)', advance = 'no') P0(istep, is) 
          end if
        end do
        write(io,*)
      end do

      close(io)

    end subroutine write_P0
!-----------------------------------------------------------------------

