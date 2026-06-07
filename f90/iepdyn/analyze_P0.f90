!-----------------------------------------------------------------------
    subroutine calc_P0_from_Rij(option, boundary, f)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(in)    :: boundary
      type(s_func),     intent(inout) :: f

      ! Local
      !
      integer :: nt_range, nstate
      real(8) :: dt

      ! Dummy
      !
      integer :: istep, is, js
      real(8) :: rval, weight

      ! Arrays
      !
      real(8), allocatable :: rsum(:, :) 


      ! Setup
      !
      nstate   = option%nstate
      nt_range = option%nt_range
      nstate   = option%nstate
      dt       = option%dt_out

      allocate(rsum(nstate, nstate))
      rsum = 0.0d0

      f%P0 = 0.0d0
      do is = 1, nstate
        if (.not. option%is_initial(is)) cycle

        weight      = option%state_weight(is) 
        f%P0(:, is) = weight
        do js = 1, nstate
          if (.not. boundary%is_connected(js, is)) cycle

          do istep = 1, nt_range
            rsum(js, is)    = rsum(js, is) - dt * f%R(istep - 1, js, is)
            f%P0(istep, is) = f%P0(istep, is) + rsum(js, is)
          end do
        end do

      end do

      deallocate(rsum)

    end subroutine calc_P0_from_Rij 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine write_P0(output, option, f)
!-----------------------------------------------------------------------
      implicit none

      type(s_output),   intent(in)    :: output
      type(s_option),   intent(in)    :: option
      type(s_func),     intent(inout) :: f

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
            write(io,'(e15.7,2x)', advance = 'no') f%P0(istep, is) 
          end if
        end do
        write(io,*)
      end do

      close(io)

    end subroutine write_P0
!-----------------------------------------------------------------------

