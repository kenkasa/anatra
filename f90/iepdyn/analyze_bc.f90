!-----------------------------------------------------------------------
    subroutine set_reflection(output, option, boundary, f)
!-----------------------------------------------------------------------
      implicit none

      type(s_output),   intent(in)    :: output
      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(inout) :: boundary
      type(s_func),     intent(inout) :: f

      ! I/O
      !
      integer :: io

      ! Local
      !
      integer                :: nstate, nt_range, nt_extend, nboundary
      real(8)                :: dt
      character(len=MaxChar) :: fname

      ! Dummy
      !
      integer :: is, js, is1, is2, ib, jb, id, istep, jstep, iref
      integer :: js1, js2

      ! Arrays
      !


      if (.not. option%use_reflection_state) return

      ! Setup
      !
      nstate    = option%nstate
      nt_range  = option%nt_range
      nt_extend = option%nt_extend
      dt        = option%dt_out
      nboundary = boundary%nboundary

      if (.not. allocated(boundary%conv_direc)) then
        allocate(boundary%conv_direc(-nboundary:nboundary))
      end if
      boundary%conv_direc = .false.

      do iref = 1, option%nreflect
        is = option%reflection_state_ids(iref)
        do js = 1, nstate
          ib = boundary%p2b(is, js)
          if (ib == 0) cycle

          ! Erase transitions from reflection state
          !
          f%K(:, :, ib) = 0.0d0
          f%M(:, ib)    = 0.0d0

          ! Invert directions
          !
          boundary%conv_direc(ib) = .true.

        end do
        
      end do
      
    end subroutine set_reflection 
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine set_product(output, option, boundary, f)
!-----------------------------------------------------------------------
      implicit none

      type(s_output),   intent(in)    :: output
      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(inout) :: boundary
      type(s_func),     intent(inout) :: f

      ! I/O
      !
      integer :: io

      ! Local
      !
      integer                :: nstate, nt_range, nt_extend, nboundary
      real(8)                :: dt
      character(len=MaxChar) :: fname

      ! Dummy
      !
      integer :: is, js, is1, is2, ib, jb, id, istep, jstep, iprod
      integer :: js1, js2

      ! Arrays
      !


      if (.not. option%use_product_state) return

      ! Setup
      !
      nstate    = option%nstate
      nt_range  = option%nt_range
      nt_extend = option%nt_extend
      dt        = option%dt_out
      nboundary = boundary%nboundary

      do iprod = 1, option%nproduct
        is = option%product_state_ids(iprod)

        do js = 1, nstate

          ib = boundary%p2b(is, js)

          if (ib == 0) cycle

          ! Erase transitions from reflection product i 
          ! (X<---i<---j, X<---j<---i)
          !
          f%K(:, :, ib)  = 0.0d0
          f%M(:, ib)     = 1.0d0
          f%K(:, :, -ib) = 0.0d0
          f%M(:, -ib)    = 0.0d0

        end do
        
      end do

    end subroutine set_product 
!-----------------------------------------------------------------------
