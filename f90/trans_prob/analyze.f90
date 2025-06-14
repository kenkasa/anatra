!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_ctrl
  use mod_bootstrap
  use mod_random

  ! constants
  !

  ! structures
  !
  type :: s_state
    integer :: nmol
    integer :: nstep
    integer, allocatable :: data(:, :) ! (nstep, nmol) 
  end type s_state

  type :: s_boundary
    integer :: nboundary
    logical, allocatable :: is_connected(:, :)
    integer, allocatable :: p2b(:, :)  ! pair        => boundary id
    integer, allocatable :: b2p(:, :)  ! boundary id => pair
    logical, allocatable :: conv_direc(:)
  end type s_boundary

  type :: s_booteach
    integer, allocatable :: rand(:, :)
    real(8), allocatable :: func(:, :)
  end type s_booteach

  type :: s_bootave
    real(8), allocatable :: ave(:)
    real(8), allocatable :: err(:) 
  end type s_bootave

  ! subroutines
  !
  public  :: analyze

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option, timegrid)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),    intent(in)    :: input
      type(s_output),   intent(in)    :: output
      type(s_option),   intent(in)    :: option
      type(s_timegrid), intent(in)    :: timegrid 

      ! I/O
      !
      integer                :: io
      character(len=MaxChar) :: fname 

      ! Local
      !
      integer :: ndim, nstep, nfile, nstate
      integer :: nt_range
      integer :: nboundary

      type(s_boundary) :: boundary

      ! Dummy
      !
      integer :: ifile, istep, is, js, is1, is2, ib, id, idir

      ! Arrays
      !
      type(s_cv),    allocatable :: cv(:)
      type(s_state), allocatable :: state(:)

      real(8),       allocatable :: Rij(:, :, :), Rij_int(:, :, :)
      real(8),       allocatable :: P0(:, :)

      real(8),       allocatable :: Kijk(:, :, :), hit_count(:)
      real(8),       allocatable :: Mij(:, :)


      ! Setup
      !
      ndim     = option%ndim * option%nmol 
      nfile    = input%ncv
      nstate   = option%nstate
      nt_range = option%nt_range  
      allocate(cv(nfile), state(nfile))

      ! Read CV files
      !
      write(iw,*)
      write(iw,'("Analyze> Read CV file")')
      do ifile = 1, nfile 
        call read_cv  (input%fcv(ifile), ndim, cv(ifile))
        call get_state(option, cv(ifile), state(ifile))
      end do
      write(iw,'(">> Done")')

      ! Check state connectivity
      !
      write(iw,*)
      write(iw,'("Analyze> Get Connectivity")')
      do ifile = 1, nfile
        call get_state_connectivity(option, state(ifile), boundary)
      end do
      call show_state_connectivity(option, boundary)
      write(iw,'(">> Done")')

      ! Define boundary
      !
      write(iw,*)
      write(iw,'("Analyze> Define Boundary")')
      call define_boundary(option, boundary)
      nboundary = boundary%nboundary
      write(iw,'(">> Done")')

      ! Main part
      !
      if (option%kinetic_mode == KineticModeTransition) then

        ! Not implmented yet


      else if (option%kinetic_mode == KineticModeReaction) then


        write(iw,*)
        write(iw,'("Analyze> KineticMode = REACTION")')
        write(iw,'(">> Calculate Rij(t) and Kijk(t)")')

        ! Calculate Rij(t) and P0(t)
        !
        write(iw,*)
        write(iw,'("o Rij(t)")')
        allocate(Rij    (0:nt_range, nstate, nstate))
        allocate(Rij_int(0:nt_range, nstate, nstate))
        allocate(P0     (0:nt_range, nstate))

        Rij = 0.0d0
        if (.not. option%read_Rij_bin) then
          do ifile = 1, nfile
            call calc_Rij_wo_normalize(option, state(ifile), Rij) 
          end do
          call normalize_Rij(option, Rij)
        end if
        call Rij_bin             (output, option, boundary, Rij)   
        call running_integral_Rij(option, Rij, Rij_int)
        call calc_P0_from_Rij    (option, boundary, Rij, P0)

        write(iw,'(">> Done")')

        ! Calculate Kijk(t) and Mij(t)
        !
        write(iw,*)
        write(iw,'("o Kijk(t)")')
        allocate(Kijk(0:nt_range, nstate, -nboundary:nboundary))
        allocate(Mij (0:nt_range, -nboundary:nboundary))
        allocate(hit_count(-nboundary:nboundary))

        Kijk       = 0.0d0
        hit_count = 0.0d0
        do ifile = 1, nfile
          call calc_Kijk_wo_normalize(option, boundary, state(ifile), Kijk, hit_count)
        end do
        call normalize_Kijk(option, boundary, Kijk, hit_count)

        if (option%read_Kijk_bin .or. option%write_Kijk_bin) then 
          call Kijk_bin    (output, option, boundary, Kijk)
          if (.not. option%write_Kijk_bin) then 
            call renormalize_Kijk (option, boundary, Kijk)
          end if
        end if 

        call calc_Mij_from_Kijk(option, boundary, Kijk, Mij)

        ! Check Kijk
        !
        call check_Kijk(option, boundary, Kijk)

        write(iw,'(">> Done")')

        ! Output
        !
        call write_Rij   (output, option, boundary, Rij, Rij_int)
        call write_P0    (output, option, P0)
        call write_Kijk   (output, option, boundary, Kijk)
        call write_Mij    (output, option, boundary, Mij)

        if (option%extrapolate) then

          write(iw,*)
          write(iw,'("Analyze> Start propagation")')

          ! Setup Boundary conditions 
          !
          call set_reflection(output, option, boundary, Kijk, Mij)
          call set_product   (output, option, boundary, Kijk, Mij)
          
          ! Extend timescale 
          !
          call reacdyn_tcf(output, option, boundary, Rij, P0, Kijk, Mij)

          write(iw,'(">> Done")')

        end if 

      end if


    end subroutine analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_state(option, cv, state) 
!-----------------------------------------------------------------------
      implicit none
      
      type(s_option) :: option
      type(s_cv)     :: cv
      type(s_state)  :: state

      ! Local
      !
      integer :: nmol, ndim, nstep 

      ! Dummy
      !
      integer :: istep, imol, idim, jdim, ia, istate

      ! Array
      !
      real(8), allocatable :: val(:)


      ! Setup
      !
      nmol  = option%nmol
      ndim  = option%ndim
      nstep = cv%nstep

      state%nmol  = nmol
      state%nstep = nstep

      allocate(state%data(nstep, nmol))
      allocate(val(ndim))

      state%data = - 1 
      do istep = 1, nstep
        jdim = 0
        do imol = 1, nmol
          do idim = 1, ndim
            jdim      = jdim + 1
            val(idim) = cv%data(jdim, istep)
          end do

          do istate = 1, option%nstate

            ia = 0
            do idim = 1, ndim
              if (val(idim) >= option%state_def(1, idim, istate) &
            .and. val(idim) <  option%state_def(2, idim, istate)) then
                ia = ia + 1
              end if 
            end do

            if (ia == ndim) then
              state%data(istep, imol) = istate  
            end if

          end do

          if (state%data(istep, imol) == - 1) then
            write(iw,'("Get_State> Error.")')
            write(iw,'("State determination failed. Please check state definition.")')
            stop
          end if

        end do
      end do 

    end subroutine get_state
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_state_connectivity(option, state, boundary)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_state),    intent(in)    :: state 
      type(s_boundary), intent(inout) :: boundary 

      ! Local
      !
      integer :: nstate, nmol, nstep

      ! Dummy
      !
      integer :: istep, imol
      integer :: snow, sprev


      ! Setup
      !
      nmol   = option%nmol
      nstate = option%nstate
      nstep  = state%nstep

      if (.not. allocated(boundary%is_connected)) then
        allocate(boundary%is_connected(nstate, nstate))
        boundary%is_connected = .false.
      end if

      do imol = 1, nmol
        do istep = 1, nstep
          snow = state%data(istep, imol)
          if (istep == 1) then
            sprev = snow
            cycle
          end if

          if (snow /= sprev) then
            boundary%is_connected(snow,  sprev) = .true. 
            boundary%is_connected(sprev, snow)  = .true. 
          end if

          sprev = snow

        end do
      end do

    end subroutine get_state_connectivity
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine show_state_connectivity(option, boundary)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in) :: option
      type(s_boundary), intent(in) :: boundary 

      ! Local
      !
      integer :: nstate, nmol

      ! Dummy
      !
      integer :: istep, imol, is, js
      integer :: snow, sprev


      ! Setup
      !
      nmol   = option%nmol
      nstate = option%nstate

      do js = 1, nstate - 1
        do is = js + 1, nstate
          if (boundary%is_connected(is, js)) then
            write(iw,'(i0,2x,i0)') is, js
          end if
        end do 
      end do

    end subroutine show_state_connectivity
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine define_boundary(option, boundary)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(inout) :: boundary 

      ! Local
      !
      integer :: nstate, nboundary

      ! Dummy
      !
      integer :: is, js 


      ! Setup
      !
      nstate = option%nstate

      allocate(boundary%p2b(nstate, nstate))

      boundary%p2b = 0
      nboundary    = 0
      do js = 1, nstate - 1
        do is = js + 1, nstate
          if (boundary%is_connected(is, js)) then
            nboundary            =  nboundary + 1
            boundary%p2b(is, js) =  nboundary
            boundary%p2b(js, is) = -nboundary
          end if
        end do
      end do

      boundary%nboundary = nboundary

      allocate(boundary%b2p(2, -nboundary:nboundary))

      nboundary    = 0
      boundary%b2p = 0
      do js = 1, nstate - 1
        do is = js + 1, nstate
          if (boundary%is_connected(is, js)) then
            nboundary                   = nboundary + 1
            boundary%b2p(1,  nboundary) = js
            boundary%b2p(2,  nboundary) = is
            boundary%b2p(1, -nboundary) = is
            boundary%b2p(2, -nboundary) = js
          end if
        end do
      end do

    end subroutine define_boundary 
!-----------------------------------------------------------------------

    include 'analyze_Rij.f90'
    include 'analyze_P0.f90'
    include 'analyze_Kijk.f90'
    include 'analyze_Mij.f90'
    include 'analyze_bc.f90'
    include 'analyze_reacdyn.f90'

end module mod_analyze
!=======================================================================
