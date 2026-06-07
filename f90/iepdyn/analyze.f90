!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_ctrl
  use mod_random

  ! Constants
  !

  ! Structures
  !
  type :: s_state
    integer :: nmol
    integer :: nstep
    integer :: unperturbed_id = -1
    integer :: use_for_Rij    = .true.
    integer, allocatable :: data(:, :) ! (nstep, nmol)
  end type s_state

  type :: s_boundary
    integer :: nboundary
    logical, allocatable :: is_connected(:, :)
    integer, allocatable :: p2b(:, :)  ! pair        => boundary id
    integer, allocatable :: b2p(:, :)  ! boundary id => pair
    integer, allocatable :: n_influx_boundary(:)
    integer, allocatable :: influx_boundary(:, :)
    logical, allocatable :: conv_direc(:)

    logical, allocatable :: is_cQij(:)
    real(8), allocatable :: cQij(:)
  end type s_boundary

  type :: s_func
    real(8), allocatable :: K(:, :, :)
    real(8), allocatable :: M(:, :)
    real(8), allocatable :: P0(:, :)
    real(8), allocatable :: R(:, :, :), Rint(:, :, :)
    real(8), allocatable :: hit_count(:) 
  end type s_func

  type :: s_fwrk
    integer :: nkmax = 50000
    integer :: nk    = 0
    integer, allocatable :: kmesh(:, :, :)
    real(8), allocatable :: K(:, :)
    real(8), allocatable :: h(:, :)
  end type s_fwrk

  type :: s_infprop
    real(8), allocatable :: prob(:)
    real(8), allocatable :: fe(:)
    real(8), allocatable :: fe_pair(:, :)
  end type s_infprop

  type :: s_inpcond
    integer, allocatable :: unperturbed_ids(:)
    integer, allocatable :: use_for_Rij(:)
    integer, allocatable :: nfile_each_state(:)
    integer, allocatable :: nfile_to_be_read(:)
    integer, allocatable :: ista(:)
    integer, allocatable :: iend(:)
  end type s_inpcond

  ! subroutines
  !
  public  :: analyze

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, einput, output, option, timegrid)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),       intent(in)    :: input
      type(s_extra_input), intent(in)    :: einput
      type(s_output),      intent(in)    :: output
      type(s_option),      intent(in)    :: option
      type(s_timegrid),    intent(in)    :: timegrid 

      ! I/O
      !

      ! Local
      !
      type(s_func)     :: f
      type(s_inpcond)  :: ic
      type(s_boundary) :: boundary
      type(s_cv)       :: cv
      type(s_state)    :: state
      type(s_fwrk)     :: fwrk
      type(s_infprop)  :: ip

      integer :: iseed


      ! Setup
      !
      fwrk%nkmax = option%nkmax

      call get_seed         (iseed)
      call initialize_random(iseed)


      ! Read Unperturbed_ID file 
      ! (if use_perturbed_traj = .true.)
      !
      call setup_perturb(input, option, ic)

      ! Calculate Kernels
      !
      call calc_kernels(input, output, option, ic, boundary, f, fwrk, &
                        set_boundary = .true., verbose = .true.)

      ! Write functions involved in the integral equations (IEs) as inputs
      !
      call write_IEfunc(output, option, boundary, f)

      ! Setup boundary conditions 
      !
      call setup_boundary_cond(output, option, boundary, f)

      ! Setup constant Qij 
      ! (if use_constant_Qij = .true.)
      !
      call setup_constant_Qij(option, boundary)

      ! Evaluate steady-state properties 
      ! (if check_Pint = .true. or check_Steady = .true.)
      !
      call reacdyn_pint(output, option, boundary, f, ip, write_steady = .true.)

      ! Calculate Block-Average of steady-state properties
      ! (if check_blockave = .true.)
      !
      call calc_blockave(input, output, option, boundary, ic, fwrk)

      ! Conduct Cumulative analysis 
      ! (if check_cumulative = .true.)
      !
      call calc_cumulative(input, output, option, boundary, ic, fwrk)

      ! Conduct Sensitivity-Error analysis
      ! (if check_senserr = .true.)
      !
      call calc_senserr(input, output, option, boundary, ip, ic, fwrk)


      ! Extend timescale of Pj by solving IEs
      ! (if extrapolate = .true.)
      call reacdyn_tcf(output, option, boundary, f)

    end subroutine analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine setup_perturb(input, option, ic) 
!-----------------------------------------------------------------------
      implicit none

      type(s_input),   intent(in)    :: input      
      type(s_option),  intent(in)    :: option
      type(s_inpcond), intent(inout) :: ic 

      ! Local
      !
      integer :: nstate, nfile 

      ! Dummy
      !
      integer :: is, ifile 


      if (.not. option%use_perturbed_traj) return

      write(iw,*)
      write(iw,'("Analyze> Read f_unperturbed_id file")')
      write(iw,'("Note: unperturbed state info.&
                & is used only if use_perturbed_traj = .true.")')

      ! Setup
      !
      nfile  = input%ncv
      nstate = option%nstate 

      ! Get unperturbed state id for each file
      !
      allocate(ic%unperturbed_ids(nfile), ic%use_for_Rij(nfile))
      call read_f_unperturbed_id(option, nfile, ic%unperturbed_ids, ic%use_for_Rij)
      do ifile = 1, nfile
        write(iw,'(3i10)') ifile, ic%unperturbed_ids(ifile), ic%use_for_Rij(ifile)
      end do
      
      ! Get # of files for each unperturbed state
      !
      allocate(ic%nfile_each_state(nstate))
      allocate(ic%nfile_to_be_read(nstate))
      allocate(ic%ista(nstate), ic%iend(nstate))

      ic%nfile_each_state = 0
      do ifile = 1, nfile
        is                      = ic%unperturbed_ids(ifile)
        ic%nfile_each_state(is) = ic%nfile_each_state(is) + 1
      end do

      ! Initial setting (will not be changed if check_senserr = .false.)
      !
      do is = 1, nstate
        ic%nfile_to_be_read(is) = ic%nfile_each_state(is)  
      end do

      ! Initial setting (will not be changed if check_blockave = .false.)
      !
      do is = 1, nstate
        ic%ista(is) = 1
        ic%iend(is) = ic%nfile_each_state(is)
      end do
!
    end subroutine setup_perturb 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine write_IEfunc(output, option, boundary, f) 
!-----------------------------------------------------------------------
      implicit none
      
      type(s_output),   intent(in)    :: output
      type(s_option),   intent(in)    :: option 
      type(s_boundary), intent(in)    :: boundary
      type(s_func),     intent(inout) :: f 


      write(iw,*)
      write(iw,'("Analyze> Write functions involved in IEs")')
      call write_Rij   (output, option, boundary, f)
      call write_P0    (output, option, f)
      call write_Kijk  (output, option, boundary, f)
      call write_Mjk   (output, option, boundary, f)
      write(iw,'(">> Done")')

!
    end subroutine write_IEfunc
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine setup_boundary_cond(output, option, boundary, f) 
!-----------------------------------------------------------------------
      implicit none
      
      type(s_output),   intent(in)    :: output
      type(s_option),   intent(in)    :: option 
      type(s_boundary), intent(inout) :: boundary
      type(s_func),     intent(inout) :: f 


      write(iw,*)
      write(iw,'("Analyze> Setup boundary conditions")')
      call set_reflection(output, option, boundary, f) 
      call set_product   (output, option, boundary, f)
      write(iw,'(">> Done")')
!
    end subroutine setup_boundary_cond 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine setup_constant_Qij(option, boundary) 
!-----------------------------------------------------------------------
      implicit none
      
      type(s_option),   intent(in)    :: option 
      type(s_boundary), intent(inout) :: boundary

      ! I/O
      !
      integer :: io

      ! Local
      !
      integer :: nb 

      ! Dummy
      !
      integer                :: is, js, ib
      real(8)                :: val
      character(len=MaxChar) :: line


      if (.not. option%use_constant_Qij) return

      write(iw,*)
      write(iw,'("Analyze> Setup constant Qij")')

      ! Setup
      !
      nb = boundary%nboundary

      ! Allocate
      !
      if (.not. allocated(boundary%is_cQij)) then
        allocate(boundary%is_cQij(nb))
        allocate(boundary%cQij(nb))
      end if
      boundary%is_cQij = .false.
      boundary%cQij    = 0.0d0

      ! Read
      !
      call open_file(option%f_cQij, io)

      do while (.true.)
        read(io, '(a)', end = 100) line
        if (line(1:1) == "#") cycle
        read(line, *) is, js, val

        ib = boundary%p2b(is, js)
        if (ib == 0) then
          write(iw,'("Setup_Constant_Qij> Error.")')
          write(iw,'("Constant-Qij value at unconnected (i, j) &
                     &has been detected")')
          stop
        end if
        boundary%is_cQij(ib) = .true.
        boundary%cQij(ib)    = val

        write(iw,'(i5,2x,i5,2x,e15.7)') is, js, boundary%cQij(ib)
      end do

 100  close(io)
      write(iw,'(">> Done")')
!
    end subroutine setup_constant_Qij
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine calc_blockave(input, output, option, boundary, ic, fwrk) 
!-----------------------------------------------------------------------
      implicit none
     
      type(s_input),    intent(in)    :: input 
      type(s_output),   intent(in)    :: output
      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(inout) :: boundary
      type(s_inpcond),  intent(inout) :: ic 
      type(s_fwrk),     intent(inout) :: fwrk 

      ! I/O
      !
      character(len=MaxChar) :: fname

      ! Local
      !
      integer      :: nstate, block_size
      type(s_func) :: fs

      ! Dummy
      !
      integer :: is, iblock

      ! Array
      !
      type(s_infprop), allocatable :: ipse(:)


      if (.not. option%check_blockave) return

      write(iw,*)
      write(iw,'("Analyze> Start Block analysis of steady-state properties")')

      ! Setup
      !
      nstate = option%nstate

      ! Allocate
      !
      allocate(ipse(option%nblock))

      do iblock = 1, option%nblock
      
        write(iw,'("Block ", i0)') iblock

        ! Set Start and End trajectories 
        !
        write(iw,'("  File block")')
        do is = 1, nstate
          block_size  = ic%nfile_each_state(is) / option%nblock
          ic%ista(is) = block_size * (iblock - 1) + 1
          ic%iend(is) = block_size * iblock
          write(iw,'(2x,"State ", i5, " : ", i0,2x,i0)') is, ic%ista(is), ic%iend(is) 
        end do
        write(iw,*)

        ! Calculate Kernels
        !
        call calc_kernels(input, output, option, ic, boundary, fs, fwrk, &
                          set_boundary = .false., verbose = .false.)

        ! Setup Boundary conditions
        !
        call set_reflection(output, option, boundary, fs) 
        call set_product   (output, option, boundary, fs)

        ! Calculate steady-state properties
        !
        write(fname,'(a,".steady.",i4.4)') trim(output%fhead), iblock
        call reacdyn_pint(output, option, boundary, fs, ipse(iblock), &
                          write_steady = .true., fname_out = fname)
      end do 

!
    end subroutine calc_blockave 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine calc_cumulative(input, output, option, boundary, ic, fwrk) 
!-----------------------------------------------------------------------
      implicit none
     
      type(s_input),    intent(in)    :: input 
      type(s_output),   intent(in)    :: output
      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(inout) :: boundary
      type(s_inpcond),  intent(inout) :: ic 
      type(s_fwrk),     intent(inout) :: fwrk 

      ! I/O
      !
      character(len=MaxChar) :: fname

      ! Local
      !
      integer      :: nstate, block_size
      type(s_func) :: fs

      ! Dummy
      !
      integer :: is, iblock

      ! Array
      !
      type(s_infprop), allocatable :: ipse(:)


      if (.not. option%check_cumulative) return 

      write(iw,*)
      write(iw,'("Analyze> Start Cumulative analysis of steady-state properties")')

      ! Setup
      !
      nstate = option%nstate

      ! Allocate
      !
      allocate(ipse(option%ncum))

      do iblock = 1, option%ncum
      
        write(iw,'("Block ", i0)') iblock

        ! Set Start and End trajectories 
        !
        write(iw,'("  File block")')
        if (option%cumdirec == CumDirecIncrease) then
          do is = 1, nstate
            block_size  = ic%nfile_each_state(is) / option%ncum
            ic%ista(is) = 1
            ic%iend(is) = block_size * iblock
            write(iw,'(2x,"State ", i5, " : ", i0,2x,i0)') is, ic%ista(is), ic%iend(is) 
          end do
        else if (option%cumdirec == CumDirecDecrease) then
          do is = 1, nstate
            block_size  = ic%nfile_each_state(is) / option%ncum
            ic%ista(is) = block_size * (iblock - 1) + 1
            ic%iend(is) = ic%nfile_each_state(is)
            write(iw,'(2x,"State ", i5, " : ", i0,2x,i0)') is, ic%ista(is), ic%iend(is) 
          end do
        end if
        write(iw,*)

        ! Calculate Kernels
        !
        call calc_kernels(input, output, option, ic, boundary, fs, fwrk, &
                          set_boundary = .false., verbose = .false.)

        ! Setup Boundary conditions
        !
        call set_reflection(output, option, boundary, fs) 
        call set_product   (output, option, boundary, fs)

        ! Calculate steady-state properties
        !
        if (option%calc_Steady) then
          write(fname,'(a,".steady.",i4.4)') trim(output%fhead), iblock
        else if (option%calc_Pint) then
          write(fname,'(a,".int.",i4.4)')    trim(output%fhead), iblock
        end if
        call reacdyn_pint(output, option, boundary, fs, ipse(iblock), &
                          write_steady = .true., fname_out = fname)
      end do 

!
    end subroutine calc_cumulative
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine calc_senserr(input, output, option, boundary, ip, ic, fwrk) 
!-----------------------------------------------------------------------
      implicit none
     
      type(s_input),    intent(in)    :: input 
      type(s_output),   intent(in)    :: output
      type(s_option),   intent(in)    :: option
      type(s_boundary), intent(inout) :: boundary
      type(s_infprop),  intent(inout) :: ip
      type(s_inpcond),  intent(inout) :: ic 
      type(s_fwrk),     intent(inout) :: fwrk 

      ! I/O
      !
      integer                :: io
      character(len=MaxChar) :: fname

      ! Local
      !
      integer      :: nstate, nfepair
      type(s_func) :: fs

      ! Dummy
      !
      integer :: is, js, ks

      ! Array
      !
      type(s_infprop), allocatable :: ipse(:)
      real(8),         allocatable :: festdev(:)


      if (.not. option%check_senserr) return

      write(iw,*)
      write(iw,'("Analyze> Start Sensitivity analysis of steady-state properties")')

      ! Setup
      !
      nstate = option%nstate

      ! Allocate
      !
      allocate(ipse(nstate))

      do is = 1, nstate

        allocate(ipse(is)%prob(nstate))
        allocate(ipse(is)%fe(nstate))
        allocate(ipse(is)%fe_pair(nstate, nstate))

        ! Initialize
        !
        do js = 1, nstate
          ic%nfile_to_be_read(js) = ic%nfile_each_state(js)
        end do

        ! Set # of files to be read for specific state
        !
        if (option%is_reflect(is) .or. option%is_product(is)) cycle
        write(iw,'("Checking state ", i0)') is
        
        ic%nfile_to_be_read(is) = ic%nfile_to_be_read(is) * 0.5d0 

        ! Calculate Kernels
        !
        call calc_kernels(input, output, option, ic, boundary, fs, fwrk, &
                          set_boundary = .false., verbose = .false.)

        ! Setup Boundary conditions
        !
        call set_reflection(output, option, boundary, fs) 
        call set_product   (output, option, boundary, fs)

        ! Calculate steady-state properties
        !
        call reacdyn_pint(output, option, boundary, fs, ipse(is), &
                          write_steady = .false.)

        if (option%calc_Steady) then

          if (.not. allocated(festdev)) then
            allocate(festdev(nstate))
          end if

          festdev(is) = 0.0d0
          nfepair     = 0
          do js = 1, nstate - 1
            if (option%is_reflect(js) .or. option%is_product(js)) cycle 
            do ks = js + 1, nstate
              if (option%is_reflect(ks) .or. option%is_product(ks)) cycle
              nfepair = nfepair + 1
              festdev(is) = festdev(is) + (ipse(is)%fe_pair(ks, js) - ip%fe_pair(ks, js))**2
            end do
          end do
          festdev(is) = sqrt(festdev(is)/dble(nfepair))
        end if

      end do

      write(fname,'(a,".senserr")') trim(output%fhead)
      call open_file(fname, io)
      do is = 1, nstate
        if (option%is_reflect(is) .or. option%is_product(is)) cycle
        write(io,'(i5, f15.7)') is, festdev(is)
      end do
      close(io)
      write(iw,'(">> Done")') 
!
    end subroutine calc_senserr
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
            write(iw,'("State determination failed. &
                       &Please check state definition.")')
            stop
          end if

        end do
      end do 

    end subroutine get_state
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine read_f_unperturbed_id(option, nfile, unperturbed_ids, use_for_Rij)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      integer,        intent(in)    :: nfile 
      integer,        intent(inout) :: unperturbed_ids(nfile)
      integer,        intent(inout) :: use_for_Rij(nfile) 

      ! I/O
      !
      integer :: io

      ! Dummy
      !
      integer :: ifile


      ! Initialize 
      !
      unperturbed_ids = 0
      use_for_Rij     = 0

      ! Read
      !
      call open_file(option%f_unperturbed_id, io, stat = 'old')
      do ifile = 1, nfile
        read(io,*) unperturbed_ids(ifile), use_for_Rij(ifile) 
      end do
      close(io)
      
!
    end subroutine read_f_unperturbed_id 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine calc_kernels(input, output, option, ic, b, f, fwrk, &
                            set_boundary, verbose)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),     intent(in)    :: input
      type(s_output),    intent(in)    :: output
      type(s_option),    intent(in)    :: option
      type(s_inpcond),   intent(in)    :: ic
      type(s_boundary),  intent(inout) :: b
      type(s_func),      intent(inout) :: f
      type(s_fwrk),      intent(inout) :: fwrk
      logical, optional, intent(in)    :: set_boundary
      logical, optional, intent(in)    :: verbose

      ! I/O
      !
      integer                :: io
      character(len=MaxChar) :: fname 

      ! Local
      !
      character(len=MaxChar) :: line
      integer                :: ndim, nstep, nfile, nstate
      integer                :: nt_range, unp_id
      integer                :: nboundary
      logical                :: is_end
      logical                :: sb, vb
      logical                :: is_gen

      type(s_cv)    :: cv
      type(s_state) :: state

      ! Dummy
      !
      integer :: ifile, istep, iseg, is, js, is1, is2, ib, id, idir, ierr

      ! Arrays
      !
      integer, allocatable :: ncount_traj(:)

      ! Setup
      ! 
      ndim     = option%ndim * option%nmol 
      nfile    = input%ncv
      nstate   = option%nstate
      nt_range = option%nt_range

      vb = .false.
      if (present(verbose)) vb = verbose

      ! Allocate some arrays in f structure
      !
      if (.not. allocated(f%R)) then
        allocate(f%R   (0:nt_range, nstate, nstate))
        allocate(f%Rint(0:nt_range, nstate, nstate))
        allocate(f%P0  (0:nt_range, nstate))
      end if

      f%R    = 0.0d0
      f%Rint = 0.0d0
      f%P0   = 0.0d0

      ! Allocate work space
      !
      if (.not. allocated(fwrk%h)) then
        allocate(fwrk%h(nstate, nstate))
        allocate(fwrk%K(0:nt_range, fwrk%nkmax))
        allocate(fwrk%kmesh(nstate, nstate, nstate))
        fwrk%nk    = 0
        fwrk%K     = 0.0d0
        fwrk%kmesh = 0
      end if
      fwrk%K = 0.0d0
      fwrk%h = 0.0d0

      if (.not. allocated(ncount_traj)) then
        allocate(ncount_traj(nstate))
      end if
      ncount_traj = 0

      if (option%input_type == InputTypeTIMESERIES) then

        do ifile = 1, nfile

          ! For sensitivity analysis
          !    
          if (option%use_perturbed_traj .and. option%check_senserr) then
            id = ic%unperturbed_ids(ifile)
            if (id > 0) then
              ncount_traj(id) = ncount_traj(id) + 1
              if (ncount_traj(id) > ic%nfile_to_be_read(id)) then
                cycle
              end if
            end if
          end if

          ! For block average or cumulative analysis
          !    
          if (option%use_perturbed_traj .and. (option%check_blockave .or. option%check_cumulative)) then
            id = ic%unperturbed_ids(ifile)
            if (id > 0) then
              ncount_traj(id) = ncount_traj(id) + 1
              if (.not. option%is_errex(id)) then
                if (ncount_traj(id) < ic%ista(id) .or. ncount_traj(id) > ic%iend(id)) then
                  cycle
                end if
              end if
            end if
          end if

          if (vb) then
            write(iw,'("Analyze> Read CV file: ", 2x,a)') trim(input%fcv(ifile))
          end if
       
          iseg   = 0
          is_end = .false.
          call open_file(input%fcv(ifile), io, stat = 'old')
          do while (.not. is_end)
       
            ! Initialized
            ! 
            if (allocated(cv%x)) then
              call deallocate_cv(cv)
            end if
       
            if (allocated(state%data)) then
              state%nmol  = 0
              state%nstep = 0
              deallocate(state%data)
            end if
       
            ! Read
            !
            call read_cv_split(io, ndim, '#SPLIT', cv, is_end)
       
            if (cv%nstep /= 0) then
              iseg = iseg + 1
              call get_state(option, cv, state)
       
              state%unperturbed_id = -1
              state%use_for_Rij    = -1
              if (option%use_perturbed_traj) then
                state%unperturbed_id = ic%unperturbed_ids(ifile)
                state%use_for_Rij    = ic%use_for_Rij(ifile)
              end if
       
              ! Update Connectivity
              !
              call get_state_connectivity(output, option, state, b)
       
              ! Update R- and K-functions
              !
              if (option%use_rsto_Rij) then
                call update_rsto_Rij(option, state, f)
              else
                call update_Rij_wo_normalize(option, state, f)
              end if

              call update_Kijk_wo_normalize(option, state, fwrk)

            end if
       
          end do
       
          close(io)
       
        end do

        if (option%output_histogram) then
          is_gen = .false.
          call output_Kijk_hist(option, output, fwrk, is_gen)
          call output_Rij_hist (option, output, f, is_gen)
        end if

      else if (option%input_type == InputTypeHISTOGRAM) then

        do ifile = 1, nfile

          ! For sensitivity analysis
          !    
          if (option%use_perturbed_traj .and. option%check_senserr) then
            id = ic%unperturbed_ids(ifile)
            if (id > 0) then
              ncount_traj(id) = ncount_traj(id) + 1

              if (ncount_traj(id) > ic%nfile_to_be_read(id)) then
                cycle
              end if
            end if
          end if

          ! For block average or cumulative analysis
          !    
          if (option%use_perturbed_traj .and. (option%check_blockave .or. option%check_cumulative)) then
            id = ic%unperturbed_ids(ifile)
            if (id > 0) then
              ncount_traj(id) = ncount_traj(id) + 1
              if (.not. option%is_errex(id)) then
                if (ncount_traj(id) < ic%ista(id) .or. ncount_traj(id) > ic%iend(id)) then
                  cycle
                end if
              end if
            end if
          end if

          if (vb) then
            write(iw,'("Analyze> Read CV file: ", 2x,a)') trim(input%fcv(ifile))
          end if

          iseg   = 0
          is_end = .false.
          call open_file(input%fcv(ifile), io, stat = 'old')

          call seek_line(io, 'KIJK', ierr)
          if (ierr == 0) call update_Kijk_from_hist(io, option, fwrk)

          call seek_line(io, 'RIJ', ierr)
          if (ierr == 0) call update_Rij_from_hist(io, option, f) 
            
          close(io)
        end do
        call get_state_connectivity_from_h(option, fwrk%h, b)
      end if

      sb = .true.
      if (present(set_boundary)) then
        sb = set_boundary
      end if

      if (sb) then

        ! Show connectivity
        !
        if (vb) then
          write(iw,*)
          write(iw,'("Analyze> Get Connectivity")')
          call show_state_connectivity(option, b)
        end if
        
        ! Define boundary
        !
        if (vb) then
          write(iw,*)
          write(iw,'("Analyze> Define Boundary")')
        end if
        call define_boundary(option, b)
        nboundary = b%nboundary
        if (vb) write(iw,'(">> Done")')

      end if

      ! Convert arrays
      !
      call convert_Kijk_arrays(option, b, fwrk, f)

      ! Normalize R- and K-functions
      ! 
      call normalize_Rij(option, f)
      call normalize_Kijk(option, b, f, verbose = vb)

      ! Check Kijk
      !
      if (vb) then
        call check_Kijk(option, b, f)
      end if

      ! Compute R-integration and P0
      !
      if (vb) then
        write(iw,*)
        write(iw,'("Analyze> Calculate P0 from Rij")')
      end if
      call running_integral_Rij(option, f)
      call calc_P0_from_Rij    (option, b, f)
      if (vb) write(iw,'(">> Done")')

      ! Compute Mij
      !
      if (vb) then
        write(iw,*)
        write(iw,'("Analyze> Calculate Mjk from Kijk")')
      end if
      call calc_Mjk_from_Kijk(option, b, f)
      if (vb) write(iw,'(">> Done")')

    end subroutine calc_kernels
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_state_connectivity(output, option, state, boundary)
!-----------------------------------------------------------------------
      implicit none

      type(s_output),   intent(in)    :: output
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

      ! I/O
      !
      integer                :: io
      character(len=MaxChar) :: fname


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
    subroutine get_state_connectivity_from_h(option, h, boundary)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option 
      real(8),          intent(in)    :: h(option%nstate, option%nstate)
      type(s_boundary), intent(inout) :: boundary 

      ! Local
      !
      integer :: nstate

      ! Dummy
      !
      integer :: js, ks 

      ! I/O
      !
      integer                :: io
      character(len=MaxChar) :: fname


      ! Setup
      !
      nstate = option%nstate

      if (.not. allocated(boundary%is_connected)) then
        allocate(boundary%is_connected(nstate, nstate))
        boundary%is_connected = .false.
      end if

      do ks = 1, nstate
        do js = 1, nstate
          if (h(js, ks) > 0.999d0) then
            boundary%is_connected(js,  ks) = .true. 
            boundary%is_connected(ks,  js) = .true. 
          end if
        end do 
      end do

    end subroutine get_state_connectivity_from_h
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
      integer :: is, js, ib, ic 


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

      ! Setup influx_boundary
      !
      nboundary = boundary%nboundary

      allocate(boundary%n_influx_boundary(nstate))
      allocate(boundary%influx_boundary(2*nboundary, nstate))

      do is = 1, nstate
        ic = 0
        do ib = -nboundary, nboundary
          if (ib == 0) cycle
          js = boundary%b2p(2, ib)
          if (js == is) then
            ic = ic + 1
            boundary%influx_boundary(ic, is) = ib 
          end if
        end do
        boundary%n_influx_boundary(is) = ic
      end do

    end subroutine define_boundary 
!-----------------------------------------------------------------------

    include 'analyze_Rij.f90'
    include 'analyze_P0.f90'
    include 'analyze_Kijk.f90'
    include 'analyze_Mij.f90'
    include 'analyze_bc.f90'
    include 'analyze_reacdyn.f90'
    include 'analyze_pint.f90'

end module mod_analyze
!=======================================================================
