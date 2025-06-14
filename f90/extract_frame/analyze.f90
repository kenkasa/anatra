!=======================================================================
module mod_analyze
!=======================================================================
  use mod_util
  use mod_const
  use mod_cv
  use mod_input
  use mod_output
  use mod_ctrl
  use mod_netcdfio
  use mod_dcdio
  use mod_xtcio
  use mod_traj
  use xdr, only: xtcfile


  ! constants
  !

  ! structures
  !
  !type :: s_cv
  !  integer :: nstep
  !  integer :: ndim
  !  real(8), allocatable :: data(:,:)
  !end type s_cv

  type :: s_state
    integer :: nstep
    integer :: ndim
    logical :: is_reacted = .false.
    integer, allocatable :: data(:) 
    real(8), allocatable :: hist(:)
  end type s_state

  ! subroutines
  !
  public :: analyze
  public :: get_state
  public :: extract_dcd2dcd
  public :: extract_xtc2xtc

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in)   :: input
      type(s_output), intent(in)   :: output
      type(s_option), intent(in)   :: option

      ! I/O
      !
      integer                :: io
      character(len=MaxChar) :: fname

      ! Local
      !
      integer :: ndim
      integer :: ncv, ntraj, nweight
      integer :: next
      logical :: present_traj
      logical :: present_weight

      ! Dummy
      !
      integer :: ifile, istep, jstep, idim
      integer :: iunit_weight
      integer :: trajtype

      ! Arrays
      !
      type(s_cv),    allocatable :: cv(:)
      type(s_state), allocatable :: state(:)


      ! Setup
      !
      ndim    = option%ndim
      ncv     = input%ncv
      ntraj   = input%ntraj
      nweight = input%nweight

      allocate(cv(ncv), state(ncv)) 

      ! Get input/output trajectory type
      !
      if (trim(input%ftraj(1)) /= "") then
        present_traj = .true.
      else
        present_traj = .false.
      end if

      if (present_traj) then

        if (ncv /= ntraj) then
          write(iw,'("Analyze> Error.")')
          write(iw,'("# of traj files should be the same as # of CV files.")')
          stop
        end if

        call get_trajtype(input%ftraj(1), trajtype)

      end if

      ! Read CV files
      !
      write(iw,*)
      write(iw,'("Analyze> Read Time-Series Data file")')

      do ifile = 1, ncv
        call read_cv(input%fcv(ifile), ndim, cv(ifile))
      end do

      ! Read weight file if present
      !
      present_weight = .false.
      if (trim(input%fweight(1)) /= "") then

        present_weight = .true.

        if (ncv /= nweight) then
          write(iw,'("Analyze> Error.")')
          write(iw,'("# of weight files should be the same as # of CV files.")')
          stop
        end if

        write(iw,*)
        write(iw,'("Analyze> Read Weight file")')
        do ifile = 1, nweight
          call read_cv_weight(input%fweight(ifile), cv(ifile))
        end do

      end if

      ! Get state
      !
      write(fname,'(a,".cv")') trim(output%fhead)
      call open_file(fname, io)

      write(iw,*)
      write(iw,'("Analyze> Get State")')

      next = 0
      do ifile = 1, ncv
        call get_state(option%ndim,      &
                       option%state_def, &
                       cv(ifile),        &
                       state(ifile))

        do istep = 1, cv(ifile)%nstep
          if (state(ifile)%data(istep) == REACTIVE) then
            next = next + 1
            write(io,'(i10,2x)', advance='no') next 
            do idim = 1, option%ndim 
              write(io,'(e20.10)', advance='no') cv(ifile)%data(idim, istep)
            end do
            write(io,*)
          end if
        end do

      end do

      ! Extract weight 
      !
      if (present_weight) then

        write(fname, '(a,".weight")') trim(output%fhead) 
        call open_file(fname, io)

        jstep = 0
        do ifile = 1, ncv
          do istep = 1, cv(ifile)%nstep

            if (state(ifile)%data(istep) == REACTIVE) then
              jstep = jstep + 1
              write(io,'(i10,2x,e20.10)') jstep, cv(ifile)%weight(istep)
            end if

          end do
        end do

        close(io)

      end if

      !  Extract snapshots
      !
      if (present_traj) then

        if (trajtype == TrajTypeDCD) then
       
          call extract_dcd2dcd(input,          &
                               output,         &
                               option,         &
                               cv,             &
                               state,          &
                               next)
       
        else if (trajtype == TrajTypeXTC) then
       
          call extract_xtc2xtc(input,          &
                               output,         &
                               option,         &
                               cv,             &
                               state,          &
                               next)
       
        end if

      end if


    end subroutine analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_state(ndim, state_def, cv, state)
!-----------------------------------------------------------------------
      implicit none

      integer,                intent(in)  :: ndim
      real(8),                intent(in)  :: state_def(2, ndim, Nstate)
      type(s_cv),             intent(in)  :: cv
      type(s_state),          intent(out) :: state

      integer :: istep, istate, icv, ia
      integer :: nstep
      logical :: is_assigned
      real(8) :: wrk(ndim) 


      nstep       = cv%nstep

      ! for global
      allocate(state%data(nstep))

      do istep = 1, cv%nstep
        wrk(:) = cv%data(:, istep) 
        is_assigned = .false.
        do istate = 1, nstate
          if (is_assigned) then
            exit
          else
            ia = 0
            do icv = 1, ndim 
              if (wrk(icv) >= state_def(1, icv, istate) &
                .and. wrk(icv) < state_def(2, icv, istate)) then
                ia = ia + 1
              end if
            end do

            if (ia == ndim) then 
              is_assigned       = .true.
              state%data(istep) = StateInfo(istate)
            end if

          end if
        end do

        if (.not. is_assigned) then
          state%data(istep) = OTHERS 
        end if
      end do
      
      state%nstep = nstep
      state%ndim  = ndim

    end subroutine get_state 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine extract_dcd2dcd(input,          &
                               output,         &
                               option,         &
                               cv,             &
                               state,          &
                               next)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in) :: input
      type(s_output), intent(in) :: output
      type(s_option), intent(in) :: option
      type(s_cv),     intent(in) :: cv(:)
      type(s_state),  intent(in) :: state(:)
      integer,        intent(in) :: next

      type(s_dcd) :: dcdin, dcdout

      ! I/O 
      !
      integer                :: io_i, io_o
      character(len=MaxChar) :: fname

      ! Local
      !
      integer :: nfile, natm, nstep
      logical :: traj_reactive = .false.

      ! Dummy
      !
      integer :: itraj, istep, istep_tot


      ! Setup
      !
      nfile = input%ntraj

      !  Check the existence of reactive state in trajectory
      !
      if (next == 0) then
        write(iw,'("Extract_Dcd2dcd> No reactive frames.")')
        return
      end if

      ! Write
      !
      write(fname, '(a,".dcd")') trim(output%fhead)
      call dcd_open(fname, io_o)

      ! - Get header info
      !
      call dcd_open(input%ftraj(1), io_i)
      call dcd_read_header(io_i, dcdin)

      ! - Copy header info 
      !
      natm = dcdin%natm

      dcdout%natm       = natm
      dcdout%nstep      = next 
      dcdout%dcdinfo    = dcdin%dcdinfo
      dcdout%dcdinfo(1) = next
      dcdout%dcdinfo(4) = next

      call dcd_write_header(io_o, dcdout)

      ! - Allocate
      !
      call alloc_dcd(natm, 1, dcdin)
      call alloc_dcd(natm, 1, dcdout)

      call dcd_close(io_i, input%ftraj(1))

      ! - Copy coord.
      !
      istep_tot = 0
      do itraj = 1, nfile 

        ! Check
        !
        nstep         = state(itraj)%nstep 
        traj_reactive = .false.

        do istep = 1, nstep
          if (state(itraj)%data(istep) == REACTIVE) then
            traj_reactive = .true.
            exit
          end if 
        end do

        ! Skip if no reactive frames are present 
        !
        if (.not. traj_reactive) then
           write(iw,'("Skip reading ", a)') trim(input%ftraj(itraj)) 
           cycle
        end if 

        ! Copy
        !
        call dcd_open(input%ftraj(itraj), io_i)
        call dcd_read_header(io_i, dcdin)

        nstep = dcdin%nstep

        do istep = 1, nstep

          istep_tot = istep_tot + 1
          call read_dcd_oneframe(io_i, dcdin) 

          if (state(itraj)%data(istep) == REACTIVE) then
            write(iw,'("File ", i8, " Step ", i8, " : Extracted")') itraj, istep

            dcdout%box(1:3, 1)           = dcdin%box(1:3, 1)
            dcdout%coord(1:3, 1:natm, 1) = dcdin%coord(1:3, 1:natm, 1)

            call write_dcd_oneframe(io_o, 1, dcdout)

          end if

        end do

        call dcd_close(io_i, input%ftraj(itraj))

      end do

      call dcd_close(io_o, fname)

      ! Deallocate
      !
      call dealloc_dcd(dcdin)
      call dealloc_dcd(dcdout)


    end subroutine extract_dcd2dcd
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine extract_xtc2xtc(input,          &
                               output,         &
                               option,         &
                               cv,             & 
                               state,          &
                               next)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in) :: input
      type(s_output), intent(in) :: output
      type(s_option), intent(in) :: option
      type(s_cv),     intent(in) :: cv(:)
      type(s_state),  intent(in) :: state(:)
      integer,        intent(in) :: next

      type(xtcfile) :: xtcin, xtcout

      ! I/O
      !
      integer                :: io
      character(len=MaxChar) :: fname

      ! Local
      !
      integer :: natm, nfile
      logical :: is_end

      ! Dummy
      !
      integer :: itraj, istep, istep_tot


      ! Setup Input XTC / Output XTC
      !
      write(fname,'(a,".xtc")') trim(output%fhead)

      call xtcin%init(trim(input%ftraj(1)))
      natm = xtcin%natoms
      call xtcin%close

      call xtcout%init(trim(fname), 'w')

      ! Start extract
      !
      istep_tot = 0
      do itraj = 1, input%ntraj

        call xtcin%init(trim(input%ftraj(itraj)))
       
        is_end = .false.
        istep  = 0
        do while (.not. is_end)
          call xtcin%read

          if (xtcin%STAT /= 0) then
            is_end = .true.
            exit
          else
            if (state(itraj)%data(istep + 1) == REACTIVE) then
              write(iw,'(" Step ",i8, ": Extracted")') istep_tot + 1

              call xtcout%write(xtcin%natoms,   &
                                istep_tot + 1,  &
                                xtcin%time,     &
                                xtcin%box,      &
                                xtcin%pos,      &
                                xtcin%prec)

            end if

            istep     = istep     + 1
            istep_tot = istep_tot + 1 

          end if

        end do

        call xtcin%close

      end do

      ! Deallocate memory
      !
      call xtcout%close


    end subroutine extract_xtc2xtc
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
