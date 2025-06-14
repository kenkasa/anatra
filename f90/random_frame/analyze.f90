!=======================================================================
module mod_analyze
!=======================================================================
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_ctrl
  use mod_traj
  use mod_random

  ! constants
  !

  ! structures
  !

  ! subroutines
  !
  public  :: analyze
  public  :: analyze_shuffle
  private :: heapsort_integer
  private :: dcd2dcd
  private :: dcd2xtc
  private :: xtc2xtc
  private :: dcd2dcd_shuffle
  private :: dcd2xtc_shuffle
  private :: xtc2xtc_shuffle

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in)   :: input
      type(s_output), intent(in)   :: output
      type(s_option), intent(in)   :: option

      ! Local
      !
      integer :: nstep, natm, next
      integer :: iseed
      integer :: trajtype_in, trajtype_out

      ! Dummy
      !
      integer :: i, id, istep, it

      ! Arrays
      !
      integer, allocatable :: rand(:), used_snap(:) 


      write(iw,'("Analyze> Start the extraction")')

      ! Get trajectory type
      !
      call get_trajtype(input%ftraj(1), trajtype_in)

      trajtype_out = option%output_trajtype

      ! Setup 
      !
      !nstep = dcdin%nstep

      next  = option%nsample
      iseed = option%iseed

      if (trajtype_in == TrajTypeDCD) then
        call get_total_step_from_dcd(input%ftraj, nstep)
      else if (trajtype_in == TrajTypeXTC) then
        call get_total_step_from_xtc(input%ftraj, nstep)
      end if

      ! memory allocation
      !
      allocate(rand(next), used_snap(nstep))

      ! Generate random seed
      !
      call get_seed(iseed)
      call initialize_random(iseed)

      ! Generate random numbers
      !
      if (.not. option%use_allsnap) then
        call get_random_integer(option%nsample, 1, nstep, option%duplicate, rand)
        call heapsort_integer(rand)

      ! Extract all snapshots
      !
      else
        do istep = 1, next
          rand(istep) = istep 
        end do
      end if

      write(iw,*)
      write(iw,'("Analyze> Selected snapshots")')
      do i = 1, next 
        write(iw,'(i8)', advance="no") rand(i)
        if (mod(i,10) == 0) &
          write(iw,*)
      end do

      used_snap = 0
      do istep = 1, next
        id            = rand(istep)
        used_snap(id) = used_snap(id) + 1 
      end do

      if (trajtype_in == TrajTypeDCD) then

        if (trajtype_out == TrajTypeDCD) then
          call dcd2dcd(input, output, option, rand, used_snap) 
        else if (trajtype_out == TrajTypeXTC) then
          call dcd2xtc(input, output, option, rand, used_snap)
        end if

      else if (trajtype_in == TrajTypeXTC) then

        if (trajtype_out == TrajTypeDCD) then
          !call xtc2dcd(input, output, option, rand, used_snap)
          write(iw,'("Analyze> Error.")')
          write(iw,'("Sorry, XTC => DCD convert is not supported.")')
        else if (trajtype_out == TrajTypeXTC) then
          call xtc2xtc(input, output, option, rand, used_snap)
        end if

      end if

      deallocate(rand, used_snap)

    end subroutine analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine analyze_shuffle(input, output, option)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in)   :: input
      type(s_output), intent(in)   :: output
      type(s_option), intent(in)   :: option

      ! Local
      !
      integer :: nstep, natm, next, iseed
      integer :: trajtype_in, trajtype_out

      ! Dummy
      !
      integer :: i, id, istep, it

      ! Arrays
      !
      integer, allocatable :: frame_index(:)


      write(iw,'("Analyze> Start the extraction")')

      ! Get trajectory type
      !
      call get_trajtype(input%ftraj(1), trajtype_in)

      trajtype_out = option%output_trajtype

      ! Setup 
      !
      iseed = option%iseed

      if (trajtype_in == TrajTypeDCD) then
        call get_total_step_from_dcd(input%ftraj, nstep)
      else if (trajtype_in == TrajTypeXTC) then
        call get_total_step_from_xtc(input%ftraj, nstep)
      end if

      allocate(frame_index(nstep))

      ! Generate random seed
      !
      call get_seed(iseed)
      call initialize_random(iseed)

      ! Fisher-Yates shuffle
      !
      write(iw,*)
      write(iw,'("Analyze_Shuffle> Get shuffle index")')

      do istep = 1, nstep
        frame_index(istep) = istep
      end do

      call shuffle_fisher_yates(frame_index)

      write(iw,*)
      write(iw,'("Analyze_shuffle> Shuffle index")')
      do i = 1, nstep
        write(iw,'(i8)', advance="no") frame_index(i)
        if (mod(i,10) == 0) &
          write(iw,*)
      end do

      ! Create shuffled trajectories
      !
      if (trajtype_in == TrajTypeDCD) then

        if (trajtype_out == TrajTypeDCD) then
          call dcd2dcd_shuffle(input, output, option, frame_index) 
        else if (trajtype_out == TrajTypeXTC) then
          write(iw,'("Analyze> Error.")')
          write(iw,'("Sorry, DCD => XTC convert is not supported for shuffle.")')
          stop
          !call dcd2xtc_shuffle(input, output, option, frame_index)
        end if

      else if (trajtype_in == TrajTypeXTC) then

        if (trajtype_out == TrajTypeDCD) then
          !call xtc2dcd(input, output, option, rand, used_snap)
          write(iw,'("Analyze> Error.")')
          write(iw,'("Sorry, XTC => DCD convert is not supported for shuffle.")')
          stop
        else if (trajtype_out == TrajTypeXTC) then
          write(iw,'("Analyze> Error.")')
          write(iw,'("Sorry, XTC => DCD convert is not supported for shuffle.")')
          stop
          !call xtc2xtc_shuffle(input, output, option, rand, used_snap)
        end if

      end if

      deallocate(frame_index)

    end subroutine analyze_shuffle
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
    subroutine dcd2dcd(input, output, option, rand, used_snap)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),   intent(in) :: input
      type(s_output),  intent(in) :: output
      type(s_option),  intent(in) :: option
      integer,         intent(in) :: rand(:)
      integer,         intent(in) :: used_snap(:)

      type(s_dcd) :: dcd_in, dcd_out

      ! I/O
      !
      integer                :: io_i, io_o
      character(len=MaxChar) :: fname
      character(len=MaxChar) :: finpcrd

      ! Local
      !
      integer                :: nstep_tot
      integer                :: natm, nrand, nstep, ndegen

      ! Dummy
      !
      integer :: id, jd, itraj, istep, jstep, istep_tot
      integer :: iwrite
      logical :: is_end


      write(iw,*)
      write(iw,'("Dcd2dcd> Start the sampling")')

      nrand = size(rand)

      call get_natm_from_dcd(input%ftraj(1), natm) 

      write(fname, '(a, ".dcd")') trim(output%fhead)
      call dcd_open(fname, io_o)

      call alloc_dcd(natm, 1, dcd_out)

      istep_tot = 0
      iwrite    = 0
      do itraj = 1, input%ntraj
        call dcd_open(input%ftraj(itraj), io_i)
        call dcd_read_header(io_i, dcd_in)

        nstep = dcd_in%nstep

        if (itraj == 1) then
          call alloc_dcd(dcd_in%natm, 1, dcd_in)

          dcd_out%natm       = natm
          dcd_out%dcdinfo    = dcd_in%dcdinfo
          dcd_out%dcdinfo(1) = nrand

          call dcd_write_header(io_o, dcd_out)
        end if

        do istep = 1, nstep
          istep_tot = istep_tot + 1
          
          if (mod(istep_tot, 100) == 0) then
            write(iw,'("Step ",i0)') istep_tot
          end if

          call read_dcd_oneframe(io_i, dcd_in)

          ndegen = used_snap(istep_tot)

          if (ndegen >= 1) then
            dcd_out%coord = dcd_in%coord
            dcd_out%box   = dcd_in%box
            do jstep = 0, ndegen - 1
              iwrite = iwrite + 1

              call write_dcd_oneframe(io_o, 1, dcd_out)

              if (option%out_rst7) then
                write(finpcrd,'(a,i5.5,".inpcrd")') trim(output%fhead), iwrite
                call write_inpcrd(finpcrd, dcd_out%coord(1:3, 1:natm, 1), dcd_out%box(1:3, 1))
              end if
            end do
          end if

        end do

        call dcd_close(io_i, input%ftraj(itraj))

      end do

      call dcd_close(io_o, fname) 

!
    end subroutine dcd2dcd
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine dcd2xtc(input, output, option, rand, used_snap)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),   intent(in) :: input
      type(s_output),  intent(in) :: output
      type(s_option),  intent(in) :: option
      integer,         intent(in) :: rand(:)
      integer,         intent(in) :: used_snap(:)

      real(8), parameter :: ang2nm = 0.1d0

      type(s_dcd)   :: dcd_in
      type(xtcfile) :: xtc_out
  
      ! I/O
      ! 
      integer                :: io_i
      character(len=MaxChar) :: fname
      character(len=MaxChar) :: finpcrd

      ! Local
      !
      integer :: natm, nrand, nstep, nstep_tot
      integer :: ndegen
      real(8) :: box(3, 3)
      integer :: trajtype

      ! Dummy
      !
      integer :: id, jd, iwrite, itraj
      integer :: istep, jstep, istep_tot

      ! Arrays
      !
      real(4), allocatable :: coord(:, :)


      write(iw,*)
      write(iw,'("Dcd2xtc> Start the sampling")')

      nrand = size(rand)
    
      write(fname,'(a,".xtc")') trim(output%fhead) 
      call xtc_out%init(fname, 'w')

      istep_tot = 0
      iwrite    = 0
      do itraj = 1, input%ntraj
        call dcd_open(input%ftraj(itraj), io_i)
        call dcd_read_header(io_i, dcd_in)
        nstep = dcd_in%nstep

        if (itraj == 1) then
          natm = dcd_in%natm
          call alloc_dcd(dcd_in%natm, 1, dcd_in)

          allocate(coord(1:3, 1:natm))

        end if

        do istep = 1, nstep
           istep_tot = istep_tot + 1

            if (mod(istep_tot, 100) == 0) then
              write(iw,'("Step ",i0)') istep_tot
            end if

            call read_dcd_oneframe(io_i, dcd_in)

            ndegen = used_snap(istep_tot)

            if (ndegen >= 1) then
              box = 0.0d0
              box(1, 1) = dcd_in%box(1, 1) * ang2nm
              box(2, 2) = dcd_in%box(2, 1) * ang2nm
              box(3, 3) = dcd_in%box(3, 1) * ang2nm

              coord(1:3, 1:natm) = dcd_in%coord(1:3, 1:natm, 1) * ang2nm

              do jstep = 0, ndegen - 1
                iwrite = iwrite + 1

                call xtc_out%write(natm,               &
                                   0,                  &
                                   0.0,                &
                                   real(box),          &
                                   coord(1:3, 1:natm), &
                                   real(1000.0d0))

                if (option%out_rst7) then
                  write(finpcrd,'(a,i5.5,".inpcrd")') trim(output%fhead), iwrite
                  call write_inpcrd(finpcrd,                      &
                                    dcd_in%coord(1:3, 1:natm, 1), & 
                                    dcd_in%box(1:3, 1))
                end if
              end do
            end if
        end do

        call dcd_close(io_i, input%ftraj(itraj))

      end do

      call xtc_out%close
      deallocate(coord)


    end subroutine dcd2xtc
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine xtc2xtc(input, output, option, rand, used_snap) 
!-----------------------------------------------------------------------
      implicit none

      type(s_input),   intent(in) :: input
      type(s_output),  intent(in) :: output
      type(s_option),  intent(in) :: option
      integer,         intent(in) :: rand(:)
      integer,         intent(in) :: used_snap(:)

      type(xtcfile) :: xtc_in, xtc_out

      ! I/O
      !
      integer                :: io_o
      character(len=MaxChar) :: fname 
      character(len=MaxChar) :: finpcrd

      ! Local
      !
      integer :: nstep_tot
      integer :: natm, nrand, nstep, ndegen
      real(8) :: box(3)
      logical :: is_end

      ! Dummy
      !
      integer :: id, jd, itraj, istep, jstep, istep_tot
      integer :: iwrite

      ! Arrays
      !
      real(8), allocatable   :: coord(:, :)


      write(iw,*)
      write(iw,'("Xtc2xtc> Start the sampling")')

      nrand = size(rand)

      write(fname, '(a,".xtc")') trim(output%fhead)
      call xtc_out%init(fname, 'w')

      istep_tot = 0
      iwrite    = 0
      do itraj = 1, input%ntraj
        call xtc_in%init(input%ftraj(itraj), 'r')

        if (itraj == 1) then
          natm = xtc_in%natoms
          allocate(coord(1:3, xtc_in%natoms))
        end if

        is_end = .false.
        do while (.not. is_end)
          call xtc_in%read

          if (xtc_in%STAT /= 0) then
            is_end = .true.
            exit
          end if

          istep_tot = istep_tot + 1

          if (mod(istep_tot, 100) == 0) then
            write(iw,'("Step ",i0)') istep_tot
          end if

          ndegen = used_snap(istep_tot)

          if (ndegen >= 1) then

            coord(1:3, 1:natm) = xtc_in%pos(1:3, 1:natm) * 10.0d0
            box(1)             = xtc_in%box(1, 1) * 10.0d0
            box(2)             = xtc_in%box(2, 2) * 10.0d0
            box(3)             = xtc_in%box(3, 3) * 10.0d0

            do jstep = 0, ndegen - 1
              iwrite = iwrite + 1
              call xtc_out%write(xtc_in%natoms,  &
                                 istep_tot,      &
                                 xtc_in%time,    &
                                 xtc_in%box,     &
                                 xtc_in%pos,     &
                                 xtc_in%prec)

              if (option%out_rst7) then
                write(finpcrd,'(a,i5.5,".inpcrd")') trim(output%fhead), iwrite
                call write_inpcrd(finpcrd, coord, box)
              end if
            end do
          end if

        end do

        call xtc_in%close
      end do

      call xtc_out%close

      deallocate(coord)


    end subroutine xtc2xtc
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine dcd2dcd_shuffle(input, output, option, frame_index)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),   intent(in) :: input
      type(s_output),  intent(in) :: output
      type(s_option),  intent(in) :: option
      integer,         intent(in) :: frame_index(:) 

      type(s_dcd) :: dcd_in, dcd_out
      type(s_dcd) :: dcd_tmp

      ! I/O
      !
      integer                :: io_i, io_o
      character(len=MaxChar) :: fname
      character(len=MaxChar) :: finpcrd

      ! Local
      !
      integer :: nstep_tot
      integer :: natm, nrand, nstep, ndegen
      logical :: is_end

      ! Dummy
      !
      integer :: id, jd, itraj, istep, jstep, istep_tot
      integer :: iwrite

      ! Arrays
      !
      real(8), allocatable   :: coord_store(:, :, :)
      real(8), allocatable   :: box_store(:, :)


      write(iw,*)
      write(iw,'("Dcd2dcd_Shuffle> Output shuffled trajectory")')

      ! Total # of steps
      !
      nstep_tot = size(frame_index)

      ! # of atoms
      !
      call get_natm_from_dcd(input%ftraj(1), natm)

      ! Allocation of output dcd
      !
      write(fname, '(a,".dcd")') trim(output%fhead)
      call dcd_open(fname, io_o)
      call alloc_dcd(natm, 1, dcd_out)

      ! Allocation of working space
      !
      allocate(coord_store(1:3, natm, nstep_tot))
      allocate(box_store(1:3, nstep_tot))

      ! Read Input dcd
      !
      istep_tot = 0
      do itraj = 1, input%ntraj
        write(iw,'("Read Traj: ", a)') trim(input%ftraj(itraj))

        call dcd_open(input%ftraj(itraj), io_i)
        call dcd_read_header(io_i, dcd_in)
        nstep = dcd_in%nstep

        if (itraj == 1) then
          call alloc_dcd(dcd_in%natm, 1, dcd_in)
          dcd_out%natm       = natm
          dcd_out%dcdinfo    = dcd_in%dcdinfo
          dcd_out%dcdinfo(1) = nstep_tot
          call dcd_write_header(io_o, dcd_out)
        end if

        do istep = 1, nstep
          istep_tot = istep_tot + 1

          if (mod(istep_tot, 100) == 0) then
            write(iw,'("Read step ", i0)') istep_tot
          end if

          call read_dcd_oneframe(io_i, dcd_in)

          coord_store(1:3, 1:natm, istep_tot) = dcd_in%coord(1:3, 1:natm, 1)
          box_store(1:3, istep_tot)           = dcd_in%box(1:3, 1)
        end do 

        call dcd_close(io_i, input%ftraj(itraj))

      end do

      ! Shuffle
      !
      do istep = 1, nstep_tot
        id = frame_index(istep)

        dcd_out%coord(1:3, 1:natm, 1) = coord_store(1:3, 1:natm, id)
        dcd_out%box(1:3, 1)           = box_store(1:3, id)

        call write_dcd_oneframe(io_o, 1, dcd_out)

      end do 

      call dcd_close(io_o, fname)

      call dealloc_dcd(dcd_in)
      call dealloc_dcd(dcd_out)
      deallocate(coord_store, box_store)
!
    end subroutine dcd2dcd_shuffle
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine write_inpcrd(fname, coord, box)
!-----------------------------------------------------------------------
      implicit none

      character(len=MaxChar), intent(in) :: fname
      real(8),                intent(in) :: coord(:, :)
      real(8),                intent(in) :: box(3)

      integer                :: iatm, io
      integer                :: natm
      character(len=MaxChar) :: title


      title = "TITLE : generated by ANATRA"


      call open_file(fname, io)

      natm = size(coord(1, :))

      write(io,'(a)')   trim(title)
      write(io,'(i5)')  natm
      do iatm = 1, natm
        write(io,'(3f12.7)', advance="no") &
          coord(1:3, iatm)

        if (mod(iatm, 2) == 0) &
          write(io,*)

      end do

      write(io,'(6f12.7)') &
        box(1:3), 90.0d0, 90.0d0, 90.0d0
      close(io)

    end subroutine write_inpcrd
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine heapsort_integer(arr)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(inout) :: arr(:)

      integer :: i, j, k, l
      integer :: tmp
      integer :: n


      n = size(arr)

      if (n <= 1) &
        return

      l = n / 2 + 1
      k = n

      do while (k /= 1)

        if (l > 1) then
          l   = l - 1
          tmp = arr(l) 
        else
          tmp      = arr(k)
          arr(k) = arr(l)
          k        = k - 1

          if (k == 1) then
            arr(l) = tmp
            exit
          end if

        end if

        i = l
        j = l + l

        do while (j <= k)
          if (j < k) then
            if (arr(j) < arr(j+1)) &
              j = j + 1
          end if

          if (tmp < arr(j)) then
            arr(i) = arr(j)
            i      = j
            j      = j + j
          else
            j      = k + 1
          end if
        end do

        arr(i) = tmp

      end do
      

    end subroutine heapsort_integer
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine shuffle_fisher_yates(ind)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(inout) :: ind(:)

      integer :: narr

      integer :: iarr, jarr, rand(1), tmp


      narr = size(ind)

      do iarr = 1, narr - 1
        jarr = narr - iarr + 1
        call get_random_integer(1, 1, jarr, .true., rand(1))

        tmp          = ind(jarr)
        ind(jarr)    = ind(rand(1))
        ind(rand(1)) = tmp 
      end do 

    end subroutine shuffle_fisher_yates
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
