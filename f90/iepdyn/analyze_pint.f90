!-----------------------------------------------------------------------
    subroutine reacdyn_pint(output, option, boundary, f, ip, &
                            write_steady, fname_out)
!-----------------------------------------------------------------------
      implicit none

      type(s_output),    intent(in)    :: output
      type(s_option),    intent(in)    :: option
      type(s_boundary),  intent(in)    :: boundary
      type(s_func),      intent(inout) :: f
      type(s_infprop),   intent(inout) :: ip
      logical, optional, intent(in)    :: write_steady
      character(len=MaxChar), optional, intent(in) :: fname_out

      ! I/O
      !
      integer :: io

      ! Local
      !
      integer                :: nstate, nt_range, nboundary, nbt
      real(8)                :: dt, kT
      character(len=MaxChar) :: fname
      logical                :: ws

      integer                :: thread_id

      ! Dummy
      !
      integer :: is, js, ks, is1, is2, js1, js2, iref
      integer :: iu, ju, ku
      integer :: ib, jb, inflx, jnflx, id
      integer :: lda, ldb, ldvl, ldvr, info, lwork
      real(8) :: psum
      real(8) :: t1, t2

      ! Arrays
      !
      real(8), allocatable :: Qint(:), Qinf(:)
      real(8), allocatable :: Kuu(:, :), Mu(:), Ru(:)
      real(8), allocatable :: Xuu(:, :), Pint(:), Pss(:)
      integer, allocatable :: map(:, :), mapb(:), mapb_inv(:)
      integer, allocatable :: ipiv(:)
      real(8), allocatable :: wr(:), wi(:)        ! for eigen
      real(8), allocatable :: vl(:, :), vr(:, :)  ! for eigen
      real(8), allocatable :: work(:)             ! for eigen 


      if (.not. (option%calc_Pint .or. option%calc_Steady)) return 

      ! Setup
      !
      nstate    = option%nstate
      nt_range  = option%nt_range
      dt        = option%dt_out
      nboundary = boundary%nboundary
      nbt       = nboundary * 2
      kT        = option%temperature * Boltz

      if (option%calc_Steady .and. .not. allocated(ip%prob)) then
        allocate(ip%prob(nstate))
        allocate(ip%fe(nstate))
        allocate(ip%fe_pair(nstate, nstate))
      end if

      if (.not. allocated(Mu)) then
        allocate(Mu(nbt))
        allocate(Kuu(nbt, nbt))
        allocate(Ru(nbt))
        allocate(Xuu(nbt, nbt))
        allocate(Pint(nstate))
        allocate(Pss(nstate))
        allocate(Qint(nbt), Qinf(nbt))
        allocate(ipiv(nbt)) 
        allocate(map(2, nbt))
        allocate(mapb(nbt))
        allocate(mapb_inv(-nboundary:nboundary))
        allocate(wr(nbt), wi(nbt), vl(nbt, nbt), vr(nbt, nbt))
      end if
      Mu       = 0.0d0
      Kuu      = 0.0d0
      Ru       = 0.0d0
      Xuu      = 0.0d0
      Pint     = 0.0d0
      Pss      = 0.0d0
      Qint     = 0.0d0
      map      = 0
      mapb     = 0
      mapb_inv = 0

      iu  = 0
      do ib = -nboundary, nboundary
        if (ib == 0) cycle
        iu  = iu + 1
        is1 = boundary%b2P(1, ib)
        is2 = boundary%b2p(2, ib)
        map(1, iu)          = is1
        map(2, iu)          = is2
        mapb(iu)            = ib
        mapb_inv(ib)        = iu

        Mu(iu) = sum(f%M(0:nt_range, ib))
        Ru(iu) = sum(f%R(0:nt_range, is2, is1))

        ju = 0
        do jb = -nboundary, nboundary
          if (jb == 0)  cycle
          ju  = ju + 1
          js1 = boundary%b2p(1, jb)
          js2 = boundary%b2p(2, jb)

          if (is1 == js2) then
            Kuu(iu, ju) = sum(f%K(0:nt_range, is2, jb))
          end if 
        end do
      end do

      ! Reflecting boundary
      !
      if (option%use_reflection_state) then
        do iu = 1, nbt
          ib = mapb(iu)
          is1 = boundary%b2p(1, ib)
          is2 = boundary%b2p(2, ib)
          if (boundary%conv_direc(ib)) then
            jb = -ib
            ju = mapb_inv(jb)
            Ru(ju) = sum(f%R(0:nt_range, is2, is1))
            Ru(iu) = 0.0d0

            Kuu(ju, :) = 0.0d0
            do ku = 1, nbt
              Kuu(ju, ku) = Kuu(iu, ku)
              Kuu(iu, ku) = 0.0d0
            end do 
          end if
        end do
      end if

      if (option%calc_Pint) then
        do ju = 1, nbt
          Xuu(ju, ju) = 1.0d0
          do iu = 1, nbt
            Xuu(iu, ju) = Xuu(iu, ju) - Kuu(iu, ju) * dt
          end do
        end do
        
        do iu = 1, nbt
          Qint(iu) = Ru(iu) * dt
        end do
        lda  = nbt
        ldb  = nbt
        ipiv = 0
        call dgesv(nbt, 1, Xuu, lda, ipiv, Qint, ldb, info) 
        
        do is = 1, nstate 
          if (option%is_dissoc(is)) then
            Pint(is) = 0.0d0
            cycle
          end if
        
          Pint(is) = sum(f%P0(0:nt_range, is)) * dt
          do iu = 1, nbt
            if (map(2, iu) /= is) cycle
            Pint(is) = Pint(is) + Mu(iu) * dt * Qint(iu)  
          end do
        end do
        
        ! Integrated values 
        !
        psum = 0.0
        do is = 1, nstate
          if (option%is_initial(is)) then
            psum = psum + Pint(is)
          end if
        end do

        ws = .true.
        if (present(write_steady)) then
          ws = write_steady
        end if

        if (ws) then
          write(fname,'(a,".int")') trim(output%fhead)
          if (present(fname_out)) then
            write(fname,'(a)') trim(fname_out)
          end if

          call open_file(fname, io)
          write(io,'("Pint InitTotal  ", e15.7)') psum
          do is = 1, nstate
            if (option%is_dissoc(is)) cycle
            write(io,'("Pint ", i9, 2x, e15.7)') is, Pint(is)
          end do
          close(io)
        end if

      end if

      if (option%calc_Steady) then

        ! Construct Transition kernel
        !
        Xuu = 0.0d0
        do ju = 1, nbt
          do iu = 1, nbt
            Xuu(iu, ju) = Kuu(iu, ju) * dt
          end do
        end do

        ! Solve Eigen-value problem of Xuu
        !
        lda   = nbt
        ldvl  = nbt
        ldvr  = nbt
        lwork = -1
        allocate(work(1))
        call dgeev('N', 'V', nbt, Xuu, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info) 
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))

        call dgeev('N', 'V', nbt, Xuu, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
        deallocate(work)

        if (info /= 0) then
          write(iw,'("Reacdyn_Pint> Error.")')
          write(iw,'("Failed in solving eigen equation of Kint. stop.")')
          stop
        end if

!        ! << DEBUG
!        write(iw,*)
!        write(iw,'("[ Eigenvalues ]")')
!        do iu = 1, nbt
!          write(iw,'(i5,2x,f20.10, f20.10)') iu, wr(iu), wi(iu) 
!        end do
!        write(iw,*)
!        write(iw,'("[ Eigenvector corresponding to steady state]")')
!        iu = maxloc(wr(:), dim = 1)
!
!        write(iw,'("Eigenvalue: ", f20.10)') wr(iu)
!        do ju = 1, nbt
!          write(iw,'(i5,2x,f20.10)') ju, vr(ju, iu)
!        end do 
!        ! >> DEBUG

        iu = maxloc(wr(:), dim = 1)
        Qinf(:) = vr(:, iu)

        Pss = 0.0d0
        do is = 1, nstate 
          do ju = 1, nbt
            if (map(2, ju) /= is) cycle 
            Pss(is) = Pss(is) + Mu(ju) * dt * Qinf(ju) 
          end do
        end do

!        ! << DEBUG
!        write(iw,*)
!        write(iw,'("[ Steady state population ]")')
!        do is = 1, nstate
!          write(iw,'(i5,2x,f20.10)') is, Pss(is)
!        end do
!        ! >> DEBUG

        psum = 0.0
        do is = 1, nstate
          if (option%is_initial(is)) then
            psum = psum + Pss(is)
          end if
        end do

        Pss(:) = Pss(:) / psum
        do is = 1, nstate
          ip%prob(is) = Pss(is)
          ip%fe(is)   = -kT * log(ip%prob(is))
        end do

        ip%fe_pair = 0.0d0
        do is = 1, nstate - 1
          if (option%is_reflect(is) .or. option%is_product(is)) cycle
          do js = is + 1, nstate
            if (option%is_reflect(js) .or. option%is_product(js)) cycle
            ip%fe_pair(js, is) = ip%fe(js) - ip%fe(is)
          end do
        end do

        ! Write steady-state properties
        !
        ws = .true.
        if (present(write_steady)) then
          ws = write_steady
        end if

        if (ws) then
          write(fname,'(a,".steady")') trim(output%fhead)
          if (present(fname_out)) then
            write(fname,'(a)') trim(fname_out)
          end if

          call open_file(fname, io)
          do is = 1, nstate
            if (option%is_reflect(is) .or. option%is_product(is)) cycle
            write(io,'(i5, 2x, 2(e15.7,2x))') is, ip%prob(is), ip%fe(is)
          end do
          close(io)
        end if

      end if

      ! Deallocate memory
      !
      
    end subroutine reacdyn_pint
!-----------------------------------------------------------------------
