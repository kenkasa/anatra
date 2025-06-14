!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_ctrl
  use mod_traj
  use mod_anaparm
  use mod_prmtop
  use mod_xtcio
  use mod_dcdio
  use mod_com
  use mod_potential
  !use mod_pme_str
  !use mod_pme
  !use mod_fftmkl

  ! constants
  !

  ! subroutines
  !
  public :: analyze

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option, traj)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in)    :: input
      type(s_output), intent(in)    :: output
      type(s_option), intent(in)    :: option
      type(s_traj),   intent(inout) :: traj(2)

      type(s_dcd)            :: dcd
      type(xtcfile)          :: xtc
      type(s_com)            :: com(2)
      type(s_prmtop)         :: prmtop
      type(s_prmtop)         :: prmtop2
      type(s_anaparm)        :: anaparm
      type(s_pot)            :: pot
      type(s_pot)            :: pot_dual
      !type(s_fftinfo)        :: fftinfo

      ! for Ewald
      !type(s_pmevar)         :: pmevar

      integer                :: istep, istep_tot, itraj, imol, jmol
      integer                :: iunit
      integer                :: iunit_ij, iunit_fsq, iunit_ijdual
      integer                :: iunit_fsqdual, iunit_ijdiff, iunit_fsqdiff
      integer                :: natm, nmol(2)
      integer                :: trajtype
      real(8)                :: rljcut2
      character(len=MaxChar) :: fout_ij,  fout_ij_dual,  fout_ij_diff
      character(len=MaxChar) :: fout_fsq, fout_fsq_dual, fout_fsq_diff
      logical                :: is_end

      real(8), allocatable :: uij(:, :)
      real(8), allocatable :: uLJ(:, :), uES(:, :), uEL(:, :)
      real(8), allocatable :: uCORR(:, :)
      real(8), allocatable :: frc(:, :, :)

      ! for dual calculation
      real(8), allocatable :: uij_dual  (:, :)
      real(8), allocatable :: uLJ_dual  (:, :)
      real(8), allocatable :: uES_dual  (:, :) 
      real(8), allocatable :: uEL_dual  (:, :)
      real(8), allocatable :: uCORR_dual(:, :)
      real(8), allocatable :: frc_dual  (:, :, :)

      real(8), allocatable :: uij_diff  (:, :)
      real(8), allocatable :: frc_diff  (:, :, :)


      ! Read parameter file
      !
      if (option%parmformat == ParmFormatANAPARM) then
        call read_anaparm(input%fanaparm, anaparm) 
      else if (option%parmformat == ParmFormatPRMTOP) then
        call read_prmtop(input%fprmtop, prmtop)

        if (option%pme_dual) then
          call read_prmtop(input%fprmtop2, prmtop2)
        end if 

      end if

      ! Setup
      !
      call get_com(option%mode(1),       &
                   traj(1),              &
                   com(1),               &
                   setup = .true.,       &
                   calc_coord = .false., &
                   myrank = 0)

      call get_com(option%mode(2),       &
                   traj(2),              &
                   com(2),               &
                   setup = .true.,       &
                   calc_coord = .false., &
                   myrank = 0)

      nmol(1) = com(1)%nmol
      nmol(2) = com(2)%nmol

      ! Prepare LJ coefficients & point charges
      !
      call alloc_pot((/traj(1)%natm, traj(2)%natm/), pot)
      call alloc_pot((/traj(1)%natm, traj(2)%natm/), pot_dual)
      if (option%parmformat == ParmFormatANAPARM) then
        call setup_pot(traj, pot, anaparm = anaparm)
      else if (option%parmformat == ParmFormatPRMTOP) then
        call setup_pot(traj, pot, prmtop = prmtop)
      
        if (option%pme_dual) then
          call setup_pot(traj, pot_dual, prmtop = prmtop2)
        else
          call setup_pot(traj, pot_dual, prmtop = prmtop)
        end if
      
      end if
      
      if (.not. option%calc_vdw) then
        pot%acoef  = 0.0d0
        pot%bcoef  = 0.0d0
        pot%rwell2 = 0.0d0
        pot%uljmin = 0.0d0
        pot%ljsgm  = 0.0d0
        pot%ljeps  = 0.0d0
      
        if (option%pme_dual) then
          pot_dual%acoef  = 0.0d0
          pot_dual%bcoef  = 0.0d0
          pot_dual%rwell2 = 0.0d0
          pot_dual%uljmin = 0.0d0
          pot_dual%ljsgm  = 0.0d0
          pot_dual%ljeps  = 0.0d0
        end if
      end if
     
      ! Calculate Square of Cutoff Distance 
      !
      rljcut2 = option%rljcut * option%rljcut

      allocate(uij(nmol(2), nmol(1)))
      allocate(frc(3, nmol(2), nmol(1)))
      
      if (option%elec == ElecTypePME) then

        allocate(uLJ  (nmol(2), nmol(1)), &
                 uES  (nmol(2), nmol(1)), &
                 uEL  (nmol(2), nmol(1)), &
                 uCORR(nmol(2), nmol(1)))
        ! for dual
        allocate(uLJ_dual  (nmol(2), nmol(1)),    &
                 uES_dual  (nmol(2), nmol(1)),    &
                 uEL_dual  (nmol(2), nmol(1)),    &
                 uCORR_dual(nmol(2), nmol(1)),    &
                 frc_dual  (3, nmol(2), nmol(1)), &
                 uij_dual  (nmol(2), nmol(1)),    &
                 uij_diff  (nmol(2), nmol(1)),    &
                 frc_diff  (3, nmol(2), nmol(1)))

      end if

      ! Setup output files
      !
      write(fout_ij,      '(a,".uij")')     trim(output%fhead)
      write(fout_fsq,     '(a,".fsq")')     trim(output%fhead)
      write(fout_ij_dual, '(a,".uijdual")') trim(output%fhead)
      write(fout_fsq_dual,'(a,".fsqdual")') trim(output%fhead)
      write(fout_ij_diff, '(a,".uijdiff")') trim(output%fhead)
      write(fout_fsq_diff,'(a,".fsqdiff")') trim(output%fhead)

      call open_file(fout_ij,  iunit_ij)
      call open_file(fout_fsq, iunit_fsq)

      !open(UnitOut, file=trim(fout_ij))
      !open(UnitOut2,file=trim(fout_fsq))

      if (option%pme_dual) then
        !open(UnitOut3,file=trim(fout_ij_dual))
        !open(UnitOut4,file=trim(fout_fsq_dual))
        !open(UnitOut5,file=trim(fout_ij_diff))
        !open(UnitOut6,file=trim(fout_fsq_diff))

        call open_file(fout_ij_dual,  iunit_ijdual)
        call open_file(fout_fsq_dual, iunit_fsqdual)
        call open_file(fout_ij_diff,  iunit_ijdiff)
        call open_file(fout_fsq_diff, iunit_fsqdiff)
      end if

      ! Get trajectory types
      !
      call get_trajtype(input%ftraj(1), trajtype)

      !
      !
      istep_tot = 0
      do itraj = 1, input%ntraj
        call open_trajfile(input%ftraj(itraj), trajtype, iunit, dcd, xtc)
        call init_trajfile(trajtype, iunit, dcd, xtc, natm)

        is_end = .false.
        istep  = 0
        do while (.not. is_end)
          istep     = istep     + 1
          istep_tot = istep_tot + 1 

          call read_trajfile_oneframe(trajtype, iunit, istep, dcd, xtc, is_end)

          if (is_end) exit

          call send_coord_to_traj(1, trajtype, dcd, xtc, traj(1))
          call send_coord_to_traj(1, trajtype, dcd, xtc, traj(2))

          ! Calculate Pair Energy & Force
          !
          uij = 0.0d0
          frc = 0.0d0
         
          uij_dual = 0.0d0
          frc_dual = 0.0d0

          if (option%elec == ElecTypePME) then
        
            write(iw,'("Analyze> Error.")')
            write(iw,'("Sorry, PME schme is not implemented in Fugaku ver.")')
            stop

          else
            
            if (option%pbc) then
         
              if (option%vdw == VdwTypeSTANDARD) then
                call calc_uij_pbc(1, option, com, traj, pot, rljcut2,  &
                                  uij, frc)

              else if (option%vdw == VdwTypeATTRACTIVE) then
                call calc_uij_pbc_attractive(1, option, com, traj, pot, &
                                             rljcut2, uij, frc)
           
              else if (option%vdw == VdwTypeREPULSIVE) then
                call calc_uij_pbc_repulsive(1, option, com, traj, pot,  &
                                            rljcut2, uij, frc)
           
              end if
            else
              if (option%vdw == VdwTypeSTANDARD) then
                call calc_uij_nobc(1, option, com, traj, pot, rljcut2,   &
                                   uij, frc)
           
              else if (option%vdw == VdwTypeATTRACTIVE) then
                call calc_uij_nobc_attractive(1, option, com, traj, pot, &
                                              rljcut2, uij, frc)
           
              else if (option%vdw == VdwTypeREPULSIVE) then
                call calc_uij_nobc_repulsive(1, option, com, traj, pot,  &
                                             rljcut2, uij, frc)
           
              end if
           
            end if
         
          end if

          uij_diff = uij_dual - uij
          frc_diff = frc_dual - frc
         
          if (istep_tot == 1) then
            write(iw,'("=========================================")')
            write(iw,'("      Step             Euv_tot (kcal/mol)")')
            write(iw,'("=========================================")')
          end if
          write(iw,'(i10,11x,f20.10)') istep_tot, sum(uij(1:nmol(2), 1)) 
        
          call out_uij(iunit_ij,  traj(1)%dt, istep_tot, nmol, uij)
          call out_fsq(iunit_fsq, traj(1)%dt, istep_tot, nmol, frc) 
         
          if (option%pme_dual) then
        
            call out_uij(iunit_ijdual,  traj(1)%dt, istep_tot, nmol, uij_dual)
            call out_fsq(iunit_fsqdual, traj(1)%dt, istep_tot, nmol, frc_dual)

            call out_uij(iunit_ijdiff,  traj(1)%dt, istep_tot, nmol, uij_diff) 
            call out_fsq(iunit_fsqdiff, traj(1)%dt, istep_tot, nmol, frc_diff)  
         
          end if

        end do ! step

        istep_tot = istep_tot - 1

        call close_trajfile(trajtype, iunit, dcd, xtc)

      end do   ! traj

      close(iunit_ij)
      close(iunit_fsq)

      if (option%pme_dual) then
        close(iunit_ijdual)
        close(iunit_fsqdual)
        close(iunit_ijdiff)
        close(iunit_fsqdiff)
      end if


      call dealloc_pot(pot)
      deallocate(uij)

      if (allocated(uLJ)) then
        deallocate(uLJ,      uES,      uEL,      uCORR)
        deallocate(uLJ_dual, uES_dual, uEL_dual, uCORR_dual)
      end if

      deallocate(frc)

      !if (option%elec == ElecTypePME) then
      !  call dealloc_pmevar(pmevar)
      !end if

    end subroutine analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine calc_uij_pbc(istep, option, com, traj, pot, rljcut2, uij, &
                            frc) 
!-----------------------------------------------------------------------
      implicit none

      integer,        intent(in)    :: istep
      type(s_option), intent(in)    :: option
      type(s_com),    intent(in)    :: com(2)
      type(s_traj),   intent(in)    :: traj(2)
      type(s_pot),    intent(in)    :: pot
      real(8),        intent(in)    :: rljcut2
      real(8),        intent(inout) :: uij(:, :)
      real(8),        intent(inout) :: frc(:, :, :)

      integer :: iatm, jatm, imol, jmol
      real(8) :: d(3), r, r2, r2inv, r4inv, r6inv
      real(8) :: aij, bij, u, ulj, uele, fcoef, fj(3)


      uij   = 0.0d0
      frc   = 0.0d0

      !$omp parallel private(iatm, jatm, imol, jmol, d, r, r2, r2inv, r4inv, r6inv, &
      !$omp                & aij, bij, u, ulj, uele, fcoef, fj), &
      !$omp          shared(istep), reduction(+:uij), reduction(+:frc)

      !$omp do
      do iatm = 1, traj(1)%natm
        imol = com(1)%molid(iatm)
        do jatm = 1, traj(2)%natm 
          jmol            = com(2)%molid(jatm)
          d(1:3)          = traj(2)%coord(1:3, jatm, istep) &
                          - traj(1)%coord(1:3, iatm, istep) 
          d(1:3)          = d(1:3) - traj(1)%box(1:3, istep) &
                          * anint(d(1:3) / traj(1)%box(1:3, istep))
      
          r2              = dot_product(d, d)
          r2inv           = 1.0d0 / r2
          r4inv           = r2inv * r2inv
          r6inv           = r2inv * r4inv
      
          aij             = pot%acoef(jatm, iatm)
          bij             = pot%bcoef(jatm, iatm)

          if (r2 < rljcut2) then 
            u       = r6inv * (aij*r6inv - bij)
            fcoef   = r6inv * (12.0d0 * aij * r6inv - 6.0d0 * bij) 
            fj(1:3) = fcoef * r2inv * d(1:3)
          else
            u       = 0.0d0
            fj(1:3) = 0.0d0
          end if

          ulj                  = u
          uij(jmol, imol)      = uij(jmol, imol)      + u
          frc(1:3, jmol, imol) = frc(1:3, jmol, imol) + fj(1:3) 
     
          ! electrostatic
          !
          if (option%calc_elec) then
            r = sqrt(r2)
            if (r < option%relcut) then
              u       = pot%qq(jatm, iatm) / r
              fj(1:3) = u * r2inv * d(1:3)
            else
              u       = 0.0d0
              fj(1:3) = 0.0d0
            end if
            uele                 = u
            uij(jmol, imol)      = uij(jmol, imol)      + u
            frc(1:3, jmol, imol) = frc(1:3, jmol, imol) + fj(1:3)
          end if

          ! tot-attractive
          !
          if (option%calc_elec .and. &
              option%tottype == TotTypeATTRACTIVE) then
            if (r2 <= pot%rtotwell2(jatm, iatm)) then
              uij(jmol, imol)      = uij(jmol, imol) - ulj - uele &
                                   + pot%utotmin(jatm, iatm)
              frc(1:3, jmol, imol) = 0.0d0 ! not implemented
            end if
          end if

        end do
      end do
      !$omp end do
      !$omp end parallel

    end subroutine calc_uij_pbc
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine calc_uij_pbc_attractive(istep, option, com, traj, pot,   &
                                       rljcut2, uij, frc) 
!-----------------------------------------------------------------------
      implicit none

      integer,        intent(in)    :: istep
      type(s_option), intent(in)    :: option
      type(s_com),    intent(in)    :: com(2)
      type(s_traj),   intent(in)    :: traj(2)
      type(s_pot),    intent(in)    :: pot 
      real(8),        intent(in)    :: rljcut2
      real(8),        intent(inout) :: uij(:, :)
      real(8),        intent(inout) :: frc(:, :, :)

      integer :: iatm, jatm, imol, jmol
      real(8) :: d(3), r, r2, r2inv, r4inv, r6inv
      real(8) :: aij, bij, u, fcoef, fj(3)


      uij = 0.0d0

      !$omp parallel private(iatm, jatm, imol, jmol, d, r, r2, r2inv, r4inv, r6inv, &
      !$omp                & aij, bij, u, fcoef, fj), &
      !$omp          shared(istep), reduction(+:uij), reduction(+:frc)

      !$omp do
      do iatm = 1, traj(1)%natm
        imol = com(1)%molid(iatm)
        do jatm = 1, traj(2)%natm 
          jmol            = com(2)%molid(jatm)
          d(1:3)          = traj(2)%coord(1:3, jatm, istep) &
                          - traj(1)%coord(1:3, iatm, istep) 
          d(1:3)          = d(1:3) - traj(1)%box(1:3, istep) &
                          * anint(d(1:3) / traj(1)%box(1:3, istep))
      
          r2              = dot_product(d, d)
          r2inv           = 1.0d0 / r2
          r4inv           = r2inv * r2inv
          r6inv           = r2inv * r4inv
      
          aij             = pot%acoef(jatm, iatm)
          bij             = pot%bcoef(jatm, iatm)

          if (r2 < rljcut2) then 
            u       = r6inv * (aij*r6inv - bij)
            fcoef   = r6inv * (12.0d0 * aij * r6inv - 6.0d0 * bij) 
            fj(1:3) = fcoef * r2inv * d(1:3)
            if (r2 < pot%rwell2(jatm, iatm)) then
              u       = pot%uljmin(jatm, iatm)
              fj(1:3) = 0.0d0 
            end if
          else
            u       = 0.0d0
            fj(1:3) = 0.0d0
          end if
      
          uij(jmol, imol)      = uij(jmol, imol)      + u
          frc(1:3, jmol, imol) = frc(1:3, jmol, imol) + fj(1:3) 

          ! electrostatic
          !
          if (option%calc_elec) then
            r = sqrt(r2)
            if (r < option%relcut) then
              u       = pot%qq(jatm, iatm) / r
              fj(1:3) = u * r2inv * d(1:3)
            else
              u       = 0.0d0
              fj(1:3) = 0.0d0
            end if
            uij(jmol, imol)      = uij(jmol, imol)      + u
            frc(1:3, jmol, imol) = frc(1:3, jmol, imol) + fj(1:3)
          end if

        end do
      end do
      !$omp end do
      !$omp end parallel

    end subroutine calc_uij_pbc_attractive
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine calc_uij_pbc_repulsive(istep, option, com, traj, pot,    &
                                      rljcut2, uij, frc) 
!-----------------------------------------------------------------------
      implicit none

      integer,        intent(in)    :: istep
      type(s_option), intent(in)    :: option
      type(s_com),    intent(in)    :: com(2)
      type(s_traj),   intent(in)    :: traj(2)
      type(s_pot),    intent(in)    :: pot 
      real(8),        intent(in)    :: rljcut2
      real(8),        intent(inout) :: uij(:, :)
      real(8),        intent(inout) :: frc(:, :, :)

      integer :: iatm, jatm, imol, jmol
      real(8) :: d(3), r, r2, r2inv, r4inv, r6inv
      real(8) :: aij, bij, u, fcoef, fj(3)


      uij = 0.0d0

      !$omp parallel private(iatm, jatm, imol, jmol, d, r, r2, r2inv, r4inv, r6inv, &
      !$omp                & aij, bij, u, fcoef, fj), &
      !$omp          shared(istep), reduction(+:uij), reduction(+:frc)

      !$omp do
      do iatm = 1, traj(1)%natm
        imol = com(1)%molid(iatm)
        do jatm = 1, traj(2)%natm 
          jmol            = com(2)%molid(jatm)
          d(1:3)          = traj(2)%coord(1:3, jatm, istep) &
                          - traj(1)%coord(1:3, iatm, istep) 
          d(1:3)          = d(1:3) - traj(1)%box(1:3, istep) &
                          * anint(d(1:3) / traj(1)%box(1:3, istep))
      
          r2              = dot_product(d, d)
          r2inv           = 1.0d0 / r2
          r4inv           = r2inv * r2inv
          r6inv           = r2inv * r4inv
      
          aij             = pot%acoef(jatm, iatm)
          bij             = pot%bcoef(jatm, iatm)

          if (r2 < rljcut2) then 
            u       = r6inv * (aij*r6inv - bij)
            fcoef   = r6inv * (12.0d0 * aij * r6inv - 6.0d0 * bij) 
            fj(1:3) = fcoef * r2inv * d(1:3)
            if (r2 < pot%rwell2(jatm, iatm)) then
              u = u - pot%uljmin(jatm, iatm)
            else
              u       = 0.0d0
              fj(1:3) = 0.0d0 
            end if
          else
            u       = 0.0d0
            fj(1:3) = 0.0d0
          end if
      
          uij(jmol, imol)      = uij(jmol, imol) + u
          frc(1:3, jmol, imol) = frc(1:3, jmol, imol) &
                               + fj(1:3)
      
          ! electrostatic
          !
          if (option%calc_elec) then
            r = sqrt(r2)
            if (r < option%relcut) then
              u       = pot%qq(jatm, iatm) / r
              fj(1:3) = u * r2inv * d(1:3) 
            else
              u       = 0.0d0
              fj(1:3) = 0.0d0
            end if
            uij(jmol, imol)      = uij(jmol, imol)      + u
            frc(1:3, jmol, imol) = frc(1:3, jmol, imol) + fj(1:3)
          end if

        end do
      end do
      !$omp end do
      !$omp end parallel

    end subroutine calc_uij_pbc_repulsive
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine calc_uij_nobc(istep, option, com, traj, pot, rljcut2,    &
                             uij, frc) 
!-----------------------------------------------------------------------
      implicit none

      integer,        intent(in)    :: istep
      type(s_option), intent(in)    :: option
      type(s_com),    intent(in)    :: com(2)
      type(s_traj),   intent(in)    :: traj(2)
      type(s_pot),    intent(in)    :: pot 
      real(8),        intent(in)    :: rljcut2
      real(8),        intent(inout) :: uij(:, :)
      real(8),        intent(inout) :: frc(:, :, :)

      integer :: iatm, jatm, imol, jmol
      real(8) :: d(3), r, r2, r2inv, r4inv, r6inv
      real(8) :: aij, bij, u, ulj, uele, fcoef, fj(1:3)


      uij = 0.0d0

      !$omp parallel private(iatm, jatm, imol, jmol, d, r, r2, r2inv, r4inv, r6inv, &
      !$omp                & aij, bij, u, ulj, uele, fcoef, fj), &
      !$omp          shared(istep), reduction(+:uij), reduction(+:frc)

      !$omp do
      do iatm = 1, traj(1)%natm
        imol = com(1)%molid(iatm)
        do jatm = 1, traj(2)%natm 
          jmol            = com(2)%molid(jatm)
          d(1:3)          = traj(2)%coord(1:3, jatm, istep) &
                          - traj(1)%coord(1:3, iatm, istep) 
      
          r2              = dot_product(d, d)
          r2inv           = 1.0d0 / r2
          r4inv           = r2inv * r2inv
          r6inv           = r2inv * r4inv
      
          aij             = pot%acoef(jatm, iatm)
          bij             = pot%bcoef(jatm, iatm)

          if (r2 < rljcut2) then 
            u       = r6inv * (aij*r6inv - bij)
            fcoef   = r6inv * (12.0d0 * aij * r6inv - 6.0d0 * bij) 
            fj(1:3) = fcoef * r2inv * d(1:3)
          else
            u       = 0.0d0
            fj(1:3) = 0.0d0 
          end if
          
          ulj                  = u
          uij(jmol, imol)      = uij(jmol, imol)      + u
          frc(1:3, jmol, imol) = frc(1:3, jmol, imol) + fj(1:3)
      
          ! electrostatic
          !
          if (option%calc_elec) then
            r = sqrt(r2)
            if (r < option%relcut) then
              u       = pot%qq(jatm, iatm) / r
              fj(1:3) = u * r2inv * d(1:3) 
            else
              u       = 0.0d0
              fj(1:3) = u * r2inv * d(1:3) 
            end if
            uele                 = u
            uij(jmol, imol)      = uij(jmol, imol)      + u
            frc(1:3, jmol, imol) = frc(1:3, jmol, imol) + fj(1:3) 
          end if

          ! tot-attractive
          !
          if (option%tottype == TotTypeATTRACTIVE) then
            if (r2 <= pot%rtotwell2(jatm, iatm)) then
              uij(jmol, imol)      = uij(jmol, imol) - ulj - uele & 
                                   + pot%utotmin(jatm, iatm)
              frc(1:3, jmol, imol) = 0.0d0 ! not implemented
            end if
          end if

        end do
      end do
      !$omp end do
      !$omp end parallel

    end subroutine calc_uij_nobc
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine calc_uij_nobc_attractive(istep, option, com, traj, pot,  &
                                        rljcut2, uij, frc) 
!-----------------------------------------------------------------------
      implicit none

      integer,        intent(in)    :: istep
      type(s_option), intent(in)    :: option
      type(s_com),    intent(in)    :: com(2)
      type(s_traj),   intent(in)    :: traj(2)
      type(s_pot),    intent(in)    :: pot 
      real(8),        intent(in)    :: rljcut2
      real(8),        intent(inout) :: uij(:, :)
      real(8),        intent(inout) :: frc(:, :, :)

      integer :: iatm, jatm, imol, jmol
      real(8) :: d(3), r, r2, r2inv, r4inv, r6inv
      real(8) :: aij, bij, u, fcoef, fj(3)


      uij = 0.0d0

      !$omp parallel private(iatm, jatm, imol, jmol, d, r, r2, r2inv, r4inv, r6inv, &
      !$omp                & aij, bij, u, fcoef, fj), &
      !$omp          shared(istep), reduction(+:uij), reduction(+:frc)

      !$omp do
      do iatm = 1, traj(1)%natm
        imol = com(1)%molid(iatm)
        do jatm = 1, traj(2)%natm 
          jmol            = com(2)%molid(jatm)
          d(1:3)          = traj(2)%coord(1:3, jatm, istep) &
                          - traj(1)%coord(1:3, iatm, istep) 
      
          r2              = dot_product(d, d)
          r2inv           = 1.0d0 / r2
          r4inv           = r2inv * r2inv
          r6inv           = r2inv * r4inv
      
          aij             = pot%acoef(jatm, iatm)
          bij             = pot%bcoef(jatm, iatm)

          if (r2 < rljcut2) then 
            u       = r6inv * (aij*r6inv - bij)
            fcoef   = r6inv * (12.0d0 * aij * r6inv - 6.0d0 * bij) 
            fj(1:3) = fcoef * r2inv * d(1:3)
            if (r2 < pot%rwell2(jatm, iatm)) then
              u       = pot%uljmin(jatm, iatm)
              fj(1:3) = 0.0d0
            end if
          else
            u       = 0.0d0
            fj(1:3) = 0.0d0
          end if
      
          uij(jmol, imol)      = uij(jmol, imol) + u
          frc(1:3, jmol, imol) = frc(1:3, jmol, imol) + fj(1:3)
      
          ! electrostatic
          !
          if (option%calc_elec) then
            r = sqrt(r2)
            if (r < option%relcut) then
              u       = pot%qq(jatm, iatm) / r
              fj(1:3) = u * r2inv * d(1:3) 
            else
              u       = 0.0d0
              fj(1:3) = 0.0d0
            end if
            uij(jmol, imol)      = uij(jmol, imol) + u
            frc(1:3, jmol, imol) = frc(1:3, jmol, imol) &
                                 + fj(1:3)
          end if

        end do
      end do
      !$omp end do
      !$omp end parallel

    end subroutine calc_uij_nobc_attractive
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine calc_uij_nobc_repulsive(istep, option, com, traj, pot,   &
                                       rljcut2, uij, frc) 
!-----------------------------------------------------------------------
      implicit none

      integer,        intent(in)    :: istep
      type(s_option), intent(in)    :: option
      type(s_com),    intent(in)    :: com(2)
      type(s_traj),   intent(in)    :: traj(2)
      type(s_pot),    intent(in)    :: pot 
      real(8),        intent(in)    :: rljcut2
      real(8),        intent(inout) :: uij(:, :)
      real(8),        intent(inout) :: frc(:, :, :)

      integer :: iatm, jatm, imol, jmol
      real(8) :: d(3), r, r2, r2inv, r4inv, r6inv
      real(8) :: aij, bij, u, fcoef, fj(3)


      uij = 0.0d0

      !$omp parallel private(iatm, jatm, imol, jmol, d, r, r2, r2inv, r4inv, r6inv, &
      !$omp                & aij, bij, u, fcoef, fj), &
      !$omp          shared(istep), reduction(+:uij), reduction(+:frc)

      !$omp do
      do iatm = 1, traj(1)%natm
        imol = com(1)%molid(iatm)
        do jatm = 1, traj(2)%natm 
          jmol            = com(2)%molid(jatm)
          d(1:3)          = traj(2)%coord(1:3, jatm, istep) &
                          - traj(1)%coord(1:3, iatm, istep) 
      
          r2              = dot_product(d, d)
          r2inv           = 1.0d0 / r2
          r4inv           = r2inv * r2inv
          r6inv           = r2inv * r4inv
      
          aij             = pot%acoef(jatm, iatm)
          bij             = pot%bcoef(jatm, iatm)

          if (r2 < rljcut2) then 
            u       = r6inv * (aij*r6inv - bij)
            fcoef   = r6inv * (12.0d0 * aij * r6inv - 6.0d0 * bij) 
            fj(1:3) = fcoef * r2inv * d(1:3)
            if (r2 < pot%rwell2(jatm, iatm)) then
              u = u - pot%uljmin(jatm, iatm)
            else
              u       = 0.0d0
              fj(1:3) = 0.0d0
            end if
          else
            u       = 0.0d0
            fj(1:3) = 0.0d0
          end if
      
          uij(jmol, imol)      = uij(jmol, imol)      + u
          frc(1:3, jmol, imol) = frc(1:3, jmol, imol) + fj(1:3)
      
          ! electrostatic
          !
          if (option%calc_elec) then
            r = sqrt(r2)
            if (r < option%relcut) then
              u       = pot%qq(jatm, iatm) / r
              fj(1:3) = u * r2inv * d(1:3) 
            else
              u       = 0.0d0
              fj(1:3) = 0.0d0
            end if
            uij(jmol, imol)      = uij(jmol, imol)      + u
            frc(1:3, jmol, imol) = frc(1:3, jmol, imol) + fj(1:3)
          end if

        end do
      end do
      !$omp end do
      !$omp end parallel

    end subroutine calc_uij_nobc_repulsive
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine out_uij(iunit, dt, istep, nmol, uij)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) :: iunit
      real(8), intent(in) :: dt
      integer, intent(in) :: istep
      integer, intent(in) :: nmol(2)
      real(8), intent(in) :: uij(:, :)

      integer :: imol, jmol


      ! Write Pair Energy at istep  
      !
      write(iunit,'(f20.10)',advance="no") dt * istep 
      do imol = 1, nmol(1)
        do jmol = 1, nmol(2)
          write(iunit,'(f20.10)',advance="no") uij(jmol, imol)
        end do
      end do
      write(iunit,*)

    end subroutine out_uij
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine out_fsq(iunit, dt, istep, nmol, frc)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) :: iunit
      real(8), intent(in) :: dt
      integer, intent(in) :: istep
      integer, intent(in) :: nmol(2)
      real(8), intent(in) :: frc(:, :, :)

      integer :: imol, jmol


      write(iunit,'(f20.10)',advance="no") dt * istep 
      do imol = 1, nmol(1)
        do jmol = 1, nmol(2)
          write(iunit,'(f20.10)',advance="no") &
            dot_product(frc(1:3, jmol, imol), frc(1:3, jmol, imol))
        end do
      end do
      write(iunit,*)

    end subroutine out_fsq
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
