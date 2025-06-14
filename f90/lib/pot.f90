!-------------------------------------------------------------------------------
!
!  Module   mod_potential 

!    Module for define potential (force field) parameters  
!
!    (c) Copyright 2024 Osaka Univ. All rights reserved.
!
!-------------------------------------------------------------------------------

module mod_potential

  use mod_util
  use mod_const
  use mod_prmtop
  use mod_anaparm
  use mod_traj

  implicit none

  ! parameters
  !
  integer,      parameter :: PotentialTypeUNDEFINED = 0
  integer,      parameter :: PotentialTypeAMBER     = 1
  integer,      parameter :: PotentialTypeCHARMM    = 2
  character(*), parameter :: PotentialTypes(2)      = (/'AMBER ', 'CHARMM'/)

  integer,      parameter :: nsearch     = 100
  real(8),      parameter :: rmin_search = 0.1d0
  real(8),      parameter :: rmax_search = 20.0d0
  real(8),      parameter :: dr_search   = 0.1
  real(8),      parameter :: rhuge       = 1.0d50 

  ! structures
  !
  type :: s_pot
    integer :: natm(2)
    integer :: potential_type = PotentialTypeUNDEFINED

    ! AMBER-FF
    !
    !   U_LJ(r) = A / r^12 - B / r^6 
    !
    real(8), allocatable :: acoef(:, :)    ! LJ A parameter
    real(8), allocatable :: bcoef(:, :)    ! LJ B parameter
    !
    !   U_LJ(r) = 4 * epsilon * {(sigma/r)^12 - (sigma/r)^6}
    !
    real(8), allocatable :: ljsgm(:, :)    ! LJ sigma
    real(8), allocatable :: ljeps(:, :)    ! LJ epsilon
    !
    !
    real(8), allocatable :: rwell2(:, :)   ! r^2 at ULJ minimum
    real(8), allocatable :: uljmin(:, :)   ! ULJ minimum
    !
    !   U_coulomb(r) = qq/r
    !
    real(8), allocatable :: q1(:), q2(:)
    real(8), allocatable :: qq(:, :)

  end type s_pot

  ! subroutines
  !
  public   :: alloc_pot
  public   :: dealloc_pot
  public   :: setup_pot
  private  :: setup_pot_amber
  private  :: setup_pot_anaparm
  private  :: ftot 
  !private  :: seek_prmtop_flag
  !private  :: read_prmtop_integer 
  !private  :: read_prmtop_real

  contains

    !---------------------------------------------------------------------------
    !
    !  Subroutine  alloc_pot 
    !  @brief      allocate memories for POT structure 
    !  @authors    KK
    !  @param[in]  natm : # of atoms 
    !  @param[out] pot  : structure of potential information 
    !
    !---------------------------------------------------------------------------

    subroutine alloc_pot(natm, pot)

      implicit none

      ! formal arguments
      integer,     intent(in)    :: natm(2)
      type(s_pot), intent(inout) :: pot 
      
      
      allocate(pot%acoef(natm(2), natm(1)))
      allocate(pot%bcoef(natm(2), natm(1)))
      allocate(pot%ljsgm(natm(2), natm(1)))
      allocate(pot%ljeps(natm(2), natm(1)))
      allocate(pot%rwell2(natm(2), natm(1)))
      allocate(pot%uljmin(natm(2), natm(1)))
      allocate(pot%q1(natm(1)))
      allocate(pot%q2(natm(2)))
      allocate(pot%qq(natm(2), natm(1)))

    end subroutine alloc_pot 

    !---------------------------------------------------------------------------
    !
    !  Subroutine  dealloc_pot 
    !  @brief      deallocate memories for POT structure 
    !  @authors    KK
    !  @param[inout] pot : sturcture of potential information 
    !
    !---------------------------------------------------------------------------

    subroutine dealloc_pot(pot)

      implicit none

      type(s_pot), intent(inout) :: pot 
      
     
      if (allocated(pot%acoef)) then
        deallocate(pot%acoef)
        deallocate(pot%bcoef)
      end if

      if (allocated(pot%ljsgm)) then
        deallocate(pot%ljsgm)
        deallocate(pot%ljeps)
      end if
      
      if (allocated(pot%uljmin)) then
        deallocate(pot%rwell2)
        deallocate(pot%uljmin)
      end if

      if (allocated(pot%q1)) then
        deallocate(pot%q1, pot%q2)
      end if

      if (allocated(pot%qq)) then
        deallocate(pot%qq)
      end if

    end subroutine dealloc_pot 


    !---------------------------------------------------------------------------
    !
    !  Subroutine  setup_pot 
    !  @brief      setup potential paramters 
    !  @authors    KK
    !  @param[in]    traj(2) : structure of ANATRA trajectories 
    !  @param[inout] pot     : structure of potential information 
    !  @param[in, optional] prmtop  : structure of AMBER parm7 file info.
    !  @param[in, optional] anaparm : structure of ANATRA parameter (ANAPARM) 
    !                                 file info.
    !
    !---------------------------------------------------------------------------

    subroutine setup_pot(traj, pot, prmtop, anaparm)

      implicit none

      ! formal arguments
      type(s_traj),    intent(in)           :: traj(2)
      type(s_pot),     intent(inout)        :: pot
      type(s_prmtop),  intent(in), optional :: prmtop
      type(s_anaparm), intent(in), optional :: anaparm 
     

      if (present(prmtop)) then
        call setup_pot_amber(traj, prmtop, pot)
      else if (present(anaparm)) then
        call setup_pot_anaparm(traj, anaparm, pot)
      end if

    end subroutine setup_pot


    !---------------------------------------------------------------------------
    !
    !  Subroutine    setup_pot_amber 
    !  @brief        setup potential paramters from AMBER parm7 file info. 
    !  @authors      KK
    !  @param[in]    traj(2) : structure of ANATRA trajectories 
    !  @param[in]    prmtop  : structure of AMBER parm7 file info.
    !  @param[inout] pot     : structure of potential information 
    !
    !---------------------------------------------------------------------------

    subroutine setup_pot_amber(traj, prmtop, pot)

      implicit none

      ! formal arguments
      type(s_traj),   intent(in)           :: traj(2)
      type(s_prmtop), intent(in)           :: prmtop
      type(s_pot),    intent(inout)        :: pot

      integer :: iatm, jatm, iac1, iac2, ico, ipair
      integer :: ntypes
      real(8) :: a, b, s, e


      pot%natm(1) = traj(1)%natm
      pot%natm(2) = traj(2)%natm

      ntypes = prmtop%ntypes
      do iatm = 1, traj(1)%natm
        iac1         = prmtop%iac(traj(1)%ind(iatm))
        pot%q1(iatm) = prmtop%charge(traj(1)%ind(iatm))
        do jatm = 1, traj(2)%natm

          iac2         = prmtop%iac(traj(2)%ind(jatm))
          pot%q2(jatm) = prmtop%charge(traj(2)%ind(jatm)) 

          ipair = ntypes * (max(iac1, iac2) - 1) + min(iac1, iac2)
          ico   = prmtop%ico(ipair)

          ! LJ parameters
          !
          a     = prmtop%cn1(ico)
          b     = prmtop%cn2(ico)

          pot%acoef(jatm, iatm)  = a
          pot%bcoef(jatm, iatm)  = b
          if (a < 1.0d-8 .or. b < 1.0d-8) then
            pot%rwell2(jatm, iatm) = 0.0d0 
            pot%uljmin(jatm, iatm) = 0.0d0 
          else
            pot%rwell2(jatm, iatm) = (2.0d0 * a / b)**(1.0d0 / 3.0d0)
            pot%uljmin(jatm, iatm) = - 0.25d0 * b * b / a
          end if


          if (b > 1.0d-8) then
            s                 = (a / b) ** (1.0d0 / 6.0d0)
            e                 = b / (4.0d0 * (a / b))
            pot%ljsgm(jatm, iatm) = s
            pot%ljeps(jatm, iatm) = e
          else
            pot%ljsgm(jatm, iatm) = 0.0d0 
            pot%ljeps(jatm, iatm) = 0.0d0 
          end if

          ! Charge parameters
          !
          pot%qq(jatm, iatm) = prmtop%charge(traj(1)%ind(iatm)) &
                             * prmtop%charge(traj(2)%ind(jatm))


        end do
      end do


    end subroutine setup_pot_amber

    !---------------------------------------------------------------------------
    !
    !  Subroutine    setup_pot_anaparm
    !  @brief        setup potential paramters from ANAPARM file info. 
    !  @authors      KK
    !  @param[in]    traj(2) : structure of ANATRA trajectories 
    !  @param[in]    anaparm : structure of ANAPARM file info.
    !  @param[inout] pot     : structure of potential information 
    !
    !---------------------------------------------------------------------------

    subroutine setup_pot_anaparm(traj, anaparm, pot)

      implicit none

      ! parameters
      real(8), parameter :: charge_const = 332.05221729d0 

      ! formal arguments
      type(s_traj),    intent(in)           :: traj(2)
      type(s_anaparm), intent(in)           :: anaparm 
      type(s_pot),     intent(inout)        :: pot

      integer :: iatm, jatm, ind1, ind2
      real(8) :: si, sj, ei, ej, qi, qj
      real(8) :: sij, eij, a, b


      pot%natm(1) = traj(1)%natm
      pot%natm(2) = traj(2)%natm

      do iatm = 1, traj(1)%natm
        ind1 = traj(1)%ind(iatm)
        do jatm = 1, traj(2)%natm
          ind2 = traj(2)%ind(jatm)

          ! LJ parameters
          !
          si    = anaparm%sgm(ind1)
          sj    = anaparm%sgm(ind2)
          ei    = anaparm%eps(ind1)
          ej    = anaparm%eps(ind2)

          sij   = 0.5d0 * (si + sj)
          eij   = sqrt(ei * ej)

          a     = 4.0d0 * eij * sij**12.0d0
          b     = 4.0d0 * eij * sij**6.0d0

          pot%acoef(jatm, iatm)  = a  
          pot%bcoef(jatm, iatm)  = b
          if (a < 1.0d-8 .or. b < 1.0d-8) then
            pot%rwell2(jatm, iatm) = 0.0d0 
            pot%uljmin(jatm, iatm) = 0.0d0 
          else
            pot%rwell2(jatm, iatm) = (2.0d0 * a / b)**(1.0d0 / 3.0d0)
            pot%uljmin(jatm, iatm) = - 0.25d0 * b * b / a
          end if

          pot%ljsgm(jatm, iatm) = sij 
          pot%ljeps(jatm, iatm) = eij 

          ! Charge parameters
          !
          qi                 = anaparm%charge(ind1)
          qj                 = anaparm%charge(ind2)
          pot%q1(iatm)       = sqrt(charge_const) * qi
          pot%q2(jatm)       = sqrt(charge_const) * qj
          pot%qq(jatm, iatm) = charge_const * qi * qj
        end do
      end do

    end subroutine setup_pot_anaparm

    !---------------------------------------------------------------------------
    !
    !  Subroutine  modify_pot_for_excl 
    !  @brief      modify potential paramters 
    !  @authors    KK
    !  @param[in]    traj(2) : structure of ANATRA trajectories 
    !  @param[inout] pot     : structure of potential information 
    !  @param[in, optional] prmtop  : structure of AMBER parm7 file info.
    !  @param[in, optional] anaparm : structure of ANATRA parameter (ANAPARM) 
    !                                 file info.
    !
    !---------------------------------------------------------------------------

    subroutine modify_pot_for_excl(traj, pot, prmtop, anaparm)

      implicit none

      ! formal arguments
      type(s_traj),    intent(in)           :: traj(2)
      type(s_pot),     intent(inout)        :: pot
      type(s_prmtop),  intent(in), optional :: prmtop
      type(s_anaparm), intent(in), optional :: anaparm 
     

      if (present(prmtop)) then
        call modify_pot_for_excl_amber(traj, prmtop, pot)
      else if (present(anaparm)) then
        write(iw,'("Setup_Pot_For_Excl> Error: sorry, anaparm is currently not supported.")')
        stop
      end if

    end subroutine modify_pot_for_excl 


    !---------------------------------------------------------------------------
    !
    !  Subroutine  modify_pot_for_excl_amber 
    !  @brief      modify potential paramters for AMBER 
    !  @authors    KK
    !  @param[in]    traj(2) : structure of ANATRA trajectories 
    !  @param[in     prmtop  : structure of AMBER parm7 file info.
    !  @param[inout] pot     : structure of potential information 
    !
    !---------------------------------------------------------------------------

    subroutine modify_pot_for_excl_amber(traj, prmtop, pot)

      implicit none

      ! formal arguments
      type(s_traj),    intent(in)           :: traj(2)
      type(s_pot),     intent(inout)        :: pot
      type(s_prmtop),  intent(in)           :: prmtop

      integer :: iatm, jatm, iac1, iac2, ico, ipair, ie
      integer :: iid, jid, kid, ista, iend, ihed
      integer :: ih(5)
      integer :: natm, npair
      real(8) :: a, b, s, e

      integer, allocatable :: excl(:, :)
      real(8), allocatable :: scnb(:, :), scee(:, :)


      natm = traj(1)%natm

      allocate(excl(natm, natm), scnb(natm, natm), scee(natm, natm))
      excl = 1

      ! excluded pairs
      !
      do iatm = 1, natm - 1
        iid = traj(1)%ind(iatm)
        ista = sum(prmtop%num_excl_atm(1:iid-1))
        iend = ista + prmtop%num_excl_atm(iid)

        do jatm = iatm + 1, natm
          jid = traj(1)%ind(jatm)

          if (ista /= iend) then
            do ie = ista + 1, iend
              kid = prmtop%excl_atm_list(ie)

              if (jid == kid) then
                excl(jatm, iatm) = 0
              end if
            end do
          end if

        end do

      end do

      ! 1-4 interaction pairs
      !
      scnb = 1.0d0
      scee = 1.0d0

      !   for dihedrals involving hydrogen
      !
      do ihed = 1, prmtop%nphih
        ih(1) = prmtop%dihed_inc_hyd(5*(ihed-1)+1)
        ih(2) = prmtop%dihed_inc_hyd(5*(ihed-1)+2)
        ih(3) = prmtop%dihed_inc_hyd(5*(ihed-1)+3)
        ih(4) = prmtop%dihed_inc_hyd(5*(ihed-1)+4)
        ih(5) = prmtop%dihed_inc_hyd(5*(ihed-1)+5)

        ih(1) = ih(1) / 3 + 1
        ih(4) = abs(ih(4)) / 3 + 1

        if (ih(3) > 0) then
          do iatm = 1, natm - 1
            iid = traj(1)%ind(iatm)

            if (ih(1) == iid) then
              do jatm = iatm + 1, natm
                jid = traj(1)%ind(jatm)

                if (ih(4) == jid) then
                  excl(jatm, iatm) = 1
                  scnb(jatm, iatm) = prmtop%scnb_scale_fact(ih(5))
                  scee(jatm, iatm) = prmtop%scee_scale_fact(ih(5))
                end if

              end do
            end if

          end do
        end if

      end do

      !   for dihedrals involving only non-hydrogen
      !
      do ihed = 1, prmtop%nphia
        ih(1) = prmtop%dihed_wo_hyd(5*(ihed-1)+1)
        ih(2) = prmtop%dihed_wo_hyd(5*(ihed-1)+2)
        ih(3) = prmtop%dihed_wo_hyd(5*(ihed-1)+3)
        ih(4) = prmtop%dihed_wo_hyd(5*(ihed-1)+4)
        ih(5) = prmtop%dihed_wo_hyd(5*(ihed-1)+5)

        ih(1) = ih(1) / 3 + 1
        ih(4) = abs(ih(4)) / 3 + 1

        if (ih(3) > 0) then
          do iatm = 1, natm - 1
            iid = traj(1)%ind(iatm)

            if (ih(1) == iid) then
              do jatm = iatm + 1, natm
                jid = traj(1)%ind(jatm)

                if (ih(4) == jid) then
                  excl(jatm, iatm) = 1
                  scnb(jatm, iatm) = prmtop%scnb_scale_fact(ih(5))
                  scee(jatm, iatm) = prmtop%scee_scale_fact(ih(5))
                end if

              end do
            end if

          end do
        end if

      end do

      npair = 0
      do iatm = 1, natm - 1
        do jatm = iatm + 1, natm
          npair = npair + excl(jatm, iatm)
        end do
      end do

      ! modify parameters
      !
      do iatm = 1, natm - 1
        do jatm = iatm + 1, natm
          if (excl(jatm, iatm) == 0) then
            pot%acoef(jatm, iatm)  = 0.0d0
            pot%bcoef(jatm, iatm)  = 0.0d0
            pot%rwell2(jatm, iatm) = huge(1.0d0) 
            pot%uljmin(jatm, iatm) = 0.0d0
            pot%ljsgm(jatm, iatm)  = 0.0d0
            pot%ljeps(jatm, iatm)  = 0.0d0
            pot%qq(jatm, iatm)     = 0.0d0
          else
            a = pot%acoef(jatm, iatm) / scnb(jatm, iatm)
            b = pot%bcoef(jatm, iatm) / scnb(jatm, iatm)

            pot%acoef(jatm, iatm) = a
            pot%bcoef(jatm, iatm) = b 

            if (a < 1.0d-8 .or. b < 1.0d-8) then
              pot%rwell2(jatm, iatm) = 0.0d0 
              pot%uljmin(jatm, iatm) = 0.0d0 
            else
              pot%rwell2(jatm, iatm) = (2.0d0 * a / b)**(1.0d0 / 3.0d0)
              pot%uljmin(jatm, iatm) = - 0.25d0 * b * b / a

            end if

            if (b > 1.0d-8) then
              s                 = (a / b) ** (1.0d0 / 6.0d0)
              e                 = b / (4.0d0 * (a / b))
              pot%ljsgm(jatm, iatm) = s
              pot%ljeps(jatm, iatm) = e
            else
              pot%ljsgm(jatm, iatm) = 0.0d0 
              pot%ljeps(jatm, iatm) = 0.0d0 
            end if

            pot%qq(jatm, iatm) = pot%qq(jatm, iatm) / scee(jatm, iatm)
            
          end if
        end do
      end do
      
!      do iatm = 1, natm1
!        iid  = traj(1)%ind(iatm)
!
!        ista = sum(prmtop%num_excl_atm(1:iid-1))
!        iend = ista + prmtop%num_excl_atm(iid)
!
!        do jatm = 1, natm2 
!          jid = traj(2)%ind(jatm)
!
!          if (iid == jid) then
!            pot%acoef(jatm, iatm)  = 0.0d0
!            pot%bcoef(jatm, iatm)  = 0.0d0
!            pot%rwell2(jatm, iatm) = huge(1.0d0) 
!            pot%uljmin(jatm, iatm) = 0.0d0
!            pot%ljsgm(jatm, iatm)  = 0.0d0
!            pot%ljeps(jatm, iatm)  = 0.0d0
!            pot%qq(jatm, iatm)     = 0.0d0
!          end if
!
!          if (ista /= iend) then
!            do ie = ista + 1, iend
!              kid = prmtop%excl_atm_list(ie)
!
!              if (jid == kid) then
!                pot%acoef(jatm, iatm)  = 0.0d0
!                pot%bcoef(jatm, iatm)  = 0.0d0
!                pot%rwell2(jatm, iatm) = huge(1.0d0) 
!                pot%uljmin(jatm, iatm) = 0.0d0
!                pot%ljsgm(jatm, iatm)  = 0.0d0
!                pot%ljeps(jatm, iatm)  = 0.0d0
!                pot%qq(jatm, iatm)     = 0.0d0
!              end if
!
!            end do
!          end if  
!        end do 
!      end do
       
    end subroutine modify_pot_for_excl_amber

    !---------------------------------------------------------------------------
    !
    !  Function      ftot 
    !  @brief        calculate force from total pair potential 
    !  @authors      KK
    !  @param[in]    aij     : LJ-Aij parameter 
    !  @param[in]    bij     : LJ-Bij parameter 
    !  @param[in]    qq      : electrostatic potential parameter 
    !  @param[in]    r       : interatomic distance 
    !
    !---------------------------------------------------------------------------

    function ftot(aij, bij, qq, r)

      implicit none

      ! formal arguments
      real(8), intent(in) :: aij
      real(8), intent(in) :: bij
      real(8), intent(in) :: qq
      real(8), intent(in) :: r

      real(8) :: ftot

      real(8) :: rinv


      ftot = 0.0d0
      rinv = 1.0d0 / r

      ftot = - qq * rinv**2 & 
             - 12.0d0 * aij * rinv**13 + 6.0d0 * bij * rinv**7 
      
!
    end function ftot

    !---------------------------------------------------------------------------
    !
    !  Function      utot 
    !  @brief        calculate total pair potential 
    !  @authors      KK
    !  @param[in]    aij     : LJ-Aij parameter 
    !  @param[in]    bij     : LJ-Bij parameter 
    !  @param[in]    qq      : electrostatic potential parameter 
    !  @param[in]    r       : interatomic distance 
    !
    !---------------------------------------------------------------------------

    function utot(aij, bij, qq, r)

      implicit none

      ! formal arguments
      real(8), intent(in) :: aij
      real(8), intent(in) :: bij
      real(8), intent(in) :: qq
      real(8), intent(in) :: r

      real(8) :: utot

      real(8) :: rinv


      utot = 0.0d0
      rinv = 1.0d0 / r

      utot = qq * rinv + (aij * rinv**12 - bij * rinv**6)

    end function utot 

end module mod_potential
!=======================================================================
