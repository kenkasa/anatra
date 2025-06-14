!=======================================================================
module mod_energy
!=======================================================================
  use mod_util
  use mod_const
  use mod_anatra_ermod
  use mod_traj
  use mod_com
  use mod_potential

  ! constants
  !

  ! subroutines
  !
  public :: anatra_energy
  public :: get_dist
  public :: get_cosz

  contains
!-----------------------------------------------------------------------
    subroutine anatra_energy(istep, option, traj, pot, com, uij, ljval, elval)
!-----------------------------------------------------------------------
      implicit none

      integer,        intent(in)    :: istep
      type(s_option), intent(in)    :: option
      type(s_traj),   intent(in)    :: traj(2)
      type(s_pot),    intent(inout) :: pot
      type(s_com),    intent(in)    :: com(2)
      real(8),        intent(out)   :: uij(:, :)
      real(8),        intent(out)   :: ljval(:, :)
      real(8),        intent(out)   :: elval(:, :)

      integer :: itraj, imol, jmol
      integer :: nmol(2)
      real(8) :: rljcut2


      nmol(1) = com(1)%nmol
      nmol(2) = com(2)%nmol

      if (.not. option%calc_vdw) then
        pot%acoef  = 0.0d0
        pot%bcoef  = 0.0d0
        pot%rwell2 = 0.0d0
        pot%uljmin = 0.0d0
        pot%ljsgm  = 0.0d0
        pot%ljeps  = 0.0d0
      end if

      ! Calculate Square of Cutoff Distance 
      !
      rljcut2 = option%rljcut * option%rljcut

      ! Calculate Pair Energy & Force
      !
      uij   = 0.0d0
      ljval = 0.0d0
      elval = 0.0d0

      if (option%vdw == VdwTypeSTANDARD) then
        call calc_uij(istep, option, com, traj, pot, rljcut2,   &
                      uij, ljval, elval)
      
      else if (option%vdw == VdwTypeATTRACTIVE) then
        call calc_uij_attractive(istep, option, com, traj, pot, &
                                 rljcut2, uij, ljval)
      
      else if (option%vdw == VdwTypeREPULSIVE) then
        call calc_uij_repulsive(istep, option, com, traj, pot,  &
                                rljcut2, uij, ljval)
      else if (option%vdw == VdwTypeMINDIST) then
        call calc_uij_mindist(istep, option, com, traj, pot,    &
                              rljcut2, uij, ljval)
      end if
         
    end subroutine anatra_energy
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine calc_uij(istep, option, com, traj, pot, rljcut2,    &
                        uij, ljval, elval) 
!-----------------------------------------------------------------------
      implicit none

      integer,        intent(in)    :: istep
      type(s_option), intent(in)    :: option
      type(s_com),    intent(in)    :: com(2)
      type(s_traj),   intent(in)    :: traj(2)
      type(s_pot),    intent(in)    :: pot 
      real(8),        intent(in)    :: rljcut2
      real(8),        intent(inout) :: uij(:, :)
      real(8),        intent(inout) :: ljval(:, :)
      real(8),        intent(inout) :: elval(:, :)

      integer :: iatm, jatm, imol, jmol
      real(8) :: d(3), r, r2, r2inv, r4inv, r6inv
      real(8) :: aij, bij, u, ulj, uele, fcoef


      uij   = 0.0d0
      ljval = 0.0d0
      elval = 0.0d0

      do iatm = 1, traj(1)%natm
        imol = com(1)%molid(iatm)
        do jatm = 1, traj(2)%natm 
          jmol            = com(2)%molid(jatm)
          d(1:3)          = traj(2)%coord(1:3, jatm, istep) &
                          - traj(1)%coord(1:3, iatm, istep)

          if (option%pbc) then
            d(1:3)          = d(1:3) - traj(1)%box(1:3, istep) &
                            * anint(d(1:3) / traj(1)%box(1:3, istep))
          end if

          r2              = dot_product(d, d)
          r2inv           = 1.0d0 / r2
          r4inv           = r2inv * r2inv
          r6inv           = r2inv * r4inv
      
          aij             = pot%acoef(jatm, iatm)
          bij             = pot%bcoef(jatm, iatm)

          if (r2 < rljcut2) then 
            u       = r6inv * (aij*r6inv - bij)
          else
            u       = 0.0d0
          end if
          
          ulj                  = u
          uij(jmol, imol)      = uij(jmol, imol)      + u
          ljval(jmol, imol)    = ljval(jmol, imol)    + u
      
          ! electrostatic
          !
          if (option%calc_elec) then
            r = sqrt(r2)
            if (r < option%relcut) then
              u       = pot%qq(jatm, iatm) / r
            else
              u       = 0.0d0
            end if
            uele                 = u
            uij(jmol, imol)      = uij(jmol, imol)      + u
            elval(jmol, imol)    = elval(jmol, imol)    + u
          end if

          ! tot-attractive
          !
          !if (option%tottype == TotTypeATTRACTIVE) then
          !  if (r2 <= pot%rtotwell2(jatm, iatm)) then
          !    uij(jmol, imol)      = uij(jmol, imol) - ulj - uele & 
          !                         + pot%utotmin(jatm, iatm)
          !  end if
          !end if

        end do
      end do

    end subroutine calc_uij
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine calc_uij_attractive(istep, option, com, traj, pot,  &
                                   rljcut2, uij, ljval)
!-----------------------------------------------------------------------
      implicit none

      integer,        intent(in)    :: istep
      type(s_option), intent(in)    :: option
      type(s_com),    intent(in)    :: com(2)
      type(s_traj),   intent(in)    :: traj(2)
      type(s_pot),    intent(in)    :: pot 
      real(8),        intent(in)    :: rljcut2
      real(8),        intent(inout) :: uij(:, :)
      real(8),        intent(inout) :: ljval(:, :)

      integer :: iatm, jatm, imol, jmol
      real(8) :: d(3), r, r2, r2inv, r4inv, r6inv
      real(8) :: aij, bij, u, ulj, fcoef


      uij   = 0.0d0
      ljval = 0.0d0

      do iatm = 1, traj(1)%natm
        imol = com(1)%molid(iatm)
        do jatm = 1, traj(2)%natm 
          jmol            = com(2)%molid(jatm)
          d(1:3)          = traj(2)%coord(1:3, jatm, istep) &
                          - traj(1)%coord(1:3, iatm, istep) 
          if (option%pbc) then
            d(1:3)          = d(1:3) - traj(1)%box(1:3, istep) &
                            * anint(d(1:3) / traj(1)%box(1:3, istep))
          end if
      
          r2              = dot_product(d, d)
          r2inv           = 1.0d0 / r2
          r4inv           = r2inv * r2inv
          r6inv           = r2inv * r4inv
      
          aij             = pot%acoef(jatm, iatm)
          bij             = pot%bcoef(jatm, iatm)

          if (r2 < rljcut2) then 
            u       = r6inv * (aij*r6inv - bij)
            ulj     = u
            if (r2 < pot%rwell2(jatm, iatm)) then
              u       = pot%uljmin(jatm, iatm)
            end if
          else
            u       = 0.0d0
            ulj     = 0.0d0
          end if
      
          uij(jmol, imol)      = uij(jmol, imol)   + u
          ljval(jmol, imol)    = ljval(jmol, imol) + ulj
      
          ! electrostatic
          !
          if (option%calc_elec) then
            r = sqrt(r2)
            if (r < option%relcut) then
              u       = pot%qq(jatm, iatm) / r
            else
              u       = 0.0d0
            end if
            uij(jmol, imol)      = uij(jmol, imol) + u
          end if

        end do
      end do

    end subroutine calc_uij_attractive
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine calc_uij_repulsive(istep, option, com, traj, pot,   &
                                  rljcut2, uij, ljval)
!-----------------------------------------------------------------------
      implicit none

      integer,        intent(in)    :: istep
      type(s_option), intent(in)    :: option
      type(s_com),    intent(in)    :: com(2)
      type(s_traj),   intent(in)    :: traj(2)
      type(s_pot),    intent(in)    :: pot 
      real(8),        intent(in)    :: rljcut2
      real(8),        intent(inout) :: uij(:, :)
      real(8),        intent(inout) :: ljval(:, :)

      integer :: iatm, jatm, imol, jmol
      real(8) :: d(3), r, r2, r2inv, r4inv, r6inv
      real(8) :: aij, bij, u, ulj, fcoef


      uij   = 0.0d0
      ljval = 0.0d0

      do iatm = 1, traj(1)%natm
        imol = com(1)%molid(iatm)
        do jatm = 1, traj(2)%natm 
          jmol            = com(2)%molid(jatm)
          d(1:3)          = traj(2)%coord(1:3, jatm, istep) &
                          - traj(1)%coord(1:3, iatm, istep) 
          if (option%pbc) then
            d(1:3)          = d(1:3) - traj(1)%box(1:3, istep) &
                            * anint(d(1:3) / traj(1)%box(1:3, istep))
          end if
      
          r2              = dot_product(d, d)
          r2inv           = 1.0d0 / r2
          r4inv           = r2inv * r2inv
          r6inv           = r2inv * r4inv
      
          aij             = pot%acoef(jatm, iatm)
          bij             = pot%bcoef(jatm, iatm)

          if (r2 < rljcut2) then 
            u       = r6inv * (aij*r6inv - bij)
            ulj     = u
            if (r2 < pot%rwell2(jatm, iatm)) then
              u = u - pot%uljmin(jatm, iatm)
            else
              u       = 0.0d0
            end if
          else
            u       = 0.0d0
            ulj     = 0.0d0
          end if
      
          uij(jmol, imol)      = uij(jmol, imol)      + u
          ljval(jmol, imol)    = ljval(jmol, imol)    + ulj
      
          ! electrostatic
          !
          if (option%calc_elec) then
            r = sqrt(r2)
            if (r < option%relcut) then
              u       = pot%qq(jatm, iatm) / r
            else
              u       = 0.0d0
            end if
            uij(jmol, imol)      = uij(jmol, imol)      + u
          end if

        end do
      end do

    end subroutine calc_uij_repulsive
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine calc_uij_mindist(istep, option, com, traj, pot,   &
                                rljcut2, uij, ljval)
!-----------------------------------------------------------------------
      implicit none

      integer,        intent(in)    :: istep
      type(s_option), intent(in)    :: option
      type(s_com),    intent(in)    :: com(2)
      type(s_traj),   intent(in)    :: traj(2)
      type(s_pot),    intent(in)    :: pot 
      real(8),        intent(in)    :: rljcut2
      real(8),        intent(inout) :: uij(:, :)
      real(8),        intent(inout) :: ljval(:, :)

      integer :: iatm, jatm, imol, jmol
      real(8) :: d(3), r, r2, r2inv, r4inv, r6inv


      uij   = 0.0d0
      ljval = 0.0d0

      uij = 1.0d10

      do iatm = 1, traj(1)%natm
        imol = com(1)%molid(iatm)
        do jatm = 1, traj(2)%natm 
          jmol   = com(2)%molid(jatm)
          d(1:3) = traj(2)%coord(1:3, jatm, istep) &
                 - traj(1)%coord(1:3, iatm, istep) 
          if (option%pbc) then
            d(1:3) = d(1:3) - traj(1)%box(1:3, istep) &
                   * anint(d(1:3) / traj(1)%box(1:3, istep))
          end if
      
          r2 = dot_product(d, d)
          r  = sqrt(r2)
     
          if (r < uij(jmol, imol)) &
            uij(jmol, imol) = r

        end do
      end do

    end subroutine calc_uij_mindist
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_dist(istep, pbc, box, nmol, com, dist)
!-----------------------------------------------------------------------
      implicit none

      integer,      intent(in)  :: istep
      logical,      intent(in)  :: pbc
      real(8),      intent(in)  :: box(3)
      integer,      intent(in)  :: nmol(2)
      type(s_com),  intent(in)  :: com(2)
      real(8),      intent(out) :: dist(:,:)

      integer :: imol, jmol
      real(8) :: d(3) 

     
      dist = 0.0d0

      if(pbc) then
        do imol = 1, nmol(1)
          do jmol = 1, nmol(2)
            d(1:3) = com(2)%coord(1:3, jmol, istep) &
                   - com(1)%coord(1:3, imol, istep)
            d(1:3) = d(1:3) - box(1:3) * nint(d(1:3) / box(1:3))

            dist(jmol, imol)  = sqrt(dot_product(d, d))
          end do
        end do
      else
        do imol = 1, nmol(1)
          do jmol = 1, nmol(2)
            d(1:3) = com(2)%coord(1:3, jmol, istep) &
                   - com(1)%coord(1:3, imol, istep)

            dist(jmol, imol) = sqrt(dot_product(d, d))
          end do
        end do
      end if

!
    end subroutine get_dist
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_cosz(istep, nmol, com, cosz)
!-----------------------------------------------------------------------
      implicit none

      integer,      intent(in)  :: istep
      integer,      intent(in)  :: nmol(2)
      type(s_com),  intent(in)  :: com(2)
      real(8),      intent(out) :: cosz(:,:)

      integer :: imol, jmol
      real(8) :: d(3), dist 

     
      cosz = 0.0d0

      do imol = 1, nmol(1)
        do jmol = 1, nmol(2)
          d(1:3) = com(2)%coord(1:3, jmol, istep) &
                 - com(1)%coord(1:3, imol, istep)

          dist             = sqrt(dot_product(d, d))
          cosz(jmol, imol) = d(3) / dist
        end do
      end do

!
    end subroutine get_cosz
!-----------------------------------------------------------------------

end module mod_energy
!=======================================================================
