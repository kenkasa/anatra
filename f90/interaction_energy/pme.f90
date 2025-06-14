!=======================================================================
module mod_pme
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_grid3d
  use mod_ctrl
  use mod_traj
  use mod_prmtop
  use mod_com
  use mod_potential
  use mod_pme_str
  use mod_fftmkl
  !use mod_fftw3i

  ! constants
  !

  ! structures
  !

  ! subroutines
  !
  public  :: calc_pme_short 
  public  :: calc_pme_long 
  public  :: setup_pme_gamma
  public  :: setup_pme_dfunc
  public  :: setup_pme_qfunc
  public  :: setup_pme_gfunc
  private :: bfunc
  private :: vec_prod

  contains
!-----------------------------------------------------------------------
    subroutine calc_pme_short(istep,      &
                              option,     &
                              com,        &
                              traj,       &
                              pot,        &
                              rljcut2,    &
                              uLJ,        &
                              uES,        &
                              uCORR,      &
                              frc,        &
                              trq,        &
                              pot_dual,   &
                              uLJ_dual,   &
                              uES_dual,   &
                              uCORR_dual, &
                              frc_dual,   &
                              trq_dual)
!-----------------------------------------------------------------------
      implicit none

      integer,        intent(in)    :: istep
      type(s_option), intent(in)    :: option
      type(s_com),    intent(in)    :: com(2)
      type(s_traj),   intent(in)    :: traj(2)
      type(s_pot),    intent(in)    :: pot
      real(8),        intent(in)    :: rljcut2
      real(8),        intent(inout) :: uLJ(:, :)
      real(8),        intent(inout) :: uES(:, :)
      real(8),        intent(inout) :: uCORR(:, :)
      real(8),        intent(inout) :: frc(:, :, :)
      real(8),        intent(inout) :: trq(:, :, :)
      ! for dual
      type(s_pot),    intent(in)    :: pot_dual
      real(8),        intent(inout) :: uLJ_dual  (:, :)
      real(8),        intent(inout) :: uES_dual  (:, :)
      real(8),        intent(inout) :: uCORR_dual(:, :)
      real(8),        intent(inout) :: frc_dual  (:, :, :)
      real(8),        intent(inout) :: trq_dual  (:, :, :)

      integer :: iatm, jatm, imol, jmol
      real(8) :: d(3), dmol(3), r, r2, r2inv, r4inv, r6inv
      real(8) :: aij, bij, u, usum, uele, fcoef, fj(3)
      real(8) :: aij_dual, bij_dual, u_dual, fj_dual(3)
      real(8) :: vol


      vol = traj(1)%box(1, istep) * traj(1)%box(2, istep) * traj(1)%box(3, istep)

      uLJ        = 0.0d0
      uES        = 0.0d0
      uCORR      = 0.0d0

      uLJ_dual   = 0.0d0
      uES_dual   = 0.0d0
      uCORR_dual = 0.0d0

      !$omp parallel private(iatm, jatm, imol, jmol, d, dmol, r, r2, r2inv), &
      !$omp          private(r4inv, r6inv, aij, bij, u, fcoef, fj),          &
      !$omp          private(aij_dual, bij_dual, u_dual, fj_dual),           &
      !$omp          shared(istep, vol), reduction(+:uLJ),                   &
      !$omp          reduction(+:uES), reduction(+:uCORR),                   &
      !$omp          reduction(+:frc), reduction(+:trq),                     &
      !$omp          reduction(+:uLJ_dual),   reduction(+:uES_dual),         &
      !$omp          reduction(+:uCORR_dual), reduction(+:frc_dual),         &
      !$omp          reduction(+:trq_dual)

      !$omp do
      do iatm = 1, traj(1)%natm
        imol = com(1)%molid(iatm)
        do jatm = 1, traj(2)%natm 
          jmol            = com(2)%molid(jatm)
          d(1:3)          = traj(2)%coord(1:3, jatm, istep) &
                          - traj(1)%coord(1:3, iatm, istep) 
          d(1:3)          = d(1:3) - traj(1)%box(1:3, istep) &
                          * anint(d(1:3) / traj(1)%box(1:3, istep))

          dmol(1:3)       = traj(2)%coord(1:3, jatm, istep) &
                          - com(2)%coord(1:3, jmol, istep)

          r2              = dot_product(d, d)

          ! Calculate LJ interaction
          !
          r2inv           = 1.0d0 / r2
          r4inv           = r2inv * r2inv
          r6inv           = r2inv * r4inv
      
          aij             = pot%acoef(jatm, iatm)
          bij             = pot%bcoef(jatm, iatm)

          aij_dual        = pot_dual%acoef(jatm, iatm)
          bij_dual        = pot_dual%bcoef(jatm, iatm)

          if (r2 < rljcut2) then 
            u            = r6inv * (aij*r6inv - bij)
            fcoef        = r6inv * (12.0d0 * aij * r6inv - 6.0d0 * bij)
            fj(1:3)      = fcoef * r2inv * d(1:3)
                         
            u_dual       = r6inv * (aij_dual * r6inv - bij_dual)
            fcoef        = r6inv * (12.0d0 * aij_dual * r6inv    &
                           - 6.0d0 * bij_dual)
            fj_dual(1:3) = fcoef * r2inv * d(1:3)
          else
            u            = 0.0d0
            u_dual       = 0.0d0
            fj(1:3)      = 0.0d0
            fj_dual(1:3) = 0.0d0
          end if

          uLJ(jmol, imol)           = uLJ(jmol, imol)      + u
          frc(1:3, jmol, imol)      = frc(1:3, jmol, imol) + fj(1:3) 
          
          uLJ_dual(jmol, imol)      = uLJ_dual(jmol, imol) + u_dual
          frc_dual(1:3, jmol, imol) = frc_dual(1:3, jmol, imol) &
                                      + fj_dual(1:3)

          trq(1:3, jmol, imol)      = trq(1:3, jmol, imol)          &
                                    + vec_prod(dmol(1:3), fj(1:3))
          trq_dual(1:3, jmol, imol) = trq_dual(1:3, jmol, imol)     &
                                    + vec_prod(dmol(1:3), fj_dual(1:3))
     
          ! Calculate short range part of electrostatic
          !
          r = sqrt(r2)
          if (r < option%relcut) then
            u       = pot%qq(jatm, iatm) / r * erfc(option%pme_alpha * r)
            fj(1:3) = u * r2inv * d(1:3)
            fj(1:3) = fj(1:3) + (2.0d0 * option%pme_alpha / sqrt(PI)) &
                              * exp(-option%pme_alpha**2 * r2)        & 
                              * pot%qq(jatm, iatm) * r2inv * d(1:3)

            u_dual       = pot_dual%qq(jatm, iatm) / r * erfc(option%pme_alpha * r)
            fj_dual(1:3) = u_dual * r2inv * d(1:3)
            fj_dual(1:3) = fj_dual(1:3) + (2.0d0 * option%pme_alpha / sqrt(PI)) &
                              * exp(-option%pme_alpha**2 * r2)                  & 
                              * pot_dual%qq(jatm, iatm) * r2inv * d(1:3)
          else
            u            = 0.0d0
            u_dual       = 0.0d0
            fj(1:3)      = 0.0d0
            fj_dual(1:3) = 0.0d0
          end if

          uES(jmol, imol)      = uES(jmol, imol) + u
          frc(1:3, jmol, imol) = frc(1:3, jmol, imol) + fj(1:3)

          uES_dual(jmol, imol)      = uES_dual(jmol, imol)      + u_dual
          frc_dual(1:3, jmol, imol) = frc_dual(1:3, jmol, imol) + fj_dual(1:3)

          trq(1:3, jmol, imol)      = trq(1:3, jmol, imol)          &
                                    + vec_prod(dmol(1:3), fj(1:3))
          trq_dual(1:3, jmol, imol) = trq_dual(1:3, jmol, imol)     &
                                    + vec_prod(dmol(1:3), fj_dual(1:3))

          ! Calculation correction part of electrostatic
          !
          u = - PI * pot%qq(jatm, iatm) / (option%pme_alpha**2 * vol)
          uCORR(jmol, imol) = uCORR(jmol, imol) + u

          u_dual = - PI * pot_dual%qq(jatm, iatm) / (option%pme_alpha**2 * vol)
          uCORR_dual(jmol, imol) = uCORR_dual(jmol, imol) + u_dual 

        end do
      end do
      !$omp end do
      !$omp end parallel

    end subroutine calc_pme_short
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine calc_pme_long_upara(istep,     &
                                   is_first,  &
                                   option,    &
                                   com,       &
                                   traj,      &
                                   pot,       &
                                   pmevar,    &
                                   fftinfo,   &
                                   uEL,       &
                                   frc,       &
                                   trq,       &
                                   pot_dual,  &
                                   uEL_dual,  &
                                   frc_dual,  &
                                   trq_dual)
!-----------------------------------------------------------------------
      implicit none

      integer,             intent(in)    :: istep
      logical,             intent(in)    :: is_first
      type(s_option),      intent(in)    :: option
      type(s_com),         intent(in)    :: com(2)
      type(s_traj),        intent(in)    :: traj(2)
      type(s_pot),         intent(in)    :: pot
      type(s_pmevar),      intent(inout) :: pmevar
      type(s_fftinfo),     intent(inout) :: fftinfo 
      real(8),             intent(inout) :: uEL(:, :) 
      real(8),             intent(inout) :: frc(:, :, :) 
      real(8),             intent(inout) :: trq(:, :, :) 
      type(s_pot),         intent(in)    :: pot_dual
      real(8),             intent(inout) :: uEL_dual(:, :) 
      real(8),             intent(inout) :: frc_dual(:, :, :) 
      real(8),             intent(inout) :: trq_dual(:, :, :) 

      integer :: igx, igy, igz
      integer :: imol, jmol, igr, ngx, ngy, ngz, ngr
      real(8) :: dmol(3), box(3), dbox(3), fj(3), fj_dual(3)
      real(8) :: db, eneval, eneval_dual
      real(8) :: ucorr_ij

      type(s_pmevar), save, allocatable ::  pmevars(:)
      real(8), allocatable :: qfunck     (:, :, :)
      real(8), allocatable :: qderiv     (:, :, :, :)


      ngx      = option%pme_grids(1)
      ngy      = option%pme_grids(2)
      ngz      = option%pme_grids(3)
      ngr      = ngx * ngy * ngz

      box(1:3) = traj(1)%box(1:3, istep)

      !if (.not. allocated(qfunck)) then
      !  allocate(qfunck(0:ngx-1, 0:ngy-1, 0:ngz-1))
      !  allocate(qderiv(1:3, 0:ngx-1, 0:ngy-1, 0:ngz-1))
      !end if


      if (.not. allocated(pmevars)) then
        allocate(pmevars(1:com(1)%nmol))
      end if
      
      ! generate reciprocal grid information
      !
      do imol = 1, com(1)%nmol
        call setup_pmevar(option%pme_grids, box, pmevars(imol))
      end do
      call setup_pmevar(option%pme_grids, box, pmevar)

      ! setup FFT variables 
      !
      fftinfo%ng3(1:3) = pmevars(1)%ng3(1:3)

      if (is_first) &
        call fftmkl_init(fftinfo, pmevar%qfunck, pmevar%qfuncm)

      ! generate PME functions 
      !
      do imol = 1, com(1)%nmol
        call setup_pme_splfc(option, pmevars(imol))
        call setup_pme_dfunc(option, pmevars(imol))
      end do

      ! calculate solute-solvent interactions
      !
      uEL      = 0.0d0
      uEL_dual = 0.0d0

     !$omp parallel private(imol, jmol, eneval, eneval_dual, igx, igy, igz,   &
     !$omp                  fj, fj_dual),                             &
     !$omp          shared (istep, pot, pot_dual, com, traj, ngx, ngy, ngz,   &
     !$omp                  pmevars, uEL, uEL_dual, trq, trq_dual),                    & 
     !$omp          default(shared)
     !$omp do
      do imol = 1, com(1)%nmol

        if (is_first) then

          ! generate Qfunc for solute 
          !
          call setup_pme_qfunc(istep,     &
                               imol,      &
                               option,    &
                               pot,       &
                               pot%q1,    &
                               com(1),    &
                               traj(1),   &
                               pmevars(imol),  &
                               pot_dual%q1)
        
          ! convert Q-function to G-function for solute molecule
          ! (variable name is unchanged : qfunc_u)
          call setup_pme_gfunc(option, fftinfo, pmevars(imol))

        else

          if (.not. option%pme_rigid) then

            ! generate Qfunc for solute 
            !
            call setup_pme_qfunc(istep,     &
                                 imol,      &
                                 option,    &
                                 pot,       &
                                 pot%q1,    &
                                 com(1),    &
                                 traj(1),   &
                                 pmevars(imol),  &
                                 pot_dual%q1)
         
            ! convert Q-function to G-function for solute molecule
            ! (variable name is unchanged : qfunc_u)
            call setup_pme_gfunc(option, fftinfo, pmevars(imol))

          end if

        end if
       
        do jmol = 1, com(2)%nmol
       
          call setup_pme_qfunc_v(istep,       &
                                 jmol,        &
                                 option,      &
                                 pot,         &
                                 pot%q2,      &
                                 com(2),      &
                                 traj(2),     &
                                 pmevars(imol),      &
                                 eneval,      &
                                 eneval_dual, &
                                 fj,          &
                                 fj_dual,     &
                                 trq,         &
                                 trq_dual)
       
          ! pair energy
          !
          uEL     (jmol, imol) = eneval
          uEL_dual(jmol, imol) = eneval_dual
       
          ! pair force
          !
          frc     (1:3, jmol, imol) = frc     (1:3, jmol, imol) + fj(1:3)
          frc_dual(1:3, jmol, imol) = frc_dual(1:3, jmol, imol) + fj_dual(1:3)
          
        end do

      end do
      !$omp end do
      !$omp end parallel

      !deallocate(qfunck, qderiv)

    end subroutine calc_pme_long_upara 
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
    subroutine calc_pme_long_vpara(istep,     &
                                   is_first,  &
                                   option,    &
                                   com,       &
                                   traj,      &
                                   pot,       &
                                   pmevar,    &
                                   fftinfo,   &
                                   uEL,       &
                                   frc,       &
                                   trq,       &
                                   pot_dual,  &
                                   uEL_dual,  &
                                   frc_dual,  &
                                   trq_dual)
!-----------------------------------------------------------------------
      implicit none

      integer,             intent(in)    :: istep
      logical,             intent(in)    :: is_first
      type(s_option),      intent(in)    :: option
      type(s_com),         intent(in)    :: com(2)
      type(s_traj),        intent(in)    :: traj(2)
      type(s_pot),         intent(in)    :: pot
      type(s_pmevar),      intent(inout) :: pmevar
      type(s_fftinfo),     intent(inout) :: fftinfo 
      real(8),             intent(inout) :: uEL(:, :) 
      real(8),             intent(inout) :: frc(:, :, :) 
      real(8),             intent(inout) :: trq(:, :, :) 
      type(s_pot),         intent(in)    :: pot_dual
      real(8),             intent(inout) :: uEL_dual(:, :) 
      real(8),             intent(inout) :: frc_dual(:, :, :) 
      real(8),             intent(inout) :: trq_dual(:, :, :) 

      integer :: igx, igy, igz
      integer :: imol, jmol, igr, ngx, ngy, ngz, ngr
      real(8) :: dmol(3), box(3), dbox(3), fj(3), fj_dual(3)
      real(8) :: db, eneval, eneval_dual
      real(8) :: ucorr_ij

      real(8), allocatable :: qfunck     (:, :, :)
      real(8), allocatable :: qderiv     (:, :, :, :)


      ngx      = option%pme_grids(1)
      ngy      = option%pme_grids(2)
      ngz      = option%pme_grids(3)
      ngr      = ngx * ngy * ngz

      box(1:3) = traj(1)%box(1:3, istep)

      !if (.not. allocated(qfunck)) then
      !  allocate(qfunck(0:ngx-1, 0:ngy-1, 0:ngz-1))
      !  allocate(qderiv(1:3, 0:ngx-1, 0:ngy-1, 0:ngz-1))
      !end if


      ! generate reciprocal grid information
      !
      call setup_pmevar(option%pme_grids, box, pmevar)

      ! setup FFT variables 
      !
      fftinfo%ng3(1:3) = pmevar%ng3(1:3)

      if (is_first) &
        call fftmkl_init(fftinfo, pmevar%qfunck, pmevar%qfuncm)

      ! generate PME functions 
      !
      call setup_pme_splfc(option, pmevar)
      call setup_pme_dfunc(option, pmevar)

      ! calculate solute-solvent interactions
      !
      uEL      = 0.0d0
      uEL_dual = 0.0d0

      do imol = 1, com(1)%nmol

        ! generate Qfunc for solute 
        !
        call setup_pme_qfunc(istep,     &
                             imol,      &
                             option,    &
                             pot,       &
                             pot%q1,    &
                             com(1),    &
                             traj(1),   &
                             pmevar,    &
                             pot_dual%q1)
       
        ! convert Q-function to G-function for solute molecule
        ! (variable name is unchanged : qfunc_u)
        call setup_pme_gfunc(option, fftinfo, pmevar)
       
        !$omp parallel private(jmol, eneval, eneval_dual, igx, igy, igz,         &
        !$omp                  fj, fj_dual),                                     &
        !$omp          shared (istep, pot, pot_dual, com, traj, ngx, ngy, ngz,   &
        !$omp                  pmevar, uEL, uEL_dual, trq, trq_dual),            & 
        !$omp          default(shared)
        !$omp do
        do jmol = 1, com(2)%nmol
       
          call setup_pme_qfunc_v(istep,       &
                                 jmol,        &
                                 option,      &
                                 pot,         &
                                 pot%q2,      &
                                 com(2),      &
                                 traj(2),     &
                                 pmevar,      &
                                 eneval,      &
                                 eneval_dual, &
                                 fj,          &
                                 fj_dual,     &
                                 trq,         &
                                 trq_dual)
       
          ! pair energy
          !
          uEL     (jmol, imol) = eneval
          uEL_dual(jmol, imol) = eneval_dual
       
          ! pair force
          !
          frc     (1:3, jmol, imol) = frc     (1:3, jmol, imol) + fj(1:3)
          frc_dual(1:3, jmol, imol) = frc_dual(1:3, jmol, imol) + fj_dual(1:3)
          
        end do
        !$omp end do
        !$omp end parallel

      end do

      !deallocate(qfunck, qderiv)

    end subroutine calc_pme_long_vpara 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine setup_pmevar(pme_grids, box, pmevar)
!-----------------------------------------------------------------------
      implicit none

      integer,             intent(in)    :: pme_grids(3)
      real(8),             intent(in)    :: box(3)
      type(s_pmevar),      intent(inout) :: pmevar

      integer :: ng3(3)


      ng3(1:3)        = pme_grids(1:3)

      ! setup variables
      !
      pmevar%ng3(1:3) = pme_grids(1:3)
      pmevar%box(1:3) = box(1:3)
      pmevar%del(1:3) = box(1:3) / pmevar%ng3(1:3)
      pmevar%dv(1:3)  = 1.0d0 / box(1:3) 
    
      ! allocate memory
      !
      if (.not. allocated(pmevar%dfunc)) then

        ! D-function (real(k) space)
        !
        allocate(pmevar%dfunc (0:ng3(1)/2    , 0:ng3(2) - 1, 0:ng3(3) - 1))

        ! Q-function (real(k) and reciprocal(m) space)
        !
        allocate(pmevar%qfunck(0:ng3(1) - 1  , 0:ng3(2) - 1, 0:ng3(3) - 1))
        allocate(pmevar%qfuncm(0:ng3(1)/2    , 0:ng3(2) - 1, 0:ng3(3) - 1))

        allocate(pmevar%qfunck_dual(0:ng3(1) - 1  , 0:ng3(2) - 1, 0:ng3(3) - 1))
        allocate(pmevar%qfuncm_dual(0:ng3(1)/2    , 0:ng3(2) - 1, 0:ng3(3) - 1))

        ! G-function (real(k) and reciprocal(m) space)
        !
        allocate(pmevar%gfunck(0:ng3(1) - 1  , 0:ng3(2) - 1, 0:ng3(3) - 1))
        allocate(pmevar%gfuncm(0:ng3(1)/2    , 0:ng3(2) - 1, 0:ng3(3) - 1))

        allocate(pmevar%gfunck_dual(0:ng3(1) - 1  , 0:ng3(2) - 1, 0:ng3(3) - 1))
        allocate(pmevar%gfuncm_dual(0:ng3(1)/2    , 0:ng3(2) - 1, 0:ng3(3) - 1))

        ! Q-derivative
        !
        allocate(pmevar%qderiv(3, 0:ng3(1) - 1, 0:ng3(2) - 1, 0:ng3(3) - 1))
        allocate(pmevar%qderiv_dual(3, 0:ng3(1) - 1, 0:ng3(2) - 1, 0:ng3(3) - 1))

        ! B-spline function
        !
        allocate(pmevar%splfc (0:ngmax, 3))
      end if

    end subroutine setup_pmevar
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine setup_pme_splfc(option, pmevar)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_pmevar),   intent(inout) :: pmevar 

      integer :: ixyz, ig, isp 
      integer :: ng, ng3(3)
      integer :: nsp
      real(8) :: twopi, bf
      real(8) :: sink, cosk, phase
      complex :: cmpval 


      twopi      = 2.0d0 * PI
      nsp        = option%pme_spline_order 
      ng3(1:3)   = pmevar%ng3(1:3)

      do ixyz = 1, 3
        ng = ng3(ixyz)
        do ig = 0, ng - 1
          cmpval = (0.0d0, 0.0d0)
          do isp = 0, nsp - 2
            bf     = bfunc(nsp, dble(isp + 1))
            phase  = twopi * isp * ig / ng
            cosk   = bf * cos(phase)
            sink   = bf * sin(phase)
            cmpval = cmpval + cmplx(cosk, sink)
          end do
          pmevar%splfc(ig, ixyz) = real(cmpval * conjg(cmpval))
        end do
      end do

    end subroutine setup_pme_splfc
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine setup_pme_dfunc(option, pmevar)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      type(s_pmevar), intent(inout) :: pmevar

      integer     :: igx, igy, igz, lg(3)
      integer     :: ng3(3), ng, ngr, ngG
      integer     :: nsp
      real(8)     :: v(3), vol, vsq, pisq, alphasq, spl_prod
      complex(16) :: p


      pisq     = PI * PI
      alphasq  = option%pme_alpha * option%pme_alpha
      ng3      = pmevar%ng3
      vol      = pmevar%box(1) * pmevar%box(2) * pmevar%box(3)

      pmevar%dfunc = 0.0d0
      do igz = 0, ng3(3) - 1
        do igy = 0, ng3(2) - 1
          do igx = 0, ng3(1) / 2
            spl_prod  =   pmevar%splfc(igx, 1)   & 
                        * pmevar%splfc(igy, 2)   &
                        * pmevar%splfc(igz, 3)

            if (igx <= ng3(1)/2) then
              lg(1) = igx
            else
              lg(1) = igx - ng3(1)
            end if

            if (igy <= ng3(2)/2) then
              lg(2) = igy
            else
              lg(2) = igy - ng3(2)
            end if

            if (igz <= ng3(3)/2) then
              lg(3) = igz
            else
              lg(3) = igz - ng3(3)
            end if

            v(1:3)    = lg(1:3) * pmevar%dv(1:3)
            vsq       = dot_product(v, v) 

            pmevar%dfunc(igx, igy, igz) &
              = exp(- pisq * vsq / alphasq) / (PI * vol * vsq * spl_prod)
          end do
        end do
      end do

      pmevar%dfunc(0, 0, 0) = 0.0d0

    end subroutine setup_pme_dfunc
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine setup_pme_qfunc(istep,       &
                               imol,        &
                               option,      &
                               pot,         &
                               charge,      &
                               com,         &
                               traj,        &
                               pmevar,      &
                               charge_dual)

!-----------------------------------------------------------------------
      implicit none

      integer,        intent(in)    :: istep
      integer,        intent(in)    :: imol
      type(s_option), intent(in)    :: option
      type(s_pot),    intent(in)    :: pot
      real(8),        intent(in)    :: charge(:) 
      type(s_com),    intent(in)    :: com 
      type(s_traj),   intent(in)    :: traj
      type(s_pmevar), intent(inout) :: pmevar
      real(8),        intent(in)    :: charge_dual(:) 

      integer :: ista, iend, iatm
      integer :: igx, igy, igz, ig
      integer :: jgx, jgy, jgz
      integer :: ispx, ispy, ispz
      integer :: ngx, ngy, ngz, ngr
      integer :: ixyz
      integer :: norder
      integer :: ufloor(3), rci(3)
      real(8) :: q, q_dual, pos(3), u(3), du(3), inm(3)
      real(8) :: bfx, bfy, bfz, cfx, cfy, cfz, dfx, dfy, dfz
      real(8) :: box(3)


      ! setup several variables
      ista     = com%molsta(imol)
      iend     = com%molend(imol)

      ngx      = pmevar%ng3(1) 
      ngy      = pmevar%ng3(2) 
      ngz      = pmevar%ng3(3)
      ngr      = ngx * ngy * ngz

      box(1:3) = pmevar%box(1:3)


      ! spline order (currently, hard coding)
      norder   = 4
  
      ! calculate Q-function and its derivative
      !
      pmevar%qfunck = 0.0d0
      pmevar%qderiv = 0.0d0

      pmevar%qfunck_dual = 0.0d0
      pmevar%qderiv_dual = 0.0d0 ! currently, not used

      do iatm = ista, iend 
        q           = charge(iatm)
        q_dual      = charge_dual(iatm)
        pos(1:3)    = traj%coord(1:3, iatm, istep)
        pos(1:3)    = pos(1:3) - nint(pos(1:3) / box(1:3)) * box(1:3)
        !pos(1:3)    = pos(1:3) + box(1:3) * 0.5d0

        u(1:3)      = pos(1:3) / pmevar%del(1:3)
        !u(1:3)      = pos(1:3) / pmevar%box(1:3)
        ufloor(1:3) = floor(u(1:3))
        du(1:3)     = u(1:3) - ufloor(1:3)
        !du(1:3)     = u(1:3) * pmevar%ng3(1:3)
        !rci(1:3)    = int(du(1:3))

        do ispz = 0, norder - 1
          bfz = bfunc(norder,     du(3) + ispz    )
          cfz = bfunc(norder - 1, du(3) + ispz    )
          dfz = bfunc(norder - 1, du(3) + ispz - 1)
          !bfz = bfunc(norder, du(3) - rci(3) + ispz)
          igz = modulo(ufloor(3) - ispz,     ngz)
          jgz = modulo(ufloor(3) - ispz - 1, ngz)
          !igz = modulo(rci(3) - ispz, ngz) 
          do ispy = 0, norder - 1
            bfy = bfunc(norder    , du(2) + ispy    )
            cfy = bfunc(norder - 1, du(2) + ispy    ) 
            dfy = bfunc(norder - 1, du(2) + ispy - 1) 
            !bfy = bfunc(norder, du(2) - rci(2) + ispy)
            igy = modulo(ufloor(2) - ispy    , ngy)
            jgy = modulo(ufloor(2) - ispy - 1, ngy)
            !igy = modulo(rci(2) - ispy, ngy) 
            do ispx = 0, norder - 1
              bfx = bfunc(norder,     du(1) + ispx    )
              cfx = bfunc(norder - 1, du(1) + ispx    )
              dfx = bfunc(norder - 1, du(1) + ispx - 1)
              !bfx = bfunc(norder, du(1) - rci(1) + ispx)
              igx = modulo(ufloor(1) - ispx,     ngx)
              jgx = modulo(ufloor(1) - ispx - 1, ngx)
              !igx = modulo(rci(1) - ispx, ngx) 

              pmevar%qfunck(igx, igy, igz)     &
                = pmevar%qfunck(igx, igy, igz) &
                  + q * bfx * bfy * bfz

              pmevar%qfunck_dual(igx, igy, igz)       &
                = pmevar%qfunck_dual(igx, igy, igz)   &
                  + q_dual * bfx * bfy * bfz
              
              ! derivative
              !
              if (ispx <= norder - 2) then
                pmevar%qderiv(1, igx, igy, igz)             &
                  = pmevar%qderiv(1, igx, igy, igz)         &
                    + q / pmevar%del(1)                     &
                      * (cfx - dfx) * bfy * bfz
              end if

              if (ispy <= norder - 2) then
                pmevar%qderiv(2, igx, igy, igz)             &
                  = pmevar%qderiv(2, igx, igy, igz)         &
                    + q / pmevar%del(2)                     &
                      * bfx * (cfy - dfy) * bfz
              end if

              if (ispz <= norder - 2) then
                pmevar%qderiv(3, igx, igy, igz)             &
                  = pmevar%qderiv(3, igx, igy, igz)         &
                    + q / pmevar%del(3)                     &
                      * bfx * bfy * (cfz - dfz)
              end if

            end do
          end do
        end do
      end do

    end subroutine setup_pme_qfunc
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine setup_pme_qfunc_v(istep,       &
                                 imol,        &
                                 option,      &
                                 pot,         &
                                 charge,      &
                                 com,         &
                                 traj,        &
                                 pmevar,      &
                                 eneval,      &
                                 eneval_dual, &
                                 fj,          &
                                 fj_dual,     &
                                 trq,         &
                                 trq_dual)
!-----------------------------------------------------------------------
      implicit none

      integer,        intent(in)    :: istep
      integer,        intent(in)    :: imol
      type(s_option), intent(in)    :: option
      type(s_pot),    intent(in)    :: pot
      real(8),        intent(in)    :: charge(:) 
      type(s_com),    intent(in)    :: com 
      type(s_traj),   intent(in)    :: traj
      type(s_pmevar), intent(in)    :: pmevar
      real(8),        intent(out)   :: eneval
      real(8),        intent(out)   :: eneval_dual
      real(8),        intent(out)   :: fj(3)
      real(8),        intent(out)   :: fj_dual(3)
      real(8),        intent(inout) :: trq(:, :, :)
      real(8),        intent(inout) :: trq_dual(:, :, :)

      integer :: ista, iend, iatm
      integer :: igx, igy, igz, ig
      integer :: jgx, jgy, jgz
      integer :: ispx, ispy, ispz
      integer :: ngx, ngy, ngz, ngr
      integer :: ixyz
      integer :: norder
      integer :: ufloor(3), rci(3)
      real(8) :: q, q_dual, pos(3), u(3), du(3), inm(3), dmol(3)
      real(8) :: bfx, bfy, bfz, cfx, cfy, cfz, dfx, dfy, dfz
      real(8) :: qf, deriv
      real(8) :: box(3)
      real(8) :: fj_atom(3), fj_atom_dual(3), fj_tmp(3), fj_dual_tmp(3)


      ! setup several variables
      ista     = com%molsta(imol)
      iend     = com%molend(imol)

      ngx      = pmevar%ng3(1) 
      ngy      = pmevar%ng3(2) 
      ngz      = pmevar%ng3(3)
      ngr      = ngx * ngy * ngz

      box(1:3) = pmevar%box(1:3)


      ! spline order (currently, hard coding)
      norder   = 4
  
      ! calculate Q-function and its derivative
      !
      eneval      = 0.0d0
      eneval_dual = 0.0d0
      fj          = 0.0d0
      fj_dual     = 0.0d0
      do iatm = ista, iend
        fj_atom      = 0.0d0
        fj_atom_dual = 0.0d0

        q           = charge(iatm)
        pos(1:3)    = traj%coord(1:3, iatm, istep)
        pos(1:3)    = pos(1:3) - nint(pos(1:3) / box(1:3)) * box(1:3)

        dmol(1:3)   = traj%coord(1:3, iatm, istep) &
                    - com%coord(1:3, imol, istep)  

        !pos(1:3)    = pos(1:3) + box(1:3) * 0.5d0

        u(1:3)      = pos(1:3) / pmevar%del(1:3)
        !u(1:3)      = pos(1:3) / pmevar%box(1:3)
        ufloor(1:3) = floor(u(1:3))
        du(1:3)     = u(1:3) - ufloor(1:3)
        !du(1:3)     = u(1:3) * pmevar%ng3(1:3)
        !rci(1:3)    = int(du(1:3))

        do ispz = 0, norder - 1
          bfz = bfunc(norder,     du(3) + ispz    )
          cfz = bfunc(norder - 1, du(3) + ispz    )
          dfz = bfunc(norder - 1, du(3) + ispz - 1)
          !bfz = bfunc(norder, du(3) - rci(3) + ispz)
          igz = modulo(ufloor(3) - ispz,     ngz)
          !igz = modulo(rci(3) - ispz, ngz) 
          do ispy = 0, norder - 1
            bfy = bfunc(norder    , du(2) + ispy    )
            cfy = bfunc(norder - 1, du(2) + ispy    ) 
            dfy = bfunc(norder - 1, du(2) + ispy - 1) 
            !bfy = bfunc(norder, du(2) - rci(2) + ispy)
            igy = modulo(ufloor(2) - ispy    , ngy)
            !igy = modulo(rci(2) - ispy, ngy) 
            do ispx = 0, norder - 1
              bfx = bfunc(norder,     du(1) + ispx    )
              cfx = bfunc(norder - 1, du(1) + ispx    )
              dfx = bfunc(norder - 1, du(1) + ispx - 1)
              !bfx = bfunc(norder, du(1) - rci(1) + ispx)
              igx = modulo(ufloor(1) - ispx,     ngx)
              !igx = modulo(rci(1) - ispx, ngx) 

              qf          = q * bfx * bfy * bfz 
              eneval      = eneval      &
                          + pmevar%gfunck(igx, igy, igz) * qf 
              eneval_dual = eneval_dual &
                          + pmevar%gfunck_dual(igx, igy, igz) * qf 
              
              ! derivative
              !
              !   x-direction 
              if (ispx <= norder - 2) then
                deriv      = q / pmevar%del(1) * cfx * bfy * bfz
                fj_tmp(1)       = - pmevar%gfunck(igx, igy, igz)      * deriv
                fj_dual_tmp(1)  = - pmevar%gfunck_dual(igx, igy, igz) * deriv

                fj_atom(1)      = fj_atom(1) + fj_tmp(1)
                fj_atom_dual(1) = fj_atom_dual(1) + fj_dual_tmp(1)

                fj(1)      = fj(1)      + fj_tmp(1)
                fj_dual(1) = fj_dual(1) + fj_dual_tmp(1)
                !fj(1)      = fj(1)      - pmevar%gfunck(igx, igy, igz)      * deriv
                !fj_dual(1) = fj_dual(1) - pmevar%gfunck_dual(igx, igy, igz) * deriv
              end if

              if (ispx > 0) then
                deriv      = - q / pmevar%del(1) * dfx * bfy * bfz

                fj_tmp(1)      = - pmevar%gfunck(igx, igy, igz)      * deriv
                fj_dual_tmp(1) = - pmevar%gfunck_dual(igx, igy, igz) * deriv

                fj_atom(1)      = fj_atom(1) + fj_tmp(1)
                fj_atom_dual(1) = fj_atom_dual(1) + fj_dual_tmp(1)

                fj(1)      = fj(1)      + fj_tmp(1)
                fj_dual(1) = fj_dual(1) + fj_dual_tmp(1)
                !fj(1)      = fj(1)      - pmevar%gfunck(igx, igy, igz)      * deriv
                !fj_dual(1) = fj_dual(1) - pmevar%gfunck_dual(igx, igy, igz) * deriv
              end if

              !   y-direction
              if (ispy <= norder - 2) then
                deriv      = q / pmevar%del(2) * bfx * cfy * bfz

                fj_tmp(2)      = - pmevar%gfunck(igx, igy, igz)      * deriv
                fj_dual_tmp(2) = - pmevar%gfunck_dual(igx, igy, igz) * deriv

                fj_atom(2)      = fj_atom(2) + fj_tmp(2)
                fj_atom_dual(2) = fj_atom_dual(2) + fj_dual_tmp(2)

                fj(2)      = fj(2)      + fj_tmp(2)
                fj_dual(2) = fj_dual(2) + fj_dual_tmp(2)
                !fj(2)      = fj(2)      - pmevar%gfunck(igx, igy, igz)      * deriv
                !fj_dual(2) = fj_dual(2) - pmevar%gfunck_dual(igx, igy, igz) * deriv
              end if

              if (ispy > 0) then
                deriv      = - q / pmevar%del(2) * bfx * dfy * bfz

                fj_tmp(2)      = - pmevar%gfunck(igx, igy, igz)      * deriv
                fj_dual_tmp(2) = - pmevar%gfunck_dual(igx, igy, igz) * deriv

                fj_atom(2)      = fj_atom(2) + fj_tmp(2)
                fj_atom_dual(2) = fj_atom_dual(2) + fj_dual_tmp(2)

                fj(2)      = fj(2)      + fj_tmp(2)
                fj_dual(2) = fj_dual(2) + fj_dual_tmp(2)
                !fj(2)      = fj(2)      - pmevar%gfunck(igx, igy, igz)      * deriv
                !fj_dual(2) = fj_dual(2) - pmevar%gfunck_dual(igx, igy, igz) * deriv
              end if

              !   z-direction
              if (ispz <= norder - 2) then
                deriv      = q / pmevar%del(3) * bfx * bfy * cfz

                fj_tmp(3)      = - pmevar%gfunck(igx, igy, igz)      * deriv
                fj_dual_tmp(3) = - pmevar%gfunck_dual(igx, igy, igz) * deriv

                fj_atom(3)      = fj_atom(3) + fj_tmp(3)
                fj_atom_dual(3) = fj_atom_dual(3) + fj_dual_tmp(3)

                fj(3)      = fj(3)      + fj_tmp(3)
                fj_dual(3) = fj_dual(3) + fj_dual_tmp(3)
                !fj(3)      = fj(3)      - pmevar%gfunck(igx, igy, igz)      * deriv
                !fj_dual(3) = fj_dual(3) - pmevar%gfunck_dual(igx, igy, igz) * deriv
              end if

              if (ispz > 0) then
                deriv      = - q / pmevar%del(3) * bfx * bfy * dfz

                fj_tmp(3)      = - pmevar%gfunck(igx, igy, igz)      * deriv
                fj_dual_tmp(3) = - pmevar%gfunck_dual(igx, igy, igz) * deriv

                fj_atom(3)      = fj_atom(3) + fj_tmp(3)
                fj_atom_dual(3) = fj_atom_dual(3) + fj_dual_tmp(3)

                fj(3)      = fj(3)      + fj_tmp(3)
                fj_dual(3) = fj_dual(3) + fj_dual_tmp(3)
                !fj(3)      = fj(3)      - pmevar%gfunck(igx, igy, igz)      * deriv
                !fj_dual(3) = fj_dual(3) - pmevar%gfunck_dual(igx, igy, igz) * deriv
              end if


              !if (ispz <= norder - 2) then
              !  deriv      = q / pmevar%del(3) * bfx * bfy * (cfz - dfz)
              !  fj(3)      = fj(3)      - pmevar%gfunck(igx, igy, igz)      * deriv
              !  fj_dual(3) = fj_dual(3) - pmevar%gfunck_dual(igx, igy, igz) * deriv
              !end if


              !if (ispx <= norder - 2) then
              !  deriv      = q / pmevar%del(1) * (cfx - dfx) * bfy * bfz
              !  fj(1)      = fj(1)      - pmevar%gfunck(igx, igy, igz)      * deriv
              !  fj_dual(1) = fj_dual(1) - pmevar%gfunck_dual(igx, igy, igz) * deriv
              !end if

              !if (ispy <= norder - 2) then
              !  deriv      = q / pmevar%del(2) * bfx * (cfy - dfy) * bfz
              !  fj(2)      = fj(2)      - pmevar%gfunck(igx, igy, igz)      * deriv
              !  fj_dual(2) = fj_dual(2) - pmevar%gfunck_dual(igx, igy, igz) * deriv
              !end if

              !if (ispz <= norder - 2) then
              !  deriv      = q / pmevar%del(3) * bfx * bfy * (cfz - dfz)
              !  fj(3)      = fj(3)      - pmevar%gfunck(igx, igy, igz)      * deriv
              !  fj_dual(3) = fj_dual(3) - pmevar%gfunck_dual(igx, igy, igz) * deriv
              !end if

            end do
          end do
        end do

        trq(1:3, imol, 1) = trq(1:3, imol, 1) &
                          + vec_prod(dmol(1:3), fj_atom(1:3))
        trq_dual(1:3, imol, 1) = trq_dual(1:3, imol, 1) &
                               + vec_prod(dmol(1:3), fj_atom_dual(1:3))

      end do

    end subroutine setup_pme_qfunc_v
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine setup_pme_gfunc(option, fftinfo, pmevar)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),  intent(in)    :: option
      type(s_fftinfo), intent(in)    :: fftinfo
      type(s_pmevar),  intent(inout) :: pmevar 

      integer :: igx, igy, igz, ig, igR, igI
      integer :: ngx, ngy, ngz, ngr, ngk


      ! setup several variables
      !
      ngx      = pmevar%ng3(1) 
      ngy      = pmevar%ng3(2) 
      ngz      = pmevar%ng3(3)

      ! FT: Q(real) => Q(recp)
      !

      !pmevar%qfuncm = 0.0d0
      call fftmkl_r2c(fftinfo, pmevar%qfunck, pmevar%qfuncm)

      if (option%pme_dual) then
        call fftmkl_r2c(fftinfo, pmevar%qfunck_dual, pmevar%qfuncm_dual)
      end if
      !call fftw3i_execute(fftinfo,                                   &
      !                    FFTW_SignForward,                          &
      !                    pmevar%qfunck,  &
      !                    pmevar%qfuncm)

      !call fftw3i_execute(fftinfo,                                   &
      !                    FFTW_SignForward,                          &
      !                    pmevar%qfunck(0:ngx-1, 0:ngy-1, 0:ngz-1),  &
      !                    pmevar%qfuncm(0:ngx/2, 0:ngy-1, 0:ngz-1))

      ! G = D * Q
      !
      pmevar%gfuncm = 0.0d0
      pmevar%gfuncm(:, :, :) = pmevar%dfunc(:, :, :)              &
                             * pmevar%qfuncm(:, :, :)
      if (option%pme_dual) then
        pmevar%gfuncm_dual = 0.0d0
        pmevar%gfuncm_dual(:, :, :) = pmevar%dfunc(:, :, :)       &
                                    * pmevar%qfuncm_dual(:, :, :)
      end if

      ! FT: G(recp) => G(real) 
      !
      pmevar%gfunck = 0.0d0
      call fftmkl_c2r(fftinfo, pmevar%gfuncm, pmevar%gfunck)

      pmevar%gfunck_dual = 0.0d0
      if (option%pme_dual) then
        call fftmkl_c2r(fftinfo, pmevar%gfuncm_dual, pmevar%gfunck_dual)
      end if

    end subroutine setup_pme_gfunc
!---------------------------------------------------------------------

!---------------------------------------------------------------------
    function bfunc(n, x)
!---------------------------------------------------------------------
      implicit none

      integer, intent(in) :: n
      real(8), intent(in) :: x

      real(8) :: bfunc

      real(8) :: x2, x3
      real(8) :: inv6


      if (x < 0.0d0 .or. x > dble(n)) then
        bfunc = 0.0d0
      else

        if (n == 2) then
          if (x < 1.0d0) then
            bfunc = x
          else if (x >= 1.0d0 .and. x <= 2.0d0) then
            bfunc = 2.0d0 - x
          end if
        else if (n == 3) then
          x2 = x * x
          if (x < 1.0d0) then
            bfunc =  0.5d0 * x2
          else if (x >= 1.0d0 .and. x <= 2.0d0) then
            bfunc = -0.5d0 * (2.0d0 * x2 - 6.0d0 * x + 3.0d0)
          else if (x >= 2.0d0 .and. x <= 3.0d0) then
            bfunc =  0.5d0 * (x2 - 6.0d0 * x + 9.0d0)
          end if
        else if (n == 4) then
          x2 = x * x
          x3 = x * x * x
          inv6 = 1.0d0 / 6.0d0
          if (x < 1.0d0) then
            bfunc =  inv6 * x3
          else if (x >= 1.0d0 .and. x <= 2.0d0) then
            bfunc = -inv6 * (3.0d0 * x3 - 12.0d0 * x2 + 12.0d0 * x - 4.0d0)
          else if (x >= 2.0d0 .and. x <= 3.0d0) then
            bfunc =  inv6 * (3.0d0 * x3 - 24.0d0 * x2 + 60.0d0 * x - 44.0d0)
          else if (x >= 3.0d0 .and. x <= 4.0d0) then
            bfunc = -inv6 * (x3 - 12.0d0 * x2 + 48.0d0 * x - 64.0d0)
          end if

        end if

      end if

!
    end function bfunc
!---------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    function vec_prod(v1, v2)
!-----------------------------------------------------------------------
      implicit none

      real(8), intent(in) :: v1(3)
      real(8), intent(in) :: v2(3)

      real(8) :: vec_prod(3)


      vec_prod(1) = v1(2) * v2(3) - v1(3) * v2(2)
      vec_prod(2) = v1(3) * v2(1) - v1(1) * v2(3)
      vec_prod(3) = v1(1) * v2(2) - v1(2) * v2(1)

    end function vec_prod
!-----------------------------------------------------------------------

end module mod_pme
!=======================================================================
