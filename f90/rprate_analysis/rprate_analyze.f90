!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_ctrl

  ! constants
  !

  ! structures
  !
  type :: s_kin
    real(8) :: kon1, kon1_err
    real(8) :: kon2, kon2_err
    real(8) :: kf1, kf1_err
    real(8) :: kf2, kf2_err
    real(8) :: kon_kf_ratio1, kon_kf_ratio1_err
    real(8) :: kon_kf_ratio2, kon_kf_ratio2_err
    real(8) :: chi1, chi1_err
    real(8) :: chi2, chi2_err
    real(8) :: dG, dG_err
  end type s_kin

  ! subroutines
  !
  public  :: analyze
  private :: calc_sum_err
  private :: calc_prod_err
  private :: calc_div_err

  contains
!-----------------------------------------------------------------------
    subroutine analyze(option)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),  intent(inout) :: option

      type(s_kin) :: kin


      ! Kstar calculation
      !
      if (option%use_sfe) then
        call get_Kstar(option, kin) 
      end if

      ! kon calculation
      !
      call get_kon(option, kin)

      ! kf calculation
      !
      call get_kf(option, kin)

      ! kinetic contribution calculation
      !
      call get_chi(option, kin)

      ! Output
      !
      write(iw,*)
      write(iw,'("Analyze> Evaluate rate constants")')
      write(iw,*)
      write(iw,'("+----- Overall rate constant -----+")')
      write(iw,'("kon   [1st-order] (s^-1 M^-1) = ", es15.7, " +- ", es15.7)') &
        kin%kon1, kin%kon1_err
      if (option%rporder == 2) then
        write(iw,'("kon   [2nd-order] (s^-1 M^-1) = ", es15.7, " +- ", es15.7)') &
          kin%kon2, kin%kon2_err
      end if

      write(iw,*)
      write(iw,'("+----- Difussion-controlled reaction rate constant -----+")')
      write(iw,'("kf    [1st-order] (s^-1 M^-1) = ", es15.7, " +- ", es15.7)') &
        kin%kf1, kin%kf1_err
      if (option%rporder == 2) then
        write(iw,'("kf    [2nd-order] (s^-1 M^-1) = ", es15.7, " +- ", es15.7)') &
          kin%kf2, kin%kf2_err
      end if

      write(iw,*)
      write(iw,'("+----- Contribution of diffusion part to overall rate -----+")')
      write(iw,'("kon / kf [1st-order]           = ", es15.7, " +- ", es15.7)') &
        kin%kon_kf_ratio1, kin%kon_kf_ratio1_err
      if (option%rporder == 2) then
        write(iw,'("kon / kf [2nd-order]           = ", es15.7, " +- ", es15.7)') &
          kin%kon_kf_ratio2, kin%kon_kf_ratio2_err
      end if
      write(iw,*)
      write(iw,'("+----- Thermodynamic (Kstar) and Kinetic (Chi) contributions -----+")')
      write(iw,'("Kstar                (M^-1)    = ", es15.7, " +- ", es15.7)') &
        option%Kstar, option%Kstar_err

      if (option%use_sfe) then
        write(iw,'("dG                 (kcal/mol)  = ", es15.7, " +- ", es15.7)') &
          kin%dG, kin%dG_err
      end if

      write(iw,'("chi      [1st-order] (s^-1)    = ", es15.7, " +- ", es15.7)') &
        kin%chi1, kin%chi1_err

      if (option%rporder == 2) then
        write(iw,'("chi      [2nd-order] (s^-1)    = ", es15.7, " +- ", es15.7)') &
          kin%chi2, kin%chi2_err
      end if

    end subroutine analyze 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_Kstar(option, kin)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(inout) :: option
      type(s_kin),    intent(inout) :: kin

      real(8) :: dG, dG1, dG_err, err1, err2
      real(8) :: Kstar, Kstar_err
      real(8) :: kT, beta


      dG1    = option%sfeb - option%sfed + option%dGcorr
      dG     = dG1 + option%shift_sol - option%shift_ref
      err1   = calc_sum_err(option%sfeb, option%sfed, option%sfeb_err, option%sfed_err)
      err2   = calc_sum_err(dG1, option%shift_sol, err1, option%shift_sol_err)
      dG_err = calc_sum_err(dG1, option%shift_ref, err2, option%shift_ref_err)
     
      kT     = Boltz * option%temperature
      beta   = 1.0d0 / kT

      Kstar     = exp(- beta * dG)
      Kstar_err = Kstar * beta * dG_err

      option%Kstar     = Kstar
      option%Kstar_err = Kstar_err

      kin%dG           = dG
      kin%dG_err       = dG_err

    end subroutine get_Kstar
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_kon(option, kin)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      type(s_kin),    intent(inout) :: kin

      real(8) :: kins_sq, kins_sq_err
      real(8) :: taud_sq, taud_sq_err

      real(8) :: first, first_err
      real(8) :: second, second_err
      real(8) :: third, third_err

      real(8) :: nume, nume_err
      real(8) :: denom1, denom1_err
      real(8) :: denom2, denom2_err


      kins_sq     = option%kins * option%kins
      kins_sq_err = calc_prod_err(option%kins, option%kins, &
                                  option%kins_err, option%kins_err)

      taud_sq     = option%taud * option%taud
      taud_sq_err = calc_prod_err(option%taud, option%taud, &
                                  option%taud_err, option%taud_err)


      first     = option%kins * option%taud
      first_err = calc_prod_err(option%kins, option%taud, &
                                option%kins_err, option%taud_err) 

      if (option%rporder == 2) then
        second     = kins_sq * option%taupa3
        second_err = calc_prod_err(kins_sq, option%taupa3, &
                                   kins_sq_err, option%taupa3_err)

        third      = kins_sq * taud_sq
        third_err  = calc_prod_err(kins_sq, taud_sq, &
                                   kins_sq_err, taud_sq_err)
      end if

      nume     = option%kins * option%Kstar
      nume_err = calc_prod_err(option%kins, option%Kstar, &
                               option%kins_err, option%Kstar_err)

      denom1     = 1.0d0 + first
      denom1_err = first_err 


      if (option%rporder == 2) then
        denom2     = 1.0d0 + first - second + third
        denom2_err = sqrt(first_err ** 2 + second_err ** 2 + third_err ** 2)
      end if

      ! 1st-order rate constant
      !
      kin%kon1     = nume / denom1
      !kin%kon1_err = sqrt((nume_err/denom1)**2 + (nume/denom1**2 * denom1_err)**2)
      kin%kon1_err = calc_div_err(nume, denom1, nume_err, denom1_err)

      kin%kon1     = kin%kon1     / option%timeunit
      kin%kon1_err = kin%kon1_err / option%timeunit

      ! 2nd-order rate constant
      !

      if (option%rporder == 2) then
        kin%kon2     = nume / denom2
        !kin%kon2_err = sqrt((nume_err/denom2)**2 + (nume/denom2**2 * denom2_err)**2)
        kin%kon2_err = calc_div_err(nume, denom2, nume_err, denom2_err)
       
        kin%kon2     = kin%kon2     / option%timeunit
        kin%kon2_err = kin%kon2_err / option%timeunit
      end if


    end subroutine get_kon
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_kf(option, kin)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      type(s_kin),    intent(inout) :: kin

      real(8) :: taud_sq, taud_sq_err
      real(8) :: dtau_sq, dtau_sq_err
      real(8) :: fact, fact_err
      real(8) :: denom, denom_err


      ! 1st-order diffusion controlled rate consntant
      !
      kin%kf1     = option%Kstar / option%taud
      kin%kf1_err = calc_div_err(option%Kstar,     option%taud, &
                                 option%Kstar_err, option%taud_err)

      kin%kf1     = kin%kf1     / option%timeunit
      kin%kf1_err = kin%kf1_err / option%timeunit

      kin%kon_kf_ratio1     = kin%kon1 / kin%kf1
      kin%kon_kf_ratio1_err = calc_div_err(kin%kon1, kin%kf1, kin%kon1_err, kin%kf1_err) 

      ! 2nd-order diffusion controlled rate constant
      !
      if (option%rporder == 2) then
        taud_sq     = option%taud * option%taud
        taud_sq_err = calc_prod_err(option%taud, option%taud, &
                                    option%taud_err, option%taud_err)

        dtau_sq     = option%taupa3 - taud_sq
        dtau_sq_err = calc_sum_err(option%taupa3, taud_sq, &
                                   option%taupa3_err, taud_sq_err)

        fact        = option%kins * dtau_sq
        fact_err    = calc_prod_err(option%kins, dtau_sq, &
                                    option%kins_err, dtau_sq_err)

        denom       = option%taud - fact
        denom_err   = calc_sum_err(option%taud, fact, &
                                   option%taud_err, fact_err)

        kin%kf2        = option%Kstar / denom
        kin%kf2_err    = calc_div_err(option%Kstar, denom, &
                                      option%Kstar_err, denom_err)

        kin%kf2        = kin%kf2     / option%timeunit
        kin%kf2_err    = kin%kf2_err / option%timeunit

        kin%kon_kf_ratio2     = kin%kon2 / kin%kf2
        kin%kon_kf_ratio2_err = calc_div_err(kin%kon2, kin%kf2, kin%kon2_err, kin%kf2_err)

      end if
      
    end subroutine get_kf
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_chi(option, kin)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      type(s_kin),    intent(inout) :: kin

      real(8) :: kins_sq, kins_sq_err
      real(8) :: taud_sq, taud_sq_err

      real(8) :: first, first_err
      real(8) :: second, second_err
      real(8) :: third, third_err

      real(8) :: nume, nume_err
      real(8) :: denom1, denom1_err
      real(8) :: denom2, denom2_err


      kins_sq     = option%kins * option%kins
      kins_sq_err = calc_prod_err(option%kins, option%kins, &
                                  option%kins_err, option%kins_err)

      taud_sq     = option%taud * option%taud
      taud_sq_err = calc_prod_err(option%taud, option%taud, &
                                  option%taud_err, option%taud_err)


      first     = option%kins * option%taud
      first_err = calc_prod_err(option%kins, option%taud, &
                                option%kins_err, option%taud_err) 

      if (option%rporder == 2) then
        second     = kins_sq * option%taupa3
        second_err = calc_prod_err(kins_sq, option%taupa3, &
                                   kins_sq_err, option%taupa3_err)

        third      = kins_sq * taud_sq
        third_err  = calc_prod_err(kins_sq, taud_sq, &
                                   kins_sq_err, taud_sq_err)
      end if

      nume     = option%kins
      nume_err = option%kins_err 

      denom1     = 1.0d0 + first
      denom1_err = first_err 


      if (option%rporder == 2) then
        denom2     = 1.0d0 + first - second + third
        denom2_err = sqrt(first_err ** 2 + second_err ** 2 + third_err ** 2)
      end if

      ! 1st-order kinetic contribution 
      !
      kin%chi1     = nume / denom1
      !kin%kon1_err = sqrt((nume_err/denom1)**2 + (nume/denom1**2 * denom1_err)**2)
      kin%chi1_err = calc_div_err(nume, denom1, nume_err, denom1_err)

      kin%chi1     = kin%chi1     / option%timeunit
      kin%chi1_err = kin%chi1_err / option%timeunit

      ! 2nd-order rate constant
      !

      if (option%rporder == 2) then
        kin%chi2     = nume / denom2
        kin%chi2_err = calc_div_err(nume, denom2, nume_err, denom2_err)
       
        kin%chi2     = kin%chi2     / option%timeunit
        kin%chi2_err = kin%chi2_err / option%timeunit
      end if


    end subroutine get_chi
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
    function calc_sum_err(val1, val2, err1, err2)
!-----------------------------------------------------------------------
      implicit none

      real(8), intent(in) :: val1, val2, err1, err2

      real(8) :: calc_sum_err, err
      real(8) :: p1, p2


      p1  = err1 
      p2  = err2
      err = sqrt(p1 * p1 + p2 * p2) 

      ! return value
      calc_sum_err = err

    end function calc_sum_err
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    function calc_prod_err(val1, val2, err1, err2)
!-----------------------------------------------------------------------
      implicit none

      real(8), intent(in) :: val1, val2, err1, err2

      real(8) :: calc_prod_err, err
      real(8) :: p1, p2


      p1  = val1 * err2
      p2  = val2 * err1
      err = sqrt(p1 * p1 + p2 * p2)

      ! return value
      calc_prod_err = err

    end function calc_prod_err
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    function calc_div_err(val1, val2, err1, err2)
!
!     val1, err1 : numerator
!     val2, err2 : denominator
!-----------------------------------------------------------------------
      implicit none

      real(8), intent(in) :: val1, val2, err1, err2

      real(8) :: calc_div_err, err
      real(8) :: p1, p2


      p1  = err1 / val2
      p2  = val1 * err2 / (val2 * val2) 
      err = sqrt(p1 * p1 + p2 * p2)

      ! return value
      calc_div_err = err

    end function calc_div_err
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
