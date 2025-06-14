!=======================================================================
!
!  Module   mod_const
!
!    Module for defining parameters commonly used in many modules of ANATRA  
!
!    (c) Copyright 2024 Osaka Univ. All rights reserved.
!
!=======================================================================

module mod_const

  use mod_util
  implicit none

  ! parameters
  !
  integer, parameter :: MaxChar     = 1000
  integer, parameter :: MaxTraj     = 5000
  integer, parameter :: ir          = 5, iw = 6
  integer, parameter :: UnitIN      = 10
  integer, parameter :: UnitOUT     = 20
  integer, parameter :: UnitOUT2    = 21
  integer, parameter :: UnitOUT3    = 22
  integer, parameter :: UnitOUT4    = 23
  integer, parameter :: UnitOUT5    = 24
  integer, parameter :: UnitOUT6    = 25
  integer, parameter :: UnitOUT7    = 26
  integer, parameter :: UnitDCD     = 50
  integer, parameter :: UnitCV      = 51
  integer, parameter :: UnitMPL     = 52
  integer, parameter :: UnitGNPLT   = 53 
  integer, parameter :: UnitDX      = 54
  integer, parameter :: UnitPRMTOP  = 55
  integer, parameter :: UnitANAPARM = 56
  integer, parameter :: UnitXTC     = 57

  ! physical constans 
  !
  real(8), parameter :: Pi         = acos(-1.0d0)
  real(8), parameter :: Avogadro   = 6.02214076d+23
  real(8), parameter :: Boltz      = 1.987204292510d-03  ! kcal mol^-1
  real(8), parameter :: ElecCharge = 1.602176634d-19     ! C

  ! constants for unit conversion
  !
  ! --- energy 
  real(8), parameter :: Ene_kJmol_to_kcalmol   = 2.390057360d-01
  real(8), parameter :: Ene_kJmol_to_hartree   = 3.808798850d-04  
  real(8), parameter :: Ene_kcalmol_to_kjmol   = 4.184000000d-00
  real(8), parameter :: Ene_kcalmol_to_hartree = 1.593601440d-03
  real(8), parameter :: Ene_kcalmol_to_J       = 6.947695457d-21
  real(8), parameter :: Ene_kcalmol_to_eV      = 0.043370000d0
  real(8), parameter :: Ene_hartree_to_kJmol   = 2.625499640d+03
  real(8), parameter :: Ene_hartree_to_kcalmol = 6.275094740d+02
  ! 
  ! --- length
  real(8), parameter :: Len_angs_to_nm         = 1.000000000d-01 
  real(8), parameter :: Len_angs_to_cm         = 1.000000000d-08
  real(8), parameter :: Len_angs_to_dm         = 1.000000000d-09
  real(8), parameter :: Len_angs_to_m          = 1.000000000d-10 
  real(8), parameter :: Len_angs_to_bohr       = 1.889726130d+00
  real(8), parameter :: Len_nm_to_angs         = 1.000000000d+01
  real(8), parameter :: Len_cm_to_angs         = 1.000000000d+08 
  real(8), parameter :: Len_dm_to_angs         = 1.000000000d+09 
  real(8), parameter :: Len_m_to_angs          = 1.000000000d+10
  real(8), parameter :: Len_bohr_to_angs       = 5.291772110d-01
  ! 
  ! --- volume
  real(8), parameter :: Vol_angs_to_nm         = 1.000000000d-03
  real(8), parameter :: Vol_angs_to_cm         = 1.000000000d-24
  real(8), parameter :: Vol_angs_to_dm         = 1.000000000d-27
  real(8), parameter :: Vol_angs_to_m          = 1.000000000d-30
  real(8), parameter :: Vol_angs_to_bohr       = 6.748334552d+00
  real(8), parameter :: Vol_nm_to_angs         = 1.000000000d+03
  real(8), parameter :: Vol_cm_to_angs         = 1.000000000d+24
  real(8), parameter :: Vol_dm_to_angs         = 1.000000000d+27
  real(8), parameter :: Vol_m_to_angs          = 1.000000000d+30
  real(8), parameter :: Vol_bohr_to_angs       = 1.481847116E-01

  ! subroutines
  ! 

end module mod_const
!=======================================================================
