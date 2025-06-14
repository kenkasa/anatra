!=======================================================================
module mod_analyze_str
!=======================================================================
  use mod_util
  use mod_const

  ! constants
  !

  ! structures
  !
  type :: s_state
    integer :: nstep
    integer :: step_final
    integer :: nmol
    integer :: ndim
    integer :: nBevents    = 0
    integer :: nUBevents   = 0
    !integer :: quench_step = 0
    logical, allocatable :: is_reacted(:)
    logical :: is_reacted_final = .false.
    integer, allocatable :: quench_step(:)
    integer, allocatable :: read_step(:)
    !integer, allocatable :: data(:)
    integer, allocatable :: data(:, :)
    !real(8), allocatable :: data_ts(:) 
    real(8), allocatable :: data_ts(:, :)
    !real(8), allocatable :: avedata(:)
    real(8), allocatable :: avedata(:, :)
    real(8), allocatable :: hist(:, :, :)
    real(8), allocatable :: norm(:, :)
    real(8), allocatable :: hist2(:, :)
  end type s_state

  type :: s_booteach
    integer              :: ntrial
    integer              :: nsample
    integer              :: nt_range
    integer, allocatable :: rand(:, :)
    real(8), allocatable :: taud(:), kd(:)
    real(8), allocatable :: tau_tpm(:, :, :), k_tpm(:, :, :)
    real(8), allocatable :: tauins(:), kins(:)
    real(8), allocatable :: hist(:, :, :, :)
    real(8), allocatable :: hist2(:, :, :), pa3(:, :)
    real(8), allocatable :: taud2(:, :), taupa3(:)
  end type s_booteach

  type :: s_bootave
    integer              :: nt_range
    real(8)              :: taud, taud_stdev, taud_err
    real(8)              :: kd, kd_stdev, kd_err
    real(8)              :: tauins, tauins_stdev, tauins_err
    real(8)              :: kins, kins_stdev, kins_err
    real(8)              :: taupa3, taupa3_err
    real(8), allocatable :: hist(:, :, :), hist_err(:, :, :), hist_stdev(:, :, :)
    real(8), allocatable :: cumm(:, :, :)
    real(8), allocatable :: tau_tpm(:, :), tau_tpm_err(:, :), k_tpm(:, :), k_tpm_err(:, :)
    real(8), allocatable :: hist2(:, :), pa3(:)
    real(8), allocatable :: taud2(:), taud2_err(:)
    real(8), allocatable :: cumm2(:, :), cummpa3(:)
  end type s_bootave

  ! subroutines
  !

end module mod_analyze_str
!=======================================================================
