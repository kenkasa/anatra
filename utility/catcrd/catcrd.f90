!===============================================================================
program main
!===============================================================================
  use mod_const
  use mod_util
  use mod_parse_arguments
  use mod_catcrd_ctrl
  use mod_catcrd_analyze

  implicit none

  type(s_parg)   :: parg
  type(s_option) :: option

  integer      :: i


  ! parse arguments
  !
  call get_arguments(parg)
  call assign_option_names(OptionIndices, OptionNames, parg)
  call setup_catcrd_options(parg, option)

  ! get trajectory information 
  !
  call get_trajectory_info(option)

  ! analyze
  !
  call analyze(option)


end program main
!===============================================================================
