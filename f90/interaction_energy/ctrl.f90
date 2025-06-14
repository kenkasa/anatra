!=======================================================================
module mod_ctrl
!=======================================================================
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_traj
  use mod_com
  implicit none

  ! constants
  !
  integer,      parameter :: ParmFormatANAPARM   = 1
  integer,      parameter :: ParmFormatPRMTOP    = 2
  character(*), parameter :: ParmFormatTypes(2)  = (/'ANAPARM   ',& 
                                                     'PRMTOP    '/)

  integer,      parameter :: VdwTypeSTANDARD   = 1
  integer,      parameter :: VdwTypeATTRACTIVE = 2
  integer,      parameter :: VdwTypeREPULSIVE  = 3
  integer,      parameter :: VdwTypeNUMBER     = 4 
  integer,      parameter :: VdwTypeUATINTRA   = 5 
  character(*), parameter :: VdwType(5)        = (/'STANDARD  ',& 
                                                   'ATTRACTIVE',&
                                                   'REPULSIVE ',&
                                                   'NUMBER    ',&
                                                   'UAT-INTRA '/)

  integer,      parameter :: ElecTypeBARE       = 1
  integer,      parameter :: ElecTypePME        = 2
  character(*), parameter :: ElecType(2)        = (/'BARE      ',& 
                                                    'PME       '/)

  ! structures
  !
  type :: s_option
    integer :: parmformat       = ParmFormatPRMTOP
    logical :: pbc              = .false.
    logical :: calc_vdw         = .true. 
    logical :: calc_elec        = .false.
    logical :: calc_torque      = .false.
    logical :: calc_siteesp     = .false.
    integer :: mode(2)          = (/ComModeRESIDUE, ComModeRESIDUE/)
    integer :: vdw_type         = VdwTypeSTANDARD
    integer :: elec_type        = ElecTypeBARE
    real(8) :: rvdwcut           = 12.0d0
    real(8) :: relcut           = 1.0d10
    real(8) :: pme_alpha        = 0.35d0
    integer :: pme_grids(3)     = 0
    integer :: pme_spline_order = 4
    logical :: pme_rigid        = .false.
    logical :: pme_dual         = .false.

    ! Hidden options
    !
    logical :: output_param     = .false.

  end type s_option

  ! subroutines
  !
  public  :: read_ctrl
  private :: read_ctrl_option

  contains

!-----------------------------------------------------------------------
    subroutine read_ctrl(input, output, option, trajopt)
!-----------------------------------------------------------------------
      implicit none

      integer, parameter           :: iunit = 10

      type(s_input),   intent(out) :: input
      type(s_output),  intent(out) :: output
      type(s_option),  intent(out) :: option
      type(s_trajopt), intent(out) :: trajopt

      character(len=MaxChar)       :: f_ctrl

      ! get control file name
      !
      call getarg(1, f_ctrl)

      write(iw,*)
      write(iw,'("Read_Ctrl> Reading parameters from ", a)') trim(f_ctrl)
      open(iunit, file=trim(f_ctrl), status='old')
        call read_ctrl_input  (iunit, input)
        call read_ctrl_output (iunit, output)
        call read_ctrl_option (iunit, option)
        call read_ctrl_trajopt(iunit, trajopt)
      close(iunit)

    end subroutine read_ctrl
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine read_ctrl_option(iunit, option)
!-----------------------------------------------------------------------
      implicit none
!
      integer,        intent(in)  :: iunit
      type(s_option), intent(out) :: option 

      character(len=MaxChar) :: parmformat       = "PRMTOP" 
      logical                :: pbc              = .false.
      logical                :: calc_vdw         = .true.
      logical                :: calc_elec        = .false.
      logical                :: calc_torque      = .false.
      logical                :: calc_siteesp     = .false.
      character(len=MaxChar) :: mode(2)          = (/"RESIDUE", "RESIDUE"/)
      character(len=MaxChar) :: vdw_type         = "STANDARD"
      character(len=MaxChar) :: elec_type        = "BARE"
      real(8)                :: rvdwcut           = 12.0d0
      real(8)                :: relcut           = 1.0d10
      real(8)                :: pme_alpha        = 0.35d0
      integer                :: pme_grids(3)     = (/64, 64, 64/)
      integer                :: pme_spline_order = 4
      logical                :: pme_rigid        = .false.
      logical                :: pme_dual         = .false.

      ! Hidden options
      logical                :: output_param     = .false.

      integer                :: i
      integer                :: iopt, ierr

      namelist /option_param/ parmformat,       &
                              pbc,              &
                              calc_vdw,         &
                              calc_elec,        &
                              calc_torque,      &
                              calc_siteesp,     &
                              mode,             &
                              vdw_type,         &
                              elec_type,        &
                              rvdwcut,           &
                              relcut,           &
                              pme_alpha,        &
                              pme_grids,        &
                              pme_spline_order, &
                              pme_rigid,        &
                              pme_dual,         &
                              output_param


      rewind iunit
      read(iunit, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("parmformat       = ", a)')       trim(parmformat) 
      write(iw,'("pbc              = ", a)')       get_tof(pbc)
      write(iw,'("calc_vdw         = ", a)')       get_tof(calc_vdw)
      write(iw,'("calc_elec        = ", a)')       get_tof(calc_elec)
      write(iw,'("calc_torque      = ", a)')       get_tof(calc_torque)
      write(iw,'("calc_siteesp     = ", a)')       get_tof(calc_siteesp)
      write(iw,'("mode             = ",2(a,2x))')  (trim(mode(i)), i = 1, 2)
      write(iw,'("vdw_type         = ", a)')       trim(vdw_type)
      write(iw,'("elec_type        = ", a)')       trim(elec_type)
      write(iw,'("rvdwcut           = ",f20.10)')   rvdwcut
      write(iw,'("relcut           = ",e20.10)')   relcut
      write(iw,'("pme_alpha        = ",f20.10)')   pme_alpha
      write(iw,'("pme_grids        = ",3(i0,2x))') (pme_grids(i), i = 1, 3)
      write(iw,'("pme_spline_order = ",i0)')       pme_spline_order 
      write(iw,'("pme_rigid        = ", a)')       get_tof(pme_rigid)
      write(iw,'("pme_dual         = ", a)')       get_tof(pme_dual)

      iopt = get_opt(parmformat, ParmFormatTypes, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("parmformat = ",a," is not available.")') trim(parmformat)
        stop
      end if
      option%parmformat    = iopt

      do i = 1, 2
        iopt = get_opt(mode(i), CoMMode, ierr)
        if (ierr /= 0) then
          write(iw,'("Read_Ctrl_Option> Error.")')
          write(iw,'("mode = ",a," is not available.")') trim(mode(i))
          stop
        end if
        option%mode(i)     = iopt
      end do

      iopt = get_opt(vdw_type, VdwType, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("vdw_type = ",a," is not available.")') trim(vdw_type)
        stop
      end if
      option%vdw_type     = iopt

      iopt = get_opt(elec_type, ElecType, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("elec_type = ",a," is not available.")') trim(elec_type)
        stop
      end if
      option%elec_type        = iopt
                         
      option%pbc              = pbc
      option%calc_vdw         = calc_vdw
      option%calc_elec        = calc_elec
      option%calc_torque      = calc_torque
      option%calc_siteesp     = calc_siteesp
      option%rvdwcut          = rvdwcut
      option%relcut           = relcut
      option%pme_alpha        = pme_alpha
      option%pme_grids        = pme_grids
      option%pme_spline_order = pme_spline_order
      option%pme_rigid        = pme_rigid
      option%pme_dual         = pme_dual
      option%output_param     = output_param

      if (option%calc_siteesp) then
        write(iw,'("Read_Ctrl_Option> REMARK")')
        write(iw,'("detect calc_siteesp = .true.")')
        write(iw,'("=> following parameters are modified")')
        write(iw,'("  calc_vdw    <= false")')
        write(iw,'("  calc_elec   <= true")')
        write(iw,'("  calc_torque <= false")')
        write(iw,'("  pme_dual    <= false")')
        write(iw,'("  mode(1)     <= ATOM")')
        write(iw,'("  mode(2)     <= WHOLE")')
        option%calc_vdw    = .false.
        option%calc_elec   = .true.
        option%calc_torque = .false.
        option%pme_dual    = .false.
        option%mode(1)     = ComModeATOM
        option%mode(2)     = ComModeWHOLE 
      end if

      ! Combination check
      !
      if (.not. option%pbc .and. option%elec_type == ElecTypePME) then
        write(iw,'("Read_Ctrl_Option> Error. pbc = .true. and elec = PME")')
        write(iw,'("                  can not be used at the same time.")')
        stop
      end if

    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end module
!=======================================================================
