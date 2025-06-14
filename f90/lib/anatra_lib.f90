!===============================================================================
module mod_anatra_ermod
!===============================================================================
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_com
  use mod_grid3d

  implicit none

  ! constants
  !
  integer,      parameter :: MaxState = 100
  integer,      parameter :: NdimMax  = 5 

  integer,      parameter :: ParmFormatANAPARM   = 1
  integer,      parameter :: ParmFormatPRMTOP    = 2
  character(*), parameter :: ParmFormatTypes(2)  = (/'ANAPARM   ',& 
                                                     'PRMTOP    '/)

  integer,      parameter :: VdwTypeSTANDARD   = 1
  integer,      parameter :: VdwTypeATTRACTIVE = 2
  integer,      parameter :: VdwTypeREPULSIVE  = 3
  integer,      parameter :: VdwTypeMINDIST    = 4 
  character(*), parameter :: VdwType(4)        = (/'STANDARD  ',& 
                                                   'ATTRACTIVE',&
                                                   'REPULSIVE ',&
                                                   'MINDIST   '/)

  integer,      parameter :: ElecTypeBARE       = 1
  integer,      parameter :: ElecTypePME        = 2
  character(*), parameter :: ElecType(2)        = (/'BARE      ',& 
                                                    'PME       '/)

  integer,      parameter :: TotTypeSTANDARD   = 1
  integer,      parameter :: TotTypeATTRACTIVE = 2
  character(*), parameter :: TotTypes(2)        = (/'STANDARD  ',& 
                                                    'ATTRACTIVE'/)


  integer,      parameter :: ReacCoordDISTANCE = 1
  integer,      parameter :: ReacCoordENERGY   = 2
  integer,      parameter :: ReacCoordCOSZ     = 3 
  character(*), parameter :: ReacCoordType(3)  = (/'DISTANCE  ',&
                                                   'ENERGY    ',&
                                                   'COSZ      '/)

  ! structures
  !
  type s_eneopt
    logical :: calc_vdw
    logical :: calc_elec
    logical :: ljcheck
    logical :: separate_comp
    integer :: uid, vid
    integer :: mode(2)
    integer :: tottype
    integer :: vdw
    integer :: elec
    real(8) :: umin(MaxState)
    real(8) :: umax(MaxState)
    real(8) :: umin_lj(MaxState)
    real(8) :: umax_lj(MaxState)
    real(8) :: umin_el(MaxState)
    real(8) :: umax_el(MaxState)
  end type s_eneopt

  type s_disopt
    integer :: uid, vid
    integer :: mode(2)
    real(8) :: rmin(MaxState), rmax(MaxState)
  end type s_disopt

  type s_oriopt
    integer :: uid1, uid2 
    integer :: mode(2)
    real(8) :: omin(MaxState), omax(MaxState)
  end type s_oriopt

  type :: s_option
    integer        :: ndim             = 1
    integer        :: nstate           = 1
    integer        :: parmformat       = ParmFormatPRMTOP
    logical        :: pbc              = .false.
    logical        :: use_sdf          = .false.
    logical        :: project_rpl      = .false.
    real(8)        :: rljcut           = 12.0d0
    real(8)        :: relcut           = 1.0d10
    integer        :: grid_resolution  = 1
    real(8)        :: box(3)           = 0.0d0
    real(8)        :: box_shrink(3)    = 0.0d0
    real(8)        :: sdf_threshold    = -1.0d10
    integer        :: reac_coords(NdimMax)   = 0
    integer        :: uid_fit          = 0
    integer        :: vid_fit          = 0
    integer        :: uid_ins          = 0
    type(s_eneopt) :: eneopt(NdimMax)
    type(s_disopt) :: disopt(NdimMax)
    type(s_oriopt) :: oriopt(NdimMax)

    ! Values in the following variables are changed during the operations 
    !
    logical        :: calc_vdw
    logical        :: calc_elec
    integer        :: tottype
    integer        :: vdw
    integer        :: elec
    integer        :: mode(2)
  end type s_option

  ! subroutines
  !
  public  :: check_ioparams
  public  :: read_ctrl_option
  private :: read_ctrl_energy
  private :: read_ctrl_distance

  contains
!
!-----------------------------------------------------------------------
    subroutine check_ioparams(input, output)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in) :: input
      type(s_output), intent(in) :: output 


      ! Check INPUT
      !
      if (trim(input%fxyz(1)) == "") then
        write(iw,'("Check_IOparams> Error.")')
        write(iw,'("fxyz should be specified in input_param.")')
        stop
      end if

      if (trim(input%fdxs(1)) == "") then
        write(iw,'("Check_IOparams> Error.")')
        write(iw,'("fdxs should be specified in input_param.")')
        stop
      end if

      if (trim(input%fprmtop) == "") then
        write(iw,'("Check_IOparams> Error.")')
        write(iw,'("fprmtop should be specified in input_param.")')
        stop
      end if

      ! Check OUTPUT
      !
      if (trim(output%fhead) == "") then
        write(iw,'("Check_IOparams> Error.")')
        write(iw,'("fhead should be specified in output_param.")')
        stop
      end if


    end subroutine check_ioparams
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine read_ctrl_option(iunit, calc_energy, option, myrank)
!-----------------------------------------------------------------------
      implicit none
!
      integer,           intent(in)  :: iunit
      logical,           intent(in)  :: calc_energy
      type(s_option),    intent(out) :: option
      integer, optional, intent(in)  :: myrank


      ! for option_param namelist
      !
      integer                :: ndim             = 2
      integer                :: nstate           = 1 
      character(len=MaxChar) :: parmformat       = "PRMTOP" 
      character(len=MaxChar) :: reac_coords(NdimMax)   = (/"ENERGY", "ENERGY", "ENERGY", "ENERGY", "ENERGY"/)
      logical                :: pbc              = .false.
      logical                :: use_sdf          = .false.
      logical                :: project_rpl      = .false.
      real(8)                :: rljcut           = 12.0d0
      real(8)                :: relcut           = 1.0d10
      integer                :: grid_resolution  = 1
      real(8)                :: box(3)           = 0.0d0
      real(8)                :: box_shrink(3)    = 0.0d0
      real(8)                :: sdf_threshold    = -1.0d10
      integer                :: uid_fit          = 0
      integer                :: vid_fit          = 0
      integer                :: uid_ins          = 0

      integer                :: i
      integer                :: iopt, ierr
      integer                :: irank

      namelist /option_param/ ndim,             &
                              nstate,           &
                              reac_coords,      &
                              pbc,              &
                              use_sdf,          &
                              project_rpl,      &
                              rljcut,           &
                              relcut,           &
                              grid_resolution,  &
                              box,              &
                              box_shrink,       &
                              sdf_threshold,    &
                              uid_fit,          &
                              vid_fit,          &
                              uid_ins


      if (present(myrank)) then
        irank = myrank
      else
        irank = 0
      end if

      ! Initialize
      !
      nstate          = 1 
      parmformat      = "PRMTOP"    ! Fixed
      grid_resolution = 1           ! Fixed
      uid_fit         = 0
      vid_fit         = 0
      uid_ins         = 0
      box             = 0.0d0
      box_shrink      = 0.0d0
      sdf_threshold   = -1.0d10

      rewind iunit
      read(iunit, option_param)

      if (sdf_threshold < -1.0d5) then
        sdf_threshold = 1.0d-5
      end if

      if (irank == 0) then
        write(iw,*)
        write(iw,'(">> Option section parameters")')
        write(iw,'("ndim             = ", i0)')      ndim
        write(iw,'("nstate           = ", i0)')      nstate
        write(iw,'("reac_coords      = ", 2(a,2x))') (trim(reac_coords(i)), i = 1, 2)
        write(iw,'("pbc              = ", a)')       get_tof(pbc)
        write(iw,'("use_sdf          = ", a)')       get_tof(use_sdf)
        write(iw,'("project_rpl      = ", a)')       get_tof(project_rpl)
        write(iw,'("rljcut           = ",f20.10)')   rljcut
        write(iw,'("relcut           = ",e20.10)')   relcut
        write(iw,'("grid_resolution  = ",i0)')       grid_resolution
        write(iw,'("box              = ", 3f20.10)') box(1:3)
        write(iw,'("box_shrink       = ", 3f20.10)') box_shrink(1:3)
        write(iw,'("sdf_threshold    = ",  f20.10)') sdf_threshold 
        write(iw,'("uid_fit          = ", i0)')      uid_fit
        write(iw,'("vid_fit          = ", i0)')      vid_fit
        write(iw,'("uid_ins          = ", i0)')      uid_ins
      end if

      do i = 1, ndim 
        iopt = get_opt(reac_coords(i), ReacCoordType, ierr)
        if (ierr /= 0) then
          write(iw,'("Read_Ctrl_Option> Error.")')
          write(iw,'("reac_coords = ",a," is not available.")') trim(reac_coords(i))
          stop
        end if
        option%reac_coords(i) = iopt
      end do

      iopt = get_opt(parmformat, ParmFormatTypes, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("parmformat = ",a," is not available.")') trim(parmformat)
        stop
      end if
      option%parmformat    = iopt

      ! Set parameters
      !
      option%ndim             = ndim
      option%nstate           = nstate
      option%pbc              = pbc
      option%use_sdf          = use_sdf
      option%project_rpl      = project_rpl 
      option%rljcut           = rljcut
      option%relcut           = relcut
      option%grid_resolution  = grid_resolution
      option%box              = box
      option%box_shrink       = box_shrink
      option%sdf_threshold    = sdf_threshold
      option%uid_fit          = uid_fit
      option%vid_fit          = vid_fit
      option%uid_ins          = uid_ins

      ! Read reacttion coordinate settings
      !
      do i = 1, option%ndim
        if (option%reac_coords(i) == ReacCoordDISTANCE) then
          call read_ctrl_distance(iunit, i, option%nstate, option%disopt(i), myrank)
        else if (option%reac_coords(i) == ReacCoordENERGY) then
          call read_ctrl_energy(iunit, i, option%nstate, option%eneopt(i), myrank)
        else if (option%reac_coords(i) == ReacCoordCOSZ) then
          call read_ctrl_orient(iunit, i, option%nstate, option%oriopt(i), myrank)
        end if
      end do

      if (.not. calc_energy) then
        if (option%project_rpl) then
          call read_ctrl_distance(iunit, 0, option%nstate, option%disopt(0), myrank)
        end if
      end if

      if (calc_energy) then
        if (nstate /= 1) then
          write(iw,'("Read_Ctrl_Option> Error.")')
          write(iw,'("nstate should be 1 if calc_energy = .true.")')
          stop
        end if
      end if

      ! Combination check
      !
      !if (option%elec == ElecTypePME) then
      !  write(iw,'("Read_Ctrl_Option> Error.")')
      !  write(iw,'("elec = PME is currently not available.")')
      !  stop
      !end if

    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine read_ctrl_energy(iunit, id, nstate, eneopt, myrank)
!-----------------------------------------------------------------------
      implicit none

      integer,           intent(in)  :: iunit
      integer,           intent(in)  :: id
      integer,           intent(in)  :: nstate
      type(s_eneopt),    intent(out) :: eneopt
      integer, optional, intent(in)  :: myrank

      integer                :: uid           = 0
      integer                :: vid           = 0
      logical                :: calc_vdw      = .true.
      logical                :: calc_elec     = .false.
      logical                :: separate_comp = .false.
      character(len=MaxChar) :: tottype       = "STANDARD"
      character(len=MaxChar) :: mode(2)       = (/"WHOLE", "WHOLE"/)
      character(len=MaxChar) :: vdw           = "STANDARD"
      character(len=MaxChar) :: elec          = "BARE"
      real(8)                :: umin(1:MaxState)    = 0.0d0
      real(8)                :: umax(1:MaxState)    = 0.0d0
      real(8)                :: umin_lj(1:MaxState) = 0.0d0   
      real(8)                :: umax_lj(1:MaxState) = 0.0d0
      real(8)                :: umin_el(1:MaxState) = 0.0d0
      real(8)                :: umax_el(1:MaxState) = 0.0d0

      integer                :: i
      integer                :: iopt, ierr
      integer                :: irank

      namelist /energy_param1/ uid,             &
                               vid,             &
                               calc_vdw,        &
                               calc_elec,       &
                               separate_comp,   &
                               vdw,             &
                               umin,            &
                               umax,            &
                               umin_lj,         &
                               umax_lj,         &
                               umin_el,         &
                               umax_el

      namelist /energy_param2/ uid,             &
                               vid,             &
                               calc_vdw,        &
                               calc_elec,       &
                               separate_comp,   &
                               vdw,             &
                               umin,            &
                               umax,            &
                               umin_lj,         &
                               umax_lj,         &
                               umin_el,         &
                               umax_el

      namelist /energy_param3/ uid,             &
                               vid,             &
                               calc_vdw,        &
                               calc_elec,       &
                               separate_comp,   &
                               vdw,             &
                               umin,            &
                               umax,            &
                               umin_lj,         &
                               umax_lj,         &
                               umin_el,         &
                               umax_el

      namelist /energy_param4/ uid,             &
                               vid,             &
                               calc_vdw,        &
                               calc_elec,       &
                               separate_comp,   &
                               vdw,             &
                               umin,            &
                               umax,            &
                               umin_lj,         &
                               umax_lj,         &
                               umin_el,         &
                               umax_el

      namelist /energy_param5/ uid,             &
                               vid,             &
                               calc_vdw,        &
                               calc_elec,       &
                               separate_comp,   &
                               vdw,             &
                               umin,            &
                               umax,            &
                               umin_lj,         &
                               umax_lj,         &
                               umin_el,         &
                               umax_el


      if (present(myrank)) then
        irank = myrank
      else
        irank = 0
      end if

      ! Initialize
      !
      uid           = 0
      vid           = 0
      calc_vdw      = .true.
      calc_elec     = .false.
      separate_comp = .false.
      mode          = "WHOLE"    ! fixed in this program
      tottype       = "STANDARD" ! fixed in this program
      vdw           = "STANDARD"
      elec          = "BARE"     ! fixed in this program
      umin          = 0.0d0
      umax          = 0.0d0
      umin_lj       = 0.0d0
      umax_el       = 0.0d0

      ! Read namelist
      !
      rewind iunit
      if (id == 1) then
        read(iunit, energy_param1) 
      else if (id == 2) then
        read(iunit, energy_param2)
      else if (id == 3) then
        read(iunit, energy_param3)
      else if (id == 4) then
        read(iunit, energy_param4)
      else if (id == 5) then
        read(iunit, energy_param5)
      end if

      ! Output parameter info
      !
      if (irank == 0) then
        write(iw,*)
        write(iw,'(">> Energy_Param",i0," parameters")') id
        write(iw,'("uid           = ", i0)')       uid
        write(iw,'("vid           = ", i0)')       vid
        write(iw,'("calc_vdw      = ", a)')        get_tof(calc_vdw)
        write(iw,'("calc_elec     = ", a)')        get_tof(calc_elec)
        write(iw,'("separate_comp = ", a)')        get_tof(separate_comp)
        write(iw,'("vdw           = ", a)')        trim(vdw)

        if (nstate == 1) then
          write(iw,'("umin          = ", f15.7)')    umin(1)
          write(iw,'("umax          = ", f15.7)')    umax(1)

          if (separate_comp) then
            write(iw,'("umin_lj       = ", f15.7)')    umin_lj(1)
            write(iw,'("umax_lj       = ", f15.7)')    umax_lj(1)
            write(iw,'("umin_el       = ", f15.7)')    umin_el(1)
            write(iw,'("umax_el       = ", f15.7)')    umax_el(1)
          end if
        else
          write(iw,'("+--- Energy range ---+")')
          do i = 1, nstate
            write(iw,'("State ", i0, " : ", 2f20.10)') i, umin(i), umax(i) 
          end do
          write(iw,*)
        end if

      end if

      ! Convert mode
      !
      do i = 1, 2
        iopt = get_opt(mode(i), CoMMode, ierr)
        if (ierr /= 0) then
          write(iw,'("Read_Ctrl_Option> Error.")')
          write(iw,'("mode = ",a," is not available.")') trim(mode(i))
          stop
        end if
        eneopt%mode(i)     = iopt
      end do

      ! Convert tottype
      !
      iopt = get_opt(tottype, TotTypes, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("tottype = ",a," is not available.")') trim(tottype)
        stop
      end if
      eneopt%tottype       = iopt

      ! Convert vdw 
      !
      iopt = get_opt(vdw, VdwType, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("vdw = ",a," is not available.")') trim(vdw)
        stop
      end if
      eneopt%vdw           = iopt

      ! Convert elec
      !
      iopt = get_opt(elec, ElecType, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("elec = ",a," is not available.")') trim(elec)
        stop
      end if
      eneopt%elec          = iopt

      ! Error check
      !
      if (uid == 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("uid should be not zero.")')
        stop
      end if

      if (vid == 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("vid should be not zero.")')
        stop
      end if

      if (eneopt%separate_comp) then
        if ((.not. calc_vdw) .or. (.not. calc_elec) .or. &
            eneopt%vdw  /= VdwTypeSTANDARD) then
          write(iw,'("Read_Ctrl_Option> Error.")')
          write(iw,'("separate_comp = .true. is available only if Standard VdW and Elec.")')
          stop
        end if
      end if

      eneopt%uid           = uid
      eneopt%vid           = vid
      eneopt%calc_vdw      = calc_vdw
      eneopt%calc_elec     = calc_elec
      eneopt%separate_comp = separate_comp 
      eneopt%umin          = umin
      eneopt%umax          = umax
      eneopt%umin_lj       = umin_lj
      eneopt%umax_lj       = umax_lj
      eneopt%umin_el       = umin_el
      eneopt%umax_el       = umax_el

    end subroutine read_ctrl_energy
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine read_ctrl_distance(iunit, id, nstate, disopt, myrank)
!-----------------------------------------------------------------------
      implicit none

      integer,           intent(in)  :: iunit
      integer,           intent(in)  :: id
      integer,           intent(in)  :: nstate 
      type(s_disopt),    intent(out) :: disopt
      integer, optional, intent(in)  :: myrank

      integer                :: uid
      integer                :: vid
      real(8)                :: rmin(1:MaxState)
      real(8)                :: rmax(1:MaxState)
      character(len=MaxChar) :: mode(2)

      integer  :: i
      integer  :: iopt, ierr
      integer  :: irank

      namelist /distance_param0/ uid,  &
                                 vid

      namelist /distance_param1/ uid,  &
                                 vid,  &
                                 rmin, &
                                 rmax

      namelist /distance_param2/ uid,  &
                                 vid,  &
                                 rmin, &
                                 rmax

      namelist /distance_param3/ uid,  &
                                 vid,  &
                                 rmin, &
                                 rmax

      namelist /distance_param4/ uid,  &
                                 vid,  &
                                 rmin, &
                                 rmax

      namelist /distance_param5/ uid,  &
                                 vid,  &
                                 rmin, &
                                 rmax


      if (present(myrank)) then
        irank = myrank
      else
        irank = 0
      end if

      ! Initialize
      !
      uid  = 0
      vid  = 0
      mode = "WHOLE"  ! fixed in this program
      rmin = 0.0d0
      rmax = 0.0d0

      ! Read namelist
      !
      rewind iunit
      if (id == 0) then
        read(iunit, distance_param0) 
      else if (id == 1) then
        read(iunit, distance_param1) 
      else if (id == 2) then
        read(iunit, distance_param2)
      else if (id == 3) then
        read(iunit, distance_param3)
      else if (id == 4) then
        read(iunit, distance_param4)
      else if (id == 5) then
        read(iunit, distance_param5)
      end if

      ! Output parameter info
      !
      if (irank == 0) then
        write(iw,*)
        write(iw,'(">> Distance_Param",i0," parameters")') id
        write(iw,'("uid  = ", i0)') uid
        write(iw,'("vid  = ", i0)') vid
        if (id /= 0) then
          if (nstate == 1) then
            write(iw,'("rmin          = ", f15.7)')    rmin(1)
            write(iw,'("rmax          = ", f15.7)')    rmax(1)
          else
            write(iw,'("+--- Distance range ---+")')
            do i = 1, nstate
              write(iw,'("State ", i0, " : ", 2f20.10)') rmin(i), rmax(i) 
            end do
            write(iw,*)
          end if
        end if
      end if

      ! Convert mode
      !
      do i = 1, 2
        iopt = get_opt(mode(i), CoMMode, ierr)
        if (ierr /= 0) then
          write(iw,'("Read_Ctrl_Option> Error.")')
          write(iw,'("mode = ",a," is not available.")') trim(mode(i))
          stop
        end if
        disopt%mode(i)     = iopt
      end do

      ! Error check
      !
      if (uid == 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("uid should be not zero.")')
        stop
      end if

      if (vid == 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("vid should be not zero.")')
        stop
      end if

      disopt%uid  = uid
      disopt%vid  = vid
      disopt%rmin = rmin
      disopt%rmax = rmax

    end subroutine read_ctrl_distance
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine read_ctrl_orient(iunit, id, nstate, oriopt, myrank)
!-----------------------------------------------------------------------
      implicit none

      integer,           intent(in)  :: iunit
      integer,           intent(in)  :: id
      integer,           intent(in)  :: nstate 
      type(s_oriopt),    intent(out) :: oriopt
      integer, optional, intent(in)  :: myrank

      integer                :: uid1
      integer                :: uid2
      real(8)                :: omin(1:MaxState)
      real(8)                :: omax(1:MaxState)
      character(len=MaxChar) :: mode(2)

      integer  :: i
      integer  :: iopt, ierr
      integer  :: irank

      namelist /orient_param1/ uid1, &
                               uid2, &
                               omin, &
                               omax

      namelist /orient_param2/ uid1, &
                               uid2, &
                               omin, &
                               omax

      namelist /orient_param3/ uid1, &
                               uid2, &
                               omin, &
                               omax

      namelist /orient_param4/ uid1, &
                               uid2, &
                               omin, &
                               omax

      namelist /orient_param5/ uid1, &
                               uid2, &
                               omin, &
                               omax


      if (present(myrank)) then
        irank = myrank
      else
        irank = 0
      end if

      ! Initialize
      !
      uid1 = 0
      uid2 = 0
      mode = "WHOLE"  ! fixed in this program
      omin = 0.0d0
      omax = 0.0d0

      ! Read namelist
      !
      rewind iunit
      if (id == 1) then
        read(iunit, orient_param1) 
      else if (id == 2) then
        read(iunit, orient_param2)
      else if (id == 3) then
        read(iunit, orient_param3)
      else if (id == 4) then
        read(iunit, orient_param4)
      else if (id == 5) then
        read(iunit, orient_param5)
      end if

      ! Output parameter info
      !
      if (irank == 0) then
        write(iw,*)
        write(iw,'(">> Distance_Param",i0," parameters")') id
        write(iw,'("uid1 = ", i0)') uid1
        write(iw,'("uid2 = ", i0)') uid2
        if (id /= 0) then
          if (nstate == 1) then
            write(iw,'("omin          = ", f15.7)')    omin(1)
            write(iw,'("omax          = ", f15.7)')    omax(1)
          else
            write(iw,'("+--- Orientation range ---+")')
            do i = 1, nstate
              write(iw,'("State ", i0, " : ", 2f20.10)') omin(i), omax(i) 
            end do
            write(iw,*)
          end if
        end if
      end if

      ! Convert mode
      !
      do i = 1, 2
        iopt = get_opt(mode(i), CoMMode, ierr)
        if (ierr /= 0) then
          write(iw,'("Read_Ctrl_Option> Error.")')
          write(iw,'("mode = ",a," is not available.")') trim(mode(i))
          stop
        end if
        oriopt%mode(i)     = iopt
      end do

      ! Error check
      !
      if (uid1 == 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("uid1 should be not zero.")')
        stop
      end if

      if (uid2 == 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("uid2 should be not zero.")')
        stop
      end if

      oriopt%uid1 = uid1
      oriopt%uid2 = uid2 
      oriopt%omin = omin
      oriopt%omax = omax

    end subroutine read_ctrl_orient
!-----------------------------------------------------------------------

end module mod_anatra_ermod
!===============================================================================

