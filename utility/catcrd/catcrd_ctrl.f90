!===============================================================================
module mod_catcrd_ctrl
!===============================================================================
  use mod_const
  use mod_util
  use mod_dcdio
  use mod_xtcio
  use mod_netcdfio
  use mod_parse_arguments

  implicit none

  ! constants
  !
  integer, parameter, public :: OptionIndexI       = 1
  integer, parameter, public :: OptionIndexO       = 2 
  integer, parameter, public :: OptionIndexStride  = 3
  integer, parameter, public :: OptionIndexFirst   = 4
  integer, parameter, public :: OptionIndexLast    = 5
  integer, parameter, public :: OptionIndexRect    = 6
  integer, parameter, public :: OptionIndexSelfile = 7

  integer, parameter, public :: OptionIndices(7) =     &
                                  (/OptionIndexI,      &
                                    OptionIndexO,      &
                                    OptionIndexStride, &
                                    OptionIndexFirst,  &
                                    OptionIndexLast,   &
                                    OptionIndexRect,   &
                                    OptionIndexSelfile &
                                    /)

  character(*), parameter, public :: OptionNames(7) =     &
                                       (/'-I      ',      &
                                         '-O      ',      &
                                         '-STRIDE ',      &
                                         '-FIRST  ',      &
                                         '-LAST   ',      &
                                         '-RECT   ',      &
                                         '-SELFILE'       &
                                        /)

  integer, parameter, public :: TrjTypeDCD    = 1
  integer, parameter, public :: TrjTypeXTC    = 2
  integer, parameter, public :: TrjTypeNCD    = 3 

  ! structures
  !
  type :: s_option
    character(len=MaxChar), allocatable :: trj_in(:)
    character(len=MaxChar)              :: trj_out
    integer                             :: stride
    integer                             :: first
    integer                             :: last
    logical                             :: rect
    character(len=MaxChar)              :: selfile
    integer                             :: ntrj_in
    logical                             :: out_exist
    integer                             :: trjtype_in
    integer                             :: trjtype_out
    logical                             :: selfile_exist

    integer                             :: natm
  end type s_option 

  ! subroutines
  !
  public :: setup_catcrd_options
  
  contains
!-------------------------------------------------------------------------------
    subroutine setup_catcrd_options(parg, option)
!-------------------------------------------------------------------------------
      implicit none

      type(s_parg),    intent(in)  :: parg
      type(s_option),  intent(out) :: option 


      integer :: nopts
      integer :: i, j, ind, ipos
      integer :: ilen
      integer :: dotloc, trjtype


      ! Setup Traj_in
      !
      ind     = parg%opt_ind(OptionIndexI)
      ipos    = parg%opt_pos(ind)
      ilen    = 0
      if (ipos /= 0) &
        ilen = parg%opt_len(ind)

      if (ipos == 0) then
        write(iw,'("Setup_Catcrd_Options> Error.")')
        write(iw,'("-i option is neccesary")')
        stop
      end if

      if (ilen >= 1) then
        allocate(option%trj_in(ilen))
        option%trj_in = ""
        do i = 1, ilen
          option%trj_in(i) = trim(parg%args(ipos + i))
        end do
        option%ntrj_in = ilen
      else
        write(iw,'("Setup_Catcrd_Options> Error.")')
        write(iw,'("input trajectories should be specified after -i option.")')
        stop
      end if

      ! Setup Trj_out 
      !
      ind     = parg%opt_ind(OptionIndexO)
      ipos    = parg%opt_pos(ind)
      ilen    = 0
      if (ipos /= 0) &
        ilen = parg%opt_len(ind)

      option%trj_out   = ""
      option%out_exist = .true.
      if (ipos == 0) then
        option%out_exist = .false. 
      else
        if (ilen >= 1) then
          option%trj_out = trim(parg%args(ipos + 1))
        else
          write(iw,'("Setup_Catcrd_Options> Error.")')
          write(iw,'("Output trajectory name should be specified after -o option.")')
          stop
        end if
      end if

      ! Setup Stride
      !
      ind     = parg%opt_ind(OptionIndexStride)
      ipos    = parg%opt_pos(ind)
      ilen    = 0
      if (ipos /= 0) &
        ilen = parg%opt_len(ind)

      option%stride = 1
      if (ipos == 0) then
        option%stride = 1
      else
        read(parg%args(ipos + 1),*) option%stride
      end if

      ! Setup First
      !
      ind     = parg%opt_ind(OptionIndexFirst)
      ipos    = parg%opt_pos(ind)
      ilen    = 0
      if (ipos /= 0) &
        ilen = parg%opt_len(ind)

      option%first = 1
      if (ipos == 0) then
        option%first = 1 
      else
        if (ilen >= 1) then
          read(parg%args(ipos + 1), *) option%first
        else
          write(iw,'("Setup_Catcrd_Options> Error.")')
          write(iw,'("First parameter should be specified after -first option.")')
          stop
        end if
      end if

      ! Setup Last 
      !
      ind     = parg%opt_ind(OptionIndexLast)
      ipos    = parg%opt_pos(ind)
      ilen    = 0
      if (ipos /= 0) &
        ilen = parg%opt_len(ind)

      option%last = 0 
      if (ipos == 0) then
        option%last = 0
      else
        if (ilen >= 1) then
          read(parg%args(ipos + 1), *) option%last
        else
          write(iw,'("Setup_Catcrd_Options> Error.")')
          write(iw,'("Last parameter should be specified after -last option.")')
          stop
        end if
      end if

      ! Setup RECT 
      !
      ind     = parg%opt_ind(OptionIndexRect)
      ipos    = parg%opt_pos(ind)

      option%rect = .true. 
      if (ipos == 0) then
        option%rect = .false.
      end if

      ! Setup SELFILE 
      !
      ind     = parg%opt_ind(OptionIndexSelfile)
      ipos    = parg%opt_pos(ind)
      ilen    = 0
      if (ipos /= 0) &
        ilen = parg%opt_len(ind)

      option%selfile       = ""
      option%selfile_exist = .true.
      if (ipos == 0) then
        option%selfile_exist = .false. 
      else
        if (ilen >= 1) then
          option%selfile = trim(parg%args(ipos + 1))
        else
          write(iw,'("Setup_Catcrd_Options> Error.")')
          write(iw,'("Selfile name should be specified after -selfile option.")')
          stop
        end if
      end if


      ! get input/output trajectory types
      !

      ! input
      dotloc = index(option%trj_in(1), ".", back = .true.)
      if (dotloc == 0) then
        write(iw,'("Setup_Catcrd_Options> Error.")')
        write(iw,'("extention (.XXX) is not found in input trajectory file name")')
        stop
      end if

      if (option%trj_in(1)(dotloc+1:dotloc+3) == "dcd") then
        option%trjtype_in = TrjTypeDCD 
      else if (option%trj_in(1)(dotloc+1:dotloc+3) == "xtc") then
        option%trjtype_in = TrjTypeXTC
      else if (option%trj_in(1)(dotloc+1:dotloc+2) == "nc") then
        option%trjtype_in = TrjTypeNCD
      else
        write(iw,'("Setup_Catcrd_Options> Error.")')
        write(iw,'("Unknown trjactory file format")')
        stop
      end if

      ! output
      if (option%out_exist) then
        dotloc = index(option%trj_out, ".", back = .true.)
        if (dotloc == 0) then
          write(iw,'("Setup_Catcrd_Options> Error.")')
          write(iw,'("extention (.XXX) is not found in output trajectory file name")')
          stop
        end if

        if (option%trj_out(dotloc+1:dotloc+3) == "dcd") then
          option%trjtype_out = TrjTypeDCD 
        else if (option%trj_out(dotloc+1:dotloc+3) == "xtc") then
          option%trjtype_out = TrjTypeXTC
        else if (option%trj_out(dotloc+1:dotloc+2) == "nc") then
          option%trjtype_out = TrjTypeNCD
        end if
      end if

      ! Print out option parameters
      !
      write(iw,'("Setup_Catcrd_Options> Option parameters")')
      write(iw,*)

      write(iw,'("input trajectories")')
      do i = 1, option%ntrj_in
        write(iw,'(2x,i5,2x,a)') i, trim(option%trj_in(i))
      end do
      write(iw,*)
      write(iw,'("output trajectory = ", a)')  trim(option%trj_out)
      write(iw,'("stride            = ", i0)') option%stride
      write(iw,'("first             = ", i0)') option%first
      write(iw,'("last              = ", i0)') option%last
      write(iw,'("rect              = ", a)')  get_tof(option%rect)
      write(iw,'("selfile           = ", a)')  trim(option%selfile)
      !write(iw,'("trjtype_in        = ", i0)') option%trjtype_in
      !write(iw,'("trjtype_out       = ", i0)') option%trjtype_out


    end subroutine setup_catcrd_options
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
    subroutine get_trajectory_info(option) 
!-------------------------------------------------------------------------------
      implicit none

      type(s_option), intent(inout) :: option


      ! get number of atoms (natm) from input trajectory file
      !
      write(iw,*)
      write(iw,'("Get_Trajectory_info> Getting number of atoms in trajectory...")')
      if (option%trjtype_in == TrjTypeDCD) then
        call get_natm_from_dcd(option%trj_in(1), option%natm)
        write(iw,'("number of atoms = ",i0)') option%natm
      else if (option%trjtype_in == TrjTypeXTC) then
        call get_natm_from_xtc(option%trj_in(1), option%natm)
      else if (option%trjtype_in == TrjTypeNCD) then
        call get_natm_from_netcdf(option%trj_in(1), option%natm)
      end if


    end subroutine get_trajectory_info
!-------------------------------------------------------------------------------

end module mod_catcrd_ctrl
!===============================================================================
