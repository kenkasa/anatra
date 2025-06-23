!===============================================================================
module mod_catcrd_analyze
!===============================================================================
  use mod_const
  use mod_util
  use mod_parse_arguments
  use mod_catcrd_ctrl
  use mod_netcdfio
  use mod_dcdio
  use mod_xtcio
  use mod_traj
  use xdr, only: xtcfile 

  implicit none

  ! subroutines
  !
  private :: dcd2xtc
  private :: dcd2dcd
  private :: dcd2netcdf
  private :: xtc2xtc

  contains
!-------------------------------------------------------------------------------
    subroutine analyze(option)
!-------------------------------------------------------------------------------
      implicit none

      type(s_option) :: option


      if (option%trjtype_in == TrjTypeDCD) then

        if (option%trjtype_out == TrjTypeDCD) then
          call dcd2dcd(option)
        else if (option%trjtype_out == TrjTypeXTC) then
          call dcd2xtc(option)
        else if (option%trjtype_out == TrjTypeNCD) then
          call dcd2netcdf(option)
        else
          call dcd2dcd(option) ! print out total # of steps
        end if

      else if (option%trjtype_in == TrjTypeXTC) then

        if (option%trjtype_out == TrjTypeXTC) then
          call xtc2xtc(option)
        else
          call xtc2xtc(option) ! print out total # of steps
        end if

      else if (option%trjtype_in == TrjTypeNCD) then

        if (option%trjtype_out == TrjTypeNCD) then
          call netcdf2netcdf(option)
        else
          write(iw,'("Analyze> Error.")')
          write(iw,'("Convert NetCDF to other formats is not suported yet.")')
          stop
        end if

      end if


    end subroutine analyze
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
    subroutine dcd2xtc(option)
!-------------------------------------------------------------------------------
      implicit none

      type(s_option) :: option

      real(8), parameter :: ang2nm = 0.1d0


      integer :: iatm, id, itrj, istep, istep_tot, nwrite 
      integer :: iunit_in, iunit_out
      integer :: natm, nstep, stride, first, last
      real(8) :: box(3, 3)


      real, allocatable :: coord(:, :)

      type(s_traj)    :: traj
      type(s_trajopt) :: trajopt 
      type(s_dcd)     :: dcd
      type(xtcfile)   :: xtc


      if (option%selfile_exist) then
        trajopt%molinfo(1) = trim(option%selfile)
        trajopt%nmolinfo   = 1
        call setup_traj_from_args(trajopt, 1, traj, trajid = 1)
        natm   = traj%natm
      else
        natm   = option%natm
      end if

      stride = option%stride

      first  = option%first
      last   = option%last

      allocate(coord(1:3, 1:natm))

      call xtc%init(option%trj_out, 'w')

      iunit_in = 10

      istep_tot = 0
      nwrite    = 0
      do itrj = 1, option%ntrj_in
        call dcd_open(option%trj_in(itrj), iunit_in)
        call dcd_read_header(iunit_in, dcd)
        nstep = dcd%nstep

        if (itrj == 1) then
          call alloc_dcd(option%natm, 1, dcd) 
        end if

        do istep = 1, nstep

          istep_tot = istep_tot + 1

          if (mod(istep_tot, 100) == 0) then
            write(iw,'("Read step ", i0)') istep_tot
          end if
          
          call read_dcd_oneframe(iunit_in, dcd)

          if (istep_tot >= first) then
            if (istep_tot <= last .or. last == 0) then

              if (mod(istep_tot, stride) == 0) then
                nwrite = nwrite + 1
                box = 0.0d0
                box(1, 1) = dcd%box(1, 1) * ang2nm 
                box(2, 2) = dcd%box(2, 1) * ang2nm
                box(3, 3) = dcd%box(3, 1) * ang2nm
             
                if (option%selfile_exist) then
                  do iatm = 1, natm
                    id = traj%ind(iatm)
                    coord(1:3, iatm) = dcd%coord(1:3, id, 1) * ang2nm
                  end do
                else
                  coord(1:3, 1:natm) = dcd%coord(1:3, 1:natm, 1) * ang2nm
                end if
                 
                call xtc%write(natm,               &
                               0,                  &
                               0.0,                &
                               real(box),          &
                               coord(1:3, 1:natm), &
                               real(1000.0d0))

              end if

            end if
          end if
        
        end do

        call dcd_close(iunit_in, option%trj_in(itrj))

      end do

      call xtc%close
      deallocate(coord)

      write(iw,'("Total number of frames in output : ", i0)') nwrite

    end subroutine dcd2xtc
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
    subroutine dcd2dcd(option)
!-------------------------------------------------------------------------------
      implicit none

      type(s_option) :: option

      integer :: iatm, id, itrj, istep, istep_tot, nwrite 
      integer :: iunit_in, iunit_out
      integer :: natm, nstep, stride, first, last
      integer :: nstep_tot
      real(8) :: box(3, 3)


      real, allocatable :: coord(:, :)

      type(s_traj)    :: traj
      type(s_trajopt) :: trajopt
      type(s_dcd)     :: dcd_in
      type(s_dcd)     :: dcd_out 


      if (option%selfile_exist) then
        trajopt%molinfo(1) = trim(option%selfile)
        trajopt%nmolinfo   = 1
        call setup_traj_from_args(trajopt, 1, traj, trajid = 1)
        natm   = traj%natm
      else
        natm   = option%natm
      end if

      stride = option%stride

      first  = option%first
      last   = option%last

      iunit_in  = 10
      iunit_out = 11

      ! get total number of steps in input trajectories
      !
      !nstep_tot = 0
      !do itrj = 1, option%ntrj_in
      !  call dcd_open(option%trj_in(itrj), iunit_in)
      !  call dcd_read_header(iunit_in, dcd_in)
      !  nstep_tot = nstep_tot + dcd_in%nstep
      !  call dcd_close(iunit_in, option%trj_in(itrj))
      !end do
      call get_total_step_from_dcd(option%trj_in, nstep_tot)

      write(iw,'("Total number of steps is ",i0)') nstep_tot

      if (.not. option%out_exist) &
        return


      nwrite = 0
      do istep = 1, nstep_tot
        if (istep >= first) then
          if (istep <= last .or. last == 0) then
            nwrite = nwrite + 1 
          end if
        end if
      end do

      call dcd_open(option%trj_out, iunit_out)
      call alloc_dcd(natm, 1, dcd_out)

      istep_tot = 0
      nwrite    = 0

      do itrj = 1, option%ntrj_in
        call dcd_open(option%trj_in(itrj), iunit_in)
        call dcd_read_header(iunit_in, dcd_in)
        nstep = dcd_in%nstep

        if (itrj == 1) then
          call alloc_dcd(option%natm, 1, dcd_in)

          dcd_out%natm       = natm 
          dcd_out%dcdinfo    = dcd_in%dcdinfo
          dcd_out%dcdinfo(1) = nwrite

          call dcd_write_header(iunit_out, dcd_out)
        end if

        do istep = 1, nstep
          istep_tot = istep_tot + 1

          if (mod(istep_tot, 100) == 0) then
            write(iw,'("Read step ", i0)') istep_tot
          end if

          call read_dcd_oneframe(iunit_in, dcd_in)

          if (istep_tot >= first) then
            if (istep_tot <= last .or. last == 0) then

              if (mod(istep_tot, stride) == 0) then
                if (option%selfile_exist) then
                  do iatm = 1, natm
                    id = traj%ind(iatm)
                    dcd_out%coord(1:3, iatm, 1) = dcd_in%coord(1:3, id, 1)
                  end do
                else
                  dcd_out%coord = dcd_in%coord
                end if
               
                dcd_out%box = dcd_in%box
                call write_dcd_oneframe(iunit_out, 1, dcd_out)
              end if

            end if
          end if
        end do

        call dcd_close(iunit_in, option%trj_in(itrj))

      end do

      call dcd_close(iunit_out, option%trj_out)


    end subroutine dcd2dcd
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
    subroutine xtc2xtc(option)
!-------------------------------------------------------------------------------
      implicit none

      type(s_option) :: option

      real(8), parameter :: ang2nm = 0.1d0

      integer :: nstep_tot
      integer :: iatm, id, itrj, istep, istep_tot, ixyz, nwrite 
      integer :: iunit_in, iunit_out
      integer :: natm, nstep, stride, first, last
      real    :: box(3, 3)
      logical :: is_end

      real, allocatable :: coord(:, :)

      type(s_traj)    :: traj
      type(s_trajopt) :: trajopt
      type(xtcfile)   :: xtc_in
      type(xtcfile)   :: xtc_out


      if (option%selfile_exist) then
        trajopt%molinfo(1) = trim(option%selfile)
        trajopt%nmolinfo   = 1
        call setup_traj_from_args(trajopt, 1, traj, trajid = 1)
        natm   = traj%natm

        allocate(coord(1:3, natm))
      else
        natm   = option%natm
      end if

      stride = option%stride

      first  = option%first
      last   = option%last


      nstep_tot = 0 
      if (.not. option%out_exist) then
        call get_total_step_from_xtc(option%trj_in, nstep_tot)         
        write(iw,'("Total number of steps is ",i0)') nstep_tot
        return
      end if

      call xtc_out%init(option%trj_out, 'w')

      iunit_in = 10

      istep_tot = 0
      nwrite    = 0
      do itrj = 1, option%ntrj_in

        call xtc_in%init(option%trj_in(itrj), 'r')

        is_end = .false.
        do while (.not. is_end)

          call xtc_in%read

          if (xtc_in%STAT /= 0) then
            is_end = .true.
            exit
          end if

          istep_tot = istep_tot + 1

          if (mod(istep_tot, 100) == 0) then
            write(iw,'("Read step ", i0)') istep_tot
          end if

          if (istep_tot >= first) then
            if (istep_tot <= last .or. last == 0) then

              if (mod(istep_tot, stride) == 0) then
                nwrite = nwrite + 1

                box = 0.0e0
                if (option%rect) then
                  do ixyz = 1, 3
                    box(ixyz, ixyz) = xtc_in%box(ixyz, ixyz)
                  end do
                else
                  box(:, :) = xtc_in%box(:, :)
                end if

                if (option%selfile_exist) then
                  do iatm = 1, natm
                    id = traj%ind(iatm)
                    coord(1:3, iatm) = xtc_in%pos(1:3, id)
                  end do

                  call xtc_out%write(natm,          &
                                     nwrite,        &
                                     xtc_in%time,   &
                                     box,           &
                                     coord,         &
                                     xtc_in%prec)
                else
                  call xtc_out%write(xtc_in%natoms, &
                                     nwrite,        &
                                     xtc_in%time,   &
                                     box,           &
                                     xtc_in%pos,    &
                                     xtc_in%prec)
                end if



              end if
            end if

          end if

          if (xtc_in%STAT /= 0) then 
            is_end = .true.
          end if

        end do

        call xtc_in%close

      end do

      write(iw,'("Total number of frames in output : ", i0)') nwrite


      if (allocated(coord)) &
        deallocate(coord)

      call xtc_out%close

    end subroutine xtc2xtc
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
    subroutine dcd2netcdf(option)
!-------------------------------------------------------------------------------
      implicit none

      type(s_option) :: option

      integer :: iatm, id, itrj, istep, istep_tot, nwrite 
      integer :: iunit_in, iunit_out
      integer :: natm, nstep, stride, first, last
      integer :: icount, nstep_tot

      ! For NetCDF
      !
      integer :: dim_frame, dim_spatial, dim_atom
      integer :: var_coords, var_box, var_angle, retval
      integer :: start(3), count(3)
      integer :: start_box(2), count_box(2)
      real(8) :: box(3, 1), angle(3, 1)

      real(8), allocatable :: coord(:, :, :)

      type(s_traj)    :: traj
      type(s_trajopt) :: trajopt
      type(s_dcd)     :: dcd_in


      if (option%selfile_exist) then
        trajopt%molinfo(1) = trim(option%selfile)
        trajopt%nmolinfo   = 1
        call setup_traj_from_args(trajopt, 1, traj, trajid = 1)
        natm   = traj%natm
      else
        natm   = option%natm
      end if

      allocate(coord(3, natm, 1))

      count     = (/3, natm, 1/)
      count_box = (/3, 1/)

      stride = option%stride
      first  = option%first
      last   = option%last

      iunit_in  = 10
      iunit_out = 11

      ! Get total number of steps in input trajectories
      !
      call get_total_step_from_dcd(option%trj_in, nstep_tot)

      write(iw,'("Total number of steps is ",i0)') nstep_tot

      if (.not. option%out_exist) &
        return

      nwrite = 0
      do istep = 1, nstep_tot
        if (istep >= first) then
          if (istep <= last .or. last == 0) then
            nwrite = nwrite + 1 
          end if
        end if
      end do

      istep_tot = 0
      nwrite    = 0

      call netcdf_open(option%trj_out, iunit_out , is_write = .true.)

      ! - Define dimensions
      !
      retval = nf90_def_dim(iunit_out, "frame",   nf90_unlimited, dim_frame)
      retval = nf90_def_dim(iunit_out, "spatial", 3,              dim_spatial)
      retval = nf90_def_dim(iunit_out, "atom",    natm,           dim_atom)

      ! - Define coordinate
      !
      retval = nf90_def_var(iunit_out, "coordinates",  nf90_real, (/dim_spatial, dim_atom, dim_frame/), var_coords) 
      retval = nf90_def_var(iunit_out, "cell_lengths", nf90_real, (/dim_spatial, dim_frame/),           var_box)
      retval = nf90_def_var(iunit_out, "cell_angles",  nf90_real, (/dim_spatial, dim_frame/),           var_angle)

      retval = nf90_put_att(iunit_out, var_coords,  "units",             "angstrom")
      retval = nf90_put_att(iunit_out, var_box,     "units",             "angstrom")
      retval = nf90_put_att(iunit_out, var_angle,   "units",             "degree")
      retval = nf90_put_att(iunit_out, nf90_global, "Conventions",       "AMBER")
      retval = nf90_put_att(iunit_out, nf90_global, "ConventionVersion", "1.0")
      retval = nf90_put_att(iunit_out, nf90_global, "program",           "ANATRA")
      retval = nf90_put_att(iunit_out, nf90_global, "programVersion",    "1.0")

      retval = nf90_enddef(iunit_out)


      icount = 0
      do itrj = 1, option%ntrj_in
        call dcd_open(option%trj_in(itrj), iunit_in)
        call dcd_read_header(iunit_in, dcd_in)
        nstep = dcd_in%nstep

        if (itrj == 1) then
          call alloc_dcd(option%natm, 1, dcd_in)
        end if

        do istep = 1, nstep
          istep_tot = istep_tot + 1

          if (mod(istep_tot, 100) == 0) then
            write(iw,'("Read step ", i0)') istep_tot
          end if

          call read_dcd_oneframe(iunit_in, dcd_in)

          if (istep_tot >= first) then
            if (istep_tot <= last .or. last == 0) then

              if (mod(istep_tot, stride) == 0) then
                if (option%selfile_exist) then
                  do iatm = 1, natm
                    id = traj%ind(iatm)
                    coord(1:3, iatm, 1) = dcd_in%coord(1:3, id, 1)
                  end do
                else
                  coord(1:3, 1:natm, 1) = dcd_in%coord(1:3, 1:natm, 1)
                end if
             
                icount    = icount + 1 
                start     = (/1, 1, icount/)
                start_box = (/1, icount/)

                box(1:3, 1)   = dcd_in%box(1:3, 1)
                angle(1:3, 1) = 90.0d0 
                retval =  nf90_put_var(iunit_out, var_coords, coord(1:3, 1:natm, 1), start = start,     count = count)
                retval =  nf90_put_var(iunit_out, var_box,    box(1:3, 1),           start = start_box, count = count_box) 
                retval =  nf90_put_var(iunit_out, var_angle,  angle(1:3, 1),         start = start_box, count = count_box) 
              end if

            end if
          end if
        end do

        call dcd_close(iunit_in, option%trj_in(itrj))

      end do

      call netcdf_close(iunit_out)


    end subroutine dcd2netcdf
!-------------------------------------------------------------------------------
    
!-------------------------------------------------------------------------------
    subroutine netcdf2netcdf(option)
!-------------------------------------------------------------------------------
      implicit none

      type(s_option) :: option

      integer :: iatm, id, itrj, istep, istep_tot, nwrite 
      integer :: iunit_in, iunit_out
      integer :: natm, nstep, stride, first, last
      integer :: icount, nstep_tot

      ! For NetCDF
      !
      integer :: dim_frame, dim_spatial, dim_atom
      integer :: var_coords, var_box, var_angle, retval
      integer :: start(3), count(3)
      integer :: start_box(2), count_box(2)
      real(8) :: box(3, 1), angle(3, 1)

      real(8), allocatable :: coord(:, :, :)

      type(s_traj)    :: traj
      type(s_trajopt) :: trajopt
      type(s_netcdf)  :: nc_in 


      if (option%selfile_exist) then
        trajopt%molinfo(1) = trim(option%selfile)
        trajopt%nmolinfo   = 1
        call setup_traj_from_args(trajopt, 1, traj, trajid = 1)
        natm   = traj%natm
      else
        natm   = option%natm
      end if

      allocate(coord(3, natm, 1))

      count     = (/3, natm, 1/)
      count_box = (/3, 1/)

      stride = option%stride

      first  = option%first
      last   = option%last

      iunit_in  = 10
      iunit_out = 11

      ! Get total number of steps in input trajectories
      !
      call get_total_step_from_netcdf(option%trj_in, nstep_tot)

      write(iw,'("Total number of steps is ",i0)') nstep_tot

      if (.not. option%out_exist) &
        return

      nwrite = 0
      do istep = 1, nstep_tot
        if (istep >= first) then
          if (istep <= last .or. last == 0) then
            nwrite = nwrite + 1 
          end if
        end if
      end do

      istep_tot = 0
      nwrite    = 0

      call netcdf_open(option%trj_out, iunit_out , is_write = .true.)

      ! - Define dimensions
      !
      retval = nf90_def_dim(iunit_out, "frame",   nf90_unlimited, dim_frame)
      retval = nf90_def_dim(iunit_out, "spatial", 3,              dim_spatial)
      retval = nf90_def_dim(iunit_out, "atom",    natm,           dim_atom)

      ! - Define coordinate
      !
      retval = nf90_def_var(iunit_out, "coordinates",  nf90_real, (/dim_spatial, dim_atom, dim_frame/), var_coords) 
      retval = nf90_def_var(iunit_out, "cell_lengths", nf90_real, (/dim_spatial, dim_frame/),           var_box)
      retval = nf90_def_var(iunit_out, "cell_angles",  nf90_real, (/dim_spatial, dim_frame/),           var_angle)

      retval = nf90_put_att(iunit_out, var_coords,  "units",             "angstrom")
      retval = nf90_put_att(iunit_out, var_box,     "units",             "angstrom")
      retval = nf90_put_att(iunit_out, var_angle,   "units",             "degree")
      retval = nf90_put_att(iunit_out, nf90_global, "Conventions",       "AMBER")
      retval = nf90_put_att(iunit_out, nf90_global, "ConventionVersion", "1.0")
      retval = nf90_put_att(iunit_out, nf90_global, "program",           "ANATRA")
      retval = nf90_put_att(iunit_out, nf90_global, "programVersion",    "1.0")

      retval = nf90_enddef(iunit_out)


      icount = 0
      do itrj = 1, option%ntrj_in

        call netcdf_open(option%trj_in(itrj), iunit_in)
        call netcdf_read_dimension(iunit_in, nc_in)
        nstep = nc_in%nstep

        do istep = 1, nstep
          istep_tot = istep_tot + 1

          if (mod(istep_tot, 100) == 0) then
            write(iw,'("Read step ", i0)') istep_tot
          end if


          if (istep_tot >= first) then
            if (istep_tot <= last .or. last == 0) then

              if (mod(istep_tot, stride) == 0) then

                call netcdf_read_oneframe(iunit_in, istep, nc_in)

                if (option%selfile_exist) then
                  do iatm = 1, natm
                    id = traj%ind(iatm)
                    coord(1:3, iatm, 1) = nc_in%coord(1:3, id, 1)
                  end do
                else
                  coord(1:3, 1:natm, 1) = nc_in%coord(1:3, 1:natm, 1)
                end if
             
                icount    = icount + 1 
                start     = (/1, 1, icount/)
                start_box = (/1, icount/)

                box(1:3, 1)   = nc_in%box(1:3, 1)
                angle(1:3, 1) = 90.0d0
                write(iw,'("WRITE")') 
                retval =  nf90_put_var(iunit_out, var_coords, coord(1:3, 1:natm, 1), start = start,     count = count)
                retval =  nf90_put_var(iunit_out, var_box,    box(1:3, 1),           start = start_box, count = count_box) 
                retval =  nf90_put_var(iunit_out, var_angle,  angle(1:3, 1),         start = start_box, count = count_box) 
              end if

            end if
          end if
        end do

        call netcdf_close(iunit_in)

      end do

      call netcdf_close(iunit_out)


    end subroutine netcdf2netcdf
!-------------------------------------------------------------------------------

end module mod_catcrd_analyze 
!===============================================================================
