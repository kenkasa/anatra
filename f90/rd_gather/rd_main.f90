!=======================================================================
program main 
!=======================================================================
  use mod_const
  use mod_ctrl
  use mod_analyze

  type(s_input)   :: input
  type(s_output)  :: output
  type(s_rdinfo)  :: rdinfo 

  call show_title
  call show_usage
  call read_ctrl(input, output, rdinfo)
  call analyze(input, output, rdinfo)

end program main 
!=======================================================================

!-----------------------------------------------------------------------
subroutine show_title
!-----------------------------------------------------------------------
  implicit none

  write(6,'("==================================================")')
  write(6,*)
  write(6,'("                 RDF Gathering")')
  write(6,*)
  write(6,'("==================================================")') 

end subroutine show_title
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine show_usage
!-----------------------------------------------------------------------
  use mod_const

  implicit none

  character(len=MaxChar) :: f_ctrl


  call getarg(1, f_ctrl)

  if (trim(f_ctrl) == "-h") then
    write(iw,'("&input_param")')
    write(iw,'(" flist = ""rdbin_filelist"" ! File that contains List of rdbin")')
    write(iw,'("/")')

    write(iw,'("&output_param")')
    write(iw,'(" frd   = ""out""            ! RDF file")')
    write(iw,'("/")')
    stop
  end if

end subroutine show_usage
!-----------------------------------------------------------------------
