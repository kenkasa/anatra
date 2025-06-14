!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_ctrl
  use mod_traj
  use mod_com

  ! subroutines
  !
  public :: analyze 

  contains
!-----------------------------------------------------------------------
    subroutine analyze(output, option, traj)
!-----------------------------------------------------------------------
      implicit none

      type(s_output), intent(in)   :: output
      type(s_option), intent(in)   :: option
      type(s_traj),   intent(in)   :: traj(2)

      type(s_com)            :: com(2)

      integer                :: ixyz, itraj
      integer                :: istep, jstep, ij, imol
      integer                :: nend
      integer                :: nmol, nstep
      real(8)                :: cost, vec(3), veclength
      character(len=MaxChar) :: fname

      real(8), allocatable :: orient_vec(:, :, :)
      real(8), allocatable :: tcf(:, :), tcfave(:)


      call get_com(option%mode, traj(1), com(1))
      call get_com(option%mode, traj(2), com(2))

      if (com(1)%nmol /= com(2)%nmol) then
        write(iw,'("Analyze> Error.")')
        write(iw,'("Number of molecules in two dcd should be the same.")')
        stop
      end if

      nmol  = com(1)%nmol
      nstep = traj(1)%nstep 

      ! allocate memory
      !
      allocate(orient_vec(3, nstep, nmol))
      allocate(tcf(0:option%tcfrange, nmol))
      allocate(tcfave(0:option%tcfrange))

      ! calculate orientation vector 
      !
      orient_vec = 0.0d0
      do istep = 1, nstep
        do imol = 1, nmol
          vec(1:3)  = com(2)%coord(1:3, imol, istep) &
                    - com(1)%coord(1:3, imol, istep)

          veclength = sqrt(dot_product(vec(1:3), vec(1:3)))
          orient_vec(1:3, istep, imol) = vec(1:3) / veclength
        end do
      end do

      ! calculate tcf
      !
      tcf    = 0.0d0
      tcfave = 0.0d0

!$omp parallel private(imol, istep, nend, ij, cost) &
!$omp        & default(shared)
!$omp do
      do imol = 1, nmol
        do istep = 1, nstep - 1

          nend = istep + option%tcfrange + 1
          if (nend > nstep) &
            nend = nstep

          do jstep = istep, nend
            ij            = jstep - istep
            cost          = dot_product(orient_vec(1:3, jstep, imol), &
                                        orient_vec(1:3, istep, imol))
            tcf(ij, imol) = tcf(ij, imol) + cost * cost 
                            
          end do
        end do 
      end do
!$omp end do
!$omp end parallel

      do imol = 1, nmol
        tcf(0, imol) = 1.0d0
        do istep = 1, option%tcfrange
          tcf(istep, imol) = tcf(istep, imol) / (nstep - istep)
          tcf(istep, imol) = 1.5d0 * tcf(istep, imol) - 0.5d0
        end do
      end do

      ! calculate averaged tcf
      !
      do istep = 0, option%tcfrange - 1
        tcfave(istep) = sum(tcf(istep, 1:nmol)) / dble(nmol)
      end do

      ! generate files
      !
      write(fname,'(a,".tcf")') trim(output%fhead)
      open(UnitOut, file=trim(fname))
      do istep = 0, option%tcfrange - 1
        write(UnitOut, '(f20.10)', advance='no') traj(1)%dt * istep
        do imol = 1, nmol
          write(UnitOut,'(f20.10)') tcf(istep, imol)
        end do
        write(UnitOut, *)
      end do
      close(UnitOut)

      write(fname,'(a,".tcfave")') trim(output%fhead)
      open(UnitOut, file=trim(fname))
      do istep = 0, option%tcfrange - 1
        write(UnitOut, '(2f20.10)') traj(1)%dt * istep, tcfave(istep)
      end do
      close(UnitOut)

      ! Deallocate memory
      !
      deallocate(orient_vec)
      deallocate(tcf, tcfave)

    end subroutine analyze
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
