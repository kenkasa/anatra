!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_ctrl

  ! structures
  type :: s_rdbin
    logical :: identical
    logical :: separate_self
    integer :: nstep
    integer :: nmol(2)
    real(8) :: vol
    real(8) :: dr
    integer :: nr
    integer, allocatable :: hist(:), hist_distinct(:), hist_self(:)
  end type s_rdbin

  ! subroutines
  !
  public :: analyze 

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, rdinfo)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in)   :: input
      type(s_output), intent(in)   :: output
      type(s_rdinfo), intent(in)   :: rdinfo

      type(s_rdbin), allocatable :: rdbin(:)

      integer :: ig
      integer :: nr, nstep, ifile
      integer :: nmol(2)
      integer :: npair, npair_distinct
      real(8) :: dr, vol
      real(8) :: r, fourpi, fact
      logical :: identical, separate_self

      real(8), allocatable :: gr(:), gr_self(:), gr_distinct(:)


      allocate(rdbin(rdinfo%nfile))
      do ifile = 1, rdinfo%nfile
        open(10,file=trim(rdinfo%frdbin(ifile)), form='unformatted')
       
        read(10) rdbin(ifile)%identical
        read(10) rdbin(ifile)%separate_self
        read(10) rdbin(ifile)%nstep
        read(10) rdbin(ifile)%nmol
        read(10) rdbin(ifile)%vol
        read(10) rdbin(ifile)%dr
        read(10) rdbin(ifile)%nr

        nr = rdbin(ifile)%nr
        allocate(rdbin(ifile)%hist(0:nr))
        allocate(rdbin(ifile)%hist_distinct(0:nr))
        allocate(rdbin(ifile)%hist_self(0:nr))

        read(10) rdbin(ifile)%hist, &
                 rdbin(ifile)%hist_distinct, &
                 rdbin(ifile)%hist_self

        close(10)
      end do

      identical     = rdbin(1)%identical
      separate_self = rdbin(1)%separate_self
      nr            = rdbin(1)%nr
      dr            = rdbin(1)%dr
      nmol          = rdbin(1)%nmol
      nstep         = sum(rdbin(:)%nstep)
      vol           = sum(rdbin(:)%vol * rdbin(:)%nstep) / nstep

      write(iw,*)
      write(iw,'(" Analyze> Show parameters read from rdbin files")')
      write(iw,'(" identical     = ",a)') get_tof(identical)
      write(iw,'(" separate_self = ",a)') get_tof(separate_self)
      write(iw,'(" nr            = ",i0)') nr
      write(iw,'(" dr            = ",f20.10)') dr
      write(iw,'(" nmol          = ",i0)') nmol
      write(iw,'(" nstep         = ",i0)') nstep
      write(iw,'(" vol           = ",f20.10)') vol
      write(iw,*)


      allocate(gr(0:nr), gr_self(0:nr), gr_distinct(0:nr))

      if (identical) then
        npair = nmol(1) * (nmol(1) - 1) / 2
      else
        npair = nmol(1) * nmol(2)
      end if

      npair_distinct = nmol(1) * (nmol(1) - 1) / 2

      gr          = 0.0d0
      gr_self     = 0.0d0
      gr_distinct = 0.0d0
      do ig = 0, nr
        do ifile = 1, rdinfo%nfile
          gr(ig)          = gr(ig) + rdbin(ifile)%hist(ig)
          if (separate_self) then
            gr_self(ig)     = gr_self(ig)     + rdbin(ifile)%hist_self(ig)
            gr_distinct(ig) = gr_distinct(ig) + rdbin(ifile)%hist_distinct(ig)
          end if
        end do
      end do

      fourpi = 4.0d0 * PI
      do ig = 1, nr
        r = dble(ig) * dr
        fact = nmol(2) / (fourpi * r * r * dr * npair * nstep)
        gr(ig) = gr(ig) * fact

        if (separate_self) then
          fact = 1.0d0 / (fourpi * r * r * dr * nmol(1) * nstep)
          gr_self(ig) = gr_self(ig) * fact
          if (npair_distinct == 0) then
            gr_distinct(ig) = 0.0d0 
          else
            fact = 1.0d0 / (fourpi * r * r * dr * npair_distinct * nstep)
            gr_distinct(ig) = gr_distinct(ig) * fact
          end if
        end if
      end do

      open(10, file=trim(output%frd))
      if (separate_self) then
        do ig = 0, nr 
          r = dr * dble(ig) 
          write(10,'(4f20.10)') r, gr(ig), gr_distinct(ig), gr_self(ig) 
        end do
      else
        do ig = 0, nr 
          r = dr * dble(ig) 
          write(10,'(2f20.10)') r, gr(ig) 
        end do
      end if

      close(10)

      deallocate(gr, gr_self, gr_distinct)

    end subroutine analyze
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
