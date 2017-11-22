
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


! ... MNRACC.F

      subroutine mnracc (mnrflo,gross,net)

      implicit none

! ... Argument declarations
      real      mnrflo, gross, net

! ... Update mineralization accumulators.
! ... written by vek 05/91

! ... Input:
! ...   mnrflo = mineralization value returned by esched.
! ...            A negative value indicates immobilization.
! ...     
! ... Transput:
! ...   gross = gross mineralization (sums only mineralization, not
! ...           immoblization)
! ...   net   = net mineralization (mineralization - immobilization)

! ... Gross mineralization
      if (mnrflo .gt. 0.0) then
        gross = gross + mnrflo
      endif

! ... Net mineralization
      net = net + mnrflo

      return
      end
