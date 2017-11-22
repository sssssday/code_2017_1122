
!               Copyright 1993 Colorado State University
!                       All Rights Reserved


      subroutine csched(cflow,
     &                  froma1,tob1,
     &                  froma2,tob2,
     &                  frac,accum)

      implicit none
      !include 'const.inc'
      !include 'zztim.inc'

! ... Argument declarations
      real    cflow, froma1, tob1, froma2, tob2,
     &        frac,accum(2)

! ... Schedule C flows for decomposition from Box A to Box B
! ... written by vek 05/91

! ... Input:
! ...   cflow  = total C flow from Box A to Box B
! ...   protop = proportion top, either the labeled pool 
! ...            OR cisofr or cisotf
! ...   probot = proportion bottom, either the total pool or 1
! ...   frac   = amount of fractionation
! ...   time   = simulation time (passed in /zztim/)
! ...
! ... Transput:
! ...   froma1   = state variable for unlabeled C in Box A
! ...   tob1     = state variable for unlabeled C in Box B
! ...   froma2   = state variable for labeled C in Box A
! ...   tob2     = state variable for labeled C in Box B
! ...   accum(1) = unlabeled accumulator
! ...   accum(2) = labeled accumulator

! ... Local variables
      real    cflow1, cflow2

! ... Determine the amount of labeled to flow
      cflow2 = cflow * (protop / probot) * frac

! ... Determine the amount of unlabeled to flow
      cflow1 = cflow - cflow2

! ... Flow the amounts
      call flow(froma1,tob1,time,cflow1)
      call flow(froma2,tob2,time,cflow2)

! ... Accumulate the amounts flowed
      !accum(UNLABL) = accum(UNLABL) + cflow1
      !accum(LABELD) = accum(LABELD) + cflow2

      return
      end
