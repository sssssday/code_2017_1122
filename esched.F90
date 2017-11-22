


      subroutine esched(cflow,tca,rcetob,anps,bnps,labile)

      implicit none
    

!! ...... Argument declarations
      real     cflow, tca, rcetob, anps, bnps, labile
      real     mnrflo

!! ...... Schedule N, P, or S flow and associated mineralization or
!! ...... immobilization flow for decomposition from Box A to Box B.
!! ...... written by vek 05/91

!! ...... Input:
!! ......   cflow  = C flow from Box A to Box B
!! ......   tca    = total C (unlabeled + labeled) in Box A
!! ......   rcetob = C/N, C/P, or C/S ratio of new material being added
!! ......            to Box B
!! ......   time   = simulation time, passed in /zztim/
!! ......
!! ...... Transput:
!! ......   anps   = N, P, or S state variable for Box A
!! ......   bnps   = N, P, or S state variable for Box B
!! ......   labile = minerl(1,iel) where iel indicates N, P, or S
!! ......
!! ...... Output:
!! ......   mnrflo = amount of N, P, or S that is mineralized or
!! ......            immobilized.
!! ......
!! ...... A positive value of mnrflow indicates mineralization;
!! ...... a negative value of mnrflow indicates immobilization.

!! ...... Local variables
      real     atob, immflo, outofa

!! ...... Compute and schedule N, P, and S flows

!! ...... N, P, or S flowing out of Box A is proportional to C flow.
      outofa = anps * (cflow/tca)

!! ...... Microcosm option can cause a 0/0 error on the pc.  This
!! ...... was added to avoid that situation (mse 2/95)
!! ...... Change from .AND. to .OR. could still get an underflow.
!! ...... -rm 1/96

      if (cflow .le. 0.0 .or. outofa .le. 0.0) then
        mnrflo = 0.0
        goto 999
      endif

!! ...... If C/E of Box A > C/E of new material entering Box B
      if (cflow/outofa .gt. rcetob) then

!! ........ IMMOBILIZATION occurs.

!! ........ Compute the amount of E immobilized.
!! ........ since  rcetob = cflow/(outofa+immflo),
!! ........ where immflo is the extra E needed from the mineral pool
        immflo = cflow/rcetob - outofa

!! ........ Remove the code that is turning of the E flow due to
!! ........ decomposition when the amount of E in the labile pool falls
!! ........ below zero, cak - 01/26/2011
!! ........ Do not allow immobilization to occur if there is not enough
!! ........ mineral E, cak - 09/10/02
!!       if ((labile - immflo) .lt. 0.0) then
!!         mnrflo = 0.0
!!         goto 999
!!       endif

!! ........ Schedule flow from Box A to Box B (outofa)
       !! call flow(anps,bnps,time,outofa)

           anps = anps -outofa
           bnps = bnps + outofa*10                     !Need to convert to kg/ha Zhang
!! ........ Schedule flow from mineral pool to Box B (immflo)
        !!
        !!   call flow(labile,bnps,time,immflo)
             labile = labile - immflo*10              !Need to convert to kg/ha Zhang !! what if the soil mineral pool could not privide needed nutrient????
             bnps = bnps + immflo                       !!Need to convert to kg/ha Zhang

!! ........ Return mineralization value.
       !! mnrflo = -immflo
          
      else
!! ........ MINERALIZATION occurs
!! ........ Schedule flow from Box A to Box B
        atob = cflow/rcetob 
        !!call flow(anps,bnps,time,atob)
        
           anps = anps - atob
           bnps = bnps + atob*10                       !!!Need to convert to kg/ha Zhang
          
!! ........ Schedule flow from Box A to mineral pool
        mnrflo = outofa - atob 
       !! call flow(anps,labile,time,mnrflo)
        anps = anps - mnrflo
        labile = labile + mnrflo*10                    !!Need to convert to kg/ha Zhang
        

      endif

999   continue

      return
      end
