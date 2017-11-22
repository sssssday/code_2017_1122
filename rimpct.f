

!!... RTIMP.F

      real function rtimp(riint, rictrl, rootc)
      use parm

      implicit none

!!... Argument declarations
      real      riint, rictrl, rootc

!!... This function calculates and returns a value between 0-1 which is the impact
!!... of root biomass on available nutrients.  It is used in the calculation of
!!... total plant production in RESTRP.

!!... Called From:     GROWTH
!!...                  CROPDYNC
!!...                  TREEDYNC
!!...                  TREEGROW

!!... Check added to handle underflow potential in exp intrinsic
      if ((rictrl * rootc * 2.5) .gt. 33) then
        rtimp = 1.0
      else
        rtimp = (1.0 - riint * exp(-rictrl * rootc * 2.5))
      endif

      return
      end
