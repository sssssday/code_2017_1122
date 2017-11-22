
c               Copyright 1993 Colorado State University
c                       All Rights Reserved


!!.... GPDF.F

      real function gpdf(x,a,b,c,d)

      implicit none

!!.... Argument declarations
      real    x, a, b, c, d

!!.... ******************** flowlib *********************
!!....
!!.... (run-time sub-set of modaid, exclusive of modctl)
!!....
!!.... Release 1.0  (first formal release of modaid)
!!....
!!....   james m. vevea
!!....   natural resource ecology lab
!!....   colorado state university
!!....   fort collins, colorado  80523
!!....
!!.... This routine is functionally equivalent to the routine of the
!!.... same name, described in the publication:
!!....
!!....   Some Graphs and their Functional Forms 
!!....   Technical Report No. 153
!!....   William Parton and George Innis (1972)
!!....   Natural Resource Ecology Lab.
!!....   Colorado State University
!!....   Fort collins, Colorado  80523

!!.... 12/90 Corrected by McKeown - exponent on frac changed from d to c

!!.... Local variables
      real    frac

      frac = (b-x) / (b-a)
      gpdf = 0.
      if (frac .gt. 0.) then
        gpdf = exp(c/d * (1. - frac**d)) * (frac**c)
      endif

      return
      end
