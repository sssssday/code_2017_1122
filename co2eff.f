
      subroutine co2eff (ttime)
      
      use parm

      implicit none


!!...Argument declarations
      integer ttime

!!...Compute the effect of atmospheric CO2 concentration.

!!...modified
!!...K. Killian 21-Mar-2003  Make co2sys consistent with co2tm input
!!...K. Killian 15-Jul-1998  Changed the step function concentration
!!...so that it finds the upper concentration.

!!...Function declarations
      real    ::  effect, line, ramp
      external  line, ramp

!!...Local variables
      integer ::   iel, mnmx, system, j, idf
      real    ::  co2conc
      !!integer :: nelem
      
      j = 0
      idf = 0

      
      j = ihru
      idf = idplt(j)

!!...Reset all effects to 1.0
      !!do 30 system = CRPSYS, FORSYS  only consider forst
      
      !!  co2cpr(2) = 1.0
       !! co2ctr(system) = 1.0
        co2crs = 1.0
        !!nelem = 2 !!  only consider N and P
        do 20 mnmx = IMIN, IMAX
          do 10 iel = 1, nelem
            co2cce(2,mnmx,iel) = 1.0
10        continue
20      continue
!!30    continue

!!...If there is no co2 effect, return
      !!if (co2sys .lt. 0) goto 999
      !! if (co2sys .le. 0) goto 999

!!...Calculate the co2 effect

!!...Calculate the new co2 concentration
 !!     if (time .le. co2tm(1)) then
c ..... use base concentration if time < start of effect

!!      else if (co2rmp .eq. 0) then
c ..... raised Constant concentration
 !!       co2conc = co2ppm(2)

 !!     else
c ..... Ramping
 !!       if (time .lt. co2tm(2)) then
 !!         co2conc = ramp(time,co2tm(1),co2ppm(1),co2tm(2),co2ppm(2))
 !!       else
 !!         co2conc = co2ppm(2)
 !!       endif
 !!     endif
            
      co2conc = co2(hru_sub(j)) 
      
      co2crs = line(co2conc,350.0,1.0,700.0,CO2IRS2f(idf))


!!...Calculate effect on production
      !! do 40 system = CRPSYS, FORSYS
       !! co2cpr(system) = effect(co2ipr(system),co2conc)
!!40    ! continue

!!!!...Calculate effect on PET
!!      do 50 system = CRPSYS, FORSYS
!!        co2ctr(system) = effect(co2itr(system),co2conc)
!!50    continue

!!...Calculate effect on C/E
      do 80 iel = 1, nelem
        do 70 mnmx = IMIN, IMAX
!!          do 60 system = CRPSYS, FORSYS
            co2cce(2,mnmx,iel) =
     &        effect(CO2ICE2f(idf,mnmx,iel),co2conc)
!!60        continue
70      continue
80    continue

!!...Calculate effect on root/shoot
!!...Reference co2 concentration = 350.0 at 1.0
!!      do 90 system = CRPSYS, FORSYS
!!        co2crs(2) = line(co2conc,350.0,1.0,700.0,co2irs(system))
!!90    continue

999   continue

      return
      end


      real function effect(co2input, co2conc)

      implicit none
      real      co2input, co2conc

!!...Reference co2 concentration = 350.0
      effect = 1 + (co2input-1) / (log10(2.0)) * (log10(co2conc/350.0))

      return
      end
