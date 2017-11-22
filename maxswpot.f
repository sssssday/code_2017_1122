


      real function maxswpot(numlayers)
      use parm

      implicit none
!!      include 'parfx.inc'
 !!     include 'plot1.inc'
 !!     include 'site.inc'

!!... Argument declarations
      integer numlayers

!!... This function calculates the soil water potential of the wettest
!!... soil layer in the plant rooting zone.

!!... Declaration explanations:
!!...   numlayers - number of soil layers in the plant rooting zone

!!... Local variable explanations:
!!...   b      - slope of retention curve
!!...   BAR2CM - conversion factor for bars to centimeters H2O
!!...   base   - base value for power function
!!...   expon  - exponent value for power function
!!...   lyr    - current soil layer
!!...   psis   - "saturation" matric potential of "ilyr" (cm H2O ?)
!!...   swptnl - soil water potential of the current layer (bars)
!!...   theta  - volumetric soil water content * 100
!!...   thetas - volumetric soil water content at saturation for layer
!!...            (% volume)

!!... Local variables
      integer          BAR2CM, lyr,j
      real             b, psis, swptnl, theta, thetas
      double precision base, expon
      real :: sol_thick(sol_nly(ihru))
      
      j = ihru
!!... 
      BAR2CM = 1024

      maxswpot = 9999
      
      do 10 lyr = 1, numlayers 
      
        
        
           if (lyr == 1) then
	              sol_thick(lyr) = sol_z(lyr,j)
	        else	
	              sol_thick(lyr) = sol_z(lyr,j) - sol_z(lyr-1,j)
	        end if
       
10      end do 

      do 130 lyr = 1, numlayers  
!!..... Calculate the soil water potential for the current soil layer

        if (sol_st(lyr,j) .gt. 0) then
          theta = (sol_thick(lyr)*0.1) / (sol_thick(lyr)*0.1) * 100  !! *0.1 to convert mm to cm
          thetas = (-14.2 * sol_sand(lyr,j)*0.01) -           
     & (3.7 * sol_clay(lyr,j)*0.01) + 50.5
          base =  theta / thetas
          b = (-0.3 * sol_sand(lyr,j)*0.01) +      
     &  (15.7 * sol_clay(lyr,j)*0.01) + 3.10
          expon = b
          psis = 10.0**((-1.58 * sol_sand(lyr,j)*0.01) - 
     &      (0.63 * sol_clay(lyr,j)*0.01) + 2.17)
          swptnl = (psis / (base**expon)) / BAR2CM
        else
          swptnl = 80.0
        endif

        if (swptnl .lt. maxswpot) then
          maxswpot = swptnl   !! get the maximum potential for the given layers
        endif
130   continue

!!... Place a limit the maximum water potential of the wettest layer
!!... because roots are not able to extract water below this level,
!!... cak - 03/20/2008

   
      if (maxswpot .gt. 30.0) then
        maxswpot = 30.0
      endif

      return
      end
