
     
      real function anerob (drain, rprpet, pet, micosm)

      use parm
      
      implicit none
      
      integer j, idf
      

!! ...... Argument declarations
      real     drain, rprpet, pet
      integer  micosm
 
!! ...... This function calculates the impact of soil anaerobic conditions
!! ...... on decomposition.  It returns a multiplier 'anerob' whose value
!! ...... is 0-1.

!! ...... If a microcosm is being simulated, return the value for aneref(3)
!! ...... which is set by the user.

!! ...... Called From:  CYCLE

!! ...... Declaration explanations:
!! ......   aneref(1) - ratio RAIN/PET with maximum impact
!! ......   aneref(2) - ratio RAIN/PET with minimum impact
!! ......   aneref(3) - minimum impact
!! ......   drain     - percentage of excess water lost by drainage
!! ......   newrat    - local var calculated new (RAIN+IRRACT+AVH2O(3))/PET ratio
!! ......   pet       - potential evapotranspiration
!! ......   rprpet    - actual (RAIN+IRRACT+AVH2O(3))/PET ratio

!! ...... Local variables
      real      newrat, tslope, xh2o
      real      aneref(3)
      
      aneref(1) = 1.5
      aneref(2) = 3.0
      aneref(3) = 1.0
      anerob = 1.0



!! ...... Determine if RAIN/PET ratio is greater than the ratio with
!! ...... maximum impact.
!! ...... If it is winter time the value for rprpet has been calculated
!! ...... using snow that is melting into the soil profile, cak - 10/21/02
      if (rprpet .gt. aneref(1)) then
        xh2o = (rprpet - aneref(1)) * pet * (1.0 - drain)
        if (xh2o .gt. 0) then
          newrat = aneref(1) + (xh2o / pet)
          tslope = (1.0 - aneref(3)) / (aneref(1) - aneref(2))
          anerob = 1.0 + tslope * (newrat - aneref(1))
        endif

        if (anerob .lt. aneref(3)) then
          anerob = aneref(3)
        endif
      endif


      return
      end