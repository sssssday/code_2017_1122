

      real function leafa(lai, rleavc, fbrchc, rlwodc, 
     &  cprodfLeft, cprodf)
      
      use parm

      implicit none
      

!!... Argument declarations
      real      rleavc, fbrchc, rlwodc, lai
      real      cprodfLeft, cprodf

!!... Compute the fraction of production going to leaves in crops
!!... and woody plants based on water and nutrient availability.
!!...   cprodf     - total potential C production for all tree parts
!!...                (gC/m^2)
!!...   cprodfLeft - amount of carbon still available for allocation
!!...                (gC/m^2)
!!...   leafprod   - amount of leaf production required to reach
!!...                optimal LAI
!!...   rleavc     - tree leaf carbon (gC/m^2)
!!...   fbrchc     - tree fine branch carbon (gC/m^2)
!!...   rlwodc     - tree large wood carbon (gC/m^2)

!!... Local Variables
      real tblai
      real leafprod
      real rleavc_opt
      integer j, idf
      
      j = 0
      idf = 0
          
      j = ihru
      idf = idplt(j)

      if (rleavc .lt. 0.0) then
        write(*,*) 'Error in leafa, rleavc < 0.0'
        write(*,*) 'rleavc = ', rleavc
        STOP
      endif
      if (rlwodc .lt. 0.0) then
        write(*,*) 'Error in leafa, rlwodc < 0.0'
        STOP
      endif
      if (cprodfLeft .lt. 0.0) then
        write(*,*) 'Error in leafa, cprodfLeft < 0.0'
        STOP
      endif
      if (cprodf .le. 0.0) then
        write(*,*) 'Error in leafa, cprodf <= 0.0'
        STOP
      endif

!!... Calculate theoretical maximum for LAI based on large wood biomass,
!!... cak - 07/24/02
!!... Include fine branch carbon as part of the woody component in the
!!... LAI calculation, cak - 10/20/2006
c      call lacalc(lai, rleavc, rlwodc, btolai, maxlai, klai)
c      call lacalc(lai, rlwodc, maxlai, klai)
      !!call lacalc(lai, fbrchc, rlwodc, maxlai, klai)
      if (lai .lt. 0.1) then
        lai = 0.1
      endif
          
      if (BTOLAIf(idf) .le. 0.0) then
        write(*,*) 'Error in leafa, btolai <= 0.0, check tree.100'
        STOP
      endif
!!... Calculate optimal rleavc
      if (lai > MAXLAIf(idf)) then
      
        leafa = 0.
        go to 999         
      
      endif
      
      rleavc_opt = 2.5 * lai / BTOLAIf(idf)
      print*, 'leafopt=', rleavc_opt

      if (rleavc_opt .gt. rleavc) then
!!..... if possible, increase leaf biomass so that optimum is reached
        leafprod = min((rleavc_opt-rleavc), cprodfLeft)
      else
!!.....  optimum leaf area has already been acheived
        leafprod = 0.0
      endif

      leafa = leafprod / cprodf

      if (leafa .lt. 0.0) then
        write(*,*) 'Error in leafa, leafa < 0.0'
        write(*,*) 'leafa = ', leafa
        STOP
      endif
      if (leafa .gt. 1.0) then
        write(*,*) 'Error in leafa, leafa > 1.0'
        write(*,*) 'leafa = ', leafa
        STOP
      endif
999   return
      end
