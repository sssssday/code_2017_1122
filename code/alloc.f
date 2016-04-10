
      subroutine alloc(gpp, tavg, swatlyr, tree_cfrac, checkgrow)
      use parm

      implicit none

!!    need to check wether should use swatlyr or all the soil layers. 
!! ... Argument declarations
      real gpp, tavg, tree_cfrac(6)
      integer checkgrow,  swatlyr

!! ... Compute carbon allocation fractions for tree parts (LEAF, FROOT,
!! ... RBRCH, LWOOD, and CROOT).
!! ...
!! ... For seasonal deciduous/conifer forest system, reapportion growth
!! ... This applies only to non-drought systems.
!! ... From Bill Parton (April 2000):
!! ... "For the drought decidious plants the main impact is to have leaves
!! ... drop in response to drought stress. You don't want leaf growth
!! ... initiated the same way as traditional decidious plants. Drought
!! ... decidious plants have their leaf growth start when it get wet enough
!! ... and is not controlled by temperature directly like the traditional
!! ... decidious plants."
!! ...

!! ... Function declarations
      integer grochk,j, idf, ii
      real       rtimp
      external   grochk, rtimp
      
      
!! ... Local variables
      integer :: greenUpCnt
      integer :: iel, lyr
      real,dimension(3) :: eavail
      real :: demand,  fracrc, leafprod,maxrootf
     &        maxNfix, rimpct, rootprod, totale, toler
      !! integer :: nelem !! how many elements considered in the limitation calculation
      !! nutrient demond based on carbon uptake and available nitrogen
      integer :: l
      
      real:: avn, avNH4, avNO3
      real:: h2ogef !! water stress on growth    
      real:: petdly !! potential daily ET (cm h2o/day)
      real :: twstress !! wetest water stress
      real :: maxrwcf 
      real, dimension(10) :: rwct !! relative water content in each soil layer
      real :: riint, rictrl
      !!real,dimension(:) :: tree_a2drat
      real,dimension(2) :: favail
      real :: maxrootf
      real, dimension(10) :: rwcf
      real, dimension(3) :: ta2drat
      
      data greenUpCnt /0/
      save greenUpCnt

      j = 0
      idf = 0
      
      tree_a2drat = 0.0
      
      j = ihru
      idf = idplt(j)
      
      h2ogef = 0.0
      favail = 0.


      toler = 1.0E-30
      twstress = 0.0

!! ... Use tree specific favail(1) value
      maxrootf = SFAVAIL2f(idf)

!! ... Estimate the fraction of carbon going to the roots
      fracrc = (TFRTCWf(idf,1)+TFRTCWf(idf,2)+TFRTCNf(idf,1)+
     &      TFRTCNf(idf,2))/4.0

!! ... Estimate the fine root and leaf production based on total
!! ... potential production
      leafprod = gpp - gpp * fracrc
      rootprod = gpp * fracrc



   !!   do 20 iel = 1, 1   !! only consider N and P currently. 1: N; 2: P
   !! only consider N and P currently. 1: N; 2: P
        iel = 1
        
        favail(1) = SFAVAIL2f(idf)
        
        !! calculation of favail(2) based on the first soil layer
        
        favail(2) = max(0.2,
     &              min(0.2 + (sol_no3(1,j)/10. + sol_nh3(1,j)/10.)*
     &              (0.4 - 0.2) / 2.0,
     &              0.4))
        
        
        
!! ..... Nutrients available to trees are in the top tlaypg layers,
!! ..... cak 01/29/03
!!        do 30 lyr = 1, nlayer
       !! do 30 lyr = 1, swatlyr    !! check if we need to use all layers
          !!if (minerl(lyr,iel) .gt. toler) then
            !availm(iel) = availm(iel) + minerl(lyr, iel)
      !!    endif
!!30      continue
!!20    continue

!!  .. Calculate impact of root biomass on available nutrients
      !! following two parameters are both from the fixed parameters of daycent
      riint = 0.8
      rictrl = 0.015

      rimpct = rtimp(riint, rictrl, frootjc(j)+frootmc(j))

!!  .. Calculate soil available nutrients, based on a maximum fraction
!!  .. (favail) and the impact of root biomass (rimpct), adding storage.
      do 45 iel = 1, nelem
        eavail(iel) = (availm(j,iel) * maxrootf * rimpct) + 
     &                 forstg(j,iel)
45    continue


       
      
      !!nelem = 1 !! here only consider N limitation 2, means consider both N and P, 3 means N,P,S
     !! Estimate the demand
      !! caluclate total available N in soils need to change later to account for P and S
     
     
      do 50  iel = 1, nelem      
      
        
     !! N FIXATION
        if (iel .eq. 1) then
          call nfix
        endif
     !! DEMAND based on the maximum E/C ratio.
     
     
        demand = 0.0
        demand = demand + leafprod * (1.0 / CERFORf(idf,IMIN,LEAF,iel))
        demand = demand + rootprod * (1.0 / CERFORf(idf,IMIN,FROOT,iel))
         if (iel .eq. 1) then
        totale = eavail(iel) + fixn
         endif

     !! New calculation -mdh 5/10/01
        tree_a2drat(j,iel) = min(1.0, totale / demand)
        tree_a2drat(j,iel) = max(0.0, tree_a2drat(j,iel))
        ta2drat(iel) = tree_a2drat(j,iel)
50    continue

!!...... calculate the soil water impacts
           petdly = pet_day/10.0
   !! calculate tree water stress  !! find the wetest soil layer
      maxrwcf = -9999
      do 140 ii = 1, swatlyr
           rwcf(ii) = (sol_st(ii,j)*0.1) / (sol_z(ii,j)*0.1)  !! relative wetness
        if (rwcf(ii) .gt. maxrwcf) then
          maxrwcf = rwcf(ii)
        endif
140   continue

      twstress = twstress + min(1.0, maxrwcf)        
      
        if (petdly .ge. .01) then
!! ....... Calculate potential growth based on the relative water content
!! ....... of the wettest soil layer, cak - 12/06/04
          h2ogef = 1.0/(1.0+exp(WSCOEFF2f(idf,2)*
     &                             (WSCOEFF2f(idf,1)-twstress)))
        else
          h2ogef = 0.01
        endif
      call froota(ta2drat,h2ogef)

!!  .. Decidious forest

      if (DECIDf(idf) .eq. 1) then
!!  .... Determine if we are within the greenup period for deciduous trees
        if (checkgrow .eq. 1) then   !! check to use greenUpCnt to control the green up period
          greenUpCnt = greenUpCnt + 1
!!  ...... All allocation can be to leaves during greenup
!!  ...... Allow some fine root growth during leaf out, cak - 03/27/2007
          tree_cfrac(LEAF) = 0.9
          tree_cfrac(FROOT) = 0.1
          if (greenUpCnt .gt. 30) then
            write(*,*) 'Error in treeDynC, greenUpCnt > 30'   !! check there is no processs to limit greenUpCnt to 30, 
            STOP
          endif
        elseif (decidgrow) then
!!  ...... If we are in the period between leaf out and leaf drop, allocate
!!  ...... to fine roots first and then to leaves
          greenUpCnt = 0
          tree_cfrac(FROOT) = ffroota
          tree_cfrac(LEAF) = 1.0 - tree_cfrac(FROOT)
        else
!!  ...... No growth occurs this time period
          gpp = 0.0   !! check how to determine this 
          goto 999
        endif
!!  .. Drought decidious forest
      else if (DECIDf(idf) .eq. 2) then
!!  .... Determine if we are within the greenup period for drought
!!  .... deciduous trees
        if (hrsinc .and. (h2ogef .gt. 0.5)) then
!!  ...... All allocation can be to leaves during greenup
!!  ...... Allow some fine root growth during leaf out, cak - 03/27/2007
          tree_cfrac(LEAF) = 0.9
          tree_cfrac(FROOT) = 0.1
        else
!!  ...... Allocate to fine roots first and then to leaves
          tree_cfrac(FROOT) = ffroota
          tree_cfrac(LEAF) = 1.0 - tree_cfrac(FROOT)
        endif
!!  .. Evergreen forest
      else if (DECIDf(idf) .eq. 0) then
!!  .... Allocate to fine roots first and then to leaves
       
        tree_cfrac(FROOT) = ffroota
        tree_cfrac(LEAF) = 1.0 - tree_cfrac(FROOT)
      else
        write(*,*) "Invalid forest type in treeDynC!"
        STOP
      endif

      tree_cfrac(FBRCH) = 0.0
      tree_cfrac(LWOOD) = 0.0
      tree_cfrac(CROOT) = 0.0

999   continue

      return
      end