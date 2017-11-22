
      subroutine wdeath (tavewk, bgwfunc, tfrac, avgstemp, thrsinc, ldrmlt,mrspweff, swatlyr)
        use parm
        implicit none
        
!!.......Argument declarations
      real      :: tavewk, bgwfunc, tfrac, avgstemp, ldrmlt
      integer   :: swatlyr
      logical   :: thrsinc
      
!!.......Death of leaves, fine branches, large wood, fine roots, and coarse roots.
!!.......Modifications:
!!.......Corrected a bug in the death of lableled fine branches, coarse wood, 
!!.......and coarse roots. Use the actual labled fraction instead of the growth 
!!.......ratio, cisotf. Prevents the removal of to much labeled material when 
!!.......non equilibrium, under ratio, wood dies.  7/2007  K. Killian

!!.......Function declarations
      real      catanf, gpdf, maxswpot, mrspweff
      external  catanf, gpdf, maxswpot

!!.......Local variables
      integer iel, drpdys, j, idf, lyr
      logical drpdlv
      real  accum(ISOS), ctodie, etodie, fr14, recres(MAXIEL), &
              tostore, srfclittr, soillittr,     &        
              rtdh, tempeff, watreff, cturn, eturn, temp
      data    drpdys /0/
      
      real :: dayhrs
      real :: sitlat
      integer :: month
      integer :: evntyp
      real, dimension (10):: sl_frootmc, sl_frootjc  !! dead fine root(mature and juvenile)that enter soil litter pool in each soil layer
      
      real :: ttrdis
      real, dimension (20):: prdis
!!    function declaration
      real daylenth   
!!      integer :: nelem
        
!!.......Saved variables
      save drpdlv               !drop dead leaves
      save drpdys               !days over which to drop leaves

      accum(LABELD) = 0.0
      accum(UNLABL) = 0.0

      
      j = 0
      idf = 0
      
      
      sl_frootmc  = 0.0
      sl_frootjc = 0.0

      
      j = ihru
      idf = idplt(j)
      month = i_mo

      dayhrs = daylenth (iida, ihru)
      sitlat = sub_lat(ihru)
      
      !if (curyr .gt. 1) drpdlv = .FALSE.  !! no death in the first year
      if (curyr == 1) drpdlv = .FALSE.  !! no death in the first year
!!.......Death of leaves
!!.......NOTE:  WOODDR(1)   - the death rate in fall for deciduous forests
!!.......       LEAFDR(MTH) - the monthly death rate for leaves in every
!!.......                     case except for fall in deciduous forests.
      if (leafc(j) .gt. 0.0001) then                !gC/m2
        if (DECIDf(idf) .ge. 1) then

!!....... Deciduous forest
!!....... If the daylight hours are increasing - it must be spring
         if ((thrsinc) .and. (tavewk .gt. TMPLFFf(idf))) drpdlv = .FALSE.
         
!!....... If daylight hours are decreasing and the temperature is low
!!....... enough drop leaves for fall season.
!!....... Add check for number of daylight hours to conditional for
!!....... determining if leaf drop should occur, cak - 06/30/03
!!....... If leaf drop has not occurred by the time the winter solstice
!!....... is reached force leaf drop to occur, cak - 10/28/04
            if (dofy ==1) drpdys = 0
        
          if (DECIDf(idf) .eq. 1) then              !Decidious forest
            if (( (tavewk .lt. TMPLFFf(idf)) .and. (.not. drpdlv) .and. &
                (.not. thrsinc) .and. (dayhrs .lt. 12.0)) .or.          &
               ((.not.drpdlv) .and. (sitlat .ge. 0) .and.               &
                (month .eq. 12)) .or.                                   & 
              ((.not.drpdlv) .and. (sitlat .lt. 0) .and.                &
                (month .eq. 6))) then
!!........... Allow leaf drop to occur over a 30 day period, dropping
!!........... all of the remaining leaves on the last day.
              if (drpdys .lt. 30) then
                ctodie = leafc(j) * WOODDRf(idf,LEAF)*tfrac   !!!! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX why leaf carbon pool is not updated here
                
                drpdys = drpdys + 1
                decidgrow = .FALSE.
               else
                ctodie = leafc(j) * WOODDRf(idf,LEAF)
                drpdlv = .TRUE.
                drpdys = 0
               !! evntyp = 2
              endif
            else
              ctodie = leafc(j) * LEAFDRf(idf,month)*tfrac
            endif
          elseif (DECIDf(idf) .eq. 2) then              !Drought decidious forest
!!......... Drought deciduous forest
!!......... Compute death for drought deciduous forests
            ctodie = leafc(j) * (1. - bgwfunc) * WOODDRf(idf,LEAF)*tfrac
          endif
        else
!!....... Continuous forest
!!....... Use leaf death rate multiplier from EACHYR
          ctodie = leafc(j) * LEAFDRf(idf,month)*tfrac * ldrmlt
        endif
       leafc(j) = leafc(j) - ctodie   
     endif  
!! Finished  calculating how much leaf died        
       
       
!!..... Compute E/C ratios
        !!nelem = 1
        do 10 iel = 1, nelem
          if (leafc(j)>0.01) then
              recres(iel) = eleaf(iel,j) / leafc(j)

    !!....... Compute flow to retranslocation storage
              tostore = recres(iel) * ctodie * FORRTFf(idf,iel)
              forstg(iel,j) = forstg(iel,j) + tostore
              eleaf(iel,j) = eleaf(iel,j) - tostore              

    !!....... Decrease E/C by the amount that is retranslocated
              recres(iel) = recres(iel) * (1 - FORRTFf(idf,iel))          
          else
          
              tostore = eleaf(iel,j)*FORRTFf(idf,iel)
              forstg(iel,j) = forstg(iel,j) + tostore
              eleaf(iel,j) = eleaf(iel,j) - tostore          
          endif
10      continue



!! currently we do not have such management activities, but may consider to add in case need to harvest leaf in the future. So keep the following codes here.
!!..... If evntyp is greater than 1 the leaves go to the source/sink
!!..... rather than to the litter, cak - 02/07/2006
      !! fr14 = rlvcis(LABELD) / rleavc
         !! evntype refers to the tree removal type
         !!         = 0 for cutting event
         !!         = 1 for fire event
         !!         = 2 leaves are removed from system
         
          !!IS THIS PART REPETITIVE, WHEN WE DO NOT COSNIDER EVNTYPE???
          do 15 iel = 1, nelem
              if (leafc(j)>0.01) then
                etodie = ctodie * (eleaf(iel,j) / leafc(j))
                eleaf(iel,j) = eleaf(iel,j) - etodie
                !!esrsnk(iel,j) = esrsnk(iel,j) + etodie
              else
                etodie = eleaf(iel,j)*(1-FORRTFf(idf,iel))
                esrsnk(iel,j) = esrsnk(iel,j) + etodie
                eleaf(iel,j) = 0.
              endif  
15        continue
        
!! following code should be used to particition leaf fall to surface structure and metbolic litter pools. Need to discuss with Xuesong
!! ctodie: dead carbon from leaf
!! recres: matix storing element:carbon ratio
!! SRFC: indicating where litter fall to, surface (SRFC =1) or belowground (SOIL = 2).

        call partit(ctodie, recres, SRFC, WDLIGf(idf,LEAF)) !! leaf fall has been substracted from leaf carbon and element pools in previous steps




!!.......Add code to age fine roots, juvenile fine roots age to the mature
!!.......fine root pool.  Modify this subroutine so that the death rate of
!!.......roots is a function of soil water potential and soil temperature,
!!.......cak - 06/28/2007
!!.......See:  A Model of Production and Turnover of Roots in Shortgrass Prairie
!!.......      Parton, Singh, and Coleman, 1978
!!.......      Journal of Applied Ecology
!!.......Cap the temperature effect on fine roots at -2 and +28 degrees C
      if (avgstemp .gt. 28.0) then
        temp = 28.0
      else if (avgstemp .lt. -2.0) then
        temp = -2.0
      else
        temp = avgstemp
      endif

!!.......Soil temperature effect on aging of juvenile roots
      if (frootjc(j) .gt. 0.0) then
        tempeff = gpdf(temp, 37.0, 0.0, 3.0, 3.0)
        cturn = TMXTURNf(idf) * tempeff * frootjc(j) * tfrac
        
        frootjc(j) = frootjc(j) - cturn   !! yound root become mature
        frootmc(j) = frootmc(j) + cturn
        

        do 60 iel = 1, nelem
          eturn = cturn * efrootj(iel,j) / frootjc(j)
          efrootj(iel,j) = efrootj(iel,j) - eturn
          efrootm(iel,j) = efrootm(iel,j) + eturn
          
60      continue
      endif

!!.......Soil temperature effect on root death rate
      tempeff = (temp - 10.0)**2 / 4.0 * 0.00175 + 0.1
      tempeff = min(tempeff, 0.5)
!!.......Soil water potential effect on root death rate
       
      
      watreff = 0.5 + (1.0 / DPI) * atan(DPI *0.05 * (mrspweff- 35))
      rtdh = max(tempeff, watreff)


      ctodie = 0.0  !! check if this variable need to be initialized
!!.......Death of juvenile fine roots
      if (frootjc(j) .gt. 0.0) then
        ctodie = frootjc(j) *WOODDRf(idf,FROOTJ) * rtdh * tfrac
        
        do 20 iel = 1, nelem
          recres(iel) = efrootj(iel,j) / frootjc(j)
20      continue
        !!fr14 = frtcisj(LABELD) / frootcj  !C14 label
!!..... A fraction of the dead roots are transferred to the surface
!!..... litter layer, the remainder goes to the soil litter layer
!!..... cak - 05/14/2007
        srfclittr = ctodie * WRDSRFCf(idf)
        soillittr = ctodie - srfclittr
        frootjc(j) = frootjc(j) - ctodie
        
!! distribute total soil litter from juvenile fine root to each soil layer based on the root fraction
     
          do lyr = 2, sol_nly(j)
         
             sl_frootjc(lyr) = rdis(lyr,j) * soillittr
             
     !! need to call the partition function here. check with xuesong later
     !!      sl_frootjc(lyr):: litter from juvenil fine root distributed to each layer
     
             call partit(sl_frootjc(lyr), recres, lyr, WDLIGf(idf,FROOTJ))        
             
          end do   
        
        
    !! check nutrient out    
        do 21 iel = 1, nelem
          etodie = ctodie * (efrootj (iel,j)  / frootjc(j))
          efrootj (iel,j) = efrootj (iel,j) - etodie 
21      continue
      endif

!!.......Death of mature fine roots
      if (frootmc(j) .gt. 0.0) then
        ctodie = frootmc(j) * WOODDRf(idf,FROOTM) * rtdh * tfrac
        do 25 iel = 1, nelem
          recres(iel) = efrootm(iel,j) / frootmc(j)
25      continue
      !!  fr14 = frtcism(LABELD) / frootcm
!!..... A fraction of the dead roots are transferred to the surface
!!..... litter layer, the remainder goes to the soil litter layer
!!..... cak - 05/14/2007
        srfclittr = ctodie * WRDSRFCf(idf)
        soillittr = ctodie - srfclittr
        frootmc(j) = frootmc(j) - ctodie


 !! distribute total soil litter from mature fine root to each soil layer based on the root fraction         
          
              do lyr = 2, sol_nly(j)
             
                 sl_frootmc(lyr) = rdis(lyr,j) * soillittr
         
 !! need to call the partition function here. check with xuesong later
 !! sl_frootmc(lyr): litter from mature fine root distributeed to each soil layer
 
                 call partit(sl_frootmc(lyr), recres, lyr,  WDLIGf(idf,FROOTJ))        
                 
              end do   

 !! check nutrient out       
        do 26  iel = 1, nelem
          etodie = ctodie * (recres(iel))
          efrootm (iel,j) = efrootm (iel,j) - etodie 
26      continue

!! XXXXX         call partit(srfclittr, recres, SRFC, frtcism, frootem,    !! not account for yet, need to discuss with xuesong
!! XXXXX      &              wdlig(FROOTM), fr14)
!! XXXXX         call partit(soillittr, recres, SOIL, frtcism, frootem,    !! not account for yet, need to discuss with xuesong
!! XXXXX      &              wdlig(FROOTM), fr14)
      endif

!!.......Fine Branches, Large Wood, and Coarse Roots go to the dead wood
!!.......compartments: WOOD1, WOOD2, WOOD3

!!.......Death of fine branches
      if (brchc(j).gt. 0.0) then
        ctodie = brchc(j) * WOODDRf(idf,FBRCH)*tfrac
!!..... remove labled fraction instead of the growth ratio, cisotf.  
!!..... KLK 7/2007
!        call csched(ctodie, cisotf, 1.0,
         brchc(j) = brchc(j) - ctodie
         woodc1(j) = woodc1(j) + ctodie
!! XXXXX         call csched(ctodie, fbrcis(LABELD), fbrchc,    !! not accounted for yet
!! XXXXX      &              fbrcis(UNLABL), wd1cis(UNLABL),
!! XXXXX      &              fbrcis(LABELD), wd1cis(LABELD),
!! XXXXX      &              1.0, accum)

        do 30 iel = 1, nelem
          etodie = ctodie * (ebrch(iel,j) / brchc(j))
          ebrch(iel,j) = ebrch(iel,j) - etodie
          woode1(iel,j) = woode1(iel,j) + etodie
          !! XXXXX           call flow(fbrche(iel), wood1e(iel), time, etodie) !!  not accounted for yet
30      continue

      endif

!!.......Death of large wood
      if (largwc(j).gt. 0.0) then
        ctodie = largwc(j) * WOODDRf(idf,LWOOD)*tfrac
!!..... remove labled fraction instead of the growth ratio, cisotf.  
!!..... KLK 7/2007
!        call csched(ctodie, cisotf, 1.0,
         largwc(j) = largwc(j) - ctodie
         woodc2(j) = woodc2(j) + ctodie
!! XXXXX         call csched(ctodie, rlwcis(LABELD), rlwodc,   !!  not accounted for yet
!! XXXXX      &              rlwcis(UNLABL), wd2cis(UNLABL),
!! XXXXX      &              rlwcis(LABELD), wd2cis(LABELD),
!! XXXXX      &              1.0, accum)

        do 40 iel = 1, nelem
          etodie = ctodie * (elargw(iel,j) / largwc(j))
          elargw(iel,j) = elargw(iel,j) - etodie
          woode2(iel,j) = woode2(iel,j) + etodie
 !! XXXXX          call flow(rlwode(iel), wood2e(iel), time, etodie) !!  not accounted for yet
40      continue
      endif

        
        !calculate root distribution
        ttrdis = 0.
        idf = idplt(j)
 !!       prdis = 0.
 !!        rdis = 0.
 !!      do lyr = 2, sol_nly(j)
  !!      
  !!       if (lyr<sol_nly(j)) then 
  !!          !prdis: percentage root distribution for a specific layer
  !!          prdis(lyr) = 0.5*(exp(-raf(idf)*sol_z(lyr-1,j)/1000.)+exp(-rbf(idf)*sol_z(lyr-1,j)/1000.)-exp(-raf(idf)*sol_z(lyr,j)/1000.) &
 !!                -exp(-rbf(idf)*sol_z(lyr,j)/1000.))     !raf and rbf read in
 !!        else
  !!          prdis(lyr) = 0.5*(exp(-raf(idf)*sol_z(lyr-1,j)/1000.)+exp(-rbf(idf)* sol_z(lyr-1,j)/1000.))       
  !!       end if        
  !!       ttrdis = ttrdis + prdis(lyr)
  !!     end do
       
 !!      do lyr = 2, sol_nly(j)
  !!          rdis(lyr,j) = prdis(lyr) / ttrdis 
  !!      end do


!!.......Death of coarse roots
      if (csrootc(j).gt. 0.0) then
        ctodie = csrootc(j) * WOODDRf(idf,CROOT)*tfrac
!!..... remove labled fraction instead of the growth ratio, cisotf.  
!!..... KLK 7/2007
!        call csched(ctodie, cisotf, 1.0,
          csrootc(j) = csrootc(j) - ctodie
          
            do lyr = 2, sol_nly(j)
                woodc3(lyr,j) = woodc3(lyr,j) + ctodie * rdis(lyr,j)
            end do 

!! XXXXX         call csched(ctodie, crtcis(LABELD), crootc, !! not account for yet
!! XXXXX      &              crtcis(UNLABL), wd3cis(UNLABL),
!! XXXXX      &              crtcis(LABELD), wd3cis(LABELD),
!! XXXXX      &              1.0, accum)

        do 50 iel = 1, nelem
          etodie = ctodie * (ecsroot(iel,j) / csrootc(j))
          
          ecsroot(iel,j) = ecsroot(iel,j) - etodie
          
          do lyr = 2, sol_nly(j)
            woode3(iel,lyr,j) = woode3(iel,lyr,j) + etodie * etodie* rdis(lyr,j)
          end do 
 !! XXXXX          call flow(croote(iel), wood3e(iel), time, etodie)  !! not account for yet !
50      continue
      endif
      
      
      
      leafn(j) =  eleaf(1,j)
      leafp(j) = eleaf(2,j)
      brchn(j)  = ebrch(1,j)
      brchp(j)  = ebrch(2,j)
      largwn(j)  = elargw(1,j)
      largwp(j)  = elargw(2,j)
      frootjn(j)  =  efrootj (1,j) 
      frootjp(j)  =  efrootj (2,j) 
      frootmn(j)  =  efrootm (1,j) 
      frootmp(j)  =  efrootm (2,j) 
      csrootn(j)  =  ecsroot(1,j)
      csrootp(j)  =  ecsroot(2,j)

      return
      end
