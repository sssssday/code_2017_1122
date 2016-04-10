      subroutine growfrst
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!   

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm
      implicit none 
      integer :: j,l
      integer:: idf !! forest id
      real :: delg, par, ruedecl, beadj, reg, f, ff, deltalai
      real :: laimax, rto
      real :: swavgfd



!! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!! photosysnthesis


!! calculate ppdf
!!~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    !!calculate daily ppfd and ppfd at noon
!!        !!Chen et al., 2005

!!       !!input:
!!        !!daylength  (second)
!!        !!zenith     (radians) solar zenith angle at noon
!!        !!swavgfd    (W/m2)  daylight average shortwave flux

!!        !!output:
!!
!!       !!ppfd_sunlai_daily         !!(umol/m2/s) daily average photosynthetic photo flux density absorbed by per sun leaf area
!!      !! ppfd_shadelai_daily       !!(umol/m2/s) daily average photosynthetic photo flux density absorbed by per shade leaf area
!!        !!ppfd_sunlai_noon          !!(umol/m2/s) photosynthetic photo flux density absorbed by per sun leaf area at noon
!!        !!ppfd_shadelai_noon        !!(umol/m2/s) photosynthetic photo flux density absorbed by per shade leaf area at noon
!!        !!swtrans                   !!(W/m2)  transmitted shortwave flux

           
        real :: sunhight, ash, cosash_diff
        real :: fdif,fdir,R
        real :: sdif,sdir, sdif_under, sscater
        real :: parw_ppfd
        real :: ppfd_sunlai_daily, ppfd_shadelai_daily
        real :: ppfd_sunlai_noon, ppfd_shadelai_noon
        real :: day_lai
        real :: dtor
        real :: tdtor
        real :: cossh, cosash, cos60
        real :: cosavz_dif
        real :: fsunppdf,fshadeppdf, ftmin,fvpdd,fco2,ft
        real :: sungs, shadegs
        real :: aco2, bco2
        real :: soildep !! change later, SWAT has multiple soil layers
        !!real :: prdis, tprdis!! root distribution of each layer
       !! real :: rdis !! summed root distribution of all layers
       ! real :: soilm_potsat, soilm_pot !! soil matix potential
        !!real :: sol_w
        real :: tmoist !! soil water impacts on stomatal conductance and photosynthesis
        !!real :: sol_wet
       !! real :: sol_sita_sat
       !! real :: sol_sita
        !!real :: fclay
        real :: photos
        real :: photo_sun
        real :: photo_shade
        real :: air_pre, pr_v1, pr_v2
        real :: po2
        real :: prfr, prfp,prfg,prft, prma, prr !! air pressure related variables
        real :: sunleafn_cn, shadeleafn_cn
        real :: biomass, sunl_cn, shal_cn  !!sun and shade leaf c/n ratio
        real :: lai, sunlai, shalai
        real :: fnsun, fnsha
        !! real :: m_respire
        real :: phe_stg !! phenology stage
        integer :: get_stage
        real :: cal_phen
        !!real :: alloc
        real :: get_nitrogen
        real :: day_stage, day_phen
         real :: dmaxleafc
        real :: ttemp
        real :: tleafc, tbrchc, tfrootjc,tfrootmc,tcsrootc, tseedc, tinic
        real :: tlargwc
        real :: tdumy
        real :: falloc
        real :: soil_sand, soil_clay
        real :: fdelg
        real :: atemp, stemp
        integer :: dumsun, dumsha
         real :: n11, n22, n33, n44, n55, n66, n77, n88, n99, n1010
        real, dimension (20):: rdis, prdis, solm_potsat, solm_pot, sol_sita_sat, sol_wet, fclay, sol_w
        real :: sshade, ssun
        real :: leaf_cn
        real :: maxleafc
        integer :: count
        integer :: iptr !! variable used to determine whether tree is juvenile or mature
        real :: mrspteff !! temperature effect on maintenance respiration
        real :: mrspweff !! soil effect on maintenance respiration
        real, dimension (10):: dsoilthickness !! soil thickness of 10 layers defined by daycent. change unit from cm to mm
        real :: tlaypgdth  !! soil depth of the TLAYPG layers (lays nutrient affect plant growth)
        integer ilyer !! count soil layers defined by daycent
        integer tlaypgswat !!  lay numbers in swat corresponds to tlaypg in dahycent
        real :: tfrac !! fraction of one day in a given month
        real :: mrspReduce !! factors to reduce maintenance respiration based on the carbon storge pool
        real :: rleavc_opt !! optimum leaf carbon
        real :: mrleaf, mrbrch, mrlgwood, mrfrootj, mrfrootm, mrcroot !! maintenance respiration
       
        real, dimension (6):: tree_cfrac !! fraction of different plant tissuse for carbon allocation
        real :: cdlength, pdlength !! day length of previous day and current day
        integer :: cday, pday !!current day and previous day
        integer :: checkgrow   !! check growing state of trees
        real :: rimpct  !! root impacts on nutrient avaibility
        integer :: forgrow   !! indicate forest start to grow ! curretly do not consider this because forest start to growth when photosynthesis is above zero
        real :: cprodfLeft !! left gpp after allocation to fine root
        real :: remCfrac  !! allocation fraction left after allocating C to fine root and leaf
        real :: totCup !! temperory variable to sum up allocation fractions
        real :: sum_cfrac !! temperory variable to sum up allocation fractions
        integer :: ipart  !! plant tissue numbers
        real :: grleaf, grbrch, grlgwood, grfrootj, grfrootm, grcroot !! growth respiration
        real :: cprodfdy, eprodfdy(3) !! carbon and nutrient production (increase) in plant
        real :: avn,avNH4,avNO3
        real, dimension (1):: fgpp  !! variable used to transfer actural gpp in the restrp process. 
        real,dimension(2) :: favail
        integer :: eli
        real, dimension(2,6,2)::ccefor 
       real ::riint
       real :: rictrl
       real :: treeNfix
       real :: relyld !! fraction of actural and potential production
       real :: elimit !! id of element that is limiting growth
       real, dimension(4,MAXIEL):: uptake
       real, dimension (FPARTS) :: euf
       real :: amt, namt !! nutrient uptake temperory variable
       real :: fsol
       real :: fno3, fnh3, tno3nh3
       real :: tavewk !! average temperature for the last seven days
       real :: lfncon
       real :: lfncmin
       real :: lfncmax
       real :: ldrmlt
       real :: tavgp7d !! average temperature of the past 7 days
       integer :: dd, cdd !! day counter
       integer :: yeardays 
       real :: avgstemp
       real :: bgwfunc
       integer :: ilyr
       real,dimension (10) :: rel_wc
       real :: av_rel_wc
       real,dimension (10) :: sol_swclimit 
       integer :: iel, lyr
       
       real,dimension (6) :: mfprd
       real :: toler
       real :: calcup
       
 !!    declare function type
       real maxswpot, line, daylenth, rtimp, leafa
       integer grochk
       integer tdays !! days in one month
      !!data growday /0/
      !!save growday
      
 !!    initializing local variables    
     
         rdis = 0.  
         prdis = 0.  
         solm_potsat = 0.  
         solm_pot = 0.   
         sol_sita_sat = 0.   
         sol_wet = 0.   
         fclay = 0.   
         sol_w = 0.  
         count = 0
         
         mrleaf = 0. 
         mrbrch = 0.
         mrlgwood = 0.
         mrfrootj = 0. 
         mrfrootm = 0.
         mrcroot = 0.
         
         ipart = 0
         
         grleaf = 0. 
         grbrch = 0.
         grlgwood = 0.
         grfrootj = 0. 
         grfrootm = 0.
         grcroot = 0.
         
         avn = 0.
         avNH4 = 0.  
         avNO3 = 0.
         fgpp = 0.
         favail = 0.
         treeNfix = 0.
         elimit = 0.
         uptake = 0.
         euf = 0.
         fsol = 0.
         
         lfncon = 0.
         lfncmin = 0.
         lfncmax = 0.
         ldrmlt = 1.
         tavgp7d = 0.
         
         eli = 0
         ilyr = 0
         rel_wc = 0.
         av_rel_wc = 0.
         sol_swclimit = 0.
         cprodfdy = 0.
         
     !! calculate how many days are there in a month    
         
        select case(mo_chk)
        case (9, 4, 6, 11)
          tdays = 30
        case (2)
          tdays = 29 - leapyr
        case default
          tdays = 31
        end select
         
         tfrac = 1.0/tdays
         tree_cfrac = 0. !! check if all elements are initialized as 0
         iel = 0
         toler = 0.000001
        
        cal_temp(3) = tfrac
        
        
         j = ihru
        
         rto = 1.
         idf = idplt(j)
        dtor = 3.1415926/180.     !!degree to radians

        tdtor = dtor*(-23.4)*cos(dtor*360.0*(dofy+10.0)/365.0)

        cossh = sin(tdtor)*sin(dtor*sub_lat(j)) +cos(tdtor)*cos(dtor*sub_lat(j))
        cossh = min(cossh, 0.99999)
        sunhight = acos(cossh)
       
        
        print*,dofy

   if ( dofy ==1)then 
      !! processes only occure once every year   
       ccefor = 0.
       growday = 0 
           call co2eff (iyr)
           
!! ... Added effect of co2 for forest; done here because not calcualted
!! ... dynamically based on biomass like grassland/crop
!! ... Direct CO2 effects only C/E ratio of leaves.
         !!nelem = 1
      do 30 iel = 1, nelem
        ccefor(IMIN,LEAF,iel) = CERFORf(idf,IMIN,LEAF,iel) * co2cce(2,IMIN,iel)
                           
        ccefor(IMAX,LEAF,iel) = CERFORf(idf,IMAX,LEAF,iel) * co2cce(2,IMAX,iel)
                             
30    continue

      do 50 ipart = 2, FPARTS-1
        do 40 iel = 1, nelem 
          ccefor(IMIN,ipart,iel) = CERFORf(idf,IMIN,ipart,iel)
          ccefor(IMAX,ipart,iel) = CERFORf(idf,IMAX,ipart,iel)
40      continue 
50    continue

 !!! variable related to death rate muultiplier
 
 !! ... Change leaf death rate multiplier if you have floating C/E ratios.
      if (ccefor(IMIN,LEAF,N) .ne. ccefor(IMAX,LEAF,N)) then
        if (leafc(j) .gt. 0) then
          lfncon = eleaf(j,N) / leafc(j)
          lfncmin = 1 / ccefor(IMIN,LEAF,N)
          lfncmax = 1 / ccefor(IMAX,LEAF,N)
          ldrmlt = 1 + (MAXLDRf(idf) - 1) *                &
                (lfncon - lfncmin) / (lfncmax - lfncmin)
        endif
      endif
 
        

   endif    
 !!         cosash = Cos(avz)
!!        float cos60  = cos(1.04)
!!        float cosash_dif = 0.537 + 0.025*flai
!!        float fdif,fdir         !! the fraction of diffuse and direct solar radiation

           
!!        sunhight = 2.0       !! assume height is 2.0, will change later
        ash = (sunhight + 3.1415926/2.0 )/2.0      !! avz in DLEM!!average zenith during the day
        cosash = cos(ash)
      !!   R = (hru_ra(j)*1.0e6*86400) / (cosash * 1367) !! convert hru_ra MJ/m2/day to W/s
        cos60 = cos(dtor*60.0)

!!      in DLEM, fdif was first calculated with R, but then assigned with this values. need to check later
        fdif = 0.3
       
        fdif = max(0.,min(1.,fdif))
        fdir = 1.0 - fdif
        swavgfd = hru_ra(j)*1.0e6/86400.0
        sdif = swavgfd * fdif    !!!! (W/m2) diffuse solar radiation
        sdir = swavgfd * fdir   !! !! (W/m2) direct solar radiation
!!        float sdif_under                !! (W/m2) diffuse solar radiation under canopy
!1        float sscater                   !! (W/m2) scattered solar radiation by sun leaf
       
        lai = leafc(j) * avg_slaf(idf) !!temporarill closed and changed to the maxlai
         !!if (lai>11.0)lai = 11.
       !!lai = MAXLAIf(idf) 
        if (lai < 0.0001) then
            lai = 0.0001
            count = count+1
            
            leafc(j) = lai / avg_slaf(idf)
            leafn(j) = leafc (j)/ fcn_leafminf(idf)
        end if

   
        cosavz_dif = 0.537 + 0.025*lai;
        sdif_under = sdif * exp(-0.5 * omgaf(idf) * lai/cosavz_dif) !! modify omgaf later
        sscater = 0.07*omgaf(idf)*sdir*(1.1-0.1*lai)*exp(-cosash)

        if(lai > 0) then
            sshade = (sdif-sdif_under)/lai + sscater
        else
            sshade = 0.
        end if
        ssun   = sshade + sdir*cos60/cosash

        
        parw_ppfd = 4.55         !! !!convert from PAR W/m2 -> umol/m2/s PPFD
                                       !! !!Ruimy et al., 1995
        
        if(leafc(j)>0.)then
                ppfd_sunlai_daily = 0.45*ssun*parw_ppfd
                ppfd_shadelai_daily = 0.45*sshade*parw_ppfd

        else
                ppfd_sunlai_daily = 0.
                ppfd_shadelai_daily = 0.
        endif


        ppfd_sunlai_noon = 3.14/2. * ppfd_sunlai_daily
        ppfd_shadelai_noon = 3.14/2. * ppfd_shadelai_daily
        
!!        float tempvar2 = exp(-pam.ext_coef*lai)

!!        swabs = swavgfd * (1.-pam.sw_alb) * (1.-tempvar2)     !!absorbed solar radiation  over biome area

!!        swtrans = swavgfd * exp(-0.5*pam.omega * lai/cosavz_dif )

!!
!! ~ ~ ~ ~ ~  ~ ~ ~ ~ ~  ~ ~ ~ ~ ~  ~ ~ ~ ~ ~  ~ ~ ~ ~ ~  ~ ~ ~ ~ ~  ~ ~ ~ ~ ~  ~ ~ ~ ~ ~  ~ ~ ~ ~ ~  ~ ~ ~ ~ ~  ~ ~ ~ ~ ~

!!
!! calculate stomatal conductance

!!      1, calculate soil wetness, fclay
!!      2, calculate matrix potential
 !!     3, calcualte root distribution
!!      4, calcualate fmoist and other f factors
!! ####### for soil properties, do not consider different lays yet
!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!      !!Stomatal conductance
!!        !!Chen J. et al., 2005
!!        !!Running & Coughlan, 1988
!!        !!Jarvis & Morison, 1981
!!        !!Chen F. et. al., 1996

!!        !!input:
!!        !!ITEMpsnparameters psnparam
!!        !!float tmin              !!foliage temperature (degree celsius) daily min temperature
!!        !!float vpd             !!(pa) vapor pressure deficit
!!        !!float btrans          !!coefficient for stomata conductance (0-1)
!!        !!float ppdf            !!(umol/m2/s) photosynthetic photo flux density

!!        !!output: stomatal conductance m/s



        soil_sand = 0.
        soil_clay = 0.
        soildep = 0.

!!          calculate wetness
        do l = 1, sol_nly(j)
            soil_sand = soil_sand + sol_sand(l,j)
            soil_clay = soil_clay + sol_clay(l,j)
            soildep   = soildep + sol_z(l,j)
            solm_potsat(l)= -10. * (10**(1.88-0.0131*sol_sand(l,j)))
            sol_sita_sat(l) = 0.489-0.00126*sol_sand(l,j)
            sol_wet(l) = sol_st(l,j)/(sol_z(l,j))/sol_sita_sat(l)
        if (sol_wet(l)<0.01) sol_wet(l)=0.01

            fclay(l)= 2.91 + 0.159 * sol_clay(l,j)
        
        end do

        !!soil_sand = soil_sand / sol_nly(j)
        !!soil_clay = soil_clay / sol_nly(j)

       
      
       !!sol_wet = sol_sita / sol_sita_sat  !! here ice in soil is not considered

   
!!         calculate fcaly
    
!!        calculate matrix potential
       
       do l = 1, sol_nly(j)
            if (sol_wet(l)>1.)sol_wet(l)=1.
       
            solm_pot(l) = solm_potsat(l) * (sol_wet(l)**(-fclay(l)))

            if(solm_pot(l) < -1.e8) solm_pot(l) = -1.e8
        
            sol_w (l) = (-150000. - solm_pot(l))/(-150000. - solm_potsat(l))
            if(sol_w(l) <= 0.) sol_w(l) = 0.
            if(sol_w(l) > 1.) sol_w(l) = 1.
        
       end do
      
!!     calcualte fmoist


        


!!     calculate root distribution
      tmoist = 0.

       do l = 2, sol_nly(j)
       
        if (l<sol_nly(j)) then 
           prdis(l) = 0.5*(exp(-raf(idf)*sol_z(l-1,j)/1000.)+exp(-rbf(idf)*sol_z(l-1,j)/1000.)-exp(-raf(idf)*sol_z(l,j)/1000.) &
                -exp(-rbf(idf)*sol_z(l,j)/1000.))
        else
           prdis(l) = 0.5*(exp(-raf(idf)*sol_z(l-1,j)/1000.)+exp(-rbf(idf)* sol_z(l-1,j)/1000.))
       
        end if 
       
            rdis(l) = prdis(l)
       end do
       
       do l = 2, sol_nly(j)
        
            tmoist = tmoist + rdis(l)* sol_w(l)
        
       end do
       
       if (tmoist < 1.e-10)tmoist = 1.e-10
     
          

         f_moist(j) = tmoist

!!
!!        !!fppdf = ppdf*psnparam.ppdfcoef/(1.0+ppdf*psnparam.ppdfcoef)

!!        !!Biome-BGC method
!!      example: g_noon_sun = biome[ipft].Stomatal_conductance(ppsn[ecotype],daytemp[day],dtmin[day],dPAir[day],dvpd[day],btrans,biome[ipft].p_leaf.ppfd_sunlai_noon,dco2[day],1)
        fsunppdf = ppfd_sunlai_noon/(75.+ppfd_sunlai_noon)
        fshadeppdf = ppfd_shadelai_noon/(75.+ppfd_shadelai_noon)
!!        !!if(tc<1.0) ft = 0.0

!!        !!else if(tc<psnparam.topt) ft = log(tc)/log(psnparam.topt)
!!        !!else if(tc<psnparam.trange) ft = cos(3.14159/2.0*(tc-psnparam.topt)/(psnparam.trange-psnparam.topt))
!!        !!else ft = 0.0

!!        !! NOAH-MP's method
        ft = 1.0- 0.0016* ( (toptdf(idf) - tmpav(j))*(toptdf(idf) - tmpav(j)))
        ft = max(ft,0.0001)

!!        /* freezing night minimum temperature multiplier */

       
        if(vpd < vpd_openf(idf)) then

            fvpdd = 1

        else if (vpd < vpd_closef(idf)) then

            fvpdd = (vpd_closef(idf)-vpd)/(vpd_closef(idf)-vpd_openf(idf))

        else

            fvpdd =0.0

        end if

 !!       !!CO2 effect on conductance is derived from Ainsworth & Long (2004)
!!        !!Base is set as 350 ppm, elevated CO2 assumed as mean 550 (+200PPM), the average response of 0.8 is used
        aco2 = -0.001
        bco2 = 1.35
        fco2 = aco2*co2(hru_sub(j)) + bco2

        sungs = max(gmaxf(idf)*fsunppdf*fvpdd*tmoist*fco2*ft,gmaxf(idf))
       
        shadegs = max(gmaxf(idf)*fshadeppdf*fvpdd*tmoist*fco2*ft,gmaxf(idf))
!!      calculate air pressure

       prfr = 0.0065        !! (m3 Pa mol-1 K-1) gas law constant */
       prfg = 9.80665         !! (m s-2) standard gravitational accel. */
       prfp= 101325.0         !! (Pa) standard pressure at 0.0 m elevation */
       prft = 288.15
       prr = 8.3143
       prma = 2.89644e-2


        pr_v1 = 1.0 - (prfr * sub_elev(j))/prft
        pr_v2 = prfg / (prfr * (prr / prma))

        air_pre = prfp * (pr_v1) ** pr_v2
        po2 = 0.209 * air_pre
        
!!     calculate N limitation

      
        sunlai =2.*cos(sunhight)*(1.-exp(-0.5*omgaf(idf)* lai/cos(sunhight)))
        shalai = MAXLAIf(idf) - sunlai
        
        
        

        biomass = leafc(j)/0.45  !! qichun modify later
        dumsun = 1
        dumsha = 2
        sunl_cn = leaf_cn (biomass, leafn(j), lai, sunlai, j, dumsun)
        shal_cn = leaf_cn (biomass, leafn(j), lai, sunlai, j, dumsha)

     !!   n33 = leafc(j) +brchc(j) + largwc(j) + lgrootc(j) + frootjc(j)+ frootmc(j) + seedc(j) 
        
         if(sunl_cn < fcn_leafminf (idf)) then
            fnsun = 1.
         else
            fnsun = exp(-0.08*(sunl_cn - fcn_leafminf(idf)))
         end if
         
          !!if (fnsun<0.5)fnsun=0.5  !! change later
         if(shal_cn < fcn_leafminf(idf)) then
            fnsha = 1.
         else
            fnsha = exp(-0.08*(shal_cn - fcn_leafminf(idf)))
         endif
         fnsun = 1.
         fnsha = 1.
          !!if (fnsha<0.5)fnsha=0.5
         !! print*, fnsun

!!        photosynthesis
       tdumy = 0.0
      

        photo_sun = photos(sungs, tdumy,tmpav(j),tmoist, po2, co2(hru_sub(j)), fnsun, ppfd_sunlai_daily, air_pre, j)
     
       
        photo_shade = photos(shadegs, tdumy,tmpav(j),tmoist, po2, co2(hru_sub(j)), fnsha, ppfd_shadelai_daily, air_pre, j)
  !!
      !! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX should calculate f_pgpp here
     
        f_pgpp(j) =  (sunlai*photo_sun + shalai*photo_shade) * 12.011e-6 * dayl(j) * 3600.  !! umol/m2 leaf/s -> gC/day, dayl: day lenghten in hour
       
        !write(*,*) leafn(j), brchn(j), largwn(j), lgrootn(j), frootjn(j),frootmn(j), seedn(j)
        





               
        print*, 'f_pgpp(j)=', f_pgpp(j)
        if (f_pgpp(j)<0.0001) go to 1099
       
       !!tleafc = leafc(j)
       !!tbrchc = brchc(j)
       !!tlargwc = largwc(j)
       !!tfrootjc = frootjc(j)
       !!tfrootmc = frootmc(j)
       !1tcsrootc = csrootc(j)
       !!tseedc = seedc(j)
     !  tinic = inicdf(idf) !! 
      !! atemp = tmpav(j)
    !!   stemp = sol_tmp(2,j) !! temperature use T of the second layer. change later
       
       !! call m_respire(atemp, stemp, tleafc, tbrchc, tlargwc,tfrootjc,tfrootmc,tcsrootc, tseedc, tinic, j)   !! here does not consider respiration from different layers, change later


  
  !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !! following processes are from Daycent 
    
    !! ... Determine nutrients available to plants for growth.
  !! add one more dimension to most variables defined in daycent, need to check how will this affect these variables when they are used as arguments of functions
    
  !! determine wether forest is yong or old~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
!! ...... Determine old or new forest
!! ...... iptr points to new forest carbon allocation fractions (iptr = 1) or
!! ...... mature forest carbon allocation fractions (iptr = 2) for each of the
!! ...... tree parts; leaves, fine roots, fine branches, large wood, and
!! ...... coarse roots.  Switch from new forest allocation fractions to old
!! ...... forest allocation fractions at time = swold
      if (curyr .le. SWOLDf(idf)) then
!! ........ Use juvenile forest C allocation fractions
        iptr = 1
      else
!! ........ Use mature forest C allocation fractions
        iptr = 2
      endif
  
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !! temperature effect on maintenance respiration
    
    mrspteff= 0.1 * exp(0.07 * tmpav(j))
    
    if (mrspteff > 1.0)mrspteff = 1.0
    if (mrspteff < 0.0)mrspteff = 0.0
    
   !! soil moisture effect on maintenance respiration 
   
        dsoilthickness (1:10) = (/100.0, 200.0, 150.0, 150.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0/) !! need to check if values are assigned correctly.
    
    
    tlaypgdth = 0.0
    
     
    

    do ilyer = 1, TLAYPGf(idf)   !! check if value should plus one
  
       tlaypgdth = tlaypgdth + dsoilthickness(ilyer)
        
    end do
    
    
    
    
    tlaypgswat = 0
    soildep = 0.
    
    do lyr = 1, sol_nly(j)

       soildep   = soildep + sol_z(lyr,j)
       
       if (soildep > tlaypgdth) exit  !! test this later
    end do
    
    
       tlaypgswat = lyr !! get the swat layer that will be considered in the following steps
       
       
        do lyr = 1, tlaypgswat
        avNO3 = avNO3 + sol_no3(lyr,j)/10.
        avNH4 = avNH4 + sol_nh3(lyr,j)/10.
        end do

        avn = avNH4 + avNO3;
        
        availm(j,1) = avn  
        
        
       
       mrspweff = maxswpot(tlaypgswat) 
                    
    if (mrspweff > 1.0)mrspteff = 1.0
    if (mrspweff < 0.0)mrspteff = 0.0
        
         !! set mimium value for the 
          if (carbostg(j) .lt. 15.0) then
          
              csrsnk(j) = csrsnk(j)- (15.0 - carbostg(j))
              carbostg(j) = 15.0
              
          end if    
              
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
    !! calculate the maintence reducing factor based on the size of carbon storage pool
  
       !! call lacalc(lai, fbrchc, rlwodc, maxlai, klai) we do not use this method since we calculate lai directly from leaf carbon
      rleavc_opt = lai / (2.5 * BTOLAIf(idf))
      if (rleavc_opt > 0.000000001) then
        if (carbostg(j) .lt. FMRSPLAIf(idf,3) * rleavc_opt) then
          mrspReduce = line(carbostg(j), FMRSPLAIf(idf,1) * rleavc_opt, FMRSPLAIf(idf,2), FMRSPLAIf(idf,3) * rleavc_opt, FMRSPLAIf(idf,4))
        elseif (carbostg(j) .gt. FMRSPLAIf(idf,5) * rleavc_opt) then
          mrspReduce = FMRSPLAIf(idf,6)
        else
          mrspReduce = line(carbostg(j), FMRSPLAIf(idf,3) * rleavc_opt, FMRSPLAIf(idf,4), FMRSPLAIf(idf,5) * rleavc_opt, FMRSPLAIf(idf,6))
        endif
      else
        mrspReduce = 0.0
      endif
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~       
  !! maintenance respiration from all tissues   
  
        !! ids for different plant tisses
         !!LEAF = 1
         !!FROOT = 2
         !!FROOTJ = 2
         !!FROOTM = 6
         !!FBRCH = 3
         !!LWOOD = 4
         !!CROOT = 5
  
         mrleaf = tfrac *leafc(j)*mrspReduce*mrspteff* FKMRSPMXf(idf,1)
         mrfrootj = tfrac *frootjc(j)*mrspReduce*mrspteff*mrspweff* FKMRSPMXf(idf,2) 
         mrfrootm = tfrac *frootmc(j)*mrspReduce*mrspteff*mrspweff* FKMRSPMXf(idf,6)
         mrbrch = tfrac *brchc(j)*mrspReduce*mrspteff* FKMRSPMXf(idf,3)
         mrlgwood = tfrac *largwc(j)*mrspReduce*mrspteff* FKMRSPMXf(idf,4)
         mrcroot = tfrac *csrootc(j)*mrspReduce*mrspteff* FKMRSPMXf(idf,5)
   !! check later. be very very careful that the tissues pools are not updated after respiration yet.  
   !! reduce this from the carbon storage pool and added to csrsnk pool
     
        f_mr(j) = mrleaf  + mrbrch + mrlgwood + mrfrootj + mrfrootm + mrcroot
        
        !!  leafc(j) =  leafc(j) - mrleaf
        !!  frootjc(j) = frootjc(j) - mrfrootj
        !!  frootmc(j) = frootmc(j) - mrfrootm
        !!  brchc(j) = brchc(j) - mrbrch
        !!  largwc(j) = largwc(j) - mrlgwood
        !!  csrootc(j) = csrootc(j) - mrcroot
          
          
                  
        if (carbostg(j)> f_mr(j)) then 
        carbostg(j) = carbostg(j) -  f_mr(j) 
        else 
        
        f_mr(j) = carbostg(j)
        
        carbostg(j) = carbostg(j) -  f_mr(j) 
        endif
        
        csrsnk(j) = csrsnk(j) + f_mr(j) 
        
       cal_temp(1)=  f_mr(j)
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

                        
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~         

  
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
    !! carbon allocation
    !! get the alloction fraction for each tissuse
       !! determine if day hours are increasing or decreasing
          cday = iida
          pday = cday -1
          
          if (pday<1)pday = 365   !! check how to account for leap year
          
          cdlength  = daylenth (cday, j)
          
          pdlength = daylenth (pday, j)
          
          
          
        if (cdlength .lt. pdlength) then
                 hrsinc = .FALSE.
        else if (cdlength .gt. pdlength) then
          hrsinc = .TRUE.
        endif
       !! get growing stage indicator
       
       checkgrow = grochk(hrsinc,pdlength,tmpav(j))
       
       
        tswatlyr(j) = tlaypgswat
        call alloc (f_pgpp(j), tmpav(j), tlaypgswat,tree_cfrac, checkgrow)
      riint = 0.8
      rictrl = 0.015

      rimpct = rtimp(riint, rictrl, frootjc(j)+frootmc(j))

    !!f_pgpp(j) = 100.
    if (decidgrow) then
    !! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXcheck should calculate f_gpp based on f_pgpp here   because in daycent actural gpp is calcualted twic, consider to further change this later
               call restrp(elimit, nelem, availm, ccefor, 2, tree_cfrac,        &
                           f_pgpp(j), rimpct, forstg, SNFXMX2f(idf), cprodfdy,     &
                           eprodfdy, uptake, tree_a2drat, treeNfix, relyld)
    endif   
    !! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
      print*, "cprodfdy0=",cprodfdy
       
       
     if (decidgrow.and.cprodfdy>0.) then
   
     
      !! count growth day
       growday = growday +1
     
     
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 !!... Compute carbon allocation fractions for each tree part
!!... Calculate how much of the carbon the roots use

 !! here I replace cprodfdy with gpp
        cprodfLeft = cprodfdy - (cprodfdy * tree_cfrac(FROOT))   !! how this fraction was calculated

!!... Calculate how much of the carbon the leaves use, we allocate leaves
!!... up to a optimal LAI
        tree_cfrac(LEAF) = leafa(lai,leafc(j), brchc(j), largwc(j), cprodfLeft,  &
                               cprodfdy)
        remCfrac = 1.0 - tree_cfrac(FROOT) - tree_cfrac(LEAF)

!!... If we have leftover carbon allocate it to the woody plant parts
!!... using a weighted average
        if (remCfrac .lt. 1.0E-05) then
!!... for FBRCH, LWOOD, and CROOT ...
          do 60 ipart = FBRCH, CROOT
            tree_cfrac(ipart) = 0.0
60        continue
        else
!!..... for FBRCH, LWOOD, and CROOT ...
          totCup = 0.0
          do 70 ipart = FBRCH, CROOT
            tree_cfrac(ipart) = FCFRACf(idf,ipart, iptr)
            totCup = totCup + tree_cfrac(ipart)
70        continue
          if (totCup .gt. 0.0) then
!!..... for FBRCH, LWOOD, and CROOT ...
            do 80 ipart = FBRCH, CROOT
              tree_cfrac(ipart) = tree_cfrac(ipart) / totCup * remCfrac
             
80          continue
          else
            write(*,*) 'Error in treegrow'
            write(*,*) 'fcfrac(FBRCH)+fcfrac(\e'
            STOP
          endif
        endif


!! ... Error checking
        sum_cfrac = 0.0
        do 90 ipart = 1, 5
          if (tree_cfrac(ipart) .lt. 0.0) then
            write(*,*) 'Error in treegrow, tree_cfrac(ipart) < 0'
            STOP
          else if (tree_cfrac(ipart) .gt. 1.0) then
            write(*,*) 'Error in treegrow, tree_cfrac(ipart) > 1'
            STOP
          else
            sum_cfrac = sum_cfrac + tree_cfrac(ipart)
          endif
90      continue
        if (abs(1.0 - sum_cfrac) .gt. 0.001) then
          write(*,*) "Error in tree carbon allocation fractions!"
          write(*,*) "sum_cfrac = ", sum_cfrac
!!          STOP
        endif
        
        print*, "cropdfy1 = ", cprodfdy
        
        cprodfdy = 0.
    !! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     !!f_pgpp(j) = 100.
     !! limit potential gpp again here
      call restrp(elimit, nelem, availm, ccefor, FPARTS-1,                &
                  tree_cfrac, f_pgpp(j), rimpct, forstg, SNFXMX2f(idf), & 
                  cprodfdy, eprodfdy, uptake, tree_a2drat, treeNfix,    &
                  relyld)

    
    !! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
      print*, "cropdfy2 = ", cprodfdy
        
         !! fixn = treeNfix    
          
           f_gpp(j)= cprodfdy
           carbostg(j) = carbostg(j) + f_gpp(j)
             csrsnk(j) = csrsnk(j) - f_gpp(j)
        
!! ..... Calculate production for each tree part
!! ..... New variable MFPRD added for gridded output - 6/96 rm
        do 95 ipart = 1, 6
          if (ipart .eq. LEAF .or. ipart .eq. FBRCH .or.      &
             ipart .eq. LWOOD .or. ipart .eq. CROOT) then
            mfprd(ipart) = tree_cfrac(ipart) * f_gpp(j)    !! change later
          else if (ipart .eq. FROOTJ) then
            mfprd(ipart) = tree_cfrac(FROOTJ) * f_gpp(j) *    & 
                         (1.0 - WMRTFRACf(idf))
          else if (ipart .eq. FROOTM) then
            mfprd(ipart) = tree_cfrac(FROOTJ) * f_gpp(j) * WMRTFRACf(idf)
          else
            write(*,*) 'Error in treegrow, ipart out of bounds = ',  ipart          
            STOP
          endif
95      continue    


!!     calculate growth respiration
         grleaf = 0. 
         grbrch = 0.
         grlgwood = 0.
         grfrootj = 0. 
         grfrootm = 0.
         grcroot = 0.
        
           
        grleaf = mfprd(LEAF) * FGRESPf(idf,LEAF)
        grfrootj = mfprd(FROOTJ) * FGRESPf(idf,FROOTJ)
        grfrootm = mfprd(FROOTM) * FGRESPf(idf,FROOTM)
  
        grbrch = mfprd(FBRCH) * FGRESPf(idf,FBRCH)
        grlgwood = mfprd(LWOOD) * FGRESPf(idf,LWOOD)
        grcroot   = mfprd(CROOT) * FGRESPf(idf,CROOT)

        
 !!      add carbon to tissuse pools a
         leafc(j)= leafc(j)+  mfprd(LEAF)
         frootjc(j)=  frootjc(j)+ mfprd(FROOTJ)
         frootmc(j)= frootmc(j)+ mfprd(FROOTM)
         if ((DECIDf(idf).eq. 1) .and. (growday .gt. FURGDYSf(idf))) then
         mfprd(FBRCH)= mfprd(FBRCH)* (1-FLSGRESf(idf))
         mfprd(LWOOD)= mfprd(LWOOD)* (1-FLSGRESf(idf))
         mfprd(CROOT)= mfprd(CROOT)* (1-FLSGRESf(idf))
         
         grbrch = mfprd(FBRCH) * FGRESPf(idf,FBRCH)
         grlgwood = mfprd(LWOOD) * FGRESPf(idf,LWOOD)
         grcroot   = mfprd(CROOT) * FGRESPf(idf,CROOT)
         endif
         
         brchc(j)= brchc(j)+ mfprd(FBRCH)
         largwc(j)= largwc(j)+ mfprd(LWOOD)
         csrootc(j)= csrootc(j)+ mfprd(CROOT)
         
         
        f_gr(j) = grleaf + grfrootj + grfrootm + grbrch + grlgwood + grcroot
        

        carbostg(j) = carbostg(j) - mfprd(LEAF) - mfprd(FROOTJ) - mfprd(FROOTM) - mfprd(FBRCH)- mfprd(LWOOD)- mfprd(CROOT)
        
        if (carbostg(j)> f_gr(j)) then 
        carbostg(j) = carbostg(j) -  f_gr(j) 
        else 
        
        f_gr(j) = carbostg(j)
        
        carbostg(j) = carbostg(j) -  f_gr(j) 
        endif
        
        csrsnk(j) = csrsnk(j) + f_gr(j) 
        f_npp(j) = f_gpp(j)- f_gr(j)-f_mr(j)





     
      
        !! update biomass carbon pools
        
         
        
   
        !! growht respiration
        cal_temp(2) = carbostg(j)
        
         
        
       
        
!!   n inputs to plant biomass        
!!    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  

     !! calculate uptake ratio of each tissue
     
     
    !! ..... Actual Uptake
        do 110 iel = 1, nelem
          if (eprodfdy(iel) .eq. 0.0) then
            write(*,*) 'Divide by zero in treegrow, eprodfdy(iel) = 0'
            STOP
          endif
          euf(LEAF) = eup(j,LEAF,iel) / eprodfdy(iel)
          euf(FROOTJ) = (eup(j,FROOT,iel) * (1.0 - WMRTFRACf(idf))) /   &
                       eprodfdy(iel)
          euf(FROOTM) = (eup(j,FROOT,iel) * WMRTFRACf(idf)) / eprodfdy(iel)
          euf(FBRCH) = eup(j,FBRCH,iel) / eprodfdy(iel)
          euf(LWOOD) = eup(j,LWOOD,iel) / eprodfdy(iel)
          euf(CROOT) = eup(j,CROOT,iel) / eprodfdy(iel)
!! ....... Reset eprodcdy(iel) to track the actual uptake which can be
!! ....... restricted late in the growing season
          eprodfdy(iel) = 0.0
      !!  up  take from the storg pool 
      !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
 !!         !!.......... Take up nutrients from internal storage pool
!! ....... Don't allow uptake from storage if forstg is negative -mdh 8/8/00
!! ....... If the tree is decidious and the time alloted to grow the
!! ....... woody components of the trees has passed use the late season
!! ....... growth restriction parameter value to determine how much
!! ....... nutrients to flow out of the forest nutrient storage pool for
!! ....... the woody components, cak - 03/11/2010


            !!if(forstg(j,iel)<0.001)  forstg(j,iel) = 0.001
            !!call nfix
           !! print*, "fixed N = ", fixn
         if (iel .eq. N ) forstg(j,iel) = forstg(j,iel) + fixn
        
          if (forstg(j,iel) .gt. 0.0) then
            amt = uptake(ESTOR,iel) * euf(LEAF)
            forstg(j,iel) = forstg(j,iel) - amt
            eleaf(j,iel) = eleaf(j,iel) + amt
      
            !!eupprt(LEAF,iel) = eupprt(LEAF,iel) + amt
            !!eupacc(iel) = eupacc(iel) + amt
            eprodfdy(iel) = eprodfdy(iel) + amt
            
            amt = uptake(ESTOR,iel) * euf(FROOTJ)
            forstg(j,iel) = forstg(j,iel) - amt
            efrootj(j,iel) = efrootj(j,iel) + amt
            !!eupprt(FROOT,iel) = eupprt(FROOT,iel) + amt
            !!eupacc(iel) = eupacc(iel) + amt
            eprodfdy(iel) = eprodfdy(iel) + amt
            
            amt = uptake(ESTOR,iel) * euf(FROOTM)
            forstg(j,iel) = forstg(j,iel) - amt
            efrootm(j,iel) = efrootm(j,iel) + amt
            !!eupprt(FROOT,iel) = eupprt(FROOT,iel) + amt
            !!eupacc(iel) = eupacc(iel) + amt
            eprodfdy(iel) = eprodfdy(iel) + amt
            
            if ((DECIDf(idf).eq. 1) .and. (growday .gt. FURGDYSf(idf))) then
              amt = uptake(ESTOR,iel) * euf(FBRCH) * (1.0 - FLSGRESf(idf))
              forstg(j,iel) = forstg(j,iel) - amt
              ebrch(j,iel) = ebrch(j,iel)+ amt
              !!eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
              !!eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
              
              amt = uptake(ESTOR,iel) * euf(LWOOD) * (1.0 - FLSGRESf(idf))
              forstg(j,iel) = forstg(j,iel) - amt
              elargw(j,iel) = elargw(j,iel) + amt
              !!eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
              !!eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
              
              amt = uptake(ESTOR,iel) * euf(CROOT) * (1.0 - FLSGRESf(idf))
              forstg(j,iel) = forstg(j,iel) - amt
              ecsroot(j,iel) = ecsroot(j,iel) + amt
              !!eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
              !!eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
            else
              amt = uptake(ESTOR,iel) * euf(FBRCH)
              forstg(j,iel) = forstg(j,iel) - amt
              ebrch(j,iel) = ebrch(j,iel)+ amt
             !! eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
              !!eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
              amt = uptake(ESTOR,iel) * euf(LWOOD)
              forstg(j,iel) = forstg(j,iel) - amt
              elargw(j,iel) = elargw(j,iel) + amt
              !! eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
              !!eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
              amt = uptake(ESTOR,iel) * euf(CROOT)
              forstg(j,iel) = forstg(j,iel) - amt
              ecsroot(j,iel) = ecsroot(j,iel) + amt
              !!eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
              !!eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
            endif
          endif

      !!  up  take from soil
      !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        
!!.......... Take up nutrients from soil
!!.......... Nutrients for uptake are available in the top tlaypg layers,
!!.......... cak 01/29/03
!!.......... If the tree is decidious and the time alloted to grow the
!!.......... woody components of the trees has passed flow nutrients to
!!.......... the forest storage pool rather than to the component nutrient
!!.......... pools based on the forest late season growth restriction
!!.......... parameter value, cak - 03/11/2010
          do 100 lyr = 1, tlaypgswat  !!   need to check wether this is the toplayer. here we only consider the first 6 layers
            if ((sol_no3(lyr,j)/10.+ sol_nh3(lyr,j)/10.) .gt. toler) then
              fsol = 1.0
!!.............. The fsol calculation for P is not needed here to compute
!!.............. the weighted average, cak - 04/05/02
!!              if (iel .eq. P) then
!!                fsol = fsfunc(minerl(SRFC,P), pslsrb, sorpmx)
!!              endif
                tno3nh3 = sol_no3(lyr,j)/10.+ sol_nh3(lyr,j)/10.
                fno3 = sol_no3(lyr,j)/10./ tno3nh3
                fnh3 = 1.0 - fno3

              calcup = uptake(ESOIL,iel)*tno3nh3 * fsol/availm(j,iel)
!!.............. Leaves
              amt = calcup * euf(LEAF)
              if (iel .eq. N) then
                sol_no3(lyr,j) = sol_no3(lyr,j) - 10. * amt*fno3  !! check if that is l or one
                sol_nh3(lyr,j) = sol_nh3(lyr,j) - 10. * amt*fnh3
              endif
              
              eleaf(j,iel) = eleaf(j,iel) + amt
              !!eupprt(LEAF,iel) = eupprt(LEAF,iel) + amt
              !!eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
              
!!.............. Juvenile fine roots
              amt = calcup * euf(FROOTJ)
              if (iel .eq. N) then
                sol_no3(lyr,j) = sol_no3(lyr,j) - 10. * amt*fno3  !! check if that is l or one
                sol_nh3(lyr,j) = sol_nh3(lyr,j) - 10. * amt*fnh3
              endif
              
              efrootj(j,iel) = efrootj(j,iel) + amt
              !!eupprt(FROOT,iel) = eupprt(FROOT,iel) + amt
              !!eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
!!.............. Mature fine roots
              amt = calcup * euf(FROOTM)
              if (iel .eq. N) then
                sol_no3(lyr,j) = sol_no3(lyr,j) - 10. * amt*fno3  !! check if that is l or one
                sol_nh3(lyr,j) = sol_nh3(lyr,j) - 10. * amt*fnh3
              endif
              
              efrootm(j,iel) = efrootm(j,iel) + amt
              !!eupprt(FROOT,iel) = eupprt(FROOT,iel) + amt
              !!eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
!!.............. Fine branch
              namt = 0.0
              if ((DECIDf(idf).eq. 1) .and. (growday .gt. FURGDYSf(idf))) then
                amt = calcup * euf(FBRCH) * FLSGRESf(idf)
                namt = namt + amt
                
                forstg(j,iel) = forstg(j,iel)+amt
                !!eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
                !!eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
                
                amt = calcup * euf(FBRCH) * (1.0 - FLSGRESf(idf))
                namt = namt + amt
                ebrch(j,iel) = ebrch(j,iel)+ amt
                
                !!eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
                !!eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
              else
                amt = calcup * euf(FBRCH)
                namt = namt + amt
                ebrch(j,iel) = ebrch(j,iel)+ amt
                !!eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
                !!eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
              endif
              if (iel .eq. N) then
                sol_no3(lyr,j) = sol_no3(lyr,j) - 10. * namt*fno3  !! check if that is l or one
                sol_nh3(lyr,j) = sol_nh3(lyr,j) - 10. * namt*fnh3
              endif
!!.............. Large wood
              namt = 0.0
              if ((DECIDf(idf).eq. 1) .and. (growday .gt. FURGDYSf(idf))) then
                amt = calcup * euf(LWOOD) * FLSGRESf(idf)
                namt = namt + amt
                !! call flow(minerl(lyr,iel),forstg(iel),time,amt)  !! check if amt was substracted from mineral pool when call flow
                !! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX be very careful here
                forstg(j,iel) = forstg(j,iel)+amt
                !!eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
                !!eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
                
                amt = calcup * euf(LWOOD) * (1.0 - FLSGRESf(idf))
                namt = namt + amt
                elargw(j,iel) = elargw(j,iel) + amt
                
                !!eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
                !!eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
              else
                amt = calcup * euf(LWOOD)
                namt = namt + amt
!!                call flow(minerl(lyr,iel),rlwode(iel),time,amt)
                elargw(j,iel) = elargw(j,iel) + amt
                !!eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
                !!eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
              endif
              if (iel .eq. N) then
                sol_no3(lyr,j) = sol_no3(lyr,j) - 10. * namt*fno3  !! check if that is l or one
                sol_nh3(lyr,j) = sol_nh3(lyr,j) - 10. * namt*fnh3
              endif
!!.............. Coarse roots
              namt = 0.0
              if ((DECIDf(idf).eq. 1) .and. (growday .gt. FURGDYSf(idf))) then
                amt = calcup * euf(CROOT) * FLSGRESf(idf)
                namt = namt + amt
                 forstg(j,iel) = forstg(j,iel)+amt  
                !!eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
                !!eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
                amt = calcup * euf(CROOT) * (1.0 - FLSGRESf(idf))
                namt = namt + amt
                ecsroot(j,iel) = ecsroot(j,iel) + amt
                
                !!eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
                !!eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
              else
                amt = calcup * euf(CROOT)
                namt = namt + amt
                !!call flow(minerl(lyr,iel),croote(iel),time,amt)
                ecsroot(j,iel) = ecsroot(j,iel) + amt
                
                !!eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
                !!eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
              endif
              if (iel .eq. N) then
                sol_no3(lyr,j) = sol_no3(lyr,j) - 10. * namt*fno3  !! check if that is l or one
                sol_nh3(lyr,j) = sol_nh3(lyr,j) - 10. * namt*fnh3
              endif
            endif
100       continue
      
      

      
      !!  up  take from N fixation
      !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
     !!....... Take up nutrients from nitrogen fixation
!!........... If the tree is decidious and the time alloted to grow the
!!........... woody components of the trees has passed flow nutrients to
!!........... the forest storage pool rather than to the component nutrient
!!........... pools based on the forest late season growth restriction
!!........... parameter value, cak - 03/11/2010
          if (iel .eq. N .and. treeNfix .gt. 0) then
!!............. Leaves
            amt = uptake(ENFIX,iel) * euf(LEAF)  !! check wether this equales fixn
            eleaf(j,iel) = eleaf(j,iel) + amt
            forstg(j,iel) = forstg(j,iel)-amt
            esrsnk(j,iel) = esrsnk(j,iel) - amt
            !!eupprt(LEAF,iel) = eupprt(LEAF,iel) + amt
            !!eupacc(iel) = eupacc(iel) + amt
            eprodfdy(iel) = eprodfdy(iel) + amt
!!............. Juvenile fine roots
            amt = uptake(ENFIX,iel) * euf(FROOTJ)
            efrootj(j,iel) = efrootj(j,iel) + amt
            forstg(j,iel) = forstg(j,iel)-amt
            esrsnk(j,iel) = esrsnk(j,iel) - amt
            !!eupprt(FROOT,iel) = eupprt(FROOT,iel) + amt
            !!eupacc(iel) = eupacc(iel) + amt
            eprodfdy(iel) = eprodfdy(iel) + amt
!!............. Mature fine roots
            amt = uptake(ENFIX,iel) * euf(FROOTM)
            efrootm(j,iel) = efrootm(j,iel) + amt
             forstg(j,iel) = forstg(j,iel)-amt
            esrsnk(j,iel) = esrsnk(j,iel) - amt
            !!eupprt(FROOT,iel) = eupprt(FROOT,iel) + amt
            !!eupacc(iel) = eupacc(iel) + amt
            eprodfdy(iel) = eprodfdy(iel) + amt
!!............. Fine branch
             if ((DECIDf(idf).eq. 1) .and. (growday .gt. FURGDYSf(idf))) then
              amt = uptake(ENFIX,iel) * euf(FBRCH) * FLSGRESf(idf)
              
              esrsnk(j,iel) = esrsnk(j,iel) - amt
              
              !!forstg(j,iel) = forstg(j,iel) + amt
              !!eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
              !!eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
              amt = uptake(ENFIX,iel) * euf(FBRCH) * (1.0 - FLSGRESf(idf))
              ebrch(j,iel) = ebrch(j,iel)+ amt
              forstg(j,iel) = forstg(j,iel) - amt
              esrsnk(j,iel) = esrsnk(j,iel) - amt
              !!eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
              !!eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
            else
              amt = uptake(ENFIX,iel) * euf(FBRCH)
              ebrch(j,iel) = ebrch(j,iel)+ amt
              forstg(j,iel) = forstg(j,iel) - amt
              esrsnk(j,iel) = esrsnk(j,iel) - amt
              !!eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
              !!eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
            endif
!!............. Large wood
             if ((DECIDf(idf).eq. 1) .and. (growday .gt. FURGDYSf(idf))) then
              amt = uptake(ENFIX,iel) * euf(LWOOD) * FLSGRESf(idf)
              
              esrsnk(j,iel) = esrsnk(j,iel) - amt
              !!forstg(j,iel) = forstg(j,iel) + amt
              !!eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
              !!eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
              amt = uptake(ENFIX,iel) * euf(LWOOD) * (1.0 - FLSGRESf(idf))
              elargw(j,iel) = elargw(j,iel) + amt
              forstg(j,iel) = forstg(j,iel) - amt
              esrsnk(j,iel) = esrsnk(j,iel) - amt
              !!eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
              !!eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
            else
              amt = uptake(ENFIX,iel) * euf(LWOOD)
              elargw(j,iel) = elargw(j,iel) + amt
              forstg(j,iel) = forstg(j,iel) - amt
              esrsnk(j,iel) = esrsnk(j,iel) - amt
              !!eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
              !!eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
            endif
!!............. Coarse roots
             if ((DECIDf(idf).eq. 1) .and. (growday .gt. FURGDYSf(idf))) then
              amt = uptake(ENFIX,iel) * euf(CROOT) * FLSGRESf(idf)
              
              esrsnk(j,iel) = esrsnk(j,iel) - amt
              !!forstg(j,iel) = forstg(j,iel) + amt
              !!eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
              !!eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
              amt = uptake(ENFIX,iel) * euf(CROOT) * (1.0 - FLSGRESf(idf))
              ecsroot(j,iel) = ecsroot(j,iel) + amt
              forstg(j,iel) = forstg(j,iel) - amt
              esrsnk(j,iel) = esrsnk(j,iel) - amt
              !!eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
              !!eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
            else
              amt = uptake(ENFIX,iel) * euf(CROOT)
              ecsroot(j,iel) = ecsroot(j,iel) + amt
              forstg(j,iel) = forstg(j,iel) - amt
              esrsnk(j,iel) = esrsnk(j,iel) - amt
              !!eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
              !!eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
            endif
          endif

     

      
!!....... Take up nutrients from automatic fertilizer
!!....... If the tree is decidious and the time alloted to grow the
!!....... woody components of the trees has passed flow nutrients to
!!....... the forest storage pool rather than to the component nutrient
!!....... pools based on the forest late season growth restriction
!!....... parameter value, cak - 03/11/2010
          if (aufert .ne. 0) then
            if (uptake(EFERT,iel) .gt. 0.) then
!!........... Automatic fertilizer added to plant pools
!!........... Leaves
              amt = uptake(EFERT,iel) * euf(LEAF)
              eleaf(j,iel) = eleaf(j,iel) + amt
              esrsnk(j,iel) = esrsnk(j,iel) - amt
              !!eupprt(LEAF,iel) = eupprt(LEAF,iel) + amt
              !!eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
!!........... Juvenile fine roots
              amt = uptake(EFERT,iel) * euf(FROOTJ)
              efrootj(j,iel) = efrootj(j,iel) + amt
              esrsnk(j,iel) = esrsnk(j,iel) - amt
              !!eupprt(FROOT,iel) = eupprt(FROOT,iel) + amt
              !!eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
!!........... Mature fine roots
              amt = uptake(EFERT,iel) * euf(FROOTM)
              efrootm(j,iel) = efrootm(j,iel) + amt
              esrsnk(j,iel) = esrsnk(j,iel) - amt
              !!eupprt(FROOT,iel) = eupprt(FROOT,iel) + amt
              !!eupacc(iel) = eupacc(iel) + amt
              eprodfdy(iel) = eprodfdy(iel) + amt
!!........... Fine branch
             if ((DECIDf(idf).eq. 1) .and. (growday .gt. FURGDYSf(idf))) then   !! check the original code to see if it is .gt. or .ge.
                amt = uptake(EFERT,iel) * euf(FBRCH) * FLSGRESf(idf)
                forstg(j,iel) = forstg(j,iel) + amt
                esrsnk(j,iel) = esrsnk(j,iel) - amt
                !!eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
                !!eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
                amt = uptake(EFERT,iel) * euf(FBRCH) * (1.0 - FLSGRESf(idf))
                ebrch(j,iel) = ebrch(j,iel)+ amt
                esrsnk(j,iel) = esrsnk(j,iel) - amt
                !!eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
                !!eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
              else
                amt = uptake(EFERT,iel) * euf(FBRCH)
                ebrch(j,iel) = ebrch(j,iel)+ amt
                esrsnk(j,iel) = esrsnk(j,iel) - amt
                !!eupprt(FBRCH,iel) = eupprt(FBRCH,iel) + amt
                !!eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
              endif
!!........... Large wood
             if ((DECIDf(idf).eq. 1) .and. (growday .gt. FURGDYSf(idf))) then
                amt = uptake(EFERT,iel) * euf(LWOOD) * FLSGRESf(idf)
                forstg(j,iel) = forstg(j,iel) + amt
                esrsnk(j,iel) = esrsnk(j,iel) - amt
                !!eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
                !!eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
                amt = uptake(EFERT,iel) * euf(LWOOD) * (1.0 - FLSGRESf(idf))
                elargw(j,iel) = elargw(j,iel) + amt
                esrsnk(j,iel) = esrsnk(j,iel) - amt
                !!eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
                !!eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
              else
                amt = uptake(EFERT,iel) * euf(LWOOD)
                elargw(j,iel) = elargw(j,iel) + amt
                esrsnk(j,iel) = esrsnk(j,iel) - amt
                !!eupprt(LWOOD,iel) = eupprt(LWOOD,iel) + amt
                !!eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
              endif
!!........... Coarse roots
             if ((DECIDf(idf).eq. 1) .and. (growday .gt. FURGDYSf(idf))) then
                amt = uptake(EFERT,iel) * euf(CROOT) * FLSGRESf(idf)
                forstg(j,iel) = forstg(j,iel) + amt
                esrsnk(j,iel) = esrsnk(j,iel) - amt
                !!eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
                !!eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
                amt = uptake(EFERT,iel) * euf(CROOT) * (1.0 - FLSGRESf(idf))
                ecsroot(j,iel) = ecsroot(j,iel) + amt
                esrsnk(j,iel) = esrsnk(j,iel) - amt
                !!eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
                !!eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
              else
                amt = uptake(EFERT,iel) * euf(CROOT)
                ecsroot(j,iel) = ecsroot(j,iel) + amt
                esrsnk(j,iel) = esrsnk(j,iel) - amt
                !!eupprt(CROOT,iel) = eupprt(CROOT,iel) + amt
                !!eupacc(iel) = eupacc(iel) + amt
                eprodfdy(iel) = eprodfdy(iel) + amt
              endif
  !!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX following processes are not checked yet         
!!........... Automatic fertilizer added to mineral pool
  !!XXXXXXX             if (favail(iel) .eq. 0.0) then
   !!XXXXXXX              write(*,*) 'Error in treegrow, favail(iel) = 0'
   !!XXXXXXX              STOP
    !!XXXXXXX           endif
    !!XXXXXXX           amt = uptake(EFERT,iel) * (1.0/favail(iel) - 1.0)
!!........... For now we are only accumulating fertilizer additions
!!........... for crop/grass, cak - 06/06/2008
!!              fertot(iel) = fertot(iel) + uptake(EFERT,iel) + amt
!!              fertac(iel) = fertac(iel) + uptake(EFERT,iel) + amt
     !!XXXXXXX          if (iel .eq. N) then
    !!XXXXXXX             lyr = SRFC
     !!XXXXXXX            call update_npool(lyr, amt, frac_nh4_fert,         &
     !!XXXXXXX                             frac_no3_fert, ammonium, nitrate, &
      !!XXXXXXX                            subname)
       !!XXXXXXX        endif
 !! XXXXX              call flow(esrsnk(j,iel),minerl(SRFC,iel),time,amt)  !! check check
         endif
         endif
110    continue


!!    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        

    else
            f_gpp(j) = 0.0
        do 140 iel = 1, MAXIEL
          eprodfdy(iel) = 0.0
          do 130 ipart = 1, FPARTS
            eup(j,ipart,iel) = 0.0
130       continue

140     continue
     
    end if        
    
    

    


!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !! calcuate the average temperature of the past 7 days (including current day)
  
         if (Mod((iyr-1),4) == 0) then 
          
          yeardays = 366
        else 
          yeardays = 365
        end if
    
    
    
     cdd = 0
     
     tavgp7d =0.
     tavgp7d = tavgp7d + tavgdate(j,curyr,(dofy))
  
  do dd = 1, 6
     if ((dofy-dd ) .gt.0) then
     
        tavgp7d = tavgp7d + tavgdate(j,curyr,(dofy-dd))
        cdd = cdd + 1
     else if (curyr > 1) then
        tavgp7d = tavgp7d + tavgdate(j,(curyr-1),(dofy-dd+yeardays))
        cdd = cdd + 1
     
       
     end if
  
  
  end do
  
  
  tavgp7d = tavgp7d / real(cdd+1)   !! need to chedk if this is right
 
         
!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
            do ilyr = 1, 2
            !! in daycent  swclimit[ilyr] = wiltpt[ilyr] - deltamin; deltamin is an input variable. here we assume it is 0.01cm h20
            
            sol_swclimit (ilyr) = sol_wpmm(ilyr,j)*0.1-0.01
            rel_wc(ilyr) = ((sol_st(ilyr, j)+sol_wpmm(ilyr,j))*0.1/(sol_z(ilyr,j)*0.1) - sol_swclimit (ilyr) ) /(sol_fc(ilyr, j)*0.1 + 0.01)     !! equation here are different from that in daycent becasue swat already minus wilting point in sol_st and sol_fc
            if (rel_wc(ilyr)<0.0) rel_wc(ilyr) = 0.
            if (rel_wc(ilyr)>1.0) rel_wc(ilyr) = 1.
            end do
            
            av_rel_wc = (rel_wc(1)*sol_z(1,j) +  rel_wc(2)*sol_z(2,j))/(sol_z(1,j) + sol_z(2,j))
            
           if (av_rel_wc > 1.0)then
           
            bgwfunc = 1.0
           else 
            bgwfunc = 1.0/(1.0 + 30.0 * exp(-9.0 * av_rel_wc))
          
           end if
!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                  

    !! weighted soil average temperature. Daycent only consider two soil layers

      avgstemp = (sol_tmp(1,j)*(sol_z(1,j))+sol_tmp(2,j)*(sol_z(2,j)))/ (sol_z(1,j)+sol_z(2,j))   !! check. layer depth in swat and daycent are different
!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX   
      
       
       call wdeath(tavgp7d, bgwfunc, tfrac, avgstemp, hrsinc,ldrmlt,mrspweff,tlaypgswat)

        
!!     return values to SWAT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       !! biomass
       bioday = 10.*f_npp(j)/0.42  !! convert unit from g c/m2 to kg/ha
       bio_ms(j) = bio_ms(j) +  bioday  !! need to double check how bio_ms(j) is calculated from different carbon pools
!!     !! plant height
       if (matyrsf(idf) > 0) then
              rto = float(curyr_mat(j)) / matyrsf(idf)
       else
              rto = 1.
       end if
        cht(j) = rto * chtmaxf(idf)
     
      !! lai
         
        laiday(j) = lai
      
      !!  dmaxleafc = maxleafc(dofy,j)!! closed temporarily change later
      !!  laimax = dmaxleafc * avg_slaf(idf)
      !!if (laiday(j) > laimax) laiday(j) = laimax   !! closed temporarily change later
        !!lai_yrmx(j) = lai
      !! if (laiday(j) > lai_yrmx(j)) laiday(j) = lai_yrmx(j) 
     
      !! phu accumulation
        fdelg = 0.
        if (phu_plt(j) > 0.1) then
          fdelg = (tmpav(j) - tbasef(idf)) / phu_plt(j)
        end if
        if (fdelg < 0.) fdelg = 0.
        phuacc(j) = phuacc(j) + fdelg
      !! ET
        if (phuacc(j) > 0.5 .and. phuacc(j) < 0.99) then
            plt_et(j) = plt_et(j) + ep_day + es_day
            plt_pet(j) = plt_pet(j) + pet_day
        end if
         
!!      laimxfr(j)  seems not used in SWAT
!!      olai(j) intermedia variable
!!       phuacc(:)
!!       plt_et(j) = plt_et(j) + ep_day + es_day !! in grow.f, this process is regulated by phenology
!!        plt_pet(j) = plt_pet(j) + pet_day
        if (bio_ms(j)>0) then
!!
           rwt(j) = (frootjc(j)+frootmc(j)+csrootc(j))/(bio_ms(j)*0.45)
        else
        
           rwt(j) = 0.
        endif   

      
      
      !!   check pool balance
      !!cal_temp(1) = n11
      !!cal_temp(2) = n22
      !!cal_temp(3) = n33
      !! cal_temp(4) = n44
      
      !! cal_temp(5) = n55
      !!cal_temp(6) = n66
     !! cal_temp(7) = n77
      !!cal_temp(8) = 0
     !! cal_temp(9) = 0
     
!! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
1099      return
      end subroutine

!! there should be Nup, Pup and nutriend allocation here



!!--------------------------------------------------------------------------------
      real Function photos(gs, gmin, tc, btrans, o2, tco2, fnf, ppf, atmp, jth)
!!       photosynthesis
       use parm
       !!real :: photos
       implicit none
       real, intent(in)::gs, gmin, tc,btrans,o2, tco2, fnf, ppf, atmp
       integer, intent(in):: jth
       real:: a,b,c,d,e,gn

       real:: a1, a2     !!a1:  Rubisco-limited net photosynthesis rates
                           !!a2:  light-limited photosynthesis rates

      !!  float kc25 = psnparam.kc25        !!co2 michaelis-menten constant at 25c (pa)
      !!  float ko25 = psnparam.ko25        !!o2 michaelis-menten constant at 25c (pa)
      !!  float ako = psnparam.ako          !!q10 for ko25
      !!  float akc = psnparam.akc          !!q10 for kc25
      !!  float vcmx25 = psnparam.vcmx25    !!maximum rate of carboxylation at 25c (umol co2/m**2/s)
      !!  float avcmx = psnparam.avcmx      !!q10 for vcmx25


        real :: kc         !!co2 michaelis-menten constant (pa)
        real :: ko         !!o2 michaelis-menten constant (pa)
        real :: k          !!enzyme kinetics
        real :: cp         !!co2 compensation point (pa)
        real :: vcmx       !!maximum rate of carboxylation (umol co2/m**2/s)
        real :: j          !!electron transport (umol co2/m**2/s)
        real :: jmax       !!light-saturated rate of electron transport
        real :: rd         !!daytime leaf dark respiration
        real :: fsw        !!soil water's effect on carbonxylation rate
        real :: ftc        !! foliar temperature's effect
        real :: sqa
        integer :: idf
        

          a = 0.
          b = 0.
          c = 0.
          d = 0.
          e = 0.
          sqa = 0.
                  
                  
        idf = idplt(jth)

        fsw = (btrans / 5.0) + 0.8  !! change later
          !!fsw = btrans
        ftc = 1.0 /(1. + exp( (-2.2e5+710.*(tc+273.16)) / (8.314*(tc+273.16))))

        gn = 1.e9 * gs/(8.3143*(tc+273.16))         !!!!m/s -> umol m-2 s1 pa-1

 !! need to check >=
        if(gn>=1.0e-6) then
                ko = ko25f(idf)*(akof(idf)**((tc-25.0)/10.))

                kc = kc25f(idf)*(akcf(idf)**((tc-25.0)/10.))
                k = kc * (1.+ o2/ko)

                cp = 0.000192*o2*(1.75**((tc-25.0)/10.))

                vcmx = vcmx25f(idf) * (avcmxf(idf)**((tc-25.0)/10.))  * ftc * fnf * fsw   !!!!maximum rate of carboxylation (umol co2/m**2/s

        
                jmax = 29.1 + 1.64*vcmx                !!!!Wullschleger, 1993

                j = jmax * ppf/(ppf+2.1*jmax)          !!!!Farquhar & von Caemmerer, 1982

                rd = 0.015 * vcmx

                if(vcmx <= 1.0e-6) then
                        a1 = 0.
                        a2 = 0.

                else

                        a = (k+tco2)*(k+tco2)
                        b = 2.*(2.*cp+k-tco2)*vcmx + 2.*(tco2+k)*rd
                        c = (vcmx-rd)*(vcmx-rd)
                        d = sqrt(a*gn*gn+b*gn+c)
                        e = sqrt(a*gmin*gmin+b*gmin+c)
                        sqa = sqrt(a)
                  
                  

                        a1 = 1.27/(2.0*(gn-gmin)) * (0.5*sqa*(gn*gn-gmin*gmin) & 
                        + sqrt(c)*(gn-gmin*gmin)- (2.*a*gn+b)*d/(4.*a)         &
                        + (2.*a*gmin+b)*sqrt(e)/(4.*a)                        &
                        + (b*b - 4.*a*c)*log((2.*a*gn+b+2.*sqa*d)/            &
                        (2.*a*gmin+b+2.*sqa*e))/(8.0*a**1.5))


                         a = (2.3*cp+tco2)*(2.3*cp+tco2)
                         b = 0.4*(4.3*cp-tco2)*j + 2.*(tco2+2.3*cp)*rd
                         c = (0.2*j-rd)*(0.2*j-rd)
                         d = sqrt(a*gn*gn+b*gn+c)
                         e = sqrt(a*gmin*gmin+b*gmin+c)
                         sqa = sqrt(a)

                         a2 = 1.27/(2.0*(gn-gmin)) * (0.5*sqa*(gn*gn-gmin*gmin)   &
                         + sqrt(c)*(gn-gmin*gmin)- (2.*a*gn+b)*d/(4.*a)           &
                         + (2.*a*gmin+b)*sqrt(e)/(4.*a)                           &
                         + (b*b - 4.*a*c)*log((2.*a*gn+b+2.*sqa*d)/               &
                         (2.*a*gmin+b+2.*sqa*e))/(8.0*a**1.5))
                    print*, "hahaha"

                 end if
         end if
         
           if(a1<=0.) a1 = 0.
           if(a2<=0.) a2 = 0.

       
           photos = min(a1,a2)
       

        if(ppf <=0. ) photos = 0.


      end function


!!--------------------------------------------------------------------------------
      real Function leaf_cn (biomass, tleafn, lai, laisun, jth, dummyn)

        use parm
        implicit none
        real, intent(in):: biomass
        real, intent(in):: tleafn
        real, intent(in):: lai
        real, intent(in):: laisun
        integer, intent(in):: jth, dummyn
        real :: fsun, sung, cn_sunleaf, cn_shaleaf
        real :: ncoef
        integer :: idf
        real :: shag


        real :: avgnconc        !1(gN/gBiomass)average n concentration in the leaf
        !! real :: leaf_cn

        idf = idplt(jth)

        ncoef = exp(-ext_nf(idf)*0.5)
        if(biomass>0.) avgnconc = tleafn/biomass

        if(lai>=1.0e-12 .and. biomass >= 1.0e-12) then
                fsun = laisun/lai
                sung = avgnconc/((fsun+(1.-fsun)*ncoef))

                cn_sunleaf = 0.45 / sung
                cn_shaleaf = 0.45 / (sung*ncoef)

        else

                cn_sunleaf = -1.
                cn_shaleaf = -1.
        end if
        
        if (dummyn<2) then

            leaf_cn = cn_sunleaf
        else 
            leaf_cn = cn_shaleaf
        end if
        
       end function
!!--------------------------------------------------------------------------------
!!       maintenance respiration
!!--------------------------------------------------------------------------------
!!      phenology stage

!!----------------------------------------------------------------------------------------


!!     carbon allocation
!!--------------------------------------------------------------------------------
!!     alloc
!!--------------------------------------------------------------------------------





